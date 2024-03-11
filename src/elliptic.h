/**
   PV related routines:

   - Invert elliptic equation: del^2 psi + Gamma psi = q
   - and compute relative vorticity: omega = del^2 psi

    TODO:
      - solve one direction with tridiagonal solver (Thomas algorithm)
 */


// FFTW in/outputs and plans
double *wrk1;
fftw_complex *wrk1_alias;
fftw_plan transfo_direct, transfo_inverse;

// normalisation factor
double norm;

#include "eigmode.h"

#define idx_fft(i,j) (j-1)*N_fft + (i-1)

#ifdef _MPI
  #define idx_fft2(i,j) (i)*(2*Nx) + 2*(j)
#else
  #define idx_fft2(i,j) (j)*Nxp2 + 2*(i)
#endif

void init_elliptic(){

  /* create FFT plans */
  if (bc_fac == -1) { // periodic

    N_fft = Nxp2;
    norm = NX*NY;

    wrk1 = fftw_alloc_real(2*alloc_local);
    wrk1_alias = (fftw_complex*) &wrk1[0];

    #ifdef _MPI
      transfo_direct = fftw_mpi_plan_dft_r2c_2d(NY, NX, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_OUT);
      transfo_inverse = fftw_mpi_plan_dft_c2r_2d(NY, NX, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
    #else
      transfo_direct  = fftw_plan_dft_r2c_2d(Nx, Ny, wrk1, wrk1_alias, FFTW_EXHAUSTIVE);
      transfo_inverse = fftw_plan_dft_c2r_2d(Nx, Ny, wrk1_alias, wrk1, FFTW_EXHAUSTIVE);
    #endif
  } else {

    N_fft = Nx;
    norm = 4*NXp1*NYp1;

    wrk1 = fftw_alloc_real(alloc_local);

    #ifdef _MPI
      transfo_direct = fftw_mpi_plan_r2r_2d(NY, NX, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_RODFT00, FFTW_RODFT00,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_OUT);
      transfo_inverse = fftw_mpi_plan_r2r_2d(NY, NX, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_RODFT00, FFTW_RODFT00, 
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
    #else
      transfo_direct  = fftw_plan_r2r_2d(Ny,Nx, wrk1, wrk1,
                                        FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
      transfo_inverse = fftw_plan_r2r_2d(Ny,Nx, wrk1, wrk1,
                                        FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    #endif
  }

  
  // Prepare vertical mode inversion
  init_eigmode();
  compute_eigmode();
}

void invert_pv(double *q, double *psi, double *omega) {

  // reset psi. Temporary: only needed if doing FFT mode by mode
  for(int k = 0; k<nl; k++){
    for(int j = 0; j<Nyp2; j++){
      for(int i = 0;i <Nxp2; i++){
        psi[idx(i,j,k)] = 0.;
      }
    }
  }
  
  // loop on modes
  for(int k = 0; k<nl; k++){
  
    for(int j = 1; j<Nyp1; j++){
      for(int i = 1; i <Nxp1; i++){
        wrk1[idx_fft(i,j)] = 0;
        // inner loop on layers: projection on modes
        for (int l = 0; l < nl ; l++) {
          wrk1[idx_fft(i,j)] += cl2m[k*nl+l]*q[idx(i,j,l)];
        }
      }
    }

    // direct transform
    fftw_execute(transfo_direct);

    if (bc_fac == -1){ // periodic boundary conditions

      for(int j = 0; j < Nx; j++){
        for (int i = 0; i < Nk; i++){
          double fact = - (sq(L[j]) + sq(K[i]) + iRd2[k]);
          if (fact != 0){
            wrk1[idx_fft2(i,j)] = wrk1[idx_fft2(i,j)]/fact;
            wrk1[idx_fft2(i,j) + 1] = wrk1[idx_fft2(i,j) + 1]/fact;
          }
        }
      }

      if (rank == 0) {
        wrk1[0] = 0;
        wrk1[1] = 0;
      }

    } else { // bounded conditions

      for(int j = 1; j < Nyp1; j++){ // f = -g / ( kx^2 + ky^2 + 1/Rd^2)
        for (int i = 1; i < Nxp1; i++){
          double fact = - (sq(L[i]) + sq(K[j]) + iRd2[k]);
          wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/fact;
        }
      }
    }

    // inverse transform
    fftw_execute(transfo_inverse);
    
    // Normalisation
    for(int j = 1; j<Nyp1; j++){
      for(int i = 1;i <Nxp1; i++){
        wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/norm;
      }
    }

    for(int j = 1;j<Nyp1; j++){
      for(int i = 1;i <Nxp1; i++){
        // loop on layer: populating psi (layer l)
        for (int l = 0; l < nl ; l++) {
          psi[idx(i,j,l)] += cm2l[l*nl+k]*wrk1[idx_fft(i,j)];
        }
      }
    }
  } // mode loop

  // adjust boundary conditions
  // in periodic, it is important to set the psi BC before q and omega
  adjust_bc(psi, 0, psi);
  adjust_bc(q,   1, psi);

  // calculate omega
  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Nyp1; j++){
      for(int i = 1; i<Nxp1; i++){
        omega[idx(i,j,k)] = laplacian(psi);
      }
    }
  }
  adjust_bc(omega, 1, psi);

}

void clean_fft(){

  fftw_free(wrk1);
  
  fftw_destroy_plan(transfo_direct); 
  fftw_destroy_plan(transfo_inverse);

  #ifdef _MPI
    fftw_mpi_cleanup();
  #else
    fftw_cleanup();
  #endif
}
