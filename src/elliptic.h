/**
   PV related routines:

   - Invert elliptic equation: del^2 psi + Gamma psi = q
   - and compute relative vorticity: omega = del^2 psi

    TODO:
      - solve one direction with tridiagonal solver (Thomas algorithm)
 */


// FFTW in/outputs and plans
double *wrk1;
fftw_complex *wrk2;
fftw_plan transfo_direct, transfo_inverse;

// normalisation factor
double norm;

#include "eigmode.h"

#define idx_fft(i,j) (j-1)*N_fft + (i-1)

#ifdef _MPI
  #define idx_fft2(i,j) (i)*(2*Nxm1) + 2*(j)
#else
  #define idx_fft2(i,j) (j)*Nxp1 + 2*(i)
#endif

void init_elliptic(){

  /* create FFT plans */
  if (bc_fac == -1) { // periodic

    N_fft = Nxp1;
    norm = NXm1*NYm1;

    wrk1 = fftw_alloc_real(2*alloc_local);
  
    #ifdef _MPI
      transfo_direct = fftw_mpi_plan_dft_r2c_2d(NYm1, NXm1, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_OUT);
      transfo_inverse = fftw_mpi_plan_dft_c2r_2d(NYm1, NXm1, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
    #else
      transfo_direct  = fftw_plan_dft_r2c_2d(Nxm1, Nym1, wrk1, wrk1, FFTW_EXHAUSTIVE);
      transfo_inverse = fftw_plan_dft_c2r_2d(Nxm1, Nym1, wrk1, wrk1, FFTW_EXHAUSTIVE);
    #endif
  } else {

    N_fft = Nxm1;
    norm = 4*NX*NY;

    wrk1 = fftw_alloc_real(alloc_local);

    #ifdef _MPI
      transfo_direct = fftw_mpi_plan_r2r_2d(NYm1, NXm1, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_RODFT00, FFTW_RODFT00,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_OUT);
      transfo_inverse = fftw_mpi_plan_r2r_2d(NYm1, NXm1, wrk1, wrk1, MPI_COMM_WORLD,
                                            FFTW_RODFT00, FFTW_RODFT00, 
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
    #else
      transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, wrk1, wrk1,
                                        FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
      transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, wrk1, wrk1,
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
    for(int j = 0; j<Nyp1; j++){
      for(int i = 0;i <Nxp1; i++){
        psi[idx(i,j,k)] = 0.;
      }
    }
  }
  
  // loop on modes
  for(int k = 0; k<nl; k++){
  
    for(int j = 1; j<Ny; j++){
      for(int i = 1; i <Nx; i++){
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

      for(int j = 0; j < Nxm1; j++){
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

      for(int j = 1; j < Ny; j++){ // f = -g / ( kx^2 + ky^2 + 1/Rd^2)
        for (int i = 1; i < Nx; i++){
          double fact = - (sq(L[i]) + sq(K[j]) + iRd2[k]);
          wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/fact;
        }
      }
    }

    // inverse transform
    fftw_execute(transfo_inverse);
    
    // Normalisation
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/norm;
      }
    }

    for(int j = 1;j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        // loop on layer: populating psi (layer l)
        for (int l = 0; l < nl ; l++) {
          psi[idx(i,j,l)] += cm2l[l*nl+k]*wrk1[idx_fft(i,j)];
        }
      }
    }
  } // mode loop

  // adjust boundary conditions
  adjust_bc(q, psi);

  // calculate omega
  calc_omega(psi, omega);

}

void clean_fft(){

  fftw_free(wrk1);
  if (bc_fac == -1){
    fftw_free(wrk2);
  }
  
  fftw_destroy_plan(transfo_direct); 
  fftw_destroy_plan(transfo_inverse);

  #ifdef _MPI
    fftw_mpi_cleanup();
  #else
    fftw_cleanup();
  #endif
}
