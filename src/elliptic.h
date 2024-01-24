/**
   Invert elliptic equation
    del^2 psi = q

    TODO:
      - solve one direction with tridiagonal solver (Thomas algorithm)
 */


// FFTW in/outputs and plans
double *wrk1;

fftw_plan transfo_direct, transfo_inverse;

// MPI FFTW
ptrdiff_t alloc_local, local_n0, local_0_start;


#include "eigmode.h"

#define idx_fft(i,j) (j-1)*Nxm1 + (i-1)

void init_elliptic(){
  
  
  /**
     Prepare local indices
  */

  Nx = NX;
  Ny = NY;

  Nxm1 = Nx - 1;
  Nym1 = Ny - 1;

  Nxp1 = Nx + 1;
  Nyp1 = Ny + 1;

  NXp1 = NX + 1;
  NYp1 = NY + 1;

#ifdef _MPI

  long int fft_dims[2];
  fft_dims[0] = Nym1;
  fft_dims[1] = Nxm1;
  
  int block0 = (int) NY/n_ranks;

  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_many(2, &fft_dims[0],1,block0,MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);

  wrk1 = fftw_alloc_real(alloc_local);

  /* create plan for out-of-place */
  transfo_direct = fftw_mpi_plan_r2r_2d(Nym1, Nxm1, wrk1, wrk1, MPI_COMM_WORLD,
                                        FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_OUT);
  transfo_inverse = fftw_mpi_plan_r2r_2d(Nym1, Nxm1, wrk1, wrk1, MPI_COMM_WORLD,
                                         FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
  
  J0 = local_0_start + 1; // member the fourier grid starts at the index 1 of the real space grid
  Ny = local_n0 + 1;
  Nym1 = local_n0;
  Nyp1 = local_n0 + 2;
  
#else

  J0 = 1;
  alloc_local = Nxm1*Nym1;

  wrk1 = fftw_alloc_real(alloc_local);

  
  transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, wrk1, wrk1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
  transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, wrk1, wrk1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
  
#endif

  
  // Prepare vertical mode inversion
  init_eigmode();
  compute_eigmode();
}

void invert_pv(double *q, double *psi) {

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

    // solve elliptic in fourrier space
    for(int j = 1; j < Ny; j++){ // f = -g / ( kx^2 + ky^2 + Rd^2)
      for (int i = 1; i < Nx; i++){
        double fact = - (sq(L[i]) + sq(K[j]) + iRd2[k]);
        wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/fact;
      }
    }

    // inverse transform
    fftw_execute(transfo_inverse);

    // Normalisation
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        wrk1[idx_fft(i,j)] = wrk1[idx_fft(i,j)]/(4*(Nxm1 + 1)*NY);
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
