/**
   Invert elliptic equation
    del^2 psi = q

    TODO:
      - solve one direction with tridiagonal solver (Thomas algorithm)
 */


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
  
  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_2d(Nym1, Nxm1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  
  in1 = fftw_alloc_real(alloc_local);
  in2 = fftw_alloc_real(alloc_local);
  out1 = fftw_alloc_real(alloc_local);
  out2 = fftw_alloc_real(alloc_local);
  
  /* create plan for out-of-place */
  transfo_direct = fftw_mpi_plan_r2r_2d(Nym1, Nxm1, in1, out1, MPI_COMM_WORLD,
                                        FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_DESTROY_INPUT|FFTW_MPI_TRANSPOSED_OUT);
  transfo_inverse = fftw_mpi_plan_r2r_2d(Nym1, Nxm1, in2, out2, MPI_COMM_WORLD,
                                         FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_DESTROY_INPUT|FFTW_MPI_TRANSPOSED_IN);
  
  J0 = local_0_start + 1; // member the fourier grid starts at the index 1 of the real space grid
  Ny = local_n0 + 1;
  Nym1 = local_n0;
  Nyp1 = local_n0 + 2;
  
#else

  J0 = 1;
  
  in1  = calloc( Nxm1*Nym1, sizeof( double ) );
  in2  = calloc( Nxm1*Nym1, sizeof( double ) );
  out1 = calloc( Nxm1*Nym1, sizeof( double ) );
  out2 = calloc( Nxm1*Nym1, sizeof( double ) );
  
  transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, in1, out1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_DESTROY_INPUT);
  transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, in2, out2, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE|FFTW_DESTROY_INPUT);
  
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
        in1[idx_fft(i,j)] = 0;
        // inner loop on layers: projection on modes
        for (int l = 0; l < nl ; l++) {
          in1[idx_fft(i,j)] += cl2m[k*nl+l]*q[idx(i,j,l)];
        }
      }
    }
    
    // direct transform
    fftw_execute(transfo_direct);

    // solve elliptic in fourrier space
    for(int j = 1; j < Ny; j++){ // f = -g / ( kx^2 + ky^2 + Rd^2)
      for (int i = 1; i < Nx; i++){
        double fact = - (sq(L[i]) + sq(K[j]) + iRd2[k]);
        in2[idx_fft(i,j)] = out1[idx_fft(i,j)]/fact;
      }
    }

    // inverse transform
    fftw_execute(transfo_inverse);

    // Normalisation (different for MPI as we need to take the global Nyt)
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        out2[idx_fft(i,j)] = out2[idx_fft(i,j)]/(4*(Nxm1 + 1)*NY);
      }
    }

    for(int j = 1;j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        // loop on layer: populating psi (layer l)
        for (int l = 0; l < nl ; l++) {
          psi[idx(i,j,l)] += cm2l[l*nl+k]*out2[idx_fft(i,j)];
        }
      }
    }
  } // mode loop

  // adjust boundary conditions
  adjust_bc(q, psi);

}

void clean_fft(){

  free(in1);
  free(in2); 
  free(out1); 
  free(out2); 
  
  fftw_destroy_plan(transfo_direct); 
  fftw_destroy_plan(transfo_inverse);

  #ifdef _MPI
    fftw_mpi_cleanup();
  #else
    fftw_cleanup();
  #endif
}
