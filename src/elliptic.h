/**
   Invert elliptic equation
    del^2 psi = q

    TODO:
      - MPI routines
      - solve one direction with tridiagonal solver (Thomas algorithm)
      - Multi layer
 */


#include "eigmode.h"

#define idx_fft(i,j) (j-1)*Nxm1 + (i-1)

void init_elliptic(){
  
  Nxm1 = Nx - 1;
  Nym1 = Ny - 1;

  Nxp1 = Nx + 1;
  Nyp1 = Ny + 1;

  #ifdef _MPI

    NY = Nym1;
    NX = Nxm1;

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(NY, NX, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    in1 = fftw_alloc_real(alloc_local);
    in2 = fftw_alloc_real(alloc_local);
    out1 = fftw_alloc_real(alloc_local);
    out2 = fftw_alloc_real(alloc_local);
    
    /* create plan for out-of-place */
    transfo_direct = fftw_mpi_plan_r2r_2d(NY, NX, in1, out1, MPI_COMM_WORLD,
                                FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    transfo_inverse = fftw_mpi_plan_r2r_2d(NY, NX, in2, out2, MPI_COMM_WORLD,
                                FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

    /* From now on Ny and all related variables will be the local values, and the global 
    values will be stored in Nt.*/

    Nyt = Ny;
    Nytp1 = Nyp1;
    Nytm1 = Nym1;
    Ny_start = local_0_start + 1;
    Ny_startm1 = local_0_start;
    Nym1 = local_n0;
    Nyp1 = local_n0+2;
    Ny = local_n0+1;
    
  #else
    
    fprintf(stdout,"Prepare fft..");

    in1  = calloc( Nxm1*Nym1, sizeof( double ) );
    in2  = calloc( Nxm1*Nym1, sizeof( double ) );
    out1 = calloc( Nxm1*Nym1, sizeof( double ) );
    out2 = calloc( Nxm1*Nym1, sizeof( double ) );

    transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, in1, out1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, in2, out2, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

    fprintf(stdout,"..done\n");

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
        double fact = - (sq(K[i]) + sq(L[j]) + iRd2[k]);
        in2[idx_fft(i,j)] = out1[idx_fft(i,j)]/fact;
      }
    }

  // inverse transform
  fftw_execute(transfo_inverse);

  // Scaling (different for MPI as we need to take the global Ny)
  #ifdef _MPI
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        out2[idx_fft(i,j)] = out2[idx_fft(i,j)]/(4*(Nxm1 + 1)*(Nytm1 + 1));
      }
    }
  #else
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        out2[idx_fft(i,j)] = out2[idx_fft(i,j)]/(4*(Nxm1 + 1)*(Nym1 + 1));
      }
    }
  #endif

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
    void fftw_mpi_cleanup(void);
  #else
    fftw_cleanup();
  #endif
}
