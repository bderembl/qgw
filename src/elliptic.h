/**
   Invert elliptic equation
    del^2 psi = q

    TODO:
      - MPI routines
      - solve one direction with tridiagonal solver (Thomas algorithm)
      - Multi layer
 */


#include <fftw3.h>

#include "eigmode.h"

double *in1;
double *in2; 
double *out1; 
double *out2; 

fftw_plan transfo_direct, transfo_inverse;

void init_elliptic(){
  fprintf(stdout,"Prepare fft..");

  in1  = calloc( Nxm1*Nym1, sizeof( double ) );
  in2  = calloc( Nxm1*Nym1, sizeof( double ) );
  out1 = calloc( Nxm1*Nym1, sizeof( double ) );
  out2 = calloc( Nxm1*Nym1, sizeof( double ) );

  transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, in1, out1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
  transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, in2, out2, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

  fprintf(stdout,"..done\n");

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
    int ll = 0;
    for(int j = 1; j<Ny; j++){
      for(int i = 1; i <Nx; i++){
        in1[ll] = 0;
        // inner loop on layers: projection on modes
        for (int l = 0; l < nl ; l++) {
          //        in1[idx_in(i,j)] = q[idx(i,j,k)];
          in1[ll] += cl2m[k*nl+l]*q[idx(i,j,l)];
        }
        ll += 1;
      }
    }
    
  // direct transform
  fftw_execute(transfo_direct);

  // solve elliptic in fourrier space
  ll = 0;
    for(int j = 1; j < Ny; j++){ // f = -g / ( kx^2 + ky^2 )
      for (int i = 1; i < Nx; i++){

        double fact = - (sq(i*pi/Lx) + sq(j*pi/Ly) + iRd2[k]);
        in2[ll] = out1[ll]/fact;
        ll += 1;
      }
    }

  // inverse transform
  fftw_execute(transfo_inverse);

  // Scaling
  ll = 0;
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
//        out2[idx_fft(i,j)] = out2[idx_fft(i,j)]/(4*(Nxm1 + 1)*(Nym1 + 1));
        out2[ll] = out2[ll]/(4*(Nxm1 + 1)*(Nym1 + 1));
        ll += 1;
      }
    }

  ll = 0;
    for(int j = 1;j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        // loop on layer: populating psi (layer l)
        for (int l = 0; l < nl ; l++) {
          psi[idx(i,j,l)] += cm2l[l*nl+k]*out2[ll];
        }
        ll += 1;
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
  fftw_cleanup();
}
