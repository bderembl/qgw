/**
   Invert elliptic equation
    del^2 psi = q

    TODO:
      - MPI routines
      - solve one direction with tridiagonal solver (Thomas algorithm)
      - Multi layer
 */


#include <fftw3.h>

double *in1;
double *in2; 
double *out1; 
double *out2; 

fftw_plan transfo_direct, transfo_inverse;

void init_fft(){
  fprintf(stdout,"Prepare fft..");

  in1  = calloc( Nxm1*Nym1, sizeof( double ) );
  in2  = calloc( Nxm1*Nym1, sizeof( double ) );
  out1 = calloc( Nxm1*Nym1, sizeof( double ) );
  out2 = calloc( Nxm1*Nym1, sizeof( double ) );

  transfo_direct  = fftw_plan_r2r_2d(Nym1,Nxm1, in1, out1, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
  transfo_inverse = fftw_plan_r2r_2d(Nym1,Nxm1, in2, out2, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

  fprintf(stdout,"..done\n");

}

void invert_pv(double *q, double *psi) {


  int l = 0;

  for(int j = 1; j<Ny; j++){
    for(int i = 1; i <Nx; i++){
      in1[l] = q[idx(i,j)];
      l=l+1;
    }
  }

  // direct transform
  fftw_execute(transfo_direct);

  // solve elliptic in fourrier space
  l = 0;
  for(int j = 1; j < Ny; j++){ // f = -g / ( kx^2 + ky^2 )
    for (int i = 1; i < Nx; i++){ 

      double fact = - (sq(i*pi/Lx) + sq(j*pi/Ly));
      in2[l] = out1[l]/fact;
      l=l+1;

    }
  }

  // inverse transform
  fftw_execute(transfo_inverse);

  // Scaling
  l = 0;
  for(int j = 1; j<Ny; j++){
    for(int i = 1;i <Nx; i++){
      out2[l] = out2[l]/(4*(Nxm1 + 1)*(Nym1 + 1));
      l=l+1;
    }
  }

  l = 0;
  for(int j = 1;j<Ny; j++){
    for(int i = 1;i <Nx; i++){
      psi[idx(i,j)] = out2[l]; 
      l=l+1;
    }
  }

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
