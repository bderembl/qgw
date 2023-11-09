/**
   QG model 

   Compile with 
     gcc -O3 -Wall qg.c -o qg.e -lm -lfftw3 -lnetcdf

   TODO
     - Documentation
     - Test cases
     - MPI
     - GPU
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "extra.h"

#define sq(x) ((x)*(x)) // alias for square function
double pi = 3.141592653589793;

// field variables
double *psi;
double *q;
double *X;
double *Y;

List *params; 

// space and time constants
int Nx, Ny;
int Nxm1, Nym1;
int Nxp1, Nyp1;
int nl = 1;
double Lx, Ly;
double Delta;
double t = 0;
double dt = 0.;
double tend = 1;
int it = 0;
int it_out = 500;

// physical constants and functions
double beta = 0.;
double nu = 0.;
double tau0 = 0.;
double bc_fac = 0.;
#define forcing_q(t) ( -tau0/Ly*pi*sin(pi*Y[j]/Ly))

#include "domain.h"
#include "elliptic.h"
#include "dynamics.h"
#include "timestep.h"
#include "netcdf_io.h"


int main(int argc,char* argv[])
{
  
  /**
     Namelist and parameters
   */

  // add here variables to be read in the input file
  params = list_append(params, &Nx, "Nx", "int");
  params = list_append(params, &Ny, "Ny", "int");
  params = list_append(params, &Lx, "Lx", "double");
  params = list_append(params, &tau0, "tau0", "double");
  params = list_append(params, &beta, "beta", "double");
  params = list_append(params, &nu, "nu", "double");
  params = list_append(params, &bc_fac, "bc_fac", "double");
  params = list_append(params, &dt, "dt", "double");
  params = list_append(params, &tend, "tend", "double");
  params = list_append(params, &it_out, "it_out", "int");

  // Search for the configuration file with a given path or read params.in 
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(params, file_param);
  create_outdir();
  backup_config(file_param);
  
  /**
     Initialization
   */

  init_domain();
  init_vars();
  init_fft();
  init_timestep();

  invert_pv(q,psi);

  list_nc = list_append(list_nc, psi,"psi", "double");
  list_nc = list_append(list_nc, q, "q", "double");
  char file_tmp[90];
  sprintf (file_tmp,"%s%s", dir_out, "vars.nc");
  //create_nc("vars.nc");
  create_nc(file_tmp);

  /**
     Main Loop
  */
  while(t<tend){
    if (it%it_out == 0){
      printf("it = %d \n",it);
      write_nc();
    }
    timestep(q);
    it ++;
  }

 
  /**
     Cleanup
  */

  clean_fft();
  clean_timestep();

  free(psi);
  free(q);
  free(X);
  free(Y);
  if (list_nc) list_free(list_nc);

}
