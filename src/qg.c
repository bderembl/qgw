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
#define min(p,q) p > q ? q : p
#define max(p,q) p < q ? q : p
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
int N_c;
int nl = 1;
double Lx, Ly;
double Delta;
double t = 0;
double dt = 0.;
double tend = 1;
double dt_out = 0.1;
double t_out = 0;
double dt_print = 0.1;
double t_print = 0;
double cfl = 0.2;
double DT_max = 0;
int it = 0;

// physical constants and functions
double beta = 0.;
double nu = 0.;
double tau0 = 0.;
double bc_fac = 0.;

// grid indices
#define idx(i,j) (j)*Nxp1 + (i)

// define forcing
double eps = 1.;
double k_forc = 0.1;


#include "domain.h"
#include "elliptic.h"
#include "forcing.h"
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
  params = list_append(params, &dt_out, "dt_out", "double");
  params = list_append(params, &eps, "eps", "double");
  params = list_append(params, &k_forc, "k_forc", "double");
  params = list_append(params, &dt_print, "dt_print", "double");
  params = list_append(params, &cfl, "cfl", "double");

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

  #ifdef _STOCHASTIC
    init_stoch_forc();
    printf("Stochastic forcing. \n");
  #else
    init_det_forc();
    printf("Large-scale forcing. \n");
  #endif

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
  while(t < tend){
    if (fabs (t - t_out) < TEPS*dt){
      printf("Write output, t = %e \n",t);
      t_out += dt_out;
      write_nc();
    }

    if (t > t_print){ // only approximate here
      fprintf(stdout, "i = %d, t = %e dt = %e \n",it, t, dt);
      t_print += dt_print;
    }

    timestep(q);
    it ++;
  }

 
  /**
     Cleanup
  */
	
  #ifdef _STOCHASTIC
    clean_stoch_forcing();
  #else
    clean_det_forcing();
  #endif

  clean_fft();
  clean_timestep();
	
  free(psi);
  free(q);
  free(X);
  free(Y);
  if (list_nc) list_free(list_nc);
  if (params) list_free(params);

}
