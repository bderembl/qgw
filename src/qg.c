/**
   QG model 

   Compile with 
     gcc -O3 -Wall qg.c -o qg.e -lm -lfftw3 -llapacke -lnetcdf 
   If run with MPI compile with 
    mpicc -D_MPI -O3 qg.c -o qg.e -lfftw3_mpi -lfftw3 -lm -llapacke -lnetcdf

   Compilation flags
     -D_STOCHASTIC : add a stochastic forcing
     -D_MPI : needed for mpi compilation
    
   Run with
     ./qg.e
   or with
     mpirun -n NPROC qg.e

     create a restart file:
     ncks -d time,-1,-1 vars.nc restart.nc


   TODO
     - Documentation
     - Test cases
     - GPU
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#ifdef _MPI
  #include <mpi.h>
  #include <fftw3-mpi.h>
#endif

#define sq(x) ((x)*(x)) // alias for square function
#define min(p,q) p > q ? q : p
#define max(p,q) p < q ? q : p
#define sign(x) ((x) > 0 ? 1 : -1)
#define HUGE 1e30
double pi = 3.141592653589793;
#define nl_max 1000

// field variables
double *psi;
double *q;
double *topo;
double *omega;

// variable for printing out intermediate initialisation info
int print = 1;

// Physical Parameters
int nl = 1;
double dh[nl_max] = {1.};
double dhc[nl_max] = {1.};
double Lx, Ly;
double Delta;
double t = 0;
double dt = 0.;
double tend = 1;
double dt_out = 0.1;
double t_out = 0;
double cfl = 0.2;
double DT_max = 0;
int it = 0;
double beta = 0.;
double nu = 0.;
double nu_kin = 0.;
double hEkb = 0.;
double tau0 = 0.;
double forc_mode = 1.0;
double f0 = 1.e-5;
double bc_fac = 0.;
double N2[nl_max] = {1.};
double Ld = 0.;
double h_topo = 0.;
double w_topo = 1.;

// Local header files

#include "mpi_utils.h"
#include "extra.h"

// declaration of list type needs to occur after extra.h import
List *params; 

#include "domain.h"
#include "elliptic.h"
#include "forcing.h"
#include "dynamics.h"
#include "timestep.h"
#include "netcdf_io.h"

int main(int argc,char* argv[])
{

  init_mpi();

  /**
     Namelist and parameters
   */

  // add here variables to be read in the input file
  params = list_append(params, &NX, "NX", "int");
  params = list_append(params, &NY, "NY", "int");
  params = list_append(params, &nl, "nl", "int");
  params = list_append(params, &Lx, "Lx", "double");
  params = list_append(params, &dh, "dh", "array");
  params = list_append(params, &tau0, "tau0", "double");
  params = list_append(params, &forc_mode, "forc_mode", "double");
  params = list_append(params, &f0, "f0", "double");
  params = list_append(params, &beta, "beta", "double");
  params = list_append(params, &nu, "nu", "double");
  params = list_append(params, &nu_kin, "nu_kin", "double");
  params = list_append(params, &hEkb, "hEkb", "double");
  params = list_append(params, &N2, "N2", "array");
  params = list_append(params, &Ld, "Ld", "double");
  params = list_append(params, &h_topo, "h_topo", "double");
  params = list_append(params, &w_topo, "w_topo", "double");
  params = list_append(params, &bc_fac, "bc_fac", "double");
  params = list_append(params, &tend, "tend", "double");
  params = list_append(params, &dt_out, "dt_out", "double");
  params = list_append(params, &sigma_f, "sigma_f", "double");
  params = list_append(params, &k_f, "k_f", "double");
  params = list_append(params, &cfl, "cfl", "double");

  // Search for the configuration file with a given path or read params.in 
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  // TODO: Only rank 0 reads, then broadcasts params
  read_params(params, file_param);
  
  create_outdir();
  backup_config(file_param);
  
  /**
     Initialization
   */

  init_domain();
  init_elliptic();
  init_vars();
  init_timestep();

#ifdef _STOCHASTIC
  init_stoch_forc();
  fprintf(stdout,"Stochastic forcing. \n");
#endif

  // read q0 (restart)
  read_nc("restart.nc");
  // First inversion
  invert_pv(q, psi, omega);

  // Initialize output
  char file_tmp[90];
  sprintf (file_tmp,"%s%s", dir_out, "vars.nc");

  list_nc = list_append(list_nc, psi,"psi", "double");
  list_nc = list_append(list_nc, q, "q", "double");
  create_nc(file_tmp);

  /**
     Main Loop
  */

  while(t < tend){

    fprintf(stdout, "i = %d, t = %e dt = %e \n",it, t, dt);

    if (fabs (t - t_out) < TEPS*dt){
      
      t_out += dt_out;
      invert_pv(q, psi, omega);

      // write output
      fprintf(stdout,"Write output, t = %e \n",t);
      write_nc();

    }

    timestep(q);
    it ++;
  }

 
  /**
     Cleanup
  */

#ifdef _STOCHASTIC
  clean_stoch_forcing();
#endif

  clean_fft();
  clean_timestep();
  clean_eigmode();
  
  free(psi);
  free(q);
  free(X);
  free(Y);
  free(K);
  free(L);
  
  if (list_nc) list_free(list_nc);
  if (params) list_free(params);
  
  finalize();

}
