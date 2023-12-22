/**
   QG model 

   Compile with 
     gcc -O3 -Wall qg.c -o qg.e -lm -lfftw3 -llapacke -lnetcdf 
   If run with MPI compile with 
    mpicc -O3 qg.c -o qg.e -lfftw3_mpi -lfftw3 -lm -llapacke -lnetcdf -D_MPI (-D_PRINT)

   Compilation flags
     -D_STOCHASTIC : add a stochastic forcing
     -D_MPI : needed for mpi compilation
     -D_PRINT : stdout printing is disabled in MPI, but can be reenabled by including this flag (for debugging)
    
   Run with
     ./qg.e
   or with
     mpirun -n NPROC qg.e

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
  #include <unistd.h> //just for sleep function
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
double *X;
double *Y;

// Fourier Coefficients
double *K;
double *L;

// output variables for MPI output
double *psi_out; 
double *q_out;
int NYp1;

// FFTW in/outputs and plans
double *in1;
double *in2; 
double *out1; 
double *out2; 

fftw_plan transfo_direct, transfo_inverse;

// space and time constants
int Nx, Ny;
int Nxm1, Nym1;
int Nxp1, Nyp1;

//local or global MPI indices
int Nyt; 
int Nytm1;
int Nytp1;
int Ny_start; 
int Ny_startm1;

int rank;
int n_ranks;

// variable for printing out intermediate initialisation info
int print = 1;

// MPI FFTW
ptrdiff_t NY, NX;
ptrdiff_t alloc_local, local_n0, local_0_start;

// MPI communication variables
int *size_gather;
int *start_gather;
int *rows_gather;
int size_gather_local;
int Ny_send_start;
int Ny_send_rows;

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
double tau0 = 0.;
double forc_mode = 1.0;
double f0 = 1.e-5;
double bc_fac = 0.;
double N2[nl_max] = {1.};

#define forcing_q(t) (-tau0/dh[0]*forc_mode*pi/Ly*sin(forc_mode*pi*Y[j]/Ly))

// Local header files

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

  // sleep loop to attach debugger
  // int ii=0;
  // while (0 == ii) sleep(5);

  #ifdef _MPI
    // Initiate MPI
    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    // find out your own rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // find out total number of ranks
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    // disable stdout printing
    print = 0;

  #endif

  #ifdef _PRINT
    print = 1;
  #endif
  /**
     Namelist and parameters
   */

  // add here variables to be read in the input file
  params = list_append(params, &Nx, "Nx", "int");
  params = list_append(params, &Ny, "Ny", "int");
  params = list_append(params, &nl, "nl", "int");
  params = list_append(params, &Lx, "Lx", "double");
  params = list_append(params, &dh, "dh", "array");
  params = list_append(params, &tau0, "tau0", "double");
  params = list_append(params, &forc_mode, "forc_mode", "double");
  params = list_append(params, &f0, "f0", "double");
  params = list_append(params, &beta, "beta", "double");
  params = list_append(params, &nu, "nu", "double");
  params = list_append(params, &N2, "N2", "array");
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

  // Only rank 0 should create outdir when run in MPI
  #ifdef _MPI
    if (rank == 0){
      create_outdir();
      backup_config(file_param);
    }
  #else   
    create_outdir();
    backup_config(file_param);
  #endif
  
  /**
     Initialization
   */

  init_elliptic();
  init_domain();
  init_vars();
  init_timestep();

  #ifdef _STOCHASTIC
    init_stoch_forc();
    printf("Stochastic forcing. \n");
  #endif

  // First inversion
  invert_pv(q,psi);

  // Initialize output
  char file_tmp[90];
  sprintf (file_tmp,"%s%s", dir_out, "vars.nc");

  #ifdef _MPI
    gather_info();
    if (rank == 0){
      psi_out = calloc( Nxp1*Nytp1*nl, sizeof( double ) );
      q_out = calloc( Nxp1*Nytp1*nl, sizeof( double ) );
      list_nc = list_append(list_nc, psi_out,"psi", "double");
      list_nc = list_append(list_nc, q_out, "q", "double");
      create_nc(file_tmp);
    }
  #else 
    list_nc = list_append(list_nc, psi,"psi", "double");
    list_nc = list_append(list_nc, q, "q", "double");
    create_nc(file_tmp);
  #endif

  /**
     Main Loop
  */

  while(t < tend){

    #ifdef _MPI
      if (rank == 0){
        // Only first rank prints
        fprintf(stdout, "i = %d, t = %e dt = %e \n",it, t, dt);
      }
    #else 
      fprintf(stdout, "i = %d, t = %e dt = %e \n",it, t, dt);
    #endif

    if (fabs (t - t_out) < TEPS*dt){
      
      t_out += dt_out;
      invert_pv(q,psi);

      // write output
      #ifdef _MPI
        gather_output();
        if (rank == 0){
          printf("Write output, t = %e \n",t);
          write_nc();
        }
      #else 
        printf("Write output, t = %e \n",t);
        write_nc();
      #endif

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
  
  if (list_nc) list_free(list_nc);
  if (params) list_free(params);
  
  #ifdef _MPI    
    free(psi_out);
    free(q_out);
    MPI_Finalize();
  #endif
}
