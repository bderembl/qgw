/**
   Stochastic forcing functions

   TODO: forcing only active in upper layer. needs to be generalized
	 
*/

// forcing parameters
double sigma_f = 1.;
double k_f = 0.1;

// global constants
ptrdiff_t alloc_forc; 
int NK_f;
int Nk_f;
int N_F;
int N_f;
int J0_f;
int K0_f;
double dk;

double *forc; 
double *forc_f; 

//4d forcing
double dt_forc = 0;
double dt_forc_period = 0;
int n_rec_forc = 0;
int irec_forc1 = -1;
int irec_forc2 = 0;

double *q_forc_3d;
double *q_forc_3d_t1;
double *q_forc_3d_t2;

// forcing lists
List *list_forc1;
List *list_forc2;

#define idx_fft_f(i,j) (j)*2*NK_f + (i)

#ifdef _MPI
  #define idx_fft2_f(i,j) (i)*(2*N_F) + 2*(j)
#else
  #define idx_fft2_f(i,j) (j)*2*NK_f + 2*(i)
#endif

fftw_plan transfo_inverse_forc;

#define noise() (((double) rand()/RAND_MAX))
#define normal_noise() (sqrt(-2.*log(( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. ) ))*cos(2*pi*rand()/(double)RAND_MAX))

// functions to calculate forcing

void  init_stoch_forc(){
	
  // The result is of size Nxp1 x Nyp1 (or Nx x Ny if periodic BCs), but the fields (and the forcing vector) are Nxp2 x Nyp2.

  forc = calloc(Nxp2*Nyp2, sizeof( double ) );
  
  if (bc_fac == -1) {
    N_F = NX;
  } else {
    N_F = NXp1;
  }

  NK_f = N_F/2 + 1;

  Nk_f = NK_f;
  N_f = N_F;

  alloc_forc = Nk_f*N_f;

#ifdef _MPI

  ptrdiff_t N[] = {N_F, NK_f};
  ptrdiff_t local_n0;
  ptrdiff_t local_0_start;
  ptrdiff_t local_n1;
  ptrdiff_t local_1_start;
  
  alloc_forc = fftw_mpi_local_size_2d(N_F, NK_f, MPI_COMM_WORLD, // call for the size of the complex array
                                        &local_n0, &local_0_start);
  J0_f = local_0_start;
  N_f = local_n0;

  alloc_forc = fftw_mpi_local_size_2d(NK_f, N_F, MPI_COMM_WORLD, // call for the starting indices of the transposed complex array
                                        &local_n1, &local_1_start);
  K0_f = local_1_start;
  Nk_f = local_n1;
  
#endif

  forc_f = fftw_alloc_real(2*alloc_forc);

#ifdef _MPI

  ptrdiff_t NN[] = {N_F, N_F};
  transfo_inverse_forc = fftw_mpi_plan_dft_c2r_2d(N_F, N_F, (fftw_complex*) forc_f, forc_f, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
#else

  transfo_inverse_forc = fftw_plan_dft_c2r_2d(N_f, N_f, (fftw_complex*) forc_f, forc_f, FFTW_EXHAUSTIVE);

#endif
  
  dk = 1./Lx;

  // upper layer only
  int k = 0;

  for(int j = 0; j<Nyp2; j++){
    for(int i = 0; i <Nxp2; i++){
      forc[idx(i,j,k)] = 0.;
    }
  }
  
  for(int i = 0; i < (2*alloc_forc); i++){
    forc_f[i] = 0.;
  }
}

void calc_forc() {
  
  // initiate forcing in spectral space
  for(int j = 0; j < N_F; j++){
    for(int i = 0; i < Nk_f; i++){
      
      double l = (((j + N_F/2) % N_F) - N_F/2)*dk;
      double k = 1./Ly*(i + K0_f);

      double K2 = sq(k) + sq(l);
      double norm = sqrt(8*pi*K2)*dk;
      double envelope = sqrt(norm*exp(-sq(sqrt(K2) - k_f)/(2*sq(dk))));

      double magnitude = sigma_f*normal_noise();
      double phase = noise()*2*pi;

      forc_f[idx_fft2_f(i,j)] = envelope*magnitude*cos(phase);
      forc_f[idx_fft2_f(i,j) + 1] = envelope*magnitude*sin(phase);

    }
  }
  

  // execute FFT
  fftw_execute(transfo_inverse_forc);
  
  // upper layer only
  int k = 0;

  for(int j = 0; j<N_f; j++){
    for(int i = 0; i<N_F; i++){
      forc[idx(i+1,j+1,k)] = forc_f[idx_fft_f(i,j)];
    }
  }
  
}

void clean_stoch_forcing(){
  free(forc);
  fftw_free(forc_f);
  fftw_destroy_plan(transfo_inverse_forc);
}


/**
   4d forcing
 */

void  init_4d_forcing(){

  if (dt_forc_period){
    q_forc_3d    = calloc(Nxp2*Nyp2*nl, sizeof( double ) );
    q_forc_3d_t1 = calloc(Nxp2*Nyp2*nl, sizeof( double ) );
    q_forc_3d_t2 = calloc(Nxp2*Nyp2*nl, sizeof( double ) );

    // same name (need two lists to load different records in 2 variables)
    list_forc1 = list_append(list_forc1, q_forc_3d_t1,"q_forcing_3d", "double");
    list_forc2 = list_append(list_forc2, q_forc_3d_t2,"q_forcing_3d", "double");

    if (dt_forc == 0){
      fprintf(stdout,"if dt_forc_period is non zero, dt_forc must be non zero\n");
      exit(1);
    }

    n_rec_forc = (int)floor(dt_forc_period/dt_forc);

    if (n_rec_forc*dt_forc != dt_forc_period){
      fprintf(stdout,"dt_forc_period must be a multiple of dt_forc \n");
      exit(1);
    }

    fprintf(stdout,"init 4D forcing. \n");
  }

}


void calc_4d_forcing(){

  if (dt_forc_period){

    double w_forc1, w_forc2;
    int irec1;

    irec1 = ((int)floor(t/dt_forc))%n_rec_forc;
    
    // reload forcing
    if (irec1 != irec_forc1) {
      irec_forc1 = irec1;
      irec_forc2 = (irec1 + 1)%n_rec_forc;

      read_nc(list_forc1, "input_vars.nc", irec_forc1);
      read_nc(list_forc2, "input_vars.nc", irec_forc2);
    }

    // recompute weights and interpolate forcing
    w_forc2 = t/dt_forc - (int)floor(t/dt_forc);
    w_forc1 = 1 - w_forc2;

    for(int k = 0; k<nl; k++){
      for (int j = 1; j < Nyp1; j++){
        for (int i = 1; i < Nxp1; i++){
          q_forc_3d[idx(i,j,k)] = w_forc1*q_forc_3d_t1[idx(i,j,k)]
            + w_forc2*q_forc_3d_t2[idx(i,j,k)];
        }
      }
    }

  } // end if
}


void clean_4d_forcing(){

  if (dt_forc_period){
    list_free(list_forc1);
    list_free(list_forc2);

    free(q_forc_3d);
    free(q_forc_3d_t1);
    free(q_forc_3d_t2);
  }
}
