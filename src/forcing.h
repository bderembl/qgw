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
int K0_f = 0;
double dk;

// thin ring forcing variables
int N_P;
int N_p;
int *ind_i;
int *ind_j;
int N_ind;

double *forc; 
fftw_complex *forc_f; 
double *forc_p;

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

#ifdef _MPI
#define idx_fft_f(i,j) (j)*2*NK_f + (i)
  #define idx_fft2_f(i,j) (i)*N_F + (j)
#else
  #define idx_fft_f(i,j) (j)*N_F + (i)
  #define idx_fft2_f(i,j) (j)*NK_f + (i)
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

  forc_f = fftw_alloc_complex(alloc_forc);
  
#ifdef _MPI

  forc_p = fftw_alloc_real(2*alloc_forc);
  transfo_inverse_forc = fftw_mpi_plan_dft_c2r_2d(N_F, N_F, forc_f, forc_p, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
#else

  alloc_forc = N_F*N_F;
  forc_p = fftw_alloc_real(alloc_forc);
  transfo_inverse_forc = fftw_plan_dft_c2r_2d(N_F, N_F, forc_f, forc_p, FFTW_EXHAUSTIVE);

#endif
  
  dk = 1./Lx;
  N_P = 0;
  N_p = 0;
  N_ind = 4*pi*k_f/dk;
  ind_i = calloc(N_ind, sizeof( int ) );
  ind_j = calloc(N_ind, sizeof( int ) );
  
  // get the indices which are to be forced
  for (int j = 0; j < N_F; j ++){
    for (int i = 0; i < NK_f; i ++){

      double l = (((j + N_F/2) % N_F) - N_F/2)*dk;
      double k = (i)*dk;

      double dist_k = fabs(k - sqrt(sq(k_f) - sq(l)));
      double dist_l = min(fabs(l - sqrt(sq(k_f) - sq(k))), fabs(l + sqrt(sq(k_f) - sq(k))));

      if (dist_k < dk/2 || dist_l < dk/2){ // check whether point is to be forced
        N_P += 1;
        if (i >= K0_f && i < K0_f + Nk_f){ // check whether point is in local domain
          ind_i[N_p] = i - K0_f;
          ind_j[N_p] = j;
          N_p += 1;
        }
      }
    }
  }

  // upper layer only
  int k_ind = 0;

  for(int j = 0; j<Nyp2; j++){
    for(int i = 0; i <Nxp2; i++){
      forc[idx(i,j,k_ind)] = 0.;
    }
  }
  
  for(int i = 0; i < N_F; i++){
    for(int j = 0; j < N_f; j++){
      forc_p[idx_fft_f(i,j)] = 0.;
    }
  }
  
  for(int i = 0; i < (alloc_forc); i++){
    forc_f[i][0] = 0.;
    forc_f[i][1] = 0.;
  }
}

void calc_forc() {

  // upper layer only
  int k_ind = 0;
  
  for(int i = 0; i < (alloc_forc); i++){
    forc_f[i][0] = 0.;
    forc_f[i][1] = 0.;
  }
  
  // initiate forcing in spectral space

  for (int ii = 0; ii < N_p; ii++){
    int i = ind_i[ii];
    int j = ind_j[ii];

    double envelope = sqrt((double) 8/((N_P-1)*2))*pi*k_f;
    double magnitude = sigma_f*normal_noise();
    double phase = noise()*2*pi;

    forc_f[idx_fft2_f(i,j)][0] = envelope*magnitude*cos(phase);
    forc_f[idx_fft2_f(i,j)][1] = envelope*magnitude*sin(phase);
  }
  
  // execute FFT
  fftw_execute(transfo_inverse_forc);
  
  // upper layer only
  int k = 0;

  for(int j = 0; j<N_f; j++){
    for(int i = 0; i<N_F; i++){
      forc[idx(i+1,j+1,k)] = forc_p[idx_fft_f(i,j)];
    }
  }
  
}

void clean_stoch_forcing(){
  free(forc);
  fftw_free(forc_f);
  fftw_free(forc_p);
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
