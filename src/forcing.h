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
	
  // The result is of size Nx x Ny (or Nxm1 x Nym1 if periodic BCs), but the fields (and the forcing vector) are Nxp1 x Nyp1.

  forc = calloc(Nxp1*Nyp1, sizeof( double ) );
  
  if (bc_fac == -1) {
    N_F = NXm1;
  } else {
    N_F = NX;
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
  transfo_inverse_forc = fftw_mpi_plan_dft_c2r_2d(N_F, N_F, forc_f, forc_f, MPI_COMM_WORLD,
                                            FFTW_EXHAUSTIVE|FFTW_MPI_TRANSPOSED_IN);
#else

  transfo_inverse_forc = fftw_plan_dft_c2r_2d(N_f, N_f, forc_f, forc_f, FFTW_EXHAUSTIVE);

#endif
  
  dk = 1./Lx;

  // upper layer only
  int k = 0;

  for(int j = 0; j<Nyp1; j++){
    for(int i = 0; i <Nxp1; i++){
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
      
      double l = fmodf(1./Lx*j + N_F/2*dk, N_F*dk) - N_F/2*dk;
      double k = 1./Ly*(i + K0_f);

      double K2 = sq(k) + sq(l);
      double norm = sqrt(2*pi*K2)*dk;
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
