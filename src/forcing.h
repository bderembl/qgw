/**
   Stochastic forcing functions

   TODO: forcing only active in upper layer. needs to be generalized
	 
*/

// forcing parameters
double sigma_f = 1.;
double k_f = 0.1;

// global constants
int N_c;
double dk;

double *forc; 
fftw_complex *forc_f; 

double *out;
fftw_complex *in;

#define f_idx(i,j) (i)*N_c + (j)
#define c_idx(i,j) (i)*Nx + (j)

fftw_plan transfo_inverse_forc;

#define noise() (((double) rand()/RAND_MAX))
#define normal_noise() (sqrt(-2.*log(( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. ) ))*cos(2*pi*rand()/(double)RAND_MAX))

// functions to calculate forcing

void  init_stoch_forc(){
	
  N_c = Nx/2 + 1;
  dk = 1./Lx;
  
  forc     = calloc(Nxp1*Nyp1, sizeof( double ) );
  out      = calloc(Nx*Ny, sizeof( double ) );
  forc_f   = fftw_alloc_complex(Nx*N_c);
  in       = fftw_alloc_complex(Nx*N_c);

  // upper layer only
  int k = 0;

  for(int j = 0; j<Nyp1; j++){
    for(int i = 0; i <Nxp1; i++){
      forc[idx(i,j,k)] = 0.;
    }
  }
  
  for(int j = 0; j<Ny; j++){
    for(int i = 0; i <Nx; i++){
      out[c_idx(i,j)] = 0.;
    }
  }
  
  for(int j = 0; j<N_c; j++){
    for(int i = 0; i < Nx; i++){
      forc_f[f_idx(i,j)][0] = 0.;
      forc_f[f_idx(i,j)][1] = 0.;
      in[f_idx(i,j)][0] = 0.;
      in[f_idx(i,j)][1] = 0.;
    }
  }
	
  transfo_inverse_forc = fftw_plan_dft_c2r_2d(Nx, Ny, in, out, FFTW_EXHAUSTIVE);
}

void calc_forc() {
  
  // initiate forcing in spectral space
  for(int j = 1; j < N_c; j++){
    for(int i = 1; i < Nx; i++){
      
      // The complex array storing the fourier representation of the forcing stores values of l along the j dimension,
      // and k along the i dimension. k can be positive and negative, with the positive values being stored in 0 < i < N/2
      // and the negative values of k in N/2 < i < N.

      double k = fmodf(1./Lx*i + Nx/2*dk, Nx*dk) - Nx/2*dk;
      double l = 1./Ly*j;

      double K2 = sq(k) + sq(l);
      double norm = sqrt(2*pi*K2)*dk;
      double envelope = sqrt(norm*exp(-sq(sqrt(K2) - k_f)/(2*sq(dk))));

      double magnitude = sigma_f*normal_noise();
      double phase = noise()*2*pi;

      forc_f[f_idx(i,j)][0] = envelope*magnitude*cos(phase);
      forc_f[f_idx(i,j)][1] = envelope*magnitude*sin(phase);
    }
  }
	
  // copy forcing into input array for FFTW operation (perform an out-of-place transform)
  for(int j = 0; j < N_c; j++){
    for(int i = 0; i < Nx; i++){
      in[f_idx(i,j)][0] = forc_f[f_idx(i,j)][0];
      in[f_idx(i,j)][1] = forc_f[f_idx(i,j)][1];
    }
  }
  
  // execute FFT
  fftw_execute(transfo_inverse_forc);
  
  // Copy results into forcing vector to be used for time-stepping
  // The result is of size Nx x Ny, but the fields (and the forcing vector) are Nxp1 x Nyp1 (as the periodic values are doubled).
  // Therefore the last values of the forc() array are manually copied to be the same as the first in the last two loops.
  // Forcing is non-zero on the boundaries, but I'm assuming the values there are calculated manually by the 
  // boundary conditions anyway so I'm hoping it won't make a difference.

  // upper layer only
  int k = 0;

  for(int j = 0; j<Ny; j++){
    for(int i = 0; i<Nx; i++){
      forc[idx(i,j,k)] = out[c_idx(i,j)];
    }
  }
  
  for(int j = 0; j<Ny; j++){
    forc[idx(Nx,j,k)] = out[c_idx(0,j)];
  }
	
  for(int i = 0; i<Ny; i++){
    forc[idx(i,Ny,k)] = out[c_idx(i,0)];
  }
	
  forc[idx(Nx,Ny,k)] = out[c_idx(0,0)];
}

void clean_stoch_forcing(){
  free(forc);
  free(out);
  fftw_free(in);
  fftw_free(forc_f);
  fftw_destroy_plan(transfo_inverse_forc);
}
