/**
   Stochastic forcing functions
	 
*/

double *forc; 
fftw_complex *forc_f; 

double *out;
fftw_complex *in;

#define f_idx(i,j) (i)*N_c + (j)
#define c_idx(i,j) (i)*Nx + (j)

fftw_plan transfo_inverse_forc;

double k;
double l;
double dk;

double Norm;
double Envelope;

double magnitude;
double phase;
double PI = 3.1415;

#define normal_noise() (sqrt(-2.*log(( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. ) ))*cos(2*pi*rand()/(double)RAND_MAX))

// functions to calculate forcing

void  init_stoch_forc(){
	
	N_c = (Nx)/2 + 1;
  
  forc     = calloc(Nxp1*Nyp1, sizeof( double ) );
  out      = calloc(Nx*Ny, sizeof( double ) );
  forc_f   = fftw_alloc_complex(Nx*N_c);
  in       = fftw_alloc_complex(Nx*N_c);
	
  for(int j = 0; j<Nyp1; j++){
    for(int i = 0; i <Nxp1; i++){
			forc[idx(i,j)] = 0.;
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
	
  printf("Prepare forcing fft.. \n");
	
  transfo_inverse_forc = fftw_plan_dft_c2r_2d(Nx, Ny, in, out, FFTW_EXHAUSTIVE);
	
  fprintf(stdout,"..done\n");
}

void calc_forc() {
	
  for(int j = 1; j < N_c; j++){
    for(int i = 1; i < Nx; i++){

      // The complex array storing the fourier representation of the forcing stores values of l along the j dimension,
      // and k along the i dimension. k can be positive and negative, with the positive values being stored in 0 < i < N/2
      // and the negative values of k in N/2 < i < N.

      dk = 1./Lx;
      k = fmodf(1./Lx*i + Nx/2*dk, Nx*dk) - Nx/2*dk; 
      l = 1./Ly*j;

      Norm = pow(sq(k) + sq(l),1./2.)*dk*pow(2*pi,1./2.);
      Envelope = pow(Norm*exp(-sq(pow(sq(k) + sq(l),1./2.) - k_forc)/(2*sq(dk))),1./2.);

      magnitude = pow(eps, 1./2.)*normal_noise();
      phase = ((float) rand()/RAND_MAX)*2*PI;

      forc_f[f_idx(i,j)][0] = Envelope*magnitude*cos(phase);
      forc_f[f_idx(i,j)][1] = Envelope*magnitude*sin(phase);
    }
  }
	
	
  for(int j = 0; j < N_c; j++){
    for(int i = 0; i < Nx; i++){
      in[f_idx(i,j)][0] = forc_f[f_idx(i,j)][0];
      in[f_idx(i,j)][1] = forc_f[f_idx(i,j)][1];
    }
  }
	
  fftw_execute(transfo_inverse_forc);

  for(int j = 0; j<Ny; j++){
    for(int i = 0; i<Nx; i++){
      forc[idx(i,j)] = out[c_idx(i,j)];
    }
  }
	
  for(int j = 0; j<Ny; j++){
    forc[idx(Nx,j)] = out[c_idx(0,j)];
  }
	
  for(int i = 0; i<Ny; i++){
    forc[idx(i,Ny)] = out[c_idx(i,0)];
  }
	
  forc[idx(Nx,Ny)] = out[c_idx(0,0)];
}

void clean_stoch_forcing(){
  free(forc);
  free(out);
  fftw_free(in);
  fftw_free(forc_f);
  fftw_destroy_plan(transfo_inverse_forc);
}

void clean_det_forcing(){
  free(forc); 
}


