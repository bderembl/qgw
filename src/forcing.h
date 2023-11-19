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

/* Function to approximate inverse cumulative function of normal distribution, which I use to transform the random numbers from C
 to a normally distributed variable. From https://web.archive.org/web/20150910081113/http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c */

#include <errno.h>

static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define LOW 0.02425
#define HIGH 0.97575

double
ltqnorm(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -1 /* minus "infinity", changed to 1 due to error in code to handle such small numbers */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return 1 /* "infinity", changed to 1 due to error in code to handle such large numbers */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

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

      Norm = pow(sq(k) + sq(l),1./2.)*dk/(pow(2*pi,3./2.));
      Envelope = pow(Norm*exp(-sq(pow(sq(k) + sq(l),1./2.) - k_forc)/(2*sq(dk))),1./2.)/pow(dt, 1./2.);

      magnitude = pow(eps, 1./2.)*ltqnorm(((float) rand()/RAND_MAX));
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


