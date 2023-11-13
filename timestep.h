
/**
   Timestep routine.

   Implemented options are 
     - Adams Bashforth

   TODO:
     - Variable time step
 */

double *f1;
double *f2;
double *f3;

double dt_max = 0;
#define min(p,q) p > q ? q : p

/*Calculate dt_max*/


void check_timestep(){
	/** Compares the timestep to viscous and beta numerical instabilities, 
			and decreases it if necessary.
	 */ 
	 
		if (beta != 0 && nu != 0) {
		dt_max = min(1/(2.*beta*Lx),0.5*sq(Lx/Nx)/nu/4.);
	} else if (beta == 0 && nu != 0) {
		dt_max = 1/(2.*beta*Lx);
	} else if (beta != 0 && nu == 0) {
		dt_max = 0.5*sq(Lx/Nx)/nu/4.;
	}
	
	if (dt_max != 0 && dt > dt_max) {
	dt = dt_max; 
	}
	
}

void init_timestep(){

  f1 = calloc( Nxp1*Nyp1, sizeof( double ) );
  f2 = calloc( Nxp1*Nyp1, sizeof( double ) );
  f3 = calloc( Nxp1*Nyp1, sizeof( double ) );
	
	check_timestep();
}


void timestep(double * q){

  rhs(q, f1);
  for (int j = 0; j < Ny; j++){
    for (int i = 0; i < Nx; i++){
    q[idx(i,j)] = q[idx(i,j)]  + (dt/12.)*(23.*f1[idx(i,j)] - 16.*f2[idx(i,j)] + 5.*f3[idx(i,j)] );
    }
  }

  double *swap = f3;
  f3 = f2;
  f2 = f1;
  f1 = swap;

  t = t + dt;

}

void clean_timestep(){

  free(f1);
  free(f2);
  free(f3);

}
