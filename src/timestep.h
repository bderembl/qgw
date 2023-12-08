
/**
   Timestep routine.

   Implemented options are 
     - Adams Bashforth

   TODO:
     - Variable time step
 */

double *f0;
double *f1;
double *f2;

double dt0, dt1;

static  double TEPS = 1e-9;

/*Calculate dt_max*/


void check_timestep(){
  /** Compares the timestep to viscous and beta numerical instabilities, and
      decreases it if necessary.
  */
  
  if (beta != 0 && nu != 0) {
    DT_max = min(1/(2.*beta*Lx),sq(Lx/Nx)/nu/20.);
  } else if (beta != 0 && nu == 0) {
    DT_max = 1/(2.*beta*Lx);
  } else if (beta == 0 && nu != 0) {
    DT_max = sq(Lx/Nx)/nu/20.;
  }
  
  if (DT_max != 0 && dt > DT_max) {
    dt = DT_max; 
  }
	
}

void init_timestep(){

  f0 = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
  f1 = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
  f2 = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
	
  check_timestep();

  // initial values
  dt0 = dt;
  dt1 = dt;
}

/**
   Compute maximum possible time step dt based on CFL criterion and next output.
 */


double adjust_timestep(double *psi) {

  // Adjust dt according to CFL
  double dt_max = 1e30;

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        double u = -(psi[idx(i+1,j,k)] - psi[idx(i,j,k)])/Delta;
        double v =  (psi[idx(i,j+1,k)] - psi[idx(i,j,k)])/Delta;
        u = max(fabs(u), fabs(v));

        if (u != 0.) {
          double dt_loc = cfl*Delta/u;
          if (dt_loc < dt_max) dt_max = dt_loc;
        }
      }
    }
  }

  //TODO  MPI reduce here

  if ((dt_max > dt) && (dt < DT_max)){
    dt = (dt + 0.1*dt_max)/1.1;
    if (dt > DT_max) dt = DT_max;
  }

  if (dt_max < dt) dt = dt_max;

  // Adjust dt to reach t_out (from Basilisk)
  if (t_out > t + TEPS) {
    unsigned int n_step = (t_out - t)/dt;
    if (n_step == 0)
      dt = t_out - t;
    else{
      double dt1 = (t_out - t)/n_step;
      if (dt1 > dt*(1. + TEPS))
        dt = (t_out - t)/(n_step + 1);
      else if (dt1 < dt)
        dt = dt1;
    }
  }
  return dt;
}



void timestep(double * q){

  // Swap rhs variables for next values
  double *swap = f2;
  f2 = f1;
  f1 = f0;
  f0 = swap;

  // compute new rhs
  rhs(q, f0);

  // adjust time step
  dt1 = dt0;
  dt0 = dt;
  dt = adjust_timestep(psi);

  
  // AB3 with variable time step coefficients
  double c0 = dt*(2*sq(dt) + 3*dt*(2*dt0 + dt1) + 6*dt0*(dt0 + dt1))/(6*dt0*(dt0 + dt1));
  double c1 = sq(dt)*(-dt/3 - dt0/2 - dt1/2)/(dt0*dt1);
  double c2 = sq(dt)*(dt/3 + dt0/2)/(dt1*(dt0 + dt1));

  for(int k = 0; k<nl; k++){
    for (int j = 0; j < Ny; j++){
      for (int i = 0; i < Nx; i++){
        q[idx(i,j,k)] = q[idx(i,j,k)]  + c0*f0[idx(i,j,k)] + c1*f1[idx(i,j,k)] + c2*f2[idx(i,j,k)];
      }
    }
  }

	
  #ifdef _STOCHASTIC
    calc_forc();
    for (int j = 0; j < Ny; j++){
      for (int i = 0; i < Nx; i++){
        q[idx(i,j)] += forc[idx(i,j)]*sqrt(dt);
      }
  }
  #endif


  t = t + dt;

}

void clean_timestep(){

  free(f0);
  free(f1);
  free(f2);

}
