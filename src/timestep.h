
/**
   Timestep routine.

   Implemented options are 
     - Adams Bashforth

 */

double *f0_ab;
double *f1_ab;
double *f2_ab;

double dt0, dt1;

static  double TEPS = 1e-9;

void check_timestep(){
  /** Compares the timestep to viscous and beta numerical instabilities, and
      decreases it if necessary.
  */
  
  double nu_eff = max(nu,nu_kin);

  if (beta != 0 && nu_eff != 0) {
    DT_max = min(1/(2.*beta*Lx),sq(Lx/Nx)/nu_eff/20.);
  } else if (beta != 0 && nu_eff == 0) {
    DT_max = 1/(2.*beta*Lx);
  } else if (beta == 0 && nu_eff != 0) {
    DT_max = sq(Lx/Nx)/nu_eff/20.;
  }

  /**
     Adjusts dt with forcing
     U ~ forcing*dt*L
     dt = sqrt(cfl*Delta/forcing*L)
     or stochastic forcing
     U ~ sqrt(sigma*dt)
     dt = (cfl^2*Delta^2/sigma)^(1/3)
   */

  if (tau0 != 0 && sigma_f != 0) {
    dt = 1e-3*min(sqrt(cfl*Delta/(tau0/dh[0]*forc_mode*pi/Ly*Ly)) , pow(sq(Delta*cfl)/sigma_f, 1./3.));
  } else if (tau0 != 0 && sigma_f == 0) {
    dt = 1e-3*sqrt(cfl*Delta/(tau0/dh[0]*forc_mode*pi/Ly*Ly));
  } else if (tau0 == 0 && sigma_f != 0) {
    dt = 1e-3*pow(sq(Delta*cfl)/sigma_f, 1./3.);
  }

  if (DT_max != 0 && dt > DT_max) {
    dt = DT_max; 
  }
  if (print) {
  fprintf(stdout,"Maximum time step: DT_MAX = %g \n", DT_max);
  fprintf(stdout,"Initial time step: dt = %g \n", dt);
  }
	
}

void init_timestep(){

  f0_ab = calloc( Nxp2*Nyp2*nl, sizeof( double ) );
  f1_ab = calloc( Nxp2*Nyp2*nl, sizeof( double ) );
  f2_ab = calloc( Nxp2*Nyp2*nl, sizeof( double ) );
	
  check_timestep();

  // initial values
  dt0 = dt;
  dt1 = dt;
}

double adjust_timestep(double *psi) {

  // Adjust dt according to CFL
  double dt_max = HUGE;

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Nyp1; j++){
      for(int i = 1; i<Nxp1; i++){
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

  #ifdef _MPI
    MPI_Allreduce(&dt_max, &dt_max, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); // Compare with dt_max from other ranks
  #endif

  if (dt < dt_max) dt = 2.*dt;
  if (dt > DT_max) dt = DT_max;
  if (dt > dt_max) dt = dt_max;

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
  double *swap = f2_ab;
  f2_ab = f1_ab;
  f1_ab = f0_ab;
  f0_ab = swap;

  // compute new rhs
  rhs(q, f0_ab);

  // adjust time step
  dt1 = dt0;
  dt0 = dt;
  dt = adjust_timestep(psi);

  
  // AB3 with variable time step coefficients
  double c0 = dt*(2*sq(dt) + 3*dt*(2*dt0 + dt1) + 6*dt0*(dt0 + dt1))/(6*dt0*(dt0 + dt1));
  double c1 = sq(dt)*(-dt/3 - dt0/2 - dt1/2)/(dt0*dt1);
  double c2 = sq(dt)*(dt/3 + dt0/2)/(dt1*(dt0 + dt1));

  for(int k = 0; k<nl; k++){
    for (int j = 1; j < Nyp1; j++){
      for (int i = 1; i < Nxp1; i++){
        q[idx(i,j,k)] = q[idx(i,j,k)]  + c0*f0_ab[idx(i,j,k)] + c1*f1_ab[idx(i,j,k)] + c2*f2_ab[idx(i,j,k)];
      }
    }
  }
	
#ifdef _STOCHASTIC
    calc_forc();
    for (int j = 1; j < Nyp1; j++){
      for (int i = 1; i < Nxp1; i++){
        q[idx(i,j,0)] += forc[idx(i,j,0)]*sqrt(dt);
      }
  }
#endif

  if (dt_forc_period){

    calc_4d_forcing();

    for(int k = 0; k<nl; k++){
      for (int j = 1; j < Nyp1; j++){
        for (int i = 1; i < Nxp1; i++){
          q[idx(i,j,k)] += q_forc_3d[idx(i,j,k)]*dt;
        }
      }
    }
  }

  t = t + dt;

}

void clean_timestep(){

  free(f0_ab);
  free(f1_ab);
  free(f2_ab);

}
