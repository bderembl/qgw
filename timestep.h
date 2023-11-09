
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

void init_timestep(){

  f1 = calloc( Nxp1*Nyp1, sizeof( double ) );
  f2 = calloc( Nxp1*Nyp1, sizeof( double ) );
  f3 = calloc( Nxp1*Nyp1, sizeof( double ) );
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
