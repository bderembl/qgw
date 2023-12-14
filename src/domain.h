
/**
   Domain related routines
   
   TODO
     -MPI routines
*/

// grid indices
#define idx(i,j,k) (k)*Nxp1*Nyp1 + (j)*Nxp1 + (i)

void init_domain() {

  Delta = Lx/Nx;
  #ifdef _MPI
    Ly = Nyt*Delta;
  #else
    Ly = Ny*Delta;
  #endif

  // Coordinates in real space
  X = calloc( Nxp1, sizeof( double ) );
  Y = calloc( Nyp1, sizeof( double ) );

  for(int i = 0; i <Nxp1; i++)
    X[i] = i*Delta;
  #ifdef _MPI
    for(int j = 0; j<Nyp1; j++)
      Y[j] = (j+Ny_startm1)*Delta;
  #else
    for(int j = 0; j<Nyp1; j++)
      Y[j] = j*Delta;
  #endif


  // Coordinates in sinus space
  K = calloc( Nxp1, sizeof( double ) );
  L = calloc( Nyp1, sizeof( double ) );

  for(int i = 1; i <Nx; i++)
    K[i] = pi*(i)/Lx;

  #ifdef _MPI
    for(int j = 1; j<Ny; j++)
      L[j] = pi*(j + Ny_start)/Ly;
  #else
    for(int j = 1; j<Ny; j++)
      L[j] = pi*(j)/Ly;
  #endif

}


void  init_vars(){


  psi = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
  q   = calloc( Nxp1*Nyp1*nl, sizeof( double ) );

  for(int k = 0; k<nl; k++){
    for(int j = 0; j<Nyp1; j++){
      for(int i = 0;i <Nxp1; i++){
        q[idx(i,j,k)] = 0.;
        psi[idx(i,j,k)] = 0.;
      }
    }
  }
}


void adjust_bc(double *q, double *psi) {
  
  double psi_bc = 0.;

  for(int k = 0; k<nl; k++){

    // South
    for(int i = 0; i <Nxp1; i++){
      int j = 0;
      psi[idx(i,j,k)] = psi_bc;
      q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);;
    }
  
    // North
    for(int i = 0; i <Nxp1; i++){
      int j = Ny;
      psi[idx(i,j,k)] = psi_bc;
      q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j-1,k)] - psi_bc);
    }
  
    // West
    for(int j = 0; j <Nyp1; j++){
      int i = 0;
      psi[idx(i,j,k)] = psi_bc;
      q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i+1,j,k)] - psi_bc);
    }
  
    // East
    for(int j = 0; j <Nyp1; j++){
      int i = Nx;
      psi[idx(i,j,k)] = psi_bc;
      q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i-1,j,k)] - psi_bc);
    }
  }
}
