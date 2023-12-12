
/**
   Domain related routines
   
   TODO
     -MPI routines
*/

// grid indices
#define idx(i,j,k) (k)*Nxp1*Nyp1 + (j)*Nxp1 + (i)

void init_domain() {

  Nxm1 = Nx - 1;
  Nym1 = Ny - 1;

  Nxp1 = Nx + 1;
  Nyp1 = Ny + 1;

  Delta = Lx/Nx;
  Ly = Ny*Delta;


  X = calloc( Nxp1, sizeof( double ) );
  Y = calloc( Nxp1, sizeof( double ) );

  for(int i = 0; i <Nxp1; i++)
    X[i] = i*Delta;
  for(int j = 0; j<Nyp1; j++)
    Y[j] = j*Delta;

  for (int l = 0; l < nl-1; l++)
    dhc[l] = 0.5*(dh[l] + dh[l+1]);
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
