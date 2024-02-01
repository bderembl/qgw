/**
   Implementation of the RHS of the QG equation

   TODO:
     - Bottom drag
*/


/**
   Operators
*/


#define laplacian(p) (p[idx(i+1,j,k)] + p[idx(i-1,j,k)] + p[idx(i,j+1,k)] + p[idx(i,j-1,k)] - 4*p[idx(i,j,k)])/(sq(Delta))

#define jacobian(p,q) jacobian_lev(p,q,k,k)

//jacobian_lev is jac(p[k1], q[k2])
#define jacobian_lev(p,q,k1,k2) ((( p[idx(i+1,j,k1)]-p[idx(i-1,j,k1)])*(q[idx(i,j+1,k2)]-q[idx(i,j-1,k2)]) \
                        +(p[idx(i,j-1,k1)]-p[idx(i,j+1,k1)])*(q[idx(i+1,j,k2)]-q[idx(i-1,j,k2)]) \
                        + p[idx(i+1,j,k1)]*( q[idx(i+1,j+1,k2)] - q[idx(i+1,j-1,k2)])         \
                        - p[idx(i-1,j,k1)]*( q[idx(i-1,j+1,k2)] - q[idx(i-1,j-1,k2)])         \
                        - p[idx(i,j+1,k1)]*( q[idx(i+1,j+1,k2)] - q[idx(i-1,j+1,k2)])         \
                        + p[idx(i,j-1,k1)]*( q[idx(i+1,j-1,k2)] - q[idx(i-1,j-1,k2)])         \
                        + q[idx(i,j+1,k2)]*( p[idx(i+1,j+1,k1)] - p[idx(i-1,j+1,k1)])         \
                        - q[idx(i,j-1,k2)]*( p[idx(i+1,j-1,k1)] - p[idx(i-1,j-1,k1)])         \
                        - q[idx(i+1,j,k2)]*( p[idx(i+1,j+1,k1)] - p[idx(i+1,j-1,k1)])         \
                        + q[idx(i-1,j,k2)]*( p[idx(i-1,j+1,k1)] - p[idx(i-1,j-1,k1)]))        \
                       /(12.*Delta*Delta))

#define beta_effect(p) (beta*(p[idx(i+1,j,k)] - p[idx(i-1,j,k)])/(2*Delta))

/**
   Compute rhs of pv equation
*/

void  rhs(double *q, double * f1){

  invert_pv(q, psi);

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        f1[idx(i,j,k)] = -jacobian(psi, q)        \
          - beta_effect(psi)                    \
          + (k == 0 ? forcing_q(t) : 0.)             \
          + nu*laplacian(q);
      }
    }
  }

  // surface and bottom terms
  for(int j = 1; j<Ny; j++){
    for(int i = 1;i <Nx; i++){
      f1[idx(i,j,nl-1)] += -jacobian_lev(psi, topo, nl-1,0)*f0/dh[nl-1];
    }
  }

}

