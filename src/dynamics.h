/**
   Implementation of the RHS of the QG equation

   TODO:
     - Multi layer
     - Bottom drag
*/


/**
   Operators
*/


#define laplacian(p) (p[idx(i+1,j)] + p[idx(i-1,j)] + p[idx(i,j+1)] + p[idx(i,j-1)] - 4*p[idx(i,j)])/(sq(Delta))

#define jacobian(p,q) ((( p[idx(i+1,j)]-p[idx(i-1,j)])*(q[idx(i,j+1)]-q[idx(i,j-1)]) \
                        +(p[idx(i,j-1)]-p[idx(i,j+1)])*(q[idx(i+1,j)]-q[idx(i-1,j)]) \
                        + p[idx(i+1,j)]*( q[idx(i+1,j+1)] - q[idx(i+1,j-1)])         \
                        - p[idx(i-1,j)]*( q[idx(i-1,j+1)] - q[idx(i-1,j-1)])         \
                        - p[idx(i,j+1)]*( q[idx(i+1,j+1)] - q[idx(i-1,j+1)])         \
                        + p[idx(i,j-1)]*( q[idx(i+1,j-1)] - q[idx(i-1,j-1)])         \
                        + q[idx(i,j+1)]*( p[idx(i+1,j+1)] - p[idx(i-1,j+1)])         \
                        - q[idx(i,j-1)]*( p[idx(i+1,j-1)] - p[idx(i-1,j-1)])         \
                        - q[idx(i+1,j)]*( p[idx(i+1,j+1)] - p[idx(i+1,j-1)])         \
                        + q[idx(i-1,j)]*( p[idx(i-1,j+1)] - p[idx(i-1,j-1)]))        \
                       /(12.*Delta*Delta))

#define beta_effect(p) (beta*(p[idx(i+1,j)] - p[idx(i-1,j)])/(2*Delta))

/**
   Compute rhs of pv equation
*/

void  rhs(double *q, double * f1){

  invert_pv(q, psi);

  for(int j = 1; j<Ny; j++){
    for(int i = 1;i <Nx; i++){
      f1[idx(i,j)] = -jacobian(psi, q)          \
        - beta_effect(psi)                      \
        + forcing_q(t)                          \
        + nu*laplacian(q);
    }
  }

}

