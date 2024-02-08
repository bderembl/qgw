/**
   Implementation of the RHS of the QG equation

   TODO:
     - Bottom drag
*/


/**
   Compute rhs of pv equation
   The operators are defined in domain.h
*/

void  rhs(double *q, double * f1){

  invert_pv(q, psi, omega);

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        f1[idx(i,j,k)] = -jacobian(psi, q)        \
          - beta_effect(psi)                    \
          + nu*laplacian(q)                     \
          + nu_kin*laplacian(omega);
      }
    }
  }

  // surface and bottom terms
  for(int j = 1; j<Ny; j++){
    for(int i = 1;i <Nx; i++){
      f1[idx(i,j,0)] += forcing_q(t);
      f1[idx(i,j,nl-1)] += -jacobian_lev(psi, topo, nl-1,0)*f0/dh[nl-1];
    }
  }

}
