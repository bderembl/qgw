/**
   Implementation of the RHS of the QG equation

*/


/**
   Compute rhs of pv equation
   The operators are defined in domain.h
*/


void calc_diff(double *psi, double *lap_n_diff) {

  // copy psi
  for(int k = 0; k<nl; k++){
    for(int j = 0; j<Nyp2; j++){
      for(int i = 0; i<Nxp2; i++){
        diff[idx(i,j,k)] = psi[idx(i,j,k)];
      }
    }
  }

  // calculate lap_n_diff until the desired hyperviscosity index (0 corresponds to linear drag, 2 to standard viscosity)
  for (int n = 0; n <= n_hyper; n += 2){
    for(int k = 0; k<nl; k++){
      for(int j = 1; j<Nyp1; j++){
        for(int i = 1; i<Nxp1; i++){
          lap_n_diff[idx(i,j,k)] = laplacian(diff);
        }
      }
    }

    adjust_bc(lap_n_diff, 1, psi); // should only be used with free-slip (or periodic) BC, as the BC for lap_n_diff is also dirichlet in the bounded case (:= super-slip)

    for(int k = 0; k<nl; k++){
      for(int j = 0; j<Nyp2; j++){
        for(int i = 0; i<Nxp2; i++){
          diff[idx(i,j,k)] = lap_n_diff[idx(i,j,k)];
        }
      }
    }
  }
}


void  rhs(double *q, double * f1){

  invert_pv(q, psi);
  calc_diff(psi, lap_n_diff);

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Nyp1; j++){
      for(int i = 1; i<Nxp1; i++){
        f1[idx(i,j,k)] = -jacobian(psi, q)      \
          - beta_effect(psi)                    \
          + nu*laplacian(q)                     \
          + pow(-1, n_hyper/2 - 1)*nu_hyper*lap_n_diff[idx(i,j,k)];
      }
    }
  }

  // surface and bottom terms
  for(int j = 1; j<Nyp1; j++){
    for(int i = 1; i<Nxp1; i++){
      f1[idx(i,j,0)] += forcing_q(t);
      int k = nl-1; // for laplacian(psi) macro
      f1[idx(i,j,nl-1)] += - hEkb*f0/(2*dh[nl-1])*laplacian(psi)        \
        - jacobian_lev(psi, topo, nl-1,0)*f0/dh[nl-1];
    }
  }

}
