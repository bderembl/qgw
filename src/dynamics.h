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

void  calc_omega(double *omega, double *psi){

  // calculation at interior points 
  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        omega[idx(i,j,k)] = laplacian(psi);
      }
    }
  }

  // Boundary values (and MPI communications)

#ifdef _MPI

  int id_so = 1;    // index to send at southern boundary
  int id_no = Nym1; // index to send at northern boundary
  int rank_m1 = rank - 1;
  int rank_p1 = rank + 1;
  
  // adjust indices and ranks for first and last processor 
  if (rank == 0) { // south
    id_so = 0;
    rank_m1 = rank_crit;
  }
  if (rank == rank_crit){ // north
    id_no = Ny;
    rank_p1 = 0;
  }

  if (rank <= rank_crit){
    for(int k = 0; k<nl; k++){
      if (rank%2 == 0){ // even inner ranks (send first)
        
        // send/receive omega
        MPI_Status  status3;
            
        // send
        MPI_Send(&omega[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South 
        MPI_Send(&omega[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // receive
        MPI_Recv(&omega[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status3); // North
        MPI_Recv(&omega[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status3); //South
          
      } else { // odd inner ranks (receive first)
            
        // send/receive omega
        MPI_Status  status3;
            
        // receive
        MPI_Recv(&omega[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status3); // North
        MPI_Recv(&omega[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status3); //South
            
        // send
        MPI_Send(&omega[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&omega[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
          
      }
    }
  }
#endif

  double psi_bc = 0.;

  for(int k = 0; k<nl; k++){
      // South
#ifdef _MPI
    if (rank == 0)
#endif
      for(int i = 0; i <Nxp1; i++){
        int j = 0;
        omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);;
      }

      // North
#ifdef _MPI
    if (rank == rank_crit)
#endif
      for(int i = 0; i <Nxp1; i++){
        int j = Ny;
        omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j-1,k)] - psi_bc);
      }
    
      // West
      for(int j = 0; j <Nyp1; j++){
        int i = 0;
        omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i+1,j,k)] - psi_bc);
      }
    
      // East
      for(int j = 0; j <Nyp1; j++){
        int i = Nx;
        omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i-1,j,k)] - psi_bc);
      }
  }
}


void  rhs(double *q, double * f1){

  invert_pv(q, psi);

  if (nu_kin != 0) {
    calc_omega(omega, psi);
  }

  for(int k = 0; k<nl; k++){
    for(int j = 1; j<Ny; j++){
      for(int i = 1;i <Nx; i++){
        f1[idx(i,j,k)] = -jacobian(psi, q)        \
          - beta_effect(psi)                    \
          + (k == 0 ? forcing_q(t) : 0.)        \
          + nu*laplacian(q)                     \
          + nu_kin*laplacian(omega);
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
