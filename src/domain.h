
/**
   Domain related routines and operators
   
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


// Forcing
#define forcing_q(t) (-tau0/dh[0]*forc_mode*pi/Ly*sin(forc_mode*pi*Y[j]/Ly))

// topographic shelf
#define shelf(x,d) (1 - exp(-sq(x)/(2*sq(d))))

// global size
int NX, NY;
int NXp1, NYp1;
int NXm1, NYm1;

// local size
int Nx, Ny;
int Nxm1, Nym1;
int Nxp1, Nyp1;

// Physical grid
double *X;
double *Y;

// Fourier Coefficients
double *K;
double *L;
int N_c;

// grid indices
#define idx(i,j,k) (k)*Nxp1*Nyp1 + (j)*Nxp1 + (i)

void init_domain() {
  
  // In the periodic case we have one more grid point
  if (bc_fac == -1) {
    // If periodic, then the zero of X and Y is at the second grid point
    Delta = Lx/(NX-1);
  } else {
    Delta = Lx/NX;
  }
  
  Ly = NY*Delta;

  /**
     Prepare local indices
  */

  Nx = NX;
  Ny = NY;

  Nxm1 = Nx - 1;
  Nym1 = Ny - 1;

  Nxp1 = Nx + 1;
  Nyp1 = Ny + 1;

  NXm1 = NX - 1;
  NYm1 = NY - 1;

  NXp1 = NX + 1;
  NYp1 = NY + 1;

  N_c = Nxm1/2 +1;

  J0 = 0;
  I0 = 0;

#ifdef _MPI
  int blocksize;
  blocksize = ceil((double)NYm1/n_ranks);
  rank_crit = ceil((double)NYm1/blocksize)-1;
  
  /* get local data size and allocate */
  /* Domain decomposition is given by FFTW */
  ptrdiff_t alloc_local, local_n0, local_0_start;
  alloc_local = fftw_mpi_local_size_2d(NYm1, NXm1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  J0 = local_0_start;
  Ny = local_n0 + 1;
  Nym1 = local_n0;
  Nyp1 = local_n0 + 2;
#endif

  if (bc_fac == -1) {
    // If periodic, then the zero of X and Y is at the second grid point
    J0 = J0-1;
    I0 = I0-1;
  }

  // Coordinates in real space
  X = calloc( Nxp1, sizeof( double ) );
  Y = calloc( Nyp1, sizeof( double ) );

  for(int i = 0; i <Nxp1; i++)
    X[i] = (i + I0)*Delta;
  for(int j = 0; j<Nyp1; j++)
    Y[j] = (j + J0)*Delta;


  // Coordinates in spectral space
  L = calloc( Nxp1, sizeof( double ) );
  K = calloc( Nyp1, sizeof( double ) );

  if (bc_fac == -1) {
    for(int j = 0; j<Ny; j++)
      K[j] = fmodf(2*pi*(j)/Lx + 2*pi*Nxm1/(2*Lx), 2*pi*Nxm1/Lx) - 2*pi*Nxm1/(2*Lx);
    for(int i = 0; i <Nx; i++)
      L[i] = 2*pi*(i)/Lx;

  } else {
    for(int j = 1; j<Ny; j++)
      K[j] = pi*(j + J0)/Lx;
    for(int i = 1; i <Nx; i++)
      L[i] = pi*(i + I0)/Lx;

  }
}

void  init_vars(){


  psi  = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
  q    = calloc( Nxp1*Nyp1*nl, sizeof( double ) );
  topo = calloc( Nxp1*Nyp1, sizeof( double ) );
  omega = calloc( Nxp1*Nyp1*nl, sizeof( double ) );

  for(int k = 0; k<nl; k++){
    for(int j = 0; j<Nyp1; j++){
      for(int i = 0;i <Nxp1; i++){
        q[idx(i,j,k)] = 0.;
        psi[idx(i,j,k)] = 0.;
        omega[idx(i,j,k)] = 0.;
      }
    }
  }

  for(int j = 0; j<Nyp1; j++){
    for(int i = 0;i <Nxp1; i++){
      topo[idx(i,j,0)] = h_topo*(1-shelf(X[i],w_topo)*
                                 shelf(Lx-X[i],w_topo)*
                                 shelf(Y[j],w_topo)*
                                 shelf(Ly-Y[j],w_topo));
    }
  }
  
}


void adjust_bc(double *q, double *psi, double *omega) {

  // first periodic boundary conditions, then MPI communication and then physical domain BC

  if (bc_fac == -1) { // periodic boundary conditions

    for(int k = 0; k<nl; k++){
      // South
      for(int i = 0; i <Nxp1; i++){
        int j1 = 0;
        int j2 = Nym1;
        psi[idx(i,j1,k)] = psi[idx(i,j2,k)];
        q[idx(i,j1,k)] = q[idx(i,j2,k)];
        omega[idx(i,j1,k)] = omega[idx(i,j2,k)];
      }

      // North
      for(int i = 0; i <Nxp1; i++){
        int j1 = Ny;
        int j2 = 1;
        psi[idx(i,j1,k)] = psi[idx(i,j2,k)];
        q[idx(i,j1,k)] = q[idx(i,j2,k)];
        omega[idx(i,j1,k)] = omega[idx(i,j2,k)];
      }
    
      // West
      for(int j = 0; j <Nyp1; j++){
        int i1 = 0;
        int i2 = Nxm1;
        psi[idx(i1,j,k)] = psi[idx(i2,j,k)];
        q[idx(i1,j,k)] = q[idx(i2,j,k)];
        omega[idx(i1,j,k)] = omega[idx(i2,j,k)];
      }
    
      // East
      for(int j = 0; j <Nyp1; j++){
        int i1 = Nx;
        int i2 = 1;
        psi[idx(i1,j,k)] = psi[idx(i2,j,k)];
        q[idx(i1,j,k)] = q[idx(i2,j,k)];
        omega[idx(i1,j,k)] = omega[idx(i2,j,k)];
      }
    }
  }

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
        
        //send/receive psi
        MPI_Status  status;
            
        // send
        MPI_Send(&psi[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&psi[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // receive
        MPI_Recv(&psi[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status); //North
        MPI_Recv(&psi[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status); //South
            
        // send/receive q
        MPI_Status  status2;
            
        // send
        MPI_Send(&q[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South 
        MPI_Send(&q[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // receive
        MPI_Recv(&q[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status2); // North
        MPI_Recv(&q[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status2); //South
         
        // send/receive omega
        MPI_Status  status3;
            
        // send
        MPI_Send(&omega[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South 
        MPI_Send(&omega[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // receive
        MPI_Recv(&omega[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status3); // North
        MPI_Recv(&omega[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status3); //South
 
 
      } else { // odd inner ranks (receive first)
            
        //send/receive psi
        MPI_Status  status;
            
        // receive
        MPI_Recv(&psi[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status); //North
        MPI_Recv(&psi[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status); //South
            
        // send
        MPI_Send(&psi[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&psi[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // send/receive q
        MPI_Status  status2;
            
        // receive
        MPI_Recv(&q[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status2); // North
        MPI_Recv(&q[idx(0,0,k)],  Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status2); //South
            
        // send
        MPI_Send(&q[idx(0,id_so,k)], Nxp1, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&q[idx(0,id_no,k)], Nxp1, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
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

  if (bc_fac != -1) {
    double psi_bc = 0.;

    for(int k = 0; k<nl; k++){
        // South
  #ifdef _MPI
      if (rank == 0){
  #endif
        for(int i = 0; i <Nxp1; i++){
          int j = 0;
          psi[idx(i,j,k)] = psi_bc;
          q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);
          omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);
        }
  #ifdef _MPI
      }
  #endif

        // North
  #ifdef _MPI
      if (rank == rank_crit){
  #endif
        for(int i = 0; i <Nxp1; i++){
          int j = Ny;
          psi[idx(i,j,k)] = psi_bc;
          q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j-1,k)] - psi_bc);
          omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j-1,k)] - psi_bc);
        }
  #ifdef _MPI
      }
  #endif
      
        // West
        for(int j = 0; j <Nyp1; j++){
          int i = 0;
          psi[idx(i,j,k)] = psi_bc;
          q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i+1,j,k)] - psi_bc);
          omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i+1,j,k)] - psi_bc);
        }
      
        // East
        for(int j = 0; j <Nyp1; j++){
          int i = Nx;
          psi[idx(i,j,k)] = psi_bc;
          q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i-1,j,k)] - psi_bc);
          omega[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i-1,j,k)] - psi_bc);
        }
    }
  }
}
