
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
#define forcing_q(t) (-tau0/dh[0]*forc_mode*pi/Ly*sin(forc_mode*pi*(Y[j])/Ly))

// topographic shelf
#define shelf(x,d) (1 - exp(-sq(x)/(2*sq(d))))

// global size
int NX, NY;
int NXp1, NYp1;
int NXp2, NYp2;

// local size
int Nx, Ny;
int Nxp1, Nyp1;
int Nxp2, Nyp2;

// Physical grid
double *X;
double *Y;

// Memory size for Fourier grids
ptrdiff_t alloc_local; 

// Fourier Coefficients
double *K;
double *L;
int N_c;
int N_fft;

// grid indices
#define idx(i,j,k) (k)*Nxp2*Nyp2 + (j)*Nxp2 + (i)

void init_domain() {
  
  Delta = Lx/NX;
  Ly = NY*Delta;
  
  // In bounded domains, points are discretized at cell corners: there is one
  // less active point
  if (bc_fac >= 0) {
    NX = NX - 1;
    NY = NY - 1;
  } 

  /**
     Prepare local indices
  */

  Nx = NX;
  Ny = NY;

  Nxp1 = Nx + 1;
  Nyp1 = Ny + 1;

  Nxp2 = Nx + 2;
  Nyp2 = Ny + 2;

  NXp1 = NX + 1;
  NYp1 = NY + 1;

  NXp2 = NX + 2;
  NYp2 = NY + 2;

  N_c = Nx/2 +1;

  J0 = 0;
  I0 = 0;
  
  K0 = 0;
  Nk = NXp2/2;

  alloc_local = Nx*Ny;

#ifdef _MPI

  int blocksize = ceil((double)NY/n_ranks);
  rank_crit = ceil((double)NY/blocksize)-1;
  

  if (bc_fac == -1) { // for periodic runs the distribution of the transposed array is different
    /* get local data size and allocate */
    /* Domain decomposition is given by FFTW */
    ptrdiff_t local_n0, local_0_start;
    alloc_local = fftw_mpi_local_size_2d(NY, NXp2/2, MPI_COMM_WORLD, // call for the size of the complex array
                                        &local_n0, &local_0_start);
    J0   = local_0_start;
    Ny   = local_n0;
    Nyp1 = local_n0 + 1;
    Nyp2 = local_n0 + 2;
    
    ptrdiff_t local_n1, local_1_start;
    alloc_local = fftw_mpi_local_size_2d(NXp2/2, NY, MPI_COMM_WORLD, // call for the starting indices of the transposed complex array
                                        &local_n1, &local_1_start);
    K0 = local_1_start;
    Nk = local_n1;
    
  } else {
    
    /* get local data size and allocate */
    /* Domain decomposition is given by FFTW */
    ptrdiff_t local_n0, local_0_start;
    alloc_local = fftw_mpi_local_size_2d(NY, NX, MPI_COMM_WORLD, // call for the size of the complex array
                                        &local_n0, &local_0_start);
    J0   = local_0_start;
    Ny   = local_n0;
    Nyp1 = local_n0 + 1;
    Nyp2 = local_n0 + 2;

    Nk = Ny;
  }
#endif

  if (bc_fac == -1) {
    // If periodic, then the zero of X and Y is at the second grid point
    J0 = J0 - 1;
    I0 = I0 - 1;
  }

  // Coordinates in real space
  X = calloc( Nxp2, sizeof( double ) );
  Y = calloc( Nyp2, sizeof( double ) );

  for(int i = 0; i <Nxp2; i++)
    X[i] = (i + I0)*Delta;
  for(int j = 0; j<Nyp2; j++)
    Y[j] = (j + J0)*Delta;

  // Coordinates in spectral space
  L = calloc( Nxp2, sizeof( double ) );
  K = calloc( Nyp2, sizeof( double ) );

  if (bc_fac == -1) {
    for(int j = 0; j < Nx; j++)
      L[j] = fmodf(2*pi*(j)/Lx + 2*pi*Nx/(2*Lx), 2*pi*Nx/Lx) - 2*pi*Nx/(2*Lx);
    for(int i = 0; i < Nk; i++)
      K[i] = 2*pi*(i + K0)/Lx;

  } else {
    for(int j = 1; j<Nyp1; j++)
      K[j] = pi*(j + J0)/Lx;
    for(int i = 1; i <Nxp1; i++)
      L[i] = pi*(i + I0)/Lx;

  }
}

void  init_vars(){

  /* 3d variables */
  psi   = calloc( Nxp2*Nyp2*nl, sizeof( double ) );
  q     = calloc( Nxp2*Nyp2*nl, sizeof( double ) );
  omega = calloc( Nxp2*Nyp2*nl, sizeof( double ) );

  /* 2d variables */
  topo = calloc( Nxp2*Nyp2, sizeof( double ) );

  for(int k = 0; k<nl; k++){
    for(int j = 0; j<Nyp2; j++){
      for(int i = 0;i <Nxp2; i++){
        q[idx(i,j,k)] = 0.;
        psi[idx(i,j,k)] = 0.;
        omega[idx(i,j,k)] = 0.;
      }
    }
  }

  for(int j = 0; j<Nyp2; j++){
    for(int i = 0;i <Nxp2; i++){
      topo[idx(i,j,0)] = h_topo*(1-shelf(X[i],w_topo)*
                                 shelf(Lx-X[i],w_topo)*
                                 shelf(Y[j],w_topo)*
                                 shelf(Ly-Y[j],w_topo));
    }
  }
  
}

/**
   Adjust boundary condition for any field
   field_type = 0: psi
   field_type = 1: q or omega
   
   we need psi to adjust the boundary condition for q and omega
 */
void adjust_bc(double *field, int field_type, double *psi) {

  // first periodic boundary conditions, then MPI communication and then physical domain BC

  if (bc_fac == -1) { // periodic boundary conditions

    for(int k = 0; k<nl; k++){
      // South
      for(int i = 0; i <Nxp2; i++){
        int j1 = 0;
        int j2 = Ny;
        field[idx(i,j1,k)] = field[idx(i,j2,k)];
      }

      // North
      for(int i = 0; i <Nxp2; i++){
        int j1 = Nyp1;
        int j2 = 1;
        field[idx(i,j1,k)] = field[idx(i,j2,k)];
      }
    
      // West
      for(int j = 0; j <Nyp2; j++){
        int i1 = 0;
        int i2 = Nx;
        field[idx(i1,j,k)] = field[idx(i2,j,k)];
      }
    
      // East
      for(int j = 0; j <Nyp2; j++){
        int i1 = Nxp1;
        int i2 = 1;
        field[idx(i1,j,k)] = field[idx(i2,j,k)];
      }
    }
  }

#ifdef _MPI

  int id_so = 1;    // index to send at southern boundary
  int id_no = Ny;   // index to send at northern boundary
  int rank_m1 = rank - 1;
  int rank_p1 = rank + 1;
  
  // adjust indices and ranks for first and last processor 
  if (rank == 0) { // south
    id_so = 1;
    rank_m1 = rank_crit;
  }
  if (rank == rank_crit){ // north
    id_no = Ny;
    rank_p1 = 0;
  }

  if (rank <= rank_crit){
    for(int k = 0; k<nl; k++){
      if (rank%2 == 0){ // even inner ranks (send first)
        
        //send/receive field
        MPI_Status  status;
            
        // send
        MPI_Send(&field[idx(0,id_so,k)], Nxp2, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&field[idx(0,id_no,k)], Nxp2, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
            
        // receive
        MPI_Recv(&field[idx(0,Nyp1,k)], Nxp2, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status); //North
        MPI_Recv(&field[idx(0,0,k)],  Nxp2, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status); //South
 
      } else { // odd inner ranks (receive first)
            
        //send/receive field
        MPI_Status  status;
            
        // receive
        MPI_Recv(&field[idx(0,Nyp1,k)], Nxp2, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD, &status); //North
        MPI_Recv(&field[idx(0,0,k)],  Nxp2, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD, &status); //South
            
        // send
        MPI_Send(&field[idx(0,id_so,k)], Nxp2, MPI_DOUBLE, rank_m1, 0, MPI_COMM_WORLD); // South
        MPI_Send(&field[idx(0,id_no,k)], Nxp2, MPI_DOUBLE, rank_p1, 0, MPI_COMM_WORLD); // North
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
        for(int i = 0; i <Nxp2; i++){
          int j = 0;
          field[idx(i,j,k)] = (field_type ? 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc): psi_bc);
        }
  #ifdef _MPI
      }
  #endif

        // North
  #ifdef _MPI
      if (rank == rank_crit){
  #endif
        for(int i = 0; i <Nxp2; i++){
          int j = Nyp1;
          field[idx(i,j,k)] = (field_type ? 2*bc_fac/sq(Delta)*(psi[idx(i,j-1,k)] - psi_bc) : psi_bc);
        }
  #ifdef _MPI
      }
  #endif
      
        // West
        for(int j = 0; j <Nyp2; j++){
          int i = 0;
          field[idx(i,j,k)] = (field_type ? 2*bc_fac/sq(Delta)*(psi[idx(i+1,j,k)] - psi_bc) : psi_bc);
        }
      
        // East
        for(int j = 0; j <Nyp2; j++){
          int i = Nxp1;
          field[idx(i,j,k)] = (field_type ? 2*bc_fac/sq(Delta)*(psi[idx(i-1,j,k)] - psi_bc) : psi_bc);
        }
    }
  }
}
