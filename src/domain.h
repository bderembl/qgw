
/**
   Domain related routines
   
*/


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

// grid indices
#define idx(i,j,k) (k)*Nxp1*Nyp1 + (j)*Nxp1 + (i)

void init_domain() {

  Delta = Lx/NX;
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

  J0 = 0;

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


  // Coordinates in real space
  X = calloc( Nxp1, sizeof( double ) );
  Y = calloc( Nyp1, sizeof( double ) );

  for(int i = 0; i <Nxp1; i++)
    X[i] = i*Delta;
  for(int j = 0; j<Nyp1; j++)
    Y[j] = (j + J0)*Delta;


  // Coordinates in spectral space
  L = calloc( Nxp1, sizeof( double ) );
  K = calloc( Nyp1, sizeof( double ) );

  for(int i = 1; i <Nx; i++)
    L[i] = pi*(i)/Lx;

  for(int j = 1; j<Ny; j++)
    K[j] = pi*(j + J0)/Ly;

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

// first MPI inner communation and then physical domain BC

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
        psi[idx(i,j,k)] = psi_bc;
        q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);;
      }

      // North
#ifdef _MPI
    if (rank == rank_crit)
#endif
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
