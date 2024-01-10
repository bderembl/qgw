
/**
   Domain related routines
   
*/

// grid indices
#define idx(i,j,k) (k)*Nxp1*Nyp1 + (j)*Nxp1 + (i)

void init_domain() {

  Delta = Lx/NX;
  Ly = NY*Delta;

  // Coordinates in real space
  X = calloc( Nxp1, sizeof( double ) );
  Y = calloc( Nyp1, sizeof( double ) );

  for(int i = 0; i <Nxp1; i++)
    X[i] = i*Delta;
  for(int j = 0; j<Nyp1; j++)
    Y[j] = (j + J0 - 1)*Delta;


  // Coordinates in spectral space
  K = calloc( Nxp1, sizeof( double ) );
  L = calloc( Nyp1, sizeof( double ) );

  for(int i = 1; i <Nx; i++)
    K[i] = pi*(i)/Lx;

  for(int j = 1; j<Ny; j++)
    L[j] = pi*(j + J0)/Ly;

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

  // Overload function to perform communication with MPI neighbours and boundary adjustment at the same time

  double psi_bc = 0.;

  for(int k = 0; k<nl; k++){
    #ifdef _MPI
      if (rank == 0 && rank == n_ranks-1){ // One core only, adjust boundaries like in the single core case
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
      } else if (rank  == 0){ // First rank, don't bc_set north, but send/receive instead
        // South
        for(int i = 0; i <Nxp1; i++){
          int j = 0;
          psi[idx(i,j,k)] = psi_bc;
          q[idx(i,j,k)] = 2*bc_fac/sq(Delta)*(psi[idx(i,j+1,k)] - psi_bc);;
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

        // North
        //send/receive psi
        MPI_Send(&psi[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        
        MPI_Status  status;
        MPI_Recv(&psi[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
        
        // send/receive q
        MPI_Send(&q[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        
        MPI_Status  status2;
        MPI_Recv(&q[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status2);
        
      } else if (rank == n_ranks-1){ // Last rank, interacts with the south

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

        // South
        if (rank%2 == 1){ //odd ranks first receive, then send

        //send/receive psi
        MPI_Status  status;
        MPI_Recv(&psi[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&psi[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        
        // send/receive q
        MPI_Status  status2;
        MPI_Recv(&q[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status2);
        MPI_Send(&q[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        } else { // even ranks first send, then receive

        //send/receive psi
        MPI_Status  status;
        
        MPI_Send(&psi[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        MPI_Recv(&psi[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
        
        // send/receive q
        MPI_Status  status2;

        MPI_Send(&q[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        MPI_Recv(&q[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status2);
        
        }
      } else { // inner ranks

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

        if (rank%2 == 0){ // even inner ranks (send first)
          
          //send/receive psi
          MPI_Status  status;
          
          // send
          MPI_Send(&psi[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); // South
          MPI_Send(&psi[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD); // North
          
          // receive
          MPI_Recv(&psi[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status); //South
          MPI_Recv(&psi[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status); //North
          
          // send/receive q
          MPI_Status  status2;
          
          // send
          MPI_Send(&q[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); // South 
          MPI_Send(&q[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD); // North
          
          // receive
          MPI_Recv(&q[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status2); //South
          MPI_Recv(&q[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status2); // North
        
        } else { // odd inner ranks (receive first)
          
          //send/receive psi
          MPI_Status  status;
          
          // receive
          MPI_Recv(&psi[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status); //South
          MPI_Recv(&psi[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status); //North
          
          // send
          MPI_Send(&psi[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); // South
          MPI_Send(&psi[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD); // North
          
          // send/receive q
          MPI_Status  status2;
          
          // receive
          MPI_Recv(&q[idx(0,0,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status2); //South
          MPI_Recv(&q[idx(0,Ny,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status2); // North
          
          // send
          MPI_Send(&q[idx(0,1,k)], Nxp1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); // South 
          MPI_Send(&q[idx(0,Nym1,k)], Nxp1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD); // North
        
        }
        
        
      }

    #else
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
    #endif
  }
}
