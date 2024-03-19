/**
   MPI related routines
 */


int rank;
int n_ranks;
int rank_crit; // last rank to obtain a meaningful data chunk

int J0; 
int I0;
int K0;
int Nk;

#define pid() rank

void init_mpi() {
  
#ifdef _MPI
  // Initiate MPI
  MPI_Init(NULL, NULL);
  fftw_mpi_init();
  
  // find out your own rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // find out total number of ranks
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

  // redirect fprint(stdout to dev/null for rank > 0
  if (rank > 0) {
    freopen ("/dev/null", "w", stdout);
  }
    
#else
  rank = 0;
#endif

}

static void finalize (void)
{
#ifdef _MPI
  if (rank > 0) {
    fclose(stdout);
  }
  MPI_Finalize();
#endif
}
