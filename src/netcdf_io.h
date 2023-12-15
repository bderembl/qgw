#include <netcdf.h>
#define nodata HUGE


#define LAYERS 1
#define NDIMS 4

#define Y_NAME "y"
#define X_NAME "x"
#define REC_NAME "time"
#define LVL_NAME "level"

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
int nc_err;
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return;}

// User global variables
char file_nc[80];
List *list_nc;

/* IDs for the netCDF file, dimensions, and variables. */
int ncid;
int t_varid;

// temporary
int nc_varid[1000];
char * nc_varname[1000];
int nc_rec = -1;

void create_nc(char* file_out)
{  
   // Output size variable is changing with MPI
   #ifdef _MPI
      NYp1 = Nytp1;
      NY = Nyt;
   #else
      NYp1 = Nyp1;
      NY = Ny;
   #endif

   sprintf (file_nc,"%s", file_out);
//  if (pid() == 0) { // master

   /* LOCAL IDs for the netCDF file, dimensions, and variables. */
   int x_dimid, y_dimid, lvl_dimid, rec_dimid;
   int y_varid, x_varid, lvl_varid;
   int dimids[NDIMS];
   
   /* Create the file. */
   if ((nc_err = nc_create(file_nc, NC_CLOBBER, &ncid)))
      ERR(nc_err);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(nc_err);
#if LAYERS
   if ((nc_err = nc_def_dim(ncid, LVL_NAME, nl, &lvl_dimid)))
      ERR(nc_err);
#endif
   if ((nc_err = nc_def_dim(ncid, Y_NAME, NYp1, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, Nxp1, &x_dimid)))
      ERR(nc_err);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */
   if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid,
        		    &t_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, Y_NAME, NC_FLOAT, 1, &y_dimid,
        		    &y_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, X_NAME, NC_FLOAT, 1, &x_dimid,
        		    &x_varid)))
      ERR(nc_err);
#if LAYERS
   if ((nc_err = nc_def_var(ncid, LVL_NAME, NC_FLOAT, 1, &lvl_dimid,
        		    &lvl_varid)))
      ERR(nc_err);
#endif
   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = rec_dimid;
#if LAYERS
   dimids[1] = lvl_dimid;
   dimids[2] = y_dimid;
   dimids[3] = x_dimid;
#else
   dimids[1] = y_dimid;
   dimids[2] = x_dimid;
#endif

   /* Define the netCDF variables */
//   char * str1;
   
   if (list_nc){
     for (int i = 0; i < list_nc[0].len; i++){
       
       if ((nc_err = nc_def_var(ncid, list_nc[i].name, NC_FLOAT, NDIMS,
                                dimids, &nc_varid[i])))
         ERR(nc_err);
     }
   }

   /* /\* Assign units attributes to the netCDF variables. *\/ */
   /* if ((nc_err = nc_put_att_text(ncid, pres_varid, UNITS,  */
   /*      			 strlen(PRES_UNITS), PRES_UNITS))) */
   /*    ERR(nc_err); */
   /* if ((nc_err = nc_put_att_text(ncid, temp_varid, UNITS,  */
   /*      			 strlen(TEMP_UNITS), TEMP_UNITS))) */
   /*    ERR(nc_err); */

   /* End define mode. */
   if ((nc_err = nc_enddef(ncid)))
      ERR(nc_err);

   /*  write coordinates*/
   float yc[NY+1], xc[Nx+1];
   for (int i = 0; i < Nx+1; i++){
      xc[i] = i*Delta;
   }
   for (int i = 0; i < NY+1; i++){
      yc[i] = i*Delta;
   }

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      ERR(nc_err);
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      ERR(nc_err);

#if LAYERS
   float zc[nl];
   for (int i = 0; i < nl; i++){
      zc[i] = i;
   }
   if ((nc_err = nc_put_var_float(ncid, lvl_varid, &zc[0])))
      ERR(nc_err);
#endif

   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   fprintf(stdout,"*** SUCCESS creating example file %s!\n", file_nc);
//  } // end master

}

/** Write NetCDF file

*/



void write_nc() {

//  if (pid() == 0) { // master
    /* open file. */
    if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid)))
      ERR(nc_err);
//  }

  // write time
  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;
  /* if (pid() == 0) { // master */
    if ((nc_err = nc_put_vara_float(ncid, t_varid, startt, countt,
                                    &loctime)))
      ERR(nc_err);
  /* } // master */



  float * field = (float *)malloc(Nxp1*NYp1*nl*sizeof(float));

//  float ** field = matrix_new (N_out, N_out, sizeof(float));
  
  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  size_t start[NDIMS], count[NDIMS];
  
  
  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  start[0] = nc_rec; //time
#if LAYERS
  start[1] = 0;     //level
  start[2] = 0;      //y
  start[3] = 0;      //x
#else
  start[1] = 0;      //y
  start[2] = 0;      //x
#endif

  
  count[0] = 1;
#if LAYERS
  count[1] = nl;
  count[2] = NYp1;
  count[3] = Nxp1;
#else
  count[1] = NYp1;
  count[2] = Nxp1;
#endif  

  for (int iv = 0; iv < list_nc[0].len; iv++){

    // TODO: FOR MPI
    for (int k = 0; k < nl; k++) {
      for (int j = 0; j < NYp1; j++) {
        for (int i = 0; i < Nxp1; i++) {
//          field[NYp1*Nxp1*k + Nxp1*j + i] = nodata; // for MPI
          field[NYp1*Nxp1*k + Nxp1*j + i] = 0.;
        }
      }
    }


    // TODO: FOR MPI
    double * data_loc = (double*)list_nc[iv].data;
    for (int k = 0; k < nl; k++) {
      for (int j = 0; j < NYp1; j++) {
        for (int i = 0; i < Nxp1; i++) {
          field[NYp1*Nxp1*k + Nxp1*j + i] = data_loc[NYp1*Nxp1*k + Nxp1*j + i];
        }
      }
    }


    
/*     if (pid() == 0) { // master */
/* @if _MPI */
/*         MPI_Reduce (MPI_IN_PLACE, &field[0], N_out*N_out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD); */
/* @endif */
  

     if ((nc_err = nc_put_vara_float(ncid, nc_varid[iv], start, count,
        			      &field[0])))
         ERR(nc_err);

  /* } // master*/
/* @if _MPI */
/*   else // slave */
/*   MPI_Reduce (&field[0], NULL, N_out*N_out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD); */
/* @endif */
  }
  free(field);


   /* Close the file. */
  /* if (pid() == 0) { // master */
    if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
  /* } // master */
//   printf("*** SUCCESS writing example file %s -- %d!\n", file_nc, nc_rec);
}

#ifdef _MPI

void gather_info(){
   // First rank gathers information about the other processes, to be able to properly handle 
   // their data upon communication.

   size_gather = malloc( n_ranks*sizeof( int ) ); // Number of data to be sent from each rank
   start_gather = malloc( n_ranks*sizeof( int ) ); // Row at which the domain of each rank starts
   rows_gather = malloc( n_ranks*sizeof( int ) ); // Number of rows to be sent from each rank
   
   if (rank == 0 && rank == n_ranks-1){ // Only one rank
      size_gather_local = Nyp1*Nxp1*nl;
      Ny_send_start = Ny_start - 1;
      Ny_send_rows = Nyp1;
   } else if( rank == 0 ){ // First rank has one row more (southern boundary)
      size_gather_local = Ny*Nxp1*nl;
      Ny_send_start = Ny_start - 1;
      Ny_send_rows = Ny;
   } else if ( rank == n_ranks - 1 ) { // Last rank has one row more (northern boundary)
      size_gather_local = Ny*Nxp1*nl;
      Ny_send_start = Ny_start;
      Ny_send_rows = Ny;
   } else { 
      size_gather_local = Nym1*Nxp1*nl;
      Ny_send_start = Ny_start;
      Ny_send_rows = Nym1;
   }

   MPI_Gather(&size_gather_local, 1, MPI_INT, size_gather, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Gather(&Ny_send_start, 1, MPI_INT, start_gather, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Gather(&Ny_send_rows, 1, MPI_INT, rows_gather, 1, MPI_INT, 0, MPI_COMM_WORLD);
   
}

void gather_output(){
   // function to gather information from other ranks before writing output
   // Rank 0 will gather information and write the nc_file
   
   if (rank ==0){ 
      // Copy own data chunk into psi and q arrays
      if (n_ranks == 1){ // No reception if rank 0 is the only ranks
         for (int k = 0; k < nl; k++){
            for (int j=0; j < Nyp1; j++){
               for (int i=0; i < Nxp1; i++){
                  psi_out[idx(i,j,k)] = psi[idx(i,j,k)];
               }
            }
         }
         for (int k = 0; k < nl; k++){
            for (int j=0; j < Nyp1; j++){
               for (int i=0; i < Nxp1; i++){
                  q_out[idx(i,j,k)] = q[idx(i,j,k)];
               }
            }
         }
      } else { // Else rank 0 also has to receive data from all other ranks
         for (int k = 0; k < nl; k++){
            for (int j=0; j < Ny; j++){
               for (int i=0; i < Nxp1; i++){
                  int idx_loc = k*Nxp1*Nyp1 + Nxp1*j + i;
                  int idx_glob = k*Nxp1*Nytp1 + j*Nxp1 + i;
                  psi_out[idx_glob] = psi[idx_loc];
               }
            }
         }
         for (int k = 0; k < nl; k++){
            for (int j=0; j < Ny; j++){
               for (int i=0; i < Nxp1; i++){
                  int idx_loc = k*Nxp1*Nyp1 + Nxp1*j + i;
                  int idx_glob = k*Nxp1*Nytp1 + j*Nxp1 + i;
                  q_out[idx_glob] = q[idx_loc];
               }
            }
         }
         
         // receive remaining chunks and copy into arrays
         for (int ii = 1;ii <n_ranks; ii++){
                     
            // We reallocate reception arrays at every communication, as the number of rows and hence data points 
            // might be different depending on the rank that we are receiving from.
            double *recv_psi;
            double *recv_q;
               
            // receive psi arrays
            MPI_Status  status1;
            recv_psi = calloc( size_gather[ii], sizeof( double ) );
            MPI_Recv(recv_psi, size_gather[ii], MPI_DOUBLE, ii, 1, MPI_COMM_WORLD, &status1);
            
            // Copy psi into output array
            int j_start = start_gather[ii];
            int j_end = start_gather[ii] + rows_gather[ii]; // starting point + number of rows

            for (int k = 0; k < nl; k++){
               for (int j= j_start; j < j_end; j++){
                  for (int i=0; i < Nxp1; i++){
                     int idx_loc = k*Nxp1*rows_gather[ii] + Nxp1*(j - j_start) + i;
                     int idx_glob = k*Nxp1*Nytp1 + j*Nxp1 + i;
                     psi_out[idx_glob] = recv_psi[idx_loc];
                  }
               }
            }

            // receive q arrays
            recv_q = calloc( size_gather[ii], sizeof( double ) );
            MPI_Status  status2;
            MPI_Recv(recv_q, size_gather[ii], MPI_DOUBLE, ii, 1, MPI_COMM_WORLD, &status2);
            
            // Copy q into output array
            for (int k = 0; k < nl; k++){
               for (int j= j_start; j < j_end; j++){
                  for (int i=0; i < Nxp1; i++){
                     int idx_loc = k*Nxp1*rows_gather[ii] + Nxp1*(j - j_start) + i;
                     int idx_glob = k*Nxp1*Nytp1 + j*Nxp1 + i;
                     q_out[idx_glob] = recv_q[idx_loc];
                  }
               }
            }

            free(recv_psi);
            free(recv_q);
         }  
      }
   } else { // send arrays from other ranks
      
      double *send_psi;
      double *send_q;
      
      send_psi = calloc( size_gather_local, sizeof( double ) );

      // copy psi arrays into sending allocation
      for (int k = 0; k < nl; k++){
         for (int j= 1; j < Ny_send_rows+1; j++){
            for (int i=0; i < Nxp1; i++){
               int idx_loc = k*(Nxp1*Ny_send_rows) + Nxp1*(j-1) + i;
               send_psi[idx_loc] = psi[idx(i,j,k)];
            }
         }
      }

      MPI_Send(send_psi, size_gather_local, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      send_q = calloc( size_gather_local, sizeof( double ) );

      // copy q arrays into sending allocation
      for (int k = 0; k < nl; k++){
         for (int j= 1; j < Ny_send_rows+1; j++){
            for (int i=0; i < Nxp1; i++){
               int idx_loc = k*(Nxp1*Ny_send_rows) + Nxp1*(j-1) + i;
               send_q[idx_loc] = q[idx(i,j,k)];
            }
         }
      }

      MPI_Send(send_q, size_gather_local, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      free(send_psi);
      free(send_q);
   }


} 

#endif