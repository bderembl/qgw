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
#define ERR(e) {fprintf(stdout,"Error: %s\n", nc_strerror(e)); return;}

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

/* output indices for periodic/bounded simulations */
int NX_Out;
int NY_Out;

// flag to remove boundary points in periodic domains
int ibc = 0;

void create_nc(char* file_out)
{

  // periodic: remove one point in all directions in output
   if (bc_fac == -1) {
     ibc = 1;
   }

   NX_Out = NXp2 - 2*ibc;
   NY_Out = NYp2 - 2*ibc;


   sprintf (file_nc,"%s", file_out);
   if (pid() == 0) { // master

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
   if ((nc_err = nc_def_dim(ncid, Y_NAME, NY_Out, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, NX_Out, &x_dimid)))
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
   float yc[NY_Out], xc[NX_Out];
   for (int i = 0; i < NX_Out; i++){
      xc[i] = i*Delta;
   }
   for (int i = 0; i < NY_Out; i++){
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
  } // end master

}

/** Write NetCDF file

*/



void write_nc() {

  if (pid() == 0) { // master
    /* open file. */
    if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid)))
      ERR(nc_err);
  }

  // write time
  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;
  if (pid() == 0) { // master */
    if ((nc_err = nc_put_vara_float(ncid, t_varid, startt, countt,
                                    &loctime)))
      ERR(nc_err);
  } // master */



  float * field = (float *)malloc(NX_Out*NY_Out*nl*sizeof(float));

  
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
  count[2] = NY_Out;
  count[3] = NX_Out;
#else
  count[1] = NY_Out;
  count[2] = NX_Out;
#endif  

  for (int iv = 0; iv < list_nc[0].len; iv++){

      for (int k = 0; k < nl; k++) {
        for (int j = 0; j < NY_Out; j++) {
          for (int i = 0; i < NX_Out; i++) {
            field[NY_Out*NX_Out*k + NX_Out*j + i] = nodata; // for MPI
  //          field[Nyp2*Nxp2*k + Nxp2*j + i] = 0.;
          }
        }
      }

      double * data_loc = (double*)list_nc[iv].data;
    // TODO: FOR MPI
#if _MPI 
    if (rank <= rank_crit) // only subcritical ranks are allowed to write their data
#endif

      for (int k = 0; k < nl; k++) {
        for (int j = ibc; j < Nyp2 - ibc; j++) {
          for (int i = ibc; i < Nxp2 - ibc; i++) {
            field[NY_Out*NX_Out*k + NX_Out*(j + J0) + (i + I0)] = data_loc[idx(i,j,k)];
          }
        }
      }
    
    if (pid() == 0) { // master

#if _MPI
      MPI_Reduce (MPI_IN_PLACE, &field[0], NX_Out*NY_Out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif

      if ((nc_err = nc_put_vara_float(ncid, nc_varid[iv], start, count,
        			       &field[0])))
          ERR(nc_err);
  } // master
#if _MPI
  else // slave
  MPI_Reduce (&field[0], NULL, NX_Out*NY_Out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif
  }
  free(field);


   /* Close the file. */
  if (pid() == 0) { // master
    if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
  } // master
//   printf("*** SUCCESS writing example file %s -- %d!\n", file_nc, nc_rec);
}


/**
   Read NC

   TODO: extend this routine to read anything
 */

void read_nc(List *list_in, char* file_in, int rec_in){

  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  nc_type type;
  char varname[NC_MAX_NAME+1];
  int *dimids = NULL;


  // periodic: remove one point in all directions in output
   if (bc_fac == -1) {
     ibc = 1;
   }

   NX_Out = NXp2 - 2*ibc;
   NY_Out = NYp2 - 2*ibc;

  float * field = (float *)malloc((NX_Out)*(NY_Out)*nl*sizeof(float));

  if ((nc_err = nc_open(file_in, NC_NOWRITE, &ncfile)))
    ERR(nc_err);

  if ((nc_err = nc_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited)))
    ERR(nc_err);

  for (int iv_list = 0; iv_list < list_in[0].len; iv_list++){


  for(int iv=0; iv<nvars; iv++) {


    if ((nc_err = nc_inq_var(ncfile, iv, varname, &type, &var_ndims, dimids,
                             &var_natts)))
      ERR(nc_err);


    if (strcmp(varname,list_in[iv_list].name) == 0) {
      fprintf(stdout,"Reading variable  %s!\n", varname);

      double * data_loc = (double*)list_in[iv_list].data;

      size_t start[4], count[4];
      start[0] = rec_in; //time
      start[1] = 0;
      start[2] = 0;
      start[3] = 0;

      count[0] = 1;
      count[1] = nl;
      count[2] = NX_Out;
      count[3] = NY_Out;
      if ((nc_err = nc_get_vara_float(ncfile, iv, start, count,
                                      &field[0])))
        ERR(nc_err);


      for (int k = 0; k < nl; k++) {
        for (int j = ibc; j < Nyp2 - ibc; j++) {
          for (int i = ibc; i < Nxp2 - ibc; i++) {
            data_loc[idx(i,j,k)] = field[NY_Out*NX_Out*k + NX_Out*(j + J0) + i + I0];
          }
        }
      }
    }

  } // end nvar loop

  } // end list loop

  free(field);

  if ((nc_err = nc_close(ncfile)))
    ERR(nc_err);
  
}
