/**
   Compute eigenvectors and eigenvalue of the vertical streching
   matrix
*/

#ifndef MKL
#include <lapacke.h>
#else
#include <mkl_lapacke.h>
#endif
#include <stdio.h>

double *cl2m;
double *cm2l;
double *iRd2;

void init_eigmode(){

  for (int l = 0; l < nl-1; l++)
    dhc[l] = 0.5*(dh[l] + dh[l+1]);

  cl2m = calloc( nl*nl, sizeof( double ) );
  cm2l = calloc( nl*nl, sizeof( double ) );
  iRd2 = calloc( nl, sizeof( double ) );

}

void clean_eigmode(){
  free(cl2m);
  free(cm2l);
  free(iRd2);
}


/* Auxiliary routine: printing a matrix */
void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
  lapack_int i, j;
  fprintf(stdout, "\n %s\n", desc );
 
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) fprintf(stdout, " %6.2e", mat[i*ldm+j] );
    fprintf(stdout, "\n" );
  }
}
 
 
/* Auxiliary routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
  lapack_int i, j;
  fprintf(stdout, "\n %s\n", desc );
 
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) fprintf(stdout, " %6.2f", mat[i+j*ldm] );
    fprintf(stdout, "\n" );
  }
}
 
/* Auxiliary routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
  lapack_int j;
  fprintf(stdout, "\n %s\n", desc );
  for( j = 0; j < n; j++ ) fprintf(stdout, " %6i", vec[j] );
  fprintf(stdout, "\n" );
}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

struct array_i {
  double data;
  int index; 
};

int compare_i (const void* a, const void* b) {
  return (((struct array_i*)a)->data > ((struct array_i*)b)->data) 
    - (((struct array_i*)a)->data < ((struct array_i*)b)->data);
}


/**
 Main program */
void compute_eigmode () {
  if (print) {fprintf(stdout,"Compute vertical modes .. ");}

  double htotal = 0;
  for (int k = 0; k < nl ; k++)
    htotal += dh[k];
    
  
  // initialize matrix
  double * Gamma_mat;
  Gamma_mat = malloc (nl*nl*sizeof(double));

  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      Gamma_mat[nl*k + m] = 0.;
    }
  }
    
  if (nl > 1){

    int l = 0;
    Gamma_mat[nl*l + l+1] = -sq(f0)/N2[l]/( dhc[l]*dh[l]);
    Gamma_mat[nl*l + l ] = - Gamma_mat[nl*l + l+1];

    for (int l = 1; l < nl-1 ; l++) {
      Gamma_mat[nl*l+l-1] = -sq(f0)/N2[l-1]/( dhc[l-1]*dh[l]);
      Gamma_mat[nl*l+l+1] = -sq(f0)/N2[l]/( dhc[l]*dh[l]);
      Gamma_mat[nl*l+ l ] = - Gamma_mat[nl*l+l-1] - Gamma_mat[nl*l+l+1];
    }
    l = nl-1;
    Gamma_mat[nl*l + l-1] = -sq(f0)/N2[l-1]/( dhc[l-1]*dh[l]);
    Gamma_mat[nl*l + l ] = - Gamma_mat[nl*l + l-1];
  }
  else{
    Gamma_mat[0] = 1.;
  }

  if (print){
    print_matrix_rowmajor( "Gamma matrix", nl,nl, Gamma_mat, nl);
  }

  /**
     Declare variables for lapack
    */
  double * wr;
  double * wi;
  double * vl;
  double * vr;
  double * tmp;
  vr = malloc (nl*nl*sizeof(double));
  vl = malloc (nl*nl*sizeof(double));
  tmp = malloc (nl*nl*sizeof(double));
  
  wr = malloc (nl*sizeof(double));
  wi = malloc (nl*sizeof(double));

  struct array_i* wr2 = (struct array_i*) malloc(nl * sizeof(struct array_i));

  int info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'V', 'V', nl, Gamma_mat, nl, wr, wi,
                            vl, nl, vr, nl );

  if (info < 0) {
    fprintf(stdout,"issue with lapack in eigmode.h\n");
    exit(0);
  }

  // sort
  for (int l = 0; l < nl ; l++) {
    wr2[l].data = wr[l]; 
    wr2[l].index = l;
  }

  qsort(wr2, nl, sizeof(struct array_i), &compare_i);

  for (int l = 0; l < nl ; l++) {
    wr[l] = wr2[l].data;
  }   

  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      tmp[k*nl+m] = vr[k*nl+m];
    }
  }
  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      vr[k*nl+m] = tmp[k*nl+wr2[m].index];
    }
  }

  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      tmp[k*nl+m] = vl[k*nl+m];
    }
  }
  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      vl[k*nl+m] = tmp[k*nl+wr2[m].index];
    }
  }
  /*
  if (print) {
  print_matrix_rowmajor( "Right eigenvectors", nl, nl, vr, nl );
  print_matrix_rowmajor( "Left eigenvectors", nl, nl, vl, nl );
  }
  */
  /**
     Apply the normalisation of Flierl (1978) to the oceanic
      right eigenvectors, and ensure they are +ve at the surface
      ----------------------------------------------------------
      Normalisation is equation (2.5) of "Models of Vertical
      Structure and the Calibration of Two-layer Models" by Glenn
      R. Flierl, Dynamics of Atmospheres and Oceans vol. 2 pp
      341-381 (1978) Surface sign convention follows the Rossby
      wave papers of P.D. Killworth & J.R. Blundell in J. Physical
      Oceanography.  This is only immediately relevant for the
      ocean case.
  */

  for (int m = 0; m < nl ; m++) {
    double dotp = 0.;
    for (int k = 0; k < nl ; k++) {
      dotp += dh[k]*vr[k*nl+m]*vr[k*nl+m];
    }
    double flfac = sign(vr[m])*sqrt(htotal/dotp);
    for (int k = 0; k < nl ; k++) {
      vr[k*nl+m] = flfac*vr[k*nl+m];
    }
  }
  for (int m = 0; m < nl ; m++) {
    double dotp = 0.;
    for (int k = 0; k < nl ; k++) {
      dotp += vr[k*nl+m]*vl[k*nl+m];
    }
    for (int k = 0; k < nl ; k++) {
      vl[k*nl+m] = vl[k*nl+m]/dotp;
    }
  }

  /**
     Store matrices and deformation radius 
  */

  for (int m = 0; m < nl ; m++) {
    for (int k = 0; k < nl ; k++) {
      cl2m[k*nl+m] = vl[m*nl+k];
      cm2l[k*nl+m] = vr[k*nl+m]; // transpose matrix

      /* l2m[] = (k == m ? 1.0 : 0.); */
      /* m2l[] = (k == m ? 1.0 : 0.); */
    }
  }

  // store 1/Rd2 
  for (int k = 1; k < nl ; k++) {
    iRd2[k] = wr[k];
  }

  // set BT mode to zero (unless it's a 1.5 layer simulation)
  iRd2[0] = 0.;
  if ((nl == 1) && (Ld != 0))
    iRd2[0] = 1/sq(Ld);

  
  if (print) {

    print_matrix_rowmajor( "Right eigenvectors (renorm)", nl, nl, cm2l, nl );
    print_matrix_rowmajor( "Left eigenvectors (renorm)", nl, nl, cl2m, nl );
    fprintf(stdout,"\n");
    fprintf(stdout,"Deformation radii: \n");
    for (int k = 0; k < nl ; k++) {
      fprintf(stdout,"iRd2: %g , def radius: %g\n", iRd2[k], iRd2[k] > 0 ? sqrt(1/iRd2[k]) : 0);
    }
  }

  free(Gamma_mat);
  free(wr);
  free(wr2);
  free(wi);
  free(vl);
  free(vr);
  free(tmp);

  if (print) {fprintf(stdout,"End eigenmode computation\n ");}
  
}

