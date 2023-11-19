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
  printf( "\n %s\n", desc );
 
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( " %6.2e", mat[i*ldm+j] );
    printf( "\n" );
  }
}
 
 
/* Auxiliary routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
  lapack_int i, j;
  printf( "\n %s\n", desc );
 
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i+j*ldm] );
    printf( "\n" );
  }
}
 
/* Auxiliary routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
  lapack_int j;
  printf( "\n %s\n", desc );
  for( j = 0; j < n; j++ ) printf( " %6i", vec[j] );
  printf( "\n" );
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
void eigmod ()
{
  
  fprintf(stdout,"Compute vertical modes .. ");
  double htotal = 0;
  for (int k = 0; k < nl ; k++)
    htotal += dh[k];
    

  int print = 0;
//  int print2 = 0;
    
    // initialize matrix
    double * amat;
    amat = malloc (nl*nl*sizeof(double));

    for (int m = 0; m < nl ; m++) {
      for (int k = 0; k < nl ; k++) {
        amat[nl*k + m] = 0.;
      }
    }
     
    if (nl > 1){

      int l = 0;
      amat[nl*l + l+1] = -sq(f_cor)/N2[l]/( dhc[l]*dh[l]);
      amat[nl*l + l ] = - amat[nl*l + l+1];

      for (int l = 1; l < nl-1 ; l++) {
        amat[nl*l+l-1] = -sq(f_cor)/N2[l-1]/( dhc[l-1]*dh[l]);
        amat[nl*l+l+1] = -sq(f_cor)/N2[l]/( dhc[l]*dh[l]);
        amat[nl*l+ l ] = - amat[nl*l+l-1] - amat[nl*l+l+1];
      }
      l = nl-1;
      amat[nl*l + l-1] = -sq(f_cor)/N2[l-1]/( dhc[l-1]*dh[l]);
      amat[nl*l + l ] = - amat[nl*l + l-1];
    }
    else{
      amat[0] = 0.;
    }

    /* if (print2){ */

    /*   printf("fo = %g\n",f_cor); */
    /*   for (int l = 0; l < nl ; l++) { */
    /*     scalar Fr1 = Frl[l]; */
    /*     printf("%i, Fr = %g\n",l, Fr1[]); */
    /*     printf("%i, dh = %g\n",l, dh[l]); */
    /*   } */
    /*   print_matrix_rowmajor( "a matrix", nl,nl, amat, nl); */
    /*   print2 = 0; */
    /* } */

/*     // store amat in streching matrix (tridiagonal matrix) */
/*     for (int l = 0; l < nl ; l++) { */
/*       scalar sm0 = strechmat[3*l]; */
/* //      scalar sm1 = strechmat[3*l+1]; */
/*       scalar sm2 = strechmat[3*l+2]; */

/*       sm0[] = nl*l+l-1 < 0 ? 0. : amat[nl*l+l-1]/sq(Ro[]); */
/* //      sm1[] = amat[nl*l+ l ]/sq(Ro[]); */
/*       sm2[] = nl*l+ l+1 > nl-1 ? 0. : amat[nl*l+l-1]/sq(Ro[]); */

    /* } */


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

   int info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'V', 'V', nl, amat, nl, wr, wi,
                              vl, nl, vr, nl );

    if (info < 0) {
      printf("issue with lapack in eigmode.h\n");
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

    if (print) {
    print_matrix_rowmajor( "Right eigenvectors", nl, nl, vr, nl );
    print_matrix_rowmajor( "Left eigenvectors", nl, nl, vl, nl );
    }

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
    // set BT mode to zero
    iRd2[0] = 0.;

    if (print)
      {
    /* for (int m = 0; m < nl ; m++) { */
    /*   for (int k = 0; k < nl ; k++) { */
    /*     scalar l2m = cl2m[m*nl+k]; */
    /*     scalar m2l = cm2l[m*nl+k]; */
    /*     vl[m*nl+k] = l2m[]; */
    /*     vr[m*nl+k] = m2l[]; */
    /*   } */
    /* } */
    print_matrix_rowmajor( "Right eigenvectors (renorm)", nl, nl, cm2l, nl );
    print_matrix_rowmajor( "Left eigenvectors (renorm)", nl, nl, cl2m, nl );
    for (int k = 0; k < nl ; k++) {
      printf("vp: %g %i\n", wr[k], wr2[k].index);
      printf("iRd2: %g , def radius: %g\n", iRd2[k], iRd2[k] > 0 ? sqrt(1/iRd2[k]) : 0);
    }
    print = 0;
      }
    

  free(amat);
  free(wr);
  free(wr2);
  free(wi);
  free(vl);
  free(vr);
  free(tmp);


  fprintf(stdout,"ok\n ");
}

