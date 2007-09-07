/*

  genb_test.c : This file is part of a test suite for the TSNNLS
  distribution. This test loads a matrix (or randomly generates one,
  if none is given), and constructs a series of test NNLS problems
  based on that matrix, using the GENB algorithm given in PJV.

  The matrix must be of full rank, since we use a QR factorization as
  part of the test generation. We are also required to have a full
  CLAPACK, so this test won't run on a machine with only the minimal
  ATLAS LAPACK. The test can also be quite slow if the matrix loaded
  is very large.

  Note that this file isn't even made if we don't have a full LAPACK.
*/

#define NUM_TESTS 512
#define Msize 152
#define Nsize 67

#include<config.h>

#include"acint32_type.h"
#include"tsnnls_blas_wrappers.h"

#ifdef HAVE_STDLIB_H
 #include<stdlib.h>
#endif

#ifdef HAVE_MATH_H
 #include<math.h>
#endif

#ifdef HAVE_STDIO_H
 #include<stdio.h>
#endif

#ifdef HAVE_TIME_H
 #include<time.h>
#endif

#ifdef HAVE_CLAPACK_H
  #include <clapack.h>
#else 
  #ifdef HAVE_ATLAS_CLAPACK_H
     #include <atlas/clapack.h>
  #else
     #ifdef HAVE_VECLIB_CLAPACK_H
       #include <vecLib/clapack.h>
     #endif
  #endif
#endif

/* We now arrange to have a copy of Lapack's DGELS to link to. 
  
   If we have a full (not just Atlas) CLAPACK, we'll use the 
   clapack interface to call dgels_. In this case no prototype
   is needed, because we just included clapack.h.

   If we have a full fortran LAPACK, but no clapack, we will have used
   our mangling scheme to define DGELS_F77 in the configure
   process. This will link to the fortran DGELS, but doesn't have a
   prototype. So we add an extern declaration here to prevent compiler
   warnings.

   Note that we always use fortran BLAS calls, even when we could use 
   cblas to avoid doing this kind of thing elsewhere.

*/


#include"tsnnls.h"

#ifdef HAVE_FULL_CLAPACK
  #define DGELS_WRAPPER dgels_  
#else
  #define DGELS_WRAPPER DGELS_F77

  extern void DGELS_F77(char *trans, ACINT32_TYPE *M, ACINT32_TYPE *N, ACINT32_TYPE *NRHS,
			double *A, ACINT32_TYPE *ldA, double *B, ACINT32_TYPE *ldB,
			double *work, ACINT32_TYPE *lwork, ACINT32_TYPE *info);
#endif


double *genb(int mdim, int ndim, double *A, double *x, double *y)
 
/* Use the method PJV to generate a rhs for a test problem with 
   matrix A and solution x and y. We expect that 

   A is an m x n matrix stored in LAPACK fortran order, 
   x is an n-vector with non-negative entries, 
   y is an n-vector with non-negative entries

   we return an m-vector b, which we allocate. */

{
  int     i;
  ACINT32_TYPE     m = mdim, n = ndim;
  
  /* We are more concerned with readability than speed, since this 
     is just debugging code. So we replace the CSNE method of PJV with
     a simpler iterative refinement algorithm:

     1) Approximately solve A^T z = y using the QR decomposition via DGELS.
     2) Compute r = A^T z - y.
     3) Again solve A^T z' = r using DGELS.
     4) Set the final z to z - z'.    
     5) Compute b = Ax - z.

  */

  char   trans = 'T';
  ACINT32_TYPE    nrhs = 1;
  double *z = calloc(m,sizeof(double));
  for (i=0;i<n;i++) { z[i] = y[i]; }

  double *Awork = calloc(m*n,sizeof(double));
  for(i=0;i<m*n;i++) { Awork[i] = A[i]; }

  double *work = calloc(128*m*n,sizeof(double));
  ACINT32_TYPE    lwork = 128*m*n;
  ACINT32_TYPE    info;

  if (work == NULL) { 

    fprintf(stderr,"genb_test failed: Couldn't allocate worksize %d.\n",
	    (int)(128*m*n));
    exit(1);

  }

  DGELS_WRAPPER(&trans,&m,&n,&nrhs,Awork,&m,z,&m,work,&lwork,&info);

  if (info != 0) { 

    fprintf(stderr,"genb_test: DGELS call 1 failed due to arg %d.\n",-(int)(info));
    exit(1);

  }

  /* The vector z now contains the putative solution. We compute the
     residual A^T z - y. */

  double *r = calloc(m,sizeof(double));
  for(i=0;i<n;i++) { r[i] = y[i]; }

  trans = 'T';
  double alpha = 1.0, beta = -1.0;
  int incX = 1, incY = 1;
  int intM = m, intN = n;
  
  DGEMV_F77(&trans,&intM,&intN,&alpha,A,&intM,z,&incX,&beta,r,&incY);

  /* Having computed the residual r, we record the norm of the residual for 
     debugging purposes. */

  double resNorm = 0;

  for(i=0;i<n;i++) { resNorm += pow(r[i],2.0); }
  resNorm = sqrt(resNorm);
  
  /* We now solve A^T delta = r to find a correction step delta. */

  for (i=0;i<m*n; i++) { Awork[i] = A[i]; } /* Reset AWork. */
  
  double *delta = calloc(m,sizeof(double));
  for(i=0;i<n;i++) { delta[i] = r[i]; }

  DGELS_WRAPPER(&trans,&m,&n,&nrhs,Awork,&m,delta,&m,work,&lwork,&info);

  if (info != 0) { 

    fprintf(stderr,"genb_test: DGELS call 2 failed due to arg %d.\n",(int)(-info));
    exit(1);

  }
  
  /* We now subtract the correction step delta from z to find a better z. */

  alpha = -1.0; beta = 1.0;
  DAXPY_F77(&intM,&alpha,delta,&incX,z,&incY);

  /* We have now found a very high-quality z so that A^T z = y. To check this, 
     we recompute the residual and its norm. */

  double *r2 = calloc(m,sizeof(double));
  for(i=0;i<m;i++) { r2[i] = y[i]; }

  trans = 'T';
  alpha = 1.0; beta = -1.0;
  incX = 1; incY = 1;
  intM = m; intN = n;
  
  DGEMV_F77(&trans,&intM,&intN,&alpha,A,&intM,z,&incX,&beta,r2,&incY);

  /* Having computed the residual r, we record the norm of the residual for 
     debugging purposes. */

  double resNorm2 = 0;

  for(i=0;i<n;i++) { resNorm2 += pow(r2[i],2.0); }
  resNorm2 = sqrt(resNorm2);

  /* Now that we have seen that the residual norm has improved, we can
     find b via Ax - z = b. */

  double *b = calloc(m,sizeof(double));
  for(i=0;i<m;i++) { b[i] = z[i]; }

  alpha = 1.0; beta = -1.0;
  trans = 'N';

  DGEMV_F77(&trans,&intM,&intN,&alpha,A,&intM,x,&incX,&beta,b,&incY);

  free(work); free(Awork); free(z); free(r); free(r2); free(delta);

  return b;

}
  
  
double *random_matrix(int m, int n)

/* Creates a random sparse-ish matrix for test purposes. In each column, we fill in 
   about 1/10 of the entries with nonzero elements. */

{
  double *A = calloc(m*n,sizeof(double));
  int i,j;

  for(j=0;j<n;j++) {

    for(i=0;i<ceil((double)(m)/10.0);i++) {

      A[(rand() % Msize) + Msize*j] = 2.0*((double)(rand())/(double)(RAND_MAX)) - 1.0;

    }

  }

  return A;

}

void random_x_y(int m, int n, double **x, double **y)

/* Generates random complementary and positive x and y vectors. */

{
  int i;
  
  *x = calloc(n,sizeof(double));
  *y = calloc(n,sizeof(double));

  for(i=0;i<n;i++) {

    if (rand() % 2 == 0) {  /* Half the time. */

      (*x)[i] = (double)(rand())/(double)(RAND_MAX);  /* x in [0,1] */

    } else {

      (*y)[i] = (double)(rand())/(double)(RAND_MAX);  /* y in [0,1] */

    }

  }

}

double *fliporder(int m, int n, double *A)

/* Changes order from the column-major ordering of LAPACK to the row-major
   ordering expected by tsnnls. */

{
  double *Aflip;
  int i, j;

  Aflip = calloc(Msize*Nsize,sizeof(double));
  
  for(i=0;i<Msize;i++) {
    
    for(j=0;j<Nsize;j++) {
      
      Aflip[j + Nsize*i] = A[i + Msize*j];  /* Flip for taucs_construct... */

    }
    
  }

  return Aflip;

}
  
int main() 

{
  double *x,*y,*A,*b,*Aflip;
  int     npass = 0;
  int     i,test;

  taucs_ccs_matrix *Accs;
  taucs_double *tsnnlsX;

  double ResNorm;
  double ErrTol = 0;
  int    PrintErrWarnings = 0;

  double err;

  fprintf(stderr,"genb_tests\n\n");
  fprintf(stderr,
	  "Creating %d random test %d x %d test problems using \n"
	  "the genb algorithm of PJV. Each problem will be given \n"
	  "to the tsnnls method, and the error printed below.\n",
	  NUM_TESTS,Msize,Nsize);

  fprintf(stderr,
	  "We require an error less than 1e-8 to pass the test.\n\n");

  fprintf(stderr,
	  "#    M    N     Error         Result \n"
	  "------------------------------------ \n");

#ifdef HAVE_TIME

  srand(time(0));

#else

  srand(536);

#endif

  for(test=0;test<NUM_TESTS;test++) {

    /* Create a random problem, with solution. */

    A = random_matrix(Msize,Nsize);
    random_x_y(Msize,Nsize,&x,&y);
    b = genb(Msize,Nsize,A,x,y);

    /* Now feed the problem to tsnnls. */

    Aflip = fliporder(Msize,Nsize,A);
    Accs = taucs_construct_sorted_ccs_matrix(Aflip,Nsize,Msize);
    /*tsnnlsX = t_snnls(Accs,b,&ResNorm,ErrTol,PrintErrWarnings);*/
    tsnnlsX = t_block3(Accs,b,&ResNorm,ErrTol,PrintErrWarnings);


    /* Now we compare the solution with our guess. */

    err = 0;
    for(i=0;i<Nsize;i++) { err += pow(tsnnlsX[i] - x[i],2.0); }
    err = sqrt(err);

    fprintf(stderr,"%3d  %-4d %-4d % 7e ",test+1,Msize,Nsize,err); 

    if (err < 1e-8) {

      npass++;
      fprintf(stderr," pass\n");

    } else {

      fprintf(stderr," FAIL\n");

    }

    free(A); free(x); free(y); free(b);
    taucs_ccs_free(Accs);

  }

  fprintf(stderr,"\n");
  fprintf(stderr,"%d (of %d) tests passed.\n",npass,NUM_TESTS);

  if (npass == NUM_TESTS) {

    exit(0);

  } else {

    exit(1);

  }

}
  

    

    

    
    
  
  
