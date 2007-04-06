/*

   tsnnls_blas_wrappers.h : This (non-installed) header contains
   prototypes for the blas functions linked to by tsnnls. If you are
   having trouble linking tsnnls with a blas library, you may need to
   redefine these symbols manually.

*/

extern void DAXPY_F77(int *N,double *alpha,double *vecX, \
		      int *incX,double *vecY,int *incY);

extern void DSYMV_F77(char *uplo,int *N,double *alpha, \
		      double *matA, int *lda, double *vecX, int *incX, \
		      double *beta, double *vecY, int *incY);

extern void DSYRK_F77(char *uplo,char *trans,int *N,int *K, \
		      double *alpha, double *matA, int *ldA, \
		      double *beta, double *matC, int *ldC);

extern double DDOT_F77(int *N, double *vecX, int *incX, double *vecY, int *incY);

extern void DSCAL_F77(int *N, double *alpha, double *vecX, int *incX);

extern double DNRM2_F77(int *N, double *vecX, int *incX);

/* Subroutine */ 
//int dpocon_(char *uplo, __CLPK_integer *n, 
//	    __CLPK_doublereal *a, __CLPK_integer *lda, 
//	    __CLPK_doublereal *anorm, __CLPK_doublereal *rcond, 
//	    __CLPK_doublereal *work, __CLPK_integer *iwork, 
//	    __CLPK_integer *info);
