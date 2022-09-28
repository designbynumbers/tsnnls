/*

  cblas_wrappers.c : Wrapper functions which translate between the old-school LAPACK
  interface that Taucs expects and the new cblas calling style implemented by OpenBlas.

*/

#include<cblas.h>
#include<lapacke.h>

CBLAS_UPLO uplo_convert(char *uplo)
{
  
  if(*uplo == 'U' || *uplo == 'u') {

    return CblasUpper;
    
  } else if (*uplo == 'L' || *uplo == 'l') {

    return CblasLower;
    
  } else {

    fprintf(stderr,"uplo_convert error: called with uplo = %c != U|u|L|l \n",*uplo);
    exit(1);

  }
}

CBLAS_TRANSPOSE trans_convert(char *trans)
{
  
  if (*trans == 'T' || *trans == 't') {
    return CblasTrans;
  } else if (*trans == 'N' || *trans == 'n') {
    return CblasNoTrans;
  } else if (*trans == 'C' || *trans == 'c') {
    return CblasConjTrans;
  } else {
    fprintf(stderr,"trans_convert error: called with trans = %c != T|t|N|n|C|c\n",*trans);
    exit(1);
  }
  
}

CBLAS_SIDE side_convert(char *side)
{
  
  if (*side == 'L' || *side == 'l') {
    return CblasLeft;
  } else if (*side == 'R' || *side == 'r') {
    return CblasRight;
  } else {
    fprintf(stderr,"side_convert error: called with side = %c != L|l|R|r\n",*side);
    exit(1);
  }
  
}

CBLAS_DIAG diag_convert(char *diag)
{

  if (*diag == 'U' || *diag == 'u') {
    return CblasUnit;
  } else if (*diag == 'N' || *diag == 'n') {
    return CblasNonUnit;
  } else {
    fprintf(stderr,"diag_convert error: called with diag = %c != U|u|N|n\n",*diag);
    exit(1);
  }
  
}

int wrapper_dsyrk(char *uplo, char *trans, 
		  int *N, int *K, 
		  double* alpha, 
		  double* A, int *lda, 
		  double* beta, 
		  double* C, int *ldc)

{

  cblas_dsyrk(CblasColMajor,uplo_convert(uplo),trans_convert(trans),*N,*K,*alpha,A,*lda,*beta,C,*ldc);

  return(0);

}

int wrapper_dgemm(char *transA, char *transB, int*M, int*N, int *K,
		  double*alpha, double*A, int *lda, double*B, int *ldb, 
		  double*beta, double*C, int*ldc)
{
  
  cblas_dgemm(CblasColMajor,trans_convert(transA),trans_convert(transB),*M,*N,*K,*alpha,A,
	      *lda,B,*ldb,
	      *beta, C,*ldc);

  return(0);
}

int wrapper_trsm(char *side, char *uplo, char *transA, char *diag, 
		 int *M, int*N, double*alpha, double*A, int *lda, 
		 double*B, int *ldb)
{
  cblas_dtrsm(CblasColMajor, side_convert(side),
              uplo_convert(uplo), trans_convert(transA),
              diag_convert(diag), *M, *N,
              *alpha, A, *lda,
              B, *ldb);

  return 0;
}

double wrapper_dnrm2(int*N, double*X, int*incX)
{
  return cblas_dnrm2(*N,X,*incX);
}

int wrapper_potrf(char* uplo, int*n, double*a, int*lda, int*info)
{
  *info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,*uplo,*n,a,*lda);
  return *info;
}

