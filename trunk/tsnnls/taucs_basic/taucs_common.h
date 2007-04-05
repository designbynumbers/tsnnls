/* Taucs_common.

   This file added to the taucs_basic distribution in TSNNLS to
   include the dmalloc header in a TAUCS compile cycle, hoping to pick
   up a subtle porting bug. 

   Jason Cantarella. 4/1/2007 */

#ifndef TAUCS_COMMON
#define TAUCS_COMMON

#include <config.h>

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#define taucs_herk  DSYRK_F77
#define taucs_gemm  DGEMM_F77
#define taucs_trsm  DTRSM_F77
#define taucs_potrf DPOTRF_F77
#define taucs_dnrm2 DNRM2_F77

extern int taucs_potrf(char*, int*, taucs_datatype*, int*, int*);
extern int taucs_trsm(char *, char *, char *, char *, 
			int*, int*, taucs_datatype*, taucs_datatype*, int *, 
			taucs_datatype*, int *);
extern int taucs_gemm(char *, char *, int*, int*, int *,
			taucs_datatype*, taucs_datatype*, int *, taucs_datatype*, int *, 
			taucs_datatype*, taucs_datatype*, int*);

extern int taucs_herk(char *, char *, 
		      int *, int *, 
		      taucs_real_datatype*, 
		      taucs_datatype*, int *, 
		      taucs_real_datatype*, 
		      taucs_datatype*, int *);

//taucs_double taucs_blas_name(dnrm2)(int*, taucs_double*, int*);
//taucs_single taucs_blas_name(snrm2)(int*, taucs_single*, int*);
//taucs_double taucs_blas_name(dznrm2)(int*, taucs_dcomplex*, int*);
//taucs_single taucs_blas_name(scnrm2)(int*, taucs_scomplex*, int*);

taucs_double taucs_dnrm2(int*, taucs_double*, int*);
taucs_single taucs_blas_name(snrm2)(int*, taucs_single*, int*);
taucs_double taucs_blas_name(dznrm2)(int*, taucs_dcomplex*, int*);
taucs_single taucs_blas_name(scnrm2)(int*, taucs_scomplex*, int*);

#endif
