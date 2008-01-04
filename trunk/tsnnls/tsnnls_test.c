

/*
 * This program is free software distributed under the GPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/gpl.txt
 */
 
/* tsnnls_test is a general command-line utility providing access to 
   some of the functionality of the tsnnls library. It is mostly useful
   for debugging tsnnls, but it can be used to quickly check the results
   of a MATLAB or Octave computation from the command line. */

/*
 * Note that for both formats, row and column indexing is 1-based.
 *
 * Sparse format is a plain text file containing:
 * 
 * SPARSE
 * row_number column_number
 * number_of_nonzeros
 * r1 c1 val11
 * r1 c2 val12
 * ...
 * r2 c1 val21
 * ...
 *
 * Dense format is a plain text file corresponding to the Octave text
 * format:
 *
 * Any number of lines starting with #, as long as they contain
 *
 * # rows: m
 * # columns: n 
 *
 * val1_1 val1_2 val1_3 val1_4 ... val1_col_number
 * val2_1 ...
 * ...
 * val_row_number_1 ...
 *
 * where valm_n is the value of the mth row and nth column.
 *
 * 
 * Representative files for both of these formats are included in the
 * directory 'test_files' in the distribution root directory.
 */

#include <config.h>

#ifdef HAVE_CTYPE_H
  #include <ctype.h>
#endif

#include "lsqr.h"
#include "tsnnls.h"

//#ifdef HAVE_CBLAS_H
//  #include <cblas.h>
//#else
//  #ifdef HAVE_VECLIB_CBLAS_H
//    #include <vecLib/cblas.h>
//  #else
//    #ifdef HAVE_ATLAS_CBLAS_H
//      #include <atlas/cblas.h>
//    #endif
//  #endif
//#endif

//#ifdef HAVE_DARWIN          /* We use the Apple BLAS/LAPACK if possible. */
// #include <vecLib/vBLAS.h>
//#else
// #include "tsnnls/gsl_cblas.h"
//#endif 

#ifdef HAVE_CLAPACK_H
  #include "clapack.h"
#else
  #ifdef HAVE_VECLIB_CLAPACK_H
    #include "vecLib/clapack.h"
  #else
    #ifdef HAVE_ATLAS_CLAPACK_H
      #include "atlas/clapack.h"
    #endif
  #endif
#endif

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef HAVE_SYS_TIME_H
  #include <sys/time.h>
#endif

#ifdef HAVE_SYS_RESOURCE_H
  #include <sys/resource.h>
#endif 

#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif 

#include "tsnnls_blas_wrappers.h"

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#ifdef HAVE_ARGTABLE2_H
  #include <argtable2.h>
#endif

int VERBOSITY = 0;
enum STYPE { tsnnls, pjv, tlsqr, fallback, SOLlsqr, spiv } solver = { tsnnls };
int nconstrained = -1;

int read_sparse( FILE *fp, double **vals, int *dim, int *cols );
int read_mat( FILE *fp, double **vals, int *dim, int *cols);

struct rusage start, end;
int tpStorage;

static double* lsqrwrapper( taucs_ccs_matrix* Af, double* b );

static void
start_clock()
{
	getrusage(RUSAGE_SELF, &start);
}

static double
end_clock()
{
    getrusage(RUSAGE_SELF, &end);

    double diff = end.ru_utime.tv_sec+end.ru_utime.tv_usec*1e-6 -
                  start.ru_utime.tv_sec+start.ru_utime.tv_usec*1e-6;
	return diff;
}

double *loadvals(const char *filename,int *dim,int *cols)

/* Attempts to load A from filename. */

{

  FILE *Af;
  double *vals;

  Af = fopen(filename,"r");

  if (Af == NULL) { 

    fprintf(stderr,"tsnnls_test: Couldn't open file %s.\n", filename);
    exit(1);

  }

  if (read_mat(Af, &vals, dim, cols) == -1) { /* If reading as dense fails. */

    if (read_sparse(Af, &vals, dim, cols)) { /* If read as sparse fails. */

      fprintf(stderr,"tsnnls_test: Couldn't parse file %s.\n", filename);
      exit(1);

    }

  }

  fclose(Af);

  return vals;

}

int
read_sparse( FILE *fp, double **vals, int *dim, int *cols )
{
  int i,j;
  double val;
  int fcode;
  int nnz;
  int checknnz;

  if (fscanf(fp, "SPARSE %d %d\n", dim, cols) != 2) { return -1; }
  if (fscanf(fp, "%d\n", &nnz) != 1) { return -1; }

  (*vals) = (double *)(calloc((*dim)*(*cols),sizeof(double)));

  if (*vals != NULL) { /* We can allocate this in one go. Do it the easy way. */

    for(checknnz=0;(fcode = fscanf(fp,"%d %d %lg",&i,&j,&val)) == 3;checknnz++) {

      i--; j--; /* Adjust for 1-based indexing */
      if (i < 0 || i > *dim-1 || j < 0 || j > *cols-1) { return -1; }
      (*vals)[i*(*cols) + j] = val;

    }

    /* We got here because the fscanf didn't work. */

    if (fcode == EOF) { 

      printf("tsnnls_test: Read matrix with %d (claimed %d) nonzeros.\n",
	     checknnz,nnz);
      return 0; 

    } else { return -1; }

  } else {

    printf("tsnnls_test: Couldn't allocate construction buffer for %d x %d matrix.\n",
	   *dim,*cols);
    exit(-1);

  }

}

/* This function reads a dense matrix from Octave as described above */

int
read_mat( FILE* fp, double **vals, int *dim, int *cols)

/* Returns -1 on failure, 0 on success. */

{
  *dim = -1;
  *cols=-1;
  int r, c;
  int startc;
  char linein[4096];
  
  /* First, we scan for lines starting with # */

  while ((startc = fgetc(fp)) != EOF) {

    if (startc == '#') { /* We've started a comment line */

      fgets(linein,sizeof(linein)-1,fp);  
      sscanf(linein," rows: %d",dim);
      sscanf(linein," columns: %d",cols);

    } else if (!isspace(startc)) { /* We've run over into actual data. */

      ungetc(startc,fp); break; 

    }

  }

  if (*dim == -1 || *cols == -1) {

    if (VERBOSITY >= 5) {

      printf("tsnnls_test: File format error. An Octave text file must contain"
	     "             two lines in the form\n"
	     "\n"
	     "             # rows: <m> \n"
	     "             # columns: <n> \n"
	     "\n"
	     "             in an initial block of lines starting with '#'.\n"
	     "\n"
	     "             This file parsed rows %d and cols %d.\n",
	   *dim,*cols);
      
    }

    return -1;

  }

  /* We have parsed the topmatter, and read the size of the matrix. */

  (*vals) = (double*)malloc(sizeof(double)*(*dim)*(*cols));
 
  for( r = 0; r<(*dim); r++ ) {
    for( c=0; c< (*cols); c++ ) {
      
      if (fscanf( fp, "%lf ", &((*vals)[r*(*cols)+c]) ) != 1) { return -1; }

    }
  }

  return 0;
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

static void
tsnnls_test(taucs_ccs_matrix *A,taucs_double *realx,taucs_double *b)
{
  int xItr;
  int pItr = 0;			
  double err = 0;
  double expectedError=0;

  fprintf(stderr,"matrixdims runtime  relerror (|x'-x|/|x|) worst expected\n");
  fprintf(stderr,"--------------------------------------------------------\n");
  
  double residual;
  double ttime;
  double *x;
	
  start_clock();
	
  /* Passing a 0.0 here means we will ALWAYS perform the error-reducing final 
   * step with LSQR. 
   */

  if (solver == pjv) {

    x = t_snnls_pjv(A, b, &residual, 0.0, 0 );

  } else if (solver == fallback) {

    x = t_snnls_fallback(A, b, &residual, 0.0, 0 ) ;

  } else if (solver == tsnnls) {

    x = t_snnls(A, b, &residual, 0.0, 0);

  } else if (solver == spiv) {

    x = t_snnls_spiv(A,b,&residual,0.0,0,nconstrained);

  } else {

    printf("tsnnls_test: Illegal solver in tsnnls_test.\n");
    exit(1);

  }
	
  ttime = end_clock();
	
  printf( "%3d x %3d  %6f ", A->m, A->n, ttime );
  if( x == NULL )
    {
      printf( "Problem %d failed! (Matrix probably not positive definite)\n", pItr );
      taucs_ccs_free(A);
      free(b);
      free(realx);
      
      exit(1);
    }
  
  // compute relative error using ||(x*-x)||/||x||
  
  double* diff = malloc(sizeof(double)*A->n);
  
  for( xItr=0; xItr<A->n; xItr++ )
    diff[xItr] = realx[xItr] - x[xItr];
  
  int incX = {1};
  //err = cblas_dnrm2( A->n, diff, 1 );
  //err /= cblas_dnrm2( A->n, realx, 1 );
  
  err = DNRM2_F77(&(A->n),diff,&incX);
  err /= DNRM2_F77(&(A->n),realx,&incX);
  
  /* Expected relative error given from theory is approx. cond^2*eps */
  expectedError = taucs_rcond(A);
  expectedError = 1.0/expectedError;
  
  expectedError *= expectedError;
  expectedError *= __DBL_EPSILON__;
  
  printf( "%8e          %8e\n", err, expectedError );
  
  if( err > expectedError ) {
    fprintf( stderr, "\t*** WARNING: relative error larger "
	     "than expected, your build may not function as expected\n" );
    exit(1);
  }

  free(diff);
  free(x);
  free(b);
  free(realx);
  taucs_ccs_free(A);

  exit(0);

}

static void
tlsqr_test(taucs_ccs_matrix *A, taucs_double *realx, taucs_double *b)
{
  int xItr;
  double* x;
  double ttime;
	
  start_clock();

  if (solver == tlsqr) { 

    x = t_lsqr(A, b );

  } else if (solver == SOLlsqr) {

    x = lsqrwrapper(A, b);

  } else {

    printf("tsnnls_test: Illegal solver in tlsqr_test.\n");
    exit(1);

  }

  ttime = end_clock();
	
  printf( "tsnnls_test: %d x %d matrix tlsqr runtime %f ", A->m, A->n, ttime );
  
  if( x == NULL ) {

    printf( "tsnnls_test: Solver failed! (A probably not positive definite)\n" );
    taucs_ccs_free(A);
    free(b);
    free(realx);
    
    exit(1);
    
  }
  
  // compute relative error using ||(x*-x)||/||x||
  
  double* diff = malloc(sizeof(double)*A->n);
  double expectedError;
  
  for( xItr=0; xItr<A->n; xItr++ )
    diff[xItr] = realx[xItr] - x[xItr];
  
  double err = 0;
  int incX = {1};
  
  //err = cblas_dnrm2( A->n, diff, 1 );
  //err /= cblas_dnrm2( A->n, realx, 1 );
  
  err = DNRM2_F77(&(A->n),diff,&incX);
  err /= DNRM2_F77(&(A->n),realx,&incX);
  
  /* Expected relative error given from theory is approx. cond^2*eps */
  expectedError = ((double)1.0/taucs_rcond(A));
  expectedError *= expectedError;
  expectedError *= __DBL_EPSILON__;
  
  printf( "relative error ||x'-x||/||x||: %e / "
	  "worst expected error: %e\n", err, expectedError );
  
  if( err > expectedError ) {

    fprintf( stderr, "\t*** WARNING: relative error larger "
	     "than expected, your build may not function as expected\n" );
    exit(1);
  }
  
  free(diff);  
  free(x);
  taucs_ccs_free(A);
  free(b);
  free(realx);

  exit(0);

}

/*
 * This function provides an example of using lsqr for the sparse matrix structure
 * WITHOUT going through our solver. 
 */
extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );

static double*
lsqrwrapper( taucs_ccs_matrix* Af, double* b )
{
  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  double		*xf_raw;
    
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, Af->m, Af->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = Af->m;
  lsqr_in->num_cols = Af->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = 0;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1e13;
  lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 1000;
  lsqr_in->lsqr_fp_out = NULL;	
  for( bItr=0; bItr<Af->m; bItr++ )
    {
      lsqr_in->rhs_vec->elements[bItr] = b[bItr];
    }
  /* Here we set the initial solution vector guess, which is 
   * a simple 1-vector. You might want to adjust this value for fine-tuning
   * t_snnls() for your application
   */
  for( bItr=0; bItr<Af->n; bItr++ )
    {
      lsqr_in->sol_vec->elements[bItr] = 1; 
    }
  
  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = sparse_lsqr_mult;
  
  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Af );
  
  /* copy into xf_raw */
  xf_raw = (double*)malloc(sizeof(double)*Af->n);
  for( bItr=0; bItr<Af->n; bItr++ ) // not really bItr here, but hey
    xf_raw[bItr] = lsqr_out->sol_vec->elements[bItr];
  
  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
  
  return xf_raw;
}

int QUIET = 0;

int main( int argc, char* argv[] ) 
{
  
  FILE *outfile;
  double *Avals = NULL, *bvals = NULL, *xvals = NULL;
  int Adim,Acols,bdim,bcols,xdim,xcols;
  taucs_ccs_matrix *A;
  int tverbosity = 0;

  const char *aname = NULL;
  const char *bname = NULL;
  const char *xname = NULL;

#ifdef WITH_ARGTABLE2

  struct arg_file *arg_Afile = arg_file1("A",NULL,"<A file>","matrix for problem in .mat or .sparse format");
  struct arg_file *arg_bfile = arg_file1("b",NULL,"<b file>","rhs vector b in .mat format");
  struct arg_file *arg_xfile = arg_file0("x",NULL,"<x file>","solution vector x in .mat format");
  struct arg_lit  *arg_tsnnls = arg_lit0(NULL,"tsnnls","solve with tsnnls");
  struct arg_lit  *arg_pjv = arg_lit0(NULL,"pjv","solve with reference solver");
  struct arg_lit  *arg_tlsqr = arg_lit0(NULL,"tlsqr","solve with tlsqr");
  struct arg_lit  *arg_lsqr = arg_lit0(NULL,"lsqr","solve with SOL lsqr");
  struct arg_lit  *arg_fallback = arg_lit0(NULL,"fallback","solve with fallback");
  struct arg_lit  *arg_spiv = arg_lit0(NULL,"spiv","solve with single pivoting solver");
  struct arg_int  *arg_constrained = arg_int0(NULL,"nc","<0-m>","number of constrained vars (spiv only)");

  struct arg_int  *arg_verb = arg_int0("v","Verbosity","<0-10>","verbosity for tsnnls solver");

  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_Afile,arg_bfile,arg_xfile,arg_tsnnls,
		      arg_pjv,arg_fallback,arg_tlsqr,arg_lsqr,arg_spiv,
		      arg_constrained,arg_verb,arg_help,end};
  
  int nerrors;

  /* Now we parse the arguments */

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (arg_help->count > 0) {

    printf("tsnnls_test solves a least-squares or constrained least-squares\n"
	   "problem in the A x = b. The files are expected to be in a sparse\n"
	   "format or in the Octave text format.\n"
	   "\n");

    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    printf("If a solution x is provided, tsnnls_test compares the results of the\n"
	   "calculation with the given solution, provides timing information and\n"
	   "returns 0 if the results match the given solution and 1 otherwise.\n"
	   "\n"
	   "If no solution is provided, tsnnls_test writes the solution to x.mat.\n");

    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  /* Now we read the arguments and move into local vars. */

  aname = arg_Afile->filename[0];
  bname = arg_bfile->filename[0];
  
  if (arg_xfile->count > 0) { xname = arg_xfile->filename[0]; }
  if (arg_lsqr->count > 0)  { solver = SOLlsqr; }
  if (arg_tsnnls->count > 0) { solver = tsnnls; }
  if (arg_pjv->count > 0)    { solver = pjv; }
  if (arg_tlsqr->count > 0)  { solver  = tlsqr; }
  if (arg_fallback->count > 0) {solver = fallback; }
  if (arg_spiv->count > 0) { solver = spiv; }
  if (arg_verb->count > 0) { tverbosity = arg_verb->ival[0]; }

  if (arg_constrained->count > 0) {nconstrained = arg_constrained->ival[0];}
  
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  
#else

  if (argc < 4 || (strcmp(argv[argc-1],"--tsnnls") && strcmp(argv[argc-1],"--tlsqr") && \
		   strcmp(argv[argc-1],"--spiv") && strcmp(argv[argc-1],"--pjv") && \
		   strcmp(argv[argc-1],"--fallback")) 
      || argc > 6) {
		     
    printf("Usage: tsnnls_test <A file> <b file> "
	   "<x file (optional)> <--fallback|--tsnnls|--tlsqr|--pjv|--spiv>\n");
    exit(1);

  }

  if (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-?")) {
    
    printf("tsnnls_test solves a least-squares or constrained least-squares\n"
	   "problem in the A x = b. The files are expected to be in a sparse\n"
	   "format or in the Octave text format.\n"
	   "\n");
    
    printf("Usage: tsnnls_test <A file> <b file> <x file (optional)> "
	   "<--tsnnls|--tlsqr|--pjv|--spiv>\n\n");
    
    printf("If a solution x is provided, tsnnls_test compares the results of the\n"
	   "calculation with the given solution, provides timing information and\n"
	   "returns 0 if the results match the given solution and 1 otherwise.\n"
	   "\n"
	   "If no solution is provided, tsnnls_test writes the solution to x.mat.\n");
	            
    exit(1);

  }
  
  /* Now we attempt to parse the arguments and load files. */

  aname = argv[1];
  bname = argv[2];

  if (argc == 5) { xname = argv[3]; }

  if (!strcmp(argv[argc-1],"--tsnnls")) { solver = tsnnls; }
  else if (!strcmp(argv[argc-1],"--pjv")) { solver = pjv; }
  else if (!strcmp(argv[argc-1],"--tlsqr")) { solver = tlsqr; }
  else if (!strcmp(argv[argc-1],"--fallback")) { solver = fallback; }
  else if (!strcmp(argv[argc-1],"--lsqr")) { solver = SOLlsqr; }
  else if (!strcmp(argv[argc-1],"--spiv")) { solver = spiv; }
  else {
    
    printf("tsnnls_test: Unknown solver %s.\n",argv[argc-1]);
    exit(1);
    
  }

#endif

  printf("tsnnls_test %s (%s %s).\n",PACKAGE_VERSION, __DATE__ , __TIME__ );
  tsnnls_version(NULL,0);
  tsnnls_verbosity(tverbosity);

  Avals = loadvals(aname, &Adim, &Acols);
  A = taucs_construct_sorted_ccs_matrix(Avals, Acols, Adim);
  free(Avals);

  printf("tsnnls_test: Loaded %d x %d matrix A from %s.\n",
	 Adim,Acols,aname);

  if (nconstrained == -1) { nconstrained = Acols; }

  /* Now load b. */

  bvals = loadvals(bname, &bdim, &bcols);
  
  if (bdim != Adim || bcols != 1) {

    printf("tsnnls_test: We expect the second (b) file %s to be %d x 1.\n"
	   "             The given file is %d x %d.\n",
	   bname,Adim,bdim,bcols);
    exit(1);

  }

  printf("tsnnls_test: Loaded %d x 1 rhs vector b from %s.\n",
	 bdim,bname);

  /* Now, if present, load x */

  if (xname != NULL) {

    xvals = loadvals(xname, &xdim, &xcols);

    if (xdim != Acols || xcols != 1) {

      printf("tsnnls_test: We expect the third (x) file %s to be %d x 1.\n"
	     "             The given file is %d x %d.\n",
	     xname,Acols,xdim,xcols);
      exit(1);

    }

    printf("tsnnls_test: Loaded %d x %d solution vector x from %s.\n",
	   xdim,xcols,xname);

    if (solver == tsnnls || solver == pjv || solver == fallback || solver == spiv) {
      
      tsnnls_test(A,xvals,bvals);
      
    } else if (solver == tlsqr || solver == SOLlsqr) {
      
      tlsqr_test(A,xvals,bvals);

    } 

  } 

  /* If we've survived this long, we're in problem mode instead of test mode. */

  double residual;

  if (solver == fallback) {

    xvals = t_snnls_fallback(A, bvals, &residual, -1, 1);
    
  } else if (solver == pjv) {

    xvals = t_snnls_pjv(A, bvals, &residual, -1, 1);

  } else if (solver == spiv) {

    xvals = t_snnls_spiv(A, bvals, &residual, -1, 1,nconstrained);

  } else if (solver == tlsqr) {
   
    xvals = t_lsqr(A, bvals);

  } else if (solver == SOLlsqr) {

    xvals = lsqrwrapper(A, bvals);

  } else {

    xvals = t_snnls(A, bvals, &residual, -1, 1);

  }

  if (xvals == NULL) {

    printf("tsnnls_test: solver failed.\n");
    exit(1);

  }

  printf("tsnnls_test: Solved problem.\n");

  outfile = fopen("x.mat","w");
  
  if (outfile == NULL) {

    printf("tsnnls_test: Could not open x.mat to write solution.\n");
    exit(1);

  }

  colvector_write_mat(outfile,xvals,Acols,"x");
  
  exit(0);

}

int old_readsparse(FILE *fp,double **x,double **b,taucs_ccs_matrix **outMatrix)

{
  int dim = 0, cols=0, r, nnz, colOffset, cItr;
  double* vals = NULL;
  taucs_ccs_matrix* result;
  int vItr=0;
  int*	colcounts;
  double** colels;
  int**	colrows;
  int*	colitrs;
  fpos_t  elementsStart;
  int elr, elc;
  double theVal;
  
  fgetpos(fp, &elementsStart);

  if (fscanf(fp, "SPARSE %d %d\n", &dim, &cols) != 2) { return -1; }
  if (fscanf(fp, "%d\n", &nnz) != 1) { return -1; }
  
  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  result->n = cols;
  result->m = dim;
  result->flags = TAUCS_DOUBLE;
  
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*nnz);
  result->values.d = (double*)malloc(sizeof(double)*nnz);  
   
  vals = (double*)calloc(nnz,sizeof(double));
  colels = (double**)calloc(result->n,sizeof(double*));
  colrows = (int**)calloc(result->n,sizeof(int*));
  colcounts = (int*)calloc(result->n,sizeof(int));
  
  /* we have to munger our file reading a bit since we're given
   * output in row1 col1, col2, ... row2 col1, col2, but it's better
   * than allocating a m*n double array when m and n are large
   */
  for( r=0; r<nnz; r++ )
    {
      /* first get column counts */
      fscanf(fp, "%d %d %lf\n", &elr, &elc, &theVal);
      elr--;
      elc--; // adjust for 0 based indexing
      colcounts[elc]++;
    }
  
  /* Now that we have the count of each column, we can allocate colels to 
   * store the values on our second pass
   */
  for( cItr=0; cItr<result->n; cItr++ )
    {
      colels[cItr] = (double*)calloc(colcounts[cItr], sizeof(double));
      colrows[cItr] = (int*)calloc(colcounts[cItr], sizeof(int));
    }
  
  fsetpos(fp, &elementsStart);
  
  colitrs = (int*)calloc(result->n, sizeof(int));
  for( r=0; r<nnz; r++ )
    {
      fscanf(fp, "%d %d %lf\n", &elr, &elc, &theVal);
      elr--;
      elc--; // adjust for 0 based indexing
      colels[elc][colitrs[elc]] = theVal;
      colrows[elc][colitrs[elc]] = elr;
      colitrs[elc]++;
    }
  free(colitrs);
  
  /* now that we have the elements ordered in their columns and we've also 
   * kept track of the row indices, we are ready to construct the ccs structure 
   * which was (comparatively easy!) to allocate
   */
  colOffset = 0;
  for( cItr=0; cItr<result->n; cItr++ )
    {
      result->colptr[cItr] = colOffset;
      for( vItr=0; vItr<colcounts[cItr]; vItr++ )
	{
	  result->values.d[colOffset] = colels[cItr][vItr];
	  result->rowind[colOffset] = colrows[cItr][vItr];
	  colOffset++;
	}
    }
  result->colptr[result->n] = nnz;
  
  free(vals);
  free(colcounts);
  for( cItr=0; cItr<result->n; cItr++ )
    free(colels[cItr]);
  for( cItr=0; cItr<result->n; cItr++ )
    free(colrows[cItr]);
  free(colels);
  free(colrows);
  
  *x = (double*)malloc(sizeof(double)*cols);
  for( r=0; r<cols; r++ )
    {
      double a;
      fscanf(fp, "%lf\n", &a);
      (*x)[r] = a;
    }
  
  *b = (double*)malloc(sizeof(double)*dim);
  for( r=0; r<dim; r++ )
    fscanf(fp, "%lf\n", &(*b)[r]);
  
  *outMatrix = result;

  return 0;
}
