

/*
 * This program is free software distributed under the GPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/gpl.txt
 */
 
/*
 * This program performs several tests on the tsnnls library,
 * providing relative error and timing tests for both constrained and
 * unconstrained solvers with the following matrix input file formats:
 *
 * Note that for both formats, row and column indexing is 1-based.
 *
 * Sparse format is a plain text file containing:
 * row_number column_number
 * number_of_nonzeros
 * r1 c1 val11
 * r1 c2 val12
 * ...
 * r2 c1 val21
 * ...
 *
 * Dense format is a plain text file containing:
 * row_number
 * col_number
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

#include "lsqr.h"
#include "tsnnls.h"

#ifdef HAVE_DARWIN          /* We use the Apple BLAS/LAPACK if possible. */
 #include <vecLib/vBLAS.h>
#else
 #include "tsnnls/gsl_cblas.h"
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

struct rusage start, end;
int tpStorage;

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

static void
read_sparsetp( FILE* fp, taucs_ccs_matrix** outMatrix, double** b, double** x )
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
			
	fscanf(fp, "%d %d\n", &dim, &cols);
	fscanf(fp, "%d\n", &nnz);
	
	fgetpos(fp, &elementsStart);
	
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
}

/* This function reads a dense test problem as described above */
static void
read_tp( FILE* fp, taucs_ccs_matrix** outMatrix, double** b, double** x )
{
	int dim = 0, cols=0, r, c;
	double* vals = NULL;
			
	fscanf(fp, "%d %d\n", &dim, &cols);

	vals = (double*)malloc(sizeof(double)*dim*cols);
	for( r = 0; r<dim; r++ )
	{
		for( c=0; c<cols; c++ )
			fscanf( fp, "%lf ", &vals[r*cols+c] );
		fscanf(fp, "\n");
	}
		
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
	
	*outMatrix = taucs_construct_sorted_ccs_matrix(vals, cols, dim);

	free(vals);
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

static void
taucs_snnls_test(FILE* fp)
{
	int xItr;
	double relErr = 0.0;
	
	taucs_ccs_matrix *Accs;
	double* b;
	double* realx;

	int ran = 0;
	int pItr = 0;			
	double err = 0;
	double expectedError=0;

	
	{
		while( !feof(fp) )
		{
			double* x;
			double residual;
			double ttime;
			
			if( tpStorage == 0 )
				read_tp( fp, &Accs, &b, &realx );
			else
				read_sparsetp( fp, &Accs, &b, &realx );
				
			start_clock();
			
			/* Passing a 0.0 here means we will ALWAYS perform the error-reducing final 
			 * step with LSQR. 
			 */
			x = t_snnls(Accs, b, &residual, 0.0, 0 );
			
			ttime = end_clock();
			
			printf( "%d: Matrix of size: (%d %d) runtime: %f ", ran++, Accs->m, Accs->n, ttime );
			if( x == NULL )
			{
				printf( "Problem %d failed! (Matrix probably not positive definite)\n", pItr );
				taucs_ccs_free(Accs);
				free(b);
				free(realx);
				return;
			}
			
			// compute relative error using ||(x*-x)||/||x||
			{
				double* diff = malloc(sizeof(double)*Accs->n);
				
				for( xItr=0; xItr<Accs->n; xItr++ )
					diff[xItr] = realx[xItr] - x[xItr];
				
				err = cblas_dnrm2( Accs->n, diff, 1 );
				err /= cblas_dnrm2( Accs->n, realx, 1 );
				
				/* Expected relative error given from theory is approx. cond^2*eps */
				expectedError = taucs_rcond(Accs);
				expectedError = 1.0/expectedError;
			
				expectedError *= expectedError;
				expectedError *= __DBL_EPSILON__;
				
				printf( "relative error ||x'-x||/||x||: %e / worst expected error: %e\n", err, expectedError );
				
				if( err > expectedError )
					fprintf( stderr, "\t*** WARNING: relative error larger than expected, your build may not function as expected\n" );
				
				if( err > relErr )
					relErr = err;
				
				free(diff);
			}
			
			free(x);
				
			taucs_ccs_free(Accs);
			free(b);
			free(realx);
		}
	}	
	printf( "max relative error: %e\n", relErr );
}

static void
taucs_tlsqr_test(FILE* fp)
{
	int xItr;
	double relErr = 0.0;
	
	taucs_ccs_matrix *Accs;
	double* b;
	double* realx;

	int pItr;
	
	pItr=0;

	{
		while( !feof(fp) )
		{
			double* x;
			double ttime;
			
			if( tpStorage == 0 )
				read_tp( fp, &Accs, &b, &realx );
			else
				read_sparsetp( fp, &Accs, &b, &realx );
				
			start_clock();
			
			x = t_lsqr(Accs, b );
			
			ttime = end_clock();
			
			printf( "Matrix of size: (%d %d) runtime: %f ", Accs->m, Accs->n, ttime );
			if( x == NULL )
			{
				printf( "Problem failed! (Matrix probably not positive definite)\n" );
				taucs_ccs_free(Accs);
				free(b);
				free(realx);
				return;
			}
			
			// compute relative error using ||(x*-x)||/||x||
			{
				double* diff = malloc(sizeof(double)*Accs->n);
				double expectedError;
				
				for( xItr=0; xItr<Accs->n; xItr++ )
					diff[xItr] = realx[xItr] - x[xItr];
				
				double err = 0;
				
				err = cblas_dnrm2( Accs->n, diff, 1 );
				err /= cblas_dnrm2( Accs->n, realx, 1 );
				
				/* Expected relative error given from theory is approx. cond^2*eps */
				expectedError = ((double)1.0/taucs_rcond(Accs));
				expectedError *= expectedError;
				expectedError *= __DBL_EPSILON__;
				
				printf( "relative error ||x'-x||/||x||: %e / worst expected error: %e\n", err, expectedError );
				
				if( err > expectedError )
					fprintf( stderr, "\t*** WARNING: relative error larger than expected, your build may not function as expected\n" );

				if( err > relErr )
					relErr = err;
				
				free(diff);
			}
			
			free(x);
			taucs_ccs_free(Accs);
			free(b);
			free(realx);
		}
	}	
	printf( "max relative error: %e\n", relErr );
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

int main( int argc, char* argv[] ) 
{

	if( argc != 4 )
	{
		printf( "Need 3 arguments\n\tfirst: a test problem file\n\tsecond: (0/1) indicating: (0): testing snnls (1): testing lsqr\n\tthird: (0/1) indicating: (0): using full format (1): using sparse format\n" );
		exit(-1);
	}
	
	FILE* fp = fopen(argv[1],"r");
	if( fp == NULL )
	{
		fprintf( stderr, "Need a problem test file that exists!\n" );
		exit(-1);
	}
	
	tpStorage = atoi(argv[3]);
	if( tpStorage == 0 )
		printf( "Expecting full format\n" );
	else if( tpStorage == 1  )
		printf( "Expecting sparse format\n" );
	else
	{
		fprintf(stderr, "Unknown matrix format specifier: %d\n", tpStorage);
		return -1;
	}
	
	switch( atoi(argv[2]) )
	{
		case 0:
			printf( "taucs snnls test\n" );
			taucs_snnls_test(fp);
			break;
		
		case 1:
			printf( "taucs lsqr test\n" );
			taucs_tlsqr_test(fp);
			break;
	}
	fclose(fp);

	return 0;
}

