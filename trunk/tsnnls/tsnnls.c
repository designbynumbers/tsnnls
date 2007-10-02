 
/*
 * This program is free software distributed under the GPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/gpl.txt
 */

#include <config.h>

#include "tsnnls_blas_wrappers.h"

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif 

#include<assert.h>

//#ifdef HAVE_DARWIN  /* We expect to use the Accelerate framework */
//  #include <vecLib/vBLAS.h>
//  #include <vecLib/clapack.h>
//#else
//  #include "gsl_cblas.h"
//#endif 

//#ifdef HAVE_CBLAS_H
//  #include <cblas.h>
//#else
//  #ifdef HAVE_ATLAS_CBLAS_H
//     #include <atlas/cblas.h>
//  #else
//     #ifdef HAVE_VECLIB_CBLAS_H
//       #include <vecLib/cblas.h>
//     #endif
//  #endif
//#endif

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

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#include "acint32_type.h"
#include "lsqr.h"
#include "tsnnls.h"

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

int gVERBOSITY;

/*

    File contains an implementation of the block-pivoting non-negative least squares
    algorithm of Portugal, Judice, and Vicente. The reference for the algorithm is their
    paper:

    @article {MR95a:90059,
    AUTHOR = {Portugal, Lu{\'{\i}}s F. and J{\'u}dice, Joaqu{\'{\i}}m J. and
              Vicente, Lu{\'{\i}}s N.},
     TITLE = {A comparison of block pivoting and interior-point algorithms
              for linear least squares problems with nonnegative variables},
   JOURNAL = {Math. Comp.},
  FJOURNAL = {Mathematics of Computation},
    VOLUME = {63},
      YEAR = {1994},
    NUMBER = {208},
     PAGES = {625--643},
      ISSN = {0025-5718},
     CODEN = {MCMPAF},
   MRCLASS = {90C20 (65K05)},
  MRNUMBER = {95a:90059},
}

   The paper is stored on JSTOR.


   Cantarella/Piatek. 5/04.

*/

void tsnnls_version( char *version, size_t strlen) {
  if (version == NULL) {
    printf("tsnnls Version: %s\n",PACKAGE_VERSION);
  } else {
    (void)snprintf(version,strlen,PACKAGE_VERSION);
  }
}

void tsnnls_verbosity(int level) {

  gVERBOSITY = level;

  if (gVERBOSITY >= 10) {

    printf("tsnnls verbosity level changed to %d.\n",level);

  }

}

void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );


/* The entries of setA are indices into the array varsA. Procedure finds the 
 * negative entries of varsA which correspond to indices in A, and returns a 
 * list of them in <infeas> with size <sizeInfeas>. 
 * 
 * We assume that infeas is allocated, and is large enough to hold the result. 
 */
void infeasible( const int* setA, const double* varsA, const int sizeA, 
		 int* infeas, int* sizeInfeas )
{
  int i, pItr=0;
  for( i=0; i< sizeA; i++ ) 
    {
      if( varsA[setA[i]] < 0 ) 
	{
	  infeas[pItr++] = setA[i];
	}
    }
  *sizeInfeas = pItr;
}


/* Removes the elements of setB from setA, assuming 
 * that both sets are sorted ascending. 
 */
void 
int_difference( int* setA, int aSize, const int* setB, const int bSize, int* outSize )
{
	int aRead = 0, aWrite = 0, bPtr = 0;

	/* Some special cases up front */
	if( aSize == 0 && bSize == 0 )
	{
		*outSize = 0;
		return;
	}

	if( bSize == 0 )
	{
		*outSize = aSize;
		return;
	}

	if( aSize == 0 && bSize != 0 )
	{
		fprintf( stderr, "int_differe(): aSize is 0 but bSize is not!\n" );
		exit(-1);
	}

	*outSize = 0;

	for( aRead=aWrite=0; aRead<aSize; aRead++ ) 
	{
		/* We should update bPtr, if we can. */
		if (setA[aRead] > setB[bPtr] && bPtr < bSize-1) 
		{ 
			bPtr++;
		}

		/* We don't have a match. */
		if (setA[aRead] != setB[bPtr]) 
		{ 
			setA[aWrite++] = setA[aRead];
			(*outSize)++;
		}
	}
}


/*
 * performs union, assumes both setA, setB are !sorted in ascending order!.
 * writes the result to set A.
 */
void 
int_union( int* setA, int aSize, const int* setB, const int bSize, int* outSize )
{
	int aItr, bItr, uItr;
	aItr = bItr = uItr = 0;
	int i;
	int*	unionSet = NULL;

	*outSize = 0;

	if( aSize == 0 && bSize == 0 )
		return;

	unionSet = (int*)malloc(sizeof(int)*(aSize+bSize));
	
	if( aSize == 0 )
	{
		memcpy( unionSet, setB, sizeof(int)*bSize);
		*outSize = bSize;
	}
	else if( bSize == 0 )
	{
		memcpy( unionSet, setA, sizeof(int)*aSize);
		*outSize = aSize;
	}
	else
	{
		while( 1==1 )
		{
			if( setA[aItr] == setB[bItr] )
			{
				unionSet[uItr++] = setA[aItr];
				aItr++;
				bItr++;
			}
			else
			{
				if(setA[aItr] < setB[bItr] )
					unionSet[uItr++] = setA[aItr++];
				else
				{
					unionSet[uItr++] = setB[bItr++];
				}
			}

			if( aItr == aSize )
			{
				/* copy remaining b elements */
				for( ; bItr<bSize; bItr++ )
					unionSet[uItr++] = setB[bItr];
				break;
			}
			if( bItr == bSize )
			{
				/* copy remaining a elements */
				for( ; aItr<aSize; aItr++ )
					unionSet[uItr++] = setA[aItr];
				break;
			}
		} /* while */

		*outSize = uItr;
	} // else 
	
	/* Now we recopy the results to A. */
	for(i=0;i<*outSize;i++) 
	{
		setA[i] = unionSet[i];
	}

	free(unionSet);
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

static void
fix_zeros( double* values, int size, double rcond, int inPrintErrorWarnings )
{
	int i;
	double eps = __DBL_EPSILON__; // as defined in float.h
	for( i=0; i<size; i++ )
	{
		if( inPrintErrorWarnings != 0 )
		{
			double cond = (double)1/rcond;
						
			if( fabs(values[i]) < (cond*cond*eps) )
			{
				fprintf( stderr, "Variable is within (condition number)^2*eps of 0, accuracy results may be unexpected!\n" );
				inPrintErrorWarnings = 0; // It should suffice to notify the user once. Otherwise this will happen MANY times
			}
		}
		
		if( fabs(values[i]) < kZeroThreshold )
		{
			values[i] = 0.0;
		}
	}
}

taucs_ccs_matrix*
selectAprimeDotA( double* apda, int cols, int* F, int sizeF )
{
	taucs_ccs_matrix* result = NULL;
	int maxSize, i, j, valItr;
	
	result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

	result->n = sizeF;
	result->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
	
	maxSize = (result->n*(result->n+1))/2;

	result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
	result->rowind = (int*)malloc(sizeof(int)*maxSize);
	result->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
	
	valItr = 0;
	for( i=0; i<sizeF; i++ )
	{
		result->colptr[i] = valItr;
		for( j=i; j<sizeF; j++ )
		{
			double val = apda[cols*F[j] + F[i]];
			if( val != 0 )
			{
				result->values.d[valItr] = val;
				result->rowind[valItr] = j;
				valItr++;
			}
		}
	}
	result->colptr[i] = valItr;

	return result;
}

void
selectAprimeDotAsparse( const taucs_ccs_matrix* apda, int* F, int sizeF, 
			taucs_ccs_matrix* inOutApda )
{
  /* This routine presume that inOutApda has already been
   * allocated (this save allocation and free time) the size
   * should be the same as apda, since this will be a submatrix,
   * it cannot possibly be larger */
  
  int fItr, valItr, cItr, rItr;
  
  if( sizeF == 0 )
    {
      inOutApda->m = inOutApda->n = 0;
      return;
    }
  
  inOutApda->n = sizeF;
  inOutApda->m = sizeF;

  inOutApda->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
  
  /* start in column F[0] and select entries that are in row F[0], F[1], F[2], ... etc */
  valItr = cItr = 0;
  for( cItr=0; cItr<sizeF; cItr++ )
    {
      inOutApda->colptr[cItr] = valItr;
      rItr=apda->colptr[F[cItr]];
      /* scan down rows until we can compare to F value */
      fItr=0;
      
      /* Since this matrix is symmetric storing the lower triangle, we need to adjust fItr to start at the main diagonal
       */
      //	while( F[fItr] < F[cItr] )
      //		fItr++;
      fItr=cItr;
      
      while( 1==1 )
	{
	  /* have we exhausted F? */
	  if( fItr == sizeF )
	    break;
	  
	  /* Continues until we can compare to F's selection or we finish this column */
	  while( apda->rowind[rItr] < F[fItr] && rItr<apda->colptr[F[cItr]+1] )
	    rItr++;
	  
	  /* Are we dont with this column? */
	  if( rItr == apda->colptr[F[cItr]+1] )
	    break;
	  
	  /* We have been been skipping along the sparse entries of
	   * apda, which may or may not correspond to row requests
	   * from F. Once we're out of the above loop, we are either
	   * AT F[fItr] in which case we want this entry, or we have
	   * overshot it. If we've overshot it, we need to increase
	   * fItr.
	   */
	  if( apda->rowind[rItr] == F[fItr] )
	    {
	      inOutApda->values.d[valItr] = apda->values.d[rItr];
	      inOutApda->rowind[valItr] = fItr; /* fItr here is the row in our submatrix (since it is sizeF by sizeF) */
	      valItr++;
	      rItr++;
	      fItr++; /* we got this one */
	    }
	  else
	    {
	      fItr++;
	    }
	} // while row scanning
    } // for F column selections
  
  /* fix up the last column start, valItr is one more than the index of the last guy
   * entered, which is what we want */
  inOutApda->colptr[sizeF] = valItr;
}

taucs_ccs_matrix*
selectAprime( double* Ap, int cols, int rows, int* F, int sizeF )
{
	taucs_ccs_matrix* result = NULL;
	int maxSize, i, j, valItr;
	
	result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

	result->m = sizeF;
	result->n = cols;
	result->flags = TAUCS_DOUBLE;
	
	maxSize = result->m*result->n;

	result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
	result->rowind = (int*)malloc(sizeof(int)*maxSize);
	result->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
	
	valItr = 0;
	for( i=0; i<result->n; i++ )
	{
		result->colptr[i] = valItr;
		for( j=0; j<result->m; j++ )
		{
			double val = Ap[cols*F[j] + i];
			if( val != 0 )
			{
				result->values.d[valItr] = val;
				result->rowind[valItr] = j;
				valItr++;
			}
		}
	}
	result->colptr[i] = valItr;

	return result;
}

/* the nature of this routine is described in lsqr.h ~line 316, primarily:
 *
 *									If MODE = 0, compute  y = y + A*x,
 *									If MODE = 1, compute  x = x + A^T*y.
 */
void 
sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod )
{
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)prod;
	int cItr, rItr;
	
	if( mode == 0 )
	{
		int i,ip,j,n, rows;
		taucs_datatype Aij;

		n = A->n;
		rows = A->m;
	
	/*	double* Ax = (double*)malloc(sizeof(double)*A->m);
		ourtaucs_ccs_times_vec( A, x->elements, Ax );
		for(i=0; i<A->m; i++ )
			y->elements[i] += Ax[i];
		free(Ax);
	*/
	//	ourtaucs_ccs_times_vec( A, x->elements, y->elements );
	
		for (j=0; j<n; j++) {
		  for (ip = (A->colptr)[j]; ip < (A->colptr[j+1]); ip++) {
			i   = (A->rowind)[ip];
			Aij = (A->values.d)[ip];
			
			y->elements[i] = y->elements[i]+(x->elements[j]*Aij);
		  }
		}
	
	}
	else if( mode == 1 )
	{
		// since A^T*y = y^T*A
	/*	double* yA = (double*)malloc(sizeof(double)*A->n);
		taucs_transpose_vec_times_matrix_nosub( y->elements, A, yA );
		for(i=0; i<A->n; i++ )
			x->elements[i] += yA[i];
		free(yA);
	*/
	//	taucs_transpose_vec_times_matrix_nosub( y->elements, A, x->elements );
		for( cItr=0; cItr<A->n; cItr++ )
		{
			//result[cItr] = 0;
			for( rItr=0; rItr<A->colptr[cItr+1]-A->colptr[cItr]; rItr++ )
			{
				int tRow = A->rowind[A->colptr[cItr] + rItr];
				x->elements[cItr] += y->elements[tRow] * A->values.d[A->colptr[cItr] + rItr];
			}
		}

	}
	else
		fprintf(stderr, "Unknown mode: %ld\n", mode );
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

taucs_double*
t_snnls_pjv( taucs_ccs_matrix *A_original_ordering, taucs_double *b, 
	 double* outResidualNorm, double inRelErrTolerance, int inPrintErrorWarnings )

/* This code, now deprecated, uses an older algorithm of Portugal,
   Judice and Vicente to compute the solution of the sparse
   nonnegative least-squares problem. The code can be more efficient
   than the updated block3 algorithm used in this version of tsnnls,
   so it is still available. */

{
  taucs_ccs_matrix  *Af;
  int               p, ninf, pbar = {3};
  int               m,n,i, maxSize;
  
  int				A_rows, A_cols;
  
  int              *F, *G, *H1, *H2, SCR[1];
  int              sizeF, sizeG, sizeH1, sizeH2, sizeSCR = {1};
  int				lsqrStep=0;
  double			rcond=1;
  
  /* These variables are subsets of the column indices of the matrix A, 
   * always stored in sorted order, and sized by the corresponding
   * "size" integers. 
   * 
   * Like the row indices themselves, they are 0-based. 
   */
  taucs_double     *x,*y, *xf_raw = NULL, *yg_raw, *residual;
  
  int AprimeDotA_cols;
  
  taucs_ccs_matrix* AprimeDotA = taucs_ccs_aprime_times_a(A_original_ordering);
  taucs_ccs_matrix*   lsqrApA;
  
  /* create a copy of AprimeDotA memory wise to store the tlsqr submatrices */
  {
    lsqrApA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
    lsqrApA->n = AprimeDotA->n;
    lsqrApA->flags = TAUCS_DOUBLE;
    lsqrApA->flags = lsqrApA->flags | TAUCS_SYMMETRIC;
    lsqrApA->flags = lsqrApA->flags | TAUCS_LOWER; // rep the lower half
    lsqrApA->colptr = (int*)malloc(sizeof(int)*(lsqrApA->n+1));
    /* This is the number of nonzeros in A'*A, which we cannot overflow with a submatrix */
    maxSize = AprimeDotA->colptr[AprimeDotA->n]; 
    lsqrApA->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
    lsqrApA->rowind = (int*)malloc(sizeof(int)*maxSize);
  }
  
  if( inRelErrTolerance <= 0 )
    lsqrStep = 1;
  
  A_rows = A_original_ordering->m;
  A_cols = A_original_ordering->n;
  
  AprimeDotA_cols = A_cols;
  
  m = A_original_ordering->m;
  n = A_original_ordering->n;
  
  // This initial values is suggested by PJV
  ninf = n+1;
  
  /* We first allocate space. */
  F   = calloc(n,sizeof(int));
  G   = calloc(n,sizeof(int));
  H1  = calloc(n,sizeof(int));
  H2  = calloc(n,sizeof(int));
  
  x   = calloc(n,sizeof(taucs_double));
  y   = calloc(m,sizeof(taucs_double));
  
  /* submatrix allocation actually takes bit of time during profiling,
   * so we reuse an allocation that cannot be overflowed by
   * submatrices. Note that
   * A_original_ordering->colptr[A_original_ordering->n] is the number
   * of nonzero entries in A
   */
  Af = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  Af->colptr = (int*)malloc(sizeof(int)*(A_cols+1));
  Af->rowind = (int*)malloc(sizeof(int)*(A_original_ordering->colptr[A_original_ordering->n]));
  Af->values.d = (double*)malloc(sizeof(double)*A_original_ordering->colptr[A_original_ordering->n]);
  
  /* Next we initialize variables. */
  for(i=0;i<n;i++) 
    {
      G[i] = i;
    }
  sizeF = 0; sizeG = n;  
  p = pbar;
  
  /* Set y = A'b, which is the same as y=b'A. We perform that computation as it is faster */
  //	transpose_vec_times_matrix(b, A, G, A_original_ordering->n, A_original_ordering->m, A_original_ordering->n, y);
  taucs_transpose_vec_times_matrix(b, A_original_ordering, G, A_cols, y);
  
  int incY = {1};
  double alpha = {-1.0};
  
  //cblas_dscal(n,-1.0,y,1);      /* Scalar multiply y *= -1. */
  DSCAL_F77(&n,&alpha,y,&incY);
  
  /* Now we enter the main loop. */
  infeasible(F,x,sizeF,H1,&sizeH1);  
  infeasible(G,y,sizeG,H2,&sizeH2);
  
  for( ; sizeH1 > 0 || sizeH2 > 0 ;
       infeasible(F,x,sizeF,H1,&sizeH1),  
	 infeasible(G,y,sizeG,H2,&sizeH2)  ) 
    {
      
      
      /* We found infeasible variables. We're going to swap them 
	 between F and G according to one of the schemes below. */
      
      if (sizeH1 + sizeH2 < ninf) 
	{  
	  /* The number of infeasibles _decreased_ since the last try. */
	  /* This is good. We reset the counters and swap all infeasibles. */
	  ninf = sizeH1 + sizeH2;
	  p = pbar;
	  
	  int_union(F,sizeF,H2,sizeH2,&sizeF);
	  int_difference(F,sizeF,H1,sizeH1,&sizeF);
	  
	  int_union(G,sizeG,H1,sizeH1,&sizeG);
	  int_difference(G,sizeG,H2,sizeH2,&sizeG);
	} 
      else 
	{  
	  if( p > 0 ) 
	    {  
	      /* Things didn't go well last time-- we _increased_ the number */
	      /* of infeasibles. But we haven't run out of chances yet, so   */
	      /* we'll decrement the counter, and try the large swap again.  */
	      p--;
	      
	      int_union(F,sizeF,H2,sizeH2,&sizeF);
	      int_difference(F,sizeF,H1,sizeH1,&sizeF);
	      
	      int_union(G,sizeG,H1,sizeH1,&sizeG);
	      int_difference(G,sizeG,H2,sizeH2,&sizeG);
			} /* (p>0) */
	  else 
	    {
	      /* Real trouble. Large swaps aren't reducing the number of     */
	      /* infeasibles, even after pbar tries. So this time, we'll     */
	      /* revert to the slower "Murty's Method", and move only one    */
	      /* guy-- the guy with the highest column index. */
	      
	      if( sizeH1 > 0 && sizeH2 > 0 )
		{
		  if ( H1[sizeH1 - 1] > H2[sizeH2 - 1] ) 
		    {  
		      /* H1 contains the last column index. */

		      SCR[0] = H1[sizeH1-1];
		      
		      int_difference(F,sizeF,SCR,sizeSCR,&sizeF);
		      int_union(G,sizeG,SCR,sizeSCR,&sizeG);
		    } 
		  else 
		    {
		      /* H2 contains the last column index. */
		      
		      SCR[0] = H2[sizeH2-1];
						
		      int_union(F,sizeF,SCR,sizeSCR,&sizeF);
		      int_difference(G,sizeG,SCR,sizeSCR,&sizeG);
		    }
		}
	      else if( sizeH1 == 0 )
		{
		  SCR[0] = H2[sizeH2-1];
		  
		  int_union(F,sizeF,SCR,sizeSCR,&sizeF);
		  int_difference(G,sizeG,SCR,sizeSCR,&sizeG);
		}
	      else
		{
		  SCR[0] = H1[sizeH1-1];
		  
		  int_difference(F,sizeF,SCR,sizeSCR,&sizeF);
		  int_union(G,sizeG,SCR,sizeSCR,&sizeG);
				}		
	      /* p is still 0, and will remain so until we start making */
	      /* progress at reducing the number of infeasibilities. */
	    } /* else (p>0) */
	} /* else (sizeH1 + sizeH2 < ninf) */
      
      /* We have now altered F and G, and are ready to recompute everything. */
      /* We first clear x and y, since we won't need them anymore. */
      for(i=0;i<n;i++) 
	{ 
	  x[i] = y[i] = 0.0; 
	}
      
      /* 
       *
       * By (7) in PJV, x_F should be the solution of the unconstrained
       * least squares problem:
       * 
       * min || A_F x_F - b ||,
       * 
       * where A_F is the submatrix of A containing the columns of A 
       * with column indices in F. 
       * 
       * And by (8) in PJV, the y_F should be (A_G)'*(A_F x_F - b). 
       * 
       * We first solve the lsqr problem.
       * 
       */
      taucs_ccs_submatrix(A_original_ordering, F, sizeF, Af);
      
      if( sizeF != 0 )
	{
	  /* we compute these values based on selections based on F
	   * since it's faster than recalculating them in lsqr. This
	   * requires the use of a custom lsqr solver that expects
	   * extra parameters from snnls, however.
	   */
	  
	  selectAprimeDotAsparse(AprimeDotA, F, sizeF, lsqrApA); 	
	  
	  xf_raw = NULL;
	  if( inRelErrTolerance > 1 || (lsqrStep != 0 && inPrintErrorWarnings == 0) )
	    xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, NULL);		
	  else
	    {
	      xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, &rcond );
	      if( (1/rcond)*(1/rcond)*__DBL_EPSILON__ < inRelErrTolerance )
		lsqrStep = 1;
	    }
	  if( xf_raw == NULL )
	    return NULL; // matrix probably not positive definite
	  
	  /* taucs_snnls requires us to handle in some
	   * way the case where values returned from
	   * lsqr that are within error tolerance of
	   * zero get continually swapped between x and
	   * y because our comparison to zero is
	   * numerically more precise in the infeasibles
	   * computation. In order to terminate, we must
	   * handle this case. For our usage of snnls,
	   * it is most efficient to simply zero values
	   * within some epsilon of zero as returned
	   * from lsqr, but this solution may not be the
	   * best for all possible applications.
	   */
	  fix_zeros(xf_raw, Af->n, rcond, inPrintErrorWarnings);
	  
	  /* Now compute the residual A_F x_F - b. This is an m-vector. */
	  residual = (taucs_double *)calloc(m,sizeof(taucs_double));
	  ourtaucs_ccs_times_vec(Af,xf_raw,residual);
	}
      else
	{	  
	  /* 
	   * if sizeF is 0, the meaning of residual changes (since
	   * there really is _no_ matrix), so we'll just set the
	   * residual to -b, but we still need a zeroed residual to do
	   * the below computation to make that happen, which calloc
	   * does here
	   */
	  residual = (taucs_double *)calloc(m,sizeof(taucs_double));
	}
      
      double alpha = {-1.0};
      int incX = {1}, incY = {1};
      
      //cblas_daxpy(m,-1.0,b, 1, residual, 1);
      DAXPY_F77(&m,&alpha,b,&incX,residual,&incY);
      
      /* We now compute (A_G)'. */
      /* And finally take (A_G)'*residual. This is a sizeG-vector. */
      yg_raw = (taucs_double *)calloc(sizeG,sizeof(taucs_double));     
      
      /* 
       * We now should compute (A_G)', and take (A_G)'*residual, but 
       * A_G'*residual = residual'*A_G, and that's faster. 
       * taucs_transpose_vec_times_matrix also incorporates the 
       * selection of columns of G from which to form A_G so that 
       * we do not have to incur the computational expense of creating
       * a submatrix.
       */
      taucs_transpose_vec_times_matrix(residual, A_original_ordering, G, sizeG, yg_raw);
      
      fix_zeros(yg_raw, sizeG, rcond, inPrintErrorWarnings);
      
      /* We're now done. It remains to scatter the entries in xf_raw
       * and yg_raw over the x and y vectors, and free everything that
       * we can manage to free. 
       */
      
#ifdef HAVE_MEMSET
      
      memset(x,0,sizeof(double)*n);
      memset(y,0,sizeof(double)*n);
      
#else  /* Work around it. */
      
      for(i=0;i<n;i++) { x[i] = 0; y[i] = 0; }
      
#endif

      /* This was done with bzero, but that was scratched for 
	 portability reasons. */
      
      /* bzero(x, sizeof(double)*n); Scratched bzero calls for portability */
      /* bzero(y, sizeof(double)*n); */
      
      for(i=0;i<sizeF;i++) 
	{ 
	  x[F[i]] = xf_raw[i]; 
	}
      for(i=0;i<sizeG;i++) 
	{ 
	  y[G[i]] = yg_raw[i]; 
	}
      
      free(yg_raw);
      
      if( sizeF != 0 )
       {
      	  free(xf_raw);
       }

       free(residual);

      sizeH1 = sizeH2 = 0;
    } // for
  
  if( lsqrStep != 0 )
    {
      lsqr_input   *lsqr_in;
      lsqr_output  *lsqr_out;
      lsqr_work    *lsqr_work;
      lsqr_func    *lsqr_func;
      int bItr;
      
      alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, Af->m, Af->n );
      
      /* we let lsqr() itself handle the 0 values in this structure */
      lsqr_in->num_rows = Af->m;
      lsqr_in->num_cols = Af->n;
      lsqr_in->damp_val = 0;
      lsqr_in->rel_mat_err = 0;
      lsqr_in->rel_rhs_err = 0;
      lsqr_in->cond_lim = 1e16;
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
      
      for( bItr=0; bItr<Af->n; bItr++ ) // not really bItr here, but hey
	x[F[bItr]] = lsqr_out->sol_vec->elements[bItr];
		
      free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
    }
  
  if( outResidualNorm != NULL )
    {
      double* finalresidual = (taucs_double *)calloc(m,sizeof(taucs_double));
      ourtaucs_ccs_times_vec(A_original_ordering,x,finalresidual);

      //cblas_daxpy(m,-1.0,b, 1, finalresidual, 1);
      int incX;
      alpha = -1; incX = 1; incY = 1; 
      DAXPY_F77(&m,&alpha,b,&incX,finalresidual,&incY);

      //*outResidualNorm = cblas_dnrm2(m, finalresidual, 1);
      *outResidualNorm = DNRM2_F77(&m,finalresidual,&incX);

      free(finalresidual);
    }
  
  /* We conclude with a little more memory freeing. */
  taucs_ccs_free(Af);
  free(F); 
  free(G);     
  free(H1); 
  free(H2);
  taucs_ccs_free(AprimeDotA);
  taucs_ccs_free(lsqrApA);
  
  free(y);
  
  return x;
}


/* This comparison function is needed for qsort, it does descending */
int
compare_taucs_doubles (const void *a, const void *b)
{
  const taucs_double *da = (const taucs_double *) a;
  const taucs_double *db = (const taucs_double *) b;
  
  return (*da < *db) - (*da > *db);
}



// This is an implementation of the block3 algorithm due to Adlers
// See http://citeseer.ist.psu.edu/385071.html
// Sparse Least Squares Problems with Box Constraints
// by  Mikael Adlers
// The supposed advantage of the block3 algorithm over the t_snnls
// algorithm is that block3 is suppose to avoid cycling behavior.

#define DEBUG 1

taucs_double*
t_snnls( taucs_ccs_matrix *A_original_ordering, taucs_double *b, 
	 double* outResidualNorm, double inRelErrTolerance, int inPrintErrorWarnings )
{
  taucs_ccs_matrix  *Af;
  int               m,n,i, maxSize;
  int              A_cols;
  int              *F, *G, *H1, *H2;
  int              sizeF, sizeG, sizeH1, sizeH2;
  int				lsqrStep=0;
  double			rcond=1;


  taucs_double *p;
  taucs_double *alpha;
  double qofx, qofxplusalphap;
  /*double bb;*/

  // these are just to speed up some stuff below ever so slightly
  taucs_double tmp;
  int itmp;

  // I suppose one could double use incTmp and alphadog, but it 
  // shouldn't make much of a difference
  double incTmp = {-1.0};
  int incX = {1}, incY = {1};
  double alphadog = {-1.0};
  int alphaItr = {0};

  /* These variables are subsets of the column indices of the matrix A, 
   * always stored in sorted order, and sized by the corresponding
   * "size" integers. 
   * 
   * Like the row indices themselves, they are 0-based. 
   */
  taucs_double     *x, *y, *xf_raw = NULL, *yg_raw = NULL, *residual = NULL;
  taucs_double *Apb, *ApAx, *xplusalphap, *ApAxplusalphap;
  
  int AprimeDotA_cols;
  
  taucs_ccs_matrix* AprimeDotA = taucs_ccs_aprime_times_a(A_original_ordering);
  taucs_ccs_matrix*   lsqrApA;

  if (gVERBOSITY >= 10) {  /* In veryverbose mode, we spit out debugging crap. */

    FILE *outmat, *outb;
    if ((outmat = fopen("tsnnlsA.sparse","w")) != NULL && 
	(outb = fopen("tsnnlsb.mat","w")) != NULL) {

      taucs_ccs_write_sparse(outmat,A_original_ordering);   
      colvector_write_mat(outb,b,A_original_ordering->m,"b");
      fclose(outmat);
      fclose(outb);

      printf("\ntsnnls: Wrote A and b to tsnnlsA.sparse and tsnnlsb.mat.\n");
    
    } 

    printf("tsnnls: Called with \n"
	   "          %d x %d matrix A.\n"
	   "          %p      buffer b.\n"
	   "          %p      outResidualNorm.\n"
	   "          %g      inRelErrTolerance\n"
	   "          %d      inPrintErrorWarnings.\n\n",
	   A_original_ordering->m,A_original_ordering->n,
	   b,outResidualNorm,inRelErrTolerance,inPrintErrorWarnings);
   
    printf("tsnnls: Created %d x %d matrix AprimeDotA.\n",AprimeDotA->m,AprimeDotA->n);

  }

  /* create a copy of AprimeDotA memory wise to store the tlsqr submatrices */

  lsqrApA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  lsqrApA->n = AprimeDotA->n;
  lsqrApA->flags = TAUCS_DOUBLE;
  lsqrApA->flags = lsqrApA->flags | TAUCS_SYMMETRIC;
  lsqrApA->flags = lsqrApA->flags | TAUCS_LOWER; // rep the lower half
  lsqrApA->colptr = (int*)malloc(sizeof(int)*(lsqrApA->n+1));
  /* This is the number of nonzeros in A'*A, 
     which we cannot overflow with a submatrix */
  maxSize = AprimeDotA->colptr[AprimeDotA->n]; 
  lsqrApA->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
  lsqrApA->rowind = (int*)malloc(sizeof(int)*maxSize);
  
  if( inRelErrTolerance <= 0.0 )
    lsqrStep = 1;
  
  // A_rows = A_original_ordering->m;
  A_cols = A_original_ordering->n;
  
  AprimeDotA_cols = A_cols;
  
  m = A_original_ordering->m;
  n = A_original_ordering->n;
  
  /* We first allocate space. */
  F   = calloc(n,sizeof(int));
  G   = calloc(n,sizeof(int));
  H1  = calloc(n,sizeof(int));
  H2  = calloc(n,sizeof(int));
  
  x    = calloc(n,sizeof(taucs_double));
  y    = calloc(m,sizeof(taucs_double));
  Apb  = calloc(n,sizeof(taucs_double));
  ApAx  = calloc(n,sizeof(taucs_double));
  ApAxplusalphap = calloc(n,sizeof(taucs_double));
  xplusalphap  = calloc(n,sizeof(taucs_double));
  p = calloc(n,sizeof(taucs_double));
  alpha = calloc(n,sizeof(taucs_double));


  /* submatrix allocation actually takes bit of time during profiling,
   * so we reuse an allocation that cannot be overflowed by
   * submatrices. Note that
   * A_original_ordering->colptr[A_original_ordering->n] is the number
   * of nonzero entries in A
   */

  Af = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  Af->colptr = (int*)malloc(sizeof(int)*(A_cols+1));
  Af->rowind = (int*)malloc(sizeof(int)*
			    (A_original_ordering->colptr[A_original_ordering->n]));
  Af->values.d = (double*)malloc(sizeof(double)*
				 A_original_ordering->colptr[A_original_ordering->n]);
  
  /* Next we initialize variables, Adlers suggests starting with everything
     in the free set.*/
#ifdef HAVE_MEMSET
      
  memset(x,0,sizeof(taucs_double)*n);
  for(i=0; i<n; i++){
    F[i] = i;
  }
      
#else  /* Work around it. */
      
  for(i=0;i<n;i++){ 
    x[i] = 0.0;
    F[i] = i;
  }
  
#endif

  sizeG = 0; sizeF = n;  
  
  /* here we'll precompute A'b since we have F filled up with
     all of the columns */
  /* Set y = A'b, which is the same as y=b'A. We perform that 
     computation as it is faster */

  taucs_transpose_vec_times_matrix(b,A_original_ordering, F, n, Apb);
  
  int gflag = {1};
  int fflag;
  
  while(gflag != 0){

    fflag = 1;
    
    if (gVERBOSITY >= 10) { printf("tsnnls: g loop\n"); }

    while(fflag != 0){

      if (gVERBOSITY >= 10) { 

	printf("tsnnls: \t f loop\n"); 
        printf("A_original_ordering is an %d x %d matrix.\n",
	       A_original_ordering->m,A_original_ordering->n);
        //taucs_ccs_write_sparse(stdout,A_original_ordering);
	printf("--------\n");

      }

      /* ***************************************** */
      /* solve for xf_raw in unconstrained problem */
      /* ***************************************** */
      
      taucs_ccs_submatrix(A_original_ordering, F, sizeF, Af);
      
      if( sizeF != 0 ){
	/* we compute these values based on selections based on F
	 * since it's faster than recalculating them in lsqr. This
	 * requires the use of a custom lsqr solver that expects
	 * extra parameters from snnls, however.
	 */
	
	selectAprimeDotAsparse(AprimeDotA, F, sizeF, lsqrApA); 	
	
	/* Now we include a little debugging code. */

	if (gVERBOSITY >= 10) {

	  printf("tsnnls: \t Checking inputs to t_lsqr.\n");
	  printf("tsnnls: b is an %d-vector.\n\n",Af->m);
	  //colvector_write_mat(stdout,b,Af->m,"b");

	  printf("tsnnls: Af in an %d x %d matrix.\n\n",Af->m,Af->n);
	  //taucs_ccs_write_sparse(stdout,Af);

	}
	  
	assert(xf_raw == NULL);

	if( inRelErrTolerance > 1 || 
	    (lsqrStep != 0 && inPrintErrorWarnings == 0) ) {

	  if (gVERBOSITY >= 10) { printf("tsnnls: \t Calling t_snnlslsqr.\n"); }
	  xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, NULL);		

	} else {

	  if (gVERBOSITY >= 10) { 
	    printf("tsnnls: \t Calling t_snnlslsqr w/rcond.\n");
	  }
	  
	  xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, &rcond );
	  if( (1/rcond)*(1/rcond)*__DBL_EPSILON__ < inRelErrTolerance )
	      lsqrStep = 1;
	}
	
	if( xf_raw == NULL )
	  return NULL; // matrix probably not positive definite

	if (gVERBOSITY >= 10) { 

	  int xrcnt;

	  printf("tsnnls: \t ptr xf_raw = %p. Dumping xf_raw. \n",xf_raw); 
	  
	}

      }
      else{	  
	/* if sizeF is 0, then we need to go to the outer loop */
	fflag = 0;
	break;
      }

      /* **************************************************** */
      /* Compute p = xf_raw - x, but all in the right places  */
      /* **************************************************** */

#ifdef HAVE_MEMSET
      
      memset(p,0,sizeof(taucs_double)*n);
      
#else  /* Work around it. */
      
      for(i=0;i<n;i++){ 
	p[i] = 0.0;
      }
#endif

      if (gVERBOSITY >= 10) { printf("tsnnls: \t Spreading x values over p.\n"); }

      /* we also compute the alpha[i] values while we are in here */
      /* alpha will be a sizeF vector */
      for(i=0; i<sizeF; i++){
	itmp = F[i];
	p[itmp] = xf_raw[i] - x[itmp];

	// Alternately, could use some tolerance
	if(p[itmp] != 0){
	  alpha[i] = xf_raw[i]/p[itmp];
	}
	else{
	  // If the old and the new are the same, then we'll set
	  // alpha[i] to 88.0 and this will just get skipped
	  alpha[i] = 88.0;
	}
      }


      /* ******************************************************* */
      /* We know the alpha_i values, so determine biggest step   */
      /* in the direction p which reduces the value of q.  We    */
      /* skip anything bigger than one since that will be picked */
      /* up by the first time through the loop.                  */
      /* ******************************************************* */

      /* we'll reuse the value of q(x), so we might as well hold on to it */

      /* So here's the deal, we want to compute q(x) which is
	 q(x) = 1/2*x'A'Ax-x'A'b = x'(0.5*A'Ax - A'b)
	 We know A'A stored as AprimeDotA
	 We'll assume that ourtaucs_ccs_times_vec is quicker for computing
	 A'Ax and then we'll dot with x manually.
	 For the right side, we've precomputed A'b, so we just need to dot
	 with x.
      */

      if (gVERBOSITY >= 10) { printf("tsnnls: \t Calling ourtaucs_ccs_times_vec\n");}

      // Note: x has not been updated, so it is the version from
      // the last time through the loop, unless we swapped stuff from F to G, 
      // in which case we zeroed some stuff in A.

      ourtaucs_ccs_times_vec(AprimeDotA,x,ApAx);
      
      /* Note: there might be some fanciness that makes this quicker
	 computing q(x) and the x+alpha*p value for the first 
	 run through the loop, where alpha=1.0
      */

      qofx = 0.0;

      for(i=0; i<n; i++){
	qofx += x[i] * (0.5*ApAx[i]-Apb[i]);
      }

      // printf("qofx %f\n",qofx);

      for(i=0; i<n; i++) {
	xplusalphap[i] = x[i] + p[i];
      }

#ifdef HAVE_MEMSET

      memset(ApAxplusalphap,0,sizeof(taucs_double)*n);

#else
      for(i=0; i<n; i++){
	ApAx[i] = 0.0;
      }
#endif

      ourtaucs_ccs_times_vec(AprimeDotA,xplusalphap,ApAxplusalphap);

      qofxplusalphap = 0.0;
      
      for(i=0; i<n; i++){
	qofxplusalphap += xplusalphap[i] * (0.5*ApAxplusalphap[i]-Apb[i]);
      }

      // printf("qofxplusalphap: %g.\n",qofxplusalphap);

      /* if we have improvement in q, we need not do the following */
      alphaItr = 0;

      if(qofxplusalphap > qofx){

	if (gVERBOSITY >= 10) { printf("tsnnls: \t Calling qsort.\n"); }
	/* darn, that didn't work, so we need to sort the alpha's */
	qsort(alpha,sizeF,sizeof(taucs_double),compare_taucs_doubles);

	int tmpitr;

	for(tmpitr = 0;tmpitr < sizeF && tmpitr < 10;tmpitr++) {printf("%.2g",alpha[i]);}
	printf("\n");

	/* burn the ones where alpha >= 1 */
	alphaItr = 0;
	while(alpha[alphaItr] >= 1.0){
	  alphaItr++;
	}

	/* now we see if the step of size alpha[i] improves q */
	for(; alphaItr<sizeF; alphaItr++){
	  tmp = alpha[alphaItr];

	  for(i=0; i<n; i++){
	    xplusalphap[i] = x[i] + tmp*p[i];
	  }

#ifdef HAVE_MEMSET

	  memset(ApAxplusalphap,0,sizeof(taucs_double)*n);

#else
	  for(i=0; i<n; i++){
	    ApAx[i] = 0.0;
	  }
#endif

	  ourtaucs_ccs_times_vec(AprimeDotA,xplusalphap,ApAxplusalphap);

	  qofxplusalphap = 0.0;
      
	  for(i=0; i<n; i++){
	    qofxplusalphap += xplusalphap[i] * (0.5*ApAxplusalphap[i]-Apb[i]);
	  }
	}

	if(alphaItr >= sizeF){
	  // I don't know if this is an error, none of the alpha_i worked
	  printf("Warning: No improvement with alphas\n");

	  // This is a total hack, but it hasn't happened ever in 
	  // trivial levels of testing.  We'll assume it happens
	  // so infrequently, that this doesn't affect performance.
	  // Basically, this is going to force a double copy.
	  memcpy(xplusalphap,x,sizeof(taucs_double)*n);
	}

      }// end (qofxplusalphap > qofx)

      /* ************************************************* */
      /* well, we got a reduction in q, so now we update x */
      /* for the next (potential) round                    */
      /* ************************************************* */

      memcpy(x,xplusalphap,sizeof(taucs_double)*n);
      

      /* ************************************************* */
      /* now we see if any of the frees want to be bounded */
      /* ************************************************* */

      if (gVERBOSITY >= 10) { 

	int vcnt;

	printf("tsnnls: \t Checking infeasibles.\n"); 
	printf("tsnnls: Dumping F of size %d.\n",sizeF);

	/*for(vcnt = 0;vcnt < sizeF;vcnt++) {

	  printf("%d ",F[vcnt]);

	}
	
	printf("\n"); */

      }

      infeasible(F,x,sizeF,H1,&sizeH1);
      
      if(sizeH1 == 0){

	fflag = 0;

	// We are going to leave the loop, so squash xf_raw.
	free(xf_raw);
	xf_raw = NULL;

	if (gVERBOSITY >= 10) { printf("tsnnls: H is empty.\n"); }

	break;
      }
      else{

	if (gVERBOSITY >= 10) { printf("tsnnls: H is nonempty.\n"); }

	int_difference(F,sizeF,H1,sizeH1,&sizeF);
	int_union(G,sizeG,H1,sizeH1,&sizeG);

	// Need to zero out the stuff that turned up negative
	for(i=0; i<sizeH1; i++){
	  x[H1[i]] = 0.0;
	}
      }

      // The first thing we're going to do up top is allocate a new xf_raw.
      // So we free the old one here.

      free(xf_raw);
      xf_raw = NULL;

    } // end inner fflag loop

    /* At this point, we should be at a stationary point for x_F, having bound
       a bunch of formerly free variables. We recompute and print qofx. */

    ourtaucs_ccs_times_vec(AprimeDotA,x,ApAx);
    qofx = 0.0;
    
    for(i=0; i<n; i++){
      qofx += x[i] * (0.5*ApAx[i]-Apb[i]);
    }
    
    printf("qofx at stationary point: %g.\n",qofx);

    /* We left the inner loop because sizeH was zero. This means that F has 
       not changed during this iteration. We confirm that the equation

       A'_F A_F x_F - A'_F b + A'_F A_B x_B = 0.

       Since the bound guys are bound to zero, x_B is zero, and we need only
       check that 

       A'_F A_F x_F - A'_F b = 0. */

    double normgF = 0.0;


    /* INSERT CODE TO COMPUTE THIS HERE!!!! */

    /* Now we need to compute y_G and see if shifting anything out
       of F has created infeasibles in G */

    if (gVERBOSITY >= 10) { printf("tsnnls: F loop terminated. Computing y_g.\n");}

    taucs_ccs_submatrix(A_original_ordering, F, sizeF, Af);

    // Note: it might be simpler to allocate residual up front
    // then we'd zero it out here
    if(sizeF != 0){
      
      // update xf_raw to include alpha*p. Since we are copying
      // into a new xf_raw, we reallocate it here.
      
      xf_raw = malloc(sizeF*sizeof(double));
      
      for(i=0; i<sizeF; i++){
	xf_raw[i] = x[F[i]];
      }      
      
      /* Now compute the residual A_F x_F - b. This is an m-vector. */
      assert(residual == NULL);      
      residual = (taucs_double *)calloc(m,sizeof(taucs_double));
      ourtaucs_ccs_times_vec(Af,xf_raw,residual);
      
      free(xf_raw);  // we won't use it again below
      xf_raw = NULL; // for safety's sake.

    } else{	  
      /* 
       * if sizeF is 0, the meaning of residual changes (since
       * there really is _no_ matrix), so we'll just set the
       * residual to -b, but we still need a zeroed residual to do
       * the below computation to make that happen, which calloc
       * does here
       */
      assert(residual == NULL);
      residual = (taucs_double *)calloc(m,sizeof(taucs_double));
    }

    DAXPY_F77(&m,&incTmp,b,&incX,residual,&incY);

    // Note: We could allocate this up front as well
    /* We now compute (A_G)'. */
    /* And finally take (A_G)'*residual. This is a sizeG-vector. */
    assert(yg_raw == NULL);
    yg_raw = (taucs_double *)calloc(sizeG,sizeof(taucs_double));     

    /* 
     * We now should compute (A_G)', and take (A_G)'*residual, but 
     * A_G'*residual = residual'*A_G, and that's faster. 
     * taucs_transpose_vec_times_matrix also incorporates the 
     * selection of columns of G from which to form A_G so that 
     * we do not have to incur the computational expense of creating
     * a submatrix.
     */
    
    if (gVERBOSITY >= 10) { printf("tsnnls: Computing residual.\n"); }

    taucs_transpose_vec_times_matrix(residual, A_original_ordering, G, sizeG, yg_raw);

    /* This was the last time we used residual, so let's kill it. */
    free(residual); residual = NULL;

    /* Note, this shouldn't change, this was a check during debugging
    infeasible(F,x,sizeF,H1,&sizeH1);  

    if(sizeH1 > 0){
      // error, error, this shouldn't change
      printf("sizeH1 > 0 after the inner loop\n");
      exit(-1);
    }
    */

    // Note: if we rewrote infeasible, we could avoid computing y at all.
    // But for the sake of ease, we'll leave it like this

    /* here we are setting up y */
    for(i=0; i<sizeF; i++){
      y[F[i]] = 0.0;
    }
    
    for(i=0; i<sizeG; i++){
      y[G[i]] = yg_raw[i];
    }

    /* From here on out, we will only use y. So we discard yg_raw. */
    free(yg_raw); yg_raw = NULL;
    
    infeasible(G,y,sizeG,H2,&sizeH2);

    if(sizeH2 == 0){
      gflag = 0;
      break;
    }
    else{
      int_union(F,sizeF,H2,sizeH2,&sizeF);
      int_difference(G,sizeG,H2,sizeH2,&sizeG);
    }
  } // while gflag

  if (gVERBOSITY >= 10) { printf("tsnnls: G loop terminated.\n"); }

  if( lsqrStep != 0 )
    {

      if (gVERBOSITY >= 10) { printf("tsnnls: Doing lsqr step.\n"); }

      lsqr_input   *lsqr_in;
      lsqr_output  *lsqr_out;
      lsqr_work    *lsqr_work;
      lsqr_func    *lsqr_func;
      int bItr;
      
      alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, Af->m, Af->n );
      
      /* we let lsqr() itself handle the 0 values in this structure */
      lsqr_in->num_rows = Af->m;
      lsqr_in->num_cols = Af->n;
      lsqr_in->damp_val = 0;
      lsqr_in->rel_mat_err = 0;
      lsqr_in->rel_rhs_err = 0;
      lsqr_in->cond_lim = 1e16;
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
      
      for( bItr=0; bItr<Af->n; bItr++ ) // not really bItr here, but hey
	x[F[bItr]] = lsqr_out->sol_vec->elements[bItr];
		
      free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );

      if (gVERBOSITY >= 10) { printf("tsnnls: Survived lsqr.\n"); }

    }
  
  if( outResidualNorm != NULL )
    {

      if (gVERBOSITY >= 10) { printf("tsnnls: Computing final residual.\n"); }

      double* finalresidual = (taucs_double *)calloc(m,sizeof(taucs_double));
      ourtaucs_ccs_times_vec(A_original_ordering,x,finalresidual);

      //cblas_daxpy(m,-1.0,b, 1, finalresidual, 1);
      // int incX, alphadog;
      alphadog = -1; incX = 1; incY = 1; 
      DAXPY_F77(&m,&alphadog,b,&incX,finalresidual,&incY);

      //*outResidualNorm = cblas_dnrm2(m, finalresidual, 1);
      *outResidualNorm = DNRM2_F77(&m,finalresidual,&incX);

      free(finalresidual);
    }
  // free memory

  free(F);
  free(G);
  free(H1);
  free(H2);
  taucs_ccs_free(AprimeDotA);
  taucs_ccs_free(lsqrApA);
  taucs_ccs_free(Af);

  free(y);
  free(Apb);
  free(ApAx);
  free(ApAxplusalphap);
  free(xplusalphap);
  free(p);
  free(alpha);

  if (gVERBOSITY >= 10) { printf("tsnnls: Done.\n"); }

  return x;

}


//#pragma mark -

#ifndef taucs_add
#define taucs_add(x,y) ((x)+(y))
#endif
#ifndef taucs_mul
#define taucs_mul(x,y) ((x)*(y))
#endif

double
taucs_rcond( taucs_ccs_matrix* A )
{
  #ifndef HAVE_ATLAS_LAPACK 

  /* We have two different versions of this function. In this version,
     we have a full LAPACK, so we use dgecon and dpocon to get the
     job done. */
  
  char	NORM = '1';
  ACINT32_TYPE	N = 0,AN = 0;     /* LAPACK expects 32 bit ints */
  ACINT32_TYPE LDA = 0;
  double  ANORM = 0;
  double  RCOND = 0;
  double* WORK = NULL;
  ACINT32_TYPE*   IWORK = NULL;
  ACINT32_TYPE	INFO;
  ACINT32_TYPE*	IPIV = NULL;
  double* lapackA = NULL;
  
  /* Construct LAPACK representation of A and compute the 1 norm of A */
  int vSize;
  int cItr, rItr;
  //long int rowCount = A->m;
  ACINT32_TYPE rowCount = A->m;
  
  double localMax = 0;
  
  if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
    {
      vSize = A->n*A->n;
      rowCount = A->n;
    }
  else
    vSize = A->m*A->n;
  
  lapackA = (double*)calloc(vSize,sizeof(double));
  
  /*lapackA = (double*)malloc(sizeof(double)*vSize);
    bzero(lapackA, sizeof(double)*vSize); */
  
  for( cItr=0; cItr<A->n; cItr++ )
    {
      localMax = 0;
      for( rItr=A->colptr[cItr]; rItr<A->colptr[cItr+1]; rItr++ )
	{
	  int index = -1;
	  index = A->rowind[rItr] + cItr*rowCount;
	  if( index > vSize )
	    {
	      fprintf( stderr, "Rcond memory error!\n" );
	      exit(-1);
	    }
	  lapackA[index] = A->values.d[rItr];
	  localMax += fabs(A->values.d[rItr]);
	}
      if( localMax > ANORM )
	ANORM = localMax;
    }
  
  NORM = '1';
  N = A->n;
  AN = A->n;
  LDA = A->m;
  RCOND = 0;
  WORK = (double*)malloc(sizeof(double)*4*N);
  IWORK = (ACINT32_TYPE*)malloc(sizeof(ACINT32_TYPE)*N);
  INFO = 0;
  
  IPIV = (ACINT32_TYPE*)malloc(sizeof(ACINT32_TYPE)*min(rowCount, A->n));
  
  dgetrf_( &rowCount, &AN, lapackA, &rowCount, IPIV, &INFO );
  dgecon_( &NORM, &N, lapackA, &LDA, &ANORM, &RCOND, WORK, IWORK, &INFO );
  
  free(IPIV);
  free(IWORK);
  free(WORK);
  free(lapackA);
  
  return RCOND;

  #else 

  /* We have only a limited ATLAS LAPACK available. In this case, we use lsqr to 
     estimate the condition number of A. */
  
  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  
  double        rcond;
  
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, A->m, A->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = A->m;
  lsqr_in->num_cols = A->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = 0;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1e16;
  lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 1000;
  lsqr_in->lsqr_fp_out = NULL;	
  
  for( bItr=0; bItr<A->m; bItr++ )
    {
      lsqr_in->rhs_vec->elements[bItr] = 1;  /* We will solve Ax = [1 1 .... 1]^T. */
    }
  /* Here we set the initial solution vector guess, which is 
   * a simple 1-vector. You might want to adjust this value for fine-tuning
   * t_snnls() for your application
   */
  for( bItr=0; bItr<A->n; bItr++ ) {
    lsqr_in->sol_vec->elements[bItr] = 1; 
  }
  
  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = sparse_lsqr_mult;
  
  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, A );  
  rcond = 1/lsqr_out->mat_cond_num;
  
  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
  
  return rcond;
  
#endif
}

void
transpose_vec_times_matrix(double* b, double* A, int* F, int A_cols, \
			   int rows, int cols, double* result)
{
  // result = b'*A
  int cItr;
  int incX = {1};
  int incY;
  int N;

  N = rows;
  incY = A_cols;

  for( cItr=0; cItr<cols; cItr++ ) {
   
    //result[cItr] = cblas_ddot( rows, b, 1, &A[F[cItr]], A_cols );
    result[cItr] = DDOT_F77(&N,b,&incX,&A[F[cItr]],&incY); 
  }
}

// here cols doubles for sizeF, since F is a selection of columns from A
void
taucs_transpose_vec_times_matrix(double* b, taucs_ccs_matrix* A, int* F, int cols, 
				 double* result)
{
  int cItr, rItr;
  for( cItr=0; cItr<cols; cItr++ ) {

    result[cItr] = 0;

    for( rItr=0; rItr<A->colptr[F[cItr]+1]-A->colptr[F[cItr]]; rItr++ ) {

      int tRow = A->rowind[A->colptr[F[cItr]] + rItr];	  
      result[cItr] += b[tRow] * A->values.d[A->colptr[F[cItr]] + rItr];
    
    }
  }
}

void
transpose_vec_times_matrix_nosub(double* b, double* A, int A_cols, int rows, double* result)
{
  // result = b'*A
  int cItr;
  int N;
  int incX = {1}, incY;

  N = rows; incY = A_cols;

  for( cItr=0; cItr<A_cols; cItr++ ) {
    //result[cItr] = cblas_ddot( rows, b, 1, &A[cItr], A_cols );
    result[cItr] = DDOT_F77(&N,b,&incX,&A[cItr],&incY);
  }
}

void
taucs_transpose_vec_times_matrix_nosub(double* b, taucs_ccs_matrix* A, double* result)
{
  int cItr, rItr;
  for( cItr=0; cItr<A->n; cItr++ )
    {
      result[cItr] = 0;
      for( rItr=0; rItr<A->colptr[cItr+1]-A->colptr[cItr]; rItr++ )
	{
	  int tRow = A->rowind[A->colptr[cItr] + rItr];
	  result[cItr] += b[tRow] * A->values.d[A->colptr[cItr] + rItr];
	}
    }
}

void
ourtaucs_ccs_times_vec( taucs_ccs_matrix* m, 
			 taucs_datatype* X,
			 taucs_datatype* B )
{
	int i,ip,j,n, rows;
	taucs_datatype Aij;

	n = m->n;
	rows = m->m;
	if( (m->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
	{
		// this is a total hack, but taucs's thingy works when 
		// things are symmetric. keep in mind, however, that otherwise, 
		// it does not
		taucs_ccs_times_vec( m, X, B );
		return;
	}
		
	for(i=0; i < rows; i++) 
		B[i] = 0;

	for (j=0; j<n; j++) {
      for (ip = (m->colptr)[j]; ip < (m->colptr[j+1]); ip++) {
		i   = (m->rowind)[ip];
		Aij = (m->values.d)[ip];
		
		B[i] = taucs_add(B[i],taucs_mul(X[j],Aij));
      }
    }
}

void
taucs_ccs_submatrix( const taucs_ccs_matrix* A, const int* keptCols, const int inColCount, taucs_ccs_matrix* result)
{
	int cItr, colOffset, c2;
	
	result->m = A->m;
	result->n = inColCount;
	//result->flags = A->flags;
	result->flags = TAUCS_DOUBLE;
			
	colOffset = 0;
	for( cItr=0; cItr<inColCount; cItr++ )
	{
		result->colptr[cItr] = colOffset;
		for( c2=A->colptr[keptCols[cItr]]; c2<A->colptr[keptCols[cItr]+1]; c2++ )
		{
			result->rowind[colOffset] = A->rowind[c2];
			result->values.d[colOffset] = A->values.d[c2];
			colOffset++;
		}
	}
	result->colptr[cItr] = colOffset;
}



taucs_ccs_matrix*
taucs_ccs_transpose( const taucs_ccs_matrix* A )
{
  taucs_ccs_matrix* result = NULL;
  double* values = NULL;
  int cItr, rItr, colOffset;
	
  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
	
  result->m = A->n;
  result->n = A->m;
	
  result->flags = A->flags;
	
  // we sacrifice some memory here so we don't have to do the number of 
  // nonzero entries computation
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*result->m*result->n);
  result->values.d = (double*)malloc(sizeof(double)*result->m*result->n);
	
  values = taucs_convert_ccs_to_doubles(A);
	
  colOffset = 0;
  for( rItr=0; rItr<A->m; rItr++ )
    {
      result->colptr[rItr] = colOffset;
      for( cItr=0; cItr<A->n; cItr++ )
	{
	  double v = values[(rItr*A->n)+cItr];
	  if( v != 0 )
	    {
	      result->rowind[colOffset] = cItr;
	      result->values.d[colOffset] = v;
	      colOffset++;
	    }
	}
    }
  result->colptr[rItr] = colOffset;
	
  free(values);
	
  return result;
}

taucs_ccs_matrix*
taucs_construct_sorted_ccs_matrix( double* values, int rowsize, int rows )

/* Converts an array of doubles in C (row-major) storage order to a taucs_ccs_matrix. */
/* That is, we expect the matrix to read 

      A(1,1) A(1,2) A(1,3) ... A(1,rowsize) A(2,1) ... A(rows, 1) .. A(rows, rowsize)
     
*/

{
  taucs_ccs_matrix* result = NULL;
  int nnz = 0;
  int i, rItr, cItr, colOffset;
  double v;
	
  for( i=0; i<rowsize*rows; i++ )
    {
      if( values[i] != 0 )
		nnz++;
    }
	
  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  result->n = rowsize;
  result->m = rows;
  result->flags = TAUCS_DOUBLE;
	
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*nnz);
  result->values.d = (double*)malloc(sizeof(taucs_double)*nnz);
	
  colOffset = 0;
	
  for( cItr=0; cItr<rowsize; cItr++ )
    {
      result->colptr[cItr] = colOffset;
      for( rItr=0; rItr<rows; rItr++ )
	{
	  v = values[rItr*rowsize+cItr];
	  if( v != 0 )
	    {
	      result->rowind[colOffset] = rItr;
	      result->values.d[colOffset++] = v;
	    }
	}
    }
  result->colptr[cItr] = colOffset;
	
  return result;
}

double*
taucs_convert_ccs_to_doubles( const taucs_ccs_matrix* A )
{
	int vSize;
	int cItr, rItr;
	double* values;
	
	if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
		vSize = A->n*A->n;
	else
		vSize = A->m*A->n;
	
	values = (double*)calloc(vSize,sizeof(double));
	
	/*values = (double*)malloc(sizeof(double)*vSize);	
	  bzero(values, sizeof(double)*vSize); */
	
	int rowSize = A->n;
	
	for( cItr=0; cItr<A->n; cItr++ )
	{
		for( rItr=A->colptr[cItr]; rItr<A->colptr[cItr+1]; rItr++ )
		{
			int index = -1;
			index = rowSize*A->rowind[rItr] + cItr;
			values[index] = A->values.d[rItr];
		}
	}
	
	return values;
}

void
taucs_print_ccs_matrix( const taucs_ccs_matrix* A )
{
  double* v = taucs_convert_ccs_to_doubles(A);
  int rowCount;
  int rItr, cItr;
	
  rowCount = A->m;
	
  if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
    {
      printf( "Matrix flagged symmetric\n" );
      rowCount = A->n;
    }
	
  for( rItr=0; rItr<rowCount; rItr++ )
    {
      for( cItr=0; cItr<A->n; cItr++ )
	printf( "%5.4g ", v[rItr*A->n + cItr] );
      printf( "\n" );
    }
	
  free(v);
}

void taucs_ccs_write_sparse( FILE *fp, taucs_ccs_matrix *A)

/* Writes a taucs_ccs matrix to the SPARSE format read by tsnnls_test. */

{
  int cpItr, rwItr;

  if (fp == NULL) {

    printf("tsnnls: Can't write to NULL file pointer.\n");
    exit(1);

  }

  fprintf(fp, "SPARSE\n %d %d\n %d\n",A->m,A->n,A->colptr[A->n]);
  
  for(cpItr=0,rwItr=0;cpItr<A->n;cpItr++) {

    for(;rwItr<A->colptr[cpItr+1];rwItr++) {

      /* Must adjust for 1-based indexing. */

      fprintf(fp,"%d %d %10.16g\n",A->rowind[rwItr]+1,cpItr+1,
	      A->values.d[rwItr]);

    }

  }

}

void taucs_ccs_write_mat(FILE *fp, taucs_ccs_matrix *A) 

{
  double *vals;
  int rows,cols;

  rows = A->m; cols = A->n;
  vals = taucs_convert_ccs_to_doubles(A);
 
  if (vals == NULL) {

    printf("taucs_ccs_write_mat: Can't convert %d x %d matrix to doubles.\n",
	   A->m,A->n);
    exit(1);

  }
  
  fprintf(fp,
	  "# Created by tsnnls\n"
	  "# name: A\n"
	  "# type: matrix\n"
	  "# rows: %d\n"
	  "# columns: %d\n",rows,cols);

  int rItr,cItr;

  for( rItr=0; rItr<rows; rItr++ ) {
    
    for( cItr=0; cItr<cols; cItr++ ) {

      fprintf( fp, "%10.16lf ", vals[rItr*cols + cItr] );

    }
    
    fprintf( fp, "\n" );
  }
  
  free(vals);
  
}
  
void colvector_write_mat(FILE *fp, double *x, int rows, char *name)

/* Writes a column vector to an octave MAT file, with given name */
/* if specified. Otherwise, the name is set to "x" if the name   */
/* pointer is NULL. */

{
  int i;
  char *varname;
  char xn[2] = "x";

  if (name == NULL) {

    varname = xn;

  } else {

    varname = name;

  }
  
  fprintf(fp,
	  "# Created by tsnnls\n"
	  "# name: %s\n"
	  "# type: matrix\n"
	  "# rows: %d\n"
	  "# columns: 1\n",varname,rows);

  for(i=0;i<rows;i++) {
    
    fprintf(fp,"%10.16lf\n",x[i]);

  }

}
