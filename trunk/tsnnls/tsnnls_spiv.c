
 
/*
 * This program is free software distributed under the GPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/gpl.txt
 */

// This file is included at the end of tsnnls.c, so should be able to call
// procedures from that file without headers.

taucs_double *compute_lagrange_multipliers(taucs_ccs_matrix *A,
					   taucs_ccs_matrix *ATA,
					   taucs_double *x,taucs_double *b,
					   int nBound,int *Bound)

/* Computes Lagrange multipliers for the bound variables in A, using the variable ATA,
   which should be A transpose X A. In A is an m x n matrix, we expect x to be an n x 1 vector. */

{
  taucs_double *ATAx, *ATb,*y;
  int N=A->n,incX=1,incY=1,i;
  double alpha=-1;

  ATAx = malloc(sizeof(taucs_double)*A->n);
  ATb  = malloc(sizeof(taucs_double)*A->n);
  assert(ATAx != NULL && ATb != NULL);

  /* Compute y = -(A^T(b - Ax))^T = -(b - Ax)^T A = -b^T A + x^T (A^T A). */ 

  taucs_transpose_vec_times_matrix_nosub(b,A,ATb);      
  taucs_transpose_vec_times_matrix_nosub(x,ATA,ATAx);
  DAXPY_F77(*N,*alpha,ATb,&incX,ATAx,*incY);

  /* Now select the values corresponding to bound variables. */

  y = malloc(sizeof(taucs_double)*nBound);
  assert(y != NULL);
  for(i=0;i<nBound;i++) { y[i] = ATAx[Bound[i]]; }

  /* Now free scratch memory and return. */

  free(ATAx); free(ATb);

  return y;
}

 
void P_spiv(int n,taucs_double *x,int nconstrained)

/* Projects the constrained variables in x to legal values. */

{
  int i;

  for(i=0;i<nconstrained;i++) { x[i] = max(x[i],0); }

}

taucs_ccs_matrix taucs_ccs_matrix_new(int m, int n,int flags,int nnz)

/* Allocates a ccs matrix. */

{
  taucs_ccs_matrix *A;

  A = (taucs_ccs_matrix *)malloc(sizeof(taucs_ccs_matrix));
  assert(A != NULL);

  A->n = n;
  A->m = m;
  A->flags = TAUCS_DOUBLE | flags;

  A->colptr = (int *)malloc(sizeof(int)*(nnz+1));
  A->rowind = (int *)malloc(sizeof(int)*nnz);
  A->values.d = (double*)malloc(sizeof(int)*nnz);
  
  assert((A->colptr != NULL) && (A->rowind != NULL) && (A->values.d != NULL));

  return A;
}
  

taucs_double *solve_unconstrained(taucs_ccs_matrix *A, taucs_ccs_matrix *ATA,
				  taucs_double *b,int nFree, int *Free)

/* Solves the unconstrained problem in the current free variables and spreads
   the result across a vector of length n to give the entire solution. */

{
  taucs_ccs_matrix *Afree, *ATAfree;
  taucs_double     *xFree, *x;
  int               i;
  double rcond;

  Afree = taucs_ccs_matrix_new(A->m,A->n,TAUCS_DOUBLE,A->colptr[A->n+1]);
  ATAfree = taucs_ccs_matrix_new(A->n,A->n,TAUCS_SYMMETRIC | TAUCS_LOWER,A->n*A->n);

  if ( nFree > 0 ) {

    taucs_ccs_submatrix(A,Free,nFree,Afree);
    selectAprimeDotAsparse(ATA,Free,nFree,ATAfree);
    
    xFree = t_snnlslsqr(Afree,b,ATAfree,Free,&rcond);
    
  }

  x = calloc(sizeof(taucs_double)*A->n);
  for(i=0;i<nFree;i++) { x[Free[i]] = xFree[i]; }

  taucs_ccs_free(ATAfree);
  taucs_ccs_free(Afree);

  return x;
}

taucs_double *computep(taucs_ccs_matrix *A, taucs_ccs_matrix *ATA, 
		       taucs_double *xn, taucs_double *b,
		       int nFree,int *Free)

/* We are trying to correct the current solution xn by stepping in a
   certain direction computed with respect to a new set of free
   variables given by nFree, Free. 

   To do so, we solve 

   min ||A(xn + p) - b || 
 = min || Ap - (b - Axn) || 
 = min || Ap - ( - (Axn - b)) ||. */

{
  taucs_double *Axn = malloc(A->n*sizeof(taucs_double)), *result;
  int N=A->n,incX=1,incY=1;
  double alpha=-1.0;

  taucs_transpose_vec_times_matrix_nosub(xn,A,Axn);
  DAXPY_F77(*N,*alpha,b,&incX,Axn,&incY); // Axn = Axn - b
  DSCAL_F77(*N,*alpha,Axn,&incX);         // Axn = -Axn

  result = solve_unconstrained(A,ATA,Axn,nFree,Free); 
  free(Axn);

  return result;
}


void bindzeros(int n,taucs_double *x,int *nFree,int *Free,int *nBound,int *Bound, int nconstrained)

/* Bind any free variables whose value is now zero. Assert that the free variables are legal. */

{
  int i;
  int nNewBound = 0;
  int *NewBound = calloc(sizeof(int),n);

  /* Search the free set for new variables to bind. */

  for(i=0;i<nFree;i++) {

    assert(x[Free[i]] >= -1e-16);

    if (x[Free[i]] < 1e-16 && Free[i] < nconstrained) {  
      
      // We can never bind something with an index >= nconstrained, regardless of value.

      newBound[nNewBound++] = Free[i]; 

    } 

  }

  /* Subtract these variables from the Free set and add them to the bound set. */

  int_difference(Free,nFree,NewBound,nNewBound,&nFree);
  int_union(Bound,nBound,NewBound,nNewBound,&nBound);

}
   
bool is_optimal_point( int n, taucs_double *y, int nBound, int *Bound)

/* We check the bound variables for the correct sign on their Lagrange multipliers. */

{
  int i;

  for (i=0;i<nBound;i++) {

    if (y[Bound[i]] < 0) { return (1 == 0); } // that is, return FALSE

  }

  return (1 == 1); // that is, TRUE

}

taucs_double *t_snnls_spiv (taucs_ccs_matrix *A, taucs_double *b,
			    double *outResidualNorm, double inRelErrTolerance, 
			    int inPrintErrorWarnings, int nconstrained) 

{

// This is an implementation of a single-pivoting algorithm for partially constrained
// problems. In these problems, the first "nconstrained" variables are subject to 
// non-negativity constraints, while the remaining variables are not constrained at all.

  taucs_ccs_matrix *ATA = taucs_ccs_aprime_times_a(A);
  taucs_double     *Axn = malloc(sizeof(taucs_double)*A->n);
 
  int              nFree,nBound;
  int              *Bound = calloc(sizeof(int),A->n),*Free = calloc(sizeof(int),A->n);
  
  double            last_stationary_q = {DBL_MAX};
  int               MAXPIVOT = A_original_ordering->n * 10;
  int               pivcount = 0;

  int i;

  taucs_double      *xn, *xnp1,*y;
  double            alpha = 0;

  nBound = 0; nFree = A->n;               // Set all variables free for starters.
  for(i=0;i<A->n;i++) { Free[i] = i; }
  
  xn = solve_unconstrained(A,ATA,b,nFree,Free);
  P_spiv(A->n,xn,nconstrained);
  
  bindzeros(A->n,xn,&nFree,Free,&nBound,Bound,nconstrained);
  y = compute_lagrange_multipliers(A,ATA,xn,b,nBound,Bound);
  
  while (!is_optimal_point(A->n,y,nBound,Bound)) {

    while (!isconstrainedpt) {

      /* Solve min ||A(xn + p) - b|| in the current free variables. */

      p = computep(A,ATA,xn,b,nFree,Free);
      
      alpha = findalpha(p,xn,nFree,Free);

      //CONTINUE FROM HERE!!!


  taucs_double *p;
  taucs_double *alpha;
  double qofx, newq;
  /*double bb;*/

  // these are just to speed up some stuff below ever so slightly
  taucs_double tmp = 0;
  int itmp;

  // I suppose one could double use incTmp and alphadog, but it 
  // shouldn't make much of a difference
  double minusOne = {-1.0};
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
  taucs_double *Apb, *ApAx, *xplusalphap, *ApAxplusalphap, *Pxplusalphap;
  
  int AprimeDotA_cols;
  
  taucs_ccs_matrix* AprimeDotA = taucs_ccs_aprime_times_a(A_original_ordering);
  taucs_ccs_matrix*   lsqrApA;

  clear_tsnnls_error();

  if (gVERBOSITY >= 10) {  /* In veryverbose mode, we spit out debugging crap. */

    FILE *outmat, *outb;
    if ((outmat = fopen("tsnnlsA.sparse","w")) != NULL && 
	(outb = fopen("tsnnlsb.mat","w")) != NULL) {

      taucs_ccs_write_sparse(outmat,A_original_ordering);   
      colvector_write_mat(outb,b,A_original_ordering->m,"b");
      fclose(outmat);
      fclose(outb);

      printf("\ntsnnls_spiv: Wrote A and b to tsnnlsA.sparse and tsnnlsb.mat.\n");
    
    } 

    printf("tsnnls_spiv: Called with \n"
	   "          %d x %d matrix A.\n"
	   "          %p      buffer b.\n"
	   "          %p      outResidualNorm.\n"
	   "          %g      inRelErrTolerance\n"
	   "          %d      inPrintErrorWarnings.\n"
	   "          %d      constrained variables.\n\n",
	   A_original_ordering->m,A_original_ordering->n,
	   b,outResidualNorm,inRelErrTolerance,inPrintErrorWarnings,nconstrained);
   
    printf("tsnnls_spiv: Created %d x %d matrix AprimeDotA.\n",AprimeDotA->m,AprimeDotA->n);

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
  Pxplusalphap = calloc(n,sizeof(taucs_double));

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
  for(i=0; i<n; i++){ F[i] = i; x[i] = 1.0; }
      
#else  /* Work around it. */
      
  for(i=0;i<n;i++){ 
    x[i] = 1.0;
    F[i] = i;
  }
  
#endif

  sizeG = 0; sizeF = n; sizeAlpha = 0;
  
  /* here we'll precompute A'b since we have F filled up with
     all of the columns */
  /* Set y = A'b, which is the same as y=b'A. We perform that 
     computation as it is faster */

  taucs_transpose_vec_times_matrix(b,A_original_ordering, F, n, Apb);

  int gflag = {1};
  int fflag;

  while(gflag != 0 && pivcount < MAXPIVOT){

    pivcount++;
    fflag = 1;
    
    if (gVERBOSITY >= 10) { printf("tsnnls_spiv: g loop\n"); }

    while(fflag != 0 && pivcount < MAXPIVOT){

      pivcount++;

      if (gVERBOSITY >= 10) { 

	printf("tsnnls_spiv: \t f loop\n"); 
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

	  printf("tsnnls_spiv: \t Checking inputs to t_lsqr.\n");
	  printf("tsnnls_spiv: b is an %d-vector.\n\n",Af->m);
	  //colvector_write_mat(stdout,b,Af->m,"b");

	  printf("tsnnls_spiv: Af in an %d x %d matrix.\n\n",Af->m,Af->n);
	  //taucs_ccs_write_sparse(stdout,Af);

	}
	  
	assert(xf_raw == NULL);

	if( inRelErrTolerance > 1 || 
	    (lsqrStep != 0 && inPrintErrorWarnings == 0) ) {

	  if (gVERBOSITY >= 10) { printf("tsnnls_spiv: \t Calling t_snnlslsqr.\n"); }
	  xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, NULL);		

	} else {

	  if (gVERBOSITY >= 10) { 
	    printf("tsnnls_spiv: \t Calling t_snnlslsqr w/rcond.\n");
	  }
	  
	  xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, &rcond );
	  if( (1/rcond)*(1/rcond)*__DBL_EPSILON__ < inRelErrTolerance )
	      lsqrStep = 1;
	}
	
	if( xf_raw == NULL )
	  return NULL; // matrix probably not positive definite

	if (gVERBOSITY >= 10) { 

	  printf("tsnnls_spiv: \t ptr xf_raw = %p. Dumping xf_raw. \n",xf_raw); 
	  
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

      if (gVERBOSITY >= 10) { printf("tsnnls_spiv: \t Spreading x values over p.\n"); }

      /* we also compute the alpha[i] values while we are in here */
      /* we will deduce the size of alpha as we go-- it depends on which xf_raw */
      /* guys have an infeasible value. */

      for(i=0,sizeAlpha=0; i<sizeF; i++){

	itmp = F[i];
	p[itmp] = xf_raw[i] - x[itmp];

	if (xf_raw[i] < 0) { /* An infeasible value: compute alpha. */
	  
	  assert(fabs(p[itmp]) > 1e-12);
	  alpha[sizeAlpha++] = -x[itmp]/p[itmp]; 

	  // This should be the stepsize which makes x(n) + alpha p = 0 in 
	  // the itmp coordinate.

	  // I'm going by the sbls2 code here. So sue me.
	  // Ok, I don't get why this isn't x[i], but Adlers seems clear
	  //   alpha[sizeAlpha++] = xf_raw[i]/p[itmp].
	 	 
	}	    

      }

      alpha[sizeAlpha++] = 0; /* A step of size 0 is possible. */


      /* ******************************************************* */
      /* We know the alpha_i values, so determine biggest step   */
      /* in the direction p which reduces the value of q.  We    */
      /* skip anything bigger than one since that will be picked */
      /* up by the first time through the loop.                  */
      /* ******************************************************* */

      /* we'll reuse the value of q(x), so we might as well hold on to it */

      qofx = q(x,AprimeDotA,Apb,b,m,n);

      if (gVERBOSITY > 5) {printf("qofx to beat: %f\n",qofx);}

      // We now assemble xplusalphap and P[x + 1p]. 

      for(i=0; i<n; i++) {
	xplusalphap[i] = x[i] + p[i];
	Pxplusalphap[i] = max(xplusalphap[i],0);
      }

      double ntest=0.0;
      for(i=0; i<n; i++) {ntest += pow(Pxplusalphap[i]-xplusalphap[i],2.0); }
      
      ntest = 0.0;
      for(i=0; i<n; i++) {ntest += pow(xplusalphap[i],2.0); }

      newq = q(Pxplusalphap,AprimeDotA,Apb,b,m,n);
      if (gVERBOSITY > 5) { printf("q(P[x + 1p]): %g.\n",newq); }

      /* if we have improvement in q, we need not do the following */
      alphaItr = 0;

      if(newq > qofx){

	if (gVERBOSITY >= 10) { printf("tsnnls_spiv: \t Calling qsort.\n"); }
	/* darn, that didn't work, so we need to sort the alpha's */
	qsort(alpha,sizeAlpha,sizeof(taucs_double),compare_taucs_doubles);

	/* burn the ones where alpha >= 1 */
	alphaItr = 0;

	while(alpha[alphaItr] >= 1.0 && alphaItr < sizeAlpha){ alphaItr++; }
	assert(alphaItr < sizeAlpha);

	/* burn ending zeros, so only the last one is zero */
	while(alpha[sizeAlpha-2] == 0) { sizeAlpha--; }

	/* now we see if the step of size alpha[i] improves q */

	for(; alphaItr<sizeAlpha+8 && newq >= qofx; alphaItr++){

	  if (alphaItr < sizeAlpha-1) {

	    tmp = alpha[alphaItr];

	  } else { /* The smallest alpha doesn't work, so try shrinking further. */

	    tmp = alpha[sizeAlpha-2] * pow(0.5,alphaItr-(sizeAlpha-2));

	  }

	  if (tmp < 0) {

	    gErrorCode = 1199;
	    sprintf(gErrorString,
		    "tsnnls_spiv: lowest alpha appears to be < 0.\n");
	    return NULL;

	  }

	  for(i=0; i<n; i++){

	    xplusalphap[i] = x[i] + tmp*p[i];
	    Pxplusalphap[i] = max(xplusalphap[i],0);
	    
	  }

	  newq = q(Pxplusalphap,AprimeDotA,Apb,b,m,n);
	  if (gVERBOSITY > 5) { printf("q(P[x + (%g)p]: %g.\n",tmp,newq); }

	}

	if (alphaItr > sizeAlpha+8) { 

	  gErrorCode = 230;
	  sprintf(gErrorString,
		  "tsnnls_spiv: Reducing stepsize to %g did not \n"
		  "        produce a local reduction in qofx.\n",
		  tmp);
	  return NULL;

	}
	  
      }// end (qofxplusalphap > qofx)

      /* ************************************************* */
      /* well, we got a reduction in q, so now we update x */
      /* for the next (potential) round                    */
      /* ************************************************* */

      /* Note that we make x the (feasible) P[x + alphap], not the */
      /* potentially infeasible x + alphap. */

      memcpy(x,Pxplusalphap,sizeof(taucs_double)*n);
      
      /* ************************************************* */
      /* now we see if any of the frees want to be bounded */
      /* ************************************************* */

      if (gVERBOSITY >= 10) { 

	/* int vcnt; */

	printf("tsnnls_spiv: \t Checking infeasibles.\n"); 
	printf("tsnnls_spiv: Dumping F of size %d.\n",sizeF);

	/*for(vcnt = 0;vcnt < sizeF;vcnt++) {

	  printf("%d ",F[vcnt]);

	}
	
	printf("\n"); */

      }

      infeasible(F,Pxplusalphap,sizeF,H1,&sizeH1);
      
      if(sizeH1 == 0){

	fflag = 0;

	// We are going to leave the loop, so squash xf_raw.
	free(xf_raw);
	xf_raw = NULL;

	if (gVERBOSITY >= 10) { printf("tsnnls_spiv: H is empty.\n"); }

	break;
      }
      else{

	if (gVERBOSITY >= 10) { printf("tsnnls_spiv: H is nonempty.\n"); }

	int_difference(F,sizeF,H1,sizeH1,&sizeF);
	int_union(G,sizeG,H1,sizeH1,&sizeG);

      }

      // The first thing we're going to do up top is allocate a new xf_raw.
      // So we free the old one here.

      free(xf_raw);
      xf_raw = NULL;

    } // end inner fflag loop

    /* At this point, we should be at a stationary point for x_F, having bound
       a bunch of formerly free variables. We recompute and print qofx. */

    qofx = q(x,AprimeDotA,Apb,b,m,n);
    if (gVERBOSITY > 5) {printf("qofx at stationary point: %g.\n",qofx);}
   
    if (!(qofx < last_stationary_q)) { 

      if (inPrintErrorWarnings) {
	
	printf("tsnnls_spiv: qofx at stationary point %16g.\n"
	       "        qofx at last stationary  %16g.\n"
	       "\n"
	       "tsnnls_spiv: Aborting run.\n",
	       qofx,last_stationary_q);

      }

      gErrorCode = 295;
      sprintf(gErrorString,
	      "tsnnls_spiv: Error! qofx has not decreased between stationary points.\n"
	      "        (%16g -> %16g). Aborting run.\n",
	      last_stationary_q,qofx);
      
      return NULL;

    }

    last_stationary_q = qofx;

    /* We left the inner loop because sizeH was zero. This means that F has 
       not changed during this iteration. We confirm that the equation

       A'_F A_F x_F - A'_F b + A'_F A_B x_B = 0.

       Since the bound guys are bound to zero, x_B is zero, and we need only
       check that 

       A'_F A_F x_F - A'_F b = 0. 

       Because we have left the loop at a point where we didn't modify F
       in the last round, lsqrApA and Af are A'_FA_F and A_F respectively. 

       The vector x_F must be assembled from F and x, though. */

    /* double normgF = 0.0; */
    double *AfpAfxf, *Afpb, *xf, *gf;
    int    gfItr;

    AfpAfxf = calloc(A_original_ordering->m,sizeof(double));
    Afpb    = calloc(A_original_ordering->m,sizeof(double));
    xf      = calloc(sizeF,sizeof(double));
    gf      = calloc(sizeF,sizeof(double));

    for(gfItr=0;gfItr<sizeF;gfItr++) { xf[gfItr] = x[F[gfItr]]; }

    ourtaucs_ccs_times_vec(lsqrApA,xf,AfpAfxf);
    taucs_transpose_vec_times_matrix(b,A_original_ordering,F,sizeF,Afpb);

    for(gfItr=0;gfItr<sizeF;gfItr++) { gf[gfItr] = AfpAfxf[gfItr] - Afpb[gfItr]; }

    /* At this point, if we're really at a stationary point, this should be 
       pretty much dead zero. We check this to make sure. */

    for(gfItr=0;gfItr<sizeF;gfItr++) { 

      if (fabs(gf[gfItr]) > 1e-8) {

	gErrorCode = 74;
	sprintf(gErrorString,
		"tsnnls_spiv: Warning! Component %d of the reduced gradient is %g at the\n"
		"        end of the f loop. This suggests that we are not at a stationary\n"
		"        point and that something has gone wrong with the run.\n",
		gfItr,gf[gfItr]);

	return NULL;

      }

    }

    free(gf); free(xf); free(AfpAfxf); free(Afpb);
    
    /* Now we need to compute y_G and see if shifting anything out
       of F has created infeasibles in G */

    if (gVERBOSITY >= 10) { printf("tsnnls_spiv: F loop terminated. Computing y_g.\n");}

    // Reading the algorithm closely, it seems that the ENTIRE residual,
    // and not the constrained residual, is what's wanted here. We try it...

    // taucs_ccs_submatrix(A_original_ordering, F, sizeF, Af);

    // Note: it might be simpler to allocate residual up front
    // then we'd zero it out here
    if(sizeF != 0){
      
      // update xf_raw to include alpha*p. Since we are copying
      // into a new xf_raw, we reallocate it here.
      
      //xf_raw = calloc(sizeF,sizeof(double));    
      //for(i=0; i<sizeF; i++){ xf_raw[i] = x[F[i]]; }      
      
      /* Now compute the residual A_F x_F - b. This is an m-vector. */
      assert(residual == NULL);      
      residual = (taucs_double *)calloc(m,sizeof(taucs_double));
      //ourtaucs_ccs_times_vec(Af,xf_raw,residual);
      ourtaucs_ccs_times_vec(A_original_ordering,x,residual);

      //free(xf_raw);  // we won't use it again below
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

    DAXPY_F77(&m,&minusOne,b,&incX,residual,&incY);

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
    
    if (gVERBOSITY >= 10) { printf("tsnnls_spiv: Computing residual.\n"); }

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

    /* here we are setting up y. We only need to zero only the guys not in G */
    /* but it's safer to zero everything. */

    for(i=0; i<n; i++){ y[i] = 0.0; }
    for(i=0; i<sizeG; i++){ y[G[i]] = yg_raw[i]; }  

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

  if (gVERBOSITY >= 10) { printf("tsnnls_spiv: G loop terminated.\n"); }

  if( lsqrStep != 0 && pivcount < MAXPIVOT)
    {

      if (gVERBOSITY >= 10) { printf("tsnnls_spiv: Doing lsqr step.\n"); }

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

      if (gVERBOSITY >= 10) { printf("tsnnls_spiv: Survived lsqr.\n"); }

    }
  
  if( outResidualNorm != NULL && pivcount < MAXPIVOT)
    {

      if (gVERBOSITY >= 10) { printf("tsnnls_spiv: Computing final residual.\n"); }

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
  free(Pxplusalphap);

  free(p);
  free(alpha);

  if (gVERBOSITY >= 10) { printf("tsnnls_spiv: Done.\n"); }

  if (pivcount < MAXPIVOT) {

    return x;

  } else {

    gErrorCode = 999;
    sprintf(gErrorString,"tsnnls_spiv: Too many pivots (%d).",MAXPIVOT);

    return NULL;
  }

}
