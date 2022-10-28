/*--------------------------------------------------------------------------------*/
/* MIT License                                                                    */
/*                                                                                */
/* Copyright (c) 1995-2022 Bogdan Kosanovic                                       */
/*                                                                                */
/* Permission is hereby granted, free of charge, to any person obtaining a copy   */
/* of this software and associated documentation files (the "Software"), to deal  */
/* in the Software without restriction, including without limitation the rights   */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      */
/* copies of the Software, and to permit persons to whom the Software is          */
/* furnished to do so, subject to the following conditions:                       */
/*                                                                                */
/* The above copyright notice and this permission notice shall be included in all */
/* copies or substantial portions of the Software.                                */
/*                                                                                */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  */
/* SOFTWARE.                                                                      */
/*--------------------------------------------------------------------------------*/

/*========================================================================*/
/* File: fcmc.c                                                           */
/*                                                                        */
/* History:                                                               */
/*          01-Apr-04 : Support for MATLAB V6.1 (R12.1) and above (V2-1)  */
/*          06-Apr-95 : Bug in calculating memberships for singularities  */
/*                      is fixed. (V1-5)                                  */
/*          20-Mar-95 : Return the matrix which defines metric (V1-4)     */
/*          19-Nov-94 : Name changed from fcm to fcmc (V1-1)              */
/*                      Spurious 'if' checking type of U0_IN removed      */
/*                      Singularity will be counted as a "change"         */
/*          12-Nov-94 : Created (V1-0)                                    */
/*                                                                        */
/* Author(s): Bogdan R. Kosanovic                                         */
/*                                                                        */
/* Description: Fuzzy c-means clustering algorithm.                       */
/*                                                                        */
/*              The calling syntax is:                                    */
/*                                                                        */
/*              [xc,U,fineps,err,xctraj,A] = fcmc (x, U0, m, maxiter, ... */
/*                    maxeps, dflag);                                     */
/*                                                                        */
/*              Input:                                                    */
/*                                                                        */
/*                x       : pattern matrix (d x n real matrix)            */
/*                          d > 0 (number of dimensions)                  */
/*                          n > 1 (number of patterns)                    */
/*                U0      : initial partition matrix (C x n real matrix)  */
/*                          C > 1 (number of clusters)                    */
/*                          C <= n (better be < n)                        */
/*                m       : fuzzy exponent (real > 1), optional,          */
/*                          default=2.00                                  */
/*                maxiter : maximum number of iterations (> 0), optional  */
/*                          default=100                                   */
/*                maxeps  : maximum norm error allowed (> 0), optional    */
/*                          default=1e-3                                  */
/*                          used as a stopping criterion                  */
/*                          suggested range (1e-3 to 5e-4)                */
/*                dflag   : distance flag (0: Euclidean, 1: diagonal,     */
/*                          2: Mahalanobis) optional, default=0           */
/*                                                                        */
/*              Output:                                                   */
/*                                                                        */
/*                xc     : cluster centers (d x C real matrix)            */
/*                         each column stands for a center                */
/*                U      : final partition matrix (C x n real matrix)     */
/*                         each row represents a membership function      */
/*                fineps : final maximum norm error (scalar)              */
/*                err    : objective function values (niter x 1 real      */
/*                         vector), where niter is the final number       */
/*                         of iterations                                  */
/*                xctraj : center trajectories (d x C*niter real matrix)  */
/*                         composed of niter submatrices (d x C) for each */
/*                         iteration. Submatrices have same structure as  */
/*                         xc matrix. The last submatrix equals xc.       */
/*                A      : matrix that was used to calculate distance     */
/*                         metric.                                        */
/*                                                                        */
/* Note: Since the FCM acronym may be mistaken with the Fuzzy Cognitive   */
/*       Maps, the external interface is changed to use the FCMC acronym. */
/*                                                                        */
/*========================================================================*/

#include <math.h>
#include "mex.h"

/* #define NULL  0 this was old version of MATLAB */

/* Input Arguments */

#define x_IN        pin[0]
#define U0_IN       pin[1]
#define m_IN        pin[2]  /* optional */
#define maxiter_IN  pin[3]  /* optional */
#define maxeps_IN   pin[4]  /* optional */
#define dflag_IN    pin[5]  /* optional */

/* Output Arguments */

#define xc_OUT      pout[0]
#define U_OUT       pout[1]
#define fineps_OUT  pout[2]  /* optional */
#define err_OUT     pout[3]  /* optional */
#define xctraj_OUT  pout[4]  /* optional */
#define A_OUT       pout[5]  /* optional */

/* distance types */
#define DIST_IP_EUCLIDEAN   0
#define DIST_IP_DIAGONAL    1
#define DIST_IP_MAHALANOBIS 2

/* default values */
#define DEF_M       2.0      /* fuzzy exponent */
#define DEF_MAXITER 100      /* maximum number of iterations */
#define DEF_MAXEPS  1e-3     /* maximum norm error allowed */
#define DEF_DFLAG   DIST_IP_EUCLIDEAN   /* distance type */

/* Matrix addressing (column-by-column) */
#define mind2(i,j,ilen)  ((i)+(j)*(ilen))

/* Static variables make this code non-reentrant! */

/* Pointers to distance function and matrix */
static double (*distfun_p)(const double*,const double*);

static double *dmatrix;
static int spacedim;      /* dimension of vector space */
/* a macro to invoke distance function */
#define CLTDistance(x,y)  (*distfun_p)((x),(y))



/*-------------------------------------------------------------------------
  dist_euclid : euclidean distance function (squared!)
  RETURN : Distance between x and y
-------------------------------------------------------------------------*/
double dist_euclid (
  const double *x,      /* first point  */
  const double *y       /* second point */
)
{
  register int i;
  register double dumm, dist;

  dist = 0.0;
  for (i = 0; i < spacedim; i++) {
    dumm = x[i] - y[i];
    dist += dumm * dumm;
  }
  return (dist);
} /* dist_euclid */

/*-------------------------------------------------------------------------
  dist_diag : diagonal distance function (squared!)
  RETURN : Distance between x and y
-------------------------------------------------------------------------*/
double dist_diag (
  const double *x,      /* first point  */
  const double *y       /* second point */
)
{
  register int i;
  register double dumm, dist;

  dist = 0.0;
  for (i = 0; i < spacedim; i++) {
    dumm = x[i] - y[i];
    dist += dumm * dumm * dmatrix[i];
  }
  return (dist);
} /* dist_diag */

/*-------------------------------------------------------------------------
  dist_mahal : mahalanobis distance function (squared!)
  RETURN : Distance between x and y
-------------------------------------------------------------------------*/
double dist_mahal (
  const double *x,      /* first point  */
  const double *y       /* second point */
)
{
  register int i, j;
  register double dumm, dist;

  dist = 0.0;
  for (i = 0; i < spacedim; i++) {
    dumm = x[i] - y[i];
    for (j = 0; j < spacedim; j++)
      dist += dumm * dmatrix[mind2(i,j,spacedim)] * (x[j] - y[j]);
  }
  return (dist);
} /* dist_mahal */




/*----------------------------------------------------------------------
  clt_GetFuzzyMeans - calculate cluster centers using fuzzy memberships
  Description: Goes through all input samples and calculates
     cluster centers according to pattern membership values (matrix U).
----------------------------------------------------------------------*/
static void clt_GetFuzzyMeans (
  const int C,           /* number of clusters (in) */
  const long N,          /* number of input patterns (in) */
  const int d,           /* number of features (dimensions) (in) */
  const double *x,       /* input patterns (in) */
  const double m,        /* fuzzy exponent (in) */
  const double *U,       /* partition matrix (in) */
  double *xc             /* cluster centers (out) */
)
/* NOTE: Will NOT check input arguments !!!!!!!!!!!!!!!!! */
{
  register int i, j;
  register long k;
  register double wik,  /* Uik^m weight for k-th pattern within i-th cluster */
                  wsum;     /* sum of weights for i-th cluster */
 
  /* Go through all clusters */
  for (i = 0; i < C; i++) {
    /* initialize center information (note: xc is d x C !!!) */
    for (j = 0; j < d; j++) xc[mind2(j,i,d)] = 0.0;
    /* go through all patterns */
    for (k = 0, wsum = 0; k < N; k++) {
      wik = pow (U[mind2(i,k,C)], m);
      for (j = 0; j < d; j++) xc[mind2(j,i,d)] += wik * x[mind2(j,k,d)];
      wsum += wik;
    }
    for (j = 0; j < d; j++) xc[mind2(j,i,d)] /= wsum;   /* scale the sum */
  }
} /* clt_GetFuzzyMeans */




/*----------------------------------------------------------------------
  clt_GetFuzzyMemb - calculate fuzzy membership functions using centers
  RETURN: Squared error (value of fuzzy objective function)
  Description: Goes through all input samples and calculates membership
     functions according to the cluster centers (xc).
     Calculates the fuzzy objective function value given centers and U.
     distfun_p _must_ be selected. (in order to get squared distance)
     Also calculates the maximum norm distance relative to previous
     value of the membership functions, and the number of changes
     in the membership values.
----------------------------------------------------------------------*/
static double clt_GetFuzzyMemb (
  const int C,         /* number of clusters (in) */
  const long N,        /* number of input patterns (in) */
  const int d,         /* number of features (dimensions) (in) */
  const double *x,     /* input patterns (in) */
  const double m,      /* fuzzy exponent (in) */
  double *U,           /* partition matrix (in/out) */
  const double *xc,    /* cluster centers (in) */
  double *maxnorm,     /* max norm distance from the previous estimate (out) */
  long *changes        /* number of changes in the mambership values (out) */
)
/* NOTE: Will NOT check input arguments !!!!!!!!!!!!!!!!!! */
{
  int i, j,
      sing;    /* number of singularities */
  long k,
       chgs;       /* number of changes */
  double sum,
         NewValue,   /* new membership value */
         maxn,       /* maximum norm */
         Dik2, Djk2, /* distance squared (vi,xk) and (vj,xk) */
         sq_error;   /* fuzzy objective function value */
 
  maxn = 0;           /* assume no changes will occur */
  chgs = 0L;          /* count the number of changed elements in U */
  for (k = 0L; k < N; k++) {  /* take the next sample (pattern) */
    for (i = 0; i < C; i++) {  /* go through all clusters */
      Dik2 = CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,i,d)]);
      if (Dik2 != 0) {  /* normal case (not a singularity) */
        for (j = 0, sum = 0; j < C; j++) {
          if (i == j)
            sum += 1.0;
          else {
            Djk2 = CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)]);
            sum += pow(Dik2/Djk2,1/(m-1));
          }
        }
        NewValue = 1 / sum;
        if (NewValue != U[mind2(i,k,C)]) {
          if (fabs(NewValue - U[mind2(i,k,C)]) > maxn)  /* update maxnorm */
            maxn = fabs(NewValue - U[mind2(i,k,C)]);
          U[mind2(i,k,C)] = NewValue; chgs++; /* only in order to count them */
        }
      }
      else {    /* singularity (very unlikely case) */
        chgs++; /* since not likely, just mark it as a change, */
        /* all previous memberships should be forced to zero */
        for (j = 0; j < i; j++) U[mind2(j,k,C)] = 0.0;
        /* Count how many singularities for this sample */
        for (j = i+1, sing = 1; j < C; j++)
          if (CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)]) == 0) sing++;
        /* Set the rest of memberships to 0 or 1/sing */
        NewValue = 1.0/sing;
        if (fabs(NewValue - U[mind2(i,k,C)]) > maxn)  /* update maxnorm */
          maxn = fabs(NewValue - U[mind2(i,k,C)]);  /* only for i-th cluster */
        U[mind2(i,k,C)] = NewValue;
        for (j = i+1; j < C; j++) {
          if (CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)]) == 0)
            U[mind2(j,k,C)] = NewValue;        /* bug fix 4/6/95 */
          else
            U[mind2(j,k,C)] = 0.0;
        }
        break; /* no reasons to go further, we finished with a singular */
               /* pattern case */
      }
    }  /* for i through clusters */
  }  /* for k through patterns */

  /* Calculate error squared (note distfun_p _must_ be preselected) */
  for (k = 0L, sq_error = 0.0; k < N; k++)  /* next sample */
    for (i = 0; i < C; i++)
      sq_error += pow(U[mind2(i,k,C)], m)
                  * CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,i,d)]);

  *maxnorm = maxn; *changes = chgs;
  return (sq_error);
} /* clt_GetFuzzyMemb */




/*-------------------------------------------------------------------------
  do_fcm : fuzzy c-means algorithm
  RETURN : Status ( 0: OK)
                  ( 1: Did not converge)
-------------------------------------------------------------------------*/
int do_fcm (
  const int C,          /* number of clusters (in) */
  const long N,         /* number of input patterns (in) */
  const int d,          /* number of features (dimensions) (in) */
  const double *x,      /* input pattern matrix (in) */
  const double m,       /* fuzzy exponent (in) */
  double *U,            /* initial/final partition matrix (in/out) */
  double *xc,           /* final cluster centers' matrix (out) */
  const int maxiter,    /* maximum number of iterations (in) */
  const double eps,     /* maximum norm distance allowed (in) */
  int *niter,           /* final number of iterations (out) */
  double *fineps,       /* final maximum norm distance (out) */
  double *err,          /* vector of size maxiter, contains niter values */
                        /* of clustering error, for each iteration */
  double *xctraj        /* center trajectories (out) */
)
{
  register int i, j,
               iter;   /* number of iterations */

  long changes;  /* number of changed elements in U */

  double sq_error, /* objective function value */
         maxnorm;  /* max|U(i,k)[iter] - U(i,k)[iter+1]|, i.e. max change */

  /* Input arguments are NOT checked !!!!!!!!!!!!!!!!!! */

  /*----------------------------------------------------------------------*/
  /* Initial partition matrix is in U. Go through input samples and       */
  /* adjust centers until error is minimized and/or partition matrix does */
  /* not change                                                           */
  /*----------------------------------------------------------------------*/

  iter = 0;
  do {
    clt_GetFuzzyMeans(C, N, d, x, m, U, xc);  /* get new centers */
    if (xctraj != NULL)  /* save centers into trajectory matrix */
      for (j = 0; j < C; j++)
        for (i = 0; i < d; i++)
          xctraj[mind2(i,iter*C+j,d)] = xc[mind2(i,j,d)];

    /* get new partition matrix */
    sq_error = clt_GetFuzzyMemb (C, N, d, x, m, U, xc, &maxnorm, &changes);
    if (err != NULL && changes > 0) err[iter] = sq_error;
    iter++;
    /* Check for termination (no changes? number of iterations?) */
  } while (changes > 0 && iter < maxiter && maxnorm > eps);
  if (changes > 0)
    *niter = iter;     /* return number of iterations */
  else
    *niter = iter-1;   /* the last iteration did not produce any change */
  if (fineps != NULL) *fineps = maxnorm; /* final maximum norm distance */
  return ((changes > 0 && iter >= maxiter) ? 1 : 0);
} /* do_fcm */




/* frequent messages */
#define MSG_INVTYPINPARG 0
#define MSG_REALONLY     1
#define MSG_FULLONLY     2

static char *msgs[3] = {
  "invalid type of input arguments",   /* MSG_INVTYPINPARG */
  "works with real matrices only",     /* MSG_REALONLY */
  "works with full matrices only"      /* MSG_FULLONLY */
}; /* msgs */

/*-------------------------------------------------------------------------
  mexFunction : entry point for fcmc()
-------------------------------------------------------------------------*/
void mexFunction(
  int           nout,     /* number of output elements   */
  mxArray       *pout[],  /* pointers to output matrices */
  int           nin,      /* number of input elements    */
  const mxArray *pin[]    /* pointers to input matrices  */
)
{
  int d,      /* number of features (dimensions) */
      C,      /* number of clusters */
      idumm, idumm2, stat, i, j,
      maxiter, /* maximum number of iterations */
      niter,   /* acutal number of iterations */
      dflag;   /* distance type */

  long n,     /* number of patterns */
       il, jl, ldumm;

  double *x,  /* input pattern matrix */
         *U0, /* initial partition matrix */
         *xc, /* cluster centers */
         *A,  /* output distance matrix */
         *U,  /* final partition matrix */
         *fineps_p, /* pointer to final maximum norm error */
         *err_big,    /* objective function values (internal vector) */
         *err,
         *xctraj_big, /* center trajectories (internal vector) */
         *xctraj,
         m,  /* fuzzy exponent */
         maxeps,  /* maximum norm allowed */
         dumm,
         *xMeans,  /* means of input patterns */
         *xVars;   /* variances of input patterns */

  mxArray *dmatrix_stc[2];    /* distance matrix structure */

  /*----------------------------------*/
  /* Check the input/output arguments */
  /*----------------------------------*/

  if (nin < 2) mexErrMsgTxt ("insufficient number of input arguments");
  if (nin > 6) mexErrMsgTxt ("too many input arguments");
  if (nout > 6) mexErrMsgTxt ("too many output arguments");

  /* fetch input arguments */
  /* x matrix (d x n) each column is one pattern vector */
  d = mxGetM(x_IN); n = (long)mxGetN(x_IN);
  if (n < 2) mexErrMsgTxt ("should have at least two patterns");
  if (d < 1) mexErrMsgTxt ("invalid number of features (pattern matrix)");
  if (d > n) mexPrintf ("WARNING: more features than patterns\n");

  /* U0 matrix (C x n) each row is one membership function */
  C = mxGetM(U0_IN); ldumm = (long)mxGetN(U0_IN);
  if (C < 2) mexErrMsgTxt ("requested number of clusters less then two");
  if (ldumm != n) mexErrMsgTxt("invalid number of columns (partition matrix)");
  if (C > n) mexErrMsgTxt ("more clusters than patterns");
  if (C > (double)n/5/(d+1))
    mexPrintf ("WARNING: number of patterns may be too small\n");

  /* data types for required arguments */
  if (!mxIsNumeric(x_IN) || !mxIsDouble(x_IN) ||
      !mxIsNumeric(U0_IN) || !mxIsDouble(U0_IN))
    mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
  if (mxIsComplex(x_IN) ||
      mxIsComplex(U0_IN))
     mexErrMsgTxt (msgs[MSG_REALONLY]);
  if (mxIsSparse(x_IN) ||
      mxIsSparse(U0_IN))
     mexErrMsgTxt (msgs[MSG_FULLONLY]);

  /* m : fuzzy exponent (optional, default=2.00) */
  if (nin < 3)
    m = DEF_M;
  else {  /* m is present */
    idumm = mxGetM(m_IN); idumm2 = mxGetN(m_IN);
    if (idumm != 1 || idumm2 != 1)
      mexErrMsgTxt ("fuzzy exponent must be a scalar");
    /* check data type */
    if (!mxIsNumeric(m_IN) || !mxIsDouble(m_IN))
      mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
    if (mxIsComplex(m_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
    if (mxIsSparse(m_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);

    m = mxGetScalar (m_IN);
  }
  if (m <= 1) mexErrMsgTxt ("fuzzy exponent out of range");

  /* maxiter : maximum number of iterations (optional, default=100) */
  if (nin < 4)
    maxiter = DEF_MAXITER;
  else {  /* maxiter is present */
    idumm = mxGetM(maxiter_IN); idumm2 = mxGetN(maxiter_IN);
    if (idumm != 1 || idumm2 != 1)
      mexErrMsgTxt ("maximum number of iterations must be a scalar");
    /* check data type */
    if (!mxIsNumeric(maxiter_IN) || !mxIsDouble(maxiter_IN))
      mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
    if (mxIsComplex(maxiter_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
    if (mxIsSparse(maxiter_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);

    maxiter = (int)mxGetScalar (maxiter_IN);
  }
  if (maxiter < 1) mexErrMsgTxt ("maxiter < 1");

  /* maxeps : maximum norm error allowed (optional, default=1e-3) */
  if (nin < 5)
    maxeps = DEF_MAXEPS;
  else {  /* maxeps is present */
    idumm = mxGetM(maxeps_IN); idumm2 = mxGetN(maxeps_IN);
    if (idumm != 1 || idumm2 != 1)
      mexErrMsgTxt ("maximum norm error must be a scalar");
    /* check data type */
    if (!mxIsNumeric(maxeps_IN) || !mxIsDouble(maxeps_IN))
      mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
    if (mxIsComplex(maxeps_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
    if (mxIsSparse(maxeps_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);

    maxeps = mxGetScalar (maxeps_IN);
  }
  if (maxeps <= 0) mexErrMsgTxt ("maxeps <= 0");

  /* dflag : distance type (optional, default=DIST_IP_EUCLIDEAN (0) ) */
  if (nin < 6)
    dflag = DEF_DFLAG;
  else {  /* dflag is present */
    idumm = mxGetM(dflag_IN); idumm2 = mxGetN(dflag_IN);
    if (idumm != 1 || idumm2 != 1)
      mexErrMsgTxt ("distance flag must be a scalar");
    /* check data type */
    if (!mxIsNumeric(dflag_IN) || !mxIsDouble(dflag_IN))
      mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
    if (mxIsComplex(dflag_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
    if (mxIsSparse(dflag_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);

    dflag = (int)mxGetScalar (dflag_IN);
  }
  if (dflag < 0 || dflag > 2) mexErrMsgTxt ("invalid distance type");

  /* fetch the pointers to vector data */
  x = mxGetPr (x_IN);
  U0 = mxGetPr (U0_IN);

  /*----------------------*/
  /* Allocate some memory */
  /*----------------------*/

  /* using MATLAB matrix structures */
  /* xc matrix (d x C) each column is one cluster center */
  xc_OUT = mxCreateDoubleMatrix(d, C, mxREAL);
  xc = mxGetPr(xc_OUT);

  /* U matrix (C x n) each row is one membership function */
  U_OUT = mxCreateDoubleMatrix(C, n, mxREAL);
  U = mxGetPr(U_OUT);
  for (jl = 0L; jl < n*C; jl++) U[jl] = U0[jl]; /* initial partition matrix */

  /* final maximum norm error */
  fineps_p = (double*)NULL;
  if (nout > 2) {
    fineps_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    fineps_p = mxGetPr(fineps_OUT);
  }

  /* using MATLAB memory pool */
  err_big = xctraj_big = (double*)NULL;
  if (nout > 3) {
    err_big = (double*)mxCalloc (maxiter, sizeof(double));
    if (nout > 4)
      xctraj_big = (double*)mxCalloc (d*C*maxiter, sizeof(double));
  }

  /*-------------------------*/
  /* Prepare for fcm routine */
  /*-------------------------*/

  /* Select distance */
  spacedim = d;       /* vector space dimension */
  switch (dflag) {
  case DIST_IP_EUCLIDEAN:
    distfun_p = dist_euclid;
    dmatrix_stc[0] = dmatrix_stc[1] = (mxArray*) NULL;
    dmatrix = (double*) NULL;
    break;
  case DIST_IP_DIAGONAL:
  case DIST_IP_MAHALANOBIS:
    /* get mean */
    xMeans = (double *)mxCalloc (d, sizeof(double));
    for (i = 0; i < d; i++) {
      xMeans[i] = 0;
      for (jl = 0L; jl < n; jl++)
        xMeans[i] += x[mind2(i,jl,d)];
      xMeans[i] /= (double)n;
    }

    if (dflag == DIST_IP_DIAGONAL) {  /* diagonal */
      distfun_p = dist_diag;
      dmatrix_stc[0] = mxCreateDoubleMatrix (d, 1, mxREAL);
      dmatrix_stc[1] = (mxArray*) NULL;
      dmatrix = mxGetPr(dmatrix_stc[0]);

      /* variances and diagonal matrix */
      xVars = (double *)mxCalloc (d, sizeof(double));
      for (i = 0; i < d; i++) {
        xVars[i] = 0;
        for (jl = 0L; jl < n; jl++) {
          dumm = x[mind2(i,jl,d)] - xMeans[i];
          xVars[i] += dumm * dumm;
        }
        xVars[i] /= (double)n;
        dmatrix[i] = 1.0 / xVars[i];
      }
      mxFree(xVars);
    }
    else {  /* mahalanobis */
      distfun_p = dist_mahal;
      dmatrix_stc[1] = mxCreateDoubleMatrix (d, d, mxREAL);  /* scratch matrix */
      dmatrix = mxGetPr(dmatrix_stc[1]);

      /* create mahalanobis distance matrix */

      /* estimate autocovariance matrix */
      for (j = 0; j < d; j++)
        for (i = 0; i < d; i++) {
          for (jl = 0; jl < n; jl++)
            dmatrix[mind2(i,j,d)] += (x[mind2(i,jl,d)] - xMeans[i]) *
                                     (x[mind2(j,jl,d)] - xMeans[j]);
          dmatrix[mind2(i,j,d)] /= (double)n;
        }

      /* find the inverse */
      mexCallMATLAB(1, &dmatrix_stc[0], 1, &dmatrix_stc[1], "inv");
      dmatrix = mxGetPr(dmatrix_stc[0]);   /* inverse matrix */
      /* free the scratch matrix */
      mxDestroyArray (dmatrix_stc[1]);
    }
    mxFree(xMeans);
    break;
  default:
    mexErrMsgTxt ("internal error no. 1");
  }

  /*-------------*/
  /* execute fcm */
  /*-------------*/

  stat = do_fcm (C, n, d, x, m, U, xc, maxiter, maxeps, &niter, fineps_p,
                 err_big, xctraj_big);

  /* check output results */
  if (stat == 1)
    mexPrintf ("WARNING: did not converge after %d iterations\n", maxiter);

  /*-----------------------------------*/
  /* check and create output variables */
  /*-----------------------------------*/

  /* err vector (niter x 1) objective function values */
  if (nout > 3) {
    err_OUT = mxCreateDoubleMatrix(niter, 1, mxREAL);
    err = mxGetPr(err_OUT);
    for (i = 0; i < niter; i++) err[i] = err_big[i];
    mxFree(err_big);

    /* xctraj matrix (niter x 1) objective function values */
    if (nout > 4) {
      xctraj_OUT = mxCreateDoubleMatrix(d, C*niter, mxREAL);
      xctraj = mxGetPr(xctraj_OUT);
      for (il = 0L; il < (long)d*C*niter; il++) xctraj[il] = xctraj_big[il];
      mxFree(xctraj_big);

      /* distance matrix */
      if (nout > 5) {
        A_OUT = mxCreateDoubleMatrix(d, d, mxREAL);
        A = mxGetPr(A_OUT);
        switch (dflag) {
        case DIST_IP_EUCLIDEAN:
          for (j = 0; j < d; j++)
            for (i = 0; i < d; i++)
              A[mind2(i,j,d)] = (i==j?1:0);
          break;
        case DIST_IP_DIAGONAL:
          for (j = 0; j < d; j++)
            for (i = 0; i < d; i++)
              A[mind2(i,j,d)] = (i==j?dmatrix[i]:0);
          break;
        case DIST_IP_MAHALANOBIS:
          for (j = 0; j < d; j++)
            for (i = 0; i < d; i++)
              A[mind2(i,j,d)] = dmatrix[mind2(i,j,d)];
          break;
        default:
          mexErrMsgTxt ("internal error no. 2");
        }
      }
    }
  }

} /* mexFunction */
/* nothing past this point */
