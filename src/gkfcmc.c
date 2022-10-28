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
/* File: gkfcmc.c                                                         */
/*                                                                        */
/* History:                                                               */
/*          01-Apr-04 : Support for MATLAB V6.1 (R12.1) and above (V2-1)  */
/*          05-May-95 : Created (V2-0)                                    */
/*                                                                        */
/* Author(s): Bogdan R. Kosanovic                                         */
/*                                                                        */
/* Description: Gustafson-Kessel fuzzy c-means clustering algorithm.      */
/*                                                                        */
/*              The calling syntax is:                                    */
/*                                                                        */
/*              [xc,U,fineps,err,xctraj] = gkfcmc (x, U0, m, maxiter, ... */
/*                    maxeps, volc) ;                                     */
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
/*                volc    : volume constraints (real vector of length C). */
/*                          Contains C positive real numbers that         */
/*                          constrain the determinants of matrices Aj     */
/*                          which are shaping the resulting clusters.     */
/*                          optional, default=ones(C,1)                   */
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
#define volc_IN     pin[5]  /* optional */

/* Output Arguments */

#define xc_OUT      pout[0]
#define U_OUT       pout[1]
#define fineps_OUT  pout[2]  /* optional */
#define err_OUT     pout[3]  /* optional */
#define xctraj_OUT  pout[4]  /* optional */

/* default values */
#define DEF_M       2.0      /* fuzzy exponent */
#define DEF_MAXITER 100      /* maximum number of iterations */
#define DEF_MAXEPS  1e-3     /* maximum norm error allowed */
#define DEF_VOLC    1.0      /* volume constraints */

/* Matrix addressing (column-by-column) 2-D and 3-D */
#define mind2(i,j,ilen)  ((i)+(j)*(ilen))
#define mind3(i,j,jlen,k,klen)  mind2(j,mind2(k,i,klen),jlen)

/* Min and Max macros */
#define Min(x,y)    (((x)<(y))?(x):(y))
#define Max(x,y)    (((x)>(y))?(x):(y))

/* Static variables make this code non-reentrant! */

/* Pointers to distance function and matrix */
static double (*distfun_p)(const double*,const double*,const double*);

static double *dmatrix;
static int spacedim;      /* dimension of vector space */
/* a macro to invoke distance function */
#define CLTDistance(x,y,i)  \
        (*distfun_p)((x),(y),&dmatrix[mind3(i,0,spacedim,0,spacedim)])


/*-------------------------------------------------------------------------
  dist_mahal : mahalanobis distance function (squared!)
  RETURN : Distance between x and y, given IP matrix
-------------------------------------------------------------------------*/
double dist_mahal (
  const double *x,      /* first point  */
  const double *y,      /* second point */
  const double *A       /* distance matrix */
)
{
  register int i, j;
  register double dumm, dist;

  dist = 0.0;
  for (i = 0; i < spacedim; i++) {
    dumm = x[i] - y[i];
    for (j = 0; j < spacedim; j++)
      dist += dumm * A[mind2(i,j,spacedim)] * (x[j] - y[j]);
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
  clt_GetFuzzyScatter - calculate fuzzy scatter and distance matrices
  RETURN: Status. (0 - OK, otherwise - singular matrix)
  Description: Estimates fuzzy scatter matrices, their determinants,
     inverse, and distance matrices.
----------------------------------------------------------------------*/
static int clt_GetFuzzyScatter (
  const int C,           /* number of clusters (in) */
  const long N,          /* number of input patterns (in) */
  const int d,           /* number of features (dimensions) (in) */
  const double *x,       /* input patterns (in) */
  const double m,        /* fuzzy exponent (in) */
  const double *U,       /* partition matrix (in) */
  const double *xc,      /* cluster centers (in) */
  const double *volc,    /* volume constraints (in) */
  mxArray *tmpSf_stc     /* temporary scatter matrix */
)
/* NOTE: Will NOT check input arguments !!!!!!!!!!!!!!!!! */
{

  register int i, j1, j2;
  register long k;

  double dumm,
         det;    /* determinant */

  double *Ai,    /* i-th distance matrix */
         *Sf,    /* scatter matrix */
         *invSf; /* inverse scatter matrix */

  mxArray *det_stc,    /* determinant structure */
          *invSf_stc;  /* inverse scatter matrix structure */

  Sf = mxGetPr (tmpSf_stc);  /* pointer to scatter matrix */
  for (i = 0; i < C; i++) {   /* for each cluster */
    Ai = &dmatrix[mind3(i,0,d,0,d)];   /* locate the distance matrix */

    /* initialize fuzzy scatter matrix */
    for (j2 = 0; j2 < d; j2++)
      for (j1 = 0; j1 < d; j1++)
        Sf[mind2(j1,j2,d)] = 0;

    /* estimate it */
    for (k = 0L; k < N; k++) {
      dumm = pow (U[mind2(i,k,C)], m);
      for (j2 = 0; j2 < d; j2++)
        for (j1 = 0; j1 < d; j1++)
          Sf[mind2(j1,j2,d)] += dumm * (x[mind2(j1,k,d)] - xc[mind2(j1,i,d)]) *
                                       (x[mind2(j2,k,d)] - xc[mind2(j2,i,d)]);
    }

    /* get determinant */
    mexCallMATLAB(1, &det_stc, 1, &tmpSf_stc, "det");
    det = mxGetScalar(det_stc);
    mxDestroyArray(det_stc);

    if (det == 0) return (1);   /* singular matrix encountered */

    /* calculate inverse */
    mexCallMATLAB(1, &invSf_stc, 1, &tmpSf_stc, "inv");
    invSf = mxGetPr(invSf_stc);

    /* get distance matrix */
    dumm = pow (volc[i] * det, 1.0/d);
    for (j2 = 0; j2 < d; j2++)
      for (j1 = 0; j1 < d; j1++)
        Ai[mind2(j1,j2,d)] = dumm * invSf[mind2(j1,j2,d)];
    mxDestroyArray(invSf_stc);
  }
  return(0);
} /* clt_GetFuzzyScatter */


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
      Dik2 = CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,i,d)], i);
      if (Dik2 != 0) {  /* normal case (not a singularity) */
        for (j = 0, sum = 0; j < C; j++) {
          if (i == j)
            sum += 1.0;
          else {
            Djk2 = CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)], j);
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
          if (CLTDistance(&x[mind2(0,k,d)],&xc[mind2(0,j,d)], j) == 0) sing++;
        /* Set the rest of memberships to 0 or 1/sing */
        NewValue = 1.0/sing;
        if (fabs(NewValue - U[mind2(i,k,C)]) > maxn)  /* update maxnorm */
          maxn = fabs(NewValue - U[mind2(i,k,C)]);  /* only for i-th cluster */
        U[mind2(i,k,C)] = NewValue;
        for (j = i+1; j < C; j++) {
          if (CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)], j) == 0)
            U[mind2(j,k,C)] = NewValue;
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
                  * CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,i,d)], i);

  *maxnorm = maxn; *changes = chgs;
  return (sq_error);
} /* clt_GetFuzzyMemb */


/*-------------------------------------------------------------------------
  do_gkfcm : Gustafson-Kessel fuzzy c-means algorithm
  RETURN : Status ( 0: OK)
                  ( 1: Did not converge)
-------------------------------------------------------------------------*/
int do_gkfcm (
  const int C,          /* number of clusters (in) */
  const long N,         /* number of input patterns (in) */
  const int d,          /* number of features (dimensions) (in) */
  const double *x,      /* input pattern matrix (in) */
  const double m,       /* fuzzy exponent (in) */
  const double *volc,   /* volume constraints vector (in) */
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
  int stat;      /* status for clg_GetFuzzyScatter, etc. */
  long changes;  /* number of changed elements in U */

  double sq_error, /* objective function value */
         maxnorm;  /* max|U(i,k)[iter] - U(i,k)[iter+1]|, i.e. max change */

  mxArray *tmpdxd_stc; /* temporary scatter matrix (d x d) */

  /* Input arguments are NOT checked !!!!!!!!!!!!!!!!!! */

  /*----------------------------------------------------------------------*/
  /* Initial partition matrix is in U. Go through input samples and       */
  /* adjust centers until error is minimized and/or partition matrix does */
  /* not change                                                           */
  /*----------------------------------------------------------------------*/

  /* allocate temporary fuzzy scatter matrix */
  tmpdxd_stc = mxCreateDoubleMatrix(d, d, mxREAL);

  iter = 0; stat = 0;
  do {
    clt_GetFuzzyMeans(C, N, d, x, m, U, xc);  /* get new centers */
    if (xctraj != NULL)  /* save centers into trajectory matrix */
      for (j = 0; j < C; j++)
        for (i = 0; i < d; i++)
          xctraj[mind2(i,iter*C+j,d)] = xc[mind2(i,j,d)];

    /* get new scatter and distance matrices */
    stat = clt_GetFuzzyScatter(C, N, d, x, m, U, xc, volc, tmpdxd_stc);

    /* get new partition matrix */
    if (stat == 0) {   /* if fuzzy scatter matrices are non-singular */
      sq_error = clt_GetFuzzyMemb (C, N, d, x, m, U, xc, &maxnorm, &changes);
      if (err != NULL && changes > 0) err[iter] = sq_error;
      iter++;
    }
    else
      changes = 1;  /* to make sure *niter == iter */

    /* Check for termination (no changes? number of iterations?) */
  } while (stat == 0 && changes > 0 && iter < maxiter && maxnorm > eps);
  if (changes > 0)
    *niter = iter;     /* return number of iterations */
  else
    *niter = iter-1;   /* the last iteration did not produce any change */
  mxDestroyArray(tmpdxd_stc);  /* free tmp scatter matrix */
  if (fineps != NULL) *fineps = maxnorm; /* final maximum norm distance */
  return ((stat != 0 || (changes > 0 && iter >= maxiter)) ? 1 : 0);
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
  mexFunction : entry point for gkfcmc()
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
      niter;   /* acutal number of iterations */

  long n,     /* number of patterns */
       il, jl, ldumm;

  double *x,  /* input pattern matrix */
         *U0, /* initial partition matrix */
         *xc, /* cluster centers */
         *U,  /* final partition matrix */
         *fineps_p, /* pointer to final maximum norm error */
         *err_big,    /* objective function values (internal vector) */
         *err,
         *xctraj_big, /* center trajectories (internal vector) */
         *xctraj,
         m,  /* fuzzy exponent */
         maxeps,  /* maximum norm allowed */
         dumm,
         *volc;    /* volume constraints */

  mxArray *volc_stc,      /* volume constraints vector structure */
          *tmpdxd_stc[2]; /* temporary matrices (d x d) */

  /*----------------------------------*/
  /* Check the input/output arguments */
  /*----------------------------------*/

  if (nin < 2) mexErrMsgTxt ("insufficient number of input arguments");
  if (nin > 6) mexErrMsgTxt ("too many input arguments");
  if (nout > 5) mexErrMsgTxt ("too many output arguments");

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

  /* volc: volume constraints (optional, default=ones(C,1) */
  if (nin < 6) { /* volc not present */
    volc_stc = mxCreateDoubleMatrix (C, 1, mxREAL);
    volc = mxGetPr(volc_stc);
    for (i = 0; i < C; i++) volc[i] = DEF_VOLC;   /* default value */    
  }
  else {  /* volc is present */
    idumm = mxGetM(volc_IN); idumm2 = mxGetN(volc_IN);
    if (Min(idumm,idumm2) != 1)
      mexErrMsgTxt ("volume constraints should be a (non-empty) vector");
    if (Max(idumm,idumm2) != C)
      mexErrMsgTxt ("invalid number of volume constraints");

    /* check data type for volume constraints */
    if (!mxIsNumeric(volc_IN) || !mxIsDouble(volc_IN))
      mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
    if (mxIsComplex(volc_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
    if (mxIsSparse(volc_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);
    volc = mxGetPr(volc_IN);
  }

  /* check if all volume constraints are positive */
  for (i = 0; i < C; i++)
    if (volc[i] <= 0) mexErrMsgTxt ("all volume constraints must be positive");

  /* fetch the pointers to required vector data */
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
  dmatrix = (double*)mxCalloc(C*d*d, sizeof(double));

  /*---------------------------*/
  /* Prepare for gkfcm routine */
  /*---------------------------*/

  /* Select distance function */
  spacedim = d;       /* vector space dimension */
  distfun_p = dist_mahal;   /* use Mahalanobis distance */

  /*---------------*/
  /* execute gkfcm */
  /*---------------*/

  stat = do_gkfcm (C, n, d, x, m, volc, U, xc, maxiter, maxeps, &niter,
                   fineps_p, err_big, xctraj_big);

  /* check output results */
  if (stat == 1)
    mexPrintf ("WARNING: did not converge after %d iterations\n", maxiter);
  mxFree(dmatrix);    /* free distance matrices */

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
    }
  }

} /* mexFunction */
