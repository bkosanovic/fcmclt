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
/* File: extfpm.c                                                         */
/*                                                                        */
/* History:                                                               */
/*          01-Apr-04 : Support for MATLAB V6.1 (R12.1) and above (V2-1)  */
/*          20-Mar-95 : Distance matrix A can be specified as             */
/*                      input. Distance matrix that was used will be      */
/*                      returned in B. (V1-4)                             */
/*          11-Dec-94 : Comment on non-Euclidean metric (V1-2)            */
/*          19-Nov-94 : Created from fcmc.c (V1-1)                        */
/*                                                                        */
/* Author(s): Bogdan R. Kosanovic                                         */
/*                                                                        */
/* Description: Extrapolate fuzzy partition matrix.                       */
/*                                                                        */
/*              The calling syntax is:                                    */
/*                                                                        */
/*              [U, err, B] = extfpm (x, xc, m, dflag, A);                */
/*                                                                        */
/*              Input:                                                    */
/*                                                                        */
/*                x       : pattern matrix (d x n real matrix)            */
/*                          d > 0 (number of dimensions)                  */
/*                          n > 1 (number of patterns)                    */
/*                xc      : presumed cluster centers (d x C real matrix)  */
/*                          each column stands for a center,              */
/*                          C > 1 (number of clusters), C < n             */
/*                m       : fuzzy exponent (real > 1), optional,          */
/*                          default=2.00                                  */
/*                dflag   : distance flag (0: Euclidean, 1: diagonal,     */
/*                          2: Mahalanobis), optional, default=0          */
/*                A       : distance matrix (default=[as estimated from   */
/*                          input data in x]) Ignored for dflag=0,        */
/*                          must be diagonal for dflag=1.                 */
/*                                                                        */
/*              Output:                                                   */
/*                                                                        */
/*                U       : partition matrix (C x n real matrix) for      */
/*                          cluster centers, C > 1 (number of clusters)   */
/*                          C <= n (better be < n), each row represents   */
/*                          a membership function                         */
/*                err     : objective function value (real scalar)        */
/*                B       : distance matrix used in calculation of U.     */
/*                          Should equal A, if A was specified as input.  */
/*                                                                        */
/* WARNING:  In case of non-Euclidean distance, it may not make much      */
/* sense to extrapolate large quantities of data, since the statistics    */
/* may be rather different from the one that was used in estimation of    */
/* the cluster prototypes. Namely, the new metric may not correspond to   */
/* the metric used in estimation.                                         */
/*                                                                        */
/* This problem may be reduced by not calculating the matrix A that       */
/* modifies the metric for new data, but providing the original matrix    */
/* used in estimation of the prototypes as returned by FCMC routine.      */
/* Then, the extrapolation will be with respect to the original metric    */
/* used within the FCMC routine.                                          */
/*                                                                        */
/*========================================================================*/

#include <math.h>
#include "mex.h"

/* #define NULL  0 this was old version of MATLAB */

/* Input Arguments */

#define x_IN        pin[0]
#define xc_IN       pin[1]
#define m_IN        pin[2]  /* optional */
#define dflag_IN    pin[3]  /* optional */
#define A_IN        pin[4]  /* optional */

/* Output Arguments */

#define U_OUT       pout[0]
#define err_OUT     pout[1] /* optional */
#define B_OUT       pout[2] /* optional */

/* distance types */
#define DIST_IP_EUCLIDEAN   0
#define DIST_IP_DIAGONAL    1
#define DIST_IP_MAHALANOBIS 2

/* default values */
#define DEF_M       2.0      /* fuzzy exponent */
#define DEF_DFLAG   DIST_IP_EUCLIDEAN   /* distance type */

/* Matrix addressing (column-by-column) */
#define mind2(i,j,ilen)  ((i)+(j)*(ilen))

/* Static variables make this code non-reentrant! */

/* Pointers to distance function and matrix */
static double (*distfun_p)(const double*,const double*);

static double *dmatrix;   /* matrix that modifies distance function */
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
  do_extfpm - calculate fuzzy partition maxtrix using centers
  RETURN: Squared error (value of fuzzy objective function)
  Description: Goes through all input samples and calculates membership
     functions according to the cluster centers (xc).
     Calculates the fuzzy objective function value given centers and U.
     distfun_p _must_ be selected. (in order to get squared distance)
----------------------------------------------------------------------*/
static double do_extfpm (
  const int C,         /* number of clusters (in) */
  const long N,        /* number of input patterns (in) */
  const int d,         /* number of features (dimensions) (in) */
  const double *x,     /* input patterns (in) */
  const double m,      /* fuzzy exponent (in) */
  double *U,           /* partition matrix (out) */
  const double *xc     /* cluster centers (in) */
)
/* NOTE: Will NOT check input arguments !!!!!!!!!!!!!!!!!! */
{
  register int i, j;
  int sing;    /* number of singularities */
  register long k;
  double sum,
         NewValue,   /* new membership value */
         Dik2, Djk2, /* distance squared (vi,xk) and (vj,xk) */
         sq_error;   /* fuzzy objective function value */
 
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
        /* don't check for changes, this routine is used for a single pass */
        U[mind2(i,k,C)] = NewValue;
      }
      else {    /* singularity (very unlikely case) */
        /* all previous memberships should be forced to zero */
        for (j = 0; j < i; j++) U[mind2(j,k,C)] = 0.0;
        /* Count how many singularities for this sample */
        for (j = i+1, sing = 1; j < C; j++)
          if (CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)]) == 0) sing++;
        /* Set the rest of memberships to 0 or 1/sing */
        U[mind2(i,k,C)] = 1.0/sing;
        for (j = i+1; j < C; j++) {
          if (CLTDistance (&x[mind2(0,k,d)], &xc[mind2(0,j,d)]) == 0)
            U[mind2(i,k,C)] = 1.0/sing;
          else
            U[mind2(i,k,C)] = 0.0;
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

  return (sq_error);
} /* do_extfpm */






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
      idumm, idumm2, i, j,
      dflag;   /* distance type */

  long n,     /* number of patterns */
       il, jl, ldumm;

  double *x,  /* input pattern matrix */
         *xc, /* cluster centers */
         *A,  /* distance matrix at input */
         *B,  /* distance matrix at output */
         *U,  /* final partition matrix */
         *err_p, /* pointer to objective function value */
         m,  /* fuzzy exponent */
         sq_error, /* value of objective function */
         dumm,
         *xMeans,  /* means of input patterns */
         *xVars;   /* variances of input patterns */

  mxArray *dmatrix_stc[2];    /* distance matrix structure */

  /*----------------------------------*/
  /* Check the input/output arguments */
  /*----------------------------------*/

  if (nin < 2) mexErrMsgTxt ("insufficient number of input arguments");
  if (nin > 5) mexErrMsgTxt ("too many input arguments");
  if (nout > 3) mexErrMsgTxt ("too many output arguments");

  /* fetch input arguments */
  /* x matrix (d x n) each column is one pattern vector */
  d = mxGetM(x_IN); n = (long)mxGetN(x_IN);
  if (n < 2) mexErrMsgTxt ("should have at least two patterns");
  if (d < 1) mexErrMsgTxt ("invalid number of features (pattern matrix)");
  if (d > n) mexPrintf ("WARNING: more features than patterns\n");

  /* xc matrix (d x C) each column is one cluster center */
  idumm = mxGetM(xc_IN); C = mxGetN(xc_IN);
  if (C < 2) mexErrMsgTxt ("requested number of clusters less then two");
  if (idumm != d) mexErrMsgTxt ("invalid number of features (center matrix)");
  if (C > n) mexErrMsgTxt ("more clusters than patterns");
  if (C > (double)n/5/(d+1))
    mexPrintf ("WARNING: number of patterns may be too small\n");

  /* data types for required arguments */
  if (!mxIsNumeric(x_IN) || !mxIsDouble(x_IN) ||
      !mxIsNumeric(xc_IN) || !mxIsDouble(xc_IN))
    mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
  if (mxIsComplex(x_IN) ||
      mxIsComplex(xc_IN))
     mexErrMsgTxt (msgs[MSG_REALONLY]);
  if (mxIsSparse(x_IN) ||
      mxIsSparse(xc_IN))
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

  /* dflag : distance type (optional, default=DIST_IP_EUCLIDEAN (0) ) */
  if (nin < 4)
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
  if (dflag < DIST_IP_EUCLIDEAN || dflag > DIST_IP_MAHALANOBIS)
    mexErrMsgTxt ("invalid distance type");

  /* A matrix (d x d) used in distance calculation */
  A = (double*)NULL;
  if (nin >= 5) {  /* matrix A is given */
    if (dflag == DIST_IP_EUCLIDEAN)   /* ignore it */
      mexPrintf ("WARNING: distance matrix ignored\n");
    else { /* use it */
      if (mxGetM(A_IN) != d || mxGetN(A_IN) != d) /* check dimensions */
        mexErrMsgTxt ("invalid number of dimensions (distance matrix)");
      /* check data type */
      if (!mxIsNumeric(A_IN) || !mxIsDouble(A_IN))
        mexErrMsgTxt (msgs[MSG_INVTYPINPARG]);
      if (mxIsComplex(A_IN)) mexErrMsgTxt (msgs[MSG_REALONLY]);
      if (mxIsSparse(A_IN)) mexErrMsgTxt (msgs[MSG_FULLONLY]);

      /* fetch the pointer to distance matrix */
      A = mxGetPr (A_IN);
      if (A == (double*)NULL) mexErrMsgTxt ("internal error no. 1");
      /* the previous line may look strange, but at one point I may
         decide to allow empty matrix to be specified. In that case
         A==NULL would signal such situation and has to be eliminated
         at this point. */

      if (dflag == DIST_IP_DIAGONAL) {  /* must be diagonal matrix */
        for (j = 0; j < d; j++)
          for (i = 0; i < d; i++)
            if (i != j && A[mind2(i,j,d)] != 0)
              mexErrMsgTxt ("distance matrix must be diagonal");
      }
    }
  }

  /* fetch the pointers to vector data */
  x = mxGetPr (x_IN);
  xc = mxGetPr (xc_IN);

  /*----------------------*/
  /* Allocate some memory */
  /*----------------------*/

  /* using MATLAB matrix structures */
  /* U matrix (C x n) each row is one membership function */
  U_OUT = mxCreateDoubleMatrix(C, n, mxREAL);
  U = mxGetPr(U_OUT);    /* initial values not needed */

  /*-------------------------------*/
  /* Prepare for do_extfpm routine */
  /*-------------------------------*/

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
    if (A == (double*)NULL) {  /* distance matrix not specified */
      /* get mean */
      xMeans = (double *)mxCalloc (d, sizeof(double));
      for (i = 0; i < d; i++) {
        xMeans[i] = 0;
        for (jl = 0L; jl < n; jl++)
          xMeans[i] += x[mind2(i,jl,d)];
        xMeans[i] /= (double)n;
      }
    }

    if (dflag == DIST_IP_DIAGONAL) {  /* diagonal */
      distfun_p = dist_diag;
      dmatrix_stc[0] = mxCreateDoubleMatrix (d, 1, mxREAL);
      dmatrix_stc[1] = (mxArray*) NULL;
      dmatrix = mxGetPr(dmatrix_stc[0]);

      if (A != (double*)NULL) /* use input matrix A if specified */
        for (i = 0; i < d; i++) dmatrix[i] = A[mind2(i,i,d)];
      else {  /* estimate distance matrix */
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
    }
    else {  /* mahalanobis */
      distfun_p = dist_mahal;

      /* create mahalanobis distance matrix */
      if (A != (double*)NULL) { /* use input matrix A if specified */
        dmatrix_stc[0] = mxCreateDoubleMatrix (d, d, mxREAL);  /* distance matrix */
        dmatrix_stc[1] = (mxArray*) NULL;
        dmatrix = mxGetPr(dmatrix_stc[0]);
        for (j = 0; j < d; j++)
          for (i = 0; i < d; i++)
            dmatrix[mind2(i,j,d)] = A[mind2(i,j,d)];
      }
      else {  /* estimate distance matrix */
        dmatrix_stc[1] = mxCreateDoubleMatrix (d, d, mxREAL);  /* scratch matrix */
        dmatrix = mxGetPr(dmatrix_stc[1]);
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
    }
    if (A == (double*)NULL) mxFree(xMeans); /* distance matrix not specified */
    break;
  default:
    mexErrMsgTxt ("internal error no. 2");
  }

  /*---------------------------------------------------------*/
  /* execute do_extfpm, i.e. estimate fuzzy partition matrix */
  /*---------------------------------------------------------*/

  sq_error = do_extfpm (C, n, d, x, m, U, xc);

  /*-----------------------------------*/
  /* check and create output variables */
  /*-----------------------------------*/

  /* objective function error */
  err_p = (double*)NULL;
  if (nout > 1) {
    err_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    err_p = mxGetPr(err_OUT);
    *err_p = sq_error;

    /* distance matrix at output */
    if (nout > 2) {
      B_OUT = mxCreateDoubleMatrix(d, d, mxREAL);
      B = mxGetPr(B_OUT);
      switch (dflag) {
      case DIST_IP_EUCLIDEAN:
        for (j = 0; j < d; j++)
          for (i = 0; i < d; i++)
            B[mind2(i,j,d)] = (i==j?1:0);
        break;
      case DIST_IP_DIAGONAL:
        for (j = 0; j < d; j++)
          for (i = 0; i < d; i++)
            B[mind2(i,j,d)] = (i==j?dmatrix[i]:0);
        break;
      case DIST_IP_MAHALANOBIS:
        for (j = 0; j < d; j++)
          for (i = 0; i < d; i++)
            B[mind2(i,j,d)] = dmatrix[mind2(i,j,d)];
        break;
      default:
        mexErrMsgTxt ("internal error no. 3");
      }
    }
  }

} /* mexFunction */
