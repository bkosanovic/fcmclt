function [U,err,B] = extfpm (x, xc, m, dflag, A)
%EXTFPM   Extrapolate fuzzy partition matrix.
%         [U,err,B]=extfpm(x,xc,m,dflag,A) extrapolates the partition matrix
%         using the presumed cluster centers 'xc' (d x C real matrix) to fit
%         the pattern sequence in 'x' (d x n real matrix). For more information
%         on center and pattern matrices see FCMC. Partition matrix is returned
%         in U (C x n real matrix), and value of objective function in 'err'.
%         Distance matrix used in estimation of U is returned in B (d x d).
%         Optional input parameters are:
%
%         m     : fuzzy exponent (real > 1), default=2.00
%         dflag : distance flag (0: Euclidean, 1: diagonal, 2: mahalanobis),
%                 default=0
%         A     : distance matrix (d x d real matrix). When dflag=1 or 2,
%                 A is used to specify the metric in feature space.
%                 For dflag=1, A must be a diagonal matrix. One should use the
%                 matrix A returned from the FCMC routine as a distance matrix.
%                 For dflag=0, A is ignored. If A is not specified it will
%                 be estimated from the input data in x.
%
%         This function is often used to extrapolate the Temporal Fuzzy Sets
%         outside the time window that was used for their estimation. It
%         can also be used to obtain the initial partition matrix for
%         the postulated cluster centers, or to classify unknown patterns (x)
%         given the estimated prototypes (xc).
%
%         WARNING:  In case of non-Euclidean distance, it may not make much
%         sense to extrapolate large quantities of data, since the statistics
%         may be rather different from the one that was used in estimation of
%         the cluster prototypes. Namely, the new metric may not correspond to
%         the metric used in estimation.
%
%         This problem may be reduced by not omitting the matrix A which
%         modifies the metric for new data, but providing the original matrix
%         used in estimation of the prototypes as returned by FCMC routine.
%         Then, the extrapolation will be with respect to the original metric
%         used within the FCMC routine.
%
%         EXTFPM (in its current version) should not be used for extrapolation
%         of results obtained with the GKFCMC algorithm since it would not use
%         the resulting fuzzy covariance matrices.
%
%         See also: FCMC

% MIT License
%
% Copyright (c) 1995-2022 Bogdan Kosanovic
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Copyright (c) 1994,1995 by Bogdan R. Kosanovic
% University of Pittsburgh, PA, USA

% History: 14-May-95 Comment relating to GKFCMC added. (V2-0)
%          20-Mar-95 Distance matrix A can be specified as input. Distance
%                    matrix that was used will be returned in B. (V1-4)
%          11-Dec-94 Warning comment added (V1-2)
%          19-Nov-94 Created (V1-1)

% source is in extfpm.c
