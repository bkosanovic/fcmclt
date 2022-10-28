%GKFCMC   Gustafson-Kessel fuzzy c-means clustering.
%
%         [xc,U,fineps,err,xctraj]=gkfcmc(x,U0,m,maxiter,maxeps,volc) performs
%         the Gustafson-Kessel modification of fuzzy c-means clustering.
%         Each of the clusters will use its own IP induced metric. Hence, the
%         resulting clusters may have arbitrary ellipsoidal shapes.
%
%         Input:
%
%         x       : pattern matrix (d x n real matrix)
%                   d > 0 (number of dimensions), n > 1 (number of patterns)
%         U0      : initial partition matrix (C x n real matrix)
%                   C > 1 (number of clusters), C <= n (better be < n)
%         m       : fuzzy exponent (real > 1), optional, default=2.00
%         maxiter : maximum number of iterations (> 0), optional, default=100
%         maxeps  : maximum norm error allowed (> 0), optional, default=1e-3,
%                   a stopping criterion, suggested range (1e-3 to 5e-4)
%         volc    : volume constraints (real vector of length C). Contains
%                   C positive real numbers that constrain the determinants of
%                   matrices Aj which are shaping the resulting clusters.
%                   optional, default=ones(C,1)
%
%         Output:
%
%         xc     : cluster centers (d x C real matrix)
%                  each column stands for a center
%         U      : final partition matrix (C x n real matrix)
%                  each row represents a membership function
%         fineps : final maximum norm error (scalar)
%         err    : objective function values (niter x 1 real vector),
%                  where niter is the final number of iterations
%         xctraj : center trajectories (d x C*niter real matrix) composed of
%                  niter submatrices (d x C) for each iteration. Submatrices
%                  have same structure as xc matrix. The last submatrix
%                  equals xc.
%
%         See also: FCMC, FSCAT, CLTVALID
%
%         Note: if U0=b*ones(C,n), i.e. all entries are the same, clustering
%               algorithm would fail, producing the centroid of degenerate
%               fuzzy c-partition space, U=1/C*ones(C,n). All final centers 
%               would be at the center of mass for 'x'.
%
%         For more information on this algorithm, please refer to
%
%             J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
%                           Function Algorithms", Plenum Press, New York, 1981.

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

% Copyright (c) 1995 by Bogdan R. Kosanovic
% University of Pittsburgh, PA, USA

% History: 05-May-95 Created

% source is in gkfcmc.c
