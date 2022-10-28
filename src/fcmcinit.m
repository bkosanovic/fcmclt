function U0 = fcmcinit (C, n, tflag, arg1, arg2, arg3, arg4)
%FCMCINIT Generate initial fuzzy partition matrix.
%         U0 = fcmcinit (C,n,tflag,arg1,arg2,arg3,arg4) initializes the fuzzy
%         partition matrix U0 (C x n) for C clusters and n samples.
%         Initialization type is tflag (optional, default=0):
%
%         0 - uniform, random, with constraint that columns sum to
%             unity (first generates whole matrix without constraint
%             and then scales each column independently).
%
%             arg1 : seed (optional, default is 0 = startup value)
%
%         1 - same as 0, but with more complicated algorithm as in
%             Bezdek's paper: Computers & Geosciences Vol. 10,
%             no. 2-3, pg. 191-203, 1984. This algorithm tends
%             to assign larger memberships to top rows.
%
%             arg1 : seed (optional, default is 0 = startup value)
%
%         2 - uses arg2, to specify initial cluster centers. A half iteration
%             of the fcmc algorithm will be performed to estimate initial
%             partition matrix U0. Additional arguments are:
%
%             arg1 : x  - pattern sequence (d x n real matrix)
%             arg2 : xc - center matrix (d x C real matrix)
%             arg3 : m  - fuzzy exponent (real > 1), optional, default=2.00
%             arg4 : dflag - distance flag, optional, (0: Euclidean,
%                            1: diagonal, 2: mahalanobis), default=0
%
%             Note: if you specify xc to have equal columns (all centers
%                   overlap), clustering algorithm would fail producing
%                   U=1/C*ones(C,n), without returning objective function
%                   value nor center trajectories. All final centers would be
%                   at the center of mass for 'x'.

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

% Copyright (c) 1994 by Bogdan R. Kosanovic
% University of Pittsburgh, PA, USA

% History: 19-Nov-94 name changed from fcminit to fcmcinit (V1-1)
%                    tflag=2 added
%          13-Nov-94 tflag=1 added from sigp package
%          11-Sep-94 Created

if nargin < 3, tflag = 0; end;
if C < 2 | C >= n, error ('invalid C or n'); end;
if tflag < 0 | tflag > 2, error ('unknown initialization type'); end;

if tflag == 0   % uniform, random
  if nargin < 4, arg1 = 0; end;
  seed = arg1;
  rand('seed',seed);
  U0 = rand(C,n);
  divs = sum(U0);
  U0 = U0 ./ (ones(C,1)*divs);
elseif tflag == 1  % bezdek type initialization
  if nargin < 4, arg1 = 0; end;
  seed = arg1;
  rand('seed',seed);
  U0 = zeros(C,n);
  for k = 1:n
    maxmem = 1.0;
    for ii = 1:C-1
      next = rand/2;
      nc = C - ii + 1;
      next = maxmem * (1 - next^(1.0/nc));
      U0(ii,k) = next;
      maxmem = maxmem - next;
    end;
    if maxmem < 0, maxmem = 0; disp ('*** negative membership'); end;
    U0(C,k) = maxmem;
  end;
elseif tflag == 2  % initialization using cluster centers
  if nargin < 5, error ('not enough input arguments'); end;
  [d,dumm]=size(arg1);
  if dumm ~= n, error ('invalid number of patterns'); end;
  if d < 1, error ('invalid number of features (pattern matrix)'); end;

  xc = arg2; clear arg2;
  [dumm,dumm2]=size(xc);
  if dumm ~= d | dumm2 ~= C, error ('invalid center matrix'); end;

  % check for m and dflag
  if nargin < 6
    m = 2; dflag = 0;
  else
    m = arg3;
    if nargin < 7, dflag = 0; else, dflag = arg4; end;
  end;

  U0 = extfpm (arg1, xc, m, dflag);   % extrapolate U0 from x and xc
end;
