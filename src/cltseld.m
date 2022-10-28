function cltseld (dtype, A)
%CLTSELD  Select distance type for CLTDIST2 function.
%         cltseld (dtype,A) selects the type of inner product for CLTDIST2
%         function.
%
%           dtype : 0 - Euclidean, 1 - diagonal, 2 - Mahalanobis
%           A     : distance matrix
%                   dtype == 0, A not required (== identity)
%                   dtype == 1, d x d diagonal matrix
%                   dtype == 2, d x d matrix
%                   
%         Note: this functions makes use of FCMC_xxx global variables.
%
%         See also: CLTDIST2

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

% History: 22-Apr-95 Created

% access global variables
global FCMC_DIST_type        % distance type
global FCMC_DIST_amat        % distance matrix


% check the input arguments
if dtype < 0 | dtype > 2, error ('invalid distance type'); end;
if dtype > 0 & nargin < 2, error ('missing input argument'); end;

if dtype == 0  % Euclidean
  FCMC_DIST_type = 0;
  FCMC_DIST_amat = [];
else
  [d1,d2]=size(A);
  if d1 ~= d2 | d1 < 1, error ('invalid distance matrix dimensions'); end;
  FCMC_DIST_amat = dtype;
  if dtype == 1                % Diagonal
    FCMC_DIST_amat = diag(A);     % creates d x 1 vector
  else                         % Mahalanobis
    FCMC_DIST_amat = A;
  end;
end;
