function dd = cltdist2 (x, y)
%CLTDIST2 Calculate distance squared between two points (inner product space).
%         dd = cltdist2 (x,y) calculates distance squared between points x and
%         y in inner product space that was selected by CLTSELD function.
%         Default space is Euclidean. Input vectors x and y must be of
%         same dimensions (d x 1), where 'd' is a dimension of corresponding
%         vector space.
%
%         Note: this functions makes use of FCMC_xxx global variables.
%
%         See also: CLTSELD

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
if ~all([min(size(x)) min(size(y))] == 1),
  error ('x and y must be (nonempty) vectors');
end;
x=x(:); y=y(:); d = size(x,1);
if size(y,1) ~= d, error ('x and y must be of same dimensions'); end;

% check distance type (initialize to Euclidean if necessary)
if isempty(FCMC_DIST_type)
  cltseld(0);
  disp ('WARNING: Euclidean distance selected');
end;

% calculate distance squared
dumm = x - y;
if FCMC_DIST_type == 0    % Euclidean
  dd = dumm.' * dumm;     % distance squared
elseif FCMC_DIST_type == 1  % Diagonal
  dd = dumm.' * (FCMC_DIST_amat .* dumm);
elseif FCMC_DIST_type == 2  % Mahalanobis
  dd = dumm.' * FCMC_DIST_amat * dumm;
else
  error ('internal error (invalid distance type)');
end;
