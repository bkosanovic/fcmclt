function [Si,Fi] = fscat (x, xc, u, m)
%FSCAT    Calculate fuzzy scatter and covariance matrices for a fuzzy cluster.
%         [Si,Fi]=fscat(x,xc,u,m) estimates the fuzzy scatter and covariance
%         matrices for a fuzzy cluster with center at xc (d x 1 real vector),
%         and membership function u (1 x n real vector). Data samples are
%         provided in x (d x n real matrix). Fuzzy exponent 'm' is optional
%         (real >= 1), default is m=2. Output is a fuzzy scatter matrix
%         Si (d x d real matrix) and a fuzzy covariance matrix Fi (d x d).
%         Please note that some authors define the fuzzy scatter and fuzzy
%         covariance matrices with m=1.

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

% History: 29-Apr-95 Note on m=1 added. (V1-7)
%          15-Mar-95 Created

% Check the input arguments
if nargin < 4, m = 2; end;
if m < 1, error ('invalid fuzzy exponent'); end;

[d,n]=size(x);
if d < 1 | n < 2, error ('invalid size of sample matrix'); end;
if ~all([min(size(xc)) min(size(u))] == 1)
  error ('xc and u must be vectors');
end;
xc = xc(:);
u = u(:); u = u';
if size(xc,1) ~= d, error ('incompatible dimensionality x,xc'); end;
if size(u,2) ~= n, error ('incompatible dimensionality x,u'); end;

% Get the scatter matrix
xcmx = xc * ones(1,n);   % copy center n times (into n columns) (d x n matrix)
dif = x - xcmx;          % get difference between samples and a center
sdif = (ones(d,1) * u.^m) .* dif;   % scale difference with u^m
Si = sdif * dif';    % get the scatter matrix (d x d)

% Get the fuzzy covariance matrix
if nargout > 1
  Fi = Si / sum(u.^m);
end;
