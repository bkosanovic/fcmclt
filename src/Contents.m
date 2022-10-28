% Fuzzy c-Means Clustering Toolbox
% Version 3-0  17-Jul-2019
% Copyright (c) 1991-2022 by Bogdan R. Kosanovic
%
% Fuzzy c-Means Clustering
%   extfpm     - Extrapolate fuzzy partition matrix
%   fcmc       - Fuzzy c-means clustering
%   fcmcinit   - Generate initial fuzzy partition matrix
%   gkfcmc     - Gustafson-Kessel fuzzy c-means clustering
%
% Clustering Validity Functionals
%   cltvalid   - Calculate clustering validity functionals
%
% Other
%   fscat      - Calculate fuzzy scatter and covariance matrices for
%                a fuzzy cluster
%   testextfpm - test the fuzzy partition matrix extrapolation
%   testfcmc   - test the fcmc algorithm
%   testgk     - test the gkfcmc using two Gaussian classes
%   testgk2    - test the gkfcmc using Gustafson's cross
%   tstvalid   - test the validity functionals (using fcmc)

% Hidden utility functions for toolbox.
%   cltdist2   - Calculate distance squared between two points (inner
%                product space)
%   cltseld    - Select distance type for CLTDIST2 function

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
