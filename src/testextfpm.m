% testextfpm.m  - test the fuzzy partition matrix extrapolation
%
% xc should be like 
%
%        11.5631  45.4449
%        17.8645  22.0322        for Euclidean distance (dflag=0)
%
% A should be like
%
%         0.0040  -0.0048
%        -0.0048   0.0372        for Mahalanobis distance (dflag=2)
%
% Initial centers (C0) are
%
%         5     30
%         5     20
%
% The test data is a bit different from testfcmc.m test data

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

clear all

C = 2; mark = 'x+';
n = 11;
d = 2;
m = 2.0;
maxiter = 100;
mineps = 5e-5;


x = [11  9;
      9 14;
      8 20;
     12 24;
     17 22;
     34 20;
     42 19;
     48 16;
     52 20;
     48 26;
     45 31];

v = [14 18;
     47 21]';   % presumed centers  (these were used in one paper ;-) ;-)

% This particular initialization of FCMC algorithm also uses extfpm() !!!
C0 = [5 5; 30 20].';            % initial center positions
U0 = fcmcinit (C,n,2,x.',C0);   % initialize the membership values (startup seed)

% Select the distance type
dflag = input (' Distance type (0,1,2) ? ');

% Estimate the centers
[xc,U,fineps,err,xctraj,A] = fcmc(x.',U0,m,maxiter,mineps,dflag);

% xc should be like 
%
%        11.3783  45.6185
%        17.6907  22.1831        for Euclidean distance (dflag=0)
%
% A should be like
%
%         0.0040  -0.0048
%        -0.0048   0.0336        for Mahalanobis distance (dflag=2)

% Use presumed centers to extrapolate the fuzzy partition matrix
[U_e,err_e,B_e]=extfpm(x.',v,m,dflag);

niter = max(size(err));

figure(1);
plot(x(:,1),x(:,2),'o');
xmax = max([x(:,1); xctraj(1,:)'; 60]);
ymax = max([x(:,2); xctraj(2,:)'; 40]);
set(gca,'xlim',[0 xmax],'ylim',[0 ymax]);
hold on;
for k = 1:C
  plot(xctraj(1,k:C:C*niter), xctraj(2,k:C:C*niter), mark(k));
  plot(xctraj(1,C*niter-C+k), xctraj(2,C*niter-C+k), mark(k),'linew',2);
  plot(v(1,k),v(2,k),[mark(k),'r'],'linew',2)
end;
hold off;
title ('samples and the center trajectories');

figure(2);
plot(err,'--');   % as provided by fcmc
hold on
  plot([1 length(err)],err_e*[1 1],'r');   % as provided by extrapolation with presumed centers
hold off
ylim=get(gca,'ylim');
set(gca,'ylim', [0 ylim(2)]);
title('objective function values');

figure(3);
t = 1:n;
plot(t,U','--');     % as provided by fcmc
hold on;
  plot(t,U_e');      % as provided by extrapolation with presumed centers
hold off;
title('membership functions');
