% testfcmc.m  - test the fcmc algorithm
%
% xc should be like 
%
%        11.3783  45.6185
%        17.6907  22.1831        for Euclidean distance (dflag=0)
%
% A should be like
%
%         0.0040  -0.0048
%        -0.0048   0.0336        for Mahalanobis distance (dflag=2)

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


x = [10  8;
      9 14;
      8 20;
     12 24;
     17 22;
     34 20;
     42 19;
     48 16;
     52 20;
     48 26;
     46 32];

v = [14 18;
     47 21]';   % presumed centers  (these were used in one paper ;-) ;-)

U0 = fcmcinit (C,n);  % initialize the membership values (startup seed)

dflag = input (' Distance type (0,1,2) ? ');
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

niter = max(size(err));

figure(1);
plot(x(:,1),x(:,2),'o');
xmax = max([x(:,1); xctraj(1,:)'; 60]);
ymax = max([x(:,2); xctraj(2,:)'; 40]);
set(gca,'xlim',[0 xmax],'ylim',[0 ymax]);
hold on;
for k = 1:C
  plot(xctraj(1,k:C:C*niter), xctraj(2,k:C:C*niter), mark(k));
end;
hold off;
title ('samples and the center trajectories');

figure(2);
plot(err);
ylim=get(gca,'ylim');
set(gca,'ylim', [0 ylim(2)]);
title('objective function values');

figure(3);
t = 1:n;
plot(t,U');
title('membership functions');
