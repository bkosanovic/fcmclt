% testgk2.m  - test the gkfcmc algorithm using Gustafson's cross
%
% xc should be like
%
%    0.2399    0.1875
%    0.0681   -1.6713
%
% F1
%
%   37.8131   -0.0199          <--- large variance in x1
%   -0.0199    0.0782
%
% F2
%
%    0.1240    1.4643
%    1.4643   23.8536          <--- large variance in x2

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

C = 2; mark = 'x+';  % number of clusters and symbols for center trajectories
colr = [0, 0, 1;     % Cluster 1 is blue
        1, 0, 0];    % Cluster 2 is red
n = 20;              % number of data points
d = 2;               % feature space dimensions
m = 2.0;             % fuzzy exponent to use
maxiter = 100;       % maximum number of iterations
mineps = 1e-3;       % minimum error to stop iterations

% Gustafson's cross (n x d)
x = [-9.75 -0.15;
     -6.44  0.34;
     -4.69 -0.3;
     -2.04  0.37;
     -1.24  0.45;
      0.33 -0.08;
      5.04 -0.21;
      5.86 -0.25;
      7.54  0.16;
      7.67  0.24;
     -0.3  -8.07;
      0.13 -7.13;
     -0.37 -5.18;
      0.03 -3.33;
      0.35 -2.63;
      0.23 -2.68;
     -0.05 -2;
      0.41  0.37;
      0.69  4.75;
      0.74  8.87];

% U0 (C x n): the initial fuzzy partition matrix, it contains initial membership values
U0 = fcmcinit (C,n,1);  % initialize the membership values (startup seed)

[xc,U,fineps,err,xctraj] = gkfcmc(x.',U0,m,maxiter,mineps);

niter = max(size(err));   % Actual number of iterations

% find NMM labels
[dumm,labs]=max(U);
k1 = find(labs==1);
k2 = find(labs==2);

% get fuzzy covariance matrices
[S1,F1]=fscat(x',xc(:,1),U(1,:));
[S2,F2]=fscat(x',xc(:,2),U(2,:));
F1
F2

figure(1);
plot(x(:,1),x(:,2),'ok');
xmax = max([x(:,1); xctraj(1,:)'; 10]);
xmin = min([x(:,1); xctraj(1,:)'; -10]);
ymax = max([x(:,2); xctraj(2,:)'; 10]);
ymin = min([x(:,2); xctraj(2,:)'; -10]);
set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax]);
hold on;
for k = 1:C
  plot(xctraj(1,k:C:C*niter), xctraj(2,k:C:C*niter), mark(k), 'color', colr(k,:));
end;
hold off;
grid;
xlabel('Feature 1');
ylabel('Feature 2');
title ('Samples and the center trajectories');

figure(2);
plot(err,'k'); grid;
ylim=get(gca,'ylim');
set(gca,'ylim', [0 ylim(2)]);
xlabel('Iterations');
title('Objective function values');

figure(3);
t = 1:n;
plot(t,U(1,:),'color',colr(1,:));
hold on;
plot(t,U(2,:),'color',colr(2,:));
hold off;
grid;
xlabel('Data points');
ylabel('Membership values');
title('Membership functions');

% Display clustering results where color of a point is a mixture based on membership values
figure(4);
for k = 1:n
  plot(x(k,1),x(k,2),'.', 'MarkerSize', 15, 'color', U(1,k)*colr(1,:)+U(2,k)*colr(2,:)); hold on;
end;
hold off;
xmax = max([x(:,1); xctraj(1,:)'; 10]);
xmin = min([x(:,1); xctraj(1,:)'; -10]);
ymax = max([x(:,2); xctraj(2,:)'; 10]);
ymin = min([x(:,2); xctraj(2,:)'; -10]);
set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax]);
hold on;
for k = 1:C
  plot(xctraj(1,C*(niter-1)+k), xctraj(2,C*(niter-1)+k), mark(k), ...
      'LineWidth', 1, 'MarkerSize', 10, 'color', colr(k,:));
end;
hold off;
grid;
xlabel('Feature 1');
ylabel('Feature 2');
title ('Clustering result with cluster centers');

% nothing past this point
