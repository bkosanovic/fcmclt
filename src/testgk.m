% testgk.m  - test the gkfcmc algorithm using two Gaussian classes
%
% class 1 is N([ 3  3], [1.7 0.5])   uncorrelated
% class 2 is N([-4 -4], [1   1.3])   uncorrelated
%
% xc should be like
%
%    3.2590   -3.9495
%    2.9463   -4.1688
%
% F1
%
%    2.7873    0.0111     ---> 1.67 and 0.48 as estimates of std's
%    0.0111    0.2316
%
% F2
%    1.0437   -0.0901     ---> 1.02 and 1.38 as estimates of std's
%   -0.0901    1.9201

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
colr = [0, 0.5, 0;   % Cluster 1 is dark green
        1, 0, 0];    % Cluster 2 is red
n1 = 200; n2 = 100;  % number of data points for each class
n = n1+n2;           % total number of data points
d = 2;               % feature space dimensions
m = 2.0;             % fuzzy exponent to use
maxiter = 100;       % maximum number of iterations
mineps = 1e-3;       % minimum error to stop iterations

% Generate data points (input) (n x d)
x = zeros(n,d);

% generate Gaussian samples of class 1
randn('seed',13);   % set seed
x(1:n1,:) = ones(n1,1)*[3 3]+(ones(n1,1)*[1.7 0.5]).*randn(n1,2);

% generate Gaussian samples of class 2
randn('seed',5431); % set seed
x(n1+1:n,:) = ones(n2,1)*[-4 -4]+(ones(n2,1)*[1 1.3]).*randn(n2,2);

% perform clustering
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
plot(x(:,1),x(:,2),'.b');
xmax = max([x(:,1); xctraj(1,:)'; 10]);
xmin = min([x(:,1); xctraj(1,:)'; -10]);
ymax = max([x(:,2); xctraj(2,:)'; 10]);
ymin = min([x(:,2); xctraj(2,:)'; -10]);
set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax]);
set(gca,'dataaspectratio',[1 1 1]);
hold on;
for k = 1:C
  plot(xctraj(1,k:C:C*niter), xctraj(2,k:C:C*niter), mark(k), 'color', colr(k,:));
end;
hold off;
xlabel('Feature 1');
ylabel('Feature 2');
title ('Samples and the center trajectories');

figure(2);
plot(err,'b');
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
xlabel('Data points');
ylabel('Membership values');
title('Membership functions');

