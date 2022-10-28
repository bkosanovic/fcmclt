% tstvalid.m  - test the validity functionals
%
% should produce tables like these (distance type 0 - euclidean)
%
% C:    F      H     NFI   minHT meanHT  minRF  maxRF minNMM
% 2: 0.9106 0.1766 0.8212 1.4262 1.2813 0.0000 0.0000   5
% 3: 0.7845 0.3924 0.6768 1.2859 0.8380 0.0000 0.3333   3
% 4: 0.7411 0.4925 0.6549 1.1731 0.8707 0.0000 0.5000   2
% 5: 0.7570 0.4988 0.6962 1.1657 0.8841 0.0000 0.5000   2
%
% C:      S       Fhv      Dpa       Pd 
% 2:   0.0413  88.5299   0.0662   0.0661
% 3:   0.3302 122.5980   0.0407   0.0383
% 4:   0.1853 113.9024   0.0672   0.0643
% 5:   0.0758 102.7258   0.0806   0.0791
%
% Fhv, Dpa and Pd are not reliable for the small number of samples!

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

Crange = (2:5)';
kloops = length(Crange);
n = 11;
d = 2;
m = 2.0;
maxiter = 130;
mineps = 1e-4;

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

dflag = input (' Distance type (0,1,2) ? ');

F=zeros(kloops,1); H=F; NFI=F; MinHT = F; MeanHT = F; MinRF = F;
MaxRF = F; MinNMMcard = F; S = F; Fhv = F; Dpa = F; Pd = F;
for k = 1:kloops
  C = Crange(k);
  disp(sprintf('...clustering into %d classes', C));
  U0 = fcmcinit (C,n);  % initialize the membership values (startup seed)
  [xc,U,fineps,err,xctraj,A] = fcmc(x.',U0,m,maxiter,mineps,dflag);

  [F(k), H(k), NFI(k), MinHT(k), MeanHT(k), MinRF(k), MaxRF(k), ...
   MinNMMcard(k)] = cltvalid(U);
  [S(k), Fhv(k), Dpa(k), Pd(k)] = cltvalid(U,x.',xc,m,dflag,A);
end;
disp(sprintf (' C: %6s %6s %6s %6s %6s %6s %6s %6s', 'F  ', 'H  ', 'NFI ', ...
              'minHT', 'meanHT', 'minRF', 'maxRF', 'minNMM'));
disp(sprintf (' %d: %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %3d\n', ...
     [Crange, F, H, NFI, MinHT, MeanHT, MinRF, MaxRF, MinNMMcard]'));

disp(sprintf (' C: %8s %8s %8s %8s', 'S  ', 'Fhv ', 'Dpa ', 'Pd '));
disp(sprintf (' %d: %8.4f %8.4f %8.4f %8.4f\n', [Crange, S, Fhv, Dpa, Pd]'));
