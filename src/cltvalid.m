function [out1, out2, out3, out4, out5, out6, out7, out8] = ...
         cltvalid (U,x,xc,m,dflag,A)
%CLTVALID Calculate clustering validity functionals.
%         [F,H,NFI,MinHT,MeanHT,MinRF,MaxRF,MinNMMcard] = cltvalid(U) returns
%         the partition coefficient (F), partition entropy (H), nonfuzzy
%         index (NFI), minimum and mean hard tendencies (MinHT, MeanHT),
%         minimum and maximum relative fuzziness (MinRF, MaxRF), and the
%         minimum nearest maximum membership (NMM) cardinality (MinNMMcard) of
%         a partition represented by matrix U (C x n real matrix, C: number of
%         clusters, n: number of samples). Optimal number of clusters should be
%         obtained when F is maximum, H is minimum, NFI, MinHT or MeanHT are
%         maximum over the range of values for C. MinRF, MaxRF, and MinNMMcard
%         are only the indicators of validity. MinRF and MaxRF measure minimum
%         and maximum relative fuzziness of NMM hard clusters, respectively.
%         Values close to one are indication of bad separation, while the
%         values close to zero suggest well separated clusters. When MinNMMcard
%         equals zero it indicates that there may be clusters that are "empty",
%         or, that the fuzzy clustering algorithm failed to separate some of
%         the cluster centers, i.e. some of the final centers may be extremely
%         close one to another.
%
%         [S,Fhv,Dpa,Pd] = cltvalid(U,x,xc,m,dflag,A) returns the compactness
%         and separation index (S), fuzzy hypervolume (Fhv), average partition
%         density (Dpa), and partition density (Pd) of a partition represented
%         by matrix U (C x n), pattern matrix x (d x n), cluster centers
%         xc (d x C).
%
%         Optional arguments:
%
%           m     : fuzzy exponent (> 1), default=2 (for algorithms other
%                   than FCMC, m should always equal 2)
%           dflag : distance flag (0: Euclidean, 1: diagonal, 2: Mahalanobis),
%                   default=0
%           A     : matrix (d x d) that was used to calculate distance metric
%                   (returned by FCMC), default=identity
%
%         Optimal number of clusters should be obtained when S or Fhv are
%         minimum, Dpa or Pd are maximum over the range of values for C.
%         Optional arguments (m, dflag, A) will not directly affect
%         the calculation of Fhv, Dpa, and Pd.
%
%         See also: FCMC, FSCAT, GKFCMC
%
%
%         Notes on F:  For U in Mfc (fuzzy partition space)
%                      1/C <= F <= 1
%                      for F = 1, U is hard (zeros and ones only)
%                      for F = 1/C, U = 1/C*ones(C,n);
%
%         Notes on H:  For U in Mfc
%                      0 <= H <= log(C)
%                      for H = 0, U is hard
%                      for H = log(C), U = 1/C*ones(C,n);
%                      0 <= 1 - F <= H (strict inequality if U not hard)
%
%         For more information on F and H, please refer to
%
%             J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
%                           Function Algorithms", Plenum Press, New York, 1981.
%
%         For more information on NFI, please refer to
%
%             M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
%             Fuzzy Sets and Systems, 1:239-253, 1978.
%
%         For more information on MinHT and MeanHT, please refer to
%
%             F.F. Rivera, E.L. Zapata, and J.M. Carazo, "Cluster validity
%             based on the hard tendency of the fuzzy classification",
%             Pattern Recognition Letters, 11:7-12, 1990.
%
%         For more information on MinRF and MaxRF, please refer to
%
%             H.L. Gordon and R.L. Somorjai, "Fuzzy Cluster Analysis of
%             Molecular Dynamics Trajectories", Proteins: Structure,
%             Function, and Genetics, 14:249-264, 1992.
%
%         For more information on MinNMMcard and related functionals, please
%         refer to
%
%             B.R. Kosanovic et al., "Fuzzy Modeling of Dynamic Processes",
%             submitted for publication, 1995. (draft available upon
%             request).
%
%         Notes on MinNMMcard: if scaled with the number of samples,
%             MinNMMcard approximates the minimum a-priori probability over all
%             identified classes (MinNMMcard / n). Upper bound for the
%             minimum a-priori probability is 1/C. The number of clusters
%             for which MinNMMcard drops close to zero is an indication
%             of an upper bound for a "good" number of clusters.
%
%         For more information on S, please refer to
%
%             X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
%             IEEE Trans. PAMI, 13(8):841-847, 1991.
%
%         For more information on Fhv, Dpa, and Pd, please refer to
%
%             I. Gath and A.B. Geva, "Fuzzy clustering for the estimation of
%             the parameters of the components of mixtures of normal
%             distributions", Pattern Recognition Letters, 9:77-86, 1989.
%
%         Notes on Fhv, Dpa, and Pd: not reliable for the small number of
%             samples because the fuzzy covariance matrices cannot be
%             properly estimated. Fuzzy covariance matrices are always 
%             estimated with m=1 when calculating Fhv, Dpa, and Pd (see
%             Gath & Geva).

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

% History: 29-Apr-95 MinHT,MeanHT,MinRF,MaxRF,MinNMMcard,Fhv,Dpa,Pd (V1-7)
%          22-Apr-95 S added (V1-6)
%          27-Feb-95 Created (F, H, NFI only)

% Calculates F - partition coefficient (out1) (max)
%            H - partition entropy (out2) (min)
%            NFI - nonfuzzy index (out3) (max)
%            MinHT - minimum hard tendency (out4) (max)
%            MeanHT - mean hard tendency (out5) (max)
%            MinRF - minimum relative fuzziness (out6) (n/a)
%            MaxRF - maximum relative fuzziness (out7) (n/a)
%            MinNMMcard - minimum NMM cardinality (out8) (n/a)
%
%            S - compactness and separation index (out1) (min)
%            Fhv - fuzzy hypervolume (out2) (min)
%            Dpa - average partition density (out3) (max)
%            Pd - partition density (out4) (max)

[C,N] = size(U);
if C < 2 | C > N, error ('invalid partition matrix'); end;

if nargin == 1    % calculate F, H, and NFI
  sumF = sum(U(:).^2);
  pin = find(U > 0); sumH = sum(U(pin) .* log(U(pin)));

  F = sumF / N;  out1 = F;
  H = - sumH / N; out2 = H;
  NFI = (C*sumF - N) / N / (C-1); out3 = NFI;

  % get hard tendencies and relative fuzziness
  % start by computing first and second maximums
  [maxm1,ind1] = max(U);
  sU = U;
  for k = 1:N, sU(ind1(k),k) = -1; end;   % remove first max
  [maxm2,ind2] = max(sU);
  clear sU
  if any(maxm2 < 0), error ('internal error'); end;
  uThresh = (1+C)/2/C;   % threshold for relative fuzziness

  sumR = zeros(C,1); NMMcard = zeros(C,1); FuzCount = zeros(C,1);
  for ii = 1:C
    idx = find(ind1==ii);
    sumR(ii) = sum(maxm2(idx)./maxm1(idx));
    NMMcard(ii) = length(idx);
    FuzCount(ii) = length(find(maxm1(idx) < uThresh));
  end;

  % calculate "min" and mean hard tendencies
  idxOK = find(NMMcard ~= 0 & sumR ~= 0);
  logTs = -log10(sumR(idxOK) ./ NMMcard(idxOK));
  MinHT = max(logTs); out4 = MinHT;
  MeanHT = sum(logTs) / C; out5 = MeanHT;

  % calculate min and max relative fuzziness
  rf = ones(C,1);   % relative fuzziness for card==0 is maximum (1.0)
  idxOK = find(NMMcard > 0);
  rf(idxOK) = FuzCount(idxOK) ./ NMMcard(idxOK);
  MinRF = min(rf); out6 = MinRF;
  MaxRF = max(rf); out7 = MaxRF;
  MinNMMcard = min(NMMcard); out8 = MinNMMcard; % minimum NMM cardinality

else              % calculate S, Fhv, Dpa, and Pd
  % check input arguments
  [d, dumm] = size(x);
  if d < 1 | dumm ~= N, error ('invalid pattern matrix x'); end;
  if ~all([d,C] == size(xc)), error ('invalid center matrix xc'); end;

  % check optional arguments
  if nargin < 4, m = 2; end;
  if m <= 1, error ('invalid fuzzy exponent m'); end;
  if nargin < 5, dflag = 0; end;
  if dflag < 0 | dflag > 2, error ('invalid distance type'); end;
  if nargin < 6, A = eye(d); end;
  if ~all([d,d] == size(A)), error ('invalid distance matrix'); end;

  % select distance type
  cltseld(dflag,A);

  % calculate S, Fhv, Dpa, and Pd (rather slow due to for loops!)
  if C*N > 1000, disp ('... calculating S, Fhv, Dpa, and Pd'); end;
  Dzx = zeros(C,N); Dzz = zeros(C,C);
  Ss = 0; DpaSum = 0; Fhv = 0; xcenter = zeros(d,1);
  for ii = 1:C

    % ------ Fhv, Dpa, Pd block#1
    xcenter = xc(:,ii);
    % get fuzzy covariance matrix for this cluster
    [dumm,Fi] = fscat(x, xcenter, U(ii,:), 1);  % Gath and Geva use m=1

    sqdFi = sqrt(det(Fi));
    Fhv = Fhv + sqdFi;
    Si = 0; Finv = inv(Fi);
    % ------ end of Fhv, Dpa, Pd block#1

    for k = 1:N
      % ------ Fhv, Dpa, Pd block#2
      dxc = x(:,k) - xcenter;
      elpd = dxc.' * Finv * dxc;
      if elpd < 1, Si = Si + U(ii,k); end;
      % ------ end of Fhv, Dpa, Pd block#2

      Dzx(ii,k) = cltdist2(xcenter, x(:,k));   % S part
    end;

    % ------ Fhv, Dpa, Pd block#3
    DpaSum = DpaSum + Si / sqdFi;
    Ss = Ss + Si;
    % ------ end of Fhv, Dpa, Pd block#3

    % S part
    if ii < C
      for jj = ii+1:C
        Dzz(ii,jj) = cltdist2(xcenter, xc(:,jj));
      end;
    end;
  end;
  Dzz = Dzz(:);
  idx = find(Dzz > 0);
  if isempty(idx)
    S = inf;
  else
    S = sum(sum(Dzx.*U.^m)) / N / min(Dzz(idx));  % use m as specified!
  end;
  out1 = S;  % the first argument is S

  out2 = Fhv;                    % fuzzy hypervolume
  Dpa = DpaSum / C; out3 = Dpa;  % average partition density
  Pd = Ss / Fhv; out4 = Pd;      % partition density
  if C*N > 1000, disp ('... S, Fhv, Dpa, Pd done.'); end;

end;
