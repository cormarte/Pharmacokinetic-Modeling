%GROWROI grows a 2D region-of-interest 
%   V = GROWROI (S[, D, T, MAXITER]) grows the 2D region-of-interest (ROI)
%      'S': 2D input region-of-interest
%      'D': 2D input data set for data driven growing (default: ones copy
%           of S)
%      'T': threshold for data driven growing (default: 0.0)
%      'MAXITER': maximum number of growing iterations (default: 3)
%
%      The indexing order of S and D is Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function v = growROI (s, d, t, maxIter)
  if (nargin < 2) || isempty (d)
    d = ones (size (s));
  end
  if nargin < 4
    maxIter = 3;
  end
  if nargin < 3
    t = 0.0;
    maxiter = 1;
  end
  [ ny, nx ] = size (d);
  v = s;
  n = 1;
  m = 1;
  while (n > 0) & (m <= maxIter)
    [ y, x ] = find (s > 0);
    n = 0;
    for i=1:size(x,1)
      if y(i) > 1
        if ~v(y(i)-1,x(i)) & (d(y(i)-1,x(i)) > t)
          v(y(i)-1,x(i)) = 1;
          n = n + 1;
        end
      end
      if x(i) > 1
        if ~v(y(i),x(i)-1) & (d(y(i),x(i)-1) > t)
          v(y(i),x(i)-1) = 1;
          n = n + 1;
        end
      end
      if y(i) < ny
        if ~v(y(i)+1,x(i)) & (d(y(i)+1,x(i)) > t)
          v(y(i)+1,x(i)) = 1;
          n = n + 1;
        end
      end
      if x(i) < nx
        if ~v(y(i),x(i)+1) & (d(y(i),x(i)+1) > t)
          v(y(i),x(i)+1) = 1;
          n = n + 1;
        end
      end
    end
    m = m + 1;
  end
