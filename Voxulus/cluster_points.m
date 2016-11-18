%CLUSTER_POINTS cluster a point list using distance based merging
%   PO = CLUSTER_POINTS (PI[,T]) tries to to cluster a list of points 'PI' to
%     by merging all points whose mutual distance is below the threshold 'T'.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 06/2005
%

function po = cluster_points (pi, t)
  if nargin < 3
    t = 2.0;
  end
  ni = size (pi, 1);
  po = zeros (size (pi));
  d = zeros (ni, ni);
  for i = 1:ni
    for j = 1:i-1
      d(i,j) = sqrt ((pi(i,1) - pi(j,1))^2 + (pi(i,2) - pi(j,2))^2);
    end
  end
  k = 0;
  for i = 1:ni
    if d(i,i) > -1
      k = k + 1;
      po(k,:) = pi(i,:);
      m = 1;
      for j = i+1:ni
        if (d(j,i) > 0) & (d(j,i) <= t)
          %fprintf ('merging point %d and %d m=%d\n', i, j, m);
          po(k,:) = (m * po(k,:) + pi(j,:)) / (m + 1);
          d(j,:) = -1;
          m = m + 1;
        end
      end
    end
  end
  po = po(1:k,:);
