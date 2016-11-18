%GET_TAC time activity curve extraction
%   [Y,S] = GET_TAC (X,VOI) extracts a TAC from a 3D or 4D data set X over the 
%   volume of interest VOI
%
%     X is ZxYxX[xT] numeric, VOI is ZxYxX logical
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 02.12.2004
%

function [ y, s ] = get_tac (x, voi)
dimX = size (x);
dimVoi = size (voi);
if dimX(1) ~= dimVoi(1) || dimX(2) ~= dimVoi(2) || dimX(3) ~= dimVoi(3)
  error ('incompatible dimensions');
end
if ndims (x) < 3
  error ('must have at least 3D data');
elseif ndims (x) == 3
  y = 0;
  s = 0;
  n = 0;
  for z=1:dimX(1)
    mask = voi(z,:,:);
    if nnz (mask)
      data = x(z,:,:);
      y = y + sum (double (data(mask)));
      s = s + sum (double (data(mask)).^2);
      n = n + nnz (double (mask));
    end
  end
  y = y / n;
  s = sqrt (s / n - y * y);
else
y = zeros (dimX(4), 1);
s = zeros (dimX(4), 1);
if nnz (double (voi(:)))
  for t=1:dimX(4)
    n = 0;
    for z=1:dimX(1)
      mask = squeeze (voi(z,:,:));
      if nnz (double (mask))
        data = x(z,:,:,t);
        y(t) = y(t) + sum (double (data(mask)));
        s(t) = s(t) + sum (double (data(mask)).^2);
      end
      n = n + nnz (double (mask));
    end
    y(t) = y(t) / n;
    s(t) = sqrt (s(t) / n - y(t) * y(t));
  end
end
end
