% Y = CROP (X,ZI,YI,XI) crops 3D dataset 'X' with intervals
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%
function y = crop (x, zi, yi, xi)
  nx = size (x, 3);
  ny = size (x, 2);
  nz = size (x, 1);
  if isempty (zi)
    zi = [ 1 nz ];
  end
  if isempty (yi)
    yi = [ 1 ny ];
  end
  if isempty (xi)
    xi = [ 1 nx ];
  end
  y = x(zi(1):zi(2),yi(1):yi(2),xi(1):xi(2));
