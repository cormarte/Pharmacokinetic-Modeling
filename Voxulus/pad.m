% [Y,D] = PAD (X,V) pads 3D dataset 'X' to homogeneous size with value 'V'
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%
function [ y, d ] = pad (x, v)
  if nargin < 2
    v = 0;
  end
  nx = size (x, 3);
  ny = size (x, 2);
  nz = size (x, 1);
  n = max (size (x));
  ix = fix (1 + (n - nx) / 2);
  iy = fix (1 + (n - ny) / 2);
  iz = fix (1 + (n - nz) / 2);
  y = v * ones (n, n, n);
  y(iz:iz+nz-1,iy:iy+ny-1,ix:ix+nx-1) = x;
  d = [ n-nz n-ny n-nx ];
