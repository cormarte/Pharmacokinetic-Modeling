%INTERPOLATE_IMAGE interpolates a 2D image
%   Y = INTERPOLATE_IMAGE (X,NY[,NX]) interpolates image X to resolution NYxNX 
%     
%   A. Fischer, Philips Research Labs, Aachen, Germany, 02/2005
%

function y = interpolate_image (x, ny, nx)
  if nargin < 2
    ny = 256;
  end
  if nargin < 3
    nx = ny;
  end
  if nx & ny
    s = size (x);
    xo = 1:s(2);
    yo = 1:s(1);
    xi = linspace (1, s(2), nx);
    yi = linspace (1, s(1), ny);
    if length (s) == 3
      y = zeros (ny, nx);
      for i=1:s(3)
        y(:,:,i) = interp2 (xo, yo', squeeze (x(:,:,i)), xi, yi', 'bicubic'); 
      end
    else
      y = interp2 (xo, yo', x, xi, yi', 'bicubic'); 
    end
  else
    y = x;
  end