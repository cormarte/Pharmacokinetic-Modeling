%ROTATE3 rotate a 3D dataset
%   Y = ROTATE3 (D,THETA,PHI,PSI,C) rotates the 3D dataset D by angles 
%     THETA, PHI, PSI around axis x, y and z. The rotation center is
%     optionally given by C (default is the center of the cube).
%     
%     D is XxYxZ
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 01/2005
%

function y = rotate3 (d, theta, phi, psi, c)
  nx = size (d, 1);
  ny = size (d, 2);
  nz = size (d, 3);
  if nx ~= ny || nx ~= nz || ny ~= nz
    error ('can only rotate homogeneous volumes');
  end
  dmin = min (d(:));
  %
  % rotation matrices, Slomka, "Medical Imaging Volume 2", p. 453
  %
  wx = pi * theta / 180.0;
  wy = pi * phi / 180.0;
  wz = -pi * psi / 180.0;
  rx = [ 1, 0, 0; 0, cos(wx), -sin(wx); 0, sin(wx), cos(wx) ];
  ry = [ cos(wy), 0, sin(wy); 0, 1, 0; -sin(wy), 0, cos(wy) ];
  rz = [ cos(wz), -sin(wz), 0; sin(wz), cos(wz), 0; 0, 0, 1 ];
  r = rz * ry * rx;
  vc = [ (1+nx)/2, (1+ny)/2, (1+nz)/2 ];
  if nargin < 5
    c = vc;
  end
  %fprintf ('rotation rx=%.0f ry=%.0f rz=%.0f around x=%.0f y=%.0f z=%.0f\n', theta, phi, psi, c(1), c(2), c(3));
  %
  % create nx*ny*nz by 3 arrary with all volume coordinates
  %
  v = zeros (nx * ny * nz, 3);
  i = 1;
  for x=1:nx
    for y=1:ny
      for z=1:nz
        v(i,1) = x;
        v(i,2) = y;
        v(i,3) = z;
        i = i + 1;
      end
    end
  end
  %
  % move the given center to (0,0,0)
  %
  v(:,1) = v(:,1) - c(1);
  v(:,2) = v(:,2) - c(2);
  v(:,3) = v(:,3) - c(3);
  %
  % rotation
  %
  v = v * r;
  %
  % move (0,0,0) to the cube center
  %
  v(:,1) = v(:,1) + vc(1);
  v(:,2) = v(:,2) + vc(2);
  v(:,3) = v(:,3) + vc(3);
   
  %
  % prepare index arrays for 3D interpolation
  %
  x = zeros (size (d));
  y = zeros (size (d));
  z = zeros (size (d));
  l = 1;
  for i=1:nx
    for j=1:ny
      for k=1:nz
        x(i,j,k) = v(l,1);
        y(i,j,k) = v(l,2);
        z(i,j,k) = v(l,3);
        l = l + 1;
      end
    end
  end
  %
  % 3D interpolation
  %
  y = interp3 (d, x, y, z, 'cubic', dmin);
  %
  % somehow x and y are flipped
  %
  y = permute (y, [ 2 1 3 ]);
