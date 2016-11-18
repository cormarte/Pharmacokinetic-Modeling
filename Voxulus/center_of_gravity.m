%CENTER_OF_GRAVITY computes the center of gravity of a 3D data set
%   CG = CENTER_OF_GRAVITY (D) computes the center of gravity of the 3D data
%   set 'D'
%
%   The indexing order of D is Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 02.12.2004
%

function cg = center_of_gravity (d)
  [ nz, ny, nx ] = size (d);
  cg = zeros (3, 1);
  for z=1:nz
    for y=1:ny
      for x=1:nx
        cg(1) = cg(1) + z * d(z,y,x);
        cg(2) = cg(2) + y * d(z,y,x);
        cg(3) = cg(3) + x * d(z,y,x);
      end
    end
  end
  cg = cg / sum (d(:));
