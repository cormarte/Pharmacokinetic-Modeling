function out = process4d (d,cmd,arg1,arg2,arg3)
% PROCESS4D(X) process 4D data set X according to command CMD with the additional
%    command parameters ARG1, ARG2, ...
%    CMD can be
%    'crop'   : arg1, arg2, arg3 define the z, y and x interval
%    'flipxy' : 
%    'flipxz' :
%    'flipyz' :
%    'interp' : arg1 contains new Z,Y,X dimension
%    'pad'    : zero padding to homogeneous volume size (=max(Z,Y,X))
%    'revx'   :
%    'revy'   :
%    'revz'   :
%    'rotate' : rotate by angles in arg1 with center in arg2
%
%    A. Fischer, (C) Philips Research, 2005
%
  fprintf ('%s on %dx%dx%dx%d\n', cmd, size(d,1), size(d,2), size(d,3), size(d,4));
  if strcmp (cmd, 'crop')
    out = zeros ([ length(arg1) length(arg2) length(arg3) size(d,4) ]);
  elseif strcmp (cmd, 'flipxz')
    out = zeros ([ size(d,3) size(d,2) size(d,1) size(d,4) ]);
  elseif strcmp (cmd, 'flipyz')
    out = zeros ([ size(d,2) size(d,3) size(d,1) size(d,4) ]);
  elseif strcmp (cmd, 'flipxy')
    out = zeros ([ size(d,1) size(d,3) size(d,2) size(d,4) ]);
  elseif strcmp (cmd, 'interp')
    out = zeros ([ arg1 size(d,4) ]);
  elseif strcmp (cmd, 'pad')
    n = max (size (squeeze (d(:,:,:,1))));
    out = zeros ([ n n n size(d,4) ]);
  else
    out = zeros (size (d));
  end
  dataMin = min (d(:));
  for t=1:size(d,4)
    fprintf ('t=%d\n', t);
    v = squeeze (d(:,:,:,t));
    if strcmp (cmd, 'crop')
      v = v(arg1,arg2,arg3);
    elseif strcmp (cmd, 'flipxy')
      v = permute (v, [ 1 3 2 ]);
    elseif strcmp (cmd, 'flipyz')
      v = permute (v, [ 2 3 1 ]);
    elseif strcmp (cmd, 'flipxz')
      v = permute (v, [ 3 2 1 ]);
    elseif strcmp (cmd, 'interp')
      nx = size(v,3);
      ny = size(v,2);
      nz = size(v,1);
      [ x, y, z ] = ndgrid (1:nx, 1:ny, 1:nz);
      [ xi, yi, zi ] = meshgrid (linspace (1, nx, arg1(3)), linspace (1, ny, arg1(2)), linspace (1, nz, arg1(1)));
      tmp = interp3 (x, y, z, permute (v, [ 3 2 1 ]), xi, yi, zi, 'cubic', dataMin);
      v = permute (tmp, [ 3 2 1 ]);
    elseif strcmp (cmd, 'pad')
      nx = size(v,3);
      ny = size(v,2);
      nz = size(v,1);
      x = fix (1 + (n - nx) / 2);
      y = fix (1 + (n - ny) / 2);
      z = fix (1 + (n - nz) / 2);
      tmp = zeros (n, n, n);
      tmp(z:z+nz-1,y:y+ny-1,x:x+nx-1) = v;
      v = tmp;
    elseif strcmp (cmd, 'revx')
      nx = size(v,3);
      for z = 1:size(v,1)
        for y = 1:size(v,2)
          v(z,y,:) = v(z,y,nx:-1:1);
        end
      end
    elseif strcmp (cmd, 'revy')
      ny = size(v,2);
      for z = 1:size(v,1)
        for x = 1:size(v,3)
          v(z,:,x) = v(z,ny:-1:1,x);
        end
      end
    elseif strcmp (cmd, 'revz')
      nz = size(v,1);
      for y = 1:size(v,2)
        for x = 1:size(v,3)
          v(:,y,x) = v(nz:-1:1,y,x);
        end
      end
    elseif strcmp (cmd, 'rotate')
      v = permute (rotate3 (permute (v, [ 3 2 1 ]), arg1(1), arg1(2), arg1(3), arg2), [ 3 2 1]);
    else
      error ('unknown command ''%s''', cmd);
    end
    out(:,:,:,t) = v;
  end
