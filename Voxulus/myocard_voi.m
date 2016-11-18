%MYOCARD_VOI create a volume of interest of the myocardium
%   [V, VS] = MYOCARD_VOI (D, VI, [Q, [VERBOSITY]) creates a myocardial VOI for 
%      the 3D dataset 'D' based on the (initial) VOI 'VI'. 'VI' must select
%      the centerlines of the two long axis views and the short axis plane
%      must be the x-y plane!
%      The VOI is constructed slice by slice along the z-axis. An ellipse
%      is constructed that minimizes the distance to the currently selected
%      pixels. Afterwards the ellipse is grown where 'Q' defines the 
%      percentage of the data maximum used as threshold.
%
%      The indexing order of D, VI and V is Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function [ v, vs ] = myocard_voi (d, vi, q, verbosity)
  if nargin < 4
    verbosity = 0;
  end
  if nargin < 3
    q = 0.6;
  end
  if size (vi) ~= size (d)
    error ('incompatible data %dx%dx%d and VOI %dx%dx%d\n', ...
      size(d,1), size(d,2), size(d,3), size(vi,1), size(vi,2), size(vi,3));
  end
  [ nz, ny, nx ] = size (vi);
  v = vi;
  thresh = q * max (d(:));
  for z=1:nz
    nv = nnz (v(z,:,:));
    if nv < 2
      v(z,:,:) = 0;
    elseif nv == 2
      s = squeeze(v(z,:,:));
      [ y, x ] = find (s > 0);
      if (x(1) == x(2)) | (y(1) == y(2))
        v(z,:,:) = 0;
      end
    end
  end
  apex = 1;
  while ~nnz (v(apex,:,:))
    apex = apex + 1;
  end
  base = nz;
  while ~nnz (v(base,:,:))
    base = base - 1;
  end
  h = base - apex;
  c = center_of_gravity (v);
  if verbosity
    fprintf ('apex = %g(%d) base = %g(%d)\n', apex, nnz(v(apex,:,:)), base, nnz(v(base,:,:)));
    fprintf ('center(y,x) = %g,%g\n', c(2), c(3));
  end
  for z=base+1:base+2
    if z <= nz
      v(z,:,:) = v(z-1,:,:);
    end
  end
  for z=nz:-1:round(apex+0.1*h)+1
    s = squeeze(v(z,:,:));
    [ y, x ] = find (s > 0);
    n = size (x, 1);
    if verbosity
      fprintf ('z=%d: %d points\n', z, n);
    end
    p = cluster_points ([ x, y ], 4);
    x = p(:,1);
    y = p(:,2);
    n = size (p, 1);
    if verbosity > 1
      figure (3);
      clf;
      scatter (x, y);
      axis ([ 1 nx 1 ny ]);
      title (sprintf ('z=%g', z));
      hold on;
      scatter (c(3), c(2), '+');
    end
    if n >= 2
      %cx = mean (x);
      %cy = mean (y);
      cx = c(3);
      cy = c(2);
      w = atan2 (y - cy, x - cx) * 180 / pi;
      iw = find(w < 0);
      w(iw) = w(iw) + 360;
      [ w, i ] = sort (w);
      p1x = x(i(n));
      p1y = y(i(n));
      w1 = w(n);
      r = [ Inf, Inf ];
      for j=1:n
        p2x = x(i(j));
        p2y = y(i(j));
        w2 = w(j);
        cp = get_center ([ p1x p1y ], [ p2x p2y ], [ cx cy ]);
        rx = max ([ abs(p1x-cp(1)) abs(p2x-cp(1)) ]);
        ry = max ([ abs(p1y-cp(2)) abs(p2y-cp(2)) ]);
        w1 = atan2 (p1y - cp(2), p1x - cp(1)) * 180 / pi;
        w2 = atan2 (p2y - cp(2), p2x - cp(1)) * 180 / pi;
        if verbosity > 1
          fprintf ('  %g,%g -> %g,%g c=%g,%g rx=%g ry=%g w=[%g %g]\n', p1x, p1y, p2x, p2y, cp(1), cp(2), rx, ry, w1, w2);
        end
        if w2 < 0 & w2 < w1
          w2 = w2 + 360;
        end
        if abs (w2 - w1) < 100
          r = [ min([r(1),rx]) min([r(2),ry]) ];
          if verbosity > 1
            pw = linspace (w1, w2, 100) * pi / 180;
            plot (cp(1) + rx * cos (pw), cp(2) + ry * sin (pw), 'b');
            pause;
          end
          v(z,:,:) = setEllipse (squeeze(v(z,:,:)), cp(1), cp(2), rx, ry, [ w1 w2 ]);
        end
        p1x = p2x;
        p1y = p2y;
        w1 = w2;
      end
    else
      v(z,:,:) = 0;
    end
    if verbosity
      fprintf ('  %d voxels\n', nnz (v(z,:,:)));
    end
  end
  %
  % now the apex
  % 
  cx = c(3);
  cy = c(2);
  cz = apex + 0.1 * h;
  rx = r(1);
  ry = r(2);
  rz = 0.1 * h;
  if verbosity
    fprintf ('apex c=%g,%g,%g r=%g,%g,%g\n', cx, cy, cz, rx, ry, rz);
  end
  for phi=0:pi/20:pi/2
    for theta=0:pi/50:2*pi
      x = round (cx + rx * cos (theta) * sin (phi));
      y = round (cy + ry * sin (theta) * sin (phi));
      z = round (cz - rz * cos (phi));
      v(z,y,x) = 1;
    end
  end
  if nargout > 1
    vs = v;
  end
  for z=apex:base+2
    if z < nz
      v(z,:,:) = growROI (squeeze(v(z,:,:)), squeeze(d(z,:,:)));
      v(z,:,:) = growROI (squeeze(v(z,:,:)), squeeze(d(z,:,:)), thresh);
    end
  end
  for z=apex-1:-1:apex-2
    if z > 0
      i = find ((v(z+1,:,:) > 0) & (d(z,:,:) >= thresh));
      s = squeeze (v(z,:,:));
      s(i) = 1;
      v(z,:,:) = s;
    end
  end
  v(smooth3 (v, 'gaussian')>=0.5) = 1;

function c = get_center (p1, p2, cr)
  c1 = p1;
  if p1(1) == p2(1)
    c = [ p1(1) mean([p1(2),p2(2)]) ];
  elseif p1(2) == p2(2)
    c = [ mean([p1(1),p2(1)]) p1(2) ];
  else
    c1(1) = p2(1);
    c = c1;
    d1 = (c1 - cr) * (c1 - cr)';
    c2 = p2;
    c2(1) = p1(1);
    d2 = (c2 - cr) * (c2 - cr)';
    if d2 < d1
      c = c2;
    end
  end

function v = setEllipse (s, cx, cy, rx, ry, w)
  if nargin < 6
    w = [ 0 360 ];
  end
  [ ny, nx ] = size (s);
  v = s;
  for phi=w(1):1:w(2)
    angle = phi * pi / 180;
    i = round (cx + rx * cos (angle));
    i = min ([ i nx ]);
    i = max ([ 1 i ]);
    j = round (cy + ry * sin (angle));
    j = min ([ j ny ]);
    j = max ([ 1 j ]);
    v(j,i) = 1;
  end
