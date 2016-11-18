%POLAR_MAP create a polar map plot
%   POLAR_MAP(X[,S,L]) creates a polar map plot of the data 'X' (17 segment
%   values according to the myocardial tomographic imaging standard). 'S' can
%   define a range for the colormap (default is data min and max). 'L' is
%   a flag to show the segment names instead of the numbers.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function polar_map (x, s, l)
  hasError = size (x,2) > 1;
  if nargin < 2
    if hasError
      s = [ min(x(:,1)-x(:,2)) max(x(:,1)+x(:,2)) ];
    else
      s = [ min(x(:)) max(x(:)) ];
    end
  end
  if nargin < 3
    l = 0;
  end
  n = 25;
  m = 17;
  xp = zeros (2 * n, m);
  yp = zeros (2 * n, m);
  tx = zeros (m, 1);
  ty = zeros (m, 1);
  xpe = zeros (2 * n, 2 * m);
  ype = zeros (2 * n, 2 * m);

  name = myocard ('segment_names');
  sector = myocard ('segment_angles');
  r = myocard ('segment_radius');

  for i=1:m
    sector{i} = sector{i} * pi / 180.0;
    [ xp(:,i), yp(:,i) ] = create_sector (r{i}, sector{i}, n);
    w = sector{i};
    if length (w) > 2
      w = [ w(1)+w(3) w(2)+w(4) ];
    end
    [ tx(i), ty(i) ] = pol2cart (mean (w), mean (r{i}));
    if hasError
      dr = r{i}(2) - r{i}(1);
      r{i}(1) = r{i}(1) + dr / 3;
      r{i}(2) = r{i}(2) - dr / 3;
      w = sector{i};
      if length (w) > 2
        dw1 = w(2) - w(1);
        dw2 = w(4) - w(3);
        w = [ w(1)+dw1/2 w(2) w(3) w(3)+dw2/2 ];                
      else
        dw = w(2) - w(1);
        w = [ w(1)+dw/4 w(1)+dw/2 w(1)+dw/2 w(2)-dw/4 ];      
      end
      [ xpe(:,i), ype(:,i) ] = create_sector (r{i}, [ w(1) w(2) ], n);
      [ xpe(:,i+m), ype(:,i+m) ] = create_sector (r{i}, [ w(3) w(4) ], n);
    end
  end
  % colors
  axis ([-1 1 -1 1]);
  if hasError
    c = x(:,1)';
    fill (xp, yp, c);      
    hold on;
    c = [ x(:,1)-x(:,2); x(:,1)+x(:,2) ]';
    fill (xpe, ype, c);
    hold off;
  else
    c = x';
    fill (xp, yp, c);      
  end
  caxis (s);
  colorbar;
  if l
    for i=1:size(name,1)
      if l == 1
        label = num2str (i);
        text (tx(i), ty(i), label, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
      else
        label = name{i};
        text (tx(i), ty(i), label, 'HorizontalAlignment', 'center');
      end
    end
  end
  axis off;
  text (1.1, 1.05, 'MBF [ml/g/min]','HorizontalAlignment','center', 'FontWeight', 'bold');
  set (gca, 'FontSize', 14);
  set (gca, 'FontWeight', 'bold');

function [ xp, yp ] = create_sector (r, w, n)
  if nargin < 3
    n = 100;
  end
  xp = zeros (2 * n, 1);
  yp = zeros (2 * n, 1);
  if (length (w) > 2)
    w = [ w(1)+w(3) w(2)+w(4) ];
  end
  t = linspace (w(1), w(2), n)';
  if (r(1) > 0)
    t = linspace (w(1), w(2), n)';
    [ xp(1:n), yp(1:n) ] = pol2cart (t, r(1));
    [ xp(n+1:2*n), yp(n+1:2*n) ] = pol2cart (flipud (t), r(2));
  else
    t = linspace (w(1), w(2), 2 * n)';
    [ xp(1:2*n), yp(1:2*n) ] = pol2cart (t, r(2));
  end
