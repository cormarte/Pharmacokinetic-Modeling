%MYOCARD_SEGMENTS create standard myocardial segments from a full VOI
%   S = MYOCARD_SEGMENTS (V,A,C,VERBOSITY) creates a segmented myocardial 
%      VOI from the full myocardial VOI 'V'.
%      'A' is the angle to the anterior-anteroseptal junction to the right
%          ventricle (default: 120.0)
%      'C' denotes the center of myocardium that serves as reference point
%          for all angles (default is the center of 'V')
%
%      The indexing order of the VOIs is Z,Y,X.
%
%   requires: myocard.m
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function [ s, name ] = myocard_segments (v, a, c, verbosity)
  nz = size (v, 1);
  ny = size (v, 2);
  nx = size (v, 3);
  if nargin < 2
    a = 120.0;
  end
  if nargin < 3
    c = [];
  end
  if isempty (c)
    c = [ (1+nz)/2, (1+ny)/2, (1+nx)/2 ];
  end
  if nargin < 4
    verbosity = 1;
  end
  s = zeros (size (v), 'uint8');

  name = myocard ('segment_names');
  sector = myocard ('segment_angles');
  %
  % rotate according to right ventricle junction and 
  % normalize to [0,360]
  %
  for i = 1:16
    sector{i} = sector{i} + a - 120.0;
    tmp = [];
    for j=1:2:length(sector{i})
      w1 = sector{i}(j);
      w2 = sector{i}(j+1);
      if w2 < 0
        w1 = w1 + 360;
        w2 = w2 + 360;
        sector{i}(j:j+1) = [ w1 w2 ];
      end
      if w1 < 0
        sector{i}(j:j+1) = [ 0 w2 ];
        tmp = [ tmp w1+360 360 ];
      end
      if w2 > 360
        sector{i}(j:j+1) = [ 0 w2-360 ];
        tmp = [ tmp w1 360 ];
      end
    end
    sector{i} = [ sector{i} tmp ];
    % check for invalid sub-sectors
    tmp = [];
    for j=1:2:length(sector{i})
      if sector{i}(j) < sector{i}(j+1)
        tmp = [ tmp sector{i}(j) sector{i}(j+1) ];
      end
    end
    sector{i} = tmp;
    if verbosity > 1
      fprintf ('%s', name{i});
      fprintf (' %.0f-%.0f', sector{i});
      fprintf ('\n');
    end
  end

  apex = 1;
  while nnz (v(apex,:,:)) < 10
    apex = apex + 1;
  end
  
  base = nz;
  while nnz (v(base,:,:)) < 10
    base = base - 1;
  end
  h = base - apex;
  
  d = [ apex round(apex+0.1*h) round(apex+0.4*h) round(apex+0.7*h) ];

  for z = apex:base
    t = squeeze(s(z,:,:));
    [ y, x ] = find(squeeze(v(z,:,:)) > 0);
    w = (atan2 (y - c(2), x - c(3))) * 180 / pi;
    iw = find(w < 0);
    w(iw) = w(iw) + 360;
    k = find(squeeze(v(z,:,:)) > 0);
    if z > d(4)
      %
      % basal
      %
      for i = 1:6
        for j = 1:2:length(sector{i})
          w1 = sector{i}(j);
          w2 = sector{i}(j+1);
          t(k(find (((w >= w1) & (w < w2))))) = i;
        end
      end
    elseif z > d(3)
      %
      % mid
      %
      w1 = 0.0;
      for i=7:12
        for j = 1:2:length(sector{i})
          w1 = sector{i}(j);
          w2 = sector{i}(j+1);
          t(k(find (((w >= w1) & (w < w2))))) = i;
        end
      end
    elseif z > d(2)
      %
      % apical
      %
      for i=13:16
        for j = 1:2:length(sector{i})
          w1 = sector{i}(j);
          w2 = sector{i}(j+1);
          t(k(find (((w >= w1) & (w < w2))))) = i;
        end
      end
    else
      %
      % apex
      %
      t(k) = 17;
    end
    s(z,:,:) = t;
  end
  
  if verbosity  
    rt = 0;
    nt = 0;
    fprintf ('segment                 N     r\n');
    fprintf ('================================\n');
    for i=1:17
      k = find (s == i);
      n = size (k, 1);
      r = 100.0 * n / nnz (v);
      fprintf ('%-20s %4d %5.1f\n', name{i}, n, r); 
      nt = nt + n;
      rt = rt + r;
    end
    fprintf ('================================\n');
    fprintf ('MYOCARD              %4d %5.1f\n', nt, rt);
  end
