%MYOCARD provide myocard related definitions
%   Y=MYOCARD(X) returns a list of myocardial definitions given by 'X'
%      'X' can be: 'segment_angles', 'segment_names', 'segment_radius'
%
%   The segment definitions stem from:
%
%   Cerqueira, M. D. et al., "Standardized Myocardial Segmentation and
%   Nomenclature for Tomographic Imaging of the Heart", Circulation, Vol.
%   105, pp. 539-542, 2002
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function y = myocard (x)
  if strcmp (x, 'segment_angles')
    y = cell (17, 1);
    y{1} = [ 60 120 ];
    y{2} = [ 120 180 ];
    y{3} = [ 180 240 ];
    y{4} = [ 240 300 ];
    y{5} = [ 300 360 ];
    y{6} = [ 0 60 ];
    y{7} = [ 60 120 ];
    y{8} = [ 120 180 ];
    y{9} = [ 180 240 ];
    y{10} = [ 240 300 ];
    y{11} = [ 300 360 ];
    y{12} = [ 0 60 ];
    y{13} = [ 45 135 ];
    y{14} = [ 135 225 ];
    y{15} = [ 225 315 ];
    y{16} = [ 0 45 315 360 ];
    y{17} = [ 0 360 ]; 
  elseif strcmp (x, 'segment_names')
    y = cell (17, 1);
    y{1} = 'basal anterior';
    y{2} = 'basal anteroseptal';
    y{3} = 'basal inferoseptal';
    y{4} = 'basal inferior';
    y{5} = 'basal inferolateral';
    y{6} = 'basal anterolateral';
    y{7} = 'mid anterior';
    y{8} = 'mid anteroseptal';
    y{9} = 'mid inferoseptal';
    y{10} = 'mid inferior';
    y{11} = 'mid inferolateral';
    y{12} = 'mid anterolateral';
    y{13} = 'apical anterior';
    y{14} = 'apical septal';
    y{15} = 'apical inferior';
    y{16} = 'apical lateral';
    y{17} = 'apex';
  elseif strcmp (x, 'segment_radius')
    y = cell (17, 1);
    y{1} = [ 0.75 1.0 ];
    y{2} = [ 0.75 1.0 ];
    y{3} = [ 0.75 1.0 ];
    y{4} = [ 0.75 1.0 ];
    y{5} = [ 0.75 1.0 ];
    y{6} = [ 0.75 1.0 ];
    y{7} = [ 0.5 0.75 ];
    y{8} = [ 0.5 0.75 ];
    y{9} = [ 0.5 0.75 ];
    y{10} = [ 0.5 0.75 ];
    y{11} = [ 0.5 0.75 ];
    y{12} = [ 0.5 0.75 ];
    y{13} = [ 0.25 0.5 ];
    y{14} = [ 0.25 0.5 ];
    y{15} = [ 0.25 0.5 ];
    y{16} = [ 0.25 0.5 ];
    y{17} = [ 0 0.25 ];
  else
    error ('unknown category ''%s''', x);
  end
