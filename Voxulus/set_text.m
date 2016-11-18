function set_text (x, y, s)
%SET_TEXT draw text in axes at normalized coordinates
%   SET_TEXT (X, Y, S) draws text 's' in the current axes at normalized 
%      ([0,1]) coordinates 'x' and 'y'
%
% A. Fischer, 12/2004
%
  xlim = get (gca, 'XLim');
  ylim = get (gca, 'YLim');
  x = xlim(1) + (xlim(2) - xlim(1)) * x;
  y = ylim(1) + (ylim(2) - ylim(1)) * y;
  text (x, y, s);
