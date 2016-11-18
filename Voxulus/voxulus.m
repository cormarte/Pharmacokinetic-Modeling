function varargout = voxulus (varargin)
%  RET=VOXULUS(ARGS) 
%    is the main interface to the Voxulus kernel
%    the global variable VOXULUS_VERSION selects different versions of the
%    engine (opt,dbg). The opt is used as default.
%
% Alexander Fischer, Philips Research, 2004
%

global VOXULUS_VERSION;

if length (VOXULUS_VERSION)
  vox = strcat ('voxulus_', VOXULUS_VERSION);
else
  vox = 'voxulus_opt';
end
cmd = '';
if nargout
  cmd = '[ varargout{1}';
  for k = 2:nargout
    tmp = sprintf (', varargout{%d}', k);
    cmd = strcat (cmd, tmp);
  end
  cmd = strcat (cmd, ' ] = ');
end
tmp = sprintf ('%s (', vox);
cmd = strcat (cmd, tmp);
for k = 1:length (varargin)
  if (k > 1)
    cmd = strcat (cmd, ', ');
  end
  tmp = sprintf ('varargin{%d}', k);
  cmd = strcat (cmd, tmp);
end
cmd = strcat (cmd, ');');
%eval (cmd);
