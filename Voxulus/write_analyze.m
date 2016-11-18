function write_analyze (fname, data, time, format, scale, pixdim, origin, verbosity)
%WRITE_ANALYZE write image data in Analyze format
%   WRITE_ANALYZE(FNAME,DATA,[TIME,FORMAT,PIXDIM,ORIGIN,VERBOSITY]) 
%      writes 3D or 4D data set 'DATA' to file 'FNAME'
%   Arguments:
%      FNAME: filename without extension
%      DATA: 3D or 4D data set
%      TIME: matrix with frame start times and lengths ([] = none)
%      SCALE: scaling factor (0 = automatic, default = 1.0)
%      FORMAT: 'uchar', 'int16', 'float32', 'float64' (int16)
%      PIXDIM: pixel dimensions (1,1,1)
%      ORIGIN: volume origin (0,0,0)
%      VERBOSITY: verbosity level (1)
%
% Frank Thiele, Alexander Fischer, Philips Research, 2004
%

if (nargin < 3)
  time = [];
end
if (nargin < 4)
  format = class (data);
end
if (nargin < 5)
  scale = 1.0;
end
if (nargin < 6)
  pixdim = ones (3, 1);
end
if (nargin < 7)
  origin = zeros (3, 1);
end
if (nargin < 8)
  verbosity = 1;
end

[ pathstr, name, ext ] = fileparts (fname);
if (~strcmp (ext, '.img') && ~strcmp(ext, '.hdr') && ~strcmp (ext, '') )
  error('unknown filename extension %s for analyze format!', ext);
end

headerFileName = fullfile(pathstr, [name '.hdr']);
imageFileName = fullfile(pathstr, [name '.img']);
timeFileName = fullfile(pathstr, [name '.time']);

% init and set variables for header
dataDims = size (data);
numDim = length (dataDims);
dim = zeros (5, 1);
pixdim = zeros (4,1);

dim(1) = numDim;
dim(2) = dataDims(3); % x
dim(3) = dataDims(2); % y
dim(4) = dataDims(1); % z
if (numDim == 4)
  dim(5) = dataDims(4); % t
else
  dim(5) = 1;
end
vox_units = 'mm';

% determine min,max:
dataMin = min (data(:));
dataMax = max (data(:));
if (verbosity)
  fprintf ('min=%f max=%f\n', dataMin, dataMax);
end
scale = 1.0;
glmin = dataMin;
glmax = dataMax;

if (strcmp (format, 'bool'))
  datatype = 2;
  bitpix = 8;
  scale = 1.0;
  glmin = 0.0;
  glmax = 1.0;
  format = 'uint8';
elseif (strcmp (format, 'uint8'))
  datatype = 2;
  bitpix = 8;
  if scale == 0
    % scale to max=127
    scale = max (abs (dataMax), abs (dataMin)) / 127.0;  
  end
  glmin = 0;
  glmax = 127;
  format = 'uint8';
elseif (strcmp (format, 'int16'))
  datatype = 4;
  bitpix = 16;
  if scale == 0
    % scale to max=32767
    scale = max (abs (dataMax), abs (dataMin)) / 32767.0;
  end
  glmin = 0;
  glmax = 32767;
elseif (strcmp (format, 'float32') || strcmp (format, 'float'))
  datatype = 16;
  bitpix = 32;
  format = 'float32';
elseif (strcmp (format, 'float64') || strcmp (format, 'double'))
  datatype = 64;
  bitpix = 64;
  format = 'float64';
else
  error ('unknown format ''%s''', format);
end

vox_offset = 0;
orient = 0;
headerSize = 348;

if (verbosity)
  fprintf ('scale=%f\n', scale);
end

% open and write actual image file:
fp = fopen (imageFileName, 'w');
if (fp == -1)
  error (['Failed to open file ' imageFileName]);
end

if (verbosity)
  dataSize = dim(5)*dim(4)*dim(3)*dim(2);
  fprintf ('writing %d %s voxels to file ''%s''\n', dataSize, class (data), imageFileName);
end

n = 0;
for t=1:dim(5)
  if (verbosity>1)
    fprintf ('writing frame %d...\n', t);
  end
  for z=1:dim(4)
    for y=1:dim(3)
      %   write all dim(2) x-values at once
      if (numDim == 4)
        values = data (z,y,:,t) / scale;
      else
        values = data (z,y,:) / scale;
      end
      n = n + fwrite (fp, values, format);
    end
  end
end
fclose(fp);

if (verbosity>0)
  fprintf('wrote %u voxel values.\n',n);
end

% open header file
if (verbosity>0)
    disp (['writing header file ' headerFileName '...']);
end
fp = fopen (headerFileName, 'w');
if (fp == -1)
  error (['Failed to open header file ' headerFileName]);
end

% write 348 byte header file
% init with zero:
fwrite(fp,zeros(headerSize,1),'int8');
% start from beginning of file again:
fseek(fp,0,'bof');
fwrite(fp,headerSize,'int32');
fseek(fp,38,'bof');
fprintf(fp,'r');
fseek(fp,40,'bof');
fwrite(fp,dim,'int16');
fseek(fp,56,'bof');
fwrite(fp,vox_units,'uchar');
fseek(fp,70,'bof');
fwrite(fp,datatype,'int16');
fwrite(fp,bitpix,'int16');
fseek(fp,76,'bof');
fwrite(fp,pixdim,'float32');
fseek(fp,108,'bof');
fwrite(fp,vox_offset,'float32');
fwrite(fp,scale,'float32');
fseek(fp,140,'bof');
fwrite(fp,glmax,'int32');
fwrite(fp,glmin,'int32');
fseek(fp,253,'bof');
fwrite(fp,origin(1),'int16');
fwrite(fp,origin(2),'int16');
fwrite(fp,origin(3),'int16');
fclose(fp);

if (~isempty (time))
  % open time file
  if (verbosity>0)
    disp (['writing time file ' timeFileName '...']);
  end
  fp = fopen (timeFileName, 'w');
  if (fp == -1)
    error (['Failed to open time file ' timeFileName]);
  end
  timeSize = size (time);
  if (timeSize(1) ~= 2)
    time = time';
  end
  fprintf (fp, '# Acquisition times (start frame length) in seconds\n');
  fprintf (fp, '%d # number of acqusitions\n', size(time,2) );
  fprintf (fp, '%g %g\n', time);
  fclose (fp);
end
