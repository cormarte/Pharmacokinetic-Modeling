function [data, frame, info] = read_analyze (fname, verbosity)
%READ_ANALYZE read image data in analyze format
%   [DATA,FRAME,INFO] = READ_ANALYZE (FNAME,VERBOSITY) read Analyze image data from
%      file 'FNAME'
%   Arguments:
%      FNAME: filename without extension
%      DATA: output 3D or 4D data set
%      FRAME: matrix with frame start times and lengths ([] = none)
%      INFO: struct with additional file information (voxel dim, ...)
%
% Frank Thiele, Alexander Fischer, Philips Research, 2004
%

if nargin == 0
  error ('no filename');
else
if nargin < 2
  verbosity = 1;
end

[ pathstr, name, ext ] = fileparts (fname);
if (~strcmp (ext, '.img') && ~strcmp(ext, '.hdr') && ~strcmp (ext, '') )
  error('unknown filename extension %s for analyze format!', ext);
end

headerFileName = fullfile(pathstr, [name '.hdr']);
imageFileName = fullfile(pathstr, [name '.img']);
timeFileName = fullfile(pathstr, [name '.time']);

% open header file
fp = fopen (headerFileName, 'r');
if fp == -1
  error (['Failed to open header file ' headerFileName]);
end
if (verbosity)
  disp (['reading header file ' headerFileName '...']);
end

info.fileName = headerFileName;

% read parts of the 348 byte header file
headerSize = fread (fp, 1, 'int32');
if (headerSize ~= 348)
  error('invalid header size %d instead of 348', headerSize);
end
% move to the image_dimension part of the header
fseek (fp, 40, 'bof');
dim = fread (fp, 5, 'int16');
if (verbosity)
  fprintf ('  %dD dimensional data\n', dim(1));
end
if ( dim(1) <0 || dim(1) >15 )
  error ('little/big endian mismatch not taken care of yet!');
end
fseek(fp,56,'bof');
vox_units=fread(fp,4,'uchar');
fseek(fp,70,'bof');
datatype=fread(fp,1,'int16');
bitpix=fread(fp,1,'int16');

if datatype == 2
  format = 'uint8';
  data = uint8 (zeros (dim(4), dim(3), dim(2), dim(5)));
elseif datatype == 4
  format = 'int16';
  data = int16 (zeros (dim(4), dim(3), dim(2), dim(5)));
elseif datatype == 8
  format = 'int32';
  data = int32 (zeros (dim(4), dim(3), dim(2), dim(5)));
elseif datatype == 16
  format = 'float32';
  data = single (zeros (dim(4), dim(3), dim(2), dim(5)));
elseif datatype == 64
  format = 'float64';
  data = zeros (dim(4), dim(3), dim(2), dim(5));
else
  error ('unsupported datatype ''%d''',datatype);
end

fseek(fp,56,'bof');
info.vox_units=num2str(fread(fp,4,'4*char=>char')');
info.cal_units=num2str(fread(fp,8,'8*char=>char')');
fseek(fp,80,'bof');
info.pixdim=fread(fp,3,'float');
fseek(fp,108,'bof');
info.vox_offset(1)=fread(fp,1,'float');
info.scale=fread(fp,1,'float');
info.vox_offset(2)=fread(fp,1,'float');
info.vox_offset(3)=fread(fp,1,'float');
info.cal_min=fread(fp,1,'float');
info.cal_max=fread(fp,1,'float');
info.compressed=fread(fp,1,'float');
info.verified=fread(fp,1,'float');
info.glmax=fread(fp,1,'int32');
info.glmin=fread(fp,1,'int32');
info.description=num2str(fread(fp,80,'80*char=>char')');
fseek(fp,253,'bof');
info.originator=num2str(fread(fp,10,'10*char=>char')');
info.generated=num2str(fread(fp,10,'10*char=>char')');
info.scannum=num2str(fread(fp,10,'10*char=>char')');
info.patient_id=num2str(fread(fp,10,'10*char=>char')');
info.exp_date=num2str(fread(fp,10,'10*char=>char')');
info.exp_time=num2str(fread(fp,10,'10*char=>char')');
fclose(fp);

if (verbosity>1)
    fprintf('sizeof_hdr: <%d> \n', headerSize);
    for i=1:5
        fprintf('dim[%d]: <%d> \n', i, dim(i));
    end
    fprintf('datatype:  <%d> \n', datatype);
    fprintf('bitpix:    <%d> \n', bitpix);
    for i=1:3
        fprintf('pixdim[%d]: <%6.4f> \n',i, pixdim(i));
    end
    fprintf('SPM scale:  <%6.20f> \n', scale);
    fprintf('glmax:      <%d> \n', glmax);
    fprintf('glmin:      <%d> \n', glmin);
end

% open and read actual image file:
fp = fopen(imageFileName,'r');
if fp==-1
    error(['Failed to open image file ' imageFileName]);
end

if (verbosity)
    disp (['reading image file ' imageFileName '...']);
end

n=0;
volume=dim(5)*dim(4)*dim(3)*dim(2);
if (verbosity>0)
  fprintf ('  %d %s voxels from ''%s''\n', volume,format,imageFileName);
end

for t=1:dim(5)
  if (verbosity>1)
    fprintf('reading time frame %d...\n',t);
  end
  for z=1:dim(4)
    for y=1:dim(3)
      data(z,y,:,t) = fread (fp, dim(2), format);
      n=n+dim(2);
    end
  end
end
fclose(fp);

if (dim(1) > 3) && (nargout > 1)
  fp = fopen(timeFileName,'r');
  if fp==-1
    error(['Failed to open time file ' timeFileName]);
  end
  disp (['reading image file ' timeFileName '...']);
  fgetl (fp);
  n = fscanf (fp, '%d', 1);
  fgetl (fp);
  if n ~= dim(5)
    error ('time dimension mismatch ''%d'' versus ''%d''', n, dim(5));
  end
  [ frame, count ] = fscanf (fp, '%f');
  if count ~= 2*n
    error ('wrong number of items (%d) in time file', count / 2);
  end
  frame = reshape (frame, 2, n)';
  fclose (fp);
else
  frame = [];
end

end
