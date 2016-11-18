% [ R, C, B ] = READ_COS_PAR (FILE) read transformation parameters from
%    Philips Gemini Cardiac Oblique Section software
%    'FILE' is filename with the exported parameters
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%
function [ r, c, b ]= read_cos_par (file)
  fp = fopen (file, 'r');
  if fp == -1
    error ('cannot open parameter file ''%s''', file);
  else
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    r = fscanf (fp, '%d', 3); % ???
    p = fscanf (fp, '%d', 3); % ???
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    r = fscanf (fp, '%d', 2); % rotation angle
    c = fscanf (fp, '%d', 3); % rotation center
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    fgetl (fp);
    b = fscanf (fp, '%d', 6); % cropping boundaries
    fclose (fp);
  end
