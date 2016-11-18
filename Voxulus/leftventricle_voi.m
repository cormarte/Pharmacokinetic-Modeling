%LEFTVENTIRCLE_VOI create a volume of interest of the left ventricle
%   V = LEFTVENTRICLE_VOI (D, VI) creates a left ventricle VOI for the 3D data 
%      set 'D' from the myocardial VOI 'VI'.
%
%      The indexing order of the data and the VOIs is Z,Y,X.
%
%   A. Fischer, Philips Research Labs, Aachen, Germany, 04/2005
%

function v = leftventricle_voi (d, vi, nIterations)
  if ( nargin < 3 )
      nIterations = 1;
  end
  if ( nIterations < 1 )
      nIterations = 1;
  end
  v = logical (zeros (size (vi)));
  [ nz, ny, nx ] = size (vi);
  cg = center_of_gravity (vi);
  %
  % get the height of the myocardial VOI
  %
  h = 0;
  for z=1:nz
    if nnz (vi(z,:,:))
      h = h + 1;
    end
  end
  %
  % move the center somewhat up to the basal plane
  %
  cg(1) = cg(1) + h / 5;
  cg = round (cg);
  %
  % set the height of the blood pool VOI
  %
  ch = round (h / 10);
  %
  % select the center voxel in each slice and grow the VOI
  %
  for z = cg(1)-ch:cg(1)+ch
    v(z,cg(2),cg(3)) = 1;
    for i=1:nIterations
        v(z,:,:) = growROI (squeeze(v(z,:,:)));
    end
  end
