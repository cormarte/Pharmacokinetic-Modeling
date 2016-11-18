function show_image (img, range, voi)
  if nargin < 2 || isempty (range)
    range = [ min(img(:)), max(img(:)) ];
  end
  range = double (range);
  if nargin == 3
    if (size (img, 1) ~= size (voi, 1)) || (size (img, 2) ~= size (voi, 2))
      warning ('image (%dx%d) and voi (%dx%d) have incompatible dimensions', size(img,1), size(img,2), size(voi,1), size(voi,2));
      voi = [];
    end
  end
  if nargin < 3 || ~nnz (voi)
    voi = [];
    alpha = 1.0;
  else
    opacity = 0.4;
    voi(voi>0) = 1;
    alpha = 1 - opacity * double (voi);
  end
  %
  % clip image to given range
  %
  s = size (img);
  img(find(img<range(1))) = range(1);
  img(find(img>range(2))) = range(2);
  img = single (reshape (img, s(1), s(2)));
  %
  % scale to colormap indices
  %
  if range(2) == range(1)
    img = int8 (1.0 + (length (colormap) - 1.0) * (img - range(1)));
  else
    img = int8 (1.0 + (length (colormap) - 1.0) * (img - range(1)) ./ (range(2) - range(1)));      
  end     
  %
  %
  %
  imgRGB = zeros (size(img,1), size(img,2), 3);
  %voiRGB = zeros (size(voi,1), size(voi,2), 3);
  imgMap = colormap;
  %voiMap = colormap (gray);
  for c = 1:3
    imgMapRows = img;
    imgMapRows(find(imgMapRows==0)) = 1;
    %voiMapRows = voi;
    imgTmp = imgMap(imgMapRows, c);
    %voiTmp = voiMap(voiMapRows, c);
    imgRGB(:,:,c) = reshape (imgTmp, size (img));
    %imgVoi(:,:,c) = reshape (voiTmp, size (voi));
  end
  image (imgRGB);
  jointRGB = zeros (size (imgRGB));
  for c = 1:3
    jointRGB(:,:,c) = alpha .* imgRGB(:,:,c);
  end
  image (jointRGB);
