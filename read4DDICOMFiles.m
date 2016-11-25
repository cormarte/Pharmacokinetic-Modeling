function image = read4DDICOMFiles()

dir = 'D:\Documents de Corentin\Images\DICOM\Dynamiques\CEREBRIX\PET PETCT_CTplusFET_LM_Brain (Adult)\dynamic recon 3x10min Volume (Corrected) - 7';
[files, path, index] = uigetfile({'*.dcm','DICOM Files'; '*.*','All Files' },'DICOM Browser', 'MultiSelect', 'on', dir);

info = dicominfo(fullfile(path, files{1}))
rows = info.Rows;
columns = info.Columns;
slices = info.NumberOfSlices;
frames = info.NumberOfTimeSlices;

image = zeros(rows, columns, slices, frames, 'uint16');

for i = 1:length(files)
    
    info = dicominfo(fullfile(path, files{i}));
    index = info.ImageIndex;
    frame = idivide(index-1, slices, 'floor') + 1;
    slice = (index - (frame-1)*slices);
    image(:,:, slice, frame) = dicomread(info);
end

%figure, imshow(image(:,:,39,1));
%imcontrast();

return;