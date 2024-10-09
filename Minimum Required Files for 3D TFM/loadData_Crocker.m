function [data]=loadData_Crocker(filestem,imageName,dimensions)
%This function loads an image stack (saved as multipage TIFF format) into Matlab

%Input argument definitions
% filestem = a string indicating the directory where the image stack is located
% imageName = the name of the image
% dimensions = a 1 by 3 vector indicating the [xdimensions, y dimensions and z dimensions] of the image stack in pixels

%Output argument definitions
% data = a 3D matrix containing the pixel values of the image stack

frs=dimensions(3);
for k=1:frs
data(:,:,k)=imread(strcat(filestem,imageName,'.tif'),k);
end
data=permute(data, [2,1,3]);