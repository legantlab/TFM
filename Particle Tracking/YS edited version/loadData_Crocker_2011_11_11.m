function [data]=loadData_Crocker_2011_11_11(filestem,imageName,dimensions)
%This function loads an image stack (saved as multipage TIFF format) into Matlab

%Input argument definitions
% filestem = a string indicating the directory where the image stack is located
% imageName = the name of the image
% dimensions = a 1 by 3 vector indicating the [xdimensions, y dimensions and z dimensions] of the image stack in pixels

%Output argument definitions
% data = a 3D matrix containing the pixel values of the image stack

frs=dimensions(6)-dimensions(5);
for k=dimensions(5):dimensions(6)
data(:,:,k-dimensions(5)+1)=imread(strcat(filestem,imageName,'.TIF'),k);
end
data=permute(data, [2,1,3]);