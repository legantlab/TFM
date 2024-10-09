path = 'Z:\Regan\2020_01_31(TwoColor_Gels_Nikon)\Yu_test\channel_twocolor\GPUdecon\MIPs\';
filename_1 = '40x_20kPa_25umChannel_560_100nm_MIP_z.tif';

beadpos = matches{1,1}(:,1:3);
%beadpos = matches{1,1}(:,4:6);
imgfile_1 = strcat(path,filename_1);
%imgfile_1 = fileInfo{1,1};
% stack_info = imfinfo(imgfile_1);
% z_length = size(stack_info,1);
img_temp = imread(imgfile_1);
beadpos_plane = beadpos(:,1:2);
vec_plane = driftDisplacements(:,1:2);
figure();
imshow(img_temp);
%imshow(img_temp);
hold on
scatter(beadpos_plane(:,1),beadpos_plane(:,2),3,'o','r');
