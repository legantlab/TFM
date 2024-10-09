path = 'Z:\Yu\2020_03_18_TFM_channel_488_lifeactRFP\bead\GPUdecon\';
filename_1 = '60x_488_lifeactRFP_03z_ROI6_decon.tif';

beadpos = x{1,1}{1,1};
%beadpos = matches{1,1}(:,1:3);
%beadpos = matches{1,1}(:,4:6);
imgfile_1 = strcat(path,filename_1);
%imgfile_1 = fileInfo{1,1};
% stack_info = imfinfo(imgfile_1);
% z_length = size(stack_info,1);
z_plane = 50;
img_temp = imread(imgfile_1,z_plane);
z_range = 3;
select_array = (beadpos(:,3)>z_plane-z_range & beadpos(:,3)<z_plane+z_range);
beadpos_plane = beadpos(select_array,:);
%vec_plane = driftDisplacements(select_array,1:2);
figure();
imshow(img_temp,[7000 10000]);
%imshow(img_temp);
hold on
scatter(beadpos_plane(:,1),beadpos_plane(:,2),3,'o','r');

 
% figure();
% imshow(img_temp,[8000 20000]);
% hold on
% quiver(beadpos_plane(:,1),beadpos_plane(:,2),vec_plane(:,1),vec_plane(:,2));
% 
% filename_2 = '40x_20kPa_flat_488_560_100nm_xshift10um_green_decon_cropped.tif';
% imgfile_2 = strcat(path,filename_2);
% img_temp = imread(imgfile_2,z_plane);
% figure();
% imshow(img_temp,[8000 20000]);
% hold on
% quiver(beadpos_plane(:,1),beadpos_plane(:,2),vec_plane(:,1),vec_plane(:,2),'AutoScale','off');







% outputFileName = 'img_overlay_stack.tif'
% for K=1:length(Illu_dist_avg(1, 1, :))
% imwrite(Illu_dist_avg(:, :, K), outputFileName, 'WriteMode', 'append');
% end