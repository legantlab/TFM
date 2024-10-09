%% set up parameters
clear
clc
close all

path = 'X:\Max\RegistrationSampleData';
filename = '\20%F127D14DayOldP.tif';
imgfile = strcat(path,filename);
stack_info = imfinfo(imgfile);
z_start = 1;            % layer with bottom of the channel
z_end = length(stack_info);              % layer with top of the channel

thresh_sensitivity = 0.9;     % threshold for creating binary image

min_area = 5E3;               %minimum area to be feed into object detection
max_area = 1E9;                %maximum area to be feed into object detection
min_solid = 0.45;              %minimum solidity for object detection
max_solid = 1;               %maximum solidity for object detection



x_size = stack_info(1).Width;
y_size = stack_info(1).Height;


%% read images and create object stats
close all
stats = [];
for i = z_start:3:143*3
    img_temp = imread(imgfile,i);
    img_temp_bin = imbinarize(img_temp,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);
    BW2 = bwpropfilt(img_temp_bin,'Area',[min_area max_area],6);
    BW2 = bwpropfilt(BW2,'Solidity',[min_solid max_solid]);
    stats = [stats ; regionprops(BW2,'area','Centroid','Orientation','Solidity')];    
end    

% image is read as y, x, z


%% create stats on orientation;
ori = arrayfun(@(x) x.Orientation, stats); % angle between x-axis and ellipse created with channel
cent_cord = arrayfun(@(x) x.Centroid, stats, 'UniformOutput', false);
[ori,TF] = rmoutliers(ori); 
cent_cord = cent_cord(~TF); % centroid coordinates of all non-outlier images

axang = mean(ori)*pi/180*[0 0 1]; % average orientation angle for rotation about Z
rot_mat = rotationVectorToMatrix(axang); 

cent_y_dist = zeros(size(ori));

for i = 1:length(ori)
    cent_y_dist(i) = (cent_cord{i}(2)-y_size/2)*cosd(ori(i)) + (cent_cord{i}(1)-x_size/2)*sind(ori(i));
    % distance between y-coordinate of the centroid and y-axis
    % post-rotation
end
    
disp_vec = [0 mean(cent_y_dist)] / rot_mat(1:2,1:2) + [x_size/2 y_size/2];
% calculated average center of all channels in stack 

% cent_cord_rot = cellfun(@(x) [x(1)-x_size/2 x(2)-y_size/2]*rot_mat_xy,cent_cord,'UniformOutput',false);
% cent_y_dist = cellfun(@(x) x(2), cent_cord_rot);

img_temp_rot = imrotate(img_temp,-mean(ori));

% %% calculate rotation matrix
% pr = [0.1625,0.1625,0.3];
% % orig_vec = [0,100,0].*pr;
% % rot_vec = [yz_vec(2),x_disp,yz_vec(1)].*pr;
% orig_vec = [x_disp,0,0].*pr;
% rot_vec = [100,yz_vec(2),yz_vec(1)].*pr;
% rot_mat = vrrotvec2mat(vrrotvec(rot_vec,orig_vec));
% 
% %% calculate translational motion   (currently assume only need to move along y direction)
% % channel_width = 25;    % in um
% % channel_height = 10;    % in um
% % standard_img = zeros(z_length,y_size);
% % standard_img(floor(z_length - channel_height/pr(3)):z_length,floor(y_size/2 - channel_width/(2*pr(1))):ceil(y_size/2 + channel_width/(2*pr(1)))) = 1;
% middle_img = reshape(img_stack(:,round(y_size/2),:),y_size,[])';
% BW2 = bwpropfilt(~imbinarize(middle_img),'convex',[2000 3000]);
% stats = regionprops(BW2,'centroid');
% cent_cord = stats(1).Centroid;
% y_disp = y_size/2 - cent_cord(1);
% 
% %% perform transform
% beadpos = x{1,1}{1,1};
% beadpos_tranformed = beadpos - [disp_vec 0];
% beadpos_tranformed = beadpos_tranformed*rot_mat;


