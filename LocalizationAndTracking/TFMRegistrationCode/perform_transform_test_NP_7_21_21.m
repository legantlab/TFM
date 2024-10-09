%Edited 7-21-21


%% set up parameters
clear
clc
close all

path = 'X:\Max\RegistrationSampleData';
filename = '\20%F127D14DayOldP.tif';
filelabel = 'testing';
imgfile = strcat(path,filename);
stack_info = imfinfo(imgfile);
z_start = 1;            % layer with bottom of the channel
z_end = length(stack_info);              % layer with top of the channel

thresh_sensitivity = 0.6;     % threshold for creating binary image

min_area = 5E3;               %minimum area to be feed into object detection
max_area = 1E9;                %maximum area to be feed into object detection
min_solid = 0.4;              %minimum solidity for object detection
max_solid = 1;               %maximum solidity for object detection

x_size = stack_info(1).Width;
y_size = stack_info(1).Height;

%% Testing Binarization


img_temp = imread(imgfile, 3);
img_temp_bin = imbinarize(img_temp,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);
BW2 = bwpropfilt(img_temp_bin,'Area',[min_area max_area],6);
BW2 = bwpropfilt(BW2,'Solidity',[min_solid max_solid]); % inverts binarization so channel is selected
reg_data = [];

[reg_data(:,2), reg_data(:,1)] = find(BW2); % finds coordinates occupied by channel
reg_data(:,2) = y_size - reg_data(:,2); % flips y-coordinate for plotting
P = polyfit(reg_data(:,1), reg_data(:,2), 1); % linear regression coefficients of channel 
angle = atand(P(1)); % angle between channel and x-axis

img_rotated = imrotate(BW2, -angle+90); % orients channel vertically
img_rotated_actual = imrotate(img_temp_bin, -angle + 90); % not binary 
figure(1)
imshow(img_rotated)
% checking
figure(2)
imshow(BW2)
hold on
plot(1:x_size, y_size - polyval(P,1:x_size), 'r', 'LineWidth', 2)

%% Unnecessarily lengthy way of doing something that's probably very straight forward
reg_data_rot = [];
[reg_data_rot(:,2), reg_data_rot(:,1)] = find(img_rotated); 
% finds coordinates of channel after rotation

top_channel = max(reg_data_rot(:,2)); % bottom-most point of channel(larger index value)
bottom_channel = min(reg_data_rot(:,2)); % top-most point of channel(smaller index value)
channel_height = top_channel - bottom_channel;
left_channel = min(reg_data_rot(:,1)); % left-most point of channel
right_channel = max(reg_data_rot(:,1)); % right-most point of channel
channel_width = right_channel - left_channel;

cent_x = channel_width / 2 + left_channel; % center of channel in x

img_cropped = imcrop(img_rotated, [left_channel bottom_channel channel_width channel_height]);
% cropping rotated image to isolate the channel

sum_vals = sum(img_cropped, 2); % finds how much of each row is occupied by channel
not_channel = sum_vals < size(img_cropped, 2)*0.1; 
% rows where channel occupies low fraction of total width cropped image
box_dims = find(not_channel); % row number corresponding to rows with no channel
box_no_outs = rmoutliers(box_dims); % removes outliers -- not the center of channel

cent_y = mean(box_no_outs) + bottom_channel; % center of channel in y

x_adj = size(img_rotated, 1)/2 - cent_x;
y_adj = size(img_rotated,2)/2 - cent_y;

img_centered = imtranslate(img_rotated_actual, [x_adj, y_adj]);

figure(3)
imshow(img_centered)
hold on 
plot([0 size(img_centered, 1)], [size(img_centered, 2)/2 size(img_centered, 2)/2], 'b',...
    [size(img_centered,1)/2 size(img_centered,1)/2], [0 size(img_centered,2)], 'b', 'LineWidth', 2)
hold off

%% Putting Everything into a loop (7-13-2021 -- NP)
k = 1;
% first loop gathers rotation angles to align channels with vertical
img_temp_bin = {};
BW3 = {};

for i = 30:3:66
    %thresholding and binarizing images
    img_temp = imread(imgfile, i); % read in DIC images
    img_temp_bin{k} = imbinarize(img_temp,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);
    
    BW2 = bwpropfilt(img_temp_bin{k},'Area',[min_area max_area],6);
    BW2 = bwpropfilt(BW2,'Solidity',[min_solid max_solid]); % inverts binarization so channel is selected
    BW3{k} = imfill(BW2, 'holes');
    CC = bwconncomp(BW3{k}); % finds number of objects of binarized image
    if CC.NumObjects ~= 2
        continue
    end
    %BW2 = bwareafilt(BW2,2); % filters out all non-channel binarized objects
    % Note: sometimes the bubble/not channel stuff will be larger than the
    % channel space -- if necessary, maybe do a "pick the best looking one
    % then match the rest of the images to this nice looking image"
    
    
    % orienting channels with vertical
    reg_data = [];
    [reg_data(:,2), reg_data(:,1)] = find(BW3{k}); % finds coordinates occupied by channel
    reg_data(:,2) = y_size - reg_data(:,2); % flips y-coordinate for plotting
    P = polyfit(reg_data(:,1), reg_data(:,2), 1); % linear regression coefficients of channel
    angle = atand(P(1)); % angle between channel and x-axis
    stored_in_channel{k} = reg_data;
    stored_angs(k) = angle;
    k = k+1;
end

avg_ang = mean(rmoutliers(stored_angs)); % average rotation value

%%
% second loop rotates all layers by the average rotation angle, then
% isolates the channel to find translation amounts
m = 1;
img_rotated_actual = {};
for j = 30:3:66
    % rotating image
    img_rotated = imrotate(BW3{m}, 90-avg_ang); % orients channel vertically
    img_rotated_actual{m} = imrotate(img_temp_bin{m}, 90-avg_ang); % not binary
    stored_imgs{m} = img_rotated;
    % after rotation
    reg_data_rot = [];
    [reg_data_rot(:,2), reg_data_rot(:,1)] = find(img_rotated);
    % finds coordinates of channel after rotation
    rotated_in_channel{m} = reg_data_rot;
    
    % isolating channel
    top_channel = max(reg_data_rot(:,2)); % bottom-most point of channel(larger index value)
    bottom_channel = min(reg_data_rot(:,2)); % top-most point of channel(smaller index value)
    channel_height = top_channel - bottom_channel;
    left_channel = min(reg_data_rot(:,1)); % left-most point of channel
    right_channel = max(reg_data_rot(:,1)); % right-most point of channel
    channel_width = right_channel - left_channel;
    
    channel_dims{m} = [channel_width, channel_height]; % storing channel dimensions, WxL
    
    cent_x = channel_width / 2 + left_channel; % center of channel in x
    
    img_cropped = imcrop(img_rotated, [left_channel bottom_channel channel_width channel_height]);
    % cropping rotated image to isolate the channel
    stored_crops{m} = img_cropped;
    
    sum_vals = sum(img_cropped, 2); % finds how much of each row is occupied by channel
    not_channel = sum_vals < size(img_cropped, 2)*0.1;
    % rows where channel occupies low fraction of total width cropped image
    box_dims = find(not_channel); % row number corresponding to rows with no channel
    box_no_outs = rmoutliers(box_dims); % removes outliers -- not the center of channel
    
    cent_y = mean(box_no_outs) + bottom_channel; % center of channel in y
    
    x_adj(m) = size(img_rotated, 1)/2 - cent_x; % adjustment in x to center of image
    y_adj(m) = size(img_rotated,2)/2 - cent_y; % adjustment in y
    while m <= k - 2
        m = m+1;
    end
end

% finding an average rotation/translation
avg_x_adj = mean(rmoutliers(x_adj));
avg_y_adj = mean(rmoutliers(y_adj));
%%
k = 1;
% third loop translates all channel layers by an average amount
for j = 30:3:66
    
    % image translation
    img_centered = imtranslate(img_rotated_actual{k}, [x_adj(k), y_adj(k)]);
    % images are translated by an average amount
    stored_centered_imgs{k} = img_centered;
    
    % saving the image
    sizes(k) = size(img_centered,1);
    k = k+1;
end

%% Loop for Entire Stack

img_binarized = {};
img_initial = {};
img_rot = {};
img_trans = {};
l = 1;
for i = z_start:z_end
% testing this out
    img_initial{l} = imread(imgfile,i);
    %img_binarized{l} = imbinarize(img_initial{l},'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);
    img_rot{l} = imrotate(img_initial{l}, 90 - avg_ang);
    img_trans{l} = imtranslate(img_rot{l}, [avg_x_adj, avg_y_adj]);
    l = l+1;
end
    
%% Saving Image Trials

t = Tiff('testing6.tif','w');
tagStruct.ImageLength = max(sizes);
tagStruct.ImageWidth = max(sizes);
tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
tagStruct.BitsPerSample = 16;
tagStruct.SamplesPerPixel = 1;
tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagStruct.Software = 'MATLAB';

setTag(t,tagStruct)
write(t,squeeze(im2uint16(img_trans{1})));

for k = 2:length(img_trans)
    writeDirectory(t);
    setTag(t,tagStruct)
    write(t,squeeze(im2uint16(img_trans{k})));
end

close(t)

% img_edge_test = imread(imgfile,36);
% img_edge = imbinarize(img_edge_test,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);

%% Edge Detection Test

[~, threshold] = edge(img_edge, 'sobel');
fudgeFactor = 0.7;
BWs = edge(img_edge, 'sobel', threshold*fudgeFactor);
imshow(BWs)
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');

BWnobord = imclearborder(BWdfill, 4);

seD = strel('diamond', 1);
BWfinal = imerode(BWnobord, seD);
BWfinal = imerode(BWfinal, seD);


BWoutline = bwperim(BWfinal);
Segout = img_edge_test;
Segout(BWoutline) = 255;
imshowpair(labeloverlay(img_edge_test, BWfinal), Segout, 'montage')

% try BWareafill or something along those lines then filter for the largest
% object identified, should take out extraneous bubbles 

%% read images and create object stats
close all
k = 1;
for i = z_start:3:143*3
    img_temp = imread(imgfile,i);
    img_temp_bin = imbinarize(img_temp,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6);
    BW2 = bwpropfilt(img_temp_bin,'Area',[min_area max_area],6);
    BW2 = bwpropfilt(BW2,'Solidity',[min_solid max_solid]);
    stats = regionprops(BW2,'area','Centroid','Orientation','Solidity'); 
    ori(k) = arrayfun(@(x) x.Orientation, stats); % angle between x-axis and ellipse created with channel
    cent_cord{k} = arrayfun(@(x) x.Centroid, stats, 'UniformOutput', false);
    k = k + 1;
end    
% image is read as y, x, z


%% create stats on orientation;
[ori,TF] = rmoutliers(ori); 
cent_cord = cent_cord(~TF); % centroid coordinates of all non-outlier images

for i = 1:length(ori)
    axang = mean(ori)*pi/180*[0 0 1]; % average orientation angle for rotation about Z
    rot_mat = rotationVectorToMatrix(axang); 
    adj_x_val = cent_cord{i}(1)-x_size/2; % shifts the centroid in relation to center of image
    adj_y_val = cent_cord{i}(2)-y_size/2; % shifts y-coord of centroid
    rot_vals = rot_mat(1:2, 1:2) * [adj_x_val; adj_y_val]; % rotates adjusted coords in relation to Z
end
    
% disp_vec_y = [0 mean(cent_y_dist)] / rot_mat(1:2,1:2) + [x_size/2 y_size/2];
% disp_vec_x = [0 mean(cent_x_dist)] / rot_mat(1:2, 1:2) + [x_size/2 y_size/2];
% disp_vec = [disp_vec_x(1) disp_vec_y(2)];
% calculated average center of all channels in stack 

% cent_cord_rot = cellfun(@(x) [x(1)-x_size/2 x(2)-y_size/2]*rot_mat_xy,cent_cord,'UniformOutput',false);
% cent_y_dist = cellfun(@(x) x(2), cent_cord_rot);



%img_temp_rot = imrotate(img_temp,-mean(ori));

% for i = 143
img_temp_rot = imrotate(img_temp_bin, ori(end));
adj_diff = cent_cord{end} - disp_vec;
img_rot_adj = imtranslate(img_temp_rot, adj_diff./[x_size y_size]);
% end

subplot(1,2,1)
imshow(img_temp_bin)
subplot(1,2,2)
imshow(img_temp_rot)

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


