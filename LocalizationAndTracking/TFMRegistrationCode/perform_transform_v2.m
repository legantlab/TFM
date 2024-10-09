%% read image stack
path = 'Z:\Regan\2020_02_20(Channel_Gel_40x_Nikon)\GPUdecon\';
filename = '40X_0.3um_24ms_Orient1_ch1_decon.tif';
imgfile = strcat(path,filename);
stack_info = imfinfo(imgfile);
z_length = size(stack_info,1);
img_stack = imread(imgfile,1);
for i = 2:z_length
    img_temp = imread(imgfile,i);
    img_stack = cat(3,img_stack,img_temp);
end    

% image is read as y, x, z


%% create a series of centroids of detected features along x direction (look at yz sections)
imgfreq = 100;            %the frequency of calculating vectors
y_size = size(img_stack,1);
x_size = size(img_stack,2);
imgsection = [];
edge_rm_start = 650;       % remove the first and last N frames to reduce edge effect
edge_rm_end = 50;         % remove last N frames to reduce edge effect
cent_list = [];
sec_count = 1;
y_prev = 800;      % initial guess of approximate y position
sign = -1;
nan_count = 0;     % keep track of the number of nans in between
for i = edge_rm:imgfreq:x_size - edge_rm
    imgsection = reshape(img_stack(:,i,:),y_size,[])';
    imgsection_adj = imadjust(imgsection,[],[],0.1);
    BW2 = bwpropfilt(~imbinarize(imgsection_adj),'convex',[1000 7000]);
    stats = regionprops(BW2,'centroid');
    cent_cord = [nan,nan];
    if ~isempty(stats)
        cent_cord = [nan,nan];
        for j = 1:length(stats)
            y = stats(j).Centroid(1);
            z = stats(j).Centroid(2);
            if z > z_length/2
                if sign*(y-y_prev) < (nan_count+1)*imgfreq && sign*(y-y_prev)>0
                    cent_cord = [y,z];
                    break
                end
            end
        end
    end
        
    if ~isnan(cent_cord(1))
        y_prev = y;
        nan_count = 0;
    else
        nan_count = nan_count+1;
    end

                

    cent_list(sec_count,:) = cent_cord;
    sec_count = sec_count+1;
    
end


%% calculate rotation matrix
pr = [0.1625,0.1625,0.3];
% orig_vec = [0,100,0].*pr;
% rot_vec = [yz_vec(2),x_disp,yz_vec(1)].*pr;
orig_vec = [x_disp,0,0].*pr;
rot_vec = [100,yz_vec(2),yz_vec(1)].*pr;
rot_mat = vrrotvec2mat(vrrotvec(rot_vec,orig_vec));

%% calculate translational motion   (currently assume only need to move along y direction)
% channel_width = 25;    % in um
% channel_height = 10;    % in um
% standard_img = zeros(z_length,y_size);
% standard_img(floor(z_length - channel_height/pr(3)):z_length,floor(y_size/2 - channel_width/(2*pr(1))):ceil(y_size/2 + channel_width/(2*pr(1)))) = 1;
middle_img = reshape(img_stack(:,round(y_size/2),:),y_size,[])';
BW2 = bwpropfilt(~imbinarize(middle_img),'convex',[2000 3000]);
stats = regionprops(BW2,'centroid');
cent_cord = stats(1).Centroid;
y_disp = y_size/2 - cent_cord(1);

%% perform transform
beadpos = x{1,1}{1,1};
beadpos_tranformed = beadpos + [-x_size/2,y_disp-y_size/2,0];
beadpos_tranformed = beadpos_tranformed*rot_mat;


