%% Run file for Topology-based Particle Tracking
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   fileInfo: string for the filename prefix to load the volumetric images
%             in the current directory.
%             Input options:
%             --- If the image is not within a cell ---
%             1) fileInfo{c}{1} = 'filename*.mat' or 'filename'
%               
%             --- If image is within a cell (Recommended approach) ---  
%             2) fileInfo{c}{1} = 'filename*.mat' or 'filename' 
%                fileInfo{c}{2} = Channel number containing images you want
%                                 to run TPT on. If the channel is not
%                                 provided, then channel = 1.
%               
%   beadParam: Parameters to detect and localize particles in the images              
%              Input options: 
%              1) beadParam{c}.thres = value between 0 & 1. (Default 0.5)  
%                                      default value = 0.5
%
%                 The threshold value used to convert input images into
%                 binary images via image thresholding operations. This
%                 operation is used to detect particles in the images.
%               
%              2) beadParam{c}.minSize = int value between 0 and Inf 
%                                        default value = 3
%
%                 The minimum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size smaller than this parameter are discarded as
%                 noise.
%                   
%              3) beadParam{c}.maxSize = int value between 1 and Inf     
%                                        default value = Inf
%
%                 The maximum size of particles in voxels as connected 
%                 components in the binary images. The connected components   
%                 with size larger than this parameter are discarded as   
%                 multiply-connected particles  
%                   
%              4) beadParam{c}.winSize = size in a three column array     
%                                        default value = [7, 7, 7]
%
%                 The image subset size used to localize particle using
%                 radial symmetry method. Select a size such that a single
%                 particle just fits in the image subset.
%                   
%   tptParam:  Parameters used by T-PT to track particles              
%              Input options:        
%                   
%              1) tptParam{c}.knnFM = int value (default = 5)     
%                   
%                 The total number of neighboring particles ('q') in the
%                 that are investigated in similarity of neighborhood test.  
%                   
%              2) tptParam{c}.fmThres = int value (default = 2)    
%                   
%                 The total number of neighboring particles ('p') out of 
%                 'q' that need to match between the reference and deformed  
%                 images to satisfy simalrity of neighborhood test.  
%                   
%                 Note: Increase p/q ratio for low spatial frequency   
%                 displacement field and decrease p/q ratio for high 
%                 spatial frequency displacement  
%                   
%              3) tptParam{c}.outlrThres = real value > 0 (default = 5)   
%                   
%                 The threshold value for the normalized residual in the  
%                 universal median test. Increase for high spatial  
%                 frequency displacement field and decrease for low spatial 
%                 frequency displacement field.   
%                   
%              4) tptParam{c}.knnFD = int value > 0 (default = 16)     
%                   
%                  The number of nearest neighboring particles used to 
%                  compute particle descriptor 
%                   
%              5) tptParam{c}.nSpheres = int value > 0 (default = 2)    
%                   
%                  The number of concentric shells used to compute particle 
%                  descriptor. 
%                   
%   NOTE- 'c' lists the order in which multi-attribute particles are                    
%         tracked. If only single attribute particles, c = 1. If           
%         multi-attribute particles, list all the properties for different           
%         c values (particles)          
%                   
%                   
%  OUTPUTS                 
%  ------------------------------------------------------------------------                 
%   x:      Particle positions found for each 'c' particle at each 'time'                   
%           frame
%           
%           x{time}{c} = particle positions at time frame in MxNxO format        
%                   
%   track:  Particles links between two consecutive image frames (t = to & 
%           t = t0 + 1) for each 'c' particle at each time frame (t = time)
%                   
%                   
%           track{time}{c} = an array of size [length(x{time}{c}, 1]. The         
%                            i_th index in the array stores the index of
%                            matched particle in x{time}{c}. If i_th index    
%                            value is 0, no particle link is found in the 
%                            next image frame
               
clear; close all;

%% Script Inputs: File Locations for Microscopy Images (CZI)
%{
isCZI = 1;
stem = 'Z:\WRL\2019_06_12(PA_hydrogels_with_200nm_580_605beads)\flipped\GPUdecon\';
postTreatmentFile = '10kPa_63xobj_postTreatment_sameZ.czi'; 
postTreatmentChannel = 0;
preTreatmentFile = {'10kPa_63xobj_preTreatment.czi'}; 
preTreatmentChannel = {0};

dataFileName = '06_21_POS14_CamA_ch0_CAM1_stack0000_405nm_0000000msec_0185328912msecAbs_000x_000y_000z_0000t_decon'; %particle tracking data will be stored under this name for future reference
%}
%% Script Inputs: File Locations for Microscopy Images (TIF)
%%{

%%Ask the User for this please dear lord with option to hardcode/saves
%%previous input for testing purposes



%%
isCZI = 0;
postStem = 'Z:\Regan\2020_04_28(TFM_troubleshooting)\GPUdecon\';
postTreatmentFile = 'Post_bead_03_06_20_ROI6_decon.tif'; 
postTreatmentChannel = NaN;
preStem = 'Z:\Regan\2020_04_28(TFM_troubleshooting)\GPUdecon\';
preTreatmentFile = {'Pre_bead_03_06_20_ROI6_decon'}; %This is cell
preTreatmentFileName = 'Pre_bead_03_06_20_ROI6_decon';  %This is string version of above
preTreatmentChannel = {NaN};
analysisStem = 'Z:\Regan\2020_04_28(TFM_troubleshooting)\Analysis\ROI6\';

dataFileName = '04_29_20_ROI6_Tracking_'; %particle tracking data will be stored under this name for future reference
%% 

%}
%% Other Inputs
%Voxel Sizes
size_x = 1;
size_y = 1;
size_z = 1;

% % Bead Parameter
% beadParam{1}.thres = 0.5;
% beadParam{1}.minSize = 5;
% beadParam{1}.maxSize = 1000;
% beadParam{1}.winSize = [7, 7, 7];
% beadParam{1}.intensity = 800; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number

% Bead Parameter
beadParam{1}.thres = 0.5;
beadParam{1}.minSize = 3;          
beadParam{1}.maxSize = 1000;
beadParam{1}.winSize = [5, 5, 5];  % size of long pass filter during image preprocessing, set approximately the bead size
%beadParam{1}.intensity = 4500; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number
beadParam{1}.thresh_multi = 1.5; %The intensity threshold set to be thresh_multi x the mode of raw image for particle detection

%%
%Checkpoint: Parameter selection menu, also for below


%%
% TPT Parameter
%tptParam{1}.knnFD = 20;     %test particle linking, default is 16
tptParam{1}.knnFM = 6;
tptParam{1}.fmThres = 2;
tptParam{1}.outlrThres = 6;



save('userInputs.mat','isCZI','postStem','postTreatmentFile','postTreatmentChannel', ...
    'preStem','preTreatmentFile','preTreatmentFileName','preTreatmentChannel','analysisStem', ...
    'dataFileName', 'size_x', 'size_y', 'size_z','beadParam', 'tptParam');

%% START
%full filenames

if ~isfile([postStem postTreatmentFile])
    error('post treatment file does not exist')
end

tiffs=dir([preStem,preTreatmentFileName,'*.tif']);

fileInfo = cell(1+length(preTreatmentFile),2);
fileInfo{1,1} = [postStem postTreatmentFile]; 
fileInfo{1,2} = postTreatmentChannel;
for i = 1:length(tiffs)
    fileInfo{i+1,1} = [tiffs(i).folder '\' tiffs(i).name];
    fileInfo{i+1,2} = preTreatmentChannel{1};
end

%%
%Checkpoint: Load in image file and check to make sure it is what you
% think it is. 


%%

% Track Particles
[x, track] = funTPT(fileInfo, beadParam, tptParam, isCZI);
% save([analysisStem dataFileName '_track.mat'], 'x', 'track', 'fileInfo');

%%
%Checkpoint: Bead tracking check; this will give us the bead locations and
%tracking information so we can display the movie again with that
%information overlayed to verify tracking is actually working. 

%% Scale Bead Locations
initLocations = x{1}{1};
matches = cell(length(track),1);
displacementsWithDrift = cell(length(track),1);
displacements = cell(length(track),1);
for i = 2:length(x)
    finLocations = x{i}{1};
    ind = find(track{i-1}{1} ~= 0);
    map = track{i-1}{1}(ind);
    matches_timet = [initLocations(ind,:) finLocations(map,:)];
%     
%     outliers = isoutlier(sqrt(sum(matches_timet(:,3).^2,2)));
%     sum(outliers)
%     matches_timet = matches_timet(~outliers,:);
    
    matches_timet(:,[1 4]) = size_x*matches_timet(:,[1 4]);
    matches_timet(:,[2 5]) = size_y*matches_timet(:,[2 5]);
    matches_timet(:,[3 6]) = size_z*matches_timet(:,[3 6]);
    
    driftDisplacements = matches_timet(:,4:6) - matches_timet(:,1:3);
    displacementsWithDrift{i-1} = driftDisplacements;
    
    %Translational drift and Quadratic surface regression method of drift correction
    %Translation
    matches_timet(:,4)=matches_timet(:,4)-mean(matches_timet(:,4)-matches_timet(:,1));
    matches_timet(:,5)=matches_timet(:,5)-mean(matches_timet(:,5)-matches_timet(:,2));
    matches_timet(:,6)=matches_timet(:,6)-mean(matches_timet(:,6)-matches_timet(:,3));
    
    displacementsNew = matches_timet(:,4:6) - matches_timet(:,1:3);
%     
%     %Quadratic
%     matchedCoordinates = [zeros(size(matches_timet,1),1) matches_timet(:,1:3) zeros(size(matches_timet,1),1) matches_timet(:,4:6)];
%     CENTROIDS = [zeros(size(matches_timet,1),1) matches_timet(:,1:3)];
%     [y,CENTROIDSR] = quadDef_match_to_ref(matchedCoordinates,CENTROIDS);
%     
%     %Rotate to compensate for sample alignment (right now this is manual
%     %and needs to be tuned, but eventually will be automated)
%     
%     matches_timet(:,4:6)=rotatePointCloudX(matches_timet(:,4:6),-.8);
%     CENTROIDSR(:,2:4)=rotatePointCloudX(CENTROIDSR(:,2:4),-.8);
%     
%     displacementsNew = matches_timet(:,4:6)-CENTROIDSR(:,2:4);
%     
    matches{i-1} = matches_timet;
    displacements{i-1} = displacementsNew;
    
end

beadLocations = matches{1}(:,[1 2 3]);
save([analysisStem dataFileName '.mat'], 'matches', 'displacements','initLocations','finLocations');

%%
%Checkpoint: This looks like where the actual bead locations are stored
%after all of the regression and fitting, it might be wise to check them
%here instead of earlier and add the ability to go back and fiddle with
%the parameter space; check with Yu
%%

disp('All the beads are contained in a box with the following ranges:')
disp(['x = [' num2str(min(beadLocations(:,1))) ', ' num2str(max(beadLocations(:,1))) ']'])
disp(['y = [' num2str(min(beadLocations(:,2))) ', ' num2str(max(beadLocations(:,2))) ']'])
disp(['z = [' num2str(min(beadLocations(:,3))) ', ' num2str(max(beadLocations(:,3))) ']'])
disp('particle tracking complete!');

%% Plotting
% Some plotting is done here; check if that is what are looking for
% earlier. We may be able to split this out to some functions and then call
% them to do the plotting previously. 

displacementVectorScaling = 1;
%%{
figure(1)
scatter3(matches{1}(:,1),matches{1}(:,2),matches{1}(:,3))
hold on
scatter3(matches{1}(:,4),matches{1}(:,5),matches{1}(:,6))
hold off
title('postTreatment (blue) and preTreatment (orange) bead locations')
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
d = displacementVectorScaling;
subplot(1,2,1)
quiver3(matches{1}(:,1),-1*matches{1}(:,2),-1*matches{1}(:,3),...
    d*displacementsWithDrift{1}(:,1),d*displacementsWithDrift{1}(:,2),-1*d*displacementsWithDrift{1}(:,3),0)
title('Raw Bead Displacement Vectors')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(1,2,2)
quiver3(matches{1}(:,1),matches{1}(:,2),matches{1}(:,3),...
    d*displacements{1}(:,1),d*displacements{1}(:,2),d*displacements{1}(:,3),0)
title('Bead Displacements with Drift Correction')
xlabel('x')
ylabel('y')
zlabel('z')

figure(3)
nearTopPointedUp = matches{1}(:,3) < 50 & displacements{1}(:,3) < 0;
quiver3(matches{1}(nearTopPointedUp,1),-1*matches{1}(nearTopPointedUp,2),-1*matches{1}(nearTopPointedUp,3),...
    d*displacements{1}(nearTopPointedUp,1),-1*d*displacements{1}(nearTopPointedUp,2),-1*d*displacements{1}(nearTopPointedUp,3),0)
hold on
nearTopPointedDown = matches{1}(:,3) < 50 & displacements{1}(:,3) > 0;
quiver3(matches{1}(nearTopPointedDown,1),-1*matches{1}(nearTopPointedDown,2),-1*matches{1}(nearTopPointedDown,3),...
    d*displacements{1}(nearTopPointedDown,1),-1*d*displacements{1}(nearTopPointedDown,2),-1*d*displacements{1}(nearTopPointedDown,3),0)
hold off
title('Bead Displacements with Drift Correction (near top of gel)')
xlabel('x')
ylabel('y')
zlabel('z')

%Plot a 2D projection of the pre and post beads with vectors linking their
%positions
figure(4)
hold on
d=displacementVectorScaling;
scatter(matches{1}(:,1), matches{1}(:,2), 'r.')
scatter(matches{1}(:,4), matches{1}(:,5), 'g.')
scatter(x{1}{1}(track{1}{1}==0,1), x{1}{1}(track{1}{1}==0,2), 'r')
scatter(x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),1), ...
    x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),2), 'g')
quiver(matches{1}(:,1), matches{1}(:,2),...
    d*displacementsWithDrift{1}(:,1), d*displacementsWithDrift{1}(:,2),0, 'g')
quiver(matches{1}(:,1),matches{1}(:,2),...
    d*displacements{1}(:,1),d*displacements{1}(:,2),0,'k')
hold off

Ipost = loadtiff(fileInfo{1});
Ipre = loadtiff(fileInfo{2});
figure(5)
Ifuse = imfuse(sum(Ipost,3)>8500, sum(Ipre,3)>8500,'blend','Scaling','independent');
imshow(Ifuse)
caxis([5 20])
hold on
scatter(matches{1}(:,1), matches{1}(:,2), 200, 'r.')
scatter(matches{1}(:,4), matches{1}(:,5), 200, 'g.')
scatter(x{1}{1}(track{1}{1}==0,1), x{1}{1}(track{1}{1}==0,2), 'r')
scatter(x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),1), ...
    x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),2), 'g')
quiver(matches{1}(:,1), matches{1}(:,2),...
    d*displacementsWithDrift{1}(:,1), d*displacementsWithDrift{1}(:,2),0, 'g')
quiver(matches{1}(:,1),matches{1}(:,2),...
    d*displacements{1}(:,1),d*displacements{1}(:,2),0,'k')
hold off
legend('matched PostSDS bead','matched PreSDS bead','unmatched PostSDS bead','unmatched PreSDS bead','location','eastoutside')

