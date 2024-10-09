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
%                 binary images via image thresh  olding operations. This
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

%% Script Inputs: File Locations for Microscopy Images (TIF)
%Loading user inputs into the base workspace for ease of use with GUI
if isfile('userInputs.mat')
    load('userInputs.mat');
else
    firstTimeSetup()
end

userGUI = userInputGUI;

waitfor(userGUI)

save('userInputs.mat','isCZI','postStem','postTreatmentFile','postTreatmentChannel', ...
    'preStem','preTreatmentFile','preTreatmentChannel','analysisStem', ...
    'dataFileName', 'size_x', 'size_y', 'size_z','beadthres','beadmin','beadmax', ...
    'beadwin', 'beadmult', 'tptknnFD','tptknnFM','tptfmThres','tptoutThres','beadint', ...
    'tptoutThres');

%rebuild beadParam and tptParam to maintain consistency
beadParam{1}.thres = beadthres;
beadParam{1}.minSize = beadmin;          
beadParam{1}.maxSize = beadmax;
beadParam{1}.winSize = [beadwin, beadwin, beadwin];  % size of long pass filter during image preprocessing, set approximately the bead size
beadParam{1}.intensity = beadint; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number
beadParam{1}.thresh_multi = beadmult; %The intensity threshold set to be thresh_multi x the mode of raw image for particle detection

tptParam{1}.knnFD = tptknnFD;     %test particle linking, default is 16
tptParam{1}.knnFM = tptknnFM;
tptParam{1}.fmThres = tptfmThres;
tptParam{1}.outlrThres = tptoutThres;
tptParam{1}.nSpheres = 2;
tptParam{1}.maxIter = 5;


%% START
%full filenames

if ~isfile([postStem postTreatmentFile])
    error('post treatment file does not exist')
end

tiffs=dir([preStem,'\*.ome.tiff']);
%tiffs=dir([preStem,'\*.tif']);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'})); %Removes the folder and subfolder reference.

fileInfo = cell(1+length(preTreatmentFile),2);
fileInfo{1,1} = [postStem postTreatmentFile]; 
fileInfo{1,2} = postTreatmentChannel;
for i = 1:length(tiffs)
    fileInfo{i+1,1} = [tiffs(i).folder '\' tiffs(i).name];
    fileInfo{i+1,2} = preTreatmentChannel{1};
end

% Track Particles
[x, track, translation] = funTPT_V2(fileInfo, beadParam, tptParam, isCZI);
% Account for pixel to unit conversion here:
for i = 1:length(x)-1
    curX = x{i}{1};
    curX(:,1) = size_x.*curX(:,1);
    curX(:,2) = size_y.*curX(:,2);
    curX(:,3) = size_z.*curX(:,3); 
    x{i}{1} = curX;
end

%% Scale Bead Locations
%Rearrange x 
matches = cell(length(track)-2,1);
displacements = cell(length(track),1);

for i = 1:length(matches)
    curtracks = track{i}{1};
    %Optional drift correction if needed
    drift = mean(curtracks(:,2:4) - curtracks(:,6:8)); 
    %curtracks(:,2:4) = curtracks(:,2:4) - drift;
    curtracks(:,6:8) = curtracks(:,6:8) + drift;
    matches{i} = [curtracks(:,2:4), curtracks(:,6:8)];
    
    displacements{i} = curtracks(:,6:8) - curtracks(:,2:4);

end

beadLocations = matches{1}(:,[1 2 3]);
save([analysisStem '\' dataFileName '.mat'], 'matches', 'displacements', 'translation');

%% Plotting

displacementVectorScaling = 1;
%%{
figure(1)
time = 1;
scatter3(matches{time}(:,1),matches{time}(:,2),matches{time}(:,3))
hold on
scatter3(matches{time}(:,4),matches{time}(:,5),matches{time}(:,6))
hold off
title('postTreatment (blue) ancd preTreatment (orange) bead locations')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Posttreatment','Pretreatment')

figure(2)
d = 1;
quiver3(matches{1}(:,1),matches{1}(:,2),matches{1}(:,3),...
    d*displacements{1}(:,1),d*displacements{1}(:,2),d*displacements{1}(:,3),0)
title('Bead Displacements with Drift Correction')
xlabel('x')
ylabel('y')
zlabel('z')

%Plot a 2D projection of the pre and post beads with vectors linking their
%positions
figure(3)
hold on
time = 10;
d=1;
scatter(matches{time}(:,1), matches{time}(:,2), 'r.')
scatter(matches{time}(:,4), matches{time}(:,5), 'g.')
% scatter(x{1}{1}(track{1}{1}==0,1), x{1}{1}(track{1}{1}==0,2), 'r')
% scatter(x{time}{1}(setdiff(1:size(x{time}{1},1),track{1}{1}),1), ...
%     x{time}{1}(setdiff(1:size(x{time}{1},1),track{1}{1}),2), 'g')
% quiver(matches{time}(:,1), matches{time}(:,2),...
%     d*displacementsWithDrift{time}(:,1), d*displacementsWithDrift{time}(:,2),0, 'g')
quiver(matches{time}(:,1),matches{time}(:,2),...
    d*displacements{time}(:,1),d*displacements{time}(:,2),0,'b')
axis equal
hold off
