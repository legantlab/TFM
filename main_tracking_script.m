%% Run file for Feature Vector-based Particle Tracking
% Generates .mat file containing linked particle tracks, particle XYZ
% positions,  and drift correction translations used for computation of the
% displacement field of 3D bead localizations. GUIs are use to
% simplify parameter input. Plotting functions are provided at the end to
% assess quality of displacement field calculation. 
%
% VARIABLES
% -------------------------------------------------------------------------
%   fileInfo: string for the filename prefix to load the volumetric images
%             in the current directory. 
%               
%   beadParam: Parameters to detect and localize particles in the images.
%   These are adjusted in a subsequent GUI which shows the resulting
%   detected and localized particles as a result of parameter selection. 
%
%              Input options: 
%              1) beadParam.thres = value between 0 & 1.
%                   default 0.5
%
%                 The threshold value used to convert input images into
%                 binary images via image thresh  olding operations. This
%                 operation is used to detect particles in the images.
%               
%              2) beadParam.minSize = int value between 0 and Inf. 
%                   default value = 3
%
%                 The minimum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size smaller than this parameter are discarded as
%                 noise.
%                   
%              3) beadParam.maxSize = int value between 1 and Inf.
%                   default value = Inf
%
%                 The maximum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size larger than this parameter are discarded as
%                 multiply-connected particles
%                   
%              4) beadParam.winSize = size in a three column array.
%                   default value = [7, 7, 7]
%
%                 The image subset size used to localize particle using
%                 radial symmetry method. Select a size such that a single
%                 particle just fits in the image subset.
%                   
%  OUTPUTS
%  ------------------------------------------------------------------------
%   matches:      Particle positions found for each particle at each 'time'
%   frame
%           
%           x{time} = particle positions at time frame in MxN format
%                   
%   displacements:  Particles links between two consecutive image frames (t
%   = to & t = t0 + 1) for each particle at each time frame (t = time)
%                                      
%           track{time} = an array of size [length(x{time}), 1]. The
%                            i_th index in the array stores the index of
%                            matched particle in x{time}. If i_th index
%                            value is empty, no particle link is found in the
%                            next image frame
%
%   translation:  Cell array containing the X, Y, and Z translations
%   applied to each frame to center the dataset. 
%                   
%           translation{time} = an array of size [length(x{time}), 1]. The
%           i_th index in the array stores the index of translation in x{time}. 
%
%   Author: Max Hockenberry
%   Last Update: 10/09/2024
clear; close all;

%% Script Inputs: File Locations for Microscopy Images (TIF)
%Loading user inputs into the base workspace for ease of use with GUI
if isfile('userInputs.mat')
    load('userInputs.mat');
else
    firstTimeSetup() %Script that generates default userInputs.mat if not present. 
end

userGUI = userInputGUI;

waitfor(userGUI)

save('userInputs.mat','postStem','postTreatmentFile','postTreatmentChannel', ...
    'preStem','preTreatmentFile','preTreatmentChannel','analysisStem', ...
    'dataFileName', 'size_x', 'size_y', 'size_z','beadthres','beadmin','beadmax', ...
    'beadwin', 'beadmult','beadint');

%assemble a structure that contains the bead parameters for easy function
%input. 
beadParam{1}.thres = beadthres; %The threshold to apply for bead selection scaled to normalized data [minIntensity (0),maxIntensity (1)].
beadParam{1}.minSize = beadmin; %The minimum size of connected voxels that will be considered one bead.
beadParam{1}.maxSize = beadmax; %The maximum size of connected voxels that will be considered one bead.
beadParam{1}.winSize = [beadwin, beadwin, beadwin];  % size of long pass filter during image preprocessing, set approximately the bead size.
beadParam{1}.intensity = beadint; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number
beadParam{1}.thresh_multi = beadmult; %The intensity threshold set to be thresh_multi x the mode of raw image for particle detection

%% START
%Retrieve all .tiff images in the treatment folder and assemble structure
%containing file locations. 

if ~isfile([postStem postTreatmentFile])
    error('post treatment file does not exist')
end

tiffs=dir([preStem,'\*.tiff']);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'})); %Removes the folder and subfolder reference.

fileInfo = cell(1+length(preTreatmentFile),2);
fileInfo{1,1} = [postStem postTreatmentFile]; 
fileInfo{1,2} = postTreatmentChannel;
for i = 1:length(tiffs)
    fileInfo{i+1,1} = [tiffs(i).folder '\' tiffs(i).name];
    fileInfo{i+1,2} = preTreatmentChannel{1};
end

% Track Particles
[x, track, translation] = trackBeads(fileInfo, beadParam);

% Convert pixel sizes into actual units using xyz voxel sizes defined
% previously. 
for i = 1:length(x)-1
    curX = x{i}{1};
    curX(:,1) = size_x.*curX(:,1);
    curX(:,2) = size_y.*curX(:,2);
    curX(:,3) = size_z.*curX(:,3); 
    x{i}{1} = curX;
end

%% Scale Bead Locations and Save Results
%Rearrange x (bead locations) into matches and displacements arrays. 
matches = cell(length(track)-2,1);
displacements = cell(length(track),1);

for i = 1:length(matches)
    curtracks = track{i}{1};
    %Optional drift correction if needed
    drift = mean(curtracks(:,2:4) - curtracks(:,6:8)); 
    curtracks(:,6:8) = curtracks(:,6:8) + drift;
    matches{i} = [curtracks(:,2:4), curtracks(:,6:8)];
    
    displacements{i} = curtracks(:,6:8) - curtracks(:,2:4);
end

%Save .mat file in analysis folder for further computation of traction
%forces. 
save([analysisStem '\' dataFileName '.mat'], 'matches', 'displacements', 'translation');

%% Plotting - Quality Control and Visualization of Displacement Fields

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
time = 1;
d=1;
scatter(matches{time}(:,1), matches{time}(:,2), 'r.')
scatter(matches{time}(:,4), matches{time}(:,5), 'g.')
quiver(matches{time}(:,1),matches{time}(:,2),...
    d*displacements{time}(:,1),d*displacements{time}(:,2),0,'b')
axis equal
hold off
