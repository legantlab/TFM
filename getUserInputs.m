function [] = getUserInputs()
%getUserInputs Queries user for inputs with GUI
%   Opens User input GUI to provide file inputs for TFM code and outputs
%   inputs to .mat file for the rest of the code to use later. Will save
%   .mat file to directory to use as default inputs for ease of use/only
%   slight delay over hard coding. 

% -- MH 2020

%%
% User inputs include:
% isCZI: Check if file input is of file format czi for different processing
% file location data: file stem, pre/post treatment files, channel data,
% data/output file name, analysis stem
% Voxel size: Size of each x,y,z voxel in image
% Bead parameters
% TPT parameters
% 

%%
%Search for .mat input file and preset values, else initialize to 0 is
%is fine. 

% if isfile('userInputs.mat')
%     
% else
%     
%     %Data variables
%     isCZI = 0;
%     postStem = 0;
%     postTreatmentFile = 0; 
%     postTreatmentChannel = NaN;
%     preStem = 0;
%     preTreatmentFile = 0;
%     preTreatmentFileName = 0; 
%     preTreatmentChannel = {NaN};
%     analysisStem = 0;
% 
%     dataFileName = 0;
%     
%     %Voxel Sizes
%     size_x = 1;
%     size_y = 1;
%     size_z = 1;
%     
%     %Bead parameters
%     beadthres = 0;
%     beadmin = 0;
%     beadmax = 0; 
%     beadwin = 0;
%     Beadint = 0;
%     beadmult = 0;
%     beadParam{1}.thres = 0.5;
%     beadParam{1}.minSize = 3;          
%     beadParam{1}.maxSize = 1000;
%     beadParam{1}.winSize = [5, 5, 5];  % size of long pass filter during image preprocessing, set approximately the bead size
%     beadParam{1}.intensity = 4500; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number
%     beadParam{1}.thresh_multi = 1.5; %The intensity threshold set to be thresh_multi x the mode of raw image for particle detection
%     
%     %TPT parameters
%     tptknnFD = 0;
%     tptknnFM = 0;
%     tptfmThres = 0;
%     tptoutThres = 0;
%     tptParam{1}.knnFD = 20;     %test particle linking, default is 16
%     tptParam{1}.knnFM = 6;
%     tptParam{1}.fmThres = 2;
%     tptParam{1}.outlrThres = 6;
% end
%Open selection GUI for data input

userGUI = userInputGUI;

%wait for the GUI continue button to be pressed
waitfor(userGUI)

%Output .mat file with input parameters
end