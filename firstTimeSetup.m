function firstTimeSetup()
%firstTimeSetup Creates .mat file with default tracking variables if not present
%   Creates .mat file containing default tracking variables if not present in
%   workplace folder. 
%
%   Inputs:  
%       None
%   Outputs:
%       Saves a .mat file into the workplace directory containing default
%       particle tracking values defined below. 
%
%   Author: Max Hockenberry
%   Last Update: 10/09/2024

%   Image Input variables
    postStem = 'test';
    postTreatmentFile = 'test'; 
    preStem = 'test';
    preTreatmentFile = {'test'};
    preTreatmentFileName = 'test'; 

    analysisStem = 'test';

    dataFileName = 'test';
    
    %Voxel Sizes in microns
    size_x = 0.2;
    size_y = 0.2;
    size_z = 0.8;
    
    %Bead tracking parameters
    beadthres = 0.5; %The threshold to apply for bead selection scaled to normalized data [minIntensity (0),maxIntensity (1)].
    beadmin = 3; %The minimum size of connected voxels that will be considered one bead.
    beadmax = 1000; %The maximum size of connected voxels that will be considered one bead.
    beadwin = 5;  % size of long pass filter during image preprocessing, set approximately the bead size.
    beadint = 4500; %The intensity threshold for a bead - If beads are dimmer (or if using a smaller bit image) decrease this number
    beadmult = 1.5; %The intensity threshold set to be thresh_multi x the mode of raw image for particle detection
    
    save('userInputs.mat','postStem','postTreatmentFile', ...
    'preStem','preTreatmentFile','preTreatmentFileName','analysisStem', ...
    'dataFileName', 'size_x', 'size_y', 'size_z','beadthres','beadmin','beadmax', ...
    'beadwin', 'beadmult','beadint');
end

