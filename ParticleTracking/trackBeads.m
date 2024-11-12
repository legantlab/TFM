function [x, track, translation] = trackBeads(fileNames, beadParameter)
% [x, track, translation] = trackBeads(fileInfo, beadParameter)
% is the main function that performs the single particle tracking on time
% increment of volumetric images.
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   fileInfo: string for the filename prefix to load the volumetric images
%             in the current directory as a cell array. 
%
%   beadParam: Parameters to detect and localize particles in the images
%              Input options: 
%              1) beadParam.thres = value between 0 & 1.
%                                        default value = 0.5
%
%                 The threshold value used to convert input images into
%                 binary images via image thresholding operations. This
%                 operation is used to detect particles in the images.
%
%              2) beadParam.minSize = int value between 0 and Inf
%                                        default value = 3
%
%                 The minimum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size smaller than this parameter are discarded as
%                 noise.
%
%              3) beadParam.maxSize = int value between 1 and Inf
%                                        default value = Inf
%
%                 The maximum size of particles in voxels as connected
%                 components in the binary images. The connected components
%                 with size larger than this parameter are discarded as
%                 multiply-connected particles
%
%              4) beadParam.winSize = size in a three column array
%                                        default value = [5, 5, 5]
%
%                 The image subset size used to localize particle using
%                 radial symmetry method. Select a size such that a single
%                 particle just fits in the image subset.
%
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   x:      Particle positions found for each particle at each 'time'
%           frame
%
%           x{time} = particle positions at time frame in MxNxO format
%
%   track:  Particles links between two consecutive image frames (t = to &
%           t = t0 + 1) for each particle at each time frame (t = time)
%
%
%           track{time} = an array of size [length(x{time}, 1]. The
%                            i_th index in the array stores the index of
%                            matched particle in x{time}. If i_th index
%                            value is 0, no particle link is found in the
%                            next image frame
%
%   track:  Particles links between two consecutive image frames (t = to &
%           t = t0 + 1) for each particle at each time frame (t = time)
%
%
%           track{time} = an array of size [length(x{time}, 1]. The
%                            i_th index in the array stores the index of
%                            matched particle in x{time}. If i_th index
%                            value is 0, no particle link is found in the
%                            next image frame
%   Author: Max Hockenberry
%   Last Update: 10/23/2024

%Parse inputs
bandpass_size = beadParameter{1}.winSize;
threshold = beadParameter{1}.thres;
threshmulti = beadParameter{1, 1}.thresh_multi;

%Define tracking variables - Variables to control the linking parameters
%for feature vector matching. These are settings that work well for the
%default data provided, but will need to be adjusted depending on input
%data. Provided near the linking functions are commented plotting figures
%that can help assist with assessing the results of changing the variables.
%The variables are in the following order: 
%   numNeighbors: The number of nearest neighbors to scan as a match
%
%   numVectorsRef: The number of feature vectors to use for the reference
%   dataset
%
%   numVectorsDef: The number of feature vectors to use for the deformed
%   dataset
%
%   tol: a cutoff value indicating how much more closely the feature
%   vectors of the ideal match need to be when compaired to all other
%   potential matches between the reference and deformed configurations in
%   order to consider it a correct match.

firstPassVars = [25, 3, 5, 2];
initialMatchingVars = [50, 3, 5, 1.5];

courseSmoothVars = [50,1.75]; %Smoothing filter uses numNeighbors and tol only. 
fineSmoothVars = [50,1.5];

trackingVars = [20,3,5,1.5];
trackingVars2 = [25,3,5,1.4];


%% Particle Tracking
% Handle first time point which is the reference state seperately; Need to
% apply drift correction and set so that later iterations can be run in
% parrallel. 
x = cell(length(fileNames),1);
track = cell(length(fileNames),1);
translation = cell(length(fileNames),1);

t = 2;
j = 1;


I = loadtiff(fileNames{t-1,1});
%Do DC/background correction
I = I - median(I,'all');
I(I<0) = 0;
        
% Parameters

%lnoise is the size of the gaussian filter, usually set to 0.5
%in all dimensions so we default to there. User input in the
%GUI can change this, so it will save and altar that after the
%GUI closes. 
lnoise = [.5 .5 .5];

%Inputb isn't used for anything in the bpass3dMB, it is
%supposed to be toggles for noclip/nopad but there is no use of
%it in the code at the moment. Just leave as [0,0]. 
inputb = [0,0];
bPassParams = {lnoise, bandpass_size, inputb};

[IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3}); %Bandpass filter - Can tune parameters for optimal fitting

bpGUI = bPassGUI(I,IBP,bPassParams,threshold,threshmulti); % GUI that allows fine tuning of bPass parameters. 
%waits for the user to press continue, then extracts value and closes app
while(bpGUI.closenum ~= 1)
    pause(1)
end

bPassParams = bpGUI.finalbpassParams; 

bandpass_size = bPassParams{2};
threshold = bpGUI.finalthresh;
threshmulti = bpGUI.finalmultithresh;

beadParameter{1}.winSize = bandpass_size;
beadParameter{1}.thres = threshold;
beadParameter{1, 1}.thresh_multi = threshmulti;
bpGUI.delete
         
% Perform particle localization. 
x{1}{j} = locateParticles(I,bandpass_size,bPassParams,threshold);

x{1}{j}(any(isnan(x{1}{j}),2),:) = [];
            
%% Correct for swelling induced by SDS detergent

disp('Correcting for SDS induced Swelling')

I = loadtiff(fileNames{t,1});
I = I - median(I,'all');
I(I<0) = 0;

x{t}{j} = locateParticles(I,bandpass_size,bPassParams,threshold);
x{t}{j}(any(isnan(x{t}{j}),2),:) = [];
%Center at 0
translation{t} = mean(x{t}{j},1);
x{t}{j} =   x{t}{j} - translation{t};

translation{1} = translation{t};
x{1}{1} = x{1}{1} - translation{t};

swellingGUI = swellingCorrectionGUI(x{1}{1},x{2}{1}); % Runs a GUI to control swelling correction. 

while(swellingGUI.closenum ~= 1)
    pause(1)
end

x{1}{1} = swellingGUI.finalBeads;
x{1}{1}(any(isnan(x{1}{1}),2),:) = [];

swellingGUI.delete

% Track beads using feature vectors. 
refBeads = [(1:length(x{1}{1}))',x{1}{1}];
strainBeads = [(1:length(x{t}{1}))',x{t}{1}];

firstPass = displacementCalc_recursive_Max_2024_02_19(refBeads,strainBeads,firstPassVars(1),...
    firstPassVars(2),firstPassVars(3),firstPassVars(4)); 

distNewMatches = firstPass(:,6:8) - firstPass(:,2:4);

% figure
% quiver(firstPass(:,2), firstPass(:,3), distNewMatches(:,1), distNewMatches(:,2),1)
% axis equal

refBeadsCor = refBeads;
strainBeadsCor = strainBeads;

dcMatches = displacementCalc_recursive_Max_2024_02_19(refBeadsCor,strainBeadsCor,initialMatchingVars(1)...
    ,initialMatchingVars(2),initialMatchingVars(3),initialMatchingVars(4));
testDCmatches = dcMatches(:,6:8) - dcMatches(:,2:4);
% figure 
% quiver(dcMatches(:,2), dcMatches(:,3), testDCmatches(:,1), testDCmatches(:,2),1)
matchedCoordinates_filt=smoothFilt_2024_02_20(dcMatches,courseSmoothVars(1), courseSmoothVars(2)); 
filt_Disps = matchedCoordinates_filt(:,6:8) - matchedCoordinates_filt(:,2:4);
disp(['Filtered Beads = ',num2str(length(dcMatches) - length(matchedCoordinates_filt))])
matchedCoordinates_filt2=smoothFilt_2024_02_20(matchedCoordinates_filt,fineSmoothVars(1),fineSmoothVars(2)); 

filt_Disps2 = matchedCoordinates_filt2(:,6:8) - matchedCoordinates_filt2(:,2:4);
disp(['Filtered Beads = ',num2str(length(matchedCoordinates_filt) - length(matchedCoordinates_filt2))])

% figure
% quiver(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3), filt_Disps2(:,1), filt_Disps2(:,2),1)
% axis equal

track{t-1}{j} = matchedCoordinates_filt2;


% Start Tracking - parallized for speed
parfor t = 3:length(fileNames) % Loop through all time points
    
    tStart = tic;
         
        % Detect and Localize Particles -----------------------------------
        tPP = tic;
     
        if t~=length(fileNames) % after first frame except last frame
            disp(['Current time point: t = ' num2str(t)])
            disp(['  Particle channel: ' num2str(fileNames{t,2})])
            disp(['  Current Filename: ' fileNames{t}])
            
       
            I = loadtiff(fileNames{t,1});
            I = I - median(I,'all');
            I(I<0) = 0;
            
            %localize particles
            x{t}{j} = locateParticles(I,bandpass_size,bPassParams,threshold);
            x{t}{j}(any(isnan(x{t}{j}),2),:) = [];
            
            translation{t} = mean(x{t}{j},1);
            x{t}{j} =   x{t}{j} - mean(x{t}{j},1);
            
            %Use Pointcloud registration to align data set to reference. 
            PC1 = pointCloud(refBeads(:,2:4));
            PC2 = pointCloud(x{t}{j});
            tform = pcregistericp(PC1, PC2);
            PC2reg = pctransform(PC2,tform);
            x{t}{j} = PC2reg.Location;      
            
            disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
            
            %Run tracking scripts

            %refBeads = [(1:length(x{1}{1}))',x{1}{1}];
            strainBeads = [(1:length(x{t}{1}))',x{t}{1}];

            firstPass = displacementCalc_recursive_Max_2024_02_19(refBeads,strainBeads,trackingVars(1),...
                trackingVars(2),trackingVars(3),trackingVars(4));

            distNewMatches = firstPass(:,6:8) - firstPass(:,2:4);

            %Implement a scheme where we compute the drift and apply that to the raw
            %displacements then retrack...
            xDrift = median(distNewMatches(:,1));
            yDrift = median(distNewMatches(:,2));
            zDrift = median(distNewMatches(:,3));
            refBeadsCor = refBeads;

            strainBeadsCor = strainBeads;
            strainBeadsCor(:,2) = strainBeadsCor(:,2) - xDrift;
            strainBeadsCor(:,3) = strainBeadsCor(:,3) - yDrift;
            strainBeadsCor(:,4) = strainBeadsCor(:,4) - zDrift;
            
            dcMatches = displacementCalc_recursive_Max_2024_02_19(refBeadsCor,strainBeadsCor,trackingVars2(1),...
                trackingVars2(2),trackingVars2(3),trackingVars2(4));
            
            % Run smoothing function twice for better smoothing (course and
            % fine grain)
            
            matchedCoordinates_filt=smoothFilt_2024_02_20(dcMatches,courseSmoothVars(1), courseSmoothVars(2)); 
            filt_Disps = matchedCoordinates_filt(:,6:8) - matchedCoordinates_filt(:,2:4);
            %disp(['Filtered Beads = ',num2str(length(dcMatches) - length(matchedCoordinates_filt))])
            matchedCoordinates_filt2=smoothFilt_2024_02_20(matchedCoordinates_filt,fineSmoothVars(1),fineSmoothVars(2));  
            filt_Disps2 = matchedCoordinates_filt2(:,6:8) - matchedCoordinates_filt2(:,2:4);
            %disp(['Filtered Beads = ',num2str(length(matchedCoordinates_filt) - length(matchedCoordinates_filt2))])
            track{t-1}{j} = matchedCoordinates_filt2;                         
                  
        end
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n'); 
end

end