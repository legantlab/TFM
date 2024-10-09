function [x, track, translation] = funTPT_V2(varargin)
%
% [x, track] = funTPT(fileInfo, beadParameter, tptParam) is the main
% function that performs the single particle tracking on time increment of
% volumetric images.
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


%Parse inputs
[fileNames, beadParameter, tptParameter, isCZI] = parseInputs(varargin{:});
nChannel = 1; %Number of multi-attribute particles
%threshold = beadParameter{1}.intensity;
bandpass_size = beadParameter{1}.winSize;
threshold = beadParameter{1}.thres;
threshmulti = beadParameter{1, 1}.thresh_multi;
%% Particle Tracking
% Handle first time point which is the reference state seperately; Need to
% apply drift correction and set so that later iterations can be run in
% parrallel. 
x = cell(length(fileNames),1);
track = cell(length(fileNames),1);
translation = cell(length(fileNames),1);

t = 2;
j = 1;

if isCZI
    data = bfopen(fileNames{t-1,1});
    omeMeta = data{1, 4};
    stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue();
    hyperstackSize = length(data{1});
    numChannels = hyperstackSize/stackSizeZ;
    I = zeros([size(data{1}{1,1}) stackSizeZ]);
    index = 1;
    for i = fileNames{t-1,2}+1:numChannels:hyperstackSize
        I(:,:,index) = data{1}{i,1};
        index = index+1;
    end
else
    I = loadtiff(fileNames{t-1,1});
    %Do DC correction
    I = I - median(I,'all');
    I(I<0) = 0;
end
        
% Parameters
ImgDim=size(I);
%Default bpassParams %%% 

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
%maybe do a DC subtraction? 
[IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3}); %Bandpass filter - Can tune parameters for optimal fitting
%for some reason this becomes a single, change back to uint16
%IBP = uint16(IBP);
bpGUI = bPassGUI(I,IBP,bPassParams,threshold,threshmulti);
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
         
x{1}{j} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold);

x{1}{j}(any(isnan(x{1}{j}),2),:) = [];


%Center at origin
% translation{1} = mean(x{1}{1},1);
% x{1}{1} = x{1}{1} - mean(x{1}{1},1);

            
%% Correct for swelling induced by SDS

disp('Correcting for SDS induced Swelling')

if isCZI
    data = bfopen(fileNames{t,1});
    omeMeta = data{1, 4};
    stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue();
    hyperstackSize = length(data{1});
    numChannels = hyperstackSize/stackSizeZ;
    I = zeros([size(data{1}{1,1}) stackSizeZ]);
    index = 1;
    for i = fileNames{t,2}+1:numChannels:hyperstackSize
        I(:,:,index) = data{1}{i,1};
        index = index+1;
    end
else
    I = loadtiff(fileNames{t,1});
    I = I - median(I,'all');
    I(I<0) = 0;
end

x{t}{j} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold);
x{t}{j}(any(isnan(x{t}{j}),2),:) = [];
%Center at 0
translation{t} = mean(x{t}{j},1);
%translation{t} = [ImgDim(2)/2, ImgDim(1)/2, ImgDim(3)/2];
x{t}{j} =   x{t}{j} - translation{t};

translation{1} = translation{t};
x{1}{1} = x{1}{1} - translation{t};

figure
scatter(x{1}{1}(:,1), x{1}{1}(:,2))
hold on 
scatter(x{2}{1}(:,1), x{2}{1}(:,2))
%Use Pointcloud registration to align data set to reference. 
% PC1 = pointCloud(x{1}{j});
% PC2 = pointCloud(x{t}{j});
% tform = pcregistericp(PC1, PC2);
% PC2reg = pctransform(PC2,tform);
% x{t}{j} = PC2reg.Location;      

swellingGUI = swellingCorrectionGUI(x{1}{1},x{2}{1});

while(swellingGUI.closenum ~= 1)
    pause(1)
end

x{1}{1} = swellingGUI.finalBeads;
x{1}{1}(any(isnan(x{1}{1}),2),:) = [];

swellingGUI.delete

% Track using new scripts
refBeads = [(1:length(x{1}{1}))',x{1}{1}];
strainBeads = [(1:length(x{t}{1}))',x{t}{1}];
%Use Pointcloud registration to align data set to reference. 
% PC1 = pointCloud(refBeads(:,2:4));
% PC2 = pointCloud(strainBeads(:,2:4));
% tform = pcregistericp(PC1, PC2);
% PC1reg = pctransform(PC1,tform);
% refBeads(:,2:4) = PC1reg.Location;   

firstPass = displacementCalc_recursive_Max_2024_02_19(refBeads,strainBeads,25,3,5,2);

distNewMatches = firstPass(:,6:8) - firstPass(:,2:4);

xDrift = median(distNewMatches(:,1));
yDrift = median(distNewMatches(:,2));
zDrift = median(distNewMatches(:,3));
figure
quiver(firstPass(:,2), firstPass(:,3), distNewMatches(:,1), distNewMatches(:,2),1)
axis equal

refBeadsCor = refBeads;
%translation{t} = {xDrift, yDrift, zDrift};
strainBeadsCor = strainBeads;
%refBeadsCor(:,2) = refBeadsCor(:,2) - xDrift;
%refBeadsCor(:,3) = refBeadsCor(:,3) - yDrift;
%refBeadsCor(:,4) = refBeadsCor(:,4) - zDrift;
% strainBeadsCor(:,2) = strainBeadsCor(:,2) - xDrift;
% strainBeadsCor(:,3) = strainBeadsCor(:,3) - yDrift;
% strainBeadsCor(:,4) = strainBeadsCor(:,4) - zDrift;
% figure
% scatter(strainBeadsCor(:,2), strainBeadsCor(:,3))
% hold on 
% scatter(refBeadsCor(:,2), refBeadsCor(:,3))
% axis equal

dcMatches = displacementCalc_recursive_Max_2024_02_19(refBeadsCor,strainBeadsCor,50,3,5,1.5);
testDCmatches = dcMatches(:,6:8) - dcMatches(:,2:4);
figure 
quiver(dcMatches(:,2), dcMatches(:,3), testDCmatches(:,1), testDCmatches(:,2),1)
matchedCoordinates_filt=smoothFilt_2024_02_20(dcMatches,50,1.75); 
filt_Disps = matchedCoordinates_filt(:,6:8) - matchedCoordinates_filt(:,2:4);
disp(['Filtered Beads = ',num2str(length(dcMatches) - length(matchedCoordinates_filt))])
matchedCoordinates_filt2=smoothFilt_2024_02_20(matchedCoordinates_filt,50,1.5); 

filt_Disps2 = matchedCoordinates_filt2(:,6:8) - matchedCoordinates_filt2(:,2:4);
disp(['Filtered Beads = ',num2str(length(matchedCoordinates_filt) - length(matchedCoordinates_filt2))])

figure
quiver(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3), filt_Disps2(:,1), filt_Disps2(:,2),1)
axis equal

track{t-1}{j} = matchedCoordinates_filt2;


% Start Tracking
parfor t = 3:length(fileNames) % Loop through all time points
    
    tStart = tic;
    
    
    for j = 1:nChannel % Loop through different types of multi-attribute particles
        
        % Detect and Localize Particles -----------------------------------
        tPP = tic;

       
        if t~=length(fileNames) % after first frame except last frame
            disp(['Current time point: t = ' num2str(t)])
            disp(['  Particle channel: ' num2str(fileNames{t,2})])
            disp(['  Current Filename: ' fileNames{t}])
            
            if isCZI
                data = bfopen(fileNames{t,1});
                omeMeta = data{1, 4};
                stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue();
                hyperstackSize = length(data{1});
                numChannels = hyperstackSize/stackSizeZ;
                I = zeros([size(data{1}{1,1}) stackSizeZ]);
                index = 1;
                for i = fileNames{t,2}+1:numChannels:hyperstackSize
                    I(:,:,index) = data{1}{i,1};
                    index = index+1;
                end
            else
                I = loadtiff(fileNames{t,1});
                I = I - median(I,'all');
                I(I<0) = 0;
            end
            
            %localize particles
            x{t}{j} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold); %may need to remove the parfor loop in here. 
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

            firstPass = displacementCalc_recursive_Max_2024_02_19(refBeads,strainBeads,20,3,5,1.5);

            distNewMatches = firstPass(:,6:8) - firstPass(:,2:4);
            
            % figure
            % quiver(firstPass(:,6), firstPass(:,7),distNewMatches(:,1), distNewMatches(:,2),0)
            % hold on 
            % scatter(firstPass(:,2), firstPass(:,3), 'r','.')
            % scatter(firstPass(:,6), firstPass(:,7), 'g','.')
            % axis equal
            
            %Implement a scheme where we compute the drift and apply that to the raw
            %displacements then retrack...
            xDrift = median(distNewMatches(:,1));
            yDrift = median(distNewMatches(:,2));
            zDrift = median(distNewMatches(:,3));
            refBeadsCor = refBeads;
            %refBeadsCor(:,2) = refBeadsCor(:,2) + xDrift;
            %refBeadsCor(:,3) = refBeadsCor(:,3) + yDrift;
            %refBeadsCor(:,4) = refBeadsCor(:,4) + zDrift;
            %translation{t} = {xDrift, yDrift, zDrift};
            strainBeadsCor = strainBeads;
            strainBeadsCor(:,2) = strainBeadsCor(:,2) - xDrift;
            strainBeadsCor(:,3) = strainBeadsCor(:,3) - yDrift;
            strainBeadsCor(:,4) = strainBeadsCor(:,4) - zDrift;
            
            % figure
            % scatter(refBeadsCor(:,2), refBeadsCor(:,3), 'r','.')
            % hold on 
            % scatter(strainBeadsCor(:,2), strainBeadsCor(:,3), 'g','.')
            % axis equal
            
            dcMatches = displacementCalc_recursive_Max_2024_02_19(refBeadsCor,strainBeadsCor,25,3,5,1.4);
            % Plot resulting matches
%             figure
%             scatter(dcMatches(:,2), dcMatches(:,3),'.','r')
%             hold on
%             scatter(dcMatches(:,6), dcMatches(:,7),'.','g')
%             dcDisps = dcMatches(:,6:8) - dcMatches(:,2:4);
%             quiver(dcMatches(:,2), dcMatches(:,3),dcDisps(:,1), dcDisps(:,2),0)
%             axis equal
            
            % Run smoothing function twice for better smoothing (course and
            % fine grain)
            
            matchedCoordinates_filt=smoothFilt_2024_02_20(dcMatches,50,1.75); 
            filt_Disps = matchedCoordinates_filt(:,6:8) - matchedCoordinates_filt(:,2:4);
            %disp(['Filtered Beads = ',num2str(length(dcMatches) - length(matchedCoordinates_filt))])
            matchedCoordinates_filt2=smoothFilt_2024_02_20(matchedCoordinates_filt,50,1.50); 
            
            filt_Disps2 = matchedCoordinates_filt2(:,6:8) - matchedCoordinates_filt2(:,2:4);
            %disp(['Filtered Beads = ',num2str(length(matchedCoordinates_filt) - length(matchedCoordinates_filt2))])

%             figure
%             scatter(refBeadsCor(:,2), refBeadsCor(:,3), 'r','o')
%             hold on 
%             scatter(strainBeadsCor(:,2), strainBeadsCor(:,3), 'g','o')
%             scatter(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3), 'r','.')
%             hold on
%             scatter(matchedCoordinates_filt2(:,6), matchedCoordinates_filt2(:,7), 'g','.')
%             quiver(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3),filt_Disps2(:,1), filt_Disps2(:,2),0)
%             axis equal

            track{t-1}{j} = matchedCoordinates_filt2;                         
%           track{t-1}{j} = TPT(x{1}{j},x{t}{j},finaltptParams{j},predictor);
          
        end 
        
    end
    
    disp(['Total Elaspsed Time = ', num2str(toc(tStart)),' seconds']);fprintf('\n'); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = parseInputs(varargin)
% varargout = parseInputs(filename,beadParameter, psptParameter)


%%% Parse filenames
filename = varargin{1};

%%% Detection and Localization Parameters
beadParameter = varargin{2};

% Define default values
thres = 0.5;
minSize = 5;
maxSize = Inf;
winSize = [7,7,7];
dccd = [1,1,1];
abc = [1,1,1];
forloop = 1;
randNoise = 1/10^7; % Something small
xy = 1; % Assume symmetrical if not given
z = 1;  % Assume symmetrical if not givenf
diameter = 5;   % Assume if not given
intensity = 600;
thresh_multi = 1.5;

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'thres',thres);
    addParameter(p,'minSize',minSize);
    addParameter(p,'maxSize',maxSize);
    addParameter(p,'winSize',winSize);
    addParameter(p,'dccd',dccd);
    addParameter(p,'abc',abc);
    addParameter(p,'forloop',forloop);
    addParameter(p,'randNoise',randNoise);
    addParameter(p,'xy',xy);
    addParameter(p,'z',z);
    addParameter(p,'diameter',diameter);
    addParameter(p,'intensity',intensity);
    addParameter(p,'thresh_multi',thresh_multi);
    
    parse(p,beadParameter{i})
    
    beadParameter{i} = p.Results;
    
end

%%% TPT Parameters
tptParameter = varargin{3};

% Define default values/We redefine these after the bpassGUI
knnFD =16;
knnFM = 4;
fmThres = 2;
maxIter = 14;
nSpheres = 2;
outlrThres = 5;

for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'knnFD',knnFD);
    addParameter(p,'knnFM',knnFM);
    addParameter(p,'fmThres',fmThres);
    addParameter(p,'maxIter',maxIter);
    addParameter(p,'nSpheres',nSpheres);
    addParameter(p,'outlrThres',outlrThres);
    
    parse(p,tptParameter{i})
    
    tptParameter{i} = p.Results;
    
end

%%% Outputs
varargout{1} = filename;
varargout{2} = beadParameter;
varargout{3} = tptParameter;
varargout{4} = varargin{4};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = loadFile(fileInfo,idx,randNoise)
I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I),
    if numel(I)==1, I = I{1};
    else
        I = I{fileInfo.datachannel};
    end
end

I = double(I);
I = I/max(I(:));
end