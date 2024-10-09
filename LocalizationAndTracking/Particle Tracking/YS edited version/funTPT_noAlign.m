function [x, track] = funTPT(varargin)
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

% Start Tracking
for t = 2:length(fileNames) % Loop through all time points
    
    tStart = tic;
    
    %On the first time point, run the bandpassGUI to assess the bpass
    %parameters we want to use for the rest of the images. Also add
    %threshold slider as well to set threshold. Perhaps run the entire set
    %of localization on the beads and do the check here for the first
    %image set. 
    
    
    for j = 1:nChannel % Loop through different types of multi-attribute particles
        
        % Detect and Localize Particles -----------------------------------
        tPP = tic;

        if t == 2 % timePoint = 1
            disp(['Current time point: t = ' num2str(t-1)])
            disp(['  Particle channel: ' num2str(fileNames{t-1,2})])
            disp(['  Current Filename: ' fileNames{t-1}])
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
            %[newRef, trans] = alignOrigin(x{1}{1},0);
%             load('averagedSwellingProfile7x40_expDerivedProfilePixel.mat')
           
                     
%             %Plot
            
% 
%   
%             hold on 
%             %align element centroids
%             [newElements, ~] = alignOrigin(elemCents2D,0);
%             scatter3(newElements(:,1), newElements(:,2), newElements(:,3))
%             %interpolate reference to uniform grid - actually this makes no
%             %sense without a displacement field
% %             elemCents2D = newCentroids;
% %             elemCents2D(:,1:2) = elemCents2D(:,1:2)/.199;
% %             elemCents2D(:,3) = elemCents2D(:,3)/.8;
% %             averageDisplacement = finalDisps;
% %             averageDisplacement(:,1:2) = averageDisplacement(:,1:2)/.199;
% %             averageDisplacement(:,3) = averageDisplacement(:,3)/.8;
%             displacedbeads = newElements + averageDisplacement; 
% %             elemCents2D(:,3) = elemCents2D(:,3)*1.025;
%             [corBeads] = interpDisps({[newElements, displacedbeads]},newRef);
%             corBeads2 = corBeads{1}(:,1:3) - corBeads{1}(:,4:6);
%             scatter3(corBeads2(:,1), corBeads2(:,2), corBeads2(:,3))
%             [corBeads3, ~] = alignOrigin(corBeads2,0); %Reallign to center after swelling correction? 


            %x{1}{1} = newRef;
            %x{1}{1} = newRef;
            
            %% Correct for swelling induced by SDS
            
            disp('Correcting for SDS induced Swelling')
            % load('T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\AlignedSwellingData.mat') %Now make this process general
            % referenceBeads = x{1}{1};
            % newReference = interpDisps({[finalElements(:,1),finalElements(:,2), finalElements(:,3), ...
            %         finalElements(:,1) + finalDisplacements(:,1)...
            %         , finalElements(:,2) + finalDisplacements(:,2)...
            %         , finalElements(:,3) + finalDisplacements(:,3)]},referenceBeads);
            % scale = 9;%%1.75; % We might want some metric to compare the correct against to fine tune this scaling factor. 
            % x{1}{1} = [newReference{1}(:,1) + (scale*newReference{1}(:,4)), newReference{1}(:,2) + ...
            %         (scale*newReference{1}(:,5)), newReference{1}(:,3) + (scale*0.35*newReference{1}(:,6))];
            
            %Now load another dataset to set tracking/swelling parameters
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

            %add some random jitter in z for flat substrates.

            x{t}{j}(:,3) = x{t}{j}(:,3) + rand(length(x{t}{j}(:,3)),1);
            
            %[newFrame, ~] = alignOrigin(x{t}{1},0);
            
            
            %newFrame(:,2) = newFrame(:,2) - 8;
            
            %x{t}{1} = newFrame;
%             figure
%             scatter3(newRef(:,1), newRef(:,2), newRef(:,3))
%             hold on 
%             scatter3(newframe(:,1)-2.5, newframe(:,2), newframe(:,3))
%             Swelling Correction GUI

%Manual swelling correction - manually correct for now then update the
%swellingCorrectionGUI to allow for dynamic scaling
% %Load average swelling profile
% load('averagedSwellingProfile7x40_expDerivedProfilePixel.mat')
% 
% %Move to origin
% [newCents, ~] = alignOrigin(elemCents2D,0);
% 
% figure
% scatter3(newCents(:,1),newCents(:,2), newCents(:,3))
% hold on 
% scatter3(newCents(:,1) + averageDisplacement(:,1)...
%                     , newCents(:,2) + averageDisplacement(:,2)...
%                     , newCents(:,3) + averageDisplacement(:,3))
% %Interpolate displacements to a bead locations
% newReference = interpDisps({[newCents(:,1),newCents(:,2), newCents(:,3), ...
%                     newCents(:,1) + averageDisplacement(:,1)...
%                     , newCents(:,2) + averageDisplacement(:,2)...
%                     , newCents(:,3) + averageDisplacement(:,3)]},newRef);
%                %Apply displacements and plot
%                corRef = [newRef(:,1) - newReference{1}(:,4), ...
%                    newRef(:,2) - newReference{1}(:,5), ...
%                    newRef(:,3) - 1.5*newReference{1}(:,6)];
%                close all
%                figure
%                scatter3(newframe(:,1)-2.5, newframe(:,2), newframe(:,3))
%                hold on 
%                scatter3(corRef(:,1), corRef(:,2), corRef(:,3))
%                axis equal % So we need scaling it seems, worth and generally applying to the GUI. 
%             x{1}{1} = corRef;
%             x{2}{1} = [newframe(:,1)-2.5, newframe(:,2), newframe(:,3)];

            
            
            % Now do tpt GUI
            % Particle Tracking -----------------------------------------------
            tptParameter{j}.sizeI = size(I); 
            
%             tptGUI = TPTGUI(x,tptParameter{1});
%             %update tptParameter from this
%             %waitfor(tptGUI)
%             while(tptGUI.closenum ~= 1)
%                 pause(1)
%             end
% 
%             finaltptParams = tptGUI.tptParamFinal;
% 
%             tptGUI.delete
            
            if j>1
                predictor.flag = true;
                predictor.x0 = x{t-1}{j-1};
                temp.x1 = x{t}{j-1};
                temp.track = track{t-1}{j-1};
                predictor.u = temp.x1(temp.track(temp.track>0),:)...
                    - predictor.x0(temp.track>0,:);
                predictor.x1 = temp.x1(temp.track(temp.track>0),:);
        
            elseif j == 1
                predictor.flag = false;
                predictor.x0 = [];
                predictor.u = [];
            
            end
            
            track{t-1}{j} = TPT_2d(x{1}{j},x{t}{j}, tptParameter{j},predictor);
        else % after first frame
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
            x{t}{j} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold);
            x{t}{j}(any(isnan(x{t}{j}),2),:) = [];

            %[x{t}{1}] = interpDisps({[newFrame, newFrame]},newCentroids);
            %x{t}{1} = newFrame;
            disp(['    Time to localize particles = ', num2str(toc(tPP)),' seconds']);
            
            if j>1
                predictor.flag = true;
                predictor.x0 = x{t-1}{j-1};
                temp.x1 = x{t}{j-1};
                temp.track = track{t-1}{j-1};
                predictor.u = temp.x1(temp.track(temp.track>0),:)...
                    - predictor.x0(temp.track>0,:);
                predictor.x1 = temp.x1(temp.track(temp.track>0),:);
        
            elseif j == 1
                predictor.flag = false;
                predictor.x0 = [];
                predictor.u = [];
            
            end
        
            %Run TPT
            tptParameter{j}.sizeI = size(I);
            
            track{t-1}{j} = TPT_2d(x{1}{j},x{t}{j},tptParameter{j},predictor);
          
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