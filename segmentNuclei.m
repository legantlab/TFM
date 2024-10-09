%Test script for segmenting nuclei
%Load data
fold = 'T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA_Sub';
maskOut = 'T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA_Sub_Masks';
tiffs = dir([fold,'\*.tiff']);
allNucAspects = zeros(length(tiffs),1);
allVolume = zeros(length(tiffs),1);
allDiam = zeros(length(tiffs),1);
masks = cell(length(tiffs),1);
axisLengths = zeros(length(tiffs),3);
allProps = cell(length(tiffs),1);
for i = 1:length(tiffs)
    
    I = readtiff([tiffs(i).folder,'\',tiffs(i).name]);
    I = I - 525; %Background subtract images
%     imagesc(I(:,:,20))
    I(I<0) = 0;
%     axis equal
    I = imresize3(I,'scale',[1,1,4]); % Resize to make volume isotropic
    
    % Segment nuclei  
    BW = imbinarize(I,'global');
%     imagesc(BW(:,:,20))
    %cleanup
    SE = strel('disk',5);
    clean = imerode(BW,SE);
    %clean = imerode(BW,SE);
    clean = imdilate(clean,SE);
    %clean = imdilate(clean,SE); %Multible smaller steps rather than one big one for improved performance
    clean = imfill(clean,"holes");
    
%     figure
%     imagesc(clean(:,:,15))
    maskProps = regionprops3(clean,'all');
    if ~isempty(maskProps)
        [~,volMax] = max(maskProps.Volume);
        nucProps = maskProps(volMax,:);
        %allNucAspects(i) = max(nucProps.PrincipalAxisLength)/min(nucProps.PrincipalAxisLength);
        allNucAspects(i) = (nucProps.PrincipalAxisLength(1))/(nucProps.PrincipalAxisLength(2));
        allVolume(i) = nucProps.Volume * (.2^3);
        masks{i} = nucProps.ConvexImage;
        allDiam(i) = nucProps.EquivDiameter;
        axisLengths(i,:) = size(nucProps.ConvexImage{1});
        allProps{i} = nucProps;
    end
    
    %Output the result

    % Define the output filename for the multi-page TIFF stack
    outputFilename = [maskOut,'\mask',sprintf('%03d',i),'.tiff']; % Specify your desired output filename
    % Create a Tiff object for writing the stack
    t = Tiff(outputFilename, 'w');
    % Get the size of the Mask - type your array's name
    [numRows, numCols, numSlices] = size(clean);
    % Loop through each slice and add it to the multi-page TIFF stack
    for slice = 1:numSlices
        % Extract the 2D slice
        sliceData = clean(:, :, slice);
        
        % Write the current slice to the multi-page TIFF stack
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('BitsPerSample', 1); % Set BitsPerSample to 32 for single precision data
        t.setTag('SamplesPerPixel', 1);
        %t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
        t.setTag('ImageLength', numRows);
        t.setTag('ImageWidth', numCols);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.write(sliceData);
        
        % Add a new page to the multi-page TIFF stack
        if slice ~= numSlices
            t.writeDirectory();
        end
    end
    % Close the TIFF file
    t.close();

end

%% Further analysis
allSA = zeros(length(allProps),1);
nucAspects = zeros(length(allProps),3);
for i = 1:40
    curProp = allProps{i};
    allSA(i) = curProp.SurfaceArea;
    axisLength(i,:) = curProp.PrincipalAxisLength;
    curImg = curProp.Image{1};
    xy = regionprops(squeeze(max(curImg,[],3)),'MajorAxisLength','MinorAxisLength');
    xyAspect = xy.MajorAxisLength/xy.MinorAxisLength;
    xz = regionprops(squeeze(max(curImg,[],2)),'MajorAxisLength','MinorAxisLength');
    xzAspect = xz.MajorAxisLength/xz.MinorAxisLength;
    yz = regionprops(squeeze(max(curImg,[],1)),'MajorAxisLength','MinorAxisLength');
    yzAspect = yz.MajorAxisLength/yz.MinorAxisLength;
    nucAspects(i,:) = [xyAspect, xzAspect, yzAspect];

    convexImg = curProp.ConvexImage{1};
    % Define the output filename for the multi-page TIFF stack
    outputFilename = [maskOut,'\mask',sprintf('%03d',i),'.tiff']; % Specify your desired output filename
    % Create a Tiff object for writing the stack
    t = Tiff(outputFilename, 'w');
    % Get the size of the Mask - type your array's name
    [numRows, numCols, numSlices] = size(convexImg);
    % Loop through each slice and add it to the multi-page TIFF stack
    for slice = 1:numSlices
        % Extract the 2D slice
        sliceData = convexImg(:, :, slice);
        
        % Write the current slice to the multi-page TIFF stack
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('BitsPerSample', 1); % Set BitsPerSample to 32 for single precision data
        t.setTag('SamplesPerPixel', 1);
        %t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
        t.setTag('ImageLength', numRows);
        t.setTag('ImageWidth', numCols);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.write(sliceData);
        
        % Add a new page to the multi-page TIFF stack
        if slice ~= numSlices
            t.writeDirectory();
        end
    end
    % Close the TIFF file
    t.close();

end
figure
plot(nucAspects(:,1))
hold on 
plot(nucAspects(:,2))
plot(nucAspects(:,3))
legend('XY','XZ','YZ')

% 
% exampleNuc = allProps{10}.ConvexImage{1};
% exProps = regionprops3(exampleNuc,'all');

%% Now try processing all of the data
% saveZStackByTime('T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\C3-F06.tif','T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\DNA')
% saveZStackByTime('T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\C3-F08.tif','T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\DNA')
% saveZStackByTime('T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\C3-F11.tif','T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\DNA')

F00 = "T:\Max\2023-10-08\Tiffs\F00\DNA";
F01 = "T:\Max\2023-10-08\Tiffs\F01\DNA";
F02 = "T:\Max\2023-10-08\Tiffs\F02\DNA";
F05 = "T:\Max\2023-10-08\Tiffs\F05\DNA";
%ON1
F06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\DNA";
F08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\DNA";
F09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\DNA";
F11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\DNA";
F15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA";

folds = [F00, F01, F02, F05, F06, F08, F09, F11, F15];
allPropsSum = cell(length(folds),1);

for j = 1:length(folds)
    curFold = folds(j);
    tiffs = dir([convertStringsToChars(curFold),'\*.tiff']);
    allProps = cell(length(tiffs),1);
        for i = 1:length(tiffs)
            
            I = readtiff([tiffs(i).folder,'\',tiffs(i).name]);
            I = I - 525; %Background subtract images
        %     imagesc(I(:,:,20))
            I(I<0) = 0;
        %     axis equal
            I = imresize3(I,'scale',[1,1,4]); % Resize to make volume isotropic
            
            % Segment nuclei  
            BW = imbinarize(I,'global');
        %     imagesc(BW(:,:,20))
            %cleanup
            SE = strel('disk',5);
            clean = imerode(BW,SE);
            %clean = imerode(BW,SE);
            clean = imdilate(clean,SE);
            %clean = imdilate(clean,SE); %Multible smaller steps rather than one big one for improved performance
            clean = imfill(clean,"holes");
            
        %     figure
        %     imagesc(clean(:,:,15))
            maskProps = regionprops3(clean,'all');
            if ~isempty(maskProps)
                [~,volMax] = max(maskProps.Volume);
                nucProps = maskProps(volMax,:);
                allProps{i} = nucProps;
            end
            %Output the result
        
     
        
        end
allPropsSum{j} = allProps;
end

%% Process Confine Data
allPropsAspect = cell(length(allPropsSum),1);
for i = 1:length(allPropsSum)
    curProps = allPropsSum{i};
    curPropsAspect = zeros(length(curProps),1);
    curFold = folds(i); 
    mkdir([convertStringsToChars(curFold),'\Automasks'])
    for j = 1:length(curProps)
        curTime = curProps{j};
        if ~isempty(curTime)
            %Retrieve properties of the mask
            curPropsAspect(j) = curTime.PrincipalAxisLength(1)/curTime.PrincipalAxisLength(3);
            curMask = curTime.Image{1};
            
%             % Define the output filename for the multi-page TIFF stack
%             outputFilename = [convertStringsToChars(curFold),'\Automasks\',sprintf('%03d',j),'.tiff']; % Specify your desired output filename
%             % Create a Tiff object for writing the stack
%             t = Tiff(outputFilename, 'w');
%             % Get the size of the Mask - type your array's name
%             [numRows, numCols, numSlices] = size(curMask);
%             % Loop through each slice and add it to the multi-page TIFF stack
%             for slice = 1:numSlices
%                 % Extract the 2D slice
%                 sliceData = curMask(:, :, slice);
%                 
%                 % Write the current slice to the multi-page TIFF stack
%                 t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
%                 t.setTag('Compression', Tiff.Compression.None);
%                 t.setTag('BitsPerSample', 1); % Set BitsPerSample to 32 for single precision data
%                 t.setTag('SamplesPerPixel', 1);
%                 %t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
%                 t.setTag('ImageLength', numRows);
%                 t.setTag('ImageWidth', numCols);
%                 t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%                 t.write(sliceData);
%                 
%                 % Add a new page to the multi-page TIFF stack
%                 if slice ~= numSlices
%                     t.writeDirectory();
%                 end
%             end
%             % Close the TIFF file
%             t.close();



        end
    end
    allPropsAspect{i} = curPropsAspect;
end


%% Repeat but for straight data
% saveZStackByTime('T:\Max\2023-12-20\Tiffs\F02\DNA-F02_Cropped_Registered.tif','T:\Max\2023-12-20\Tiffs\F02\DNA')
% saveZStackByTime('T:\Max\2023-12-20\Tiffs\F09\DNA_F09.tif','T:\Max\2023-12-20\Tiffs\F09\DNA')
% saveZStackByTime('T:\Max\2023-12-20\Tiffs\F13\DNA-F13_Registered_Croppped.tif','T:\Max\2023-12-20\Tiffs\F13\DNA')

DNAF01 = "T:\Max\2023-12-20\Tiffs\F01\DNA";
DNAF02 = "T:\Max\2023-12-20\Tiffs\F02\DNA";
DNAF09 = "T:\Max\2023-12-20\Tiffs\F09\DNA";
DNAF13 = "T:\Max\2023-12-20\Tiffs\F13\DNA";
DNAF14 = "T:\Max\2023-12-20\Tiffs\F14\DNA";

DNAF00_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\DNA";
DNAF13_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F13\DNA";
DNAF14_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F14\DNA";

folds = [DNAF01; DNAF02; DNAF09; DNAF13; DNAF14; DNAF00_01; DNAF13_01; DNAF14_01];
allPropsSumStraight = cell(length(folds),1);

for j = 1:length(folds)
    curFold = folds(j);
    tiffs = dir([convertStringsToChars(curFold),'\*.tiff']);
    allProps = cell(length(tiffs),1);
        for i = 1:length(tiffs)
            
            I = readtiff([tiffs(i).folder,'\',tiffs(i).name]);
            I = I - 525; %Background subtract images
        %     imagesc(I(:,:,20))
            I(I<0) = 0;
        %     axis equal
            I = imresize3(I,'scale',[1,1,4]); % Resize to make volume isotropic
            
            % Segment nuclei  
            BW = imbinarize(I,'global');
        %     imagesc(BW(:,:,20))
            %cleanup
            SE = strel('disk',5);
            clean = imerode(BW,SE);
            %clean = imerode(BW,SE);
            clean = imdilate(clean,SE);
            %clean = imdilate(clean,SE); %Multible smaller steps rather than one big one for improved performance
            clean = imfill(clean,"holes");
            
        %     figure
        %     imagesc(clean(:,:,15))
            maskProps = regionprops3(clean,'all');
            if ~isempty(maskProps)
                [~,volMax] = max(maskProps.Volume);
                nucProps = maskProps(volMax,:);
                allProps{i} = nucProps;
            end
            %Output the result
        
     
        
        end
allPropsSumStraight{j} = allProps;
end

%% Process Straight Channel data
allPropsAspectStraight = cell(length(allPropsSumStraight),1);
for i = 1:length(allPropsSumStraight)
    curProps = allPropsSumStraight{i};
    curPropsAspect = zeros(length(curProps),1);
    curFold = folds(i); 
    mkdir([convertStringsToChars(curFold),'\Automasks'])
    for j = 1:length(curProps)
        curTime = curProps{j};
        if ~isempty(curTime)
            %Retrieve properties of the mask
            curPropsAspect(j) = curTime.PrincipalAxisLength(1)/curTime.PrincipalAxisLength(3);
            curMask = curTime.Image{1};
            
%             % Define the output filename for the multi-page TIFF stack
%             outputFilename = [convertStringsToChars(curFold),'\Automasks\',sprintf('%03d',j),'.tiff']; % Specify your desired output filename
%             % Create a Tiff object for writing the stack
%             t = Tiff(outputFilename, 'w');
%             % Get the size of the Mask - type your array's name
%             [numRows, numCols, numSlices] = size(curMask);
%             % Loop through each slice and add it to the multi-page TIFF stack
%             for slice = 1:numSlices
%                 % Extract the 2D slice
%                 sliceData = curMask(:, :, slice);
%                 
%                 % Write the current slice to the multi-page TIFF stack
%                 t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
%                 t.setTag('Compression', Tiff.Compression.None);
%                 t.setTag('BitsPerSample', 1); % Set BitsPerSample to 32 for single precision data
%                 t.setTag('SamplesPerPixel', 1);
%                 %t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
%                 t.setTag('ImageLength', numRows);
%                 t.setTag('ImageWidth', numCols);
%                 t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%                 t.write(sliceData);
%                 
%                 % Add a new page to the multi-page TIFF stack
%                 if slice ~= numSlices
%                     t.writeDirectory();
%                 end
%             end
%             % Close the TIFF file
%             t.close();



        end
    end
    allPropsAspectStraight{i} = curPropsAspect;
end

