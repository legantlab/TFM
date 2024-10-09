%Script to read in traction data from multible data sets and compile into
%cell arrays for analysis. 
%ON2
F00 = "T:\Max\2023-12-20\Tiffs\F00\2024-01-16_ComputedTractionsAlignedValuesCropped.mat";
F01 = "T:\Max\2023-12-20\Tiffs\F01\2024-04-28_ComputedTractionsAlignedValues.mat";
F02 = "T:\Max\2023-12-20\Tiffs\F02\2024-04-28_ComputedTractionsAlignedValues.mat";
F03 = "T:\Max\2023-12-20\Tiffs\F03\2024-01-22_ComputedTractionsAlignedValuesCropped.mat";
F09 = "T:\Max\2023-12-20\Tiffs\F09\2024-04-28_ComputedTractionsAlignedValues.mat";
F12 = "T:\Max\2023-12-20\Tiffs\F12\2024-01-22_ComputedTractionsAlignedValuesCropped.mat";
F13 = "T:\Max\2023-12-20\Tiffs\F13\2024-04-28_ComputedTractionsAlignedValues.mat";
F14 = "T:\Max\2023-12-20\Tiffs\F14\2024-04-28_ComputedTractionsAlignedValues.mat";
F15 = "T:\Max\2023-12-20\Tiffs\F15\2024-01-15_ComputedTractionsAlignedValuesCropped.mat";
F16 = "T:\Max\2023-12-20\Tiffs\F16\2024-01-15_ComputedTractionsAlignedValuesCropped.mat";

F00_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\2024-04-28_ComputedTractionsAlignedValues.mat";
F13_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F13\2024-04-28_ComputedTractionsAlignedValues.mat";
F14_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F14\2024-04-28_ComputedTractionsAlignedValues.mat";
%files = [F00; F01; F02; F03; F09; F12; F13; F14; F15; F16];
%files = [F00; F01; F02; F03; F09; F12; F13; F14;F15;F16];
files = [F01; F02; F09; F13; F14; F00_01; F13_01; F14_01];
%% Load Masks and apply to traction fields
compiledMasks = cell(length(files),1);
MasksF00 = "T:\Max\2023-12-20\Tiffs\F00\LA_Masks";
MasksF01 = "T:\Max\2023-12-20\Tiffs\F01\LA_Masks";
MasksF02 = "T:\Max\2023-12-20\Tiffs\F02\LA_Masks";
MasksF03 = "T:\Max\2023-12-20\Tiffs\F03\LA_02_Masks";
MasksF09 = "T:\Max\2023-12-20\Tiffs\F09\LA_Masks";
MasksF12 = "T:\Max\2023-12-20\Tiffs\F12\LA_Masks";
MasksF13 = "T:\Max\2023-12-20\Tiffs\F13\LA_Masks";
MasksF14 = "T:\Max\2023-12-20\Tiffs\F14\LA_Masks";
MasksF15 = "T:\Max\2023-12-20\Tiffs\F15\LA_Masks";
MasksF16 = "T:\Max\2023-12-20\Tiffs\F16\LA_Masks";

MasksF00_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\LA_Masks";
MasksF13_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F13\LA_Masks";
MasksF14_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F14\LA_Masks";
%maskDir = [MasksF00;MasksF01;MasksF02;MasksF03;MasksF09;MasksF12;MasksF13;MasksF14; MasksF15; MasksF16];
maskDir = [MasksF01;MasksF02;MasksF09;MasksF13;MasksF14; MasksF00_01; MasksF13_01; MasksF14_01];
for j = 1:length(maskDir)
    curDir = char(maskDir(j));
    
    tiffs=dir([curDir,'\*.tif']);
    tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
    currentAllMasks = cell(length(tiffs),1);
    for i = 1:length(tiffs(:,1))
           curMask = [tiffs(i).folder,'\',tiffs(i).name];
           currentAllMasks{i} = loadtiff(curMask);
    
    end
    compiledMasks{j} = currentAllMasks;
end

% Do the same for DNA masks
compiledDNA = cell(length(files),1);
DNAF00 = "T:\Max\2023-12-20\Tiffs\F00\DNA_Masks";
DNAF01 = "T:\Max\2023-12-20\Tiffs\F01\DNA_Masks";
DNAF02 = "T:\Max\2023-12-20\Tiffs\F02\DNA_Masks";
DNAF03 = "T:\Max\2023-12-20\Tiffs\F03\DNA_02_Masks";
DNAF09 = "T:\Max\2023-12-20\Tiffs\F09\DNA_Masks";
DNAF12 = "T:\Max\2023-12-20\Tiffs\F12\DNA_Masks";
DNAF13 = "T:\Max\2023-12-20\Tiffs\F13\DNA_Masks";
DNAF14 = "T:\Max\2023-12-20\Tiffs\F14\DNA_Masks";
DNAF15 = "T:\Max\2023-12-20\Tiffs\F15\DNA_Masks";
DNAF16 = "T:\Max\2023-12-20\Tiffs\F16\DNA_Masks";

DNAF00_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\DNA_Masks";
DNAF13_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\DNA_Masks";
DNAF14_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\DNA_Masks";

%DNADir = [DNAF00; DNAF01; DNAF02; DNAF03; DNAF09; DNAF12; DNAF13; DNAF14; DNAF15; DNAF16];
DNADir = [DNAF01; DNAF02; DNAF09; DNAF13; DNAF14; DNAF00_01; DNAF13_01; DNAF14_01];
translations = cell(length(DNADir),1);
NucCentroids = cell(length(DNADir),1);
for j = 1:length(DNADir)
    curDir = char(DNADir(j));
    
    tiffs=dir([curDir,'\*.tif']);
    tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
    DNAAllMasks = cell(length(tiffs),1);
    curDNACentroids = zeros(length(tiffs),2);
    for i = 1:length(tiffs(:,1))
           curMask = [tiffs(i).folder,'\',tiffs(i).name];
           DNAAllMasks{i} = loadtiff(curMask);
           MaskStats = regionprops(bwareafilt(logical(DNAAllMasks{i}),1),'all');
           if isempty(MaskStats)
               curDNACentroids(i,:) = nan;
           else
               curDNACentroids(i,:) = MaskStats.Centroid;
           end
    
    end
    compiledDNA{j} = DNAAllMasks;
    ImgSize = size(curMask);
    translations{j} = [ImgSize(2)/2, ImgSize(1)/2];
    NucCentroids{j} = curDNACentroids;
end
 % We will need to alter this as we are no longer doing translations the
 % same way, instead we can find the translation we need just from the size
 % of the mask. 
%Load localizations to get the translation for each data set
% DispF00 = "T:\Max\2023-10-08\Tiffs\F00\5x40_ON2_F00_Full.mat";
% DispF01 = "T:\Max\2023-10-08\Tiffs\F01\5x40_ON2_F01_Full.mat";
% DispF02 = "T:\Max\2023-10-08\Tiffs\F02\5x40_ON2_F02_Full.mat";
% DispF05 = "T:\Max\2023-10-08\Tiffs\F05\5x40_ON2_F05_Full.mat";
% DispF12 = "T:\Max\2023-10-08\Tiffs\F12\5x40_ON2_F12_Full.mat";
% 
% DispF06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\5x40_ON1_F06_Swell_Corrected.mat";
% DispF08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\5x40_ON1_F08_Swell_Corrected.mat";
% DispF09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\5x40_ON1_F09_Swell_Corrected.mat";
% DispF11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\5x40_ON1_F11_Swell_Corrected.mat";
% DispF15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\5x40_ON1_F15_Cropped.mat";
% 
% dispfiles = [DispF00, DispF01, DispF02, DispF05, DispF12, DispF06, DispF08, DispF09, DispF11, DispF15];
% translations = cell(length(dispfiles),1);
% 
% for i = 1:length(dispfiles)
%     load(dispfiles(i))
%     translations{i} = trans;
% end


%% Load data
F00 = "T:\Max\2023-12-20\Tiffs\F00\2024-01-16_ComputedTractionsAlignedValuesCropped.mat";
F01 = "T:\Max\2023-12-20\Tiffs\F01\2024-04-28_ComputedTractionsAlignedValues.mat";
F02 = "T:\Max\2023-12-20\Tiffs\F02\2024-04-28_ComputedTractionsAlignedValues.mat";
F03 = "T:\Max\2023-12-20\Tiffs\F03\2024-01-22_ComputedTractionsAlignedValuesCropped.mat";
F09 = "T:\Max\2023-12-20\Tiffs\F09\2024-04-28_ComputedTractionsAlignedValues.mat";
F12 = "T:\Max\2023-12-20\Tiffs\F12\2024-01-22_ComputedTractionsAlignedValuesCropped.mat";
F13 = "T:\Max\2023-12-20\Tiffs\F13\2024-04-28_ComputedTractionsAlignedValues.mat";
F14 = "T:\Max\2023-12-20\Tiffs\F14\2024-04-28_ComputedTractionsAlignedValues.mat";
F15 = "T:\Max\2023-12-20\Tiffs\F15\2024-01-15_ComputedTractionsAlignedValuesCropped.mat";
F16 = "T:\Max\2023-12-20\Tiffs\F16\2024-01-15_ComputedTractionsAlignedValuesCropped.mat";

F00_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F00\2024-04-28_ComputedTractionsAlignedValues.mat";
F13_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F13\2024-04-28_ComputedTractionsAlignedValues.mat";
F14_01 = "S:\Max\2024_03_21_MicrochannelTFM\StraightOvernight\2024-03-21\tiffs\F14\2024-04-28_ComputedTractionsAlignedValues.mat";
%files = [F00; F01; F02; F03; F09; F12; F13; F14; F15; F16];
%files = [F00; F01; F02; F03; F09; F12; F13; F14;F15;F16];
files = [F01; F02; F09; F13; F14; F00_01; F13_01; F14_01];

displacements = cell(length(files),1);
tractions = cell(length(files),1);
%rawLocalizations = cell(length(files),1);
for i = 1:length(files(:,1))
    load(files(i))
    displacements{i} = beadDisps;
    tractions{i} = traction_vector;
    %rawLocalizations{i} = matches;
end

%% Compute Traction quantifications
timeArray = [1,1,1,1];

% Get the area of each element for computing force balance
fileName = 'T:\Max\2023-11-30\Tiffs\HR_Straight\StraightChannel_Coarsened_Set_Loaded_SolverDeck.inp';
data = abaqusInpRead(fileName);
[newIds, TR2,newNodeCoords, surfNodes,elemCents2Dnew] = organizeGeometry(data,1,1);
elemList = TR2.ConnectivityList;
elemCoords = TR2.Points;
elemAreas = zeros(length(elemList),1);
for i = 1:length(elemList)
    curElem = elemList(i,:);
    p1 = elemCoords(curElem(1),:);
    p2 = elemCoords(curElem(2),:);
    p3 = elemCoords(curElem(3),:);
    
    V1 = p2 - p1;
    V2 = p3 - p1;
    elemAreas(i) = 0.5*norm(cross(V1,V2));
end

%Strain Energy
E = 20000; %Elastic modulus of the substrate, should then give "force" in Pa
compiledSE = cell(length(files),1);
compiledSEDen = cell(length(files),1);
compiledForces = cell(length(files),1);
compiledRelSE = cell(length(files),1);
%alignedTimes = cell(length(files),1);
dilation = strel('diamond',50);
compiledDistances = cell(length(files),1);
compiledNucAR = cell(length(files),1);
compiledSol = cell(length(files),1);
compiledOri = cell(length(files),1);
compiledPeri = cell(length(files),1);
compiledDistHull = cell(length(files),1);
compiledLADistHull = cell(length(files),1);
compiledDistancesUnsort = cell(length(files),1);
compiledCellArea = cell(length(files),1);
compiledNucArea = cell(length(files),1);
compiledNucCellDist = cell(length(files),1);
compiledCellCirc = cell(length(files),1);
compiledNCM = cell(length(files),1);


for j = 1:length(files)
traction_vector = tractions{j};
correctedbeadDisps = displacements{j};

%Apply Masks here
curMasks = compiledMasks{j}; 
curDNAMasks = compiledDNA{j}; 

strainEnergy = zeros(length(traction_vector),1);
vecForce = zeros(length(traction_vector{1}(:,1)),1);
strainEnergyDen = zeros(length(traction_vector),1);
relStrainEnergy = zeros(length(traction_vector),1);
distanceCenter = zeros(length(traction_vector),1);
currentTranslation = translations{j};
nucAspect = zeros(length(traction_vector),1);
solidity = zeros(length(traction_vector),1);
orientation = zeros(length(traction_vector),1);
perimeter = zeros(length(traction_vector),1);
NCM = zeros(length(traction_vector),1);
distanceHullExtent = zeros(length(traction_vector),1);
distanceLAHullExtent = zeros(length(traction_vector),1);
nucArea = zeros(length(traction_vector),1);
cellArea = zeros(length(traction_vector),1);
cellNucDist = zeros(length(traction_vector),1);
cellCirc = zeros(length(traction_vector),1);

NCMx = zeros(length(traction_vector),1);
NCMy = zeros(length(traction_vector),1);
%NCMz = zeros(length(traction_vector),1);
    %Compute hull of the original mask and compute distance to center of
    %channel. 
    curNucCentroids = NucCentroids{j};
    curNucCentroids = reshape(curNucCentroids(~isnan(curNucCentroids)),[],2); %remove nans?
    %We can just assume the center is just the center of the path?
    channelCenter = [(min(curNucCentroids(:,1)) + max(curNucCentroids(:,1)))/2 ,...
        (min(curNucCentroids(:,2)) + max(curNucCentroids(:,2)))/2];

    % curNucCentroids = curNucCentroids(1:end-1,:);
    % centDists = squareform(pdist(curNucCentroids));
    % [centMax, centMaxIdxX] = max(centDists);
    % [~, centMaxIdxY] = max(centMax);
    % centMaxIdxX = centMaxIdxX(centMaxIdxY);
    % 
    % channelCenter = [(curNucCentroids(centMaxIdxX,1)+curNucCentroids(centMaxIdxY,1))/2, ...
    %     (curNucCentroids(centMaxIdxX,2)+curNucCentroids(centMaxIdxY,2))/2];%Is this the most reasonable approximation of the centroid of the channel?

% 
%     firstMask = curDNAMasks{1};
%     lastMask = curDNAMasks{end-1};
% 
%     firstMaskStats = regionprops(bwareafilt(logical(firstMask),1),'all');
%     lastMaskStats = regionprops(bwareafilt(logical(lastMask),1),'all');
%     firstMaskCent = firstMaskStats.Centroid;
%     lastMaskCent = lastMaskStats.Centroid;
%     channelCenter = [(lastMaskCent(1)+firstMaskCent(1))/2,(firstMaskCent(2)+lastMaskCent(2))/2];

for i = 1:length(strainEnergy)
    curMask = curMasks{i};
    
    curMaskdil = imdilate(curMask, dilation);
    curDNAmask = curDNAMasks{i};
    imgSize = size(curMaskdil);
    elemCentroids = elemCents2Dnew(:,1:2);
    elemCentroids(:,1) = round(elemCentroids(:,1)/.199 + currentTranslation(1) - 1);
    elemCentroids(:,2) = round(elemCentroids(:,2)/.199 + currentTranslation(2) - 1);

    [ycoords, xcoords] = find(curMaskdil);
    coordinateMask = ismember(elemCentroids,[xcoords,ycoords],'rows');
    maskedCoords = elemCentroids(:,:);
    
    %channelCenter = [0,0];%mean(elemCentroids); We really should define the center of the channel as the mid point of the migration

    LAmaskStats = regionprops(bwareafilt(logical(curMask),1),'all');

    if isempty(LAmaskStats) %this isn't quite working because its computing the outline which isn't necessarilly correct. 
        distCenterLAHullmin = 100000;
        cellArea(i) = 0;
        cellCent = [0,0];
        cellCirc(i) = 0;
         
    else
        cellArea(i) = LAmaskStats.Area;
        cellCirc(i) = LAmaskStats.Circularity;
        LAmaskHull = LAmaskStats.ConvexHull;
        cellCent = LAmaskStats.Centroid;
        distCenterLAHullmin = inf;
        for k = 1:length(LAmaskHull)
                curPoint = LAmaskHull(k,:);
                distCenterLAHull = (sqrt((curPoint(1) - channelCenter(1))^2 + (curPoint(2) - channelCenter(2))^2))*.199;
                if distCenterLAHull < distCenterLAHullmin
                        distCenterLAHullmin = distCenterLAHull;
                end
        end
    end

    distanceLAHullExtent(i) = distCenterLAHullmin;
    

    maskStats = regionprops(bwareafilt(logical(curDNAmask),1),'all');
    

    if isempty(maskStats)
        maskCentroid = [0,0];
        nucAspect(i) = 0;
        solidity(i) = 0;
        orientation(i) = 0;
        perimeter(i) = 0;
        distCenterHullmin = 100000;
        nucArea(i) =0;
    else
        maskCentroid = maskStats.Centroid;
        nucArea(i) = maskStats.Area;
        nucAspect(i) = maskStats.MajorAxisLength/maskStats.MinorAxisLength;
        solidity(i) = maskStats.Solidity;
        orientation(i) = maskStats.Orientation;
        perimeter(i) = maskStats.Perimeter;


        maskHull = maskStats.ConvexHull;
        distCenterHullmin = inf;
        for k = 1:length(maskHull)
                curPoint = maskHull(k,:);
                distCenterHull = (sqrt((curPoint(1) - channelCenter(1))^2 + (curPoint(2) - channelCenter(2))^2))*.199;
                if distCenterHull < distCenterHullmin
                        distCenterHullmin = distCenterHull;
                end
        end
    end
        
    distCenter = (sqrt((maskCentroid(1) - channelCenter(1))^2 + (maskCentroid(2) - channelCenter(2))^2))*.199;
    
    distanceHullExtent(i) = distCenterHullmin;

    distanceCenter(i) = distCenter;

    curDisps = correctedbeadDisps{i}(:,4:6);
    %curDisps = correctedbeadDisps{i}(coordinateMask,4:6);
    magDisps = vecnorm(curDisps,2,2);
    curForce = traction_vector{i}(:,4:6) * E; 
    %curForce = traction_vector{i}(coordinateMask,4:6) * E;
    magForce = vecnorm(curForce,2,2);
    
    SE = 0.5 * sum(magDisps.*magForce,"all");
    strainEnergy(i) = SE;
    strainEnergyDen(i) = SE/cellArea(i);
    vecForce(:,i) = magForce;
    %vecForce(coordinateMask,i) = magForce;

    %Compute distance between nucleus and cell centroid
    cellNucDist(i) = sqrt((maskCentroid(1) - cellCent(1))^2 + (maskCentroid(2) - cellCent(2))^2);
    %Compute Net Contractile Moment
    %Compute Net Contractile Moment - Redo as a function of cell centroid
    %in 2D
    centerOfMesh = mean(elemCents2Dnew);
    %centeredMesh = elemCents2Dnew(:,1:2) - cellCent*.199;

    centeredMesh = sqrt(((elemCents2Dnew(:,1:2) - (maskCentroid*.199)).^2));
        %Compute forces from current traction using element areas
    curForceArea = curForce.*elemAreas*10^-6;
%     conMoments = [sum(centeredMesh(:,1).*curForceArea(:,1)), sum(centeredMesh(:,1).*curForceArea(:,2)), sum(centeredMesh(:,1).*curForceArea(:,3));...
%         sum(centeredMesh(:,2).*curForceArea(:,1)),sum(centeredMesh(:,2).*curForceArea(:,2)),sum(centeredMesh(:,2).*curForceArea(:,3));...
%         sum(centeredMesh(:,3).*curForceArea(:,1)),sum(centeredMesh(:,3).*curForceArea(:,2)),sum(centeredMesh(:,3).*curForceArea(:,3))];

    conMoments = [sum(centeredMesh(:,1).*curForceArea(:,1)), sum(centeredMesh(:,1).*curForceArea(:,2));...
        sum(centeredMesh(:,2).*curForceArea(:,1)),sum(centeredMesh(:,2).*curForceArea(:,2))];
    diagConMoments = diag(conMoments);
    NCMx(i) = diagConMoments(1);
    NCMy(i) = diagConMoments(2);
%Compute xyz distances from cell center, will need to use the whole set of
%of the image? We will need the XYZ centroid from each image which is
%nontrivial without some sort of segmentation. 
% LAmask = regionprops(bwareafilt(logical(curMask),1));
% 
% curCentroid = LAmask.Centroid;
% %convert to units
% curCentroid = (curCentroid -1 - currentTranslation (1:2))*.199;
% curDisps(:,1:2) = curDisps(:,1:2) - curCentroid;
% curNCMmat = [sum(curDisps(:,1).*curForce(:,1)), sum(curDisps(:,1).*curForce(:,2));...
%     sum(curDisps(:,2).*curForce(:,1)), sum(curDisps(:,2).*curForce(:,2))]; 
% curNCMmat = [sum(curDisps(:,1).*curForce(:,1)), sum(curDisps(:,1).*curForce(:,2)), sum(curDisps(:,1).*curForce(:,3));...
%     sum(curDisps(:,2).*curForce(:,1)), sum(curDisps(:,2).*curForce(:,2)), sum(curDisps(:,2).*curForce(:,3));...
%     sum(curDisps(:,3).*curForce(:,1)), sum(curDisps(:,3).*curForce(:,2)), sum(curDisps(:,3).*curForce(:,3))];

end
%relStrainEnergy = (strainEnergy - min(strainEnergy(strainEnergy>0)))/(max(strainEnergy) - min(strainEnergy(strainEnergy>0)));
relStrainEnergy = strainEnergy/(max(strainEnergy));
strainEnergyDen(strainEnergyDen==0) = 0;
compiledSEDen{j} = strainEnergyDen;
compiledRelSE{j} = relStrainEnergy;
compiledSE{j} = strainEnergy;
compiledForces{j} = vecForce;
compiledNucAR{j} = nucAspect;
compiledSol{j} = solidity;
compiledOri{j} = orientation;
compiledPeri{j} = perimeter;
compiledDistancesUnsort{j} = distanceCenter;
compiledCellArea{j} = cellArea;
compiledNucArea{j} = nucArea;
compiledNucCellDist{j} = cellNucDist;
compiledCellCirc{j} = cellCirc;
compiledNCM{j} = [NCMx,NCMy];
%Create linspace to align times
%alignedTimes{j} = linspace(1,length(strainEnergy), length(strainEnergy)) - timeArray(j);
%set all values before the minimum to negative, redo for local min? Works
%for now with our existing data but may want to check if we have multible
%transit events, etc. 
[~,minIdx] = min(distanceCenter);
distanceCenter(1:minIdx) = distanceCenter(1:minIdx)*-1;

compiledDistances{j} = distanceCenter;

%Do similar trick for hull distances
[~, minHullIdx] = min(distanceHullExtent);
distanceHullExtent(1:minHullIdx) = distanceHullExtent(1:minHullIdx)*-1;
compiledDistHull{j} = distanceHullExtent;

[~, minLAHullIdx] = min(distanceLAHullExtent);
distanceLAHullExtent(1:minLAHullIdx) = distanceLAHullExtent(1:minLAHullIdx)*-1;
compiledLADistHull{j} = distanceLAHullExtent;
end

%% Plotting
close all
scale = 1%10^-6;
%Relative Strain Energies from strain in middle of channel
sampleGrid = linspace(-50,50,11);
figure
hold on 
compiledInterpSE = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpRelSE = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpSEDen = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpNAR = zeros(length(sampleGrid),length(compiledRelSE));
compiledHullInterpSE = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpNARHull = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpOriHull = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpOri = zeros(length(sampleGrid),length(compiledRelSE));
compiledTransitTimes = zeros(1,length(compiledRelSE));
compiledAvgStrainTransit = zeros(1,length(compiledRelSE));
compiledMaxStrainTransit = zeros(1,length(compiledRelSE));
compiledInterpNucCellDist = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpCellCirc = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpCellArea = zeros(length(sampleGrid),length(compiledRelSE));
compiledNucVel = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpCellVel = zeros(length(sampleGrid),length(compiledRelSE));

compiledInterpNCMx = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpNCMy = zeros(length(sampleGrid),length(compiledRelSE));
%compiledInterpNCMz = zeros(length(sampleGrid),length(compiledRelSE));

for i = 1:length(compiledSE)
    plot(compiledDistances{i},smooth(compiledRelSE{i})*scale,'MarkerSize',10,'Marker','.')
    
    
    %Relative Peak Strain Energy
    [sortDisps, sortIdx] = sort(compiledDistances{i});
    sortDisps = sortDisps+(rand(size(sortDisps,1),1)*.000001);
   
    interpolatedRelSE = interp1(sortDisps,(compiledRelSE{i}(sortIdx)),sampleGrid);
    compiledInterpRelSE(:,i) = interpolatedRelSE;
    %Strain Energy
    interpolatedSE = interp1(sortDisps,(compiledSE{i}(sortIdx)),sampleGrid);
    compiledInterpSE(:,i) = interpolatedSE;

    %Strain Energy Density
    interpolatedSEDen = interp1(sortDisps,(compiledSEDen{i}(sortIdx)),sampleGrid);
    compiledInterpSEDen(:,i) = interpolatedSEDen;

    %Nuclear Aspect Ratio
    interpolatedNAR = interp1(sortDisps, compiledNucAR{i}(sortIdx),sampleGrid);
    compiledInterpNAR(:,i) = interpolatedNAR;

    %SE vs Minimum Hull Distance
    [sortHullDisps,sortHullIdx] = sort(compiledDistHull{i});
    sortHullDisps = sortHullDisps+(rand(size(sortHullDisps,1),1)*.000001);

    interpolatedSE2 = interp1(sortHullDisps,(compiledRelSE{i}(sortHullIdx)),sampleGrid);
    compiledHullInterpSE(:,i) = interpolatedSE2;

    %Nuclear Aspect Ratio vs min Hull Distance
    interpolatedNARHull = interp1(sortHullDisps, compiledNucAR{i}(sortHullIdx),sampleGrid);
    compiledInterpNARHull(:,i) = interpolatedNARHull;

    %Interpolated orientation vs min Hull Distance
    interpolatedOriHull = interp1(sortHullDisps, abs(compiledOri{i}(sortHullIdx)),sampleGrid);
    compiledInterpOriHull(:,i) = interpolatedOriHull;

    %Interpolated orientation vs centroid distance
    interpolatedOri = interp1(sortDisps, abs(compiledOri{i}(sortIdx)),sampleGrid);
    compiledInterpOri(:,i) = interpolatedOri;

    %Transit Time computation, when centroid passes -20 to 20 for each time
    %point
    curDistances = compiledDistances{i};
    timeVals = (linspace(1,length(curDistances) ,length(curDistances)) - 1)*5;
    transitTime = interp1(curDistances + (rand(size(curDistances,1),1)*.000001),timeVals,[-20,20]);
    %transitTime = transitTime(2) - transitTime(1);
    compiledTransitTimes(i) = transitTime(2) - transitTime(1);
    %Now compute average/max strain energy between this period
    transitPeriod = linspace(transitTime(1), transitTime(2),100);
    interpolatedStrainTransit = interp1(compiledDistances{i}+ (rand(size(compiledDistances{i},1),1)*.000001),compiledRelSE{i},transitPeriod);
    compiledAvgStrainTransit(i) = mean(interpolatedStrainTransit,'omitnan');
    compiledMaxStrainTransit(i) = max(interpolatedStrainTransit,[],'omitnan');

    %Also size of the cells during transit? Or nucleus? 
    %Need to do interpolation to determine areas, etc. 
    interpolatedCellArea = interp1(sortHullDisps,compiledCellArea{i},sampleGrid);
    compiledInterpCellArea(:,i) = (interpolatedCellArea-min(interpolatedCellArea))/(max(interpolatedCellArea) - min(interpolatedCellArea));
    %Compute distance between nuclear and cell centroids? 
    interpolatedCellNucDist = interp1(sortHullDisps, compiledNucCellDist{i}(sortIdx),sampleGrid);
    %normalize
    %interpolatedCellNucDist = (interpolatedCellNucDist-min(interpolatedCellNucDist))/(max(interpolatedCellNucDist) - min(interpolatedCellNucDist));
    compiledInterpNucCellDist(:,i) = interpolatedCellNucDist;

    %Compiled Cell Eccentricity
    interpolatedCellCirc = interp1(sortHullDisps, compiledCellCirc{i}(sortIdx),sampleGrid);
    compiledInterpCellCirc(:,i) = interpolatedCellCirc;

    %Cell momement? Nuclear moments? 

    %Compute velocity of cell/nucleus as a function of strain? 
    curHullDists = compiledDistances{i};
    nucVel = diff(curHullDists);
    compiledInterpCellVel(:,i) = interp1(sortHullDisps(1:end-1), nucVel,sampleGrid);

     %NCM interpolation
    interpolatedNCMx = interp1(sortHullDisps, compiledNCM{i}(sortIdx,1),sampleGrid);
    compiledInterpNCMx(:,i) = interpolatedNCMx;

    interpolatedNCMy = interp1(sortHullDisps, compiledNCM{i}(sortIdx,2),sampleGrid);
    compiledInterpNCMy(:,i) = interpolatedNCMy;

%     interpolatedNCMz = interp1(sortHullDisps, compiledNCM{i}(sortIdx,3),sampleGrid);
%     compiledInterpNCMz(:,i) = interpolatedNCMz;
end

xline(-20)
xline(20)
ylabel('Relative Strain Energy');
xlabel('Nuclear Distance to Center of Confinement (um)')
%xlabel('Time Relative to Center of Channel (mins)')
legend('F00','F01','F02','F05','F12','F06','F08','F09','F11','F15', 'F00_01', "F13_01", "F14_01")


figure %Histogram bins of data not interpolated
allDisps = [];
allSE = [];
for j = 1:length(compiledSE)
    curDists = round(compiledDistances{j}/10);
    curSE = compiledSE{j};

    uniqueVals = unique(curDists);
    
    curAvg = [];
    for k = 1:length(uniqueVals)
        curUnique = uniqueVals(k);
        curAvg(k) = mean(curSE(curDists==curUnique));
    end

%     allDisps = [allDisps;curDists];
%     allSE = [allSE;curSE];
    
    allDisps = [allDisps;uniqueVals];
    allSE = [allSE;curAvg'];

end

boxplot(allSE, (allDisps)*10)
hold on 
swarmchart(((allDisps+7)),allSE,10,'r','filled')
ylim([0, 5*10^5])
%xlim([2.5 13.5])
xlabel('Nuclear Distance to Center of Confinement (um)')
ylabel('Strain Energy (pJ)')

avgRelSE = mean(compiledInterpRelSE,2,'omitnan');
stdRelSE = std(compiledInterpRelSE,[],2,'omitnan')/sqrt(size(compiledInterpRelSE(1,:),2));
figure
errorbar(sampleGrid,avgRelSE,stdRelSE,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgRelSE,'r','linewidth',1)
ylabel('Relative Peak Strain Energy')
xlabel('Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

avgSE = mean(compiledInterpSE,2,'omitnan');
stdSE = std(compiledInterpSE,[],2,'omitnan')/sqrt(size(compiledInterpSE(1,:),2));
figure
errorbar(sampleGrid,avgSE,stdSE,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgSE,'r','linewidth',1)
ylabel('Strain Energy')
xlabel('Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

avgSEDen = mean(compiledInterpSEDen,2,'omitnan');
stdSEDen = std(compiledInterpSEDen,[],2,'omitnan')/sqrt(size(compiledInterpSEDen(1,:),2));
figure
errorbar(sampleGrid,avgSEDen,stdSEDen,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgSEDen,'r','linewidth',1)
ylabel('Strain Energy Density')
xlabel('Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

avgNAR = mean(compiledInterpNAR,2,'omitnan');
errNUC = std(compiledInterpNAR,[],2,'omitnan')/sqrt(size(compiledInterpNAR(1,:),2));
figure
errorbar(sampleGrid,avgNAR,errNUC,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNAR, 'r','linewidth',1)
ylabel('Average Nuclear Aspect Ratio')
xlabel('Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Interpolated Hull distances plot
avgHullSE = mean(compiledHullInterpSE,2,'omitnan');
errHullSE = std(compiledHullInterpSE,[],2,'omitnan')/sqrt(size(compiledHullInterpSE(1,:),2));
figure
errorbar(sampleGrid,avgHullSE,errHullSE,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgHullSE,'r','linewidth',1)
ylabel('Relative Peak Strain Energy')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Interpolated Hull distances plot NAR
avgHullNAR = mean(compiledInterpNARHull,2,'omitnan');
errHullNAR = std(compiledInterpNARHull,[],2,'omitnan')/sqrt(size(compiledInterpNARHull(1,:),2));
figure
errorbar(sampleGrid,avgHullNAR,errHullNAR,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgHullNAR,'r','linewidth',1)
ylabel('Nuclear Aspect Ratio')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Orientation of nucleus vs distance to confinement
avgHullOri = mean(compiledInterpOriHull,2,'omitnan');
errHullOri = std(compiledInterpOriHull,[],2,'omitnan')/sqrt(size(compiledInterpOriHull(1,:),2));
figure
errorbar(sampleGrid,avgHullOri,errHullOri,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgHullOri,'r','linewidth',1)
ylabel('Nuclear Orientation')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Orientation of nucleus vs Centroid distance to confinement

avgOri = mean(compiledInterpOri,2,'omitnan');
errOri = std(compiledInterpOri,[],2,'omitnan')/sqrt(size(compiledInterpOri(1,:),2));
figure
errorbar(sampleGrid,avgOri,errOri,'b','LineWidth',0.5)
hold on 
plot(sampleGrid,avgOri,'r','linewidth',1)
ylabel('Nuclear Orientation')
xlabel('Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

% Transit time versus maximum SE in channel? Or average SE in channel
figure
scatter(compiledTransitTimes, compiledAvgStrainTransit)
xlabel('Confinement Transit Times (mins)')
ylabel('Average Relative Strain Energy during Transit')

figure
scatter(compiledTransitTimes, compiledMaxStrainTransit)
xlabel('Confinement Transit Times (mins)')
ylabel('Max Relative Strain Energy during Transit')

%Nuclear Cell Centroid Distances
avgNucCellDist = mean(compiledInterpNucCellDist,2,'omitnan');
errNucCellDist = std(compiledInterpNucCellDist,[],2,'omitnan')/sqrt(size(compiledInterpNucCellDist(1,:),2));
figure
errorbar(sampleGrid,avgNucCellDist,errNucCellDist,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNucCellDist, 'r', 'LineWidth',1)
ylabel('Nuclear/Cell Centroid Distance (um)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Cell eccentricity
avgCellEccent = mean(compiledInterpCellCirc,2,'omitnan');
errCellEccent = std(compiledInterpCellCirc,[],2,'omitnan')/sqrt(size(compiledInterpCellCirc(1,:),2));
figure
errorbar(sampleGrid,avgCellEccent,errCellEccent,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgCellEccent, 'r', 'LineWidth',1)
ylabel('Cell Circularity')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Cell area during transit? 
avgCellArea = mean(compiledInterpCellArea,2,'omitnan');
errCellArea = std(compiledInterpCellArea,[],2,'omitnan')/sqrt(size(compiledInterpCellArea(1,:),2));
figure
errorbar(sampleGrid,avgCellArea,errCellArea,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgCellArea, 'r', 'LineWidth',1)
ylabel('Normalized Cell Area')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Cell velocity versus nuclear distance to center
avgNucVel = mean(compiledInterpCellVel,2,'omitnan');
errCellVel = std(compiledInterpCellVel,[],2,'omitnan')/sqrt(size(compiledInterpCellVel(1,:),2));
figure
errorbar(sampleGrid,avgNucVel,errCellVel,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNucVel, 'r', 'LineWidth',1)
ylabel('Average Cell Velocity (um/min)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%Net Contractile Moments
avgNCMx = mean(compiledInterpNCMx,2,'omitnan');
errNCMx = std(compiledInterpNCMx,[],2,'omitnan')/sqrt(size(compiledInterpNCMx(1,:),2));
figure
errorbar(sampleGrid,avgNCMx,errNCMx,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNCMx, 'r', 'LineWidth',1)
ylabel('Average NCMx (pNm)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

avgNCMy = mean(compiledInterpNCMy,2,'omitnan');
errNCMy = std(compiledInterpNCMy,[],2,'omitnan')/sqrt(size(compiledInterpNCMy(1,:),2));
figure
errorbar(sampleGrid,avgNCMy,errNCMy,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNCMy, 'r', 'LineWidth',1)
ylabel('Average NCMy (pNm)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

% avgNCMz = mean(compiledInterpNCMz,2,'omitnan');
% errNCMz = std(compiledInterpNCMz,[],2,'omitnan')/sqrt(size(compiledInterpNCMz(1,:),2));
% figure
% errorbar(sampleGrid,avgNCMz,errNCMz,'b','LineWidth',0.5)
% hold on 
% plot(sampleGrid, avgNCMz, 'r', 'LineWidth',1)
% ylabel('Average NCMz (pNm)')
% xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
% xline(-20)
% xline(20)

%Combined NCM graph
figure
errorbar(sampleGrid,avgNCMx,errNCMx,'LineWidth',0.5)
hold on 
plot(sampleGrid, avgNCMx, 'LineWidth',1)

errorbar(sampleGrid,avgNCMy,errNCMy,'LineWidth',0.5)
hold on 
plot(sampleGrid, avgNCMy, 'LineWidth',1)

% errorbar(sampleGrid,avgNCMz,errNCMz,'LineWidth',0.5)
% hold on 
% plot(sampleGrid, avgNCMz, 'LineWidth',1)

ylabel('Average NCM (pNm)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%% Compile measurements into a cell array for storage and comparison
% 
% save('T:\Max\2023-12-20\Tiffs\IA32LAStraightValues.mat',"sampleGrid","compiledInterpSE","compiledInterpRelSE","compiledInterpSEDen","compiledInterpNAR", "compiledHullInterpSE", "compiledInterpNARHull", "compiledInterpOriHull",...
%     "compiledInterpOri", "compiledTransitTimes", "compiledAvgStrainTransit", "compiledMaxStrainTransit", "compiledInterpNucCellDist",...
%     "compiledInterpCellCirc", "compiledInterpCellArea", "compiledNucVel" ,"compiledInterpCellVel")

 %% Plot individual traces if needed
% 
% figure
% for i = 6:10%length(compiledSE)
%     plot(compiledDistances{i}, compiledRelSE{i})
%     hold on 
% 
% end
% xlim([-25,25])
% legend('F00','F01','F02','F03','F09','F12','F13','F14','F15','F16')
%% Cross quantifications

% %Average NAR vs Average Relative Strain Energy
% figure
% scatter(avgNAR,avgRelSE)

%%     Compute distributions of force
% 
% %loop through force/traction vectors
% for i = 1:length(tractions)
%     curForce = tractions{i};
%     for j = 1:length(curForce)
%         curForceMap = curForce{j};
%         figure
%         scatter(curForceMap(:,2),curForceMap(:,4),'.')
%         hold on 
%         scatter(curForceMap(:,2),curForceMap(:,5),'.')
%         scatter(curForceMap(:,2),curForceMap(:,6),'.')
%         legend('X Traction', 'Y Traction', 'Z Traction ')
%         close()
% 
%         XInterp = interp1(curForceMap(:,2), curForceMap(:,4),sampleGrid);
%         YInterp = interp1(curForceMap(:,2), curForceMap(:,5),sampleGrid);
%         ZInterp = interp1(curForceMap(:,2), curForceMap(:,6),sampleGrid);
%         plot(sampleGrid,XInterp)
%         hold on 
%         plot(sampleGrid,YInterp)
%         plot(sampleGrid,ZInterp)
%         legend('X Traction', 'Y Traction', 'Z Traction ')
%     end
%     
% 
% end