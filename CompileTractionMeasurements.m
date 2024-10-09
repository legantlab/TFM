%Script to read in traction data from multible data sets and compile into
%cell arrays for analysis. 
%ON2
F00 = "T:\Max\2023-10-08\Tiffs\F00\2024-04-28_ComputedTractionsAlignedValues.mat";
F01 = "T:\Max\2023-10-08\Tiffs\F01\2024-04-28_ComputedTractionsAlignedValues.mat";
F02 = "T:\Max\2023-10-08\Tiffs\F02\2024-04-28_ComputedTractionsAlignedValues.mat";
F05 = "T:\Max\2023-10-08\Tiffs\F05\2024-04-28_ComputedTractionsAlignedValues.mat";
F12 = "T:\Max\2023-10-08\Tiffs\F12\2024-04-28_ComputedTractionsAlignedValuces.mat";

%ON1
F06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\2024-04-27_ComputedTractionsAlignedValues.mat";
F08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\2024-04-27_ComputedTractionsAlignedValues.mat";
F09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\2024-04-28_ComputedTractionsAlignedValues.mat";
F11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\2024-04-28_ComputedTractionsAlignedValues.mat";
F15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\2024-04-28_ComputedTractionsAlignedValues.mat";
%files = [F00; F01; F02; F05; F06; F08; F09; F11; F15];
files = [F01; F02; F05; F06; F08; F09; F11; F15];
%% Load Masks and apply to traction fields
compiledMasks = cell(length(files),1);

MasksF00 = "T:\Max\2023-10-08\Tiffs\F00\LA_Masks";
MasksF01 = "T:\Max\2023-10-08\Tiffs\F01\LA_Masks";
MasksF02 = "T:\Max\2023-10-08\Tiffs\F02\LA_Masks1";
MasksF05 = "T:\Max\2023-10-08\Tiffs\F05\LA_Masks";
MasksF12 = "T:\Max\2023-10-08\Tiffs\F12\LA_Masks";
MasksF06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\LA_masks";
MasksF08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\LA_masks";
MasksF09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\LA_masks";
MasksF11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\LA_masks";
MasksF15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\LA_masksSub";

%maskDir = [MasksF00; MasksF01; MasksF02; MasksF05; MasksF06; MasksF08; MasksF09; MasksF11; MasksF15];
maskDir = [MasksF01; MasksF02; MasksF05; MasksF06; MasksF08; MasksF09; MasksF11; MasksF15];
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

DNAF00 = "T:\Max\2023-10-08\Tiffs\F00\DNA_Masks";
DNAF01 = "T:\Max\2023-10-08\Tiffs\F01\DNA_Masks";
DNAF02 = "T:\Max\2023-10-08\Tiffs\F02\DNA_Masks1";
DNAF05 = "T:\Max\2023-10-08\Tiffs\F05\DNA_Masks";
DNAF12 = "T:\Max\2023-10-08\Tiffs\F12\DNA_Masks";

DNAF06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\DNA_masks";
DNAF08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\DNA_masks"; 
DNAF09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\DNA_masks";
DNAF11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\DNA_masks";
DNAF15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA_masksSub";

%DNADir = [DNAF00; DNAF01; DNAF02; DNAF05; DNAF06; DNAF08; DNAF09; DNAF11; DNAF15];
DNADir = [DNAF01; DNAF02; DNAF05; DNAF06; DNAF08; DNAF09; DNAF11; DNAF15];
for j = 1:length(DNADir)
    curDir = char(DNADir(j));
    
    tiffs=dir([curDir,'\*.tif']);
    tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
    DNAAllMasks = cell(length(tiffs),1);
    for i = 1:length(tiffs(:,1))
           curMask = [tiffs(i).folder,'\',tiffs(i).name];
           DNAAllMasks{i} = loadtiff(curMask);
    
    end
    compiledDNA{j} = DNAAllMasks;
end

%Load localizations to get the translation for each data set
DispF00 = "T:\Max\2023-10-08\Tiffs\F00\F00_2024Reprocess.mat";
DispF01 = "T:\Max\2023-10-08\Tiffs\F01\F01_2024Reprocess.mat";
DispF02 = "T:\Max\2023-10-08\Tiffs\F02\F02_2024Reprocess.mat";
DispF05 = "T:\Max\2023-10-08\Tiffs\F05\F05_2024Reprocess.mat";
DispF12 = "T:\Max\2023-10-08\Tiffs\F12\F12_2024Reprocess.mat";

DispF06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\F06_2024Reprocess_Smooth.mat";
DispF08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\F08_2024Reprocess_Smooth.mat";
DispF09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\F09_2024Reprocess_Smooth.mat";
DispF11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\F11_2024Reprocess_Smooth.mat";
DispF15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\F15_2024Reprocess_Smooth.mat";

%dispfiles = [DispF00, DispF01, DispF02, DispF05, DispF06, DispF08, DispF09, DispF11, DispF15];
dispfiles = [DispF01, DispF02, DispF05, DispF06, DispF08, DispF09, DispF11, DispF15];
translations = cell(length(dispfiles),1);

for i = 1:length(dispfiles)
    load(dispfiles(i))
    translations{i} = translation;
end


%% Load data
%ON2
F00 = "T:\Max\2023-10-08\Tiffs\F00\2024-04-28_ComputedTractionsAlignedValues.mat";
F01 = "T:\Max\2023-10-08\Tiffs\F01\2024-04-28_ComputedTractionsAlignedValues.mat";
F02 = "T:\Max\2023-10-08\Tiffs\F02\2024-04-28_ComputedTractionsAlignedValues.mat";
F05 = "T:\Max\2023-10-08\Tiffs\F05\2024-04-28_ComputedTractionsAlignedValues.mat";
F12 = "T:\Max\2023-10-08\Tiffs\F12\2024-04-28_ComputedTractionsAlignedValues.mat";

%ON1
F06 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\2024-04-27_ComputedTractionsAlignedValues.mat";
F08 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F08\2024-04-27_ComputedTractionsAlignedValues.mat";
F09 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F09\2024-04-28_ComputedTractionsAlignedValues.mat";
F11 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F11\2024-04-28_ComputedTractionsAlignedValues.mat";
F15 = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\2024-04-28_ComputedTractionsAlignedValues.mat";

%files = [F00; F01; F02; F05; F06; F08; F09; F11; F15];
files = [F01; F02; F05; F06; F08; F09; F11; F15];
displacements = cell(length(files),1);
tractions = cell(length(files),1);
rawLocalizations = cell(length(files),1);
for i = 1:length(files(:,1))
    load(files(i))
    displacements{i} = beadDisps;
    tractions{i} = traction_vector;
    rawLocalizations{i} = matches;
end

%% Compute Traction quantifications
%timeArray = [13, 20,20, 10,  23, 13, 15, 25, 29];
timeArray = [20,20, 10,  23, 13, 15, 25, 29];
% Get the area of each element for computing force balance
fileName = 'T:\Max\2023-10-08\Tiffs\F02\TestData\5x40_Raw_SmoothMesh2p5_s10Course_QA_Tetra_Properties_Loads_SPCs_Loaded_SolverDeck.inp';
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
alignedTimes = cell(length(files),1);
dilation = strel('diamond',10);
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
compiledtractionProj = cell(length(files),1);
compiledxContractile = cell(length(files),1);
compiledxForces5micSum = cell(length(files),1);
compiledxForces5micAvg = cell(length(files),1);
compiledxForces5micMax = cell(length(files),1);
compiledxDisps5micSum = cell(length(files),1);
compiledxDisps5micAvg = cell(length(files),1);
compiledxDisps5micMax = cell(length(files),1);
compiledMaxForce = cell(length(files),1);

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
currentTranslations = translations{j};
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
tractionProj = zeros(length(traction_vector),3);
xContractileSum = zeros(length(traction_vector),1);
force5micSum= zeros(length(traction_vector),1);
force5micAvg = zeros(length(traction_vector),1);
force5micMax = zeros(length(traction_vector),1);
disps5micSum = zeros(length(traction_vector),1);
disps5micAvg = zeros(length(traction_vector),1);
disps5micMax = zeros(length(traction_vector),1);
maxForce = zeros(length(traction_vector),1);

for i = 1:length(strainEnergy)
    curMask = curMasks{i};
    currentTrans = currentTranslations{i};
    

    curMaskdil = imdilate(curMask, dilation);
    curDNAmask = curDNAMasks{i};
    imgSize = size(curMaskdil);
    elemCentroids = elemCents2Dnew(:,1:2);
    elemCentroids(:,1) = round(elemCentroids(:,1)/.199 + currentTrans(1) - 1);
    elemCentroids(:,2) = round(elemCentroids(:,2)/.199 + currentTrans(2) - 1);

    [ycoords, xcoords] = find(curMaskdil);
    coordinateMask = ismember(elemCentroids,[xcoords,ycoords],'rows');
    maskedCoords = elemCentroids(coordinateMask,:);

    %Compute hull of the original mask and compute distance to center of
    %channel. 
    channelCenter = mean(elemCentroids);

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
        nucAspect(i) = nan;
        solidity(i) = 0;
        orientation(i) = 0;
        perimeter(i) = 0;
        distCenterHullmin = 100000;
        nucArea(i) =0;
        maskHull = nan;
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

    %curDisps = correctedbeadDisps{i}(coordinateMask,4:6);
    curDisps = correctedbeadDisps{i}(:,4:6);
    magDisps = vecnorm(curDisps,2,2);
    %curForce = traction_vector{i}(coordinateMask,4:6) * E; 
    curTraction = traction_vector{i}(:,4:6) * E; 
    magForce = vecnorm(curTraction,2,2);
    maxForce(i) = max(magForce);
    SE = 0.5 * sum(magDisps.*magForce,"all");
    strainEnergy(i) = SE;
    strainEnergyDen(i) = SE/cellArea(i);
    %vecForce(coordinateMask,i) = magForce;
    vecForce(:,i) = magForce;
    %Compute distance between nucleus and cell centroid
    cellNucDist(i) = sqrt((maskCentroid(1) - cellCent(1))^2 + (maskCentroid(2) - cellCent(2))^2);
    %Compute Net Contractile Moment - Redo as a function of cell centroid
    %in 2D
    centerOfMesh = mean(elemCents2Dnew);
    %centeredMesh = elemCents2Dnew(:,1:2) - cellCent*.199;

    centeredMesh = sqrt(((elemCents2Dnew(:,1:2) - (maskCentroid*.199)).^2));
        %Compute forces from current traction using element areas
    curForceArea = curTraction.*elemAreas*10^-6;
%     conMoments = [sum(centeredMesh(:,1).*curForceArea(:,1)), sum(centeredMesh(:,1).*curForceArea(:,2)), sum(centeredMesh(:,1).*curForceArea(:,3));...
%         sum(centeredMesh(:,2).*curForceArea(:,1)),sum(centeredMesh(:,2).*curForceArea(:,2)),sum(centeredMesh(:,2).*curForceArea(:,3));...
%         sum(centeredMesh(:,3).*curForceArea(:,1)),sum(centeredMesh(:,3).*curForceArea(:,2)),sum(centeredMesh(:,3).*curForceArea(:,3))];

    conMoments = [sum(centeredMesh(:,1).*curForceArea(:,1)), sum(centeredMesh(:,1).*curForceArea(:,2));...
        sum(centeredMesh(:,2).*curForceArea(:,1)),sum(centeredMesh(:,2).*curForceArea(:,2))];
    diagConMoments = diag(conMoments);
    NCMx(i) = diagConMoments(1);
    NCMy(i) = diagConMoments(2);
    %NCMz(i) = diagConMoments(3);

    %Project force vectors? This is in Pa, maybe only do in significant
    %tractions (>2 STD)
    xProjSTD = std(abs(curTraction(:,1)));
    yProjSTD = std(abs(curTraction(:,2)));
    zProjSTD = std(abs(curTraction(:,3)));
    absTractions = abs(curTraction);

    tractionProj(i,1) = sum(absTractions(:,1));
    tractionProj(i,2) = sum(absTractions(:,2));
    tractionProj(i,3) = sum(absTractions(:,3));

    %Compute contractile tractions based on geometry, neg is towards cell,
    %pos is away
    xNeg = elemCents2Dnew(:,1)<0;
    xTractionsContractile = curTraction(:,1);
    xTractionsContractile(xNeg) = xTractionsContractile(xNeg)*-1;
    xContractileSum(i) = sum(xTractionsContractile);
    
    %Compute forces with 5 microns of the nucleus for summation
    if ~isnan(maskHull)
    [~, distvals] = dsearchn(maskHull, elemCentroids(:,1:2));
    %forceDists = vecnorm(elemCentroids(:,1:2) - maskHull,2,2);
    force5micSum(i) = sum(xTractionsContractile(distvals<=(5/0.199)));
    force5micAvg(i) = mean(xTractionsContractile(distvals<=(5/0.199)));
    force5micMax(i) = max(xTractionsContractile(distvals<=(5/0.199)));

    %Do the same but for displacements. 
    xDispsContractile = curDisps(:,1);
    xDispsContractile(xNeg) = xDispsContractile(xNeg)*-1;
    disps5micSum(i) = sum(xDispsContractile(distvals<=(5/0.199)));
    disps5micAvg(i) = mean(xDispsContractile(distvals<=(5/0.199)));
    disps5micMax(i) = max(xDispsContractile(distvals<=(5/0.199)));
    end
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
strainEnergyDen(strainEnergyDen==0) = 0;
compiledSEDen{j} = strainEnergyDen;
relStrainEnergy = strainEnergy/(max(strainEnergy));
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
compiledtractionProj{j} = tractionProj;
compiledxContractile{j} = xContractileSum;
compiledxForces5micSum{j} = force5micSum;
compiledxForces5micAvg{j} = force5micAvg;
compiledxForces5micMax{j} = force5micMax;
compiledxDisps5micSum{j} = disps5micSum;
compiledxDisps5micAvg{j} = disps5micAvg;
compiledxDisps5micMax{j} = disps5micMax;
compiledMaxForce{j} = maxForce;
    
%Create linspace to align times
alignedTimes{j} = linspace(1,length(strainEnergy), length(strainEnergy)) - timeArray(j);
%set all values before the minimum to negative, redo for local min? Works
%for now with our existing data but may want to check if we have multible
%transit events, etc. 
%Changed this to be the value defined by timeArray - MH 4/24/2024
[~,minIdx] = min(distanceCenter);
distanceCenter(1:timeArray(j)) = distanceCenter(1:timeArray(j))*-1;

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
sampleGrid = [50];%linspace(-50,50,11); %
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
compiledInterpNCM = zeros(length(sampleGrid),length(compiledRelSE));

compiledInterpxProj = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpyProj = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpzProj = zeros(length(sampleGrid),length(compiledRelSE));

compiledInterpxContractile = zeros(length(sampleGrid),length(compiledRelSE));
compiledInterpforces5micSum = zeros(length(sampleGrid), length(compiledRelSE));
compiledInterpforces5micAvg = zeros(length(sampleGrid), length(compiledRelSE));
compiledInterpforces5micMax = zeros(length(sampleGrid), length(compiledRelSE));

compiledInterpdisp5micSum = zeros(length(sampleGrid), length(compiledRelSE));
compiledInterpdisp5micAvg = zeros(length(sampleGrid), length(compiledRelSE));
compiledInterpdisp5micMax = zeros(length(sampleGrid), length(compiledRelSE));

compiledAvgMaxForce = zeros(length(sampleGrid), length(compiledRelSE));

for i = 1:length(compiledSE)
    plot(compiledDistances{i},compiledRelSE{i}*scale,'MarkerSize',10,'Marker','.')
     
    %Relative Peak Strain Energy
    [sortDisps, sortIdx] = sort(compiledDistances{i});
    sortDisps = sortDisps+(rand(size(sortDisps,1),1)*.000001);
   
    interpolatedRelSE = interp1(sortDisps,(compiledRelSE{i}(sortIdx)),sampleGrid);
    compiledInterpRelSE(:,i) = interpolatedRelSE;

    %plot(sampleGrid,interpolatedRelSE)
    %hold on 
    %Strain Energy
    interpolatedSE = interp1(sortDisps,(compiledSE{i}(sortIdx)),sampleGrid);
    %Instead just average the values for each point within each bin.
    %outside the confinement, inside the confinment, and out again
%     interpolatedSE(1) = mean(compiledSE{i}(sortDisps<-20),'omitnan');
%     interpolatedSE(2) = mean(compiledSE{i}(sortDisps>-20 & sortDisps<20),'omitnan');
%     interpolatedSE(3) = mean(compiledSE{i}(sortDisps>20),'omitnan');
    compiledInterpSE(:,i) = avgBins(compiledSE{i},sortDisps,sampleGrid);
    %compiledInterpSE(:,i) = interpolatedSE;

    %Strain Energy Density
    interpolatedSEDen = interp1(sortDisps,(compiledSEDen{i}(sortIdx)),sampleGrid);
    compiledInterpSEDen(:,i) = interpolatedSEDen;

    %Nuclear Aspect Ratio
    interpolatedNAR = interp1(sortDisps, compiledNucAR{i}(sortIdx),sampleGrid);
    %compiledInterpNAR(:,i) = interpolatedNAR;
    compiledInterpNAR(:,i) = avgBins(compiledNucAR{i},sortDisps,sampleGrid);
    

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
    %compiledInterpCellCirc(:,i) = interpolatedCellCirc;
    compiledInterpCellCirc(:,i) = avgBins(compiledCellCirc{i},sortDisps,sampleGrid);

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

    interpolatedNCM = interp1(sortHullDisps, compiledNCM{i}(sortIdx,1) + compiledNCM{i}(sortIdx,2),sampleGrid);
    compiledInterpNCM(:,i) = interpolatedNCM;
%     interpolatedNCMz = interp1(sortHullDisps, compiledNCM{i}(sortIdx,3),sampleGrid);
%     compiledInterpNCMz(:,i) = interpolatedNCMz;

    %XYZ Force projections
    interpxProj = interp1(sortHullDisps, compiledtractionProj{i}(sortIdx,1),sampleGrid);
    %compiledInterpxProj(:,i) = interpxProj;
    compiledInterpxProj(:,i) = avgBins(compiledtractionProj{i}(:,1),sortDisps,sampleGrid);
    
    interpyProj = interp1(sortHullDisps, compiledtractionProj{i}(sortIdx,2),sampleGrid);
    %compiledInterpyProj(:,i) = interpyProj;
    compiledInterpyProj(:,i) = avgBins(compiledtractionProj{i}(:,2),sortDisps,sampleGrid);

    interpzProj = interp1(sortHullDisps, compiledtractionProj{i}(sortIdx,3),sampleGrid);
    %compiledInterpzProj(:,i) = interpzProj;
    compiledInterpzProj(:,i) = avgBins(compiledtractionProj{i}(:,3),sortDisps,sampleGrid);
    
    %Total contractile sum
    interpxContract = interp1(sortHullDisps, compiledxContractile{i}(sortIdx),sampleGrid);
    %compiledInterpxContractile(:,i) = interpxContract;
    compiledInterpxContractile(:,i) = avgBins(compiledxContractile{i},sortDisps,sampleGrid);

    %Forces within 5 mircons
    interpxContract5micSum = interp1(sortHullDisps, compiledxForces5micSum{i}(sortIdx),sampleGrid);
    compiledInterpforces5micSum(:,i) = interpxContract5micSum;

    interpxContract5micAvg = interp1(sortHullDisps, compiledxForces5micAvg{i}(sortIdx),sampleGrid);
    %compiledInterpforces5micAvg(:,i) = interpxContract5micAvg;
    compiledInterpforces5micAvg(:,i) = avgBins(compiledxForces5micAvg{i},sortDisps,sampleGrid);

    interpxContract5micMax = interp1(sortHullDisps, compiledxForces5micMax{i}(sortIdx),sampleGrid);
    compiledInterpforces5micMax(:,i) = interpxContract5micMax;
    
    %Disps within 5 microns
    interpxDisps5micSum = interp1(sortHullDisps, compiledxDisps5micSum{i}(sortIdx),sampleGrid);
    compiledInterpdisp5micSum(:,i) = interpxDisps5micSum;
    
    interpxDisps5micAvg = interp1(sortHullDisps, compiledxDisps5micAvg{i}(sortIdx),sampleGrid);
    %compiledInterpdisp5micAvg(:,i) = interpxDisps5micAvg;
    compiledInterpdisp5micAvg(:,i) = avgBins(compiledxDisps5micAvg{i},sortDisps,sampleGrid);

    interpxDisps5micMax = interp1(sortHullDisps, compiledxDisps5micMax{i}(sortIdx),sampleGrid);
    compiledInterpdisp5micMax(:,i) = interpxDisps5micMax;

    compiledAvgMaxForce(:,i) = avgBins(compiledMaxForce{i},sortDisps,sampleGrid);

end

xline(-20)
xline(20)
ylabel('Relative Strain Energy');
xlabel('Nuclear Distance to Center of Confinement (um)')
%xlabel('Time Relative to Center of Channel (mins)')
legend('F00','F01','F02','F05','F12','F06','F08','F09','F11','F15')

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
swarmchart(((allDisps+8)),allSE,10,'r','filled')
ylim([0, 5*10^5])
%xlim([2.5 13.5])
xlabel('Nuclear Distance to Center of Confinement (um)')
ylabel('Strain Energy (pJ)')
 


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
ylabel('Strain Energy (pJ)')
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
avgNCM = mean(compiledInterpNCM,2,'omitnan');
errNCM = std(compiledInterpNCM,[],2,'omitnan')/sqrt(size(compiledInterpNCM(1,:),2));
errorbar(sampleGrid,avgNCM,errNCM,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgNCM,'r', 'LineWidth',1)


% errorbar(sampleGrid,avgNCMz,errNCMz,'LineWidth',0.5)
% hold on 
% plot(sampleGrid, avgNCMz, 'LineWidth',1)

ylabel('Average NCM (pNm)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

%XYZ projections
figure
avgxProj = mean(compiledInterpxProj,2,'omitnan');
errxProj = std(compiledInterpxProj,[],2,'omitnan')/sqrt(size(compiledInterpxProj(1,:),2));
errorbar(sampleGrid,avgxProj,errxProj,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxProj,'r', 'LineWidth',1)
ylabel('Sum xTractions (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

figure
avgyProj = mean(compiledInterpyProj,2,'omitnan');
erryProj = std(compiledInterpyProj,[],2,'omitnan')/sqrt(size(compiledInterpyProj(1,:),2));
errorbar(sampleGrid,avgyProj,erryProj,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgyProj,'r', 'LineWidth',1)
ylabel('Sum yTractions (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

figure
avgzProj = mean(compiledInterpzProj,2,'omitnan');
errzProj = std(compiledInterpzProj,[],2,'omitnan')/sqrt(size(compiledInterpzProj(1,:),2));
errorbar(sampleGrid,avgzProj,errzProj,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgzProj,'r', 'LineWidth',1)
ylabel('Sum zTractions (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20)

% Plot xContractile Sums
figure
avgxCont = mean(compiledInterpxContractile  ,2,'omitnan');
errxCont = std(compiledInterpxContractile,[],2,'omitnan')/sqrt(size(compiledInterpxContractile(1,:),2));
errorbar(sampleGrid,avgxCont,errxCont,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxCont,'r', 'LineWidth',1)
ylabel('X Contractile Tractions (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

% Plot contractile forces within 5 mic nucleus
figure
avgxCont5micSum = mean(compiledInterpforces5micSum  ,2,'omitnan');
errxCont5micSum = std(compiledInterpforces5micSum,[],2,'omitnan')/sqrt(size(compiledInterpforces5micSum(1,:),2));
errorbar(sampleGrid,avgxCont5micSum,errxCont5micSum,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxCont5micSum,'r', 'LineWidth',1)
ylabel('Sum X Contractile Tractions 5micNuc (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

figure
avgxCont5micAvg = mean(compiledInterpforces5micAvg  ,2,'omitnan');
errxCont5micAvg = std(compiledInterpforces5micAvg,[],2,'omitnan')/sqrt(size(compiledInterpforces5micAvg(1,:),2));
errorbar(sampleGrid,avgxCont5micAvg,errxCont5micAvg,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxCont5micAvg,'r', 'LineWidth',1)
ylabel('Average X Contractile Tractions 5micNuc (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

figure
avgxCont5micMax = mean(compiledInterpforces5micMax  ,2,'omitnan');
errxCont5micMax = std(compiledInterpforces5micMax,[],2,'omitnan')/sqrt(size(compiledInterpforces5micMax(1,:),2));
errorbar(sampleGrid,avgxCont5micMax,errxCont5micMax,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxCont5micMax,'r', 'LineWidth',1)
ylabel('Max X Contractile Tractions 5micNuc (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

% Disps
figure
avgxdisp5micSum = mean(compiledInterpdisp5micSum,2,'omitnan');
errxdisp5micSum = std(compiledInterpdisp5micSum,[],2,'omitnan')/sqrt(size(compiledInterpdisp5micSum(1,:),2));
errorbar(sampleGrid,avgxdisp5micSum,errxdisp5micSum,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxdisp5micSum,'r', 'LineWidth',1)
ylabel('Sum x Disps (um)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

figure
avgxdisp5micAvg = mean(compiledInterpdisp5micAvg,2,'omitnan');
errxdisp5micAvg = std(compiledInterpdisp5micAvg,[],2,'omitnan')/sqrt(size(compiledInterpdisp5micAvg(1,:),2));
errorbar(sampleGrid,avgxdisp5micAvg,errxdisp5micAvg,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxdisp5micAvg,'r', 'LineWidth',1)
ylabel('Avg x Disps (um)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

figure
avgxdisp5micMax = mean(compiledInterpdisp5micMax,2,'omitnan');
errxdisp5micMax = std(compiledInterpdisp5micMax,[],2,'omitnan')/sqrt(size(compiledInterpdisp5micMax(1,:),2));
errorbar(sampleGrid,avgxdisp5micMax,errxdisp5micMax,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgxdisp5micMax,'r', 'LineWidth',1)
ylabel('Max x Disps (um)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

%Plot Max Force
figure
avgMaxForce = mean(compiledAvgMaxForce,2,'omitnan');
errMaxForce = std(compiledAvgMaxForce,[],2,'omitnan')/sqrt(size(compiledAvgMaxForce(1,:),2));
errorbar(sampleGrid,avgMaxForce,errMaxForce,'b','LineWidth',0.5)
hold on 
plot(sampleGrid, avgMaxForce,'r', 'LineWidth',1)
ylabel('Max Force (Pa)')
xlabel('Minimum Nuclear Distance to Center of Confinement (um)')
xline(-20)
xline(20) 

%% Save Resulting Measurements
% save('T:\Max\2023-12-20\Tiffs\IA32LAConfinementValues.mat',"sampleGrid","compiledInterpSE","compiledInterpSEDen","compiledInterpRelSE","compiledInterpNAR", "compiledHullInterpSE", "compiledInterpNARHull", "compiledInterpOriHull",...
%     "compiledInterpOri", "compiledTransitTimes", "compiledAvgStrainTransit", "compiledMaxStrainTransit", "compiledInterpNucCellDist",...
%     "compiledInterpCellCirc", "compiledInterpCellArea", "compiledNucVel" ,"compiledInterpCellVel")
%% Cross quantifications

% %Average NAR vs Average Relative Strain Energy
% figure
% scatter(avgNAR,avgRelSE)

%%     Compute distributions of force

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