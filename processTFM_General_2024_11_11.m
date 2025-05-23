%% Run file for Processing Tractions from microchannel TFM data
% Takes input displacement field and localization data from
% main_tracking_script.m and computes tractions using inputted geometry
% specific Green's function. Results are saved to an output folder along
% with plots of the reconstructions compared to input displacements,
% fitting of interpolated displacements to geometry, projections of the
% tractions, and Amira mesh files for the displacements and tractions for
% 3D visualization.
%
%  INPUTS
% -------------------------------------------------------------------------
%   outputFold: folder path to save output files and plots
%   voxelSize: array containing xyz voxel dimensions in microns.
%   fileName: path to the solver deck and geometry files outputted from
%   hypermesh.
%   disps: path to the displacement data computed from
%   main_tracking_script.m.
%   Masks: path to the masks of the cell used to mask the displacement
%   field.
%   nodalSolutionsFile: path to the solved FEA displacements produced by
%   optistruct which are used to construct a geometry specific Green's
%   function.
%
%  OUTPUTS
%  ------------------------------------------------------------------------
%   traction_vector: cell array containing a double array for each time
%   point with computed traction vectors in the form of an nx6 matrix where
%   each row is the xyz coordinate position of the centroid of the model
%   geometry elements and then the xyz components of the traction vector.
%   The resulting traction vector components can be multiplied by the
%   Young's modulus (in Pa) for your specific substrate stiffness. 
%
%   Several plots, visualization figures and output files for 3D
%   visualization in Amira. See end of script for details. 
%
%   Author: Max Hockenberry
%   Last Update: 11/11/2024

outputFold = '.\ExampleDataSet'; %Folder to store output tractions and amira maps.
plotFlag = 1; %If flag = 1, outputs a series of QA plots to the output folder.
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots'])
end

voxelSize = [0.199, 0.199, 0.8]; % VoxelSize XYZ in microns
%% Open solverdeck using abaqus reader
fileName = '.\ExampleDataSet\5x40SolverDeck.inp';
data = abaqusInpRead(fileName);
%% Reorganize geometry and inputs
[newIds, TR2,newNodeCoords, surfNodes,elemCents2Dnew] = organizeGeometry(data,1,1);
%% Get the area of each element for computing force balance
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
%% Load displacements and mask
disps = ".\ExampleDataSet\ExampleDataSetDisplacements.mat";
load(disps)

%Write a function to do all the below.
for i = 1:length(matches) % Fix pixel scaling.
    curMatch = matches{i};
    curMatch(:,1:2) = curMatch(:,1:2)*voxelSize(1);
    curMatch(:,3) = curMatch(:,3)*voxelSize(3);
    curMatch(:,4:5) = curMatch(:,4:5)*voxelSize(1);
    curMatch(:,6) = curMatch(:,6)*voxelSize(3);
    matches{i} = curMatch; %#ok<SAGROW>
end

%only consider elements within 40 microns of the center of the confinement
regElems = elemCents2Dnew(elemCents2Dnew(:,1)>-40 & elemCents2Dnew(:,1)<40,:); 

%Align elemCents2Dnew to matches
PC1 = pointCloud(regElems);
PC2 = pointCloud(matches{1}(:,1:3));

tform = pcregistericp(PC1,PC2);
PC1reg = pctransform(pointCloud(elemCents2Dnew), tform);
elemCents2Dnew = PC1reg.Location;

%Mask the raw displacements by cell mask
Masks = ".\ExampleDataSet\Cell_Mask";
curDir = char(Masks);

tiffs=dir([curDir,'\*.tif']);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
currentAllMasks = cell(length(tiffs),1);
for i = 1:length(tiffs(:,1))
    curMask = [tiffs(i).folder,'\',tiffs(i).name];
    currentAllMasks{i} = loadtiff(curMask);
end

dilation = strel('diamond',100); %Dilate the mask by 20 microns
maskedCoords = matches;
for j = 1:length(currentAllMasks)
    curMask = rot90(currentAllMasks{j},3);
    %curMask = (currentAllMasks{j});
    ImgSize = size(curMask);
    curMask = imdilate(curMask, dilation);
    matchCoords = matches{j}(:,:);
    matchCoords(:,1) = round(matchCoords(:,1)/voxelSize(2) + ImgSize(2)/2 - 1);
    matchCoords(:,2) = round(matchCoords(:,2)/voxelSize(1) + ImgSize(1)/2 - 1);

    [ycoords, xcoords] = find(curMask);
    coordinateMask = ismember(matchCoords(:,1:2),[xcoords,ycoords],'rows');
    maskedCoords{j} = matches{j}(coordinateMask,:);
    %rather set the second set of matches of the not coordinate mask to be
    %equal to the first set of corrdinates
    matches{j}(~coordinateMask,4:6) = matches{j}(~coordinateMask,1:3); %#ok<SAGROW>
    maskedCoords{j} = matches{j}(coordinateMask,:);
end

time = 1;
newDisps = matches{time}(:,4:6) - matches{time}(:,1:3);

if plotFlag == 1
    %Plot elements and matches to ensure they align
    figure
    scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))
    hold on
    scatter3(matches{time}(:,4), matches{time}(:,5), matches{time}(:,6))
    axis equal
    legend('Elements','Data')
    saveas(gcf,[outputFold, '\', 'QAplots\ModelAlignment.fig'])
    close()

    %Plot Masked Coordinates
    figure
    scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))
    hold on
    scatter3(matches{time}(:,4), matches{time}(:,5), matches{time}(:,6))
    scatter3(maskedCoords{time}(:,1), maskedCoords{time}(:,2),maskedCoords{time}(:,3),'b')
    axis equal
    legend('Elements','Data', 'Masked Data')
    saveas(gcf,[outputFold, '\', 'QAplots\MaskedCoordinates.fig'])
    close()
end
%% Interpolate matches to surface Nodes
beadDisps = interpDispsDisplacements(matches,elemCents2Dnew);

%Plot inteperpolated displacements versus raw
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots\InterpDisps'])
    for i = 1:length(beadDisps)
        curBeadDisps = beadDisps{i};
        curMatch = matches{i};
        curDisps = displacements{i};
        figure
        quiver(curBeadDisps(:,1), curBeadDisps(:,2), ...
            curBeadDisps(:,4),curBeadDisps(:,5))
        hold on
        quiver(curMatch(:,1), curMatch(:,2), ...
            curDisps(:,1),curDisps(:,2))
        legend('Interpolated Displacements','Raw Displacements')
        axis equal
        set(gcf,'units','normalized','position',[0 0 1 1])
        saveas(gcf,[outputFold, '\', 'QAplots\InterpDisps\t',sprintf( '%04d', i ),'_MaskedCoordinates.png'])
        close()
    end
end
%% Read nodal solutions and parse
%nodalSolutionsFile = ".\5x40_Raw_SmoothMesh2p5_s10Course_QA_Tetra_Properties_Loads_SPCs_Loaded.h5";
%nodeDisps = readNodalSolutions(nodalSolutionsFile,length(data.ElementSets.Data),newIds);

%As the solved FEA geometry is large (~several GB), the resulting nodal
%displacements have been provided as a .mat file. The code to read the
%nodal solutions file output from optistruct has been provided for
%demonstration. - MH 11/11/2024
load('.\ExampleDataSet\nodalDisplacements1.mat') 
load('.\ExampleDataSet\nodalDisplacements2.mat') 
nodeDisps = [nodes1;nodes2];
%% Construct Green's matrix
[beadu] = composeGreens(nodeDisps,elemCents2Dnew, TR2);
%% Compute CSVD and Tractions
%Matrix decomposition
[U2,s2,V2]=csvd(beadu);

%% L Curve
%Solve for first time point
temp=beadDisps{1}(:,4:6)';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U2,s2,y,'Tikh');
if plotFlag == 1
    saveas(gcf,[outputFold, '\', 'QAplots\Lcurve.png'])
    close()
end
%% Plotting
close all
traction_vector=cell(length(beadDisps),1);
reg_corner_timelapse = 10;%reg_corner;
newElemPositions=cell(length(beadDisps),1);
for kk =1:length(beadDisps)
    temp=beadDisps{kk}(:,4:6)';
    y=temp(:);
    %reg_corner_timelapse=reg_corner; %Uncomment to use L-curve determined
    %regularization value.
    [x_lambda,rho,eta] = tikhonov2(U2,s2,V2,y,reg_corner_timelapse);

    %Compute updated positions for displacements
    computedDisp = beadu*x_lambda;
    computedDispArray = [computedDisp(1:3:end), computedDisp(2:3:end), computedDisp(3:3:end)];
    newElemPositions{kk} = elemCents2Dnew + computedDispArray;
    traction_vector{kk}=[newElemPositions{kk},x_lambda(1:3:end),...
        x_lambda(2:3:end),x_lambda(3:3:end)];
end

date = string(datetime('today','Format','yyyy-MM-dd'));
scale = 20000; %Young's Modulus.
figure
kk = 1; %time to plot
quiver(traction_vector{kk}(:,1),traction_vector{kk}(:,2),traction_vector{kk}(:,4)*scale,traction_vector{kk}(:,5)*scale,1)
axis equal
figure
quiver3(traction_vector{kk}(:,1),traction_vector{kk}(:,2),traction_vector{kk}(:,3),traction_vector{kk}(:,4)...
    *scale,traction_vector{kk}(:,5)*scale,traction_vector{kk}(:,6)*scale,1)
axis equal

save(strcat(outputFold, '\', date,'_ComputedTractionsAlignedValues.mat'),'traction_vector','elemCents2Dnew','beadDisps')
%% Save 2D QA Plots of tractions
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots\Tractions'])
    for i = 1:length(traction_vector)
        figure
        quiver(traction_vector{i}(:,1),traction_vector{i}(:,2),traction_vector{i}(:,4),traction_vector{i}(:,5),1)
        axis equal
        set(gcf,'units','normalized','position',[0 0 1 1])
        saveas(gcf,[outputFold, '\', 'QAplots\Tractions\t',sprintf( '%04d', i ),'_tractions.png'])
        close()
    end
end

%% Output vectors in Amira formatting
mkdir([outputFold, '\',convertStringsToChars(date),'_Tractions'])
for i = 1:length(traction_vector)
    disp(i)
    translatePos = traction_vector{i}(:,1:3);
    curTranslate = translation{i};
    translatePos(:,1) = (translatePos(:,1)) + (curTranslate(1))*voxelSize(1);
    translatePos(:,2) = (translatePos(:,2)) + (curTranslate(2))*voxelSize(2);
    translatePos(:,3) = (translatePos(:,3)) + curTranslate(3)*voxelSize(3);
    text=buildAmiraMesh(1*traction_vector{i}(:,4:6),translatePos);
    fid=fopen([outputFold, '\',convertStringsToChars(date),'_Tractions', '\t', sprintf( '%04d', i ), 'Forces.am'],'wt');
    fprintf(fid,text);
    fclose(fid);
end

%% Output ScalarMagnitude of vectors in Amira formatting
mkdir([outputFold, '\',convertStringsToChars(date),'_TractionsContour'])
for i = 1:length(traction_vector)
    translatePos = traction_vector{i}(:,1:3);
    curTranslate = translation{i};
    translatePos(:,1) = (translatePos(:,1)) + (curTranslate(1))*voxelSize(1);
    translatePos(:,2) = (translatePos(:,2)) + (curTranslate(2))*voxelSize(2);
    translatePos(:,3) = (translatePos(:,3)) + curTranslate(3)*voxelSize(3);
    text=buildAmiraMeshMagnitude(1*traction_vector{i}(:,4:6),translatePos);
    fid=fopen([outputFold, '\',convertStringsToChars(date),'_TractionsContour', '\t', sprintf( '%04d', i ), 'Forces.am'],'wt');
    fprintf(fid,text);
    fclose(fid);
end
%% Output Displacement vectors in Amira formatting
load(disps)
for i = 1:length(matches) % Fix pixel scaling.
    curMatch = matches{i};
    curMatch(:,1:2) = curMatch(:,1:2)*voxelSize(1);
    curMatch(:,3) = curMatch(:,3)*voxelSize(3);
    curMatch(:,4:5) = curMatch(:,4:5)*voxelSize(1);
    curMatch(:,6) = curMatch(:,6)*voxelSize(3);
    matches{i} = curMatch; %#ok<SAGROW>
end

mkdir([outputFold, '\',convertStringsToChars(date),'_Disps'])
for i = 1:length(matches)
    disp(i)
    translatePos = matches{i}(:,1:3);
    curDisp = matches{i}(:,4:6) - matches{i}(:,1:3);
    curTranslate = translation{i};
    translatePos(:,1) = (translatePos(:,1)) + (curTranslate(1))*voxelSize(1);
    translatePos(:,2) = (translatePos(:,2)) + (curTranslate(2))*voxelSize(1);
    translatePos(:,3) = (translatePos(:,3)) + curTranslate(3)*voxelSize(3);
    text=buildAmiraMesh(1*curDisp,translatePos);
    fid=fopen([outputFold, '\',convertStringsToChars(date),'_Disps', '\t', sprintf( '%04d', i ), 'Disps.am'],'wt');
    fprintf(fid,text);
    fclose(fid);
end

%% Output the mesh with translation and Displacement from cell
mkdir([outputFold,'\AlignedStls'])
MeshPoints = TR2.Points;

for i = 1:length(matches)
    trans = translation{i};
    curMeshPoints = MeshPoints;
    curMeshPoints(:,1) = curMeshPoints(:,1) + trans(1)*voxelSize(1);
    curMeshPoints(:,2) = curMeshPoints(:,2) + trans(2)*voxelSize(1);
    curMeshPoints(:,3) = curMeshPoints(:,3) + trans(3)*voxelSize(3);

    TR3 = triangulation(TR2.ConnectivityList,curMeshPoints);
    stlwrite([outputFold,'\AlignedStls', '\t', sprintf( '%04d', i ),'AlignedSTL.stl'],struct('faces',TR3.ConnectivityList,'vertices',TR3.Points))

end