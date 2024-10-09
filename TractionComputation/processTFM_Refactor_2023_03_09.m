%Refactored TFM Script

close all
clear
clc

%% Open solverdeck using abaqus reader
fileName = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrected.inp';
data = abaqusInpRead(fileName);

%% Reorganize geometry and inputs
[newIds, TR2,newNodeCoords, surfNodes] = organizeGeometry(data,1,1);

%% Load displacements
load("U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlCorrectedUnits.mat")

%% Optional step to strip out SDS derived substrate expansion
%matches = polyfitThermalExpansion(matches,displacements,11, 5, 5, 1);

%% Interpolate matches to surface Nodes
beadDisps = interpDisps(matches,newNodeCoords(surfNodes,2:4));

%% Filter node lists by those we report forces from
appliedNodeList = zeros(length(data.NodeSets.Data),4);
appliedNodeList(:,1) = data.NodeSets.Data;
for i = 1:length(appliedNodeList)
    curNode = appliedNodeList(i,1);
    nodeID = find(curNode == newIds(:,1));
    appliedNodeList(i,1) = nodeID;
    appliedNodeList(i,2:4) = newNodeCoords(nodeID,2:4);
end
%% Interpolate Bead displacements to reference frame
%beadDisps = interpDisps(matches,appliedNodeList);

%% Identify mesh elements for each observed bead displacement
[newElem3DConn,beadDispsCleaned, beadElems] = organizeGeometry3D(data,newIds,newNodeCoords, matches,TR2);

%% Read nodal solutions and parse
nodalSolutionsFile = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrectedRetest.h5';
nodeDisps = readNodalSolutions(nodalSolutionsFile,elemCents2D,newIds);

%% Compute CSVD and Tractions
%Reorganize into the structure we want
[beadu]=beadDispCalcLinear(beadDispsCleaned{1}(:,1:3),beadElems,appliedNodeList,...
    [1:1:size(appliedNodeList,1)]',nodeDisps,newElem3DConn,newNodeCoords);

%Matrix decomposition 
[U,s,V]=csvd(beadu);

%Solve for first time point
temp=beadDispsCleaned{1}(:,1:3)';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U,s,y,'Tikh');

%Compute forces at different time points
reg_corner_timelapse=0.0001; %L-curve picks out the wrong corner right now

%% Plotting
force_vector={};
for kk=1:length(beadDispsCleaned)
    kk
    temp=beadDispsCleaned{kk}(:,4:6)';
    y=temp(:);
    reg_corner_timelapse=1; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov(U,s,V,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[beadDispsCleaned{1}(:,1:3),x_lambda(1:3:end-1),...
        x_lambda(2:3:end-1),x_lambda(3:3:end-1)];
end


%Save force_vector and elemCents2D for later plotting/to return to this
%point...
date = string(datetime('today','Format','yyyy-MM-dd'));
figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4),force_vector{1}(:,5),1)
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),force_vector{1}(:,6),1)
axis equal
%save(strcat('U:\Max\2023_01_18_CellsInChannels\S2\', date,'_ComputedTractionsVarsNice.mat'),'force_vector','elemCents2D','U','s','V','beadDispsCleaned')

%% More visualization

figure %Vector map on mesh with colored mesh elements by force
%We should probably try to compute pressure/force to account for element
%size
triaArea =TriaElementArea(TR2);
forceVectorSum = vecnorm(force_vector{1}(:,4:6),2,2);
colorMapTest = normalize(forceVectorSum,'range',[0, 1]);
colorMapTest = colorMapTest./triaArea;
patch('Faces',TR2.ConnectivityList,'Vertices',TR2.Points,'FaceColor','flat','CData',colorMapTest)
hold on
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),...
    force_vector{1}(:,6),1,'r')

%% Visualize Displacements and Computed forces on top of each other

figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),...
    force_vector{1}(:,6),1,'r')
hold on 
quiver3(beadDispsCleaned{1}(:,1),beadDispsCleaned{1}(:,2),beadDispsCleaned{1}(:,3),...
    beadDispsCleaned{1}(:,4),beadDispsCleaned{1}(:,5),...
    beadDispsCleaned{1}(:,6),1,'b')
title('Interpolated Displacements and Computed Tractions')
axis equal