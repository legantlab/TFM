%Refactored TFM Script

close all
clear
clc

%% Open solverdeck using abaqus reader
fileName = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrected.inp';
data = abaqusInpRead(fileName);

%% Reorganize geometry and inputs
[newIds, TR2,newNodeCoords, surfNodes,elemCents2D] = organizeGeometry(data,1,1);

%% Load displacements
load("U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlCorrectedUnits.mat")

%% Optional step to strip out SDS derived substrate expansion
%matches = polyfitThermalExpansion(matches,displacements,11, 5, 5, 1);
%% Interpolate matches to surface Nodes
beadDisps = interpDisps(matches,elemCents2D);

%% Load swellling control and subtract from computed displacements
load('U:\Max\2023_03_08_ChannelSwelling\AveragedSwellingDisplacementDataSet.mat')
scale = 2.5;
for i = 1:length(beadDisps)
    curBeadDisp = beadDisps{i};
    curBeadDisp(:,4:6) = curBeadDisp(:,4:6) - (scale*averagedSwellingInterp);
    correctedbeadDisps{i} = curBeadDisp;
end

%Plot
figure
subplot(1,2,1)
quiver3(beadDisps{1}(:,1), beadDisps{1}(:,2),beadDisps{1}(:,3),beadDisps{1}(:,4),beadDisps{1}(:,5),beadDisps{1}(:,6),0)
subplot(1,2,2)
quiver3(correctedbeadDisps{1}(:,1), correctedbeadDisps{1}(:,2),correctedbeadDisps{1}(:,3),correctedbeadDisps{1}(:,4),correctedbeadDisps{1}(:,5),correctedbeadDisps{1}(:,6),0)
%% Read nodal solutions and parse
nodalSolutionsFile = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrected.h5';
nodeDisps = readNodalSolutions(nodalSolutionsFile,length(data.ElementSets.Data),newIds);
%% Construct Green's matrix
[beadu] = composeGreens(nodeDisps, beadDisps,elemCents2D, TR2);

%% Compute CSVD and Tractions

%Matrix decomposition 
[U,s,V]=csvd(beadu);
%%
%Solve for first time point
temp=beadDisps{1}(:,4:6)';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U,s,y,'Tikh');

%Compute forces at different time points
%reg_corner_timelapse=0.0001; %L-curve picks out the wrong corner right now

%% Plotting
force_vector={};
for kk=1:length(beadDisps)
    kk
    temp=beadDisps{kk}(:,4:6)';
    y=temp(:);
    reg_corner_timelapse=reg_corner; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U,s,V,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[elemCents2D,x_lambda(1:3:end-1),...
        x_lambda(2:3:end-1),x_lambda(3:3:end-1)];
end


%Save force_vector and elemCents2D for later plotting/to return to this
%point...
date = string(datetime('today','Format','yyyy-MM-dd'));
figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4),force_vector{1}(:,5),0)
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),force_vector{1}(:,6),0)
axis equal
%save(strcat('U:\Max\2023_01_18_CellsInChannels\S2\', date,'_ComputedTractionsVarsNice.mat'),'force_vector','elemCents2D','U','s','V','beadDispsCleaned')
%% Plot swelling
recon_swell = beadu(:,end)*x_lambda(end);
quiver(elemCents2D(:,1), elemCents2D(:,2), recon_swell(1:3:end), recon_swell(2:3:end),1)
hold on
quiver(elemCents2D(:,1), elemCents2D(:,2), beadDisps{1}(:,4),beadDisps{1}(:,5),1 )
axis equal
legend('Computed Swelling Displacement','Observed Displacement')
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