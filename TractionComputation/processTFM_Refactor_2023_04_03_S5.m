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
load("U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\S5_Full_New.mat")

%% Manual correction for registration errors
%Plot element centroids as compared to matches to see if correction
%is needed
xCor = -2.5;
yCor = 0;
zCor = 1;

figure
scatter3(matches{1}(:,1) + xCor, matches{1}(:,2) + yCor,matches{1}(:,3) + zCor,'b')
hold on 
scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3),'r')
legend('Bead localizations', 'Element Centroids')

for i = 1:length(matches)
    curMatches = matches{i};
    curMatches(:,1) = curMatches(:,1) + xCor;
    curMatches(:,2) = curMatches(:,2) + yCor;
    curMatches(:,3) = curMatches(:,3) + zCor;
    curMatches(:,4) = curMatches(:,4) + xCor;
    curMatches(:,5) = curMatches(:,5) + yCor;
    curMatches(:,6) = curMatches(:,6) + zCor;
end


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
title('Raw Displacements')
subplot(1,2,2)
quiver3(correctedbeadDisps{1}(:,1), correctedbeadDisps{1}(:,2),correctedbeadDisps{1}(:,3),correctedbeadDisps{1}(:,4),correctedbeadDisps{1}(:,5),correctedbeadDisps{1}(:,6),0)
title('Corrected Displacements')
%% Read nodal solutions and parse
nodalSolutionsFile = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrected.h5';
nodeDisps = readNodalSolutions(nodalSolutionsFile,length(data.ElementSets.Data),newIds);
%% Construct Green's matrix
[beadu] = composeGreens(nodeDisps,elemCents2D, TR2);

%% Replace last column in beadu with experimental derived swelling subcase
newbeadu = beadu;
temp = averagedSwellingInterp';
newbeadu(:,end) = temp(:);

%% Compute CSVD and Tractions
%Matrix decomposition 
%[U,s,V]=csvd(beadu);

[U2,s2,V2]=csvd(newbeadu);
%%
%Solve for first time point
temp=correctedbeadDisps{1}(:,4:6)';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U2,s2,y,'Tikh');

%Compute forces at different time points
%reg_corner_timelapse=0.0001; %L-curve picks out the wrong corner right now

%% Plotting
force_vector={};
reg_corner = 50;
for kk=1:length(correctedbeadDisps)
    kk
    temp=correctedbeadDisps{kk}(:,4:6)';
    y=temp(:);
    reg_corner_timelapse=reg_corner; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U2,s2,V2,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[elemCents2D,x_lambda(1:3:end-1),...
        x_lambda(2:3:end-1),x_lambda(3:3:end-1)];
end


%Save force_vector and elemCents2D for later plotting/to return to this
%point...
date = string(datetime('today','Format','yyyy-MM-dd'));
scale = 200;
figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4)*scale,force_vector{1}(:,5)*scale,0)
axis equal
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4)*scale,force_vector{1}(:,5)*scale,force_vector{1}(:,6)*scale,0)
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

%% Plot Timeframes
%Plot Raw displacements
for i = 1:length(beadDisps)
    curMatches = matches{i};
    figure
    quiver(matches{i}(:,1), matches{i}(:,2), displacements{i}(:,1), displacements{i}(:,2))
    axis equal
    axis off
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\Displacements\' ...
        , sprintf( '%02d', i ), '.png'])
    close
end

%Plot interpolated displacements
for i = 1:length(beadDisps)
    figure
    subplot(1,2,1)
    quiver3(beadDisps{i}(:,1), beadDisps{i}(:,2), beadDisps{i}(:,3), beadDisps{i}(:,4), ...
        beadDisps{i}(:,5),beadDisps{i}(:,6),1)
    axis equal
    title('Interpolated Displacements')
    hold on
    subplot(1,2,2)
    quiver3(matches{i}(:,1), matches{i}(:,2), matches{i}(:,3), displacements{i}(:,1), ...
        displacements{i}(:,2),displacements{i}(:,3),1)
    title('Raw Displacements')
    axis equal
    %axis off
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\InterpolatedDisps\' ...
        , sprintf( '%02d', i ), '.png'])
    close
end

%Plot Comparison between displacements and interpolated displacements and
%the centroids of the elements
for i = 1:length(beadDisps)
    figure
    quiver3(beadDisps{i}(:,1), beadDisps{i}(:,2), beadDisps{i}(:,3),beadDisps{i}(:,4), beadDisps{i}(:,5),beadDisps{i}(:,6),0,'r')
    hold on 
    quiver3(matches{i}(:,1), matches{i}(:,2), matches{i}(:,3), displacements{i}(:,1), displacements{i}(:,2),displacements{i}(:,3),'b')
    scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3))
    axis equal
    axis off
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\InterpolatedDisps\' ...
        , sprintf( '%02d', i ), '.png'])
    close
end

%Plot computed tractions vs displacements in 2D graph for timelapse
for i = 1:length(force_vector)
    curForce = force_vector{i};
    figure
    scale = 200;
    subplot(1,2,1)
    quiver(elemCents2D(:,1), elemCents2D(:,2), curForce(:,4)*scale, curForce(:,5)*scale,0)
    title('Computed Traction Force')
    axis tight manual
    hold on
    subplot(1,2,2)
    quiver(correctedbeadDisps{i}(:,1), correctedbeadDisps{i}(:,2),correctedbeadDisps{i}(:,4),correctedbeadDisps{i}(:,5),0)
    axis equal
    %legend('Computed Traction Force', 'Interpolated Displacements')
    title('Interpolated Displacements')
    %axis off
    axis tight manual
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\correcteForcesVec\TractionMaps_' ...
        , sprintf( '%02d', i ), '.svg',])
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\correctedForces\TractionMaps_' ...
        , sprintf( '%02d', i ), '.png',])
    close
end

%% Visualize Displacements and Computed forces on top of each other
% A better visualization would be computed displacements without the
% thermal drift compared to the traction forces
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),...
    force_vector{1}(:,6),1,'r')
hold on 
quiver3(beadDisps{1}(:,1),beadDisps{1}(:,2),beadDisps{1}(:,3),...
    beadu(1:3:end,end),beadu(2:3:end,end),...
    beadu(3:3:end,end),1,'b')
title('Interpolated Displacements and Computed Tractions')
legend('Computed Tractions','Displacements')
axis equal

%% Compute predicted displacments from x_lambda and assess compared to interpolated displacements
test_disps = newbeadu*x_lambda;
test_dispsX = test_disps(1:3:end);
test_dispsY = test_disps(2:3:end);
test_dispsZ = test_disps(3:3:end);

%Plot
figure
quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),test_dispsX, ...
    test_dispsY, test_dispsZ,0,'r')
hold on 
quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),correctedbeadDisps{1}(:,4), ...
    correctedbeadDisps{1}(:,5), correctedbeadDisps{1}(:,6),0,'b')
legend('Predicted Displacments','Interpolated Displacements')
axis equal tight

%Plot residuals
residuals = [correctedbeadDisps{1}(:,4) - test_dispsX, correctedbeadDisps{1}(:,5) - test_dispsY, ...
    correctedbeadDisps{1}(:,6) - test_dispsZ];
numBins = 20;
figure
subplot(1,3,1)
hist(vecnorm(residuals,1,2),numBins)
title('Residual error between computed and observed displacements')
allDisps = beadDisps{1}(:,4:6);
ylim([0,1800])
xlabel('RE (microns)')
ylabel('Count')

subplot(1,3,2)
hist(vecnorm(allDisps,1,2),numBins)
title('Histogram of Raw Displacements')
xlabel('RE (microns)')
ylabel('Count')
ylim([0,1800])

subplot(1,3,3)
corallDisps = correctedbeadDisps{1}(:,4:6);
hist(vecnorm(corallDisps,1,2),numBins)
title('Histogram of Corrected Displacements')
ylim([0,1800])
xlabel('RE (microns)')
ylabel('Count')
%Plot residual errors onto elements
figure
quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),residuals(:,1), ...
    residuals(:,2), residuals(:,3),0,'r')
title('Residual error between computed and observed displacements')


%It would also be good to look at displacments as a function of distance
%from cell, but we need a cell mask for that. 