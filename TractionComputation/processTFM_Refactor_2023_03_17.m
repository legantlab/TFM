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
%Apply any registration correction
for i = 1:length(matches)
    curMatches = matches{i};
    %curMatches(:,1:2) = curMatches(:,1:2)*0.1844;
    %curMatches(:,4:5) = curMatches(:,4:5)*0.1844;
    %curMatches(:,3) = curMatches(:,3)*0.8;
    %curMatches(:,6) = curMatches(:,6)*0.8;

    curMatches(:,1) = curMatches(:,1) - 3;
    curMatches(:,3) = curMatches(:,3) + 2;
    curMatches(:,4) = curMatches(:,4) - 3;
    %curMatches(:,3) = curMatches(:,3) + 1.5;
    curMatches(:,6) = curMatches(:,6) + 2;
    matches{i} = curMatches;
    curDisplacement = curMatches(:,4:6) - curMatches(:,1:3);
    displacements{i} = curDisplacement;
end

%Plotting
figure
% quiver3(matches{1}(:,1), matches{1}(:,2),matches{1}(:,3), displacements{1}(:,1),...
%     displacements{1}(:,2),displacements{1}(:,3),3)
scatter3(matches{1}(:,1), matches{1}(:,2),matches{1}(:,3))
title('Corrected Displacements')
hold on 
scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))
legend('Bead Localizations', 'Element Centroids')
%% Interpolate matches to surface Nodes
beadDisps = interpDisps(matches,elemCents2D);

%% Optional step to strip out SDS derived substrate expansion
%matches = polyfitThermalExpansion(matches,displacements,11, 5, 5, 1);
newDisps = interpThermalExpansion(beadDisps, 0);
beadDisps = newDisps;
%% Read nodal solutions and parse
nodalSolutionsFile = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\SwellingControl\SwellingControlExpandedNodeSetCorrected.h5';
nodeDisps = readNodalSolutions(nodalSolutionsFile,length(data.ElementSets.Data),newIds);
%% Construct Green's matrix
[beadu] = composeGreens(nodeDisps, elemCents2D, TR2);

%% Common checks on beadu
%Reverse z component of last thermal subcase?
%beadu(3:3:end, end) = beadu(3:3:end, end)*-1;
%beadu(1:3:end, end) = beadu(1:3:end, end) * -1;
figure
subcase = 15550; 
quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), beadu(1:3:end, subcase),...
    beadu(2:3:end, subcase), beadu(3:3:end,subcase),1);
hold on 
quiver3(matches{end}(:,1), matches{end}(:,2),matches{end}(:,3),...
    displacements{end}(:,1), displacements{end}(:,2),displacements{end}(:,3),1)
legend('Simulated Swelling', 'Computed Displacements')
%% Compute CSVD and Tractions

%Matrix decomposition 
[U,s,V]=csvd(beadu);
%%
%Solve for first time point
temp=beadDisps{1}(:,4:6)';
%temp(3,:) = temp(3,:)*-1;
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U,s,y,'Tikh');

%Compute forces at different time points
reg_corner_timelapse=10; %L-curve picks out the wrong corner right now

%% Plotting
force_vector={};
for kk=1:length(beadDisps)
    kk
    temp=beadDisps{kk}(:,4:6)';
    %temp(3,:) = temp(3,:)*-1; %Account for swelling reversal
    y=temp(:);
    reg_corner_timelapse=10; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U,s,V,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[elemCents2D,x_lambda(1:3:end-1),...
        x_lambda(2:3:end-1),x_lambda(3:3:end-1)];
end

%%
%Save force_vector and elemCents2D for later plotting/to return to this
%point...
date = string(datetime('today','Format','yyyy-MM-dd'));
plotTime = 1;
figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{plotTime}(:,4),force_vector{plotTime}(:,5),0)
hold on 
scatter(elemCents2D(:,1), elemCents2D(:,2))
quiver(matches{plotTime}(:,1),matches{plotTime}(:,2),displacements{plotTime}(:,1),displacements{plotTime}(:,2),1)
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{plotTime}(:,4),...
    force_vector{plotTime}(:,5),force_vector{plotTime}(:,6),3)
axis equal
%save(strcat('U:\Max\2023_01_18_CellsInChannels\S2\', date,'_ComputedTractionsVarsNice.mat'),'force_vector','elemCents2D','U','s','V','beadDispsCleaned')
%% Plot swelling
recon_swell = beadu(:,end)*x_lambda(end)*5;
figure
subplot(1,2,1)
quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), ...
    recon_swell(1:3:end), recon_swell(2:3:end), recon_swell(3:3:end), 0)
axis equal
hold on
subplot(1,2,2)
quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), beadu(1:3:end,end),...
    beadu(2:3:end,end), beadu(3:3:end,end),0)
%quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), ...
%    beadDisps{1}(:,4),beadDisps{1}(:,5),beadDisps{1}(:,6),1 )
axis equal
%legend('Calculated Swelling Displacement','Simulated Swelling Displacement')

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
    subplot(1,2,1)
    quiver(elemCents2D(:,1), elemCents2D(:,2), curForce(:,4), curForce(:,5),2)
    title('Computed Traction Force')
    hold on
    subplot(1,2,2)
    quiver(beadDisps{i}(:,1), beadDisps{i}(:,2),beadDisps{i}(:,4),beadDisps{i}(:,5),1)
    axis equal
    %legend('Computed Traction Force', 'Interpolated Displacements')
    title('Interpolated Displacements')
    %axis off
    saveas(gcf, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\Forces_Interp\TractionMaps_' ...
        , sprintf( '%02d', i ), '.png'])
    close
end


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