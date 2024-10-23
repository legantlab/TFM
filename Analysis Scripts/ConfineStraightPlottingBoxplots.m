%ConfineVsStraightPlotting Scripts
close all
clear
clc
%straightFile = 'T:\Max\2023-12-20\Tiffs\IA32LAStraightValues.mat';
%straightFile = 'straightIA32CompiledResultsBoxChart3.mat';
straightFile = 'C:\Users\LegantLab\Documents\decode_master_repo\TFM\Datasets\2024_09_19_StraightData_10micAvgBins_80mic.mat';
%confineFile = 'T:\Max\2023-12-20\Tiffs\IA32LAConfinementValues.mat';
%confineFile = 'IA32ConfineValues.mat';
%confineFile = 'confineIA32CompiledResultsBarCharts.mat';
confineFile = 'C:\Users\LegantLab\Documents\decode_master_repo\TFM\Datasets\2024_09_19_ConfineData_10micAvgBins_80mic.mat';
%% Plot Strain Energy
figure(1)
hold on 
load(confineFile,'compiledInterpSE','compiledDistances','sampleGrid')
confineSE = compiledSE;
cellSE = cell(length(compiledSE),1);
for i = 1:length(cellSE)
    curSE = compiledSE{i};
    curSE(curSE==0) = nan;
    cellSE{i} = curSE;
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confineSEVals = plotBoxCharts(compiledDistances,cellSE, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledSE','compiledDistances','sampleGrid')
straightSE = compiledSE;
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
cellSE = cell(length(straightSE),1);
for i = 1:length(cellSE)
    curSE = straightSE{i};
    curSE(curSE==0) = nan;
    cellSE{i} = curSE;
end
straightSEVals = plotBoxCharts(compiledDistances,cellSE, sampleGrid, 1,options);

xlabel('Distance to Center of Confinement (microns)')
ylabel('Strain Energy (pJ)')
% xline(-20)
% xline(20)
curYlim = ylim;
pbaspect([2,1,1])
ylim([0 6*10^5])
% pVals = zeros(length(straightSE(:,1)),1);
% for i = 1:length(straightSE(:,1))
%     curStraightSE = straightSE(i,:);
%     curConfineSE = confineSE(i,:);
%     [~,pVals(i)] = ttest2(curStraightSE, curConfineSE);
% end

%T testing
sePvals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = confineSEVals(i,:);
    curStraightSE = straightSEVals(i,:);
    [~,sePvals(i)] = ttest2(curConfineSE,curStraightSE);
end
%% Plot Aspect ratio
figure(2)

load(confineFile,'compiledInterpNAR','compiledDistances','sampleGrid')
nucAR = cell(length(compiledInterpNAR),1);
for i = 1:length(nucAR)
    curNAR = compiledInterpNAR{i};
    curNAR(curNAR==0) = nan;
    nucAR{i} = curNAR;
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,nucAR, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpNAR','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
nucAR = cell(length(compiledInterpNAR),1);
for i = 1:length(nucAR)
    curNAR = compiledInterpNAR{i};
    curNAR(curNAR==0) = nan;
    nucAR{i} = curNAR;
end
plotBoxCharts(compiledDistances,nucAR, sampleGrid, 1,options)


% xline(-20);
% xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Nuclear Aspect Ratio')
pbaspect([2,1,1])
%% Plot SE Density
figure(3)

load(confineFile,'compiledSEDen','compiledDistances','sampleGrid')

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledSEDen, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledSEDen','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledSEDen, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Strain Energy Density (pJ/um^2)')
%% Plot Relative Strain Energy
figure(4)

load(confineFile,'compiledRelSE','compiledDistances','sampleGrid')

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledRelSE, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledRelSE','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledRelSE, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Relative Peak Strain Energy')

%% Plot Speed - FIX THIS
figure(5)

load(confineFile,'compiledNucVel','compiledDistances','sampleGrid')

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledRelSE, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledRelSE','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledRelSE, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Speed (um/min)')
legend('','Straight Channels', '','Confinements')
%% Plot maximum relative strain during transit
% load(straightFile)
% figure(6)
% maxStrainValsStraight = [mean(compiledInterpSE(30:70,:)), nan];
% 
% load(confineFile)
% boxchart([mean(compiledInterpSE(30:70,:));maxStrainValsStraight])

%% Plot Nuclear-Cytoplasm centroid distances
figure(7)
load(confineFile,'compiledNucCellDist','compiledDistances','sampleGrid')

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledNucCellDist, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledNucCellDist','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledNucCellDist, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell/Nucleus centroid distance (um)')
legend('','Straight Channels', '','Confinements')

%% Plot Interpolated Cell Circularity

figure(10)
load(confineFile,'compiledCellCirc','compiledDistances','sampleGrid')
cellCirc = cell(length(compiledCellCirc),1);
for i = 1:length(cellCirc)
    curCirc = compiledCellCirc{i};
    curCirc(curCirc==0) = nan;
    cellCirc{i} = curCirc;
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,cellCirc, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledCellCirc','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
cellCirc = cell(length(compiledCellCirc),1);
for i = 1:length(cellCirc)
    curCirc = compiledCellCirc{i};
    curCirc(curCirc==0) = nan;
    cellCirc{i} = curCirc;
end

plotBoxCharts(compiledDistances,cellCirc, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Circularity')
pbaspect([2,1,1])

%% Plot Interpolated Cell Areas

figure(11)
load(confineFile,'compiledCellArea','compiledDistances','sampleGrid')

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledCellArea, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledCellArea','compiledDistances','sampleGrid')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledCellArea, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Area (um^2)')

%% Plot NCMs
figure(12)
load(confineFile,'compiledNCM','compiledDistances','sampleGrid')
projx = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     projx{i} = curNCM(:,1);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledNCM','compiledDistances','sampleGrid')
projx = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     projx{i} = curNCM(:,1);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('NCMy (pNm)')

figure(13)
load(confineFile,'compiledNCM','compiledDistances','sampleGrid')
NCMy = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     NCMy{i} = curNCM(:,2);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,NCMy, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledNCM','compiledDistances','sampleGrid')
NCMy = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     NCMy{i} = curNCM(:,2);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,NCMy, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('NCMy (pNm)')

figure(14)

load(confineFile,'compiledNCM','compiledDistances','sampleGrid')
NCM = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     NCM{i} = curNCM(:,1) + curNCM(:,2);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,NCM, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledNCM','compiledDistances','sampleGrid')
NCM = cell(length(compiledNCM),1);
for i = 1:length(compiledNCM)
     curNCM = compiledNCM{i};
     NCM{i} = curNCM(:,1) + curNCM(:,2);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,NCM, sampleGrid, 1,options)

xlabel('Distance to Center of Trajectory (um)')
ylabel('NCM (pNm)')

%% Plot Strain Energy - Relative to Initial value
figure(20)
hold on 
load(confineFile,'compiledSE','compiledDistances','sampleGrid')
confineSE = compiledSE;
for i = 1:length(compiledSE)
    curSE = compiledSE{i};
    curSE = curSE/max(curSE(curSE>0),[],'omitnan');
    compiledSE{i} = curSE;
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,compiledSE, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledSE','compiledDistances','sampleGrid')
straightSE = compiledSE;
for i = 1:length(compiledSE)
    curSE = compiledSE{i};
    curSE = curSE/max(curSE(curSE>0),[],"omitnan");
    compiledSE{i} = curSE;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,compiledSE, sampleGrid, 1,options)

xlabel('Distance to Center of Confinement (microns)')
ylabel('Strain Energy (pJ)')

%% Traction Projections

figure(21) %xProj
load(confineFile,'compiledtractionProj','compiledDistances','sampleGrid')
projx = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,1)==0,1) = nan;
     projx{i} = curProj(:,1);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledtractionProj','compiledDistances','sampleGrid')
projx = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,1)==0,1) = nan;
     projx{i} = curProj(:,1);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)
% 
% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum X Forces (Pa)')
pbaspect([2,1,1])

figure(22) %yProj
load(confineFile,'compiledtractionProj','compiledDistances','sampleGrid')
projy = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,2)==0,2) = nan;
     projy{i} = curProj(:,2);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,projy, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledtractionProj','compiledDistances','sampleGrid')
projy = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,2)==0,2) = nan;
     projy{i} = curProj(:,2);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,projy, sampleGrid, 1,options)

% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum Y Forces (Pa)')
ylim([0,5* 10^5]) %Clip data to show trends
pbaspect([2,1,1])

figure(23) %zProj
load(confineFile,'compiledtractionProj','compiledDistances','sampleGrid')
projz = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,3)==0,3) = nan;
     projz{i} = curProj(:,3);
end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
plotBoxCharts(compiledDistances,projz, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledtractionProj','compiledDistances','sampleGrid')
projz = cell(length(compiledtractionProj),1);
for i = 1:length(compiledtractionProj)
     curProj = compiledtractionProj{i};
     curProj(curProj(:,3)==0,3) = nan;
     projz{i} = curProj(:,3);
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
plotBoxCharts(compiledDistances,projz, sampleGrid, 1,options)

% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum Z Forces (Pa)')
pbaspect([2,1,1])

%% Plot xContractile Forces
figure(24)
hold on 
load(confineFile,'compiledxContractile','compiledDistances','sampleGrid')
confinexCon = compiledxContractile;
% for i = 1:length(confinexCon)
%     curxCon = confinexCon{i};
%     curxCon = curxCon/max(curxCon(curxCon>0),[],'omitnan');
%     confinexCon{i} = curxCon;
% end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexCon = plotBoxCharts(compiledDistances,confinexCon, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxContractile','compiledDistances','sampleGrid')
straightSE = compiledxContractile;
% for i = 1:length(straightSE)
%     curxCon = straightSE{i};
%     curxCon = curxCon/max(curxCon(curxCon>0),[],"omitnan");
%     straightSE{i} = curxCon;
% end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightxCon = plotBoxCharts(compiledDistances,straightSE, sampleGrid, 1,options);

xlabel('Distance to Center of Confinement (microns)')
ylabel('xContractile Tractions (Pa)')
pbaspect([2,1,1])

[~,xConpVals] = ttest2(confinexCon', straightxCon');

%% Tractions within 5 microns - Sum
figure(25)
hold on 
load(confineFile,'compiledxForces5micSum','compiledDistances','sampleGrid')
confinexCon5Sum = compiledxForces5micSum;
% for i = 1:length(confinexCon)
%     curxCon = confinexCon{i};
%     curxCon = curxCon/max(curxCon(curxCon>0),[],'omitnan');
%     confinexCon{i} = curxCon;
% end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexCon5Sum = plotBoxCharts(compiledDistances,confinexCon5Sum, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxForces5micSum','compiledDistances','sampleGrid')
straightSE5Sum = compiledxForces5micSum;
for i = 1:length(straightSE5Sum)
    curxCon = straightSE5Sum{i};
    curxCon(curxCon == 0) = nan;
    straightSE5Sum{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightxCon5Sum = plotBoxCharts(compiledDistances,straightSE5Sum, sampleGrid, 1,options);

xlabel('Distance to Center of Confinement (microns)')
ylabel('Sum Contractile Tractions within 5 microns(Pa)')
pbaspect([2,1,1])

[~,xConp5SumVals] = ttest2(confinexCon5Sum', straightxCon5Sum');

%% Tractions within 5 microns - Avg
figure(26)
hold on 
load(confineFile,'compiledxForces5micAvg','compiledDistances','sampleGrid')
confinexCon5Avg = compiledxForces5micAvg;
% for i = 1:length(confinexCon)
%     curxCon = confinexCon{i};
%     curxCon = curxCon/max(curxCon(curxCon>0),[],'omitnan');
%     confinexCon{i} = curxCon;
% end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexCon5Avg = plotBoxCharts(compiledDistances,confinexCon5Avg, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxForces5micAvg','compiledDistances','sampleGrid')
straightSE5Avg = compiledxForces5micAvg;
for i = 1:length(straightSE5Avg)
    curxCon = straightSE5Avg{i};
    curxCon(curxCon == 0) = nan;
    straightSE5Avg{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightxCon5Avg = plotBoxCharts(compiledDistances,straightSE5Avg, sampleGrid, 1,options);

xlabel('Distance to Center of Confinement (microns)')
ylabel('Avg Contractile Tractions within 5 microns(Pa)')
pbaspect([2,1,1])

[~,xConp5AvgVals] = ttest2(confinexCon5Avg', straightxCon5Avg');

%% Tractions within 5 microns - Max
figure(27)
hold on 
load(confineFile,'compiledxForces5micMax','compiledDistances','sampleGrid')
confinexCon5Max = compiledxForces5micMax;
% for i = 1:length(confinexCon)
%     curxCon = confinexCon{i};
%     curxCon = curxCon/max(curxCon(curxCon>0),[],'omitnan');
%     confinexCon{i} = curxCon;
% end
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexCon5Max = plotBoxCharts(compiledDistances,confinexCon5Max, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxForces5micMax','compiledDistances','sampleGrid')
straightSE5Max = compiledxForces5micMax;
for i = 1:length(straightSE5Max)
    curxCon = straightSE5Max{i};
    curxCon(curxCon == 0) = nan;
    straightSE5Max{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightSE5Max = plotBoxCharts(compiledDistances,straightSE5Max, sampleGrid, 1,options);

xlabel('Distance to Center of Confinement (microns)')
ylabel('Max Contractile Tractions within 5 microns(Pa)')
pbaspect([2,1,1])

[~,xConp5MaxVals] = ttest2(confinexCon5Max', straightSE5Max');

%% Displacements - Sum
figure(28)
hold on 
load(confineFile,'compiledxDisps5micSum','compiledDistances','sampleGrid')
confinexDisp5Sum = compiledxDisps5micSum;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexDisp5Sum = plotBoxCharts(compiledDistances,confinexDisp5Sum, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxDisps5micSum','compiledDistances','sampleGrid')
straightDisp5Sum = compiledxDisps5micSum;
for i = 1:length(straightDisp5Sum)
    curxCon = straightDisp5Sum{i};
    curxCon(curxCon == 0) = nan;
    straightDisp5Sum{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightDisp5Sum = plotBoxCharts(compiledDistances,straightDisp5Sum, sampleGrid, 1,options);
xlabel('Distance to Center of Confinement (microns)')
ylabel('Sum Displacements within 5 microns (um)')
pbaspect([2,1,1])

[~,xDisp5SumVals] = ttest2(confinexDisp5Sum', straightDisp5Sum');

%% Displacements - Avg
figure(29)
hold on 
load(confineFile,'compiledxDisps5micAvg','compiledDistances','sampleGrid')
confinexDisp5Avg = compiledxDisps5micAvg;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexDisp5Avg = plotBoxCharts(compiledDistances,confinexDisp5Avg, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxDisps5micAvg','compiledDistances','sampleGrid')
straightDisp5Avg = compiledxDisps5micAvg;
for i = 1:length(straightDisp5Avg)
    curxCon = straightDisp5Avg{i};
    curxCon(curxCon == 0) = nan;
    straightDisp5Avg{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightDisp5Avg = plotBoxCharts(compiledDistances,straightDisp5Avg, sampleGrid, 1,options);
xlabel('Distance to Center of Confinement (microns)')
ylabel('Avg Displacements within 5 microns (um)')
pbaspect([2,1,1])

[~,xDisp5AvgVals] = ttest2(confinexDisp5Avg', straightDisp5Avg');

%% Displacements - Max
figure(30)
hold on 
load(confineFile,'compiledxDisps5micMax','compiledDistances','sampleGrid')
confinexDisp5Max = compiledxDisps5micMax;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
confinexDisp5Max = plotBoxCharts(compiledDistances,confinexDisp5Max, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledxDisps5micMax','compiledDistances','sampleGrid')
straightDisp5Max = compiledxDisps5micMax;
for i = 1:length(straightDisp5Max)
    curxCon = straightDisp5Max{i};
    curxCon(curxCon == 0) = nan;
    straightDisp5Max{i} = curxCon;
end
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
straightDisp5Max = plotBoxCharts(compiledDistances,straightDisp5Max, sampleGrid, 1,options);
xlabel('Distance to Center of Confinement (microns)')
ylabel('Max Displacements within 5 microns (um)')
pbaspect([2,1,1])

[~,xDisp5MaxVals] = ttest2(confinexDisp5Max', straightDisp5Max');