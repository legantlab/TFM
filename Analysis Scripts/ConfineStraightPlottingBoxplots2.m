%ConfineVsStraightPlotting Scripts
close all
clear
clc
%straightFile = 'T:\Max\2023-12-20\Tiffs\IA32LAStraightValues.mat';
%straightFile = 'straightIA32CompiledResultsBoxChart3.mat';
straightFile = 'C:\Users\LegantLab\Documents\decode_master_repo\TFM\Datasets\2024_09_24_StraightData_AllValues.mat';
%confineFile = 'T:\Max\2023-12-20\Tiffs\IA32LAConfinementValues.mat';
%confineFile = 'IA32ConfineValues.mat';
%confineFile = 'confineIA32CompiledResultsBarCharts.mat';
confineFile = 'C:\Users\LegantLab\Documents\decode_master_repo\TFM\Datasets\2024_09_24_ConfineData_AllValues.mat';
resampleGrid = [0];%[-30,-15, -5, 5, 15, 30];
sampleGrid = resampleGrid;
%% Plot Strain Energy
figure(1)
hold on 
load(confineFile,'compiledInterpSE')

confineSE = compiledInterpSE;
confineSE(confineSE==0) = nan;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = confineSE(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

hold on 

load(straightFile,'compiledInterpSE')
straightSE = compiledInterpSE;
straightSE(straightSE==0) = nan;
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = straightSE(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',10,'XJitter','rand')
end
%straightSEVals = plotBoxCharts(compiledDistances,cellSE, sampleGrid, 1,options);

%xlabel('Distance to Center of Confinement (microns)')
ylabel('Strain Energy (pJ)')
% xline(-20)
% xline(20)
curYlim = ylim;
xlim([-2,2])
pbaspect([1,1,1])
%ylim([0 6*10^5])
% pVals = zeros(length(straightSE(:,1)),1);
% for i = 1:length(straightSE(:,1))
%     curStraightSE = straightSE(i,:);
%     curConfineSE = confineSE(i,:);
%     [~,pVals(i)] = ttest2(curStraightSE, curConfineSE);
% end

%T testing
sePvals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = confineSE(i,:);
    curStraightSE = straightSE(i,:);
    [~,sePvals(i)] = ttest2(curConfineSE,curStraightSE);
end
%% Plot Aspect ratio
figure(2)
hold on 
load(confineFile,'compiledInterpNAR')
nucARConfine = compiledInterpNAR;
nucARConfine(nucARConfine==0) = nan;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = nucARConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

%plotBoxCharts(compiledDistances,nucAR, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpNAR')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
nucARStraight = compiledInterpNAR;
nucARStraight(nucARStraight==0) = nan;
for i = 1:length(sampleGrid)
    curVals = nucARStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%plotBoxCharts(compiledDistances,nucAR, sampleGrid, 1,options)


% xline(-20);
% xline(20);
%xlabel('Distance to Center of Trajectory (um)')
ylabel('Nuclear Aspect Ratio')
pbaspect([1,1,1])
xlim([-2,2])
%T testing
nucPvals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = nucARConfine(i,:);
    curStraightSE = nucARStraight(i,:);
    [~,nucPvals(i)] = ttest2(curConfineSE,curStraightSE);
end
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
hold on
load(confineFile,'compiledInterpCellCirc')
cellCircConfine = compiledInterpCellCirc;
cellCircConfine(cellCircConfine==0) = nan;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = cellCircConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

%plotBoxCharts(compiledDistances,cellCirc, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpCellCirc')
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
cellCircStraight = compiledInterpCellCirc;
cellCircStraight(cellCircStraight==0) = nan;
for i = 1:length(sampleGrid)
    curVals = cellCircStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

%plotBoxCharts(compiledDistances,cellCirc, sampleGrid, 1,options)

%xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Circularity')
pbaspect([1,1,1])
xlim([-2,2])
%T testing
cirPvals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = cellCircConfine(i,:);
    curStraightSE = cellCircStraight(i,:);
    [~,cirPvals(i)] = ttest2(curConfineSE,curStraightSE);
end

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
hold on 
load(confineFile,'compiledInterpxProj')
projxConfine = compiledInterpxProj;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = projxConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpxProj')
projxStraight = compiledInterpxProj;

options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = projxStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%plotBoxCharts(compiledDistances,projx, sampleGrid, 1,options)
% 
% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum X Forces (Pa)')
pbaspect([1,1,1])
xlim([-2,2])
ylim([0.5*10^5, 4.5*10^5])
xTractionProjpVals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = projxConfine(i,:);
    curStraightSE = projxStraight(i,:);
    [~,xTractionProjpVals(i)] = ttest2(curConfineSE,curStraightSE);
end

%%
figure(22) %yProj
hold on
load(confineFile,'compiledInterpyProj')
projyConfine = compiledInterpyProj;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = projyConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end


%plotBoxCharts(compiledDistances,projy, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpyProj')
projyStraight = compiledInterpyProj;

options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = projyStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

%plotBoxCharts(compiledDistances,projy, sampleGrid, 1,options)

% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum Y Forces (Pa)')
% ylim([0,5* 10^5]) %Clip data to show trends
pbaspect([1,1,1])
xlim([-2,2])
ylim([0.5*10^5, 4.5*10^5])
yTractionProjpVals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = projyConfine(i,:);
    curStraightSE = projyStraight(i,:);
    [~,yTractionProjpVals(i)] = ttest2(curConfineSE,curStraightSE);
end
%%
figure(23) %zProj
hold on 
load(confineFile,'compiledInterpzProj')
projzConfine = compiledInterpzProj;

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = projzConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%plotBoxCharts(compiledDistances,projz, sampleGrid, 1,options)

hold on 

load(straightFile,'compiledInterpzProj')
projzStraight = compiledInterpzProj;

options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = projzStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%plotBoxCharts(compiledDistances,projz, sampleGrid, 1,options)

% xlabel('Distance to Center of Trajectory (um)')
% ylabel('Sum Z Forces (Pa)')
pbaspect([1,1,1])
xlim([-2,2])
ylim([0.5*10^5, 4.5*10^5])
zTractionProjpVals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = projzConfine(i,:);
    curStraightSE = projzStraight(i,:);
    [~,zTractionProjpVals(i)] = ttest2(curConfineSE,curStraightSE);
end
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
load(confineFile,'compiledInterpforces5micAvg')
confinexCon5Avg = compiledInterpforces5micAvg;

options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = confinexCon5Avg(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%confinexCon5Avg = plotBoxCharts(compiledDistances,confinexCon5Avg, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledInterpforces5micAvg')
straightSE5Avg = compiledInterpforces5micAvg;
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = straightSE5Avg(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%straightxCon5Avg = plotBoxCharts(compiledDistances,straightSE5Avg, sampleGrid, 1,options);

%xlabel('Distance to Center of Confinement (microns)')
ylabel('Avg Contractile Tractions within 5 microns(Pa)')
pbaspect([1,1,1])
xlim([-2,2])

%[~,xConp5AvgVals] = ttest2(confinexCon5Avg', straightxCon5Avg');
%T testing
xTractionAvgpVals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = confinexCon5Avg(i,:);
    curStraightSE = straightSE5Avg(i,:);
    [~,xTractionAvgpVals(i)] = ttest2(curConfineSE,curStraightSE);
end
%% Tractions within 5 microns - Max
figure(27)
hold on 
load(confineFile,'compiledI','compiledDistances','sampleGrid')
confinexCon5Max = compiledxForces5micMax;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = confinexDisp5Avg(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',2)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%confinexCon5Max = plotBoxCharts(compiledDistances,confinexCon5Max, sampleGrid, 1,options);

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
load(confineFile,'compiledInterpdisp5micAvg')
confinexDisp5Avg = compiledInterpdisp5micAvg;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = confinexDisp5Avg(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%confinexDisp5Avg = plotBoxCharts(compiledDistances,confinexDisp5Avg, sampleGrid, 1,options);

hold on 

load(straightFile,'compiledInterpdisp5micAvg')
straightDisp5Avg = compiledInterpdisp5micAvg;
% for i = 1:length(straightDisp5Avg)
%     curxCon = straightDisp5Avg{i};
%     curxCon(curxCon == 0) = nan;
%     %straightDisp5Avg{i} = curxCon;
% end

options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = straightDisp5Avg(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%straightDisp5Avg = plotBoxCharts(compiledDistances,straightDisp5Avg, sampleGrid, 1,options);
%xlabel('Distance to Center of Confinement (microns)')
ylabel('Avg Displacements within 5 microns (um)')
pbaspect([1,1,1])
xlim([-2,2])

%[~,xDisp5AvgVals] = ttest2(confinexDisp5Avg, straightDisp5Avg);
%T testing
xDispAvgpVals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = confinexDisp5Avg(i,:);
    curStraightSE = straightDisp5Avg(i,:);
    [~,xDispAvgpVals(i)] = ttest2(curConfineSE,curStraightSE);
end

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

%% Max Force
figure(31)
hold on 
load(confineFile,'compiledAvgMaxForce')
maxForceConfine = compiledAvgMaxForce;
maxForceConfine(maxForceConfine==0) = nan;
options = struct('boxColor','r','swarmColor','k','swarmOpacity',.5,'offset',-1);
for i = 1:length(sampleGrid)
    curVals = maxForceConfine(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end

hold on 

load(straightFile,'compiledAvgMaxForce')
maxForceStraight = compiledAvgMaxForce;
maxForceStraight(maxForceStraight==0) = nan;
options = struct('boxColor','b','swarmColor','k','swarmOpacity',.5,'offset',1);
for i = 1:length(sampleGrid)
    curVals = maxForceStraight(i,:);
    curVals(curVals==0) = nan;
    curX = sampleGrid(i)+ options.offset;
    boxchart(curX*ones(length(curVals),1), curVals,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',1)
    swarmchart(curX, curVals,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')
end
%straightSEVals = plotBoxCharts(compiledDistances,cellSE, sampleGrid, 1,options);

%xlabel('Distance to Center of Confinement (microns)')
ylabel('Max Traction (Pa)')
% xline(-20)
% xline(20)
%curYlim = ylim;
pbaspect([1,1,1])
xlim([-2,2])
%ylim([0 6*10^5])
% pVals = zeros(length(straightSE(:,1)),1);
% for i = 1:length(straightSE(:,1))
%     curStraightSE = straightSE(i,:);
%     curConfineSE = confineSE(i,:);
%     [~,pVals(i)] = ttest2(curStraightSE, curConfineSE);
% end

%T testing
maxForcePvals = zeros(length(sampleGrid),1);
for i = 1:length(sampleGrid)
    curConfineSE = maxForceConfine(i,:);
    curStraightSE = maxForceStraight(i,:);
    [~,maxForcePvals(i)] = ttest2(curConfineSE,curStraightSE);
end