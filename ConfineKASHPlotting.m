%ConfineVsStraightPlotting Scripts
close all
clear
clc
%straightFile = 'T:\Max\2023-12-20\Tiffs\IA32LAStraightValues.mat';
straightFile = 'ConfineDataControl_5mic.mat';
%confineFile = 'T:\Max\2023-12-20\Tiffs\IA32LAConfinementValues.mat';
%confineFile = 'IA32ConfineValues.mat';
confineFile = 'ConfineDataDox_5mic.mat';
%% Plot Strain Energy
load(straightFile)

options = struct();
figure(1)
options.handle = 1;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpSE',options)   
hold on 
straightSE = compiledInterpSE;

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpSE',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Strain Energy (pJ)')
legend('','Control', '','DN KASH')
xlim([-40,40])
confineSE = compiledInterpSE;

pVals = zeros(length(straightSE(:,1)),1);
for i = 1:length(straightSE(:,1))
    curStraightSE = straightSE(i,:);
    curConfineSE = confineSE(i,:);
    [~,pVals(i)] = ttest2(curStraightSE, curConfineSE);
end

%% Plot Aspect ratio
load(straightFile)

options = struct();
figure(2)
options.handle = 2;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpNARHull',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpNARHull',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Nuclear Aspect Ratio')
legend('','Control', '','DN KASH')
xlim([-40,40])
%% Plot SE Density
load(straightFile)

options = struct();
figure(3)
options.handle = 3;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpSEDen',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpSEDen',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Strain Energy Density (pJ/um^2)')
legend('','Control', '','DN KASH')
xlim([-40,40])
%% Plot Relative Strain Energy
load(straightFile)

options = struct();
figure(4)
options.handle = 4;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledHullInterpSE',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledHullInterpSE',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Relative Peak Strain Energy')
legend('','Control', '','DN KASH')
xlim([-40,40])

%% Plot Speed
load(straightFile)

options = struct();
figure(5)
options.handle = 5;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpCellVel',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpCellVel',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Speed (um/min)')
legend('','Control', '','DN KASH')
xlim([-40,40])
ylim([0,10])

%% Plot maximum relative strain during transit
% load(straightFile)
% figure(6)
% maxStrainValsStraight = [mean(compiledInterpSE(30:70,:)),nan];
% 
% load(confineFile)
% boxchart([[mean(compiledInterpSE(30:70,:))];maxStrainValsStraight])

%% Plot Nuclear-Cytoplasm centroid distances
load(straightFile)

options = struct();
figure(7)
options.handle = 7;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpNucCellDist',options)   
hold on 
figure(8)
for i = 1:5
    plot(sampleGrid, compiledInterpNucCellDist(:,i))
    hold on 
end

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpNucCellDist',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell/Nucleus centroid distance (um)')
legend('','Control', '','DN KASH')
xlim([-40,40])

%Individual Plots
figure(9)
for i = 1:6
    plot(sampleGrid, compiledInterpNucCellDist(:,i))
    hold on 
end
hold off
%% Plot Interpolated Cell Circularity

load(straightFile)

options = struct();
figure(10)
options.handle = 10;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpCellCirc',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpCellCirc',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Circularity')
legend('','Control', '','DN KASH')
xlim([-40,40])
ylim([0,1])

%% Plot Interpolated Cell Areas

load(straightFile)

options = struct();
figure(11)
options.handle = 11;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpCellArea',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpCellArea',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('Cell Area (um^2)')
legend('','Control', '','DN KASH')
xlim([-40,40])
%ylim([0,1])

%% Plot NCMs
load(straightFile)

options = struct();
figure(12)
options.handle = 12;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpNCMx',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpNCMx',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('NCMx (pNm)')
legend('','Control', '','DN KASH')
xlim([-40,40])
%ylim([0,1])

load(straightFile)

options = struct();
figure(13)
options.handle = 13;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar(compiledInterpNCMy',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(compiledInterpNCMy',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('NCMy (pNm)')
legend('','Control', '','DN KASH')
xlim([-40,40])
%ylim([0,1])

load(straightFile)

options = struct();
figure(14)
options.handle = 14;
options.error = 'sem';
options.x_axis = sampleGrid;
plot_areaerrorbar((compiledInterpNCMx + compiledInterpNCMy)',options)   
hold on 

load(confineFile)
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar((compiledInterpNCMx + compiledInterpNCMy)',options)   

xline(-20);
xline(20);
xlabel('Distance to Center of Trajectory (um)')
ylabel('NCM (pNm)')
legend('','Control', '','DN KASH')
xlim([-40,40])
ylim([0,1])

%% Plot Strain Energy bar charts
load(straightFile)
straightDisps = (allDisps*10);
straightSE = allSE;
load(confineFile)
confineDisps = (allDisps*10);
confineSE = allSE;

% Run statistics on bins
pvals = zeros(11,1);
for i = 1:length(pvals)
    curBin = -50 + 10*(i-1);
    curSEStraight = straightSE(straightDisps==curBin);
    curSEConfine = confineSE(confineDisps==curBin);
    [~,pvals(i)] = ttest2(curSEStraight,curSEConfine);

end


figure(15)
%multiple_boxplot({straightSE, confineSE})
boxplot([straightSE;confineSE],[straightDisps; confineDisps+5],'ColorGroup',[zeros(length(straightSE),1);ones(length(confineSE),1)])
ylim([0, 5*10^5])
xlim([0.5 27])
xlabel('Nuclear Distance to Center of Confinement (um)')
ylabel('Strain Energy (pJ)')
hold on 
swarmchart((((straightDisps/10) +7)),straightSE,10,'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
swarmchart((((confineDisps/10) +8)),confineSE,10,'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

allDisps = categorical(allDisps*10);
boxplot(allSE, allDisps)
hold on 
swarmchart(((allDisps+7)),allSE,10,'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

load(confineFile)
allDisps = categorical(allDisps*10);
boxplot(allSE, (allDisps)*10)
hold on 
swarmchart(((allDisps+8)),allSE,10,'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

ylim([0, 5*10^5])
%xlim([2.5 13.5])
xlabel('Nuclear Distance to Center of Confinement (um)')
ylabel('Strain Energy (pJ)')
