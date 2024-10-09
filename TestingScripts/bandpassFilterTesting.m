%Code to assess bandpass filter settings
close all
%Load in image
I = loadtiff('T:\Max\2023_07_01_IA32_SpinningDisk\Tiffs\F13\DMSORefTest\F13_Registered_Rotated_Beads-4.tif');
I2 = loadtiff('T:\Max\2023_07_01_IA32_SpinningDisk\Tiffs\F13\DMSORefTest\F13_Registered_Rotated_Beads-2.tif');

% %Try darkcurrent subtraction
% I = I - 480;
% I(I<0) = 0;
% I2 = I - 480;
% I2(I2<0) = 0;
%Initial Parameters
threshmulti = 1;
bPassParams = cell(3,1);
threshold = 0.90;
%Sweep various bandpass filters
xybPass = 5;%linspace(20,40,5);
zbPass = 10;%linspace(30,70, 9);
bPassParams{1} = [.5,.5,.5];
bPassParams{3} = [0,0];
%Store images
results = cell(length(xybPass),length(zbPass));
particles = cell(length(xybPass),length(zbPass));
for i = 1:length(xybPass)
    curXYbPass = xybPass(i);
    for j = 1:length(zbPass)
        curzbPass = zbPass(j);
        bPassParams{2} = [curXYbPass, curXYbPass, curzbPass];
        [IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3});
        results{i,j} = IBP;

        %localizeParticles
        particles{i,j} = locateParticles(I,threshmulti,bPassParams{2},bPassParams,threshold);

    end


end

%Display results as maximum intensities, both XY and XZ planes to look at
%xy and z reconstructions. 
figure
count = 1
for i = 1:length(xybPass)
    for j = 1:length(zbPass)
        count
        subplot(length(xybPass),length(zbPass),count)
        imagesc(max(results{i,j},[],3))
        hold on 
        scatter(particles{i,j}(:,1)+1, particles{i,j}(:,2)+1,'r')
        title(['XY = ', num2str(xybPass(i)), 'Z = ',num2str(zbPass(j))])
        count = count + 1;
    end
end

figure
count = 1
for i = 1:length(xybPass)
    for j = 1:length(zbPass)
        count
        subplot(length(xybPass),length(zbPass),count)
        imagesc(squeeze(max(results{i,j},[],2)))
        hold on 
        scatter(particles{i,j}(:,3)+1, particles{i,j}(:,2)+1,'r')
        title(['XY = ', num2str(xybPass(i)), 'Z = ',num2str(zbPass(j))])
        count = count + 1;
    end
end


%% Compare two images
%Once settings are selected, do two images and compare localizations
%between them. 
close all
x = cell(2,1);
x{1} = locateParticles(I,threshmulti,bPassParams{2},bPassParams,threshold);
x{2} = locateParticles(I2,threshmulti,bPassParams{2},bPassParams,threshold);

figure
imagesc(max(results{i,j},[],3))
hold on 
scatter(x{1}(:,1)+1, x{1}(:,2)+1, 'r')
scatter(x{2}(:,1)+1, x{2}(:,2)+1, 'y')

% Link particles
%Initial Parameter setting
tptParam{1}.knnFD = 16;
tptParam{1}.fmThres = 3;
tptParam{1}.knnFM = 5;
tptParam{1}.nSpheres = 3;
tptParam{1}.outlrThres = 20;
tptParam{1}.maxIter = 5;
tptParam{1}.sizeI = size(I);

predictor.flag = false;
predictor.x0 = [];
predictor.u = [];

track{j} = TPT(x{1},x{2},tptParam{1},predictor);

initLocations = x{1};
finLocations = x{2};
ind = find(track{1} ~= 0);
map = track{1}(ind);
matches = [initLocations(ind,:) finLocations(map,:)];
%Subtract drift
matches(:,4)=matches(:,4)-mean(matches(:,4)-matches(:,1));
matches(:,5)=matches(:,5)-mean(matches(:,5)-matches(:,2));
matches(:,6)=matches(:,6)-mean(matches(:,6)-matches(:,3));
displacements = [matches(:,4)-matches(:,1), matches(:,5)-matches(:,2), matches(:,6)-matches(:,3)];
% Plot
figure
%imagesc(max(I,[],3))
imshowpair(max(I,[],3),max(I2,[],3))

hold on 
scatter((matches(:,1)+1), (matches(:,2)+1), 'r')
scatter((matches(:,4)+1), (matches(:,5)+1), 'y')
title('Tracked Particles')

%Look in Z
figure
imshowpair(squeeze(max(I,[],2)),squeeze(max(I2,[],2)))

figure
h1 = histogram(displacements(:,1));
hold on 
h2 = histogram(displacements(:,2));
h3 = histogram(displacements(:,3));
h1.Normalization = 'probability';
h1.BinWidth = 0.05;
h2.Normalization = 'probability';
h2.BinWidth = 0.05;
h3.Normalization = 'probability';
h3.BinWidth = 0.05;
xlabel('Displacement (microns)')
ylabel('Probability')
legend('X', 'Y', 'Z')

%Now do histogram fitting for each
figure
hxfit = histfit(displacements(:,1));
xfit = fitdist(displacements(:,1),'normal');
hold on 
hyfit = histfit(displacements(:,2));
yfit = fitdist(displacements(:,2),'normal');

hzfit = histfit(displacements(:,3));
zfit = fitdist(displacements(:,3),'normal');

xlabel('Displacement (microns)')
ylabel('Probability')

title('Histogram Fits')

display(['Xfit = ', num2str(xfit.sigma*1000), ' nm'])
display(['Yfit = ', num2str(yfit.sigma*1000), ' nm'])
display(['Zfit = ', num2str(zfit.sigma*1000), ' nm'])