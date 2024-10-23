%% Test Script 2024/01/08
load("straightChannelElements.mat")
load('T:\Max\2023-12-20\Tiffs\F18\F18_Control_01.mat')

elemCents2Dnew(:,3) = elemCents2Dnew(:,3)*1;
%% 
figure
scatter3(finLocations(:,1), finLocations(:,2), finLocations(:,3))
hold on 
scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))

%% Point cloud alignment
PC1 = pointCloud(initLocations);
PC2 = pointCloud(elemCents2Dnew);

tform = pcregistericp(PC1,PC2);

%Use alignment to move PC1 to PC2
PC1reg = pctransform(PC1,tform);
figure
pcshow(PC1)
hold on 
pcshow(PC2)
pcshow(PC1reg)
legend('Initial Locations','Elements','Final Locations')

%% Now try on something harder
% Register a messier data set? With swelling to non swelling? 
originalData = initLocations;
originalData(:,1) = originalData(:,1) + trans(1)*.199;
originalData(:,2) = originalData(:,2) + trans(2)*.199;
originalData(:,3) = originalData(:,3) + trans(3)*.8;

originalData =  rotate_3D(originalData','z',trans(6))';

figure
scatter3(originalData(:,1), originalData(:,2), originalData(:,3))
hold on
scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))

PC1 = pointCloud(originalData);
tform = pcregistericp(PC1,PC2);

PC1reg = pctransform(PC1,tform);
figure
%pcshow(PC1)
hold on 
pcshow(PC2)
pcshow(PC1reg)
legend('Initial Locations','Elements','Final Locations')

%% Try this with the confinements
load('T:\Max\2023-11-30\Tiffs\F09\F09_Full')
load('averagedSwellingProfile5x40_2023_12_07_Pixel.mat')

figure
scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))
hold on 
scatter3(initLocations(:,1), initLocations(:,2), initLocations(:,3))

%%
originalData = initLocations;
originalData(:,1) = originalData(:,1)+trans(1)*.199;
originalData(:,2) = originalData(:,2)+trans(2)*.199;
originalData(:,3) = originalData(:,3)+trans(3)*.8;

originalData =  rotate_3D(originalData','z',trans(6))';

PC1 = pointCloud(elemCents2D);
PC2 = pointCloud(originalData);
tform = pcregistericp(PC2, PC1);

PC2reg = pctransform(PC2,tform);

figure
scatter3(PC2reg.Location(:,1), PC2reg.Location(:,2), PC2reg.Location(:,3))
hold on 
scatter3(PC1.Location(:,1), PC1.Location(:,2), PC1.Location(:,3))


