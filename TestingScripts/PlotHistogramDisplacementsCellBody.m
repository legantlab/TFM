%Script that computes displacement histogram away from cell body
close all
%Load cell body mask in 3D
cellMask = fliplr(rot90(loadtiff('U:\Max\2023_03_28_IA32ChannelsTryp\S5\LA_Segment_Flip\Time_01.ome.tiff'),2));
cellMask = imtranslate(cellMask,[xCor/.156-45, yCor/.156, zCor/.8]);
%Plot MIP of cellMask
cellMaskMIP = (max(cellMask,[],3));

imagesc(cellMaskMIP)
hold on 
%Need to make sure are using corrected displacements

quiver(correctedbeadDisps{1}(:,1)/.156 + 1, correctedbeadDisps{1}(:,2)/.156 + 1,correctedbeadDisps{1}(:,4), correctedbeadDisps{1}(:,5),'r')

%%
%Get list of voxel coordinates of 3D object
voxList = regionprops3(cellMask,"VoxelList");
voxList = voxList.VoxelList{1};

%Compute distances between each point in mask and the localizations
matchDisps = zeros(length(correctedbeadDisps{1}(:,1)),1);

for i = 1:numel(correctedbeadDisps{1}(:,1))
    curPoint = correctedbeadDisps{1}(i,1:3);
    curPoint(1) = curPoint(1)/.156 + 1;
    curPoint(2) = curPoint(2)/.156 + 1;
    curPoint(3) = curPoint(3)/.8 + 1;
    temp = inf;
    for j = 1:numel(voxList(:,1))
       x = voxList(j,1) - curPoint(1);
       y = voxList(j,2) - curPoint(2);
       z = voxList(j,3) - curPoint(3);
       r = sqrt(x^2 + y^2 + z^2);

       if r < temp
           temp = r;
       end
    matchDisps(i) = temp;

    end
    
end

%% Histograms
displacementsMag = correctedbeadDisps{1}(:,4:6);
displacementsMag = vecnorm(displacementsMag,2,2);
figure
subplot(2,1,1)
scatter(matchDisps,displacementsMag)
xlabel('Distance from Cell (microns)')
ylabel('Bead Displacements (microns)')
hold on 
expfit = fit(matchDisps, displacementsMag,'power1');
plot(expfit,matchDisps,displacementsMag)
subplot(2,1,2)
boxplot(displacementsMag,'Orientation','horizontal','PlotStyle','traditional')
axis off
axis tight
%%
%Also Compute for X, Y, Z seperately
figure
scatter(matchDisps,abs(correctedbeadDisps{1}(:,4)))
title('X Displacement Decay')
xlabel('Distance from Cell (microns)')
ylabel('Bead Displacements (microns)')

figure
scatter(matchDisps,abs(correctedbeadDisps{1}(:,5)))
title('Y Displacement Decay')
xlabel('Distance from Cell (microns)')
ylabel('Bead Displacements (microns)')

figure
scatter(matchDisps,abs(correctedbeadDisps{1}(:,6)))
title('Z Displacement Decay')
xlabel('Distance from Cell (microns)')
ylabel('Bead Displacements (microns)')