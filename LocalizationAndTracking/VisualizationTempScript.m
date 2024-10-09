%VisualizationTempScript
close all

Ipost = loadtiff(fileInfo{1});
Ipre = loadtiff(fileInfo{3});
test = sum(Ipre,3);
XYScaleCorrect = 1;
figure
imagesc(sum(Ipost,3))
hold on
scatter((matches{2}(:,1)/XYScaleCorrect)+1, matches{2}(:,2)/XYScaleCorrect + 1, 200, 'r.')
axis equal
hold off

figure
imagesc(sum(Ipre,3))
hold on
scatter((matches{1}(:,4)/XYScaleCorrect) + 1, matches{1}(:,5)/XYScaleCorrect + 1, 200, 'g.')
axis equal

figure
hold on 
XYScaleCorrect = 1;
scatter((matches{1}(:,4)/XYScaleCorrect) + 1, matches{1}(:,5)/XYScaleCorrect + 1, 200, 'g.')
scatter((matches{2}(:,1)/XYScaleCorrect) + 1, matches{2}(:,2)/XYScaleCorrect + 1, 200, 'r.')
quiver((matches{1}(:,1)/XYScaleCorrect) + 1, matches{1}(:,2)/XYScaleCorrect + 1, displacements{1}(:,1), ...
    displacements{1}(:,2), 0)
axis equal
set(gca, 'YDir','reverse')
% scatter(x{1}{1}(track{1}{1}==0,1), x{1}{1}(track{1}{1}==0,2), 'r')
% scatter(x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),1), ...
%     x{2}{1}(setdiff(1:size(x{2}{1},1),track{1}{1}),2), 'g')
% quiver(matches{1}(:,1)/XYScaleCorrect, matches{1}(:,2)/XYScaleCorrect,...
%     d*displacementsWithDrift{1}(:,1), d*displacementsWithDrift{1}(:,2),0, 'g')
% quiver(matches{1}(:,1)/XYScaleCorrect,matches{1}(:,2)/XYScaleCorrect,...
%     d*displacements{1}(:,1),d*displacements{1}(:,2),0,'k')

%We also need a figure to look at the Z localizations
%%
figure
imagesc(squeeze(sum(Ipost,2)))
hold on
scatter(matches{2}(:,2)/XYScaleCorrect + 1, matches{2}(:,3)/XYScaleCorrect + 1, 200)

%% Output 3D scatter as tiff?
%[xq, yq, zq] = meshgrid(1:.156:227,1:.156:227,1:1.09:56);
[xq, yq] = meshgrid(1:.156:995,1:.156:750);
vq = griddata(matches{2}(:,1), matches{2}(:,2), matches{2}(:,3), xq, yq);
figure
mesh(xq,yq,vq) %Interestingly this makes a pretty convincing mesh of the channel, we might have use for this later but not exactly what I need. 

%% Generate convolved image from 3D points
%We should turn this into a general function and add to workflow
%Get 3d points into voxel space
[sizeX, sizeY, sizeZ] = size(Ipre);
convolvedImage = zeros(sizeX,sizeY,sizeZ);

%Loop through matches list and add points
for i = 1:length(matches{1}(:,1))
    curPoint = matches{1}(i,4:6);
    
    %Convert dimensions to voxels
    
    curPoint = [curPoint(1)/0.156, curPoint(2)/0.156 , curPoint(3)/1.09];
    %Round dimensions to closest voxel
    roundPoint = round(curPoint);
    
    %display(roundPoint)
    convolvedImage(roundPoint(2),roundPoint(1),roundPoint(3)) = 1;

end

%Apply Gaussian with approximate sigma to image
gaussSigmaEst = [1,1,2]; %Rough guestimation based on current imaging settings, we could derive values but I got these
% From 1d fits of the beads in fiji and then refining from there until it
% matched roughly what the beads looked like in the image. 
convolvedImage = imgaussfilt3(convolvedImage,gaussSigmaEst);

figure
imagesc(sum(convolvedImage,3))

%Output image as .tiff
matrix = convolvedImage;
for i=1:sizeZ
  tiff = matrix(:, :, i);
  %convert to 8bit
  tiff = uint16(rescale(tiff,0,65535));
  outputFileName = ['U:\Max\2023_01_18_CellsInChannels\S1\TestFrame\ConvolvedImage', num2str(i,3), '.tiff'];
  imwrite(tiff,outputFileName,'WriteMode', 'overwrite')
end

%This works great!
