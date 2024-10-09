%Code to compute average displacement field for swelling corrections.
%Load swelling data set nx2 cell array. Use raw localizations and align
%each point cloud using align origin. Then run tracking again using TPT GUI
%and get swelling correction interpolated to a uniform grid? 

fold = 'T:\Max\2023-12-20\Tiffs\F18\Compiled';

fileList = dir(fold);
fileList = fileList(~ismember({fileList.name},{'.','..'}));

%Load field to be interpolated to
% load('C:\Users\Max\Documents\GitHub\TFM\averagedSwellingProfile7x40Final.mat');

%Loop through each file and redo localization and tracking? Should just
%need to move each match and then redo computation of displacements...
newMatches = cell(2,1);
bigMatch = [];
newDisps = cell(2,1);
bigDisp = [];
for i = 1:length(fileList)
    curfile = [fileList(i).folder,'\',fileList(i).name];
    load(curfile)

    %align matches{1} to origin
    [curReference, ~] = alignOriginStraight(matches{1}(:,1:3),0);
    [curSwell, ~] = alignOriginStraight(matches{1}(:,4:6),0);

    %Plot to confirm
    figure
    scatter3(curReference(:,1), curReference(:,2), curReference(:,3))
    hold on 
    scatter3(curSwell(:,1), curSwell(:,2), curSwell(:,3))
    
    %compute displacements for averaging
    newDisp = curReference - curSwell;
    %Interpolate to uniform grid for averaging? 

    quiver3(curSwell(:,1), curSwell(:,2), curSwell(:,3),...
        newDisp(:,1), newDisp(:,2), newDisp(:,3),0)
    legend('Reference Beads','Data','Swelling')
    hold off
    newMatches{i} = [curReference, curSwell];
    newDisps{i} = newDisp;
    bigMatch = [bigMatch; [curReference, curSwell]];
    bigDisp = [bigDisp,; newDisp];
end

%% Construct large data set and interpolate to 3D model
%Load elements data
load('straightChannelElements.mat')
elemCents2D = elemCents2Dnew;
elemCents2D(:,3) = elemCents2D(:,3)*0.9;
%element centroids to origin
[newCentroids, ~] = alignOriginStraight(elemCents2D,0);
hold on 
scatter3(newCentroids(:,1), newCentroids(:,2), newCentroids(:,3))
%Scale centroids z to account for height differences? 


%Perform interpolation?
[beadDisps] = interpDisps({bigMatch},newCentroids);
finalDisps = beadDisps{1}(:,4:6);
figure
quiver3(beadDisps{1}(:,1), beadDisps{1}(:,2),beadDisps{1}(:,3),beadDisps{1}(:,4),...
    beadDisps{1}(:,5),beadDisps{1}(:,6),0)
hold on 
quiver3(curSwell(:,1), curSwell(:,2), curSwell(:,3),...
        newDisp(:,1), newDisp(:,2), newDisp(:,3),0)
scatter3(beadDisps{1}(:,1), beadDisps{1}(:,2),beadDisps{1}(:,3))
scatter3(bigMatch(:,1), bigMatch(:,2), bigMatch(:,3))
%[finalElements,finalDisplacements,trans] = alignChannels([bigMatch(:,1:3), bigDisp],newCentroids,1);
%save('averagedSwellingProfile7x40_2023_09_01_summedData.mat','newCentroids','finalDisps')

%Compute pixelwise Swelling correction
newCentroids(:,1) = newCentroids(:,1)/.199;
newCentroids(:,2) = newCentroids(:,2)/.199;
newCentroids(:,3) = newCentroids(:,3)/.8;

finalDisps(:,1) = finalDisps(:,1)/.199;
finalDisps(:,2) = finalDisps(:,2)/.199;
finalDisps(:,3) = finalDisps(:,3)/.8;
figure
scatter3(newCentroids(:,1), newCentroids(:,2), newCentroids(:,3))
%save('averagedSwellingProfileStraight_2024_01_04_Pixel.mat','newCentroids','finalDisps')