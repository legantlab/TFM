%Script to compute average displacement field on interpolated grid field
close all
%Get .mat files of each displacement field

fold = 'T:\Max\2023-11-30\Tiffs\HR_Straight\F12\Compiled';

fileList = dir(fold);
fileList = fileList(~ismember({fileList.name},{'.','..'}));

%Load field to be interpolated to
%load('C:\Users\Max\Documents\GitHub\TFM\averagedSwellingProfile7x40Final.mat');
fileName = 'T:\Max\2023-11-30\Tiffs\HR_Straight\StraightChannel_Coarsened_Set_Loaded_SolverDeck.inp';
data = abaqusInpRead(fileName);
[newIds, TR2,newNodeCoords, surfNodes,elemCents2Dnew] = organizeGeometry(data,1,1);

elemCents2D = elemCents2Dnew;
%% Loop through and add each set of matches to a cell array
allMatches = cell(length(fileList(:,1)),1);
xDisps = zeros(length(elemCents2D),length(allMatches));
yDisps = zeros(length(elemCents2D),length(allMatches));
zDisps = zeros(length(elemCents2D),length(allMatches));

trans = zeros(6,length(allMatches));

for i = 1:length(allMatches)   
    curFile = [fileList(i).folder, '\', fileList(i).name];
    load(curFile)

    %Need a registration/Alignment step
    %Loop through each set of matches and replace each with a rotated and
    %translated set using rotate_3D and simple transformations. 
    
    % Once interpolations are determined, save to a variable so we can call
    % later. 

    % Find optimum shift for each swelling data set

    for j = 1:length(matches)
        %%
        close all
        curMatch = matches{j};
        curDisp = curMatch(:,4:6) - curMatch(:,1:3);
        
        figure
        scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))
        hold on 

        %Perform manual Alignment
        %set translations
        xT = 0;
        yT =-80;
        zT = 4.5;

        xRot = 0;
        yRot = 0;
        zRot = 0;  

        %Now do translation

        curMatchPre = rotate_3D(curMatch(:,1:3)','z',zRot)';
        curMatchPost = rotate_3D(curMatch(:,4:6)','z',zRot)';
             
        curMatchPre = [curMatchPre(:,1) + xT, curMatchPre(:,2) + yT, curMatchPre(:,3) + zT];
        curMatchPost = [curMatchPost(:,1) + xT, curMatchPost(:,2) + yT, curMatchPost(:,3) + zT];
        

        % Save final translations
        trans(1,i) = xT;
        trans(2,i) = yT;
        trans(3,i) = zT;
        trans(4,i) = xRot;
        trans(5,i) = yRot;
        trans(6,i) = zRot;

        scatter3(curMatchPost(:,1), curMatchPost(:,2), curMatchPost(:,3))
        legend('Element Centroids', 'Translated Matched Beads')
        axis equal
        
        %% Set matches to new pre and post
        matches{j} = [curMatchPre,curMatchPost];
        
        

    end

    %% Now repeat pulling from the saved rotations
    curInterp = interpDisps(matches,elemCents2D);
    allMatches{i} = curInterp{1};

    xDisps(:,i) = curInterp{1}(:,4);
    yDisps(:,i) = curInterp{1}(:,5);
    zDisps(:,i) = curInterp{1}(:,6);

end

%% Plot without any translational corrections
set = 1;
figure
scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
hold on 
scatter3(matches{set}(:,4), matches{set}(:,5),matches{set}(:,6),'r')
legend('Element Centroids', 'Matched Beads')
% Compute Average Displacement at each point and plot
avgXDisp = mean(xDisps,2)/.199;
avgYDisp = mean(yDisps,2)/.199;
avgZDisp = mean(zDisps,2)/.209;
averageDisplacement = [avgXDisp, avgYDisp, avgZDisp];

figure
quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),avgXDisp,avgYDisp, avgZDisp,1)
hold on
scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
scatter3(matches{set}(:,4), matches{set}(:,5),matches{set}(:,6),'r')
legend('Averaged Swelling', 'Element Centroids', 'Beads')
%Sometimes you need to remove outlier data from interp plots which is
%straightforward to do, just plot interpolation and remove those rows from
%newXDisps that have obvious interpolation/tracking/localization failures. 
newCentroids = elemCents2D;
finalDisps = averageDisplacement;

%% Save
save('averagedSwellingProfileStraightPixel.mat','finalDisps','newCentroids','trans')
%% Plot element centroids to matches to compute transformation
match1 = matches{1};
rotmatch1 = rotate_3D(match1(:,1:3)','z',.1);
rotmatch1 = rotmatch1';

figure
scatter3(rotmatch1(:,1)-32.5, rotmatch1(:,2) - 5, rotmatch1(:,3) +1,'r')
hold on 
scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
legend('Rotated Matches','Element Centroids')
%% Plot without any translational corrections
figure
scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
hold on 
scatter3(matches{1}(:,1), matches{1}(:,2),matches{1}(:,3)+2,'r')
legend('Element Centroids', 'Matched Beads')
%% Plot all bead positions to look for deviations
figure
for i = 1:length(rotMatches)
    curMatches = rotMatches{i};
    scatter3(curMatches(:,4), curMatches(:,5),curMatches(:,6))

end
