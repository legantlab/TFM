%Script to compute average displacement field on interpolated grid field
close all
%Get .mat files of each displacement field

fold = 'C:\Users\maxho\Desktop\CompiledDisplacements';

fileList = dir(fold);
fileList = fileList(~ismember({fileList.name},{'.','..'}));

%Load field to be interpolated to
load('C:\Users\Max\Documents\GitHub\TFM\averagedSwellingProfile7x40Final.mat');
elemCents2D = elemCents2Dnew;
%Loop through and add each set of matches to a cell array
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

    %% Find optimum shift for each swelling data set

    for j = 1:length(matches)
        curMatch = matches{j};
        curDisp = curMatch(:,4:6) - curMatch(:,1:3);
        %Perform automated Alignment
        [finalElements,finalDisplacements] = alignChannels(swellingData,referenceLocs,0);

        %set translations
        xT = -9;
        yT = 0;
        zT = 1.5;

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
        

    end
        % %Plot for checking
        % figure
        % scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
        % hold on 
        % scatter3(curMatchPre(:,1), curMatchPre(:,2),curMatchPre(:,3),'r')
        % legend('Element Centroids', 'Matched Beads')

    %% Now repeat pulling from the saved rotations
    for j = 1:length(matches)
        curMatch = matches{j};

        %set translations
        xT = trans(1,i);
        yT = trans(2,i);
        zT = trans(3,i);

        xRot = trans(4,i);
        yRot = trans(5,i);
        zRot = trans(6,i);

        %Now do translation

        curMatchPre = rotate_3D(curMatch(:,1:3)','z',zRot)';
        curMatchPost = rotate_3D(curMatch(:,4:6)','z',zRot)';
       
        curMatchPre = [curMatchPre(:,1) + xT, curMatchPre(:,2) + yT, curMatchPre(:,3) + zT];
        curMatchPost = [curMatchPost(:,1) + xT, curMatchPost(:,2) + yT, curMatchPost(:,3) + zT];
        matches{j} = [curMatchPre, curMatchPost];
    end


    rotMatches{i} = matches{1};
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
scatter3(rotMatches{set}(:,4), rotMatches{set}(:,5),rotMatches{set}(:,6),'r')
legend('Element Centroids', 'Matched Beads')
% Compute Average Displacement at each point and plot
avgXDisp = mean(xDisps,2);
avgYDisp = mean(yDisps,2);
avgZDisp = mean(zDisps,2);
averageDisplacement = [avgXDisp, avgYDisp, avgZDisp];

figure
quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),avgXDisp,avgYDisp, avgZDisp,1)
hold on
scatter3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),'b')
scatter3(allMatches{set}(:,1), allMatches{set}(:,2),allMatches{set}(:,3),'r')
legend('Averaged Swelling', 'Element Centroids', 'Beads')
%Sometimes you need to remove outlier data from interp plots which is
%straightforward to do, just plot interpolation and remove those rows from
%newXDisps that have obvious interpolation/tracking/localization failures. 

%% Save
save('averagedSwellingProfile7x40Final.mat','averageDisplacement','elemCents2D','trans')
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
