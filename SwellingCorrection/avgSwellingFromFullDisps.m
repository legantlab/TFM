%Averaging displacement fields to get average swelling profile testing
load('averagedSwellingProfileStraight.mat')
elemCents2D(:,3) = elemCents2D(:,3)*.9;
%loop through each match and interpolate swelling to uniform grid
allInterpX = zeros(length(elemCents2D), length(matches));
allInterpY = zeros(length(elemCents2D), length(matches));
allInterpZ = zeros(length(elemCents2D), length(matches));
for i = 1:length(matches)
    curMatch = matches{i};

    [finalLocalizations, trans] = alignOriginStraight(curMatch(:,4:6),0);

%     figure
%     scatter3(finalLocalizations(:,1), finalLocalizations(:,2), finalLocalizations(:,3))
%     hold on 
%     scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))

    curMatch(:,4:6) = finalLocalizations;
    curMatch(:,1:3) = finalLocalizations-displacements{i};

    curInterp = interpDisps({curMatch},elemCents2D);    
    curInterp{1}(isnan(curInterp{1})) = 0;
    allInterpX(:,i) = curInterp{1}(:,4);
    allInterpY(:,i) = curInterp{1}(:,5);
    allInterpZ(:,i) = curInterp{1}(:,6);
% 
%     figure
%     quiver3(curInterp{i}(:,1), curInterp{i}(:,2), curInterp{i}(:,3), curInterp{i}(:,4),curInterp{i}(:,5),curInterp{i}(:,6))
%     hold on 
%     quiver3(finalLocalizations(:,1), finalLocalizations(:,2), finalLocalizations(:,3),displacements{i}(:,1), displacements{i}(:,2), ...
%         displacements{i}(:,3))


end
avgInterpX = median(allInterpX,2);
avgInterpY = median(allInterpY,2);
avgInterpZ = median(allInterpZ,2);

quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), avgInterpX, avgInterpY, avgInterpZ)

figure
%quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), allInterpX(:,73), allInterpY(:,73), allInterpZ(:,73))
quiver3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3), averageDisplacement(:,1), averageDisplacement(:,2),...
    averageDisplacement(:,3))