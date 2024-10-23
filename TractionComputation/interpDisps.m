function [beadDisps] = interpDisps(matches,appliedNodeList)
%interpDisps Performs scattered interpolant of beads onto the reference
%state grid. 
%{ 
Interpolates the localized beads and displacements onto the reference
%frame to simplify the traction solution. 
Inputs:
- matches: Cell array containing matched bead positions from localization
and tracking code
Outputs:
- beadDisps: Cell array containing the interpolated bead positions (4:6) and
matched reference position (1:3)

%   Author: Max Hockenberry
%   Last Update: 10/23/2024

%}   
    beadDisps = {};
    
    for ii=1:length(matches)
        ii
        curDisp = matches{ii}(:,4:6) - matches{ii}(:,1:3);
        sI1 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,1),'nearest', 'nearest');%'linear','none'%,'nearest', 'nearest');
        sI2 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,2),'nearest', 'nearest');%,'nearest', 'nearest');
        sI3 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,3),'nearest', 'nearest');%,'nearest', 'nearest');
        beadDisps{ii,1}=[appliedNodeList(:,1:3),sI1(appliedNodeList(:,1:3)),sI2(appliedNodeList(:,1:3)),sI3(appliedNodeList(:,1:3))];
        beadDisps{ii,1}(isnan(beadDisps{ii,1})) = 0;

        %Plot if needed
%         figure
%         subplot(1,2,1)
%         quiver3(matches{ii}(:,1), matches{ii}(:,2),matches{ii}(:,3),curDisp(:,1), curDisp(:,2),curDisp(:,3),1,'b')
%         title('Raw Displacements')
%         hold on 
%         subplot(1,2,2)
%         quiver3(beadDisps{ii}(:,1),beadDisps{ii}(:,2),beadDisps{ii}(:,3),beadDisps{ii}(:,4),beadDisps{ii}(:,5),beadDisps{ii}(:,6),1,'b')
%         title('Interpolated Displacements')   
% 
% %         Plot element centroids as compared to matches to see if correction
% %         is needed
%         figure
%         scatter3(matches{ii}(:,1) , matches{ii}(:,2),matches{ii}(:,3),'b')
%         hold on 
%         scatter3(appliedNodeList(:,1), appliedNodeList(:,2), appliedNodeList(:,3),'r')
%         legend('Bead localizations', 'Element Centroids')

    end


end