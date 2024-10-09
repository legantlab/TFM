%Random QC Figures

%% Load image and scatter plot beads on MIPs
% image = loadtiff('T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\Beads\Time_01.ome.tiff');
% 
% %% 
% figure
% imagesc((max(image,[],3)))
% hold on 
% scatter((matches{1}(:,4))/.199 - 8.25, (matches{1}(:,5))/.199 -8,'x','r')


%% Compute displacements and make pngs

for i = 1:length(matches)
    i
    figure
    hold on
    time = i;
    d=displacementVectorScaling;
    scatter(matches{time}(:,1), matches{time}(:,2), 'r.')
    scatter(matches{time}(:,4), matches{time}(:,5), 'g.')
    scatter(x{1}{1}(track{1}{1}==0,1), x{1}{1}(track{1}{1}==0,2), 'r')
    scatter(x{time}{1}(setdiff(1:size(x{time}{1},1),track{1}{1}),1), ...
        x{time}{1}(setdiff(1:size(x{time}{1},1),track{1}{1}),2), 'g')
    quiver(matches{time}(:,1), matches{time}(:,2),...
        d*displacementsWithDrift{time}(:,1), d*displacementsWithDrift{time}(:,2),0, 'g')
    quiver(matches{time}(:,1),matches{time}(:,2),...
        d*displacements{time}(:,1),d*displacements{time}(:,2),0,'b')
    hold off
    saveas(gcf,['T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\displacements\',sprintf('%02d', i ),'.png'])
    close all
end