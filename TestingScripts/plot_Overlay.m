%Overlay Plot script
%Script that overlays inputted images of cells with the traction map
%generated from TFM to compare traction localization with actual cell
%position. 


%Import Cell image
cellImage = loadtiff('T:\Max\2023-06-22\Tiffs\GreatMovies\ON_F09\MIPs\Beads\Time_01.ome.tiff');


%Threshold quiver values
threshForces = force_vector{1}(:,4:6);
maxF = norm(max(threshForces));
for i = 1:length(threshForces(:,1))
    curRow = threshForces(i,:);
    mag = norm(curRow);
    if mag<.01*maxF
        threshForces(i,:) = 0;
    end

end
% 
% %Plot MIP overlay
% figure
% imagesc(imrotate(max(cellImage,[],3),180))
% hold on 
% quiverC2D(force_vector{1}(:,1)/.199,force_vector{1}(:,2)/.199,threshForces(:,1),threshForces(:,2),1)
% axis equal
% axis off

%Plot image
figure
%subplot(1,2,1)
imagesc(flipud(cellImage))
set(gca,'YDir','normal')
hold on 
%subplot(1,2,2)
%scatter((force_vector{1}(:,1)-9)/0.1933,(force_vector{1}(:,2))/0.1993,'r','.')
quiverC2D((force_vector{1}(:,1)-9)/0.1933,(force_vector{1}(:,2))/0.1993,threshForces(:,1),threshForces(:,2),1)
axis equal


%Plot image
figure
subplot(1,2,1)
imagesc(imrotate(max(cellImage,[],3),180))
subplot(1,2,2)
quiverC3D(force_vector{5}(:,1),force_vector{5}(:,2),force_vector{5}(:,3),threshForces(:,1),threshForces(:,2),threshForces(:,3),0)
axis tight
axis off

%% Loop through directory of images

fold = 'U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S4\bead_MIPS\';

for i = 0:1

    curImg = imread([fold, 'Bead_MIPS_',num2str(i,'%04.f'),'.tif']);
    figure
    imagesc(curImg)

    hold on 
    %quiver(elemCents2D(:,1)/.156 -35,elemCents2D(:,2)/.156 - 75,force_vector{i+1}(:,4)/.156,force_vector{i+1}(:,5)/.156,1,'r')
    %axis equal
    scatter(matches{i+1}(:,1)/0.1844 + 1, matches{i+1}(:,2)/0.1844 + 1,'r')
    set(gca,'visible','off')
    set(gca, "XDir",'reverse')
    axis equal
    
    %set(gcf, 'PaperUnits', 'inches')
    %set(gcf, 'PaperPosition', [0 0 7.25, 9.125])
    ax= gca;
    %xlim([0 1324]) 
    %ylim([0 2136])
   
    %figsave(gcf,['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S4\LA_Forces_Overlay\', 'TractionForcesAllMIP_',num2str(i,'%04.f'),'.png'],[2100 2580])
    
    %exportgraphics(ax,['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S4\LA_Forces_Interp_Overlay\', 'TractionForcesAllMIP_',num2str(i,'%04.f'),'.png'],'Resolution',1000)
    %^Most common save option but sometimes does not have matching img
    %dimensions. 

    %output_size = [1324, 2136];%Size in pixels
    %resolution = 1000;%Resolution in DPI
    %set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
    %axis off
    
    %print(['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S4\LA_Forces_Overlay\', 'TractionForcesAllMIP_',num2str(i,'%04.f'),'.tif'],'-dtiff',['-r' num2str(resolution)]);
    
    %exportgraphics(ax,['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S4\LA_Forces_Overlay_Vec\', ...
        %'TractionForcesAllMIP_',num2str(i,'%04.f'),'.pdf'],'ContentType','vector')
      
    %Export as gif
%     frame = getframe(ax);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);
% 
%     if i == 0
%         imwrite(imind, cm, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\', ...
%         'TractionForcesLA_corrected.gif'],'LoopCount',inf)
% 
%     else
%         imwrite(imind, cm, ['U:\Max\2023_03_03_IA32CellsInChannels\Timelapse_01\S5\', ...
%         'TractionForcesLA_corrected.gif'],'WriteMode','append')
%     end
%     close
end

%% 
figure
imshowpair(curImg, reference)
hold on 
quiver(matches{1}(:,1)/0.1844 + 1, matches{1}(:,2)/0.1844 + 1, displacements{1}(:,1)/0.1844,displacements{1}(:,2)/0.1844,0)