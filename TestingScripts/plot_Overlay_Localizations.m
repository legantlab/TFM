%Overlay Plot script
%Script that overlays inputted images of cells with the traction map
%generated from TFM to compare traction localization with actual cell
%position. 


%Import Cell image
cellImage = imread('U:\Max\2023_01_18_CellsInChannels\S4\TestFrame2D\OutputGFP\MAX_40x_U2OS_Dox_LA_GFP_Straight_SDS_TimeLapse_006_S4_Registered_Cropped-2Full-1.tif');
beadImg1 = imread('U:\Max\2023_01_18_CellsInChannels\S4\ExampleDataForFigures\BeadsT1.tif');
beadImg2 = imread('U:\Max\2023_01_18_CellsInChannels\S4\ExampleDataForFigures\BeadsT2.tif');
%Plot image
figure
imagesc(cellImage)
hold on 
quiver(elemCents2D(:,1)/.156,elemCents2D(:,2)/.156,force_vector{1}(:,4)/.156,force_vector{1}(:,5)/.156,maxF,'r')


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
% figure
% imshowpair(curImg, reference)
% hold on 
% quiver(matches{1}(:,1)/0.1844 + 1, matches{1}(:,2)/0.1844 + 1, displacements{1}(:,1)/0.1844,displacements{1}(:,2)/0.1844,0)