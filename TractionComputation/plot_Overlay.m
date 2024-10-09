%Overlay Plot script
%Script that overlays inputted images of cells with the traction map
%generated from TFM to compare traction localization with actual cell
%position. 


%Import Cell image
cellImage = imread('U:\Max\2023_01_18_CellsInChannels\S4\TestFrame2D\OutputGFP\MAX_40x_U2OS_Dox_LA_GFP_Straight_SDS_TimeLapse_006_S4_Registered_Cropped-2Full-1.tif');

%Plot image
figure
imagesc(cellImage)
hold on 
quiver(elemCents2D(:,1)/.156,elemCents2D(:,2)/.156,force_vector{1}(:,4)/.156,force_vector{1}(:,5)/.156,maxF,'r')


%% Loop through directory of images

fold = 'U:\Max\2023_01_18_CellsInChannels\S4\AllFramesGFP';

for i = 0:22

    curImg = imread([fold, '\MAX_40x_U2OS_Dox_LA_GFP_Straight_SDS_TimeLapse_006_S4_Registered_Cropped-2Full_',num2str(i,'%04.f'),'.tif']);
    figure
    imagesc(curImg)
    hold on 
    quiver(elemCents2D(:,1)/.156,elemCents2D(:,2)/.156,force_vector{i+1}(:,4)/.156,force_vector{i+1}(:,5)/.156,1,'r')
    axis equal
    set(gca,'visible','off')
    ax= gca;
    exportgraphics(ax,['U:\Max\2023_01_18_CellsInChannels\S4\TestFrame2D\OutputImages\', 'TractionForcesAllMIP_',num2str(i),'.png'],'Resolution',600)

end