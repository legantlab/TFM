%Script that takes a multidimensional .tif (z stack, one channel, time) and
%splits it by time and saves as z stack. Refactored to take multible
%channels at once. 
close all 
%Image Location

I = bfopen('T:\Max\2024_01_24_CommitteeMeetingPrep\ExampleCells\Straight\F02\LA-F02_Cropped_Registered.tif');
OutputFold = 'T:\Max\2024_01_24_CommitteeMeetingPrep\ExampleCells\Straight\F02\LA';
%%
%Read Image specific metadata on dimensions
metaData = I{1,4};

zSize = metaData.getPixelsSizeZ(0).getValue();
tSize = metaData.getPixelsSizeT(0).getValue();
xSize = metaData.getPixelsSizeX(0).getValue();
ySize = metaData.getPixelsSizeY(0).getValue();
% 
xVox = metaData.getPixelsPhysicalSizeX(0).value();    
xVox = xVox.doubleValue;
yVox = metaData.getPixelsPhysicalSizeY(0).value(); 
yVox = yVox.doubleValue;
zVox = metaData.getPixelsPhysicalSizeZ(0).value(); 
zVox = zVox.doubleValue;


%%
%Now loop through the tiff and save each image to an output folder
counter = 1;
for i = 1:tSize
    
   %Data structure goes through z stacks, then through time. We can do this operation with a counter
   
   %Preallocate stack 
   curI = zeros(ySize,xSize,zSize,'uint16');
   for k = 1:zSize
       
       %Use counter to pull current image and save to stack
       curIZ = I{1,1}{counter,1};
       
       curI(:,:,k) = curIZ;
       
       %Increment Counter
       counter = counter + 1;
   end
    metaDataOut = createMinimalOMEXMLMetadata(curI);

    pixelSizeX = ome.units.quantity.Length(java.lang.Double(xVox), ome.units.UNITS.MICROMETER);
    
    metaDataOut.setPixelsPhysicalSizeX(pixelSizeX,0);

    pixelSizeY = ome.units.quantity.Length(java.lang.Double(yVox), ome.units.UNITS.MICROMETER);
    metaDataOut.setPixelsPhysicalSizeY(pixelSizeY,0);
    
    pixelSizeZ = ome.units.quantity.Length(java.lang.Double(zVox), ome.units.UNITS.MICROMETER);
    metaDataOut.setPixelsPhysicalSizeZ(pixelSizeZ,0);
    
    
   disp(['Saving time point ' num2str(i) ' of ' num2str(tSize)])
   
   bfsave(curI, [OutputFold '\Time_' num2str(i,'%04.f') '.ome.tiff'],'metadata',metaDataOut)
   disp('Starting Next Image')
end


%%
% 
% %Get total number of positions
% totPos = size(I);
% totPos = totPos(1);
% 
% %Initiaize cell container the same size as the image
% 
% %Initialize a cell array to hold each position
% allPos = cell(totPos, 1);
% 
% %Loop through positions and process each
% for k = 1:length(allPos)
%     %Cur Position
%     curPos = I{k,1};
%     
%     %Get current image sizes/dimensions
%     
%     %Get size of images from first image? 
%     imgLW = (size(curPos{1})./resizeFac);
%     imgLW = ceil(imgLW);
%     
%     store = zeros(imgLW(1), imgLW(2), length(curPos));
%     
%     for i = 1:length(I{1})
%        curI = I{1}{i};
% 
%        curResize = imresize(curI, 1/resizeFac);
% 
%        store(:,:,i) = curResize;
%    
%     end
%     
%     allPos{k} = store;
% end
% 
% 
% 
% 
% %%
% %Reconstruct
% 
% %Save output, each position seperately
% for j = 1:length(allPos)
%    curPos = allPos{j};
%    bfsave(curPos, ['Downsampled_' num2str(resizeFac) '_Pos' num2str(j) '.tiff'])
%     
% end
% 
% %Visualize
% %imagesc(curPos(:,:,1))
% 
