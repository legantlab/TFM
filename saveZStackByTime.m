function saveZStackByTime(ImagePath,OutputFold)
%saveZStackByTime Saves a .tiff file z stack and seperates into individual
%stacks by time. 
%   Detailed explanation goes here
I = bfopen(ImagePath);
%%
%Read Image specific metadata on dimensions
metaData = I{1,4};

zSize = metaData.getPixelsSizeZ(0).getValue();
tSize = metaData.getPixelsSizeT(0).getValue();
xSize = metaData.getPixelsSizeX(0).getValue();
ySize = metaData.getPixelsSizeY(0).getValue();
% 
% xVox = metaData.getPixelsPhysicalSizeX(0).value();    
% xVox = xVox.doubleValue;
% yVox = metaData.getPixelsPhysicalSizeY(0).value(); 
% yVox = yVox.doubleValue;
% zVox = metaData.getPixelsPhysicalSizeZ(0).value(); 
% zVox = zVox.doubleValue;
% Sometimes need to manually do the pixel sizes, not sure why. 
xVox = 0.199;
yVox = 0.199;
zVox = 0.8;

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

end