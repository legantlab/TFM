%Threshold 3D Volume

%Load z-stack images into seperate frames
I = bfopen(' U:\Max\2023_03_28_IA32ChannelsTryp\S5\C2-S5_Registered.tif');
OutputFold = 'U:\Max\2023_03_28_IA32ChannelsTryp\S5\LA_Segment2';
%% Thresholding

%Threshold each and output as masks
metaData = I{1,4};

zSize = metaData.getPixelsSizeZ(0).getValue();
tSize = metaData.getPixelsSizeT(0).getValue();
xSize = metaData.getPixelsSizeX(0).getValue();
ySize = metaData.getPixelsSizeY(0).getValue();

xVox = metaData.getPixelsPhysicalSizeX(0).value();    
xVox = xVox.doubleValue;
yVox = metaData.getPixelsPhysicalSizeY(0).value(); 
yVox = yVox.doubleValue;
zVox = metaData.getPixelsPhysicalSizeZ(0).value(); 
zVox = zVox.doubleValue;

%Initialize storage array
allImages = cell(tSize, 1);

counter = 1;
for i = 1:tSize
    
   %Data structure goes through z stacks, then through time. We can do this operation with a counter
   
   %Preallocate stack 
   curI = zeros(ySize,xSize,zSize,'uint8');
   for k = 1:zSize
       
       %Use counter to pull current image and save to stack
       curIZ = I{1,1}{counter,1};
       
       curI(:,:,k) = curIZ;
       
       %Increment Counter
       counter = counter + 1;
   end    
   %thresh = graythresh(curI);
   BW = edge3(curI,"sobel",0.05);
      
   %Erode then dilate
   SE = strel('disk',3);
   
   BW = imdilate(BW,SE);
   SE = strel('disk',2);
   BW = imerode(BW,SE);
   
   %Fill holes
   BW = imfill(BW, 'holes');

   %Smooth object
   seD = strel('diamond',1);
   BW = imerode(BW,seD);
   BW = imerode(BW,seD);

   %Filter only single object
   

   CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    filtered_vol = false(size(BW));
    filtered_vol(CC.PixelIdxList{idx}) = true;
   
   BW = filtered_vol;
   allImages{i} = BW;
    
   %Save output
   metaDataOut = createMinimalOMEXMLMetadata(int8(BW));

    pixelSizeX = ome.units.quantity.Length(java.lang.Double(xVox), ome.units.UNITS.MICROMETER);
    
    metaDataOut.setPixelsPhysicalSizeX(pixelSizeX,0);
    
    pixelSizeY = ome.units.quantity.Length(java.lang.Double(yVox), ome.units.UNITS.MICROMETER);
    metaDataOut.setPixelsPhysicalSizeY(pixelSizeY,0);
    
    pixelSizeZ = ome.units.quantity.Length(java.lang.Double(zVox), ome.units.UNITS.MICROMETER);
    metaDataOut.setPixelsPhysicalSizeZ(pixelSizeZ,0);

   disp(['Saving time point ' num2str(i) ' of ' num2str(tSize)])
   
   bfsave(int8(BW), [OutputFold '\Time_' num2str(i,'%02.f') '.ome.tiff'],'metadata',metaDataOut)
   disp('Starting Next Image')

end


%% Visualizations
imshow(allImages{1}(:,:,27))
