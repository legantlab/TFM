%Script to load masks from nice movie and make 3D STLs for Amira

LA = dir('T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\LA\Masks\*.tif');
DNA = dir('T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA\Masks\*.tif');

for i = 1:length(LA)
    curLA = loadtiff([LA(i).folder,'\',LA(i).name]);

    LAmask = regionprops3(logical(curLA),'VoxelList','Volume');
    [~,LAmaskIDX] = max(LAmask.Volume);
    LAPoints = LAmask(LAmaskIDX,:).VoxelList{1};
    
    shp = alphaShape(LAPoints(:,1)*.199, LAPoints(:,2)*.199, LAPoints(:,3)*.8,1);
    [faces, vertices] = boundaryFacets(shp);
    fv.faces = faces;
    fv.vertices = vertices;
    fv = reducepatch(fv,0.01,'fast');
    stlwrite(['T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\LA\MaskSTLs_Smooth\Mask',...
        sprintf('%03d',i),'.stl'],fv)

    curDNA = loadtiff([DNA(i).folder,'\',DNA(i).name]);
    DNAmask = regionprops3(logical(curDNA),'VoxelList','Volume');
    [~,DNAmaskIDX] = max(DNAmask.Volume);
    DNAPoints = DNAmask(DNAmaskIDX,:).VoxelList{1};
    
    shp = alphaShape(DNAPoints(:,1)*.199, DNAPoints(:,2)*.199, DNAPoints(:,3)*.8,1);
    [faces, vertices] = boundaryFacets(shp);
    fv.faces = faces;
    fv.vertices = vertices;
    fv = reducepatch(fv,0.01,'fast');
    stlwrite(['T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F15\DNA\MasksSTLs_Smooth\Mask',...
        sprintf('%03d',i),'.stl'],fv)
end