%Script to read .nd2s in a folder and output as .tiff stacks
fold = 'T:\Max\2023_07_18_WF_JR20_PDMS\Stiff';

imageDir = dir([fold, '\*.nd2']);

for k = 1:length(imageDir)
    data = bfopen([imageDir(k).folder, '\', imageDir(k).name]);
    xy = size(data{1}{1});
    temp = zeros(xy(1),xy(2),length(data{1}));
    for i = 1:length(data{1})
        curFrame = data{1}{i};
        temp(:,:,i) = curFrame;      
    end

    % Export in TIFF 
    write3Dtiff(uint16(temp), [imageDir(k).folder, '\', imageDir(k).name,'.tiff'])

end
