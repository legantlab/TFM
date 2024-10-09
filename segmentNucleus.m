%Script to segment nucleus from z stacks of SIR-DNA then compute
%confinement metric. 

%% Load tiffs

folder = 'T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\DNA';
outFold = 'T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\DNAseg';
tiffs=dir([folder,'\*.ome.tiff']);
% tiffs=dir([preStem,'\*.tif']);
%tiffs=dir([preStem,preTreatmentFileName]);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));

storeBWs = cell(length(tiffs),1);
storeProps = cell(length(tiffs),1);

for i = 1:length(tiffs)
    i
    curI = loadtiff([tiffs(i).folder,'\',tiffs(i).name]);
    I2D = max(curI,[],3);
    [~,threshold] = edge(I2D,'sobel');
    fudgeFactor = 1;
    BWs = edge(I2D,'sobel',threshold * fudgeFactor);
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    BWsdil = imdilate(BWs,[se90 se0]);
    BWdfill = imfill(BWsdil,'holes');
    
    BWnobord = imclearborder(BWdfill,4);
    
    seD = strel('diamond',2);
    BWfinal = imerode(BWnobord,seD);
    %BWfinal = imerode(BWfinal,seD);
    nucSeg = bwareafilt(BWfinal,1);
    storeBWs{i} = nucSeg;
    imwrite(nucSeg,[outFold,'\','frame',num2str(i,3),'.tiff'])

    nucProp = regionprops(nucSeg,'all');
    storeProps{i} = nucProp;
end

%% Plot Variables
axisRatioStore = zeros(length(storeProps),1);
areaStore = zeros(length(storeProps),1);
for j = 1:length(axisRatioStore)-1
    curProp = storeProps{j};
    axisRatioStore(j) = curProp.MajorAxisLength/curProp.MinorAxisLength;
    areaStore(j) = curProp.Area;
end
figure
plot(linspace(1,length(tiffs),length(tiffs)),axisRatioStore)
xlim([0 15]);
xlabel('Frame')
ylabel('Aspect Ratio')

figure
plot(linspace(1,length(tiffs),length(tiffs)),areaStore)
xlim([0 15]);
xlabel('Frame')
ylabel('Area (pix)')

%%
I = loadtiff('T:\Max\2023-06-22\Tiffs\GreatMovies\ON_F09\DNA\Time_01.ome.tiff');

%Display MIP
I2D = max(I,[],3);

imagesc(I2D)
axis equal
%% Threshold/Segment Nucleus
[~,threshold] = edge(I2D,'sobel');
fudgeFactor = 1;
BWs = edge(I2D,'sobel',threshold * fudgeFactor);
se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');

BWnobord = imclearborder(BWdfill,4);

seD = strel('diamond',2);
BWfinal = imerode(BWnobord,seD);
%BWfinal = imerode(BWfinal,seD);
nucSeg = bwareafilt(BWfinal,1);

imshow(nucSeg)


% level = graythresh(double(I2D));
% 
% BW = imbinarize(I2D,.01);
% 
% imshow(BW)

%% Compute degree of confinement via imageprops
nucProp = regionprops(nucSeg,'all');