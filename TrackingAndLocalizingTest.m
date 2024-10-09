%Test Script for localizing and tracking particles effeciently

%Load in data
postStem = '\\LEGANTLABSERVER\LegantLab_DG1_VD1\Max\2021_09_24_F127GelsJR20\New folder\Processed Data\Series5\TestData4ReallySmall\';
postTreatmentFile = 'TestData3zFixedcropReallySmall-3.tif';

preStem = postStem;
if ~isfile([postStem postTreatmentFile])
    error('post treatment file does not exist')
end

tiffs=dir([preStem,'\*.tif']);
%tiffs=dir([preStem,preTreatmentFileName]);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
%Removes the folder and subfolder reference.

fileInfo = cell(1+length(tiffs),2);
fileInfo{1,1} = [postStem postTreatmentFile]; 
fileInfo{1,2} = 0;
for i = 1:length(tiffs)
    fileInfo{i+1,1} = [tiffs(i).folder '\' tiffs(i).name];
    fileInfo{i+1,2} = 0;
end

%Hardcode bPass Parameters and TPT parameters 



%Do bPass filtering and show results


[x, track] = funTPT(fileInfo, beadParam, tptParam, isCZI);

%Do localization and show results



%Do TpT and show results


