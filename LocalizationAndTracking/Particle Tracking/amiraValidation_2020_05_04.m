clear all;
close all;
clc;

trackingDataStem='Z:\Regan\2020_04_28(TFM_troubleshooting)\Analysis\ROI6\';
filename='04_29_20_ROI6_Tracking_.mat';
load([trackingDataStem filename]);
analysisStem='Z:\Regan\2020_04_28(TFM_troubleshooting)\Analysis\ROI6\';

postFilename='postBeadLocations';
preFilename='preBeadLocations';
matchedFilename='matchedBeads';

matched=matches{1,1};

numPost=size(initLocations);
numPre=size(finLocations);

vPost=zeros(numPost(1),numPost(2));
vPre=zeros(numPre(1),numPre(2));

%write post bead locations to a file readable by Amira
% text=buildAmiraMesh(vPost,initLocations);
% fid=fopen([analysisStem postFilename '.am'],'wt');
% fprintf(fid,text);
% fclose(fid);

%write post bead locations to a file readable by Amira
text=buildAmiraMesh(vPre,finLocations);
fid=fopen([analysisStem preFilename '.am'],'wt');
fprintf(fid,text);
fclose(fid);

%write post bead locations to a file readable by Amira
text=buildAmiraMesh(displacements{1,1},matched(:,1:3));
fid=fopen([analysisStem matchedFilename '.am'],'wt');
fprintf(fid,text);
fclose(fid);