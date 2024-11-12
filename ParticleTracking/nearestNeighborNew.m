function [reference]=nearestNeighborNew(reference,map,numneighb,~)
%This function computes the distance from each bead in the reference
%configuration to every bead in the deformed configuration. Note, reference and
%deformed configuration can be the same dataset.

%Input argument definitions
% reference = a cell array containing the centroids for the beads in the
% reference configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

% map = a cell array containing the centroids for the beads in the
% reference configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

%numneighb = an integer corresponding to the number of nearest neighbors to
%map or a given bead

%refequalmap = irrelevant, include for legacy purposes

%Output argument definitions
% reference = a cell array containing the original coordinates of the input
% dataset in the first column.  Each consecutive column is the 1st, 2nd,
% 3rd,... up to numneighb coordinates of the nearest neighbors.
%   Author: Max Hockenberry
%   Last Update: 10/23/2024

%Initialize some variables 
numspotsref=length(reference);
%numspotsmap=length(map);
mapMAT=cell2mat(map);
refMAT=cell2mat(reference);
%dist=zeros(numspotsmap,1);

Idx = knnsearch(mapMAT(:,2:4), refMAT(:,2:4),'k',numneighb+2);
reference2 = cell(length(Idx),numneighb+2);
for k = 1:length(Idx)
    curMap = mapMAT(Idx(k,:),:);
    reference2{k,1} = refMAT(k,:);
    for l = 1:numneighb+1   
        reference2{k,l+1} = curMap(l,:);
    end
end

reference = reference2;