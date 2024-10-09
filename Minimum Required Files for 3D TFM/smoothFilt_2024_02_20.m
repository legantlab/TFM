function [matchedCoordinates_filt]=smoothFilt_2024_02_20(matchedCoordinates,numNeighbors,tol)
% matchedCoordinates={zeros(200,4),zeros(200,4)};
%This function finds matches beads between reference and deformed
%conditions and outputs the vectors cooresponding to each bead displacement

%Input argument definitions
% reference = a cell array containing the centroids for the beads in the
% reference configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

% deformed = a cell array containing the centroids for the beads in the deformed
% configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

%Note, the length of reference and deformed need not be the same

% numNeighbors=50; %The number of neighbors to scan for initial feature vector generation
% numVectorsRef=3; %The number of feature vectors to match between reference and deformed configurations
% numVectorsDef=5; %The number of feature vectors to search in the deformed configuration

reference=mat2cell(cat(2,[1:1:length(matchedCoordinates)]',matchedCoordinates(:,2:4)),ones(length(matchedCoordinates),1),4);
neighborRef_Ref = nearestNeighborNew(reference,reference,numNeighbors,1);

% 
% 
% tic
% 'mapping reference configuration'
% neighborRef_Ref=nearestNeighborNew(reference,reference,numNeighbors,1); %A mapping of each bead and its neighbors in the reference configuration
% toc
% 'mapping deformed configuration'
% tic
% neighborDef_Def=nearestNeighborNew(deformed,deformed,numVectorsDef,1); %A mapping of each bead and its neighbors in the deformed configuration
% toc
% 'mapping reference to deformed configuration'
% tic
% neighborRef_Def=nearestNeighborNew(reference,deformed,numNeighbors,0); %A mapping of each bead in the reference configuration and its neighbors in the deformed configuration
% toc
t=1; %A counter to keep track of how many beads are matched

%For each bead in the reference configuration, re-sort the feature vectors of each of the
%neighbors in the deformed configuration to give the best match

t=1;
for i=1:length(reference)
    refbead=reference{i}; %index and centroid of the reference bead in the material
    refVector=[matchedCoordinates(refbead(1,1),6:8)-matchedCoordinates(refbead(1,1),2:4)];
    for j=2:numNeighbors+1
        neighbor_bead=neighborRef_Ref{i,j}; %index and centroid of the jth nearest neighbor to the reference bead in the deformed configuration
        neighborVector(j-1,:)=[matchedCoordinates(neighbor_bead(1,1),6:8)-matchedCoordinates(neighbor_bead(1,1),2:4)];% Essentially the vector between a bead and its n neighbors
    end
    dist=sqrt(sum((bsxfun(@minus,neighborVector,refVector)).^2,2));
%     wAvg=sum(bsxfun(@times,1./dist,neighborVector))./sum(1./dist);
    if abs(mean(neighborVector(:,1)-refVector(1)))<(tol*mean(std(neighborVector(:,1))))&&abs(mean(neighborVector(:,2)-refVector(2)))<(tol*mean(std(neighborVector(:,2))))&&abs(mean(neighborVector(:,3)-refVector(3)))<(tol*mean(std(neighborVector(:,3))))
%     if mean(abs(mean(neighborVector)-refVector))<(tol*mean(std(neighborVector)))
        matchedCoordinates_filt(t,:)=matchedCoordinates(refbead(1,1),:);
        t=t+1;
    end
end



