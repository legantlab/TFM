function [matchIndex]=matchFeatureVector_2024_02_19(refVectors,defVectors,numVectorsRef,numVectorsDef,numNeighbors,tol)
%This function find the optimal permutations of matches for feature vectors
%in the reference and deformed configuration (i.e. determine the best
%combination/ordering of feature vectors in the deformed configuration to
%match to the reference configuration).

%Input argument definitions
%refVectors = A n by 3*nVR array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z,...,u3NVR_x, u3NVR_y, u3NVR_z].  Each row is a bead (up to n reference beads) and u's are Cartesian components of each feature vector.
%defVectors = A n by 3*nVD*numNeighbors array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, ...,u3NVD*numNeighbors_x, u3NVD*numNeighbors_y, u3NVD*numNeighbors_z]. Each row is vectors for the potential matches to a given bead in the reference configuration.
%   Author: Max Hockenberry
%   Last Update: 10/23/2024
C=genPerms(numVectorsDef,numVectorsRef); %Generate perumutations for matching vectors
numPerms=length(C);
PV=zeros(numPerms,numVectorsRef*3);
matchIndex=zeros(length(refVectors),1);

for i=1:numPerms %Convert permutation matrix to match indices of matrices storing refVectors and defVectors  (eg taking 1 --> 1,2,3 and 2 --> 4,5,6 and 4 --> 10,11,12 etc)
    t=1;
    for j=1:numVectorsRef
        PV(i,t:t+2)=C(i,j)*3-2:C(i,j)*3;
        t=t+3;
    end
end

residual=zeros(numPerms,numNeighbors);
for i=1:length(refVectors)
    refVec = refVectors(i,:);
    for j=1:numNeighbors
        
        for k=1:numPerms
            %this is really slow, but we do call it a million times so maybe no way to get around that.
            residual(k,j)=(norm(refVec-defVectors(i,(j-1)*3*numVectorsDef+PV(k,:))));

        end
    end
    [val,ind]=sort(min(residual));
    if tol*val(1)<val(2)   %Threshold for how more closely the vectors need to match when compared to all other possible matches
        matchIndex(i)=4*ind(1)+1;
    end

end