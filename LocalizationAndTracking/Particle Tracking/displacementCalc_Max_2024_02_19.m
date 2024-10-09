function [matchedCoordinates]=displacementCalc_Max_2024_02_19(reference,deformed,numNeighbors,numVectorsRef,numVectorsDef,tol)
%This function matches beads between reference and deformed
%conditions and outputs the matched centroids.

%Input argument definitions
% reference = a matrix containing the centroids for the beads in the
% reference configuration of the form [bead index, xcoord, ycoord, zcoord].  Each row
% corresponds to one bead.  "bead index" is a numerical ordering to keep
% track of each bead (based on the order in which the centroids were
% found).

% deformed = a matrix containing the centroids for the beads in the deformed
% configuration of the form [bead index, xcoord, ycoord, zcoord].  Each row
% corresponds to one bead. Note, the length of reference and deformed need
% not be the same.

%numNeighbors = the number of neighbors to scan as a potential match.
%Usually between 20 and 100.  Should be larger if drift and swelling are
%present between the datasets.

%numVectorsRef = the number of feature vectors to generate for the
%reference dataset.  A good estimate is 2 or 3.

%numVectorsDef = the number of feature vectors to generate for the deformed
%dataset.  Should be 1 or 2 larger than numVectorsRef.

%tol = a cutoff value indicating how much more closely the feature vectors of the ideal match need to
%be when compaired to all other potential matches between the reference and deformed configurations in order to
%consider it a correct match. Typically a cutoff between 1.25 and 2 is
%appropriate.

%Output arguement definitions
%matchedCoordinates = an array containing the matched beads in the reference and deformed configurations
%where each column is [bead_indexref, xcoordref ycoordref, zcoordref, bead_indexdef,
%xcoorddef ycoorddef, zcoorddef]


%Initialize variables and convert to cell arrays for nearest neighbor
%matching
matchedCoordinates={};
reference=mat2cell(cat(2,[1:1:length(reference)]',reference(:,2:4)),ones(length(reference),1),4);
deformed=mat2cell(cat(2,[1:1:length(deformed)]',deformed(:,2:4)),ones(length(deformed),1),4);

%Parallize this?
neighborRef_Ref=nearestNeighborNew(reference,reference,numVectorsRef,1);
neighborDef_Def=nearestNeighborNew(deformed,deformed,numVectorsDef,1);
neighborRef_Def=nearestNeighborNew(reference,deformed,numNeighbors,1);

% t=1; %A counter to keep track of how many beads are matched

%For each bead in the reference configuration, re-sort the feature vectors of each of the
%neighbors in the deformed configuration to give the best match

%Convert cell arrays back into matrices
reference=cell2mat(reference);
deformed=cell2mat(deformed);
neighborRef_Ref=cell2mat(neighborRef_Ref);
neighborDef_Def=cell2mat(neighborDef_Def);
neighborRef_Def=cell2mat(neighborRef_Def);

%Initialize some variables
refVectors=zeros(length(reference),3*numVectorsRef); %A n by 3*nVR array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, ...,u3NVR_x, u3NVR_y, u3NVR_z].  Each row is a bead (up to n reference beads) and u's are Cartesian components of each feature vector.
%matchIndex=zeros(length(reference),1); %A vector to store the indicies of the matched beads

%Generate feature vectors for each bead of the reference configuration
t=1;
for i=4:4:4*numVectorsRef
    refVectors(:,t:t+2)=neighborRef_Ref(:,i+2:i+4)-neighborRef_Ref(:,2:4); 
    t=t+3;
end

%Generate feature vectors for each potential match in the deformed configuration
defVectors=zeros(length(reference),4*numNeighbors);%A n by 3*nVD*numNeighbors array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, ...,u3NVD*numNeighbors_x, u3NVD*numNeighbors_y, u3NVD*numNeighbors_z].
                                                   %Each row is vectors for the potential matches to a given bead in the reference configuration.
t=1;
for k=4:4:4*numNeighbors+1
    for i=4:4:4*numVectorsDef
        defVectors(:,t:t+2)=neighborDef_Def(neighborRef_Def(:,k+1),i+2:i+4)-neighborDef_Def(neighborRef_Def(:,k+1),2:4);
        t=t+3;
    end
end

%This is what should be parallized I believe
index=matchFeatureVector_2024_02_19(refVectors,defVectors,numVectorsRef,numVectorsDef,numNeighbors,tol);

%Some comment here! Thanks Wes - MH 2024
matchIndex=zeros(length(refVectors),1);
for i=1:length(matchIndex)
    if index(i)>0
        matchIndex(i)=neighborRef_Def(i,index(i));
    end
end
    matchedCoordinates=[reference(matchIndex>0,1:4),deformed(matchIndex(matchIndex>0),1:4)];
end




