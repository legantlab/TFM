function [matchedCoordinates]=displacementCalc_recursive_Max_2024_02_19(reference,deformed,numNeighbors,numVectorsRef,numVectorsDef,tol)
%This function recursively calls displacemenCalc to match beads in the
%reference and deformed configurations.  After each round, matched beads
%are removed from the datasets and the function is called again until fewer
%than a specified percentage of possible matches are confirmed.

%Input argument definitions
% reference = a cell array containing the centroids for the beads in the
% reference configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

% deformed = a cell array containing the centroids for the beads in the deformed
% configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

%numNeighbors = the number of neighbors to scan as a potential match.
%Usually between 20 and 100.  Should be larger if drift and swelling are
%between the datasets.

%numVectorsRef = the number of feature vectors to generate for the
%reference dataset.

%numVectorsDef = the number of feature vectors to generate for the deformed
%dataset.  Should be 1 or 2 larger than numVectorsRef.

%tol = a cutoff value indicating how much more closely the feature vectors of the ideal match need to
%be when compaired to all other potential matches between the reference and deformed configurations in order to
%consider it a correct match. See displacementCalc for more details.

%Output arguement definitions
%matchedCoordinates = an array containing the matched beads in the reference and deformed configurations
%where each column is [indexref, xcoordref ycoordref, zcoordref, indexdef, xcoorddef ycoorddef, zcoorddef]

%   Author: Max Hockenberry
%   Last Update: 10/23/2024
[matchedCoordinates]=displacementCalc_Max_2024_02_19(reference,deformed,numNeighbors,numVectorsRef,numVectorsDef,tol); %initial matching
matchedCoordinates_New=matchedCoordinates; %store all matches

while(length(matchedCoordinates_New)/length(matchedCoordinates)>.0001) %As long as the number of new matches found is greater than .0001 * the number of total matches
    %remove matches from initial dataset
    reference_new=setdiff(reference(:,2:4),matchedCoordinates(:,2:4),'rows');
    deformed_new=setdiff(deformed(:,2:4),matchedCoordinates(:,6:8),'rows');
    %perform matching
    if ~isempty(reference_new)&&~isempty(deformed_new)
        [matchedCoordinates_New]=displacementCalc_Max_2024_02_19([[1:1:length(reference_new)]',reference_new],[[1:1:length(deformed_new)]',deformed_new],numNeighbors,numVectorsRef,numVectorsDef,tol);
    else 
        matchedCoordinates_New=[];
    end

    if ~isempty(matchedCoordinates_New) %check to see if any matches were made
        matchedCoordinates=[matchedCoordinates;matchedCoordinates_New];
    end
end

matchedCoordinates=[[1:1:length(matchedCoordinates)]',matchedCoordinates(:,2:4),[1:1:length(matchedCoordinates)]',matchedCoordinates(:,6:8)]; %assemble output