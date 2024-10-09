function [nodeDisps] = readNodalSolutions(nodalSolutionsFile,numElems,newIds)
%readNodalSolutions ReadsNodalSolutionsH5 file and parses into nodal
%displacements
%{
Uses h5read to parse outputed nodal solutions file from hypermesh and
arrange data for nodal displacements array. 
Input: 
- nodalSolutionsFile: String of path to nodal solutions h5 file. 
- elemCents2D: An array containing the centroids of each 2D face element
- newIds: An array containing the renumberd Ids of the nodes corresponding
to the elements
Output:
- nodeDisps: Array containing subcase ID and nodal x/y/z displacements for
each node. 
%}

    nodalSolutionsData = h5read(nodalSolutionsFile,'/OPTISTRUCT/RESULT/NODAL/DISPLACEMENT'); 
    %Get the specific data we need
    %outputStart = find(~isnan(nodalSolutions),1); %Find the starting index of the matrix
    nodeDisps = [double(nodalSolutionsData.ID),zeros(length(nodalSolutionsData.X),1),nodalSolutionsData.X,nodalSolutionsData.Y,nodalSolutionsData.Z]; 
    nodeDisps = double(nodeDisps);
    
    %add subcase info with repelem
    subcaseIDs = 1:(numElems*3); %Add +1 to num elems*3 if have thermal swelling correction
    nodeDisps(:,2) = repelem(subcaseIDs,length(nodeDisps(:,1))/length(subcaseIDs));
    numNodes = length(nodeDisps(:,1))/length(subcaseIDs);

    %Loop through first set of elements and get nodeID list
    newNodeIDs = zeros(numNodes,1);

    for i = 1:length(newNodeIDs)
        newNodeIDs(i) = newIds(nodeDisps(i,1) == newIds(:,1),2);

    end
    %Now just repeat elements a number of times equal to the number of
    %subcases
    nodeDisps(:,1) = repmat(newNodeIDs,length(subcaseIDs),1);
end