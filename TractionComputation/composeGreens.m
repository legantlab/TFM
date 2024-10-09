function [beadu] = composeGreens(nodeDisps,elemCents2D,TR2)
%composeGreens Arranges nodal displacement and beadDisps into green's
%matrix for CSVD 
%   Composes green's matrix from nodal displacments and bead displacements
%   based on interpolating displacement data to 2D node locations. Each
%   column is subcases while rows are specific nodes displacements i.e.
%   beadu is mxn matrix where m is number of nodes*3 (x,y,z disp) while n
%   is equal to number of subcases. 
    
    numSub = length(elemCents2D(:,1))*3; %Add +1 if thermal swelling correction added. 
    numNodes = length(nodeDisps(:,1))/numSub;
    %% Use shape functions to compute displacements for centroids of beads
    %Loop through each subcase (column) of nodalData
    numElems = length(elemCents2D);
    beadu = zeros(numElems*3,numSub);
    %Reshape connectivity matrix into a column vector 
    connectMatrix = reshape(TR2.ConnectivityList.',length(TR2.ConnectivityList(:,1))*3,[]);

    for i = 1:numSub    
        curSub = nodeDisps(numNodes*(i-1)+1:numNodes*i,1:5);
        [~, locb] = ismember(connectMatrix, curSub(:,1));
        test = curSub(locb,3:5);
        test2 = reshape(test,length(test)*3,[]);
        test3 = mean(reshape(test2,3,[]))';
        %This is right, but needs to be reordered again
        test6 = reshape(test3,[length(test3)/3,3]); %Test functions are just breaking up algebra to work through
        beadu(:,i) = reshape(test6', length(test3),[]);
        
    end


end