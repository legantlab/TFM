function [newIds, TR2,newNodeCoords, surfNodes,elemCents2D] = organizeGeometry(data,triaElemNodeID, triaElemSetID)
%organizeGeometry Inputs abaqus INP reader data structure reorganizes data
%{   
Inputs: 
- data: structure variable output from abaqus reader containing node and
element IDs/positions as sets. 
- triaElemNodeID: the id corresponding to the nodes of the triangular
face elements
- triaElemSetID: the id corresponding to the elements comprising the
triangular face elements
Outputs:
- newIds: An array containing the renumberd Ids of the nodes corresponding
to the 2D face elements
- TR2: A triangulation object of the 2D face elements
- newNodeCoords: An array containing the IDs and coordinates of each node
on the 2D triangular faces
- elemCents2D: An array containing the centroids of each 2D face element
%}
    nodeCoords = [data.Nodes.ID, data.Nodes.Coordinates];
    %Only report nodes from node set we output
    %nodeSet = data.NodeSets.Data;
    %nodeCoords = nodeCoords(ismember(nodeCoords(:,1),nodeSet),:);
    newIds = [nodeCoords(:,1),transpose(linspace(1,length(nodeCoords(:,1)),length(nodeCoords(:,1))))];
    
    newFaceElemConn = data.ElementSets(1).Data;
    
    elemSetConn = zeros([length(newFaceElemConn),4]);
    elemSetConn(:,1) = newFaceElemConn;
    
    triaElements = [data.Elements(1).ID, data.Elements(1).NodeIDList];
    
    for i = 1:length(newFaceElemConn)
        curConn = newFaceElemConn(i);
        %This is the element ID, now we need to get coordinates from triaIDs
        curTriaID = find(curConn == triaElements(:,1));
        curTriaNodes = triaElements(curTriaID,2:4);
        curConn = curTriaNodes;
    
    
        curConnX = curConn(1);
        curConnY = curConn(2);
        curConnZ = curConn(3);
        ConnX = find(curConnX == newIds(:,1));
        ConnY = find(curConnY == newIds(:,1));
        ConnZ = find(curConnZ == newIds(:,1));
        elemSetConn(i,2:4) = [newIds(ConnX,2),newIds(ConnY,2),newIds(ConnZ,2)];
    
    end
    
    TR2=triangulation(elemSetConn(:,2:4),nodeCoords(:,2:4)); %Setup a matlab triangulation object to hold the mesh data
    %Also need to reset nodeCoords
    newNodeCoords = [newIds(:,2),nodeCoords(:,2:4)];
    %Get centers of each 2D element
    [elemCents2D]=elemCentFind2D(newNodeCoords,elemSetConn);
   
    %Get only nodes on 2D surface from elemset
    surfaceNodes = data.NodeSets.Data;

    %relabel nodeIDs based on newIds
    for i = 1:length(surfaceNodes)
        surfaceNodes(i) = newIds(surfaceNodes(i)==newIds(:,1),2);
    end

    surfNodes = surfaceNodes;
    

end