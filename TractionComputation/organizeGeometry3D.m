function [newElem3DConn,beadDispsCleaned, beadElems] = organizeGeometry3D(data,newIds, newNodeCoords, matches,TR2)
%organizeGeometry3D Organizes 3D geometry from abaqus reader data variable
%{
Organizes the abaqus reader data 3D elements and extracts the element ID
that each displacement data point is located in. 
Input:
- data: structure variable output from abaqus reader containing node and
element IDs/positions as sets. 
- newIDs: An array containing the renumberd Ids of the nodes corresponding
to the 2D face elements
- matches: Cell array containing matched bead positions from localization
and tracking code
- newNodeCoords: An array containing the IDs and coordinates of each node
on the 2D triangular faces
- beadDisps: 
Output: 
- newElem3DConn: Array containing 3D element ID and the nodal IDs for each
node that makes up the tetrahedral element by row. 
- beadDispsCleaned: cell array containing bead positions and displacements
that are within the 3D mesh. 
- beadElems: Array containing the IDs for the elements each bead is
localized within. 
%}
    elem3DConn = [data.Elements(2).ID,data.Elements(2).NodeIDList];
    %We may need to remap these as well using the newNodeCoords...
    newElem3DConn = elem3DConn;
    for i = 1:length(newElem3DConn(:,1))
        curElem = newElem3DConn(i,:);
        cur1 = curElem(2);
        cur2 = curElem(3);
        cur3 = curElem(4);
        cur4 = curElem(5);
        Conn1 = find(cur1 == newIds(:,1));
        Conn2 = find(cur2 == newIds(:,1));
        Conn3 = find(cur3 == newIds(:,1));
        Conn4 = find(cur4 == newIds(:,1));
        newElem3DConn(i,:) = [i,newIds(Conn1,2),newIds(Conn2,2),newIds(Conn3,2),newIds(Conn4,2)];
    
    end
    

    %beadPos = beadDisps{end}(:,1:3);
    
    % elem3DConn(:,2:5) = elem3DConn(:,2:5) - nodeMin + 1;
    TR=triangulation(newElem3DConn(:,2:5),newNodeCoords(:,2:4)); %Setup a matlab triangulation object to hold the mesh data
    %Do 3D interpolation here but interpolate to grids of 2D nodes? 
    centroids2D = incenter(TR2);


    beadDisps = interpDisps(matches,centroids2D);


    elemIDindex = pointLocation(TR,beadPos(:,1:3)); %Identify which 3D element contains each bead measurement
    beadsInMeshID = [1:size(beadPos,1)]'; %Keep only beads that are found in the mesh
    beadsInMeshID = beadsInMeshID(~isnan(elemIDindex),:);
    beadsInMesh = beadPos(beadsInMeshID,:);
    
    beadDispsCleaned={};
    for kk=1:length(beadDisps)
        curDisp = beadDisps{kk}(:,:);
        %need to scale coordinates but otherwise they are all surface so we
        %should be able to use them. 
        curDisp = curDisp(beadsInMeshID,:);
        beadDispsCleaned{kk}=curDisp;
    end
    
    %% Remap elemIDindex to new system for elem3DConn
    elemIDindex = elemIDindex(~isnan(elemIDindex)); %This should already be in the new element ID system given we defined TR with it. 
    beadElems = elemIDindex; %These are the 3D elements where bead displacements are located
end