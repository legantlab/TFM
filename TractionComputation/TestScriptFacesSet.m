% Script to figure out how to seperate face elements we applied forces from
% all tria elements. 

nodeCoords = [data.Nodes.ID, data.Nodes.Coordinates];
newIds = [nodeCoords(:,1),transpose(linspace(1,length(nodeCoords(:,1)),length(nodeCoords(:,1))))];

newFaceElemConn = data.ElementSets(3).Data;

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