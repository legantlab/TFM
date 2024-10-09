%% Open data and sort into lists of nodes and elements
%fileName = 'HyperMeshFileS4_Corrected_CoarseMeshedLoadedSPCsLoaded.inp';
%data = abaqusInpRead(fileName);

%% Sanity check on the meshes
close all
%Plot node coordinates
nodeCoords = data.Nodes.Coordinates;
figure
scatter3(nodeCoords(:,1),nodeCoords(:,2),nodeCoords(:,3))

%Plot mesh

triaElements = data.Elements(1).NodeIDList;  
%Reorder triaElements node IDs to start at 1
% triaElements = triaElements - min(triaElements,[],'all') + 1;
% 
% tr = triangulation(triaElements,nodeCoords);
% figure
% trimesh(tr)

%Setup nodeIDs + coords
nodeConn = [data.Nodes.ID, nodeCoords]; % ID, X, Y, Z
nodeMin = min(nodeConn(:,1),[],'all');
nodeConn(:,1) = nodeConn(:,1) - nodeMin + 1;

%remap triaElements
triaElements = triaElements - nodeMin + 1;
triaIDs = [data.Elements(1).ID, triaElements];
tr = triangulation(triaElements,nodeCoords);
figure
trimesh(tr)


%% Alternative node ID/Coords selections
%We have an issue when the elements don't start at 1, then the index the
%triangulation is looking up won't work. We can try to get around this by
%resetting the index in order (it also skips around annoyingly). 
nodeCoords = [data.Nodes.ID, data.Nodes.Coordinates];
newIds = [nodeCoords(:,1),transpose(linspace(1,length(nodeCoords(:,1)),length(nodeCoords(:,1))))];
%Now we can use this as a lookup table to redo the node coords
%newFaceElemConn = data.Elements(1).NodeIDList;
% newFaceElemConn = data.ElementSets(3).Data;
% elemSetConn = zeros([length(newFaceElemConn),4]);
% elemSetConn(:,1) = newFaceElemConn;
% 
% for i = 1:length(newFaceElemConn)
%     curConn = newFaceElemConn(i);
%     %This is the element ID, now we need to get coordinates from triaIDs
%     curTriaID = find(curConn == triaIDs(:,1));
%     curTriaNodes = triaIDs(curTriaID,:);
%     curConn = curTriaNodes(:,2:4);
% 
%     curConnX = curConn(1);
%     curConnY = curConn(2);
%     curConnZ = curConn(3);
%     ConnX = find(curConnX == newIds(:,1));
%     ConnY = find(curConnY == newIds(:,1));
%     ConnZ = find(curConnZ == newIds(:,1));
%     elemSetConn(i,2:4) = [newIds(ConnX,2),newIds(ConnY,2),newIds(ConnZ,2)];
% 
% end

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

TR2=triangulation(elemSetConn(:,2:4),nodeCoords(:,2:4)); %Setup a matlab triangulation object to hold the mesh data
%Also need to reset nodeCoords
newNodeCoords = [newIds(:,2),nodeCoords(:,2:4)];
[elemCents2D]=elemCentFind2D(newNodeCoords,elemSetConn);
% % Find centers of points
% [elemCents2D]=elemCentFind2D(nodeConn,triaIds);

%% Load displacements
%load("TestAnalysis.mat")

%Because all displacements are near the surface we shouldn't need to do any
%masking
figure
scatter3(matches{end}(:,1)*.156,matches{end}(:,2)*.156,matches{end}(:,3)*1.09)
beadPos = matches{end}(:,1:3);
beadPos(:,1:2) = beadPos(:,1:2)*.156;
beadPos(:,3) = beadPos(:,3)*1.09;

%% Identify mesh elements for each observed bead displacement
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

% elem3DConn(:,2:5) = elem3DConn(:,2:5) - nodeMin + 1;
TR=triangulation(newElem3DConn(:,2:5),newNodeCoords(:,2:4)); %Setup a matlab triangulation object to hold the mesh data

%% 
elemIDindex = pointLocation(TR,beadPos(:,1:3)); %Identify which element contains each bead measurement
beadsInMeshID = [1:size(beadPos,1)]'; %Keep only beads that are found in the mesh
beadsInMeshID = beadsInMeshID(~isnan(elemIDindex),:);
beadsInMesh = beadPos(beadsInMeshID,:);

% 
% 
beadDispsCleaned={};
for kk=1:length(matches)
    curDisp = matches{kk}(:,:);
    curDisp(:,1:2) = curDisp(:,1:2)*.156;
    curDisp(:,4:5) = curDisp(:,4:5)*.156;
    curDisp(:,3) = curDisp(:,3)*1.09;
    curDisp(:,6) = curDisp(:,6)*1.09;
    %need to scale coordinates but otherwise they are all surface so we
    %should be able to use them. 
    beadDispsCleaned{kk}=curDisp;
end
% 
% cleanedBeadCoordinates=scaledBeadCoordinates(nearSurfIndex,:);
%%
figure
scatter3(beadsInMesh(:,1),beadsInMesh(:,2),beadsInMesh(:,3))

%% 
%Remap elemIDindex to new system for elem3DConn
elemIDindex = elemIDindex(~isnan(elemIDindex)); %This should already be in the new element ID system given we defined TR with it. 
beadElems = newElem3DConn(elemIDindex(~isnan(elemIDindex)),1); %Get the element ID for each bead
%Lists the row of elem3DConn corresponding to each bead (this can be different than the element index)

beadElemNodes = newElem3DConn(beadElems,2:5);

repNodes = unique(beadElemNodes); %Get unique beadElemNodes

%% Read in nodal solutions
nodalSolutions = readmatrix('DisplacementsOutputCourseMesh.csv'); %This takes a few minutes at least. Can't do string as it becomes enormous. 

%%

%Get the specific data we need
outputStart = find(~isnan(nodalSolutions),1); %Find the starting index of the matrix
nodeDisps = nodalSolutions(outputStart:end,[1,3,5:7]); 
%%
%Reorganize into the structure we want
[beadu]=beadDispCalcLinear(beadsInMesh,beadElems,newNodeCoords,[1:1:size(newNodeCoords,1)]',nodeDisps,newElem3DConn);

%% 
%Matrix decomposition 
[U,s,V]=csvd(beadu);

%%
%Solve for first time point
temp=beadsInMesh';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U,s,y,'Tikh');

%Compute forces at different time points
reg_corner_timelapse=0.01; %L-curve picks out the wrong corner right now
%force_vector={};
% for kk=1:51
%     kk
% temp=beadDispsCleanedInMesh{kk,1}(:,4:6)';
% y=temp(:);
[x_lambda,rho,eta] = tikhonov(U,s,V,y,reg_corner_timelapse);
force_vector=[elemCents2D,x_lambda(1:3:end),x_lambda(2:3:end),x_lambda(3:3:end)];
% end

%% Plotting!

figure
trisurf(TR2,'FaceColor','#696969','FaceAlpha',0.5,'EdgeAlpha',0.7)
hold on
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector(:,4),force_vector(:,5),-force_vector(:,6),10)
axis equal