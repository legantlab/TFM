%TFM Processing Script
close all
clear
clc
%% Open data and sort into lists of nodes and elements
fileName = 'U:\Max\2023_01_18_CellsInChannels\S1\MeshedHMFileMeshedH5solverDeck.inp';
data = abaqusInpRead(fileName);

%% Reorganize geometry and inputs
nodeCoords = [data.Nodes.ID, data.Nodes.Coordinates];
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
[elemCents2D]=elemCentFind2D(newNodeCoords,elemSetConn);

%% Load displacements
load("U:\Max\2023_01_18_CellsInChannels\S2\S2Test.mat")

%% Optional step to strip out SDS derived substrate expansion

allCents=cell2mat(matches);
allDisps=cell2mat(displacements);
[numBeads,~]=size(allCents);
CENTROIDSR=[[1:numBeads]',allCents(:,4:6)];


%Use one out of every 10 points to save time
polymodelx = polyfitn(CENTROIDSR(1:10:end,2:4),allDisps(1:10:end,1),11);
polymodely = polyfitn(CENTROIDSR(1:10:end,2:4),allDisps(1:10:end,2),[5]);
polymodelz = polyfitn(CENTROIDSR(1:10:end,2:4),allDisps(1:10:end,3),[3]);

dispsC={};
for ii=1:length(matches)
    
    [numBeads,~]=size(matches{ii}(:,1));
    CENTROIDSR=[[1:numBeads]',matches{ii}(:,4:6)];
    dispsC{ii}(:,1)=polyvaln(polymodelx,CENTROIDSR(:,2:4));
    dispsC{ii}(:,2)=polyvaln(polymodely,CENTROIDSR(:,2:4));
    dispsC{ii}(:,3)=polyvaln(polymodelz,CENTROIDSR(:,2:4));
    displacementsR_dedrift{ii}=[displacements{ii}(:,1)-dispsC{ii}(:,1),...
        displacements{ii}(:,2)-dispsC{ii}(:,2),displacements{ii}(:,3)-dispsC{ii}(:,3)];
end

%% Plot dedrift
scale = 1;
[numBeads,~]=size(matches{1}(:,1));
CENTROIDSR=[[1:numBeads]',matches{1}(:,4:6)];
figure
quiver3(CENTROIDSR(:,2)*.156,CENTROIDSR(:,3)*.156,CENTROIDSR(:,4)*1.09,...
    displacementsR_dedrift{1}(:,1)*scale,displacementsR_dedrift{1}(:,2)*scale,...
    displacementsR_dedrift{1}(:,3)*scale,1)
figure
quiver3(matches{1}(:,4),matches{1}(:,5),matches{1}(:,6),displacements{1}(:,1),...
    displacements{1}(:,2),displacements{1}(:,3),1)

%% Replace matches displacements
%If the displacements look better, then we can set the original matrices to
%the dedrifted ones and carry on. --MH 2023/02/14
for j = 1:length(matches)
    curDriftCorrectedDisp = displacementsR_dedrift{j};
    matches{j}(:,1) = matches{j}(:,1) + curDriftCorrectedDisp(:,1);
    matches{j}(:,2) = matches{j}(:,2) + curDriftCorrectedDisp(:,2);
    matches{j}(:,3) = matches{j}(:,3) + curDriftCorrectedDisp(:,3);

end

%%
beadPos = matches{end}(:,1:3);
beadPos(:,1:2) = beadPos(:,1:2)*.156;
beadPos(:,3) = beadPos(:,3)*1.09;

beadDisps = {};

for ii=1:length(matches)
    ii
    curDisp = matches{ii}(:,4:6) - matches{ii}(:,1:3);
    sI1 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,1));
    sI2 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,2));
    sI3 = scatteredInterpolant(matches{ii}(:,1),matches{ii}(:,2),matches{ii}(:,3),curDisp(:,3));
    beadDisps{ii,1}=[matches{end,1}(:,1:3),sI1(matches{end,1}(:,1:3)),sI2(matches{end,1}(:,1:3)),sI3(matches{end,1}(:,1:3))];
end

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
elemIDindex = pointLocation(TR,beadPos(:,1:3)); %Identify which element contains each bead measurement
beadsInMeshID = [1:size(beadPos,1)]'; %Keep only beads that are found in the mesh
beadsInMeshID = beadsInMeshID(~isnan(elemIDindex),:);
beadsInMesh = beadPos(beadsInMeshID,:);

beadDispsCleaned={};
for kk=1:length(matches)
    curDisp = beadDisps{kk}(:,:);
    curDisp(:,1:2) = curDisp(:,1:2)*.156;
    curDisp(:,4:5) = curDisp(:,4:5)*.156;
    curDisp(:,3) = curDisp(:,3)*1.09;
    curDisp(:,6) = curDisp(:,6)*1.09;
    %need to scale coordinates but otherwise they are all surface so we
    %should be able to use them. 
    curDisp = curDisp(beadsInMeshID,:);
    beadDispsCleaned{kk}=curDisp;
end

%% Remap elemIDindex to new system for elem3DConn
elemIDindex = elemIDindex(~isnan(elemIDindex)); %This should already be in the new element ID system given we defined TR with it. 
beadElems = newElem3DConn(elemIDindex(~isnan(elemIDindex)),1); %Get the element ID for each bead

%% Read in nodal solutions
%nodalSolutionsData = h5read('U:\Max\2023_01_18_CellsInChannels\S1\MeshedHMFileMeshedH5.h5','/OPTISTRUCT/RESULT/NODAL/DISPLACEMENT'); %This takes a few minutes at least. Can't do string as it becomes enormous. 
%Get the specific data we need
%outputStart = find(~isnan(nodalSolutions),1); %Find the starting index of the matrix
nodeDisps = [nodalSolutionsData.ID,zeros(length(nodalSolutionsData.X),1),nodalSolutionsData.X,nodalSolutionsData.Y,nodalSolutionsData.Z]; 
nodeDisps = double(nodeDisps);

%add subcase info with repelem
subcaseIDs = 1:17607;
nodeDisps(:,2) = repelem(subcaseIDs,17938);

%% Compute SVD
%Reorganize into the structure we want
[beadu]=beadDispCalcLinear(beadDispsCleaned{1}(:,1:3),beadElems,newNodeCoords,[1:1:size(newNodeCoords,1)]',nodeDisps,newElem3DConn);

%Matrix decomposition 
[U,s,V]=csvd(beadu);

%Solve for first time point
temp=beadDispsCleaned{1}(:,1:3)';
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U,s,y,'Tikh');

%Compute forces at different time points
reg_corner_timelapse=0.0001; %L-curve picks out the wrong corner right now

%% Compute tractions

force_vector={};
for kk=1:length(beadDispsCleaned)
    kk
    temp=beadDispsCleaned{kk}(:,4:6)';
    y=temp(:);
    reg_corner_timelapse=5; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov(U,s,V,y,reg_corner_timelapse);
    force_vector{kk}=[elemCents2D,-1.*x_lambda(1:3:end),-1.*x_lambda(2:3:end),-1.*x_lambda(3:3:end)];
end

%Save force_vector and elemCents2D for later plotting/to return to this
%point...
date = string(datetime('today','Format','yyyy-MM-dd'));
figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4),force_vector{1}(:,5),1)
figure
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),force_vector{1}(:,6),10)
axis equal
%save(strcat('U:\Max\2023_01_18_CellsInChannels\S2\', date,'_ComputedTractionsVarsNice.mat'),'force_vector','elemCents2D','U','s','V','beadDispsCleaned')
%% Plotting!
maxF=0;
for kk=1:length(beadDispsCleaned)
    maxF=max(maxF,max(sqrt(sum(force_vector{kk}(:,4:6).^2,2))));
end
maxF=maxF*.6;

close all

figure
trisurf(TR2,'FaceColor','#696969','FaceAlpha',0.5,'EdgeAlpha',0.7)
hold on
quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4),force_vector{1}(:,5),force_vector{1}(:,6),maxF)
axis equal

figure
quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4),force_vector{1}(:,5),maxF)
% for kk=1:length(beadDispsCleaned)
% trisurf(TR2,'FaceColor','#696969','FaceAlpha',0.5,'EdgeAlpha',0.7)
% hold on
% quiverC3D(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{kk,1}(:,4)*7,force_vector{kk,1}(:,5)*7,-force_vector{kk,1}(:,6)*7,'scale',0,'LineWidth',1,'MaxColorScale',maxF*7)
% axis equal
% view(180,90)
% saveas(gcf,['Y:\Max\2022_04_07_MatlabFEMCode\MatlabHyperMeshCode\ForceImages\topDown_',num2str(kk,'%02.f'),'.png']);
% view(180,15)
% saveas(gcf,['Y:\Max\2022_04_07_MatlabFEMCode\MatlabHyperMeshCode\ForceImages\sideOn_',num2str(kk,'%02.f'),'.png']);
% hold off
% end
for i = 1:length(force_vector)
    %curCellImage = imread
    figure
    quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{i}(:,4),force_vector{i}(:,5),0)
    set(gca,'visible','off')
    ax= gca;
    exportgraphics(ax,['U:\Max\2023_01_18_CellsInChannels\S2\OutputForceMaps\', 'TractionForces_',num2str(i),'.png'],'Resolution',600)
    %saveas(gcf,['U:\Max\2023_01_18_CellsInChannels\S4\TestFrame2D\OutputImages\', 'TractionForces_',num2str(i),'.png']);

end


