%General script for processing microchannel TFM Data - MH 03/2024
outputFold = 'T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06'; %Folder to store output tractions and amira maps. 
plotFlag = 1; %If flag = 1, outputs a series of QA plots to the output folder. 
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots'])
end
%% Open solverdeck using abaqus reader
fileName = 'T:\Max\2023-10-08\Tiffs\F02\TestData\5x40_Raw_SmoothMesh2p5_s10Course_QA_Tetra_Properties_Loads_SPCs_Loaded_SolverDeck.inp';
data = abaqusInpRead(fileName);
%% Reorganize geometry and inputs
[newIds, TR2,newNodeCoords, surfNodes,elemCents2Dnew] = organizeGeometry(data,1,1);
%% Get the area of each element for computing force balance
elemList = TR2.ConnectivityList;
elemCoords = TR2.Points;
elemAreas = zeros(length(elemList),1);
for i = 1:length(elemList)
    curElem = elemList(i,:);
    p1 = elemCoords(curElem(1),:);
    p2 = elemCoords(curElem(2),:);
    p3 = elemCoords(curElem(3),:);
    
    V1 = p2 - p1;
    V2 = p3 - p1;
    elemAreas(i) = 0.5*norm(cross(V1,V2));
end
%% Load displacements
load("T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\2023_10_07_F06_NewTracking.mat")

%Write a function to do all the below. 
for i = 1:length(matches) % Fix pixel scaling. 
    curMatch = matches{i};
    curMatch(:,1:2) = curMatch(:,1:2)*.199;
    curMatch(:,3) = curMatch(:,3)*.8;
    curMatch(:,4:5) = curMatch(:,4:5)*.199;
    curMatch(:,6) = curMatch(:,6)*.8;
    matches{i} = curMatch; %#ok<SAGROW> 
end

regElems = elemCents2Dnew(elemCents2Dnew(:,1)>-40 & elemCents2Dnew(:,1)<40,:);

%Align elemCents2Dnew to matches
PC1 = pointCloud(regElems);
PC2 = pointCloud(matches{1}(:,1:3));

tform = pcregistericp(PC1,PC2);
PC1reg = pctransform(pointCloud(elemCents2Dnew), tform);
elemCents2Dnew = PC1reg.Location;
%Try masking the raw matches displacements by cell mask
Masks = "T:\Max\2023-10-07\5x40_Overnight1\Tiffs\F06\LA_masks";
curDir = char(Masks);
    
tiffs=dir([curDir,'\*.tif']);
tiffs=tiffs(~ismember({tiffs.name},{'.','..'}));
currentAllMasks = cell(length(tiffs),1);
for i = 1:length(tiffs(:,1))
       curMask = [tiffs(i).folder,'\',tiffs(i).name];
       currentAllMasks{i} = loadtiff(curMask);
end

dilation = strel('diamond',100);
maskedCoords = matches;
for j = 1:length(currentAllMasks)-1
    %curMask = rot90(currentAllMasks{j},2);
    curMask = (currentAllMasks{j});
    ImgSize = size(curMask);
    curMask = imdilate(curMask, dilation);
    matchCoords = matches{j}(:,:);
    matchCoords(:,1) = round(matchCoords(:,1)/.199 + ImgSize(2)/2 - 1);
    matchCoords(:,2) = round(matchCoords(:,2)/.199 + ImgSize(1)/2 - 1);

    [ycoords, xcoords] = find(curMask);
    coordinateMask = ismember(matchCoords(:,1:2),[xcoords,ycoords],'rows');
    maskedCoords{j} = matches{j}(coordinateMask,:);
    %rather set the second set of matches of the not coordinate mask to be
    %equal to the first set of corrdinates
    matches{j}(~coordinateMask,4:6) = matches{j}(~coordinateMask,1:3); %#ok<SAGROW> 
    maskedCoords{j} = matches{j}(coordinateMask,:);
end

time = 6;
newDisps = matches{time}(:,4:6) - matches{time}(:,1:3);

if plotFlag == 1
    %Plot elements and matches to ensure they align
    figure
    scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))
    hold on 
    scatter3(matches{time}(:,4), matches{time}(:,5), matches{time}(:,6))
    axis equal
    legend('Elements','Data')
    saveas(gcf,[outputFold, '\', 'QAplots\ModelAlignment.fig'])
    close()
    
    %Plot Masked Coordinates
    figure
    scatter3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3))
    hold on
    scatter3(matches{time}(:,4), matches{time}(:,5), matches{time}(:,6))
    scatter3(maskedCoords{time}(:,1), maskedCoords{time}(:,2),maskedCoords{time}(:,3),'b')
    axis equal
    legend('Elements','Data', 'Masked Data')
    saveas(gcf,[outputFold, '\', 'QAplots\MaskedCoordinates.fig'])
    close()
end
%% Interpolate matches to surface Nodes
beadDisps = interpDisps(matches,elemCents2Dnew);

%Plot inteperpolated displacements versus raw
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots\InterpDisps'])
    for i = 1:length(beadDisps)
        curBeadDisps = beadDisps{i};
        curMatch = matches{i};
        curDisps = displacements{i};
        figure
        quiver(curBeadDisps(:,1), curBeadDisps(:,2), ...
            curBeadDisps(:,4),curBeadDisps(:,5))
        hold on 
        quiver(curMatch(:,1), curMatch(:,2), ...
            curDisps(:,1),curDisps(:,2))
        legend('Interpolated Displacements','Raw Displacements')
        axis equal
        set(gcf,'units','normalized','position',[0 0 1 1])
        saveas(gcf,[outputFold, '\', 'QAplots\InterpDisps\t',sprintf( '%04d', i ),'_MaskedCoordinates.png'])
        close()
    end
end
%% Read nodal solutions and parse
nodalSolutionsFile = 'T:\Max\2023-10-08\Tiffs\F02\TestData\5x40_Raw_SmoothMesh2p5_s10Course_QA_Tetra_Properties_Loads_SPCs_Loaded.h5';
nodeDisps = readNodalSolutions(nodalSolutionsFile,length(data.ElementSets.Data),newIds);
%% Construct Green's matrix
[beadu] = composeGreens(nodeDisps,elemCents2Dnew, TR2);
%% Compute CSVD and Tractions
%Matrix decomposition 
[U2,s2,V2]=csvd(beadu);

%% L Curve
%Solve for first time point
temp=beadDisps{1}(:,4:6)';
%temp(:,~coordinateMask) = 0;
y=temp(:);

%Compute L-curve
[reg_corner,rho,eta,reg_param] = l_curve(U2,s2,y,'Tikh');
if plotFlag == 1
    saveas(gcf,[outputFold, '\', 'QAplots\Lcurve.png'])
    close()
end
%% Plotting
close all
traction_vector=cell(length(beadDisps),1);
reg_corner_timelapse = 10;%reg_corner;%reg_corner*2;
newElemPositions=cell(length(beadDisps),1);
for kk =1:length(beadDisps)
    temp=beadDisps{kk}(:,4:6)';
    y=temp(:);
    %reg_corner_timelapse=reg_corner; %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U2,s2,V2,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
%Compute updated positions for displacements
    computedDisp = beadu*x_lambda;
    computedDispArray = [computedDisp(1:3:end), computedDisp(2:3:end), computedDisp(3:3:end)];
    newElemPositions{kk} = elemCents2Dnew + computedDispArray;
    traction_vector{kk}=[newElemPositions{kk},x_lambda(1:3:end),...
        x_lambda(2:3:end),x_lambda(3:3:end)];
end

date = string(datetime('today','Format','yyyy-MM-dd'));
scale = 20000;
figure
kk = 1; %time to plot
quiver(traction_vector{kk}(:,1),traction_vector{kk}(:,2),traction_vector{kk}(:,4)*scale,traction_vector{kk}(:,5)*scale,1)
axis equal
figure
quiver3(traction_vector{kk}(:,1),traction_vector{kk}(:,2),traction_vector{kk}(:,3),traction_vector{kk}(:,4)...
    *scale,traction_vector{kk}(:,5)*scale,traction_vector{kk}(:,6)*scale,1)
axis equal

%save(strcat('T:\Max\2023-10-08\Tiffs\F00\', date,'_ComputedTractionsAligned.mat'),'force_vector','elemCents2D','U2','s2','V2','correctedbeadDisps')
save(strcat(outputFold, '\', date,'_ComputedTractionsAlignedValues.mat'),'traction_vector','elemCents2Dnew','beadDisps')
%% Save 2D QA Plots of tractions
if plotFlag == 1
    mkdir([outputFold, '\', 'QAplots\Tractions'])
    for i = 1:length(traction_vector)
        figure
        quiver(traction_vector{i}(:,1),traction_vector{i}(:,2),traction_vector{i}(:,4),traction_vector{i}(:,5),1)
        axis equal
        set(gcf,'units','normalized','position',[0 0 1 1])
        saveas(gcf,[outputFold, '\', 'QAplots\Tractions\t',sprintf( '%04d', i ),'_tractions.png'])
        close()
    end
end

%% Output vectors in Amira formatting
%Need to update this to loop in each set of positions and apply translation
%for each time point
mkdir([outputFold, '\',convertStringsToChars(date),'_Forces'])
for i = 1:length(traction_vector)
    i
    translatePos = traction_vector{i}(:,1:3);
    curTranslate = translation{i}; 
    translatePos(:,1) = (translatePos(:,1)) + (curTranslate(1))*.199;
    translatePos(:,2) = (translatePos(:,2)) + (curTranslate(2))*.199;
    translatePos(:,3) = (translatePos(:,3)) + curTranslate(3)*.8;
    text=buildAmiraMesh(1*traction_vector{i}(:,4:6),translatePos);
    fid=fopen([outputFold, '\',convertStringsToChars(date),'_Forces', '\t', sprintf( '%04d', i ), 'Forces.am'],'wt');
    fprintf(fid,text);
    fclose(fid);
end