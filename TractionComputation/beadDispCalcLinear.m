function [beadu]=beadDispCalcLinear(beadCoords,beadElem,nodeCoords,repNodeList,repNodeDisp,elemList,allNodesList)
% beadCoords = an Nx3 matrix listing the coordinates of the beads in the reference configuration
% beadElem = an Nx1 vector listing the tetrahedral element that contains each bead
% nodeDisp = a matrix containing the computed nodal displacements
% elemList = element connectivity data where the first column is the element ID
% allNodesList = nx4 matrix listing the nodeID and coordinates for all nodes in model
% nodeCoords = nx4 list of nodes and coordinates where forces were applied
% Generate node index

    test = [1:1:size(allNodesList,1)]';
    nodeIndex=zeros(max(test),1);
    for i=1:length(test)
        nodeIndex(test(i))=i;
    end

numBeads=length(beadCoords);
numNodes=length(repNodeList);
numSteps=length(repNodeDisp)/numNodes;

% Initialize node coordinate and displacement matrices
    N=zeros(3,4);
    U=zeros(3*numSteps,4);
    beadu=zeros(3*numBeads,numSteps);
    beadrowind=1;
for i=1:numBeads
    %i/numBeads
    N=[allNodesList(elemList(beadElem(i),2),2:4)',allNodesList(elemList(beadElem(i),3),2:4)',...
        allNodesList(elemList(beadElem(i),4),2:4)',allNodesList(elemList(beadElem(i),5),2:4)'];
    rowind=1;
    for j=1:numSteps
        U(rowind:rowind+2,1:4)=[repNodeDisp(nodeIndex(elemList(beadElem(i),2))+(j-1)*numNodes,3:5)',...
            repNodeDisp(nodeIndex(elemList(beadElem(i),3))+(j-1)*numNodes,3:5)',...
            repNodeDisp(nodeIndex(elemList(beadElem(i),4))+(j-1)*numNodes,3:5)',...
            repNodeDisp(nodeIndex(elemList(beadElem(i),5))+(j-1)*numNodes,3:5)'];
        rowind=rowind+3;
    end
    beadu(beadrowind:beadrowind+2,:)=reshape(dispCalcLin(N,U,beadCoords(i,:)),[3,numSteps]);
    beadrowind=beadrowind+3;
end
end