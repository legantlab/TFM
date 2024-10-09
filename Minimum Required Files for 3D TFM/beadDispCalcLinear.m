function [beadu]=beadDispCalcLinear(beadCoords,beadElem,nodeCoords,repNodeList,repNodeDisp,elemList)
% beadCoords = a matrix listing the coordinates of the beads in the reference configuration
% beadElem = a vector listing the element that contains each bead
% nodeDisp = a matrix containing the computed nodal displacements
% elemList = element connectivity data
% steps = an integer listing the number of computed solutions
    
% Generate node index
    nodeIndex=zeros(max(repNodeList),1);
    for i=1:length(repNodeList)
        nodeIndex(repNodeList(i))=i;
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
    i
    N=[nodeCoords(elemList(beadElem(i),2),2:4)',nodeCoords(elemList(beadElem(i),3),2:4)',nodeCoords(elemList(beadElem(i),4),2:4)',nodeCoords(elemList(beadElem(i),5),2:4)'];
    rowind=1;
    for j=1:numSteps
        U(rowind:rowind+2,1:4)=[repNodeDisp(nodeIndex(elemList(beadElem(i),2))+(j-1)*numNodes,2:4)',repNodeDisp(nodeIndex(elemList(beadElem(i),3))+(j-1)*numNodes,2:4)',repNodeDisp(nodeIndex(elemList(beadElem(i),4))+(j-1)*numNodes,2:4)',repNodeDisp(nodeIndex(elemList(beadElem(i),5))+(j-1)*numNodes,2:4)'];
        rowind=rowind+3;
    end
    beadu(beadrowind:beadrowind+2,:)=reshape(dispCalcLin(N,U,beadCoords(i,:)),[3,numSteps]);
    beadrowind=beadrowind+3;
end
end