function [elemCents]=elemCentFind3D(nodeList,elemList)

%This function computes the centroids of triangular elements defined by
%nodeList and voidElemConn

%Input argument definitions
% nodeList = a NX4 matrix containing the node coordinates of the triangular
%mesh in the following format [index,xcoord,ycoord,zcoord] where N is the
%number of nodes comprising the mesh

% voidElemConn = a MX5 matrix containing the element connectivity for the
% tetrahedral mesh in the following format [index,node1,node2,node3,node4] where M
% is the total number of facets on the cell surface

%Output argument definitions
% elemCents = a Mx3 matrix containing the centroids of each facet on the cell
% surface in the following format [xcoord,ycoord,zcoord] where M is the
% number of facets on the cell surface

elemCents=zeros(length(elemList),3); %Initialize output arguments

%Compute element centroids
for i=1:length(elemList)
    N1=nodeList(elemList(i,2),2:4);
    N2=nodeList(elemList(i,3),2:4);
    N3=nodeList(elemList(i,4),2:4);
    N4=nodeList(elemList(i,5),2:4);
    elemCents(i,:)=[(N1(1)+N2(1)+N3(1)+N4(1))/4,(N1(2)+N2(2)+N3(2)+N4(2))/4,(N1(3)+N2(3)+N3(3)+N4(3))/4];
end