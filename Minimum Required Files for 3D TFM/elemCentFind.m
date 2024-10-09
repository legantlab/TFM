function [elemCents]=elemCentFind(nodeList,voidElemConn)

%This function computes the centroids of triangular elements defined by
%nodeList and voidElemConn

%Input argument definitions
% nodeList = a NX4 matrix containing the node coordinates of the triangular
%mesh in the following format [index,xcoord,ycoord,zcoord] where N is the
%number of nodes comprising the mesh

% voidElemConn = a MX4 matrix containing the element connectivity for the
% triangular mesh in the following format [index,node1,node2,node3] where M
% is the total number of facets on the cell surface

%Output argument definitions
% elemCents = a Mx3 matrix containing the centroids of each facet on the cell
% surface in the following format [xcoord,ycoord,zcoord] where M is the
% number of facets on the cell surface

elemCents=zeros(length(voidElemConn),3); %Initialize output arguments

%Compute element centroids
for i=1:length(voidElemConn)
    N1=nodeList(voidElemConn(i,2),2:4);
    N2=nodeList(voidElemConn(i,3),2:4);
    N3=nodeList(voidElemConn(i,4),2:4);
    elemCents(i,:)=[(N1(1)+N2(1)+N3(1))/3,(N1(2)+N2(2)+N3(2))/3,(N1(3)+N2(3)+N3(3))/3];
end