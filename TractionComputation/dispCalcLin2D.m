function [newNodalSolutions] = dispCalcLin2D(elemCents, nodeIds, nodalData)
%dispCalcLin2D Uses 2D Triangular shape function to calculate displacement
%at centroid of listed nodes. 
%   Uses 2D shape functions in 3D space to compute the displacements of the
%   elements centered on the centroid rather than the nodal displacements. 
%{
Inputs:
- elemCents: array containing element IDs and their node Ids for the 3
nodes that make up each element
- nodeIds: array containing node Ids and the xyz positions for each node
- nodalData: nodalsolution data 
Outputs:
- newNodalSolutions: Nodal solutions array where the rows are now computed
as the centroids of each element
%}

%Loop through each subcase (column) of nodalData


x1 = 1;
y1 = 1;
z1 = 1;

x2 = 4;
y2 = 5;
z2 = 6;

x3 = 5;
y3 = 8;
z3 = 9;

x = (x1 + x2 + x3)/3;
y = (y1 + y2 + y3)/3;
z = (z1 + z2 + z3)/3;

b = (y3-y + (((y1-y3)*x) - (y1-y3)*x3)/(x1-x3))/(((y1-y3)*(x2-x3))/(x1-x3) - (y2-y3));
a = (x-x3 - (x2-x3)*b)/(x1-x3);

xtest = (a*(x1) + b*(x2) + (1-a-b)*(x3));
ytest = (a*(y1) + b*(y2) + (1-a-b)*(y3));
ztest = (a*(z1) + b*(z2) + (1-a-b)*(z3));



end