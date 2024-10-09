%Generate STL mesh from inputted Tif image/bead localizations. 
testPoints = matches{23}(:,4:6);
%%Add 3D volume random points to test points volume
% newPointX = min(testPoints(:,1)) + (max(testPoints(:,1) - min(testPoints(:,1)))*rand(1000,1));
% newPointY = min(testPoints(:,2)) + (max(testPoints(:,2) - min(testPoints(:,2)))*rand(1000,1));
% newPointZ = min(testPoints(:,3)) + (max(testPoints(:,3) - min(testPoints(:,3)))*rand(1000,1));
% newPoints = [newPointX,newPointY,newPointZ];
% testPoints = [testPoints; newPoints];
%% Covert pixels to actual units
testPoints(:,1:2) = testPoints(:,1:2).*0.1844;
testPoints(:,3) = testPoints(:,3).*0.8;
%% Mask to remove features above the channel structure

maskPoints = testPoints(:,3) > 2;
testPoints = testPoints.*maskPoints;
testPoints( all(~testPoints,2), : ) = [];

%% Add points along boundary to make stl smoother surface
testPoints = [testPoints; [36, 153, 32]; [78, 0, 32]];


%%
% shp = alphaShape(testPoints(:,1),testPoints(:,2),testPoints(:,3));
% 
% [bf, P] = boundaryFacets(shp);
% 
% stlwrite(triangulation(bf,P),'testSTL.stl')

%% 
figure
xlin = linspace(min(testPoints(:,1)), max(testPoints(:,1)), 200);
ylin = linspace(min(testPoints(:,2)), max(testPoints(:,2)), 200);
[X,Y] = meshgrid(xlin, ylin);
% Z = griddata(x,y,z,X,Y,'natural');
% Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(testPoints(:,1),testPoints(:,2),testPoints(:,3),X,Y,'natural'); %Once we have this we can use the surf2stl 
Z(isnan(Z)) = -20;
mesh(X,Y,Z)
figure
test = surf2solid(X,Y,Z,'Elevation',-20);


% testMesh = mesh(X,Y,Z);
% axis tight; hold on
% plot3(testPoints(:,1),testPoints(:,2),testPoints(:,3),'.','MarkerSize',15)