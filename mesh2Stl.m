%Generate STL mesh from inputted Tif image/bead localizations. 
testPoints = matches{1}(:,4:6);
%%Add 3D volume random points to test points volume
% newPointX = min(testPoints(:,1)) + (max(testPoints(:,1) - min(testPoints(:,1)))*rand(1000,1));
% newPointY = min(testPoints(:,2)) + (max(testPoints(:,2) - min(testPoints(:,2)))*rand(1000,1));
% newPointZ = min(testPoints(:,3)) + (max(testPoints(:,3) - min(testPoints(:,3)))*rand(1000,1));
% newPoints = [newPointX,newPointY,newPointZ];
% testPoints = [testPoints; newPoints];
%% Covert pixels to actual units
%load('U:\Max\2023_03_28_IA32ChannelsTryp\S5\S5.mat')
% testPoints(:,1:2) = testPoints(:,1:2).*0.1844;
% testPoints(:,3) = testPoints(:,3).*0.8;
testPoints = matches{1}(:,4:6);
%Visualize
xCor = 0;
yCor = 0;
zCor = 0;
rotx = -0;
rotz = 0;

for i = 1:length(matches)
    curMatches = matches{i};
    curMatches(:,1:3) = rotate_3D(curMatches(:,1:3)', 'x', rotx)';
    curMatches(:,4:6) = rotate_3D(curMatches(:,4:6)', 'x', rotx)';

    curMatches(:,1:3) = rotate_3D(curMatches(:,1:3)', 'z', rotz)';
    curMatches(:,4:6) = rotate_3D(curMatches(:,4:6)', 'z', rotz)';

    curMatches(:,1) = curMatches(:,1) + xCor;
    curMatches(:,2) = curMatches(:,2) + yCor;
    curMatches(:,3) = curMatches(:,3) + zCor;
    curMatches(:,4) = curMatches(:,4) + xCor;
    curMatches(:,5) = curMatches(:,5) + yCor;
    curMatches(:,6) = curMatches(:,6) + zCor;
    matches{i} = curMatches;
end
% figure
% scatter3(matches{1}(:,1) + xCor, matches{1}(:,2) + yCor,matches{1}(:,3) + zCor,'b')

figure
scatter3(matches{1}(:,1), matches{1}(:,2),matches{1}(:,3),'b')

testPoints = matches{1}(:,4:6);

%% Mask to remove features above the channel structure

maskPoints = testPoints(:,3) < 3.5;
testPoints = testPoints.*maskPoints;
testPoints( all(~testPoints,2), : ) = [];



aboveTest = testPoints;
aboveTest(:,2) = aboveTest(:,2) + max(aboveTest(:,2));
belowTest = testPoints;
belowTest(:,2) = belowTest(:,2) + min(belowTest(:,2));

testPoints = [testPoints; aboveTest; belowTest];

figure
scatter3(testPoints(:,1), testPoints(:,2), testPoints(:,3))

%% For squiggles
%translate points to above and below in y to maintain squiggle shape

%% Add points along boundary to make stl smoother surface
testPoints = [testPoints; [-150, -150, 3.5]; [150, 150, 3.5]];%;[-9,-133,-38];[15,133,-38]];

%% Add points on extremes to build bigger area around target mesh
[minX, maxX] = bounds(testPoints(:,1));
[minY, maxY] = bounds(testPoints(:,2));
[minZ, maxZ] = bounds(testPoints(:,3));

%Round each for linspace
minX = round(minX);
minY = round(minY);
minZ = round(minZ);
maxX = round(maxX);
maxY = round(maxY);
maxZ = round(maxZ);

%Add points 50 microns in distance x,y,z from min/max point
boundPoints = [linspace([minX-50, minY-50, minZ], [minX-50,maxY+50,minZ],100); ...
    [maxX+50,minY-50,minZ]; [maxX+50,maxY+50,minZ];...
    [minX-50, minY-50, maxZ]; [minX-50,maxY+50,maxZ]; ...
    [maxX+50,minY-50,maxZ]; [maxX+50,maxY+50,maxZ]];

leftbot = [linspace(minX-50,minX-50,100);linspace(minY-50,maxY+50,100);repelem(minZ,100)]';
rightbot = [linspace(maxX+50,maxX+50,100);linspace(minY-50,maxY+50,100);repelem(minZ,100)]';
lefttop = [linspace(minX-50,minX-50,100);linspace(minY-50,maxY+50,100);repelem(maxZ,100)]';
righttop = [linspace(maxX+50,maxX+50,100);linspace(minY-50,maxY+50,100);repelem(maxZ,100)]';

newPoints = [testPoints; leftbot; rightbot; lefttop; righttop];

%% Create random points on surface defined by bounds of surface
[minX, maxX] = bounds(testPoints(:,1));
[minY, maxY] = bounds(testPoints(:,2));
[minZ, maxZ] = bounds(testPoints(:,3));
minY = minY + 25;
maxY = maxY - 25;
areaAdd = 50;
maxZ = maxZ;

leftx = (minX-areaAdd) + (minX-(minX-areaAdd)).*rand(10000,1);
lefty = minY + (maxY-minY).*rand(10000,1);

rightx = maxX + ((maxX+areaAdd)-maxX).*rand(10000,1);
righty = minY + (maxY-minY).*rand(10000,1);

topx = (minX-areaAdd) + (maxX+areaAdd-(minX-areaAdd)).*rand(10000,1);
topy = (maxY) + (maxY+areaAdd-(maxY)).*rand(10000,1);

botx = (minX-areaAdd) + (maxX+areaAdd-(minX-areaAdd)).*rand(10000,1);
boty = (minY-areaAdd) + (minY-(minY-areaAdd)).*rand(10000,1);

%Identify coordiantes of channel and remove xy coords that are in that
%space
channelxStart = -9;
channelxEnd = 9;
topx = topx(topx<= channelxStart | topx>= channelxEnd);
topy = topy(1:length(topx));
botx = botx(botx<=channelxStart | botx>= channelxEnd);
boty = boty(1:length(botx));
newPoints = [testPoints; [leftx,lefty,repelem(maxZ,10000)']; [rightx,righty,repelem(maxZ,10000)'];...
    [topx,topy,repelem(maxZ,length(topx))']; [botx,boty,repelem(maxZ,length(botx))']];

%Now add elements at bottom of channel in the channelx ranges
channelbottomTopx = (channelxStart) + (channelxEnd-channelxStart).*rand(10000,1);
channelbottomTopy = (maxY) + (maxY+areaAdd-(maxY)).*rand(10000,1);

channelxStart = -9 ;
channelxEnd = 9;
channelbottomBotx = (channelxStart) + (channelxEnd-channelxStart).*rand(10000,1);
channelbottomBoty = (minY-areaAdd) + (minY-(minY-areaAdd)).*rand(10000,1);
newPoints = [newPoints; [channelbottomTopx, channelbottomTopy, repelem(minZ, 10000)'];...
    [channelbottomBotx, channelbottomBoty, repelem(minZ+1, 10000)']];

%repelem(32,length(topx))
figure
scatter3(newPoints(:,1), newPoints(:,2), newPoints(:,3))





%%
% shp = alphaShape(testPoints(:,1),testPoints(:,2),testPoints(:,3));
% 
% [bf, P] = boundaryFacets(shp);
% 
% stlwrite(triangulation(bf,P),'testSTL.stl')

%% 
figure
xlin = linspace(min(newPoints(:,1)), max(newPoints(:,1)), 100);
ylin = linspace(min(newPoints(:,2)), max(newPoints(:,2)), 100);
[X,Y] = meshgrid(xlin, ylin);
% Z = griddata(x,y,z,X,Y,'natural');
% Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(newPoints(:,1),newPoints(:,2),newPoints(:,3),X,Y,'linear'); 
Z(isnan(Z)) = -80;
mesh(X,Y,Z)

test = surf2solid(X,Y,Z,'Elevation',-80); %Once we have this we can use the stlwrite 

%stlwrite('5x40Test.stl',test)

% testMesh = mesh(X,Y,Z);
% axis tight; hold on
% plot3(testPoints(:,1),testPoints(:,2),testPoints(:,3),'.','MarkerSize',15)