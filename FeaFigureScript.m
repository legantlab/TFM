%Script to plot Example FEA and Green's matrix results for figures. 
close all
%Assumes you have a traction dataset loaded. 

figure
%colormap("jet")
trimesh(TR2,'FaceColor','#808080','FaceAlpha',1,'EdgeColor','black')

hold on 

SampleCase = 7; 
elemNum = ceil(SampleCase/3);
elemPoints = TR2.ConnectivityList(elemNum,:);
x1 = TR2.Points(elemPoints(1),:);
x2 = TR2.Points(elemPoints(2),:);
x3 = TR2.Points(elemPoints(3),:);
SampleElemCoords = [x1; x2; x3];
elemCentroid = [sum(SampleElemCoords(:,1))/3 , sum(SampleElemCoords(:,2))/3, sum(SampleElemCoords(:,3))/3];
%Compute sample green's response
response = beadu(:,SampleCase);
responseCoords = [response(1:3:end), response(2:3:end), response(3:3:end)];

%Plot the reponse to the applied load
% quiver3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3)+1,...
%     responseCoords(:,1), responseCoords(:,2), responseCoords(:,3),'blue','LineWidth',1)

%plot3(SampleElemCoords(:,1), SampleElemCoords(:,2), SampleElemCoords(:,3),'r','LineWidth',3)
%Could also change the face color of that specific element
testConnect = [1 2 3];
DT = triangulation(testConnect,SampleElemCoords);
trimesh(DT,'facecolor','#F80')

%Now plot quiver showing applied pressures
quiver3(elemCentroid(1), elemCentroid(2),elemCentroid(3),...
    10, 0 ,0,'m','LineWidth',1,'ShowArrowHead','on','AutoScale','off','MaxHeadSize',10)
axis equal
axis off
ax = gca;
ax.Clipping = "off";
view(90,90)
zoom(2)
% 
zoom(4)
xlim([-31 -8])
ylim([-6 17])
set(gcf,'position',[560   420   420   420])


%zlim([-12.5730 -9.4730])
% set(gca, 'Color', 'none')
% plot2svg('test.svg');
% 
% hold off
%% Plot individual Element
figure
testConnect = [1 2 3];
DT = triangulation(testConnect,SampleElemCoords);
elemCentroid = [sum(SampleElemCoords(:,1))/3 , sum(SampleElemCoords(:,2))/3, sum(SampleElemCoords(:,3))/3];

trimesh(TR2)
hold on
quiver3(elemCents2Dnew(:,1), elemCents2Dnew(:,2), elemCents2Dnew(:,3),...
    responseCoords(:,1), responseCoords(:,2), responseCoords(:,3))

trimesh(DT,'facecolor','r')

%Now plot quiver showing applied pressures
quiver3(elemCentroid(1), elemCentroid(2),elemCentroid(3)+1,...
    10, 0 , 0,'g','LineWidth',1,'ShowArrowHead','on','AutoScale','off','MaxHeadSize',10)

% quiver3(elemCentroid(1), elemCentroid(2),elemCentroid(3)+1,...
%     0, 10 , 0,'g','LineWidth',1,'ShowArrowHead','on','AutoScale','off','MaxHeadSize',10)
% 
% quiver3(elemCentroid(1), elemCentroid(2),elemCentroid(3)+1,...
%     0, 0 , -10,'g','LineWidth',1,'ShowArrowHead','on','AutoScale','off','MaxHeadSize',10)

axis equal