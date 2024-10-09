%% Plots
close all
% figure
% subplot(1,3,1)
% hist(displacements{1}(:,1))
% subplot(1,3,2)
% hist(displacements{1}(:,2))
% subplot(1,3,3)
% hist(displacements{1}(:,3))
figure
image = loadtiff(fileInfo{1});
%imagesc(max(image,[],3))

image2 = loadtiff(fileInfo{2});
%imagesc(max(image2,[],3))
imshowpair(max(image,[],3),max(image2,[],3))
hold on 
scatter(matches{1}(:,1)/size_x + 1, matches{1}(:,2)/size_y + 1)
quiver(matches{1}(:,1)/size_x + 1, matches{1}(:,2)/size_y + 1, displacementsWithDrift{1}(:,1), displacementsWithDrift{1}(:,2),0)
title('Raw Images Comparison')
%scatter(matches{1}(:,4)/size_x + 1, matches{1}(:,5)/size_y + 1)
%Figure for translated image
figure
imshowpair(max(image,[],3),imtranslate(max(image2,[],3),-1*driftStore{1}(1:2)/size_x))
hold on 
scatter(matches{1}(:,1)/size_x + 1, matches{1}(:,2)/size_y + 1)
quiver(matches{1}(:,1)/size_x + 1, matches{1}(:,2)/size_y + 1, displacements{1}(:,1), displacements{1}(:,2),0)
title('Drift Corrected Comparison')

%Create Histograms of XYZ displacements
figure
h1 = histogram(displacements{1}(:,1));
hold on 
h2 = histogram(displacements{1}(:,2));
h3 = histogram(displacements{1}(:,3));
h1.Normalization = 'probability';
h1.BinWidth = 0.05;
h2.Normalization = 'probability';
h2.BinWidth = 0.05;
h3.Normalization = 'probability';
h3.BinWidth = 0.05;
xlabel('Displacement (microns)')
ylabel('Probability')
legend('X', 'Y', 'Z')

%Now do histogram fitting for each
figure
hxfit = histfit(displacements{1}(:,1));
xfit = fitdist(displacements{1}(:,1),'normal');
hold on 
hyfit = histfit(displacements{1}(:,2));
yfit = fitdist(displacements{1}(:,2),'normal');
hzfit = histfit(displacements{1}(:,3));
zfit = fitdist(displacements{1}(:,3),'normal');

xlabel('Displacement (microns)')
ylabel('Probability')

title('Histogram Fits')

display(['Xfit = ', num2str(xfit.sigma)])
display(['Yfit = ', num2str(yfit.sigma)])
display(['Zfit = ', num2str(zfit.sigma)])
%%
figure
left = 20;
right = 80;
xzImg = squeeze(max(image(:,left:right,:),[],2));
xzImg2 = squeeze(max(image2(:,left:right,:),[],2));
%imagesc(xzImg)
imshowpair(xzImg,xzImg2)
hold on 
%Get relevant matches
matchesRel = matches{1};
matchesRel = matchesRel(matchesRel(:,1)/size_y + 1 >= left-2 & matchesRel(:,1)/size_y + 1 <= right+2,:);
scatter(matchesRel(:,3)/size_z + 1, matchesRel(:,2)/size_x + 1,'r')
scatter(matchesRel(:,6)/size_z + 1, matchesRel(:,5)/size_x + 1,'b')
axis image
%scatter(matchesRel(:,4)/size_y + 1, matchesRel(:,6)/size_z + 1)
% %XZ view
% figure
% INPUT = image;
% ZZ = padarray(INPUT,[1 1],0,'post');
% [XX,YY] = meshgrid((1:size(INPUT,2)+1)-0.5,(1:size(INPUT,1)+1)-0.5);
% surface(XX,0*ZZ,YY,ZZ,'EdgeColor','none','FaceColor','flat');
% view([-50 50]); xlabel('x'); ylabel('y'); zlabel('z'); axis ij; box on; grid on;
% title('X-Z surface'); caxis([min(INPUT(:)),max(INPUT(:))]);