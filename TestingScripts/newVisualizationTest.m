%Different visualization scheme using 3D scatter plots. 

%% Load image
I = loadtiff(fileInfo{tVal + 1,1});

%%
tVal = 1;
test_RelaxBeads = matches{tVal}(:,4:6);

test_StrainedBeads = matches{tVal}(:,1:3);

curIz = imtranslate(I,driftStore{tVal}); 
% figure
% [X,Y] = meshgrid(1:size(curIz,2),1:size(curIz,1));
% warp(X,Y,50*ones(size(X)),curIz)

figure
%imshow(curIz(:,:,50), [1500 4000])

scatter3(test_RelaxBeads(:,1),test_RelaxBeads(:,2),test_RelaxBeads(:,3),'b')
hold on
scatter3(test_StrainedBeads(:,1),test_StrainedBeads(:,2),test_StrainedBeads(:,3),'r')

hold off