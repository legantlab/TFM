%Second pass at alignment, do XYZ translation first, then rotation
%Script to automatically align two 3D channel point clouds. 

%Assume that two channels are rotated approximately the same, I.E channel
%runs up and down in Y for the image (prerotation is applied) so we can
%compute the YZ profile for top and bottom, then fit to a guassian and
%align using those functions. 

%The goal is to align the swelling data with the experimental data and
%return a new swelling data coordinate system to do swelling correction. 

%Import swelling data
load('T:\Max\2023_07_08_SpinningDiskTesting\Swelling Data\averagedSwellingProfile7x40Final.mat');
%Scale swelling data set by pixel sizes
elemCents2D(:,1) = elemCents2D(:,1)/.199;
elemCents2D(:,2) = elemCents2D(:,2)/.199;
elemCents2D(:,3) = elemCents2D(:,3)/0.8;
elementsDim = [max(elemCents2D(:,1)) - (min(elemCents2D(:,1))), max(elemCents2D(:,2)) - (min(elemCents2D(:,2)))...
    ,max(elemCents2D(:,3)) - (min(elemCents2D(:,3)))]; %X,Y,Z, 
%Import experimental data
I = loadtiff('T:\Max\2023-06-22\Tiffs\GreatMovies\ON2_F18\MoreCropped\Beads\Time_53.ome.tiff');

%% Localize beads in experimental data
bandpass_size = [3,3,5];
ImgDim=size(I); %Y,X,Z
lnoise = [.5 .5 .5];
inputb = [0,0];
bPassParams = {lnoise, bandpass_size, inputb};


[IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3});

threshmulti = 1.5;
threshold = 0.99;
x{1} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold);
x{1}(any(isnan(x{1}),2),:) = [];

figure
scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
hold on 
scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))
title('Raw Localizations and Element Centroids')
%% Perform X,Z Alignment

datafiltUp = x{1}(:,2)>=ImgDim(:,1)*.75;
elemfiltUp = elemCents2D(:,2) >=elementsDim(:,2)*.75;
%Invert for guassian fit. 
dataupFitVal = [x{1}(datafiltUp,1),-1*x{1}(datafiltUp,3)];
elemupFitVal = [elemCents2D(elemfiltUp,1), -1*elemCents2D(elemfiltUp,3)];
%fix z to zero
dataupFitVal(:,2) = dataupFitVal(:,2) - median(dataupFitVal(:,2));
elemupFitVal(:,2) = elemupFitVal(:,2) - min(elemupFitVal(:,2));
%Fit to Gaussian profile
dataupFit = fit(dataupFitVal(:,1),dataupFitVal(:,2), 'gauss1');
elemupFit = fit(elemupFitVal(:,1),elemupFitVal(:,2), 'gauss1');

%plot fits
figure
scatter(dataupFitVal(:,1),dataupFitVal(:,2))
hold on 
scatter(elemupFitVal(:,1), elemupFitVal(:,2))
plot(dataupFit)
plot(elemupFit)
title('Upper Quartile fits')

%Perform alignment using fits
datayupCent = dataupFit.b1;
elemyupCent = elemupFit.b1;
ydiffUp = elemyupCent - datayupCent;


datafiltDown = x{1}(:,2)<=ImgDim(:,1)*.25;
elemfiltDown = elemCents2D(:,2) <=elementsDim(:,2)*.25;

%Invert for guassian fit. 
datadownFitVal = [x{1}(datafiltDown,1),-1*x{1}(datafiltDown,3)];
elemdownFitVal = [elemCents2D(elemfiltDown,1), -1*elemCents2D(elemfiltDown,3)];
%fix z to zero
datadownFitVal(:,2) = datadownFitVal(:,2) - median(datadownFitVal(:,2));
elemdownFitVal(:,2) = elemdownFitVal(:,2) - min(elemdownFitVal(:,2));
%Fit to Gaussian profile
datadownFit = fit(datadownFitVal(:,1),datadownFitVal(:,2), 'gauss1');
elemdownFit = fit(elemdownFitVal(:,1),elemdownFitVal(:,2), 'gauss1');

%Perform alignment using fits
dataydownCent = datadownFit.b1;
elemydownCent = elemdownFit.b1;
ydiffDown = elemydownCent - dataydownCent;

%Plot fits
figure
scatter(datadownFitVal(:,1),datadownFitVal(:,2))
hold on 
scatter(elemdownFitVal(:,1), elemdownFitVal(:,2))
plot(datadownFit)
plot(elemdownFit)
title('Lower Quartile Fits')


% Perform X Alignment
%use average difference to roughly align channels in x. 
avgxshift = (ydiffDown + ydiffUp)/2;

%subtract from swelling dataset
newElements = elemCents2D;
newElements(:,1) = newElements(:,1) - avgxshift;

%plot to confirm
figure
scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
hold on 
scatter3(newElements(:,1), newElements(:,2), newElements(:,3))
title('X Aligned Elements and Data')

%% Perform Z Alignment (sensitive to outliers?)
elemFilt = newElements(:,3) >= max(newElements(:,3)) -1; 
elemzFit = fit(newElements(elemFilt, 2),newElements(elemFilt,3),'poly1');
dataFilt = x{1}(:,3) >= max(x{1}(:,3)) - 10;

zCorr = elemzFit.p2 - median(x{1}(dataFilt,3));
newElements(:,3) = newElements(:,3) - zCorr;
figure
scatter(newElements(:,2), newElements(:,3))
hold on 
scatter(x{1}(:,2), x{1}(:,3))

%newElements(:,3) = newElements(:,3) - max(newElements(:,3));
%x{1}(:,3) = x{1}(:,3) - median(x{1}(:,3)); DONT CHANGE RAW DATA

figure
scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
hold on 
scatter3(newElements(:,1), newElements(:,2), newElements(:,3))
title('XZ Aligned Elements and Data')


%% Perform Y Alignment
%filter from center +/- 20 pixel(~10 microns)
elementsDim = [max(newElements(:,1)) - abs(min(newElements(:,1))), max(newElements(:,2)) - abs(min(newElements(:,2)))...
    ,max(newElements(:,3)) - abs(min(newElements(:,3)))]; %X,Y,Z, 

datayFilt = x{1}(:,1) >= ImgDim(:,2)/2 - 30 & x{1}(:,1) <= ImgDim(:,2)/2 + 30;
elemyFilt = newElements(:,1) >= elementsDim(:,1)/2 - 30 & newElements(:,1) <= elementsDim(:,1)/2 + 30;
%fit y components from filter with guassian (histfit?)


datayFit = fitdist(x{1}(datayFilt,2),'Normal');
elemyFit = fitdist(newElements(elemyFilt,2),'Normal');

% figure
% scatter(x{1}(datayFilt,2), -1*x{1}(datayFilt,3))
% hold on 
% xline(datayFit.mu)
% scatter(newElements(elemyFilt,2), -1*newElements(elemyFilt,3))
% xline(elemyFit.mu)

avgyshift = datayFit.mu - elemyFit.mu;

transElems = newElements;
transElems(:,2) = transElems(:,2) + avgyshift;

figure
scatter3(x{1}(:,1),x{1}(:,2), x{1}(:,3))
hold on 
scatter3(transElems(:,1), transElems(:,2), transElems(:,3))
title('XYZ aligned Elements and Data')
legend('Raw Data', 'Aligned Element Centroids')

%% Now do rotational Alignments
%Before this, we may want to move the origin to 0,0,0 given that the
%rotations will otherwise be offset (in theory). 
%Rotate X by fitting top surface
dataTopfilt = x{1}(:,3) >= max(x{1}(:,3)) -5;
elemTopfilt = transElems(:,3) >= max(transElems(:,3)) - 1;

%fit to line
dataTopfit = fit(x{1}(dataTopfilt,2), x{1}(dataTopfilt,3),'poly1');
elemTopfit = fit(transElems(elemTopfilt,2), transElems(elemTopfilt,3),'poly1');

figure
scatter(x{1}(dataTopfilt,2), x{1}(dataTopfilt,3))
hold on 
scatter(transElems(elemTopfilt,2), transElems(elemTopfilt,3))
plot(dataTopfit)
plot(elemTopfit)
title('Top Surface Fits')

xRotdata = 0; %atan(dataTopfit.p1);
xRotelem = atan(elemTopfit.p1) - atan(dataTopfit.p1);

rotData = rotate_3D(x{1}','x',-1*xRotdata)';
rotElems = rotate_3D(transElems','x',-1*xRotelem)';
rotDisps = rotate_3D(averageDisplacement','x',-1*xRotelem)';

figure
scatter3(rotData(:,1), rotData(:,2), rotData(:,3))
hold on 
scatter3(rotElems(:,1), rotElems(:,2), rotElems(:,3))
title('XYZ Translation, X rotation Corrected')

%Rotate Z by computing tilt from yz profiles. 
% This isn't working for some reason as we might expect. We may need to
% recompute the yz profiles!

%Plot yz profiles
datayrange = max(rotData(:,2)) - abs(min(rotData(:,2)));
elemyrange = max(rotElems(:,2)) - abs(min(rotElems(:,2)));

%Data
dataupFilt = rotData(:,2) >= datayrange*.75;
datadownFilt = rotData(:,2) <= datayrange*.25;
dataupVal = [rotData(dataupFilt,1), -1*rotData(dataupFilt,3) - median(-1*rotData(dataupFilt,3))];
datadownVal = [rotData(datadownFilt,1), -1*rotData(datadownFilt,3) - median(-1*rotData(dataupFilt,3))];

dataupFit = fit(dataupVal(:,1), dataupVal(:,2),'gauss1');
datadownFit = fit(datadownVal(:,1), datadownVal(:,2),'gauss1');

figure
scatter(dataupVal(:,1),dataupVal(:,2))
hold on 
scatter(datadownVal(:,1),datadownVal(:,2))
plot(dataupFit)
plot(datadownFit)
title('Data Y fits')

%Elements
elemupFilt = rotElems(:,2) >= datayrange*.75;
elemdownFilt = rotElems(:,2) <= datayrange*.25;
elemupVal = [rotElems(elemupFilt,1), -1*rotElems(elemupFilt,3) - min(-1*rotElems(elemupFilt,3))];
elemdownVal = [rotElems(elemdownFilt,1), -1*rotElems(elemdownFilt,3) - min(-1*rotElems(elemupFilt,3))];

elemupFit = fit(elemupVal(:,1), elemupVal(:,2),'gauss1');
elemdownFit = fit(elemdownVal(:,1), elemdownVal(:,2),'gauss1');

figure
scatter(elemupVal(:,1),elemupVal(:,2))
hold on 
scatter(elemdownVal(:,1),elemdownVal(:,2))
plot(elemupFit)
plot(elemdownFit)
title('Elem Y fits')

figure
scatter(dataupVal(:,1),dataupVal(:,2))
hold on
scatter(elemupVal(:,1),elemupVal(:,2))
plot(dataupFit)
plot(elemupFit)
title('Up YZ ranges')

figure
scatter(datadownVal(:,1),datadownVal(:,2))
hold on
scatter(elemdownVal(:,1),elemdownVal(:,2))
plot(datadownFit)
plot(elemdownFit)
title('Down YZ ranges')


% Z Rotation Correction
avgzrot = (abs((elemupFit.b1 - dataupFit.b1)) + abs((elemdownFit.b1 - datadownFit.b1)))/2;

distxElem = elemyrange; %(elemyrange*.75 - elemyrange*.25);
rotzElem = atan(distxElem/avgzrot) - pi/2; 
finalElements = rotate_3D(rotElems','z',-1*rotzElem)';
finalDisplacements = rotate_3D(rotDisps','z',-1*rotzElem)';

distxData = (datayrange*.75 - elemyrange*.25);
rotzData = 0; %atan(distxData/abs((elemdownFit.b1 - datadownFit.b1))) - pi/2;
finalExpData = rotate_3D(rotData','z',-1*rotzData)';

%% Replot xy profiles to make sure we have done the correction
%Upper profiles
finalExpupVal= [finalExpData(dataupFilt,1),-1*finalExpData(dataupFilt,3) - median(-1*finalExpData(dataupFilt,3))];
finalElemupVal = [finalElements(elemupFilt,1),-1*finalElements(elemupFilt,3) - min(-1*finalElements(elemupFilt,3))];
%finalExpdownVal= [finalExpData(datadownFilt,1),-1*finalExpData(datadownFilt,3)];

finaldataupFit = fit(finalExpupVal(:,1), finalExpupVal(:,2),'gauss1');
finalelemupFit = fit(finalElemupVal(:,1), finalElemupVal(:,2),'gauss1');
%finaldatadownFit = fit(finalExpdownVal(:,1), finalExpdownVal(:,2),'gauss1');

figure
scatter(finalExpupVal(:,1), finalExpupVal(:,2))
hold on 
scatter(finalElemupVal(:,1), finalElemupVal(:,2))
plot(finaldataupFit)
plot(finalelemupFit)
title('Upper Profiles')

%Lower profiles
finalExpdownVal= [finalExpData(datadownFilt,1),-1*finalExpData(datadownFilt,3) - median(-1*finalExpData(datadownFilt,3))];
finalElemdownVal = [finalElements(elemdownFilt,1),-1*finalElements(elemdownFilt,3) - min(-1*finalElements(elemdownFilt,3))];
%finalExpdownVal= [finalExpData(datadownFilt,1),-1*finalExpData(datadownFilt,3)];

finaldatadownFit = fit(finalExpdownVal(:,1), finalExpdownVal(:,2),'gauss1');
finalelemdownFit = fit(finalElemdownVal(:,1), finalElemdownVal(:,2),'gauss1');
%finaldatadownFit = fit(finalExpdownVal(:,1), finalExpdownVal(:,2),'gauss1');

figure
scatter(finalExpdownVal(:,1), finalExpdownVal(:,2))
hold on 
scatter(finalElemdownVal(:,1), finalElemdownVal(:,2))
plot(finaldatadownFit)
plot(finalelemdownFit)
title('Lower Profiles')



%% Final Plots

figure
scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
hold on 
scatter3(finalElements(:,1), finalElements(:,2), finalElements(:,3))
quiver3(finalElements(:,1), finalElements(:,2), finalElements(:,3),finalDisplacements(:,1), finalDisplacements(:,2),...
    finalDisplacements(:,3),1)
title('Translation and Rotated Data')
legend('Elements', 'Data')
