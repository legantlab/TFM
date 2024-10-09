%Script to automatically align two 3D channel point clouds. 

%Assume that two channels are rotated approximately the same, I.E channel
%runs up and down in Y for the image (prerotation is applied) so we can
%compute the YZ profile for top and bottom, then fit to a guassian and
%align using those functions. 

%Import swelling data
load('T:\Max\2023_07_08_SpinningDiskTesting\Swelling Data\averagedSwellingProfile7x40Final.mat');
%Scale swelling data set by pixel sizes
elemCents2D(:,1) = elemCents2D(:,1)/.199;
elemCents2D(:,2) = elemCents2D(:,2)/.199;
elemCents2D(:,3) = elemCents2D(:,3)/0.8;
elementsDim = [max(elemCents2D(:,1)), max(elemCents2D(:,2)),max(elemCents2D(:,3))]; %X,Y,Z
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
threshold = 0.995;
x{1} = locateParticles(I,threshmulti,bandpass_size,bPassParams,threshold);
x{1}(any(isnan(x{1}),2),:) = [];

figure
scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
hold on 
scatter3(elemCents2D(:,1), elemCents2D(:,2), elemCents2D(:,3))
%% Now make slices along y in upper and lower quartile of data (exluding constriction)
%First pass, turns out we need to correct for rotation first then we can do
%xyz translation correction using a variety of fit methods. 
filterUp = x{1}(:,2)>=ImgDim(:,1)*.75;
elementsUp = elemCents2D(:,2) >=elementsDim(:,2)*.75;
%Invert for guassian fit. 
dataupFit = [x{1}(filterUp,1),-1*x{1}(filterUp,3)];
elemupFit = [elemCents2D(elementsUp,1), -1*elemCents2D(elementsUp,3)];
%fix z to zero
dataupFit(:,2) = dataupFit(:,2) - min(dataupFit(:,2));
elemupFit(:,2) = elemupFit(:,2) - min(elemupFit(:,2));
%Fit to Gaussian profile
dataupFitval = fit(dataupFit(:,1),dataupFit(:,2), 'gauss1');
elemupFitval = fit(elemupFit(:,1),elemupFit(:,2), 'gauss1');

%plot fits
figure
scatter(dataupFit(:,1),dataupFit(:,2))
hold on 
scatter(elemupFit(:,1), elemupFit(:,2))
plot(dataupFitval)
plot(elemupFitval)

%Perform alignment using fits
datayupCent = dataupFitval.b1;
elemyupCent = elemupFitval.b1;
ydiffUp = elemyupCent - datayupCent;


filterDown = x{1}(:,2)<=ImgDim(:,1)*.25;
elementsDown = elemCents2D(:,2) <=elementsDim(:,2)*.25;

%Invert for guassian fit. 
datadownFit = [x{1}(filterDown,1),-1*x{1}(filterDown,3)];
elemdownFit = [elemCents2D(elementsDown,1), -1*elemCents2D(elementsDown,3)];
%fix z to zero
datadownFit(:,2) = datadownFit(:,2) - min(datadownFit(:,2));
elemdownFit(:,2) = elemdownFit(:,2) - min(elemdownFit(:,2));
%Fit to Gaussian profile
datadownFitval = fit(datadownFit(:,1),datadownFit(:,2), 'gauss1');
elemdownFitval = fit(elemdownFit(:,1),elemdownFit(:,2), 'gauss1');

%Perform alignment using fits
dataydownCent = datadownFitval.b1;
elemydownCent = elemdownFitval.b1;
ydiffDown = elemydownCent - dataydownCent;

%Plot fits
figure
scatter(datadownFit(:,1),datadownFit(:,2))
hold on 
scatter(elemdownFit(:,1), elemdownFit(:,2))
plot(datadownFitval)
plot(elemdownFitval)

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

%% Now correct in x using a similar strategy

%Use top 5 microns of substrate to fit a line and correct for slope in
%experimental and swelling data
expData = rmoutliers(x{1},'median');

dataTopfilt = expData(:,3) >= mean(expData(:,3)) -5;
elemTopfilt = newElements(:,3) >= max(newElements(:,3)) - 1;

%fit to line
dataTopfit = fit(expData(dataTopfilt,2), expData(dataTopfilt,3),'poly1');
elemTopfit = fit(newElements(elemTopfilt,2), newElements(elemTopfilt,3),'poly1');

figure
scatter(expData(dataTopfilt,2), expData(dataTopfilt,3))
hold on 
scatter(newElements(elemTopfilt,2), newElements(elemTopfilt,3))
plot(dataTopfit)
plot(elemTopfit)

%Rotate using slope, theta we want is just tan^-1 of slope, thanks
%pythagoras
xRotdata = atan(dataTopfit.p1);
xRotelem = atan(elemTopfit.p1);

expData = rotate_3D(x{1}','x',-1*xRotdata)';
newElements = rotate_3D(newElements','x',-1*xRotelem)';
figure
scatter(expData(:,2), expData(:,3))
hold on 
scatter(newElements(:,2), newElements(:,3))

%% Fix z rotation 
%This may be working now. 

distxElem = elementsDim(:,2)*.75 - elementsDim(:,2)*.25;
rotzElem = atan(distxElem/(elemupFitval.b1 - elemdownFitval.b1)); 
newElementsRotX = rotate_3D(newElements','z',-1*rotzElem)';

distxData = ImgDim(:,1)*.75 - ImgDim(:,1)*.25;
rotzData = atan(distxData/(dataupFitval.b1 - datadownFitval.b1));
expDataRotX = rotate_3D(expData','z',-1*rotzData)';

figure
scatter3(expData(:,1), expData(:,2),expData(:,3))
hold on 
scatter3(newElements(:,1), newElements(:,2),newElements(:,3))

figure
scatter3(newElementsRotX(:,1), newElementsRotX(:,2),newElementsRotX(:,3))
hold on
scatter3(expDataRotX(:,1), expDataRotX(:,2),expDataRotX(:,3))

%% Correct translational shifts now. 

%Correct for X via fitting gaussian profiles
elemfiltUp = newElementsRotX(:,2) >= max(newElementsRotX(:,2))*.75;
datafiltUp = expDataRotX(:,2) >= max(expDataRotX(:,2))*.75;
elemfiltDown = newElementsRotX(:,2) <= max(newElementsRotX(:,2))*.25;
datafiltDown = expDataRotX(:,2) <= max(expDataRotX(:,2))*.25;

elemfiltUpVals = newElementsRotX(elemfiltUp,:);
elemfiltDownVals = newElementsRotX(elemfiltDown,:);

datafiltUpVals = expDataRotX(datafiltUp,:);
datafiltDownVals = expDataRotX(datafiltDown,:);



