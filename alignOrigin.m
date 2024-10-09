function [finalLocalizations, trans] = alignOrigin(referenceLocs,plotFlag)
%Function that applies a xyz translation and rotation to a point cloud of
%channel localizations to bring the origin to the center of the channel xy
%and z is the top of the substrate. 
localizations = referenceLocs;
ImgDim = [max(localizations(:,1)) - min(localizations(:,1)),max(localizations(:,2)) - min(localizations(:,2)), ...
    max(localizations(:,3)) - min(localizations(:,3))];%this might need to be flipped
%Align X
datafiltUp = localizations(:,2)>=quantile(localizations(:,2),.75);

%Invert for guassian fit. 
dataupFitVal = [localizations(datafiltUp,1),-1*localizations(datafiltUp,3)];

%fix z to zero
dataupFitVal(:,2) = dataupFitVal(:,2) - median(dataupFitVal(:,2));

%Fit to Gaussian profile
dataupFit = fit(dataupFitVal(:,1),dataupFitVal(:,2), 'gauss1');

%plot fits
if plotFlag == 1
    figure
    scatter(dataupFitVal(:,1),dataupFitVal(:,2))
    hold on 
    plot(dataupFit)
    title('Upper Quartile fit')
    saveas(gcf,'alignmentPlots\UpperQuartileFit.fig')
    %close(gcf)
end
%Perform alignment using fits
datayupCent = dataupFit.b1;

ydiffUp = datayupCent;

datafiltDown = localizations(:,2) <= quantile(localizations(:,2),.25);


%Invert for guassian fit. 
datadownFitVal = [localizations(datafiltDown,1),-1*localizations(datafiltDown,3)];
%fix z to zero
datadownFitVal(:,2) = datadownFitVal(:,2) - median(datadownFitVal(:,2));

%Fit to Gaussian profile
datadownFit = fit(datadownFitVal(:,1),datadownFitVal(:,2), 'gauss1');

%Perform alignment using fits
dataydownCent = datadownFit.b1;
ydiffDown = dataydownCent;

%Plot fits
if plotFlag == 1
    figure
    scatter(datadownFitVal(:,1),datadownFitVal(:,2))
    hold on 
    plot(datadownFit)
    title('Lower Quartile Fits')
    saveas(gcf,'alignmentPlots\LowerQuartileFit.fig')
    %close(gcf)
end

% Perform X Alignment
%use average difference to roughly align channels in x. 
avgxshift = (ydiffDown + ydiffUp)/2;

%subtract from dataset
xCorLocalizations = localizations;
xCorLocalizations(:,1) = xCorLocalizations(:,1) - avgxshift;

%plot to confirm
if plotFlag == 1
    figure
    scatter3(xCorLocalizations(:,1), xCorLocalizations(:,2), xCorLocalizations(:,3))
    title('X Aligned Elements and Data')
    hold on 
    xline(0)
    saveas(gcf,'alignmentPlots\Xalignment.fig')
    %close(gcf)
end

%% Perform Z Alignment
%Compute the z value for the top of the channel by using filtering and a
%polynomial fitting? 
zSort = sort(xCorLocalizations(:,3));
zMax = zSort(end-25); %Assuming there is less than 25 outliers!
%First compute general scaling from min/maxs
dataFilt = xCorLocalizations(:,3) >= zMax - 3;
%Peform a fitting of the data using the filter
dataZFilt = fit(xCorLocalizations(dataFilt, 2),xCorLocalizations(dataFilt,3),'poly1');
%Correct z displacements
xzCorLocalizations = xCorLocalizations;
xzCorLocalizations(:,3) = xzCorLocalizations(:,3) - dataZFilt(quantile(xzCorLocalizations(:,2),.5));

if plotFlag == 1
    figure
    scatter(xCorLocalizations(:,2), xCorLocalizations(:,3))
    hold on 
    plot(dataZFilt)
    title('Z Fit')
    saveas(gcf,'alignmentPlots\zFit.fig')
    %close(gcf)

    figure
    scatter3(xzCorLocalizations(:,1),xzCorLocalizations(:,2), xzCorLocalizations(:,3))
    title('XZ Aligned Data')
    saveas(gcf,'alignmentPlots\XZAlignment.fig')
    %close(gcf)
end

%% Perform Y Alignment
%filter from center +/- 20 pixel(~10 microns)
datayFiltleft = xzCorLocalizations(:,1) >=  - 75 & xzCorLocalizations(:,1) <= 0 & xzCorLocalizations(:,3) > min(xzCorLocalizations(:,3)) + 5 & xzCorLocalizations(:,3) < - 5 ;%quantile(xzCorLocalizations(:,1),.5) = 0;
%fit y components from filter with guassian (histfit?) also filter out z
%localizations right at the bottom and top of the channels to remove some
%of the noise

datayFiltright = xzCorLocalizations(:,1) >= 0 & xzCorLocalizations(:,1) <=  75 & xzCorLocalizations(:,3) > min(xzCorLocalizations(:,3)) + 5 & xzCorLocalizations(:,3) < - 5;

datayFitleft = fit(xzCorLocalizations(datayFiltleft,2),xzCorLocalizations(datayFiltleft,1) + abs(min(xzCorLocalizations(datayFiltleft,1))),'gauss1');
datayFitright = fit(xzCorLocalizations(datayFiltright,2),-1*xzCorLocalizations(datayFiltright,1) + abs(min(-1*xzCorLocalizations(datayFiltright,1))),'gauss1');

avgyshift = (datayFitleft.b1 + datayFitright.b1)/2;

xzyCorLocalizations = xzCorLocalizations;
xzyCorLocalizations(:,2) = xzyCorLocalizations(:,2) - avgyshift;

if plotFlag == 1
    figure
    scatter(xzCorLocalizations(datayFiltleft,2), xzCorLocalizations(datayFiltleft,1) + abs(min(xzCorLocalizations(datayFiltleft,1))))
    title('Y Fit Data')
    hold on 
    scatter(xzCorLocalizations(datayFiltright,2), -1*xzCorLocalizations(datayFiltright,1) + abs(min(-1*xzCorLocalizations(datayFiltright,1))))
    plot(datayFitleft)
    plot(datayFitright)
    xlabel('y')
    ylabel('x')
    
    figure
    scatter3(xzyCorLocalizations(:,1), xzyCorLocalizations(:,2), xzyCorLocalizations(:,3))
    title('XYZ aligned Elements and Data')
    hold on 
    yline(0)
    xline(0)
    saveas(gcf,'alignmentPlots\XYZAlignment.fig')

    %close(gcf)
end

%% Perform rotational corrections to cartesian system now
sortZdata = sort(xzyCorLocalizations(:,3));
zMax = sortZdata(end-25); %Assuming there is less than 25 outliers!
dataTopfilt = xzyCorLocalizations(:,3) >= zMax - 3;

%fit to line
dataTopfit = fit(xzyCorLocalizations(dataTopfilt,2), xzyCorLocalizations(dataTopfilt,3),'poly1');

if plotFlag == 1
    figure
    scatter(xzyCorLocalizations(dataTopfilt,2), xzyCorLocalizations(dataTopfilt,3))
    hold on 
    plot(dataTopfit)
    title('Top Surface Fits')
    saveas(gcf,'alignmentPlots\TopSurfaceFits.fig')
    %close(gcf)
end

xRotdata = atan(dataTopfit.p1);

xzyXCorLocalizations = rotate_3D(xzyCorLocalizations','x',-1*xRotdata)';

if plotFlag == 1
    figure
    scatter3(xzyXCorLocalizations(:,1), xzyXCorLocalizations(:,2), xzyXCorLocalizations(:,3))
    title('XYZ Translation, X rotation Corrected')
    saveas(gcf,'alignmentPlots\XYZ_XrotAlignment.fig')
    %close(gcf)
end

%Plot yz profiles
%Data
dataupFilt = xzyXCorLocalizations(:,2) >= quantile(xzyXCorLocalizations(:,2),.75);
datadownFilt = xzyXCorLocalizations(:,2) <= quantile(xzyXCorLocalizations(:,2),.25);
dataupVal = [xzyXCorLocalizations(dataupFilt,1), -1*xzyXCorLocalizations(dataupFilt,3)...
    - median(-1*xzyXCorLocalizations(dataupFilt,3))];
datadownVal = [xzyXCorLocalizations(datadownFilt,1), -1*xzyXCorLocalizations(datadownFilt,3)...
    - median(-1*xzyXCorLocalizations(dataupFilt,3))];

dataupFit = fit(dataupVal(:,1), dataupVal(:,2),'gauss1');
datadownFit = fit(datadownVal(:,1), datadownVal(:,2),'gauss1');

if plotFlag == 1
    figure
    scatter(dataupVal(:,1),dataupVal(:,2))
    hold on 
    scatter(datadownVal(:,1),datadownVal(:,2))
    plot(dataupFit)
    plot(datadownFit)
    title('Data Y fits')
    saveas(gcf,'alignmentPlots\Yfits.fig')
    %close(gcf)
end

% Z Rotation Correction
rotDist = dataupFit.b1 - datadownFit.b1;
distx = quantile(xzyXCorLocalizations(:,2),.75) - quantile(xzyXCorLocalizations(:,2),.25); 
rotz = atan(distx/rotDist) + pi/2; %Note, revisit this to double check trig/algebra - MH 8/24/23
finalLocalizations = rotate_3D(xzyXCorLocalizations','z',-1*rotz)';

%Plot
if plotFlag == 1
    figure
    scatter3(finalLocalizations(:,1), finalLocalizations(:,2), finalLocalizations(:,3))
    title('Aligned Data Set')
    saveas(gcf,'alignmentPlots\finalAlignment.fig')
    %close(gcf)
end

trans = [avgxshift,avgyshift,dataZFilt(quantile(finalLocalizations(:,2),.5)),xRotdata,0,rotz];
end