function [finalElements,finalDisplacements, trans] = alignChannels(swellingData,referenceLocs,plotFlag)
%alignChannels Corrects reference localizations using swelling dataset
%   Uses computed average swelling profile to correct reference
%   localizations for swelling correction in 3D TFM datasets after addition
%   of SDS. Assumes swelling is computed from an averaged dataset and is a
%   3D matrix of positions (centroids in model) and displacements at those
%   positions. Reference localizations are computed from locateparticle. 
%Input: 
% --------------------------------------------------------------------
% swellingData: nx6 matrix with centroid positions (n,1:3) and displacements
% (n,4:6) to be translated and rotated to referenceLocs. Assumes that
% channel feature runs in Y up the field of view. 
% 
% referenceLocs: nx3 matrix of localized beads from reference state of
% computed data. Assumes that channel feature runs in Y up the field of
% view. 
%
% plotFlag: a flag for outputting plots of each correction step for quality
% control. 0 for no plot, 1 for plots which are saved in a new folder in
% the current directory titled alignmentPlots. 
%
%Output:
% --------------------------------------------------------------------
% finalElements: nx3 matrix with x,y,z coordinates of corrected element
% centroids that are translationally and rotationally aligned with
% reference localizations. 
%
% finalDisplacements: nx3 matrix with x,y,z displacements with rotational
% correction applied. 
%
% trans: 1x6 vector containing linear transformations to applied to swelling
% data: x,y,z,xrot,yrot,zrot
% Notes: This can be generally used to align two 3D point clouds with this
% shape, just use a matrix of zeros for the displacement vectors in
% swelling data. 

%% Perform X,Z Alignment
if plotFlag ==1
    %create directory for saving figures
    mkdir alignmentPlots
end

x = {referenceLocs};

elemCents2D = swellingData(:,1:3);
averageDisplacement = swellingData(:,4:6);

% elemCents2D(:,1) = elemCents2D(:,1)/.199;
% elemCents2D(:,2) = elemCents2D(:,2)/.199;
% elemCents2D(:,3) = elemCents2D(:,3)/0.8;
elementsDim = [max(elemCents2D(:,1)) - (min(elemCents2D(:,1))), max(elemCents2D(:,2)) - (min(elemCents2D(:,2)))...
    ,max(elemCents2D(:,3)) - (min(elemCents2D(:,3)))]; %X,Y,Z, 

datafiltUp = x{1}(:,2)>=quantile(x{1}(:,2),.75);
elemfiltUp = elemCents2D(:,2) >=quantile(elemCents2D(:,2),.75);
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
if plotFlag == 1
    figure
    scatter(dataupFitVal(:,1),dataupFitVal(:,2))
    hold on 
    scatter(elemupFitVal(:,1), elemupFitVal(:,2))
    plot(dataupFit)
    plot(elemupFit)
    title('Upper Quartile fits')
    saveas(gcf,'alignmentPlots\UpperQuartileFits.png')
    close(gcf)
end
%Perform alignment using fits
datayupCent = dataupFit.b1;
elemyupCent = elemupFit.b1;
ydiffUp = elemyupCent - datayupCent;


datafiltDown = x{1}(:,2)<=quantile(x{1}(:,1),.25);
elemfiltDown = elemCents2D(:,2) <= quantile(elemCents2D(:,2),.25);

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
if plotFlag == 1
    figure
    scatter(datadownFitVal(:,1),datadownFitVal(:,2))
    hold on 
    scatter(elemdownFitVal(:,1), elemdownFitVal(:,2))
    plot(datadownFit)
    plot(elemdownFit)
    title('Lower Quartile Fits')
    saveas(gcf,'alignmentPlots\LowerQuartileFits.png')
    close(gcf)
end

% Perform X Alignment
%use average difference to roughly align channels in x. 
avgxshift = (ydiffDown + ydiffUp)/2;

%subtract from swelling dataset
newElements = elemCents2D;
newElements(:,1) = newElements(:,1) - avgxshift;

%plot to confirm
if plotFlag == 1
    figure
    scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
    hold on 
    scatter3(newElements(:,1), newElements(:,2), newElements(:,3))
    title('X Aligned Elements and Data')
    saveas(gcf,'alignmentPlots\Xalignment.png')
    close(gcf)
end
%% Perform Z Alignment (sensitive to outliers?)
%First compute general scaling from min/maxs
elemZRange = max(newElements(:,3)) - min(newElements(:,3));
sortZ = sort(newElements(:,3));
elemFilt = newElements(:,3) >= sortZ(end-25) *.9; %Assume there is less than 25 outliers

dataFilt = x{1}(:,3) >= max(x{1}(:,3)) - 10;
dataZRange = median(x{1}(dataFilt,3)) - min(x{1}(:,3));
elemZScale = dataZRange/elemZRange; 
newElements(:,3) = newElements(:,3)*elemZScale;
%Correct z displacements by new scaling as well
averageDisplacement(:,3) = averageDisplacement(:,3) * elemZScale;
elemzFit = fit(newElements(elemFilt, 2),newElements(elemFilt,3),'poly1');
zCorr = elemzFit.p2 - median(x{1}(dataFilt,3));
newElements(:,3) = newElements(:,3) - zCorr;
if plotFlag == 1
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
    saveas(gcf,'alignmentPlots\XZAlignment.png')
    close(gcf)
end

%% Perform Y Alignment
%filter from center +/- 20 pixel(~10 microns)
elementsDim = [max(newElements(:,1)) - abs(min(newElements(:,1))), max(newElements(:,2)) - abs(min(newElements(:,2)))...
    ,max(newElements(:,3)) - abs(min(newElements(:,3)))]; %X,Y,Z, 

datayFilt = x{1}(:,1) >= quantile(x{1}(:,2),.25) & x{1}(:,1) <= quantile(x{1}(:,2),.75);
elemyFilt = newElements(:,1) >= quantile(newElements(:,1),.25) & newElements(:,1) <= quantile(newElements(:,1),.75);
%fit y components from filter with guassian (histfit?)


datayFit = fitdist(x{1}(datayFilt,2),'Normal');
elemyFit = fitdist(newElements(elemyFilt,2),'Normal');

avgyshift = datayFit.mu - elemyFit.mu;

transElems = newElements;
transElems(:,2) = transElems(:,2) + avgyshift;

if plotFlag == 1
    figure
    scatter3(x{1}(:,1),x{1}(:,2), x{1}(:,3))
    hold on 
    scatter3(transElems(:,1), transElems(:,2), transElems(:,3))
    title('XYZ aligned Elements and Data')
    legend('Raw Data', 'Aligned Element Centroids')
    saveas(gcf,'alignmentPlots\XYZAlignment.png')
    close(gcf)
end

%% Now do rotational Alignments
%Before this, we may want to move the origin to 0,0,0 given that the
%rotations will otherwise be offset (in theory). 
%Rotate X by fitting top surface
sortZElem = sort(transElems(:,3));
sortZdata = sort(x{1}(:,3));
dataTopfilt = x{1}(:,3) >= sortZdata(end-25)*.9;
elemTopfilt = transElems(:,3) >= sortZElem(end-25)*.9;

%fit to line
dataTopfit = fit(x{1}(dataTopfilt,2), x{1}(dataTopfilt,3),'poly1');
elemTopfit = fit(transElems(elemTopfilt,2), transElems(elemTopfilt,3),'poly1');

if plotFlag == 1
    figure
    scatter(x{1}(dataTopfilt,2), x{1}(dataTopfilt,3))
    hold on 
    scatter(transElems(elemTopfilt,2), transElems(elemTopfilt,3))
    plot(dataTopfit)
    plot(elemTopfit)
    title('Top Surface Fits')
    saveas(gcf,'alignmentPlots\TopSurfaceFits.png')
    close(gcf)
end

xRotdata = 0; %atan(dataTopfit.p1);
xRotelem = atan(elemTopfit.p1) - atan(dataTopfit.p1);

rotData = rotate_3D(x{1}','x',-1*xRotdata)';
rotElems = rotate_3D(transElems','x',-1*xRotelem)';
rotDisps = rotate_3D(averageDisplacement','x',-1*xRotelem)';

if plotFlag == 1
    figure
    scatter3(rotData(:,1), rotData(:,2), rotData(:,3))
    hold on 
    scatter3(rotElems(:,1), rotElems(:,2), rotElems(:,3))
    title('XYZ Translation, X rotation Corrected')
    saveas(gcf,'alignmentPlots\XYZ_XrotAlignment.png')
    close(gcf)
end
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

if plotFlag == 1
    figure
    scatter(dataupVal(:,1),dataupVal(:,2))
    hold on 
    scatter(datadownVal(:,1),datadownVal(:,2))
    plot(dataupFit)
    plot(datadownFit)
    title('Data Y fits')
    saveas(gcf,'alignmentPlots\Yfits.png')
    close(gcf)
end

%Elements
elemupFilt = rotElems(:,2) >= datayrange*.75;
elemdownFilt = rotElems(:,2) <= datayrange*.25;
elemupVal = [rotElems(elemupFilt,1), -1*rotElems(elemupFilt,3) - min(-1*rotElems(elemupFilt,3))];
elemdownVal = [rotElems(elemdownFilt,1), -1*rotElems(elemdownFilt,3) - min(-1*rotElems(elemupFilt,3))];

elemupFit = fit(elemupVal(:,1), elemupVal(:,2),'gauss1');
elemdownFit = fit(elemdownVal(:,1), elemdownVal(:,2),'gauss1');

if plotFlag == 1
    figure
    scatter(elemupVal(:,1),elemupVal(:,2))
    hold on 
    scatter(elemdownVal(:,1),elemdownVal(:,2))
    plot(elemupFit)
    plot(elemdownFit)
    title('Elem Y fits')
    saveas(gcf,'alignmentPlots\ElemYFit.png')
    close(gcf)
    
    figure
    scatter(dataupVal(:,1),dataupVal(:,2))
    hold on
    scatter(elemupVal(:,1),elemupVal(:,2))
    plot(dataupFit)
    plot(elemupFit)
    title('Up YZ ranges')
    saveas(gcf,'alignmentPlots\UpYZRange.png')
    close(gcf)
    
    figure
    scatter(datadownVal(:,1),datadownVal(:,2))
    hold on
    scatter(elemdownVal(:,1),elemdownVal(:,2))
    plot(datadownFit)
    plot(elemdownFit)
    title('Down YZ ranges')
    saveas(gcf,'alignmentPlots\DownYZRange.png')
    close(gcf)
end

% Z Rotation Correction
avgzrot = ((elemupFit.b1 - dataupFit.b1));% + ((elemdownFit.b1 - datadownFit.b1)))/2; % I don't think this is the quantity we want actually. 
avgzrot2 = ((elemdownFit.b1 - datadownFit.b1));
distxElem = elemyrange; %(elemyrange*.75 - elemyrange*.25);
rotzElem = (atan(distxElem/avgzrot) - atan(distxElem/avgzrot2))/2 + pi/2; %Note, revisit this to double check trig/algebra - MH 8/24/23
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

if plotFlag == 1
    figure
    scatter(finalExpupVal(:,1), finalExpupVal(:,2))
    hold on 
    scatter(finalElemupVal(:,1), finalElemupVal(:,2))
    plot(finaldataupFit)
    plot(finalelemupFit)
    title('Upper Profiles')
    saveas(gcf,'alignmentPlots\UpperProfiles.png')
    close(gcf)
end
%Lower profiles
finalExpdownVal= [finalExpData(datadownFilt,1),-1*finalExpData(datadownFilt,3) - median(-1*finalExpData(datadownFilt,3))];
finalElemdownVal = [finalElements(elemdownFilt,1),-1*finalElements(elemdownFilt,3) - min(-1*finalElements(elemdownFilt,3))];
%finalExpdownVal= [finalExpData(datadownFilt,1),-1*finalExpData(datadownFilt,3)];

finaldatadownFit = fit(finalExpdownVal(:,1), finalExpdownVal(:,2),'gauss1');
finalelemdownFit = fit(finalElemdownVal(:,1), finalElemdownVal(:,2),'gauss1');
%finaldatadownFit = fit(finalExpdownVal(:,1), finalExpdownVal(:,2),'gauss1');

if plotFlag == 1
    figure
    scatter(finalExpdownVal(:,1), finalExpdownVal(:,2))
    hold on 
    scatter(finalElemdownVal(:,1), finalElemdownVal(:,2))
    plot(finaldatadownFit)
    plot(finalelemdownFit)
    title('Lower Profiles')
    saveas(gcf,'alignmentPlots\LowerProfiles.png')
    close(gcf)
    
    
    
    %% Final Plots
    
    figure
    scatter3(x{1}(:,1), x{1}(:,2), x{1}(:,3))
    hold on 
    scatter3(finalElements(:,1), finalElements(:,2), finalElements(:,3))
    quiver3(finalElements(:,1), finalElements(:,2), finalElements(:,3),finalDisplacements(:,1), finalDisplacements(:,2),...
        finalDisplacements(:,3),1)
    title('Translation and Rotated Data')
    legend('Elements', 'Data')
    saveas(gcf,'alignmentPlots\finalAlignment.png')
    close(gcf)
end


trans = [-1*avgxshift, avgyshift, zCorr, -1*xRotelem, 0, -1*rotzElem];


end