function [x] = locateParticles(I,bandpass_size,bPassParams,thresh)
% [x] = locateParticles(I, beadParameter) locates particles in the image
%
% INPUTS
% -------------------------------------------------------------------------
%   I:              Input volumetric image
%   beadParameter:  Parameters to detect particle position in images
%   bandpasssize:   Size of the long pass filter
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              Voxel-level estimate of particle center in MxNxO format
%
%   Author: Max Hockenberry
%   Last Update: 10/23/2024

% Parameters
ImgDim=size(I);

[IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3}); %Bandpass filter - Can tune parameters for optimal fitting
% %These parameters are hard coded; it might be advantegous to make another
% %GUI here for ease of use and selection of optimal parameters. At least
% %display images from the start overlated with bandpass with settings shown
% %and be able to redo/edit them from there. 

minDist = mean(bandpass_size);           % minimum distance between two peaks for local maximum detection
% Bpass filter sets large parts of the array to zero which means the
% prctile returns 0 most of the time which obviously is a problem

threshold = prctile(IBP(IBP~=0),thresh*100); %*threshold_multi messes this up for some reason. 

[PeakX,PeakY,PeakZ,~]=localMaximum(IBP,minDist,1,threshold);

%Half size of the sub-image for MLE fitting of centroid (e.g. if subx = 3, then subimage extent will be +- 3 pixels from the local maxima)
subx=floor((bandpass_size(1)-1)/2);
suby=floor((bandpass_size(2)-1)/2);
subz=floor((bandpass_size(3)-1)/2);

%Avoid beads at the edges of the volume
filter=(PeakX>subx)&(PeakX<(ImgDim(1)-subx))& (PeakY>suby)&(PeakY<(ImgDim(2)-suby))& (PeakZ>subz)&(PeakZ<(ImgDim(3)-subz));

PeakX=PeakX(filter);
PeakY=PeakY(filter);
PeakZ=PeakZ(filter);

[nPeaks,~]=size(PeakX);

I_2=zeros(nPeaks,1);
xc=I_2;
yc=I_2;
zc=I_2;
residual=I_2;

 %Perform MLE fitting of the bead centroids. Implemented parallel computing
 %here for a much needed speed boost, will consider GPU parallalization but
 %CPU is fine for now. 
 % --MH 2020
 parfor ii=1:nPeaks
     subImage=double(I(PeakX(ii)-subx:PeakX(ii)+subx,PeakY(ii)-suby:PeakY(ii)+suby,PeakZ(ii)-subz:PeakZ(ii)+subz));
      [pn,res,delta,u]=MLE3D_5Param_1Emitter(subImage,[max(subImage(:)),subx+1,suby+1,subz+1,1]);
                I_2(ii)=pn(1);
                xc(ii)=PeakX(ii)+(pn(3)-(subx+1));
                yc(ii)=PeakY(ii)+(pn(2)-(suby+1));
                zc(ii)=PeakZ(ii)+(pn(4)-(subz+1));
                residual(ii)=res;
 end
x = [yc-1,xc-1,zc-1];
x(~any(~isnan(x),2),:) = []; %remove any nans
end

