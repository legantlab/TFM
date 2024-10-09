%Load bead image
close all
I = loadtiff('U:\Max\2023_05_24_DragonflySpinningDiscTesting\2023-05-24\Processed\300nmZStep\300nmZStepSingleBead1.tif');
ImgDim=size(I);
bandpass_size = bPassParams{2};

[IBP]=bpass3dMB(I, bPassParams{1},bPassParams{2},bPassParams{3}); %Bandpass filter - Can tune parameters for optimal fitting
% %These parameters are hard coded; it might be advantegous to make another
% %GUI here for ease of use and selection of optimal parameters. At least
% %display images from the start overlated with bandpass with settings shown
% %and be able to redo/edit them from there. 
% 
% figure
% imagesc(max(IBP,[],3))
% figure
% imagesc(squeeze(max(IBP,[],2)))
%%
minDist = mean(bandpass_size);           % minimum distance between two peaks for local maximum detection
% Bpass filter sets large parts of the array to zero which means the
% prctile returns 0 most of the time which obviously is a problem

%threshold = prctile(IBP(IBP~=0),threshold*100); %*threshold_multi messes this up for some reason. 

[PeakX,PeakY,PeakZ,PeakVal]=localMaximum(IBP,minDist,1,threshold);
%Insert test cases here


[nPeaks,~]=size(PeakX);

%Half size of the sub-image for MLE fitting of centroid (e.g. if subx = 3, then subimage extent will be +- 3 pixels from the local maxima)
subx=floor((bandpass_size(1)-1)/2);
suby=floor((bandpass_size(2)-1)/2);
subz=floor((bandpass_size(3)-1)/2);

%Avoid beads at the edges of the volume
filter=(PeakX>subx)&(PeakX<(ImgDim(1)-subx))& (PeakY>suby)&(PeakY<(ImgDim(2)-suby))& (PeakZ>subz)&(PeakZ<(ImgDim(3)-subz));

PeakX=PeakX(filter);
PeakY=PeakY(filter);
PeakZ=PeakZ(filter);
PeakVal=PeakVal(filter);

[nPeaks,~]=size(PeakX);

I_2=zeros(nPeaks,1);
xc=I_2;
yc=I_2;
zc=I_2;
residual=I_2;


%% Optional 8 parameter fitting for sigmoid determination:
%p_initial = [max(I(:)),subx+ 1,suby+1,subz+1,1,.5,.5,.5]; 
sigX0 = 2;
sigY0 = 2;
sigZ0 = 5;

[pn2,r,delta,u]=MLE3D_8Param_1Emitter(I,[max(I(:)),ImgDim(1)/2,ImgDim(2)/2,ImgDim(3)/2,1,sigX0,sigY0,sigZ0]);

disp(pn2(5:8))

figure
imagesc(max(I,[],3))
hold on 
scatter(pn2(2), pn2(3), 'r')
figure
imagesc(squeeze(max(I,[],2)))
hold on 
scatter(pn2(4), pn2(2), 'r')
 %% Perform MLE fitting of the bead centroids. Implemented parallel computing
 %here for a much needed speed boost, will consider GPU parallalization but
 %CPU is fine for now. 
 % --MH 2020
 for ii=1:nPeaks
     subImage=double(I(PeakX(ii)-subx:PeakX(ii)+subx,PeakY(ii)-suby:PeakY(ii)+suby,PeakZ(ii)-subz:PeakZ(ii)+subz));
      [pn,res,delta,u]=MLE3D_5Param_1Emitter(subImage,[max(subImage(:)),subx+1,suby+1,subz+1,1]);
                I_2(ii)=pn(1);
                xc(ii)=PeakX(ii)+(pn(3)-(subx+1));
                yc(ii)=PeakY(ii)+(pn(2)-(suby+1));
                zc(ii)=PeakZ(ii)+(pn(4)-(subz+1));
                residual(ii)=res;
 end

x = [xc,yc,zc];
x(~any(~isnan(x),2),:) = [];
close all
figure
subplot(1,2,1)
imagesc((max(I,[],3)))
hold on 
scatter(x(1), x(2), 'r')
%Plot guassian using fit parameters
subplot(1,2,2)
fitGauss = gauss2d(zeros(ImgDim(1),ImgDim(2)),pn2(6),[xc,yc]);
imagesc(fitGauss)

figure
subplot(1,2,1)
imagesc(squeeze(max(I,[],2)))
hold on 
scatter(x(3), x(1), 'r')
subplot(1,2,2)
fitGauss = gauss2d(zeros(ImgDim(3),ImgDim(1)),pn2(8),[zc,xc]);
imagesc(fitGauss)
