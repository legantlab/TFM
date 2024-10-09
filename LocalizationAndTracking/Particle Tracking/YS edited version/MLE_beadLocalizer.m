beads=loadtiff('..\Cell2_RFP_PreSDS.TIF');
ImgDim=size(beads);

[beadsBP]=bpass3dMB(beads, [.5,.5,.5], [2,2,2],[0,0]); %Bandpass filter - Can tune parameters for optimal fitting

threshold = 5000; %Minimum intensity in camera counts for a bead in the bandpass filtered image - can tune for optimal fitting
[PeakX,PeakY,PeakZ,PeakVal]=localMaximum(beadsBP,3,1,threshold);

[nPeaks,~]=size(PeakX);


%Half size of the sub-image for MLE fitting of centroid (e.g. if subx = 3, then subimage extent will be +- 3 pixels from the local maxima)
subx=2;
suby=2;
subz=2;

%Avoid beads at the edges of the volume
filter=(PeakX>subx)&(PeakX<(ImgDim(1)-subx))& (PeakY>suby)&(PeakY<(ImgDim(2)-suby))& (PeakZ>subz)&(PeakZ<(ImgDim(3)-subz));

PeakX=PeakX(filter);
PeakY=PeakY(filter);
PeakZ=PeakZ(filter);
PeakVal=PeakVal(filter);

[nPeaks,~]=size(PeakX);

I=zeros(nPeaks,1);
xc=I;
yc=I;
zc=I;
residual=I;

 %Perform MLE fitting of the bead centroids
 for ii=1:nPeaks
     subImage=double(beads(PeakX(ii)-subx:PeakX(ii)+subx,PeakY(ii)-suby:PeakY(ii)+suby,PeakZ(ii)-subz:PeakZ(ii)+subz));
      [pn,res,delta,u]=MLE3D_5Param_1Emitter(subImage,[max(subImage(:)),subx+1,suby+1,subz+1,1]);
                I(ii)=pn(1);
                xc(ii)=PeakX(ii)+(pn(3)-(subx+1));
                yc(ii)=PeakY(ii)+(pn(2)-(suby+1));
                zc(ii)=PeakZ(ii)+(pn(4)-(subz+1));
                residual(ii)=res;
ii
 end
 
 %Write file into format that Amira can read
fileID = fopen('..\Cell2_RFP_PreSDS_Coords.am','w');
fclose(fileID);
dlmwrite('..\Cell2_RFP_PreSDS_Coords.am',[yc-1,xc-1,zc-1],'-append','delimiter',' ');

x = [yc-1,xc-1,zc-1];
x_0test = coords_0;
x_0 = num2cell(coords_0);
save('x_0.mat', 'x_0')
save('x_0test.mat','x_0test');

