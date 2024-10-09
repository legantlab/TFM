%Script to do bead centroid fitting
close all
%Load Data
PSF1 = double(loadtiff('T:\Max\2023_07_01_IA32_SpinningDisk\Tiffs\F13\SingleBeads\Bead1.tif'));
PSF2 = double(loadtiff('T:\Max\2023_07_01_IA32_SpinningDisk\Tiffs\F13\SingleBeads\Bead2.tif'));
PSF3 = double(loadtiff('T:\Max\2023_07_01_IA32_SpinningDisk\Tiffs\F13\SingleBeads\Bead3.tif'));
ImgDim = size(PSF1);
p_initial = [max(PSF1(:))*.99,ImgDim(1)/2,ImgDim(2)/2,ImgDim(3)/2,750,0.2,0.2,0.2];
display_fit = 1;

[pn1,~]=Lsq_GaussFit_8Param(PSF1,p_initial,display_fit);
[pn2,~]=Lsq_GaussFit_8Param(PSF2,p_initial,display_fit);
[pn3,~]=Lsq_GaussFit_8Param(PSF3,p_initial,display_fit);

% Io=p_initial(1);
% xo=p_initial(2);
% yo=p_initial(3);
% zo=p_initial(4);
% bg=p_initial(5);
% sigx=p_initial(6);
% sigy=p_initial(7);
% sigz=p_initial(8);

% [pn1,~]=MLE3D_8Param_1Emitter(PSF1,p_initial);
% [pn2,~]=MLE3D_8Param_1Emitter(PSF2,p_initial);
% [pn3,~]=MLE3D_8Param_1Emitter(PSF3,p_initial);

fitPos = [pn1(2), pn1(3), pn1(4); pn2(2), pn2(3), pn2(4); pn3(2), pn3(3), pn3(4)];

dispX = (fitPos(:,1) - fitPos(1,1))*199;

dispY = (fitPos(:,2) - fitPos(1,2))*199;

dispZ = (fitPos(:,3) - fitPos(1,3))*100;

%Display average displacements in position
display(['AvgXDisp = ', num2str(mean(dispX(2:3)))])
display(['AvgYDisp = ', num2str(mean(dispY(2:3)))])
display(['AvgZDisp = ', num2str(mean(dispZ(2:3)))])
%% 

% %Display average sigmas
% disp(['Average SigX = ', num2str(mean([pn1(5),pn2(5), pn3(5)]))])
% 
% disp(['Average SigY = ', num2str(mean([pn1(6),pn2(6), pn3(6)]))])
% 
% disp(['Average SigZ = ', num2str(mean([pn1(7),pn2(7), pn3(7)]))])

%Display average sigmas
disp(['Average SigX = ', num2str(mean([pn1(6),pn2(6), pn3(6)]))])

disp(['Average SigY = ', num2str(mean([pn1(7),pn2(7), pn3(7)]))])

disp(['Average SigZ = ', num2str(mean([pn1(8),pn2(8), pn3(8)]))])