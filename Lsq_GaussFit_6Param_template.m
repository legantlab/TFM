function [pn,residuals]=Lsq_GaussFit_8Param_complete(PSF1,p_initial,display_fit)
% Lsq_Fit_6Param - Performs least square fitting to the image of a single
% molecule using a Gaussian point spread function. Will fit 6 components of
% the PSF as indicated below.
%
% Arguments:  PSF1              - an image of a single emitter
%             p_initial         - a 1x6 vector of initial guesses for the
%                                 different parameters of the Gaussian model

%                               p_initial(1); Intensity initial guess
%                               p_initial(2); x-centroid initial guess (in pixels)
%                               p_initial(3); y-centroid initial guess (in pixels)
%                               p_initial(4); sig_x initial guess (in pixels)
%                               p_initial(5); sig_y initial guess (in pixels)
%                               p_initial(6); background initial guess
%
%             display_fit       - a boolean that can be set to plot the fit image for debugging or display.
% 
% Returns:    pn                - a 1x6 vector of fit values for the
%                                 different parameters of the Gaussian model

%             residuals         - the residual value of the objective
%                                 function at the fit solution


[dimy,dimx]=size(PSF1); %Size of the image to be fit

 %Generate arrays of the x and y pixel coordinates using meshgrid
[X,Y]=meshgrid([1:1:dimx],[1:1:dimy]);

%Set options for Matlab's Lsq fitter. You can change these (especially the 'Display') for debugging.
options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 2000);

%Define the parameter limits for lsqnonlin
params_lb=[0,0,0,0,0,0]; %A 6x1 vector defining the lower bounds for the parameters to be fit
params_ub=[inf,dimx,dimy,inf,inf,inf]; %A 6x1 vector defining the upper bounds for the parameters to be fit (can be inf)

[pn, ~, residuals] = lsqnonlin(@Elliptical_Gaussian_Fit, p_initial,params_lb, params_ub, options);

%Plot the input image, the fit distribution, and the centroid estimate    
if display_fit==1
    Io= pn(1); %Intensity parameter
    xc= pn(2); %x-centroid parameter (in pixels)
    yc= pn(3); %y-centroid parameter (in pixels)
    sigx= pn(4); %Gaussian x-sigma parameter (in pixels)
    sigy= pn(5); %Gaussian y-sigma parameter (in pixels)
    bg= pn(6); %offset parameter
    Ex= 0.5*erf(((X) - xc + 0.5)/(sqrt(2)*sigx)) - 0.5*erf(((X) - xc - 0.5)/(sqrt(2)*sigx));
    Ey= 0.5*erf(((Y) - yc + 0.5)/(sqrt(2)*sigy)) - 0.5*erf(((Y) - yc - 0.5)/(sqrt(2)*sigy));
    PSF_IG= Io.*Ex.*Ey;
    subplot(1,2,1)
    imagesc(PSF1)
    hold on
    scatter(pn(2),pn(3),500,'mx');
    axis equal
    title('Input Image')
    xlabel('column index')
    ylabel('row index')
    subplot(1,2,2)
    imagesc(PSF_IG)
    hold on
    scatter(pn(2),pn(3),500,'mx');
    axis equal
    title('Fit Image')
    xlabel('column index')
    ylabel('row index')
end

function F=Elliptical_Gaussian_Fit(x)
    Io= x(1); %Intensity parameter
    xc= x(2); %x-centroid parameter (in pixels)
    yc= x(3); %y-centroid parameter (in pixels)
    sigx= x(4); %Gaussian x-sigma parameter (in pixels)
    sigy= x(5); %Gaussian y-sigma parameter (in pixels)
    bg= x(6); %offset parameter
    
    %Define the 2D Gaussian Distribution to use for fitting
    Ex= 0.5*erf(((X) - xc + 0.5)/(sqrt(2)*sigx)) - 0.5*erf(((X) - xc - 0.5)/(sqrt(2)*sigx));
    Ey= 0.5*erf(((Y) - yc + 0.5)/(sqrt(2)*sigy)) - 0.5*erf(((Y) - yc - 0.5)/(sqrt(2)*sigy));
    PSF_IG= Io.*Ex.*Ey + bg;
    
    %Compute the difference between the input image and the fit
    %distribution
    F= PSF_IG - PSF1; %A vector containing the pixel-by-pixel residuals between the fit and the input data
end
end