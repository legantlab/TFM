function [pn,CRLB,u,rChi,delta,convergeflag,pval,chi2MLE]=MLE3D_LevMarq_8Param(PSF1,X,Y,Z,p_initial,maxIter)
%A function to perform maximum likelihood fitting of single emitter spots.  This function fits the integrated intensity, the x-centroid, y-centroid, x-sigma,
%y-sigma, and per pixel background under a symmetric 2D Gaussian PSF model.

% Input Parameters
% PSF1 = a N X N X N 3D array of the photon-converted image data to fit
% p_initial = a 1 X 8 vector containing the initial guesses for the fit parameters
% maxIter = the maximum number of iterations to take before exiting

%Output Parameters
% pn = a 1 X 6 vector of the fitted parameters
%     pn(1) = the integrated intensity of the spot
%     pn(2) = the x-centroid of the spot
%     pn(3) = the y-centroid of the spot
%     pn(4) = the z-centroid of the spot
%     pn(5) = the x-sigma of the spot
%     pn(6) = the y-sigma of the spot
%     pn(7) = the y-sigma of the spot
%     pn(8) = the per pixel background of the spot
% r = the sum of the squared per-pixel differences between the image data and the model fit
% delta = an iters X 5 matrix displaying the amount that each variable is incrementated by in each iteration of the Newton solver
% flag = a flag to indicate algorithm convergence. If change in the value of all parameters is less than 1%, then flag = 1, else flag = 0

Io=p_initial(1); % initial guess for the integrated intensity of the spot
xo=p_initial(2); % initial guess for the x-centroid of the spot
yo=p_initial(3); % initial guess for the y-centroid of the spot
zo=p_initial(4); % initial guess for the y-centroid of the spot
sigx=p_initial(5); % initial guess for the x-sigma of the spot
sigy=p_initial(6); % initial guess for the y-sigma of the spot
sigz=p_initial(7); % initial guess for the y-sigma of the spot
bg=p_initial(8); % initial guess for the per pixel background of the spot

%Initialize variables
rChi=0;
delta=zeros(8,maxIter);
diagnostics=1; %Turn this on to plot the fit data
funTol = 1e-5; %The minimum threshold to define how much the Chi-square value needs to decrease in order to take a step
convergeflag=0;
FSHR=zeros(8,8);
CRLB=zeros(8,1);

%Define image model
Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sqrt(sigx^2))) + erf((0.5 + X - xo)./(sqrt(2).*sqrt(sigx^2))));
Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sqrt(sigy^2))) + erf((0.5 + Y - yo)./(sqrt(2).*sqrt(sigy^2))));
Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sqrt(sigz^2))) + erf((0.5 + Z - zo)./(sqrt(2).*sqrt(sigz^2))));
u=sqrt(bg^2)+sqrt(Io^2)*Ex.*Ey.*Ez;
chi2MLE=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));

for i=1:maxIter
    i
    %Calculate 1st derivatives
    du1_dIo=Io/sqrt(Io^2).*Ex.*Ey.*Ez;
    du1_dxo=(1./2).*sqrt(Io^2).*(sqrt(2./pi)./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx) - sqrt(2./pi)./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx)).*Ey.*Ez;
    du1_dyo=(1./2).*sqrt(Io^2).*(sqrt(2./pi)./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy) - sqrt(2./pi)./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy)).*Ex.*Ez;
    du1_dzo=(1./2).*sqrt(Io^2).*(sqrt(2./pi)./(exp((-0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz) - sqrt(2./pi)./(exp((0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz)).*Ex.*Ey;
    du1_dsigx=(1./2).*sqrt(Io^2).*((sqrt(2./pi).*(-0.5 + X - xo))./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^2) - (sqrt(2./pi).*(0.5 + X - xo))./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^2)).*Ey.*Ez;
    du1_dsigy=(1./2).*sqrt(Io^2).*((sqrt(2./pi).*(-0.5 + Y - yo))./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^2) - (sqrt(2./pi).*(0.5 + Y - yo))./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^2)).*Ex.*Ez;
    du1_dsigz=(1./2).*sqrt(Io^2).*((sqrt(2./pi).*(-0.5 + Z - zo))./(exp((-0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz.^2) - (sqrt(2./pi).*(0.5 + Z - zo))./(exp((0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz.^2)).*Ex.*Ey;
    du1_dbg=bg/sqrt(bg^2);
    
    %Calculate Jacobian
    B_k(1,1)=-sum((1-PSF1(:)./u(:)).*du1_dIo(:));
    B_k(2,1)=-sum((1-PSF1(:)./u(:)).*du1_dxo(:));
    B_k(3,1)=-sum((1-PSF1(:)./u(:)).*du1_dyo(:));
    B_k(4,1)=-sum((1-PSF1(:)./u(:)).*du1_dzo(:));
    B_k(5,1)=-sum((1-PSF1(:)./u(:)).*du1_dsigx(:));
    B_k(6,1)=-sum((1-PSF1(:)./u(:)).*du1_dsigy(:));
    B_k(7,1)=-sum((1-PSF1(:)./u(:)).*du1_dsigz(:));
    B_k(8,1)=-sum((1-PSF1(:)./u(:)).*du1_dbg(:));
    
    %Calculate Hessian - setting second derivatives to zeros
    lbda=0;
    A_kl(1,1)=sum(du1_dIo(:).*du1_dIo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(2,1)=sum(du1_dIo(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,1)=sum(du1_dIo(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,1)=sum(du1_dIo(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,1)=sum(du1_dIo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,1)=sum(du1_dIo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,1)=sum(du1_dIo(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,1)=sum(du1_dIo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    
    A_kl(1,2)=sum(du1_dxo(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,2)=sum(du1_dxo(:).*du1_dxo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(3,2)=sum(du1_dxo(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,2)=sum(du1_dxo(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,2)=sum(du1_dxo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,2)=sum(du1_dxo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,2)=sum(du1_dxo(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,2)=sum(du1_dxo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,3)=sum(du1_dyo(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,3)=sum(du1_dyo(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,3)=sum(du1_dyo(:).*du1_dyo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(4,3)=sum(du1_dyo(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,3)=sum(du1_dyo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,3)=sum(du1_dyo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,3)=sum(du1_dyo(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,3)=sum(du1_dyo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,4)=sum(du1_dzo(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,4)=sum(du1_dzo(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,4)=sum(du1_dzo(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,4)=sum(du1_dzo(:).*du1_dzo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(5,4)=sum(du1_dzo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,4)=sum(du1_dzo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,4)=sum(du1_dzo(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,4)=sum(du1_dzo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,5)=sum(du1_dsigx(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,5)=sum(du1_dsigx(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,5)=sum(du1_dsigx(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,5)=sum(du1_dsigx(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,5)=sum(du1_dsigx(:).*du1_dsigx(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(6,5)=sum(du1_dsigx(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,5)=sum(du1_dsigx(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,5)=sum(du1_dsigx(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,6)=sum(du1_dsigy(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,6)=sum(du1_dsigy(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,6)=sum(du1_dsigy(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,6)=sum(du1_dsigy(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,6)=sum(du1_dsigy(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,6)=sum(du1_dsigy(:).*du1_dsigy(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(7,6)=sum(du1_dsigy(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,6)=sum(du1_dsigy(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,7)=sum(du1_dsigz(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,7)=sum(du1_dsigz(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,7)=sum(du1_dsigz(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,7)=sum(du1_dsigz(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,7)=sum(du1_dsigz(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,7)=sum(du1_dsigz(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,7)=sum(du1_dsigz(:).*du1_dsigz(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(8,7)=sum(du1_dsigz(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,8)=sum(du1_dbg(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,8)=sum(du1_dbg(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,8)=sum(du1_dbg(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,8)=sum(du1_dbg(:).*du1_dzo(:).*PSF1(:)./u(:).^2);
    A_kl(5,8)=sum(du1_dbg(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(6,8)=sum(du1_dbg(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(7,8)=sum(du1_dbg(:).*du1_dsigz(:).*PSF1(:)./u(:).^2);
    A_kl(8,8)=sum(du1_dbg(:).*du1_dbg(:).*PSF1(:)./u(:).^2)*(1+lbda);
    
    delta(:,i)=pinv(A_kl)*B_k;
    
    Io_test=Io+delta(1,i);
    xo_test=xo+delta(2,i);
    yo_test=yo+delta(3,i);
    zo_test=zo+delta(4,i);
    sigx_test=sigx+delta(5,i);
    sigy_test=sigy+delta(6,i);
    sigz_test=sigz+delta(7,i);
    bg_test=bg+delta(8,i);
    
    %     Test whether update gives a more favorable solution as measured by a greater likelihood
    passed=0;
    t=1;
    while (~passed && t<=50 &&~convergeflag)
        %Recalculate the objective function at the proposed parameter location
        Ex_test=(1./2).*(-erf((-0.5 + X - xo_test)./(sqrt(2).*sigx_test)) + erf((0.5 + X - xo_test)./(sqrt(2).*sigx_test)));
        Ey_test=(1./2).*(-erf((-0.5 + Y - yo_test)./(sqrt(2).*sigy_test)) + erf((0.5 + Y - yo_test)./(sqrt(2).*sigy_test)));
        Ez_test=(1./2).*(-erf((-0.5 + Z - zo_test)./(sqrt(2).*sigz_test)) + erf((0.5 + Z - zo_test)./(sqrt(2).*sigz_test)));
        u_test=sqrt(bg_test^2)+sqrt(Io_test^2).*Ex_test.*Ey_test.*Ez_test;
        chi2MLE_test=2*sum(u_test(:)-PSF1(:))-2.*sum(PSF1(:).*log(u_test(:)./PSF1(:)));
        passed=chi2MLE_test<=(chi2MLE+funTol);%Compare with the previous value
        %If the value doesn't decrease (or doesn't stay within funTol of the previous value - suggesting convergence), try a smaller step
        if (~passed && t<50)
            lbda=2^t;
            t=t+1;
            A_kl(1,1)=sum(du1_dIo(:).*du1_dIo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(2,2)=sum(du1_dxo(:).*du1_dxo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(3,3)=sum(du1_dyo(:).*du1_dyo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(4,4)=sum(du1_dzo(:).*du1_dzo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(5,5)=sum(du1_dsigx(:).*du1_dsigx(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(6,6)=sum(du1_dsigy(:).*du1_dsigy(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(7,7)=sum(du1_dsigz(:).*du1_dsigz(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(8,8)=sum(du1_dbg(:).*du1_dbg(:).*PSF1(:)./u(:).^2)*(1+lbda);
            delta(:,i)=pinv(A_kl)*B_k;
            Io_test=Io+delta(1,i);
            xo_test=xo+delta(2,i);
            yo_test=yo+delta(3,i);
            zo_test=zo+delta(4,i);
            sigx_test=sigx+delta(5,i);
            sigy_test=sigy+delta(6,i);
            sigz_test=sigz+delta(7,i);
            bg_test=bg+delta(8,i);
        else if t>=50
                warning('The smallest step size has been reached, but algorithm has not converged')
                pn=[0,0,0,0,0,0,0];
                pval = -1;
                return
            else %Update parameter values
                Io=Io_test;
                xo=xo_test;
                yo=yo_test;
                zo=zo_test;
                sigx=sigx_test;
                sigy=sigy_test;
                sigz=sigz_test;
                bg=bg_test;
                Ex=Ex_test;
                Ey=Ey_test;
                Ez=Ez_test;
                u=u_test;
                convergeflag=((chi2MLE-chi2MLE_test)/chi2MLE)<funTol;
            end
        end
    end
    chi2MLE=chi2MLE_test;
end
pn=[sqrt(Io^2),xo,yo,zo,sigx,sigy,sigz,sqrt(bg^2)];
%Calculate Fisher information matrix
FSHR(1,1)=sum(du1_dIo(:).*du1_dIo(:)./u(:));
FSHR(1,2)=sum(du1_dIo(:).*du1_dxo(:)./u(:));
FSHR(1,3)=sum(du1_dIo(:).*du1_dyo(:)./u(:));
FSHR(1,4)=sum(du1_dIo(:).*du1_dzo(:)./u(:));
FSHR(1,5)=sum(du1_dIo(:).*du1_dsigx(:)./u(:));
FSHR(1,6)=sum(du1_dIo(:).*du1_dsigy(:)./u(:));
FSHR(1,7)=sum(du1_dIo(:).*du1_dsigz(:)./u(:));
FSHR(1,8)=sum(du1_dIo(:).*du1_dbg(:)./u(:));


FSHR(2,1)=FSHR(1,2);
FSHR(2,2)=sum(du1_dxo(:).*du1_dxo(:)./u(:));
FSHR(2,3)=sum(du1_dxo(:).*du1_dyo(:)./u(:));
FSHR(2,4)=sum(du1_dxo(:).*du1_dzo(:)./u(:));
FSHR(2,5)=sum(du1_dxo(:).*du1_dsigx(:)./u(:));
FSHR(2,6)=sum(du1_dxo(:).*du1_dsigy(:)./u(:));
FSHR(2,7)=sum(du1_dxo(:).*du1_dsigz(:)./u(:));
FSHR(2,8)=sum(du1_dxo(:).*du1_dbg(:)./u(:));


FSHR(3,1)=FSHR(1,3);
FSHR(3,2)=FSHR(2,3);
FSHR(3,3)=sum(du1_dyo(:).*du1_dyo(:)./u(:));
FSHR(3,4)=sum(du1_dyo(:).*du1_dzo(:)./u(:));
FSHR(3,5)=sum(du1_dyo(:).*du1_dsigx(:)./u(:));
FSHR(3,6)=sum(du1_dyo(:).*du1_dsigy(:)./u(:));
FSHR(3,7)=sum(du1_dyo(:).*du1_dsigz(:)./u(:));
FSHR(3,8)=sum(du1_dyo(:).*du1_dbg(:)./u(:));

FSHR(4,1)=FSHR(1,4);
FSHR(4,2)=FSHR(2,4);
FSHR(4,3)=FSHR(3,4);
FSHR(4,4)=sum(du1_dzo(:).*du1_dzo(:)./u(:));
FSHR(4,5)=sum(du1_dzo(:).*du1_dsigx(:)./u(:));
FSHR(4,6)=sum(du1_dzo(:).*du1_dsigy(:)./u(:));
FSHR(4,7)=sum(du1_dzo(:).*du1_dsigz(:)./u(:));
FSHR(4,8)=sum(du1_dzo(:).*du1_dbg(:)./u(:));

FSHR(5,1)=FSHR(1,5);
FSHR(5,2)=FSHR(2,5);
FSHR(5,3)=FSHR(3,5);
FSHR(5,4)=FSHR(4,5);
FSHR(5,5)=sum(du1_dsigx(:).*du1_dsigx(:)./u(:));
FSHR(5,6)=sum(du1_dsigx(:).*du1_dsigy(:)./u(:));
FSHR(5,7)=sum(du1_dsigx(:).*du1_dsigz(:)./u(:));
FSHR(5,8)=sum(du1_dsigx(:).*du1_dbg(:)./u(:));

FSHR(6,1)=FSHR(1,6);
FSHR(6,2)=FSHR(2,6);
FSHR(6,3)=FSHR(3,6);
FSHR(6,4)=FSHR(4,6);
FSHR(6,5)=FSHR(5,6);
FSHR(6,6)=sum(du1_dsigy(:).*du1_dsigy(:)./u(:));
FSHR(6,7)=sum(du1_dsigy(:).*du1_dsigz(:)./u(:));
FSHR(6,8)=sum(du1_dsigy(:).*du1_dbg(:)./u(:));

FSHR(7,1)=FSHR(1,7);
FSHR(7,2)=FSHR(2,7);
FSHR(7,3)=FSHR(3,7);
FSHR(7,4)=FSHR(4,7);
FSHR(7,5)=FSHR(5,7);
FSHR(7,6)=FSHR(6,7);
FSHR(7,7)=sum(du1_dsigz(:).*du1_dsigz(:)./u(:));
FSHR(7,8)=sum(du1_dsigz(:).*du1_dbg(:)./u(:));

FSHR(8,1)=FSHR(1,8);
FSHR(8,2)=FSHR(2,8);
FSHR(8,3)=FSHR(3,8);
FSHR(8,4)=FSHR(4,8);
FSHR(8,5)=FSHR(5,8);
FSHR(8,6)=FSHR(6,8);
FSHR(8,7)=FSHR(7,8);
FSHR(8,8)=sum(du1_dsigz(:).*du1_dbg(:)./u(:));

%Calculate CRLB
FSHRinv=inv(FSHR);
CRLB=diag(FSHRinv);

rChi=sum((PSF1(:)-u(:)).^2./PSF1(:));
pval = 0;
%

%diagnostics provides graphical plots of the function behavior in the neighborhood of the final parameter values
if diagnostics==1
%     figure
%     %Check molecule fit compared to raw data
%     Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
%     Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
%     Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
%     u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
%     subplot(1,2,1)
%     hold on
%     imagesc(u)
%     scatter(pn(3),pn(2),'bo','filled')
%     axis equal
%     title('fit image')
%     subplot(1,2,2)
%     hold on
%     imagesc(PSF1)
%     scatter(pn(3),pn(2),'bo','filled')
%     axis equal
%     title('original image')

    Irange=[sqrt(Io^2)-1000:1:sqrt(Io^2)+1000]';
    figure
    subplot(2,4,1)
    hold on
    for ii=1:length(Irange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Irange(ii)^2).*Ex.*Ey.*Ez;
        chi2MLE_I(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(Irange,chi2MLE_I,'b')
    plot(sqrt(Io^2),chi2MLE,'bo')
    xlabel('Io')
    ylabel('Log-Likelihood')
    
    subplot(2,4,2)
    hold on
    xrange=[xo-2:.01:xo+2]';
    for ii=1:length(xrange)
        Ex=(1./2).*(-erf((-0.5 + X - xrange(ii))./(sqrt(2).*sigx)) + erf((0.5 + X - xrange(ii))./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_x(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(xrange,chi2MLE_x,'r')
    plot(xo,chi2MLE,'ro')
    xlabel('xo')
    ylabel('Log-Likelihood')
    
    subplot(2,4,3)
    hold on
    yrange=[yo-2:.01:yo+2]';
    for ii=1:length(yrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yrange(ii))./(sqrt(2).*sigy)) + erf((0.5 + Y - yrange(ii))./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_y(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(yrange,chi2MLE_y,'g')
    plot(yo,chi2MLE,'go')
    xlabel('yo')
    ylabel('Log-Likelihood')
    
    subplot(2,4,4)
    hold on
    zrange=[zo-2:.01:zo+2]';
    for ii=1:length(yrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zrange(ii))./(sqrt(2).*sigz)) + erf((0.5 + Z - zrange(ii))./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_z(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(zrange,chi2MLE_z,'y')
    plot(zo,chi2MLE,'yo')
    xlabel('zo')
    ylabel('Log-Likelihood')
  
    subplot(2,4,5)
    hold on
    sigxrange=[sigx-.2:.01:sigx+.2]';
    for ii=1:length(sigxrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigxrange(ii))) + erf((0.5 + X - xo)./(sqrt(2).*sigxrange(ii))));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_sigx(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(sigxrange,chi2MLE_sigx,'k')
    plot(sigx,chi2MLE,'ko')
    xlabel('sigx')
    ylabel('Log-Likelihood')
    
    subplot(2,4,6)
    hold on
    sigyrange=[sigy-.2:.01:sigy+.2]';
    for ii=1:length(sigyrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigyrange(ii))) + erf((0.5 + Y - yo)./(sqrt(2).*sigyrange(ii))));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_sigy(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(sigyrange,chi2MLE_sigy,'k')
    plot(sigy,chi2MLE,'ko')
    xlabel('sigy')
    ylabel('Log-Likelihood')
    
    subplot(2,4,7)
    hold on
    sigzrange=[sigz-.2:.01:sigz+.2]';
    for ii=1:length(sigyrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigzrange(ii))) + erf((0.5 + Z - zo)./(sqrt(2).*sigzrange(ii))));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_sigz(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(sigzrange,chi2MLE_sigz,'k')
    plot(sigz,chi2MLE,'ko')
    xlabel('sigz')
    ylabel('Log-Likelihood')
    
    subplot(2,4,8)
    hold on
    bgrange=[sqrt(bg^2)-10:.01:sqrt(bg^2)+10]';
    for ii=1:length(bgrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
        u=sqrt(bgrange(ii)^2)+sqrt(Io^2).*Ex.*Ey.*Ez;
        chi2MLE_bg(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(bgrange,chi2MLE_bg,'c')
    plot(sqrt(bg^2),chi2MLE,'co')
    xlabel('bg')
    ylabel('Log-Likelihood')
end
end