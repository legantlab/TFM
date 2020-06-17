function [pn,CRLB,u,rChi,delta,convergeflag,pval,chi2MLE]=MLE_LevMarq_6Param_2013_11_12(PSF1,X,Y,p_initial,maxIter)
%A function to perform maximum likelihood fitting of single emitter spots.  This function fits the integrated intensity, the x-centroid, y-centroid, x-sigma,
%y-sigma, and per pixel background under a symmetric 2D Gaussian PSF model.

% Input Parameters
% PSF1 = a N X N X 1 array of the photon-converted image data to fit
% p_initial = a 1 X 6 vector containing the initial guesses for the fit parameters
% iters = the number of iterations for the Newton solver to take
% stepsize = a factor to scale the steps taken by the Newton solver

%Output Parameters
% pn = a 1 X 6 vector of the fitted parameters
%     pn(1) = the integrated intensity of the spot
%     pn(2) = the x-centroid of the spot
%     pn(3) = the y-centroid of the spot
%     pn(4) = the x-sigma of the spot
%     pn(5) = the y-sigma of the spot
%     pn(6) = the per pixel background of the spot
% r = the sum of the squared per-pixel differences between the image data and the model fit
% delta = an iters X 5 matrix displaying the amount that each variable is incrementated by in each iteration of the Newton solver
% flag = a flag to indicate algorithm convergence. If change in the value of all parameters is less than 1%, then flag = 1, else flag = 0

Io=p_initial(1); % initial guess for the integrated intensity of the spot
xo=p_initial(2); % initial guess for the x-centroid of the spot
yo=p_initial(3); % initial guess for the y-centroid of the spot
sigx=p_initial(4); % initial guess for the x-sigma of the spot
sigy=p_initial(5); % initial guess for the y-sigma of the spot
bg=p_initial(6); % initial guess for the per pixel background of the spot

%Initialize variables
rChi=0;
delta=zeros(6,maxIter);
diagnostics=1; %Turn this on to plot the fit data
funTol = 1e-5; %The minimum threshold to define how much the Chi-square value needs to decrease in order to take a step
convergeflag=0;
FSHR=zeros(6,6);
CRLB=zeros(6,1);

%Define image model
Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
u=sqrt(bg^2)+sqrt(Io^2)*Ex.*Ey;
chi2MLE=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));

for i=1:maxIter
    %Calculate 1st derivatives
    du1_dIo=Io/sqrt(Io^2).*Ex.*Ey;
    du1_dxo=(1./2).*sqrt(Io^2).*(sqrt(2./pi)./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx) - sqrt(2./pi)./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx)).*Ey;
    du1_dyo=(1./2).*sqrt(Io^2).*(sqrt(2./pi)./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy) - sqrt(2./pi)./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy)).*Ex;
    du1_dsigx=(1./2).*sqrt(Io^2).*((sqrt(2./pi).*(-0.5 + X - xo))./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^2) - (sqrt(2./pi).*(0.5 + X - xo))./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^2)).*Ey;
    du1_dsigy=(1./2).*sqrt(Io^2).*((sqrt(2./pi).*(-0.5 + Y - yo))./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^2) - (sqrt(2./pi).*(0.5 + Y - yo))./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^2)).*Ex;
    du1_dbg=bg/sqrt(bg^2);
    
    %Calculate Jacobian
    B_k(1,1)=-sum((1-PSF1(:)./u(:)).*du1_dIo(:));
    B_k(2,1)=-sum((1-PSF1(:)./u(:)).*du1_dxo(:));
    B_k(3,1)=-sum((1-PSF1(:)./u(:)).*du1_dyo(:));
    B_k(4,1)=-sum((1-PSF1(:)./u(:)).*du1_dsigx(:));
    B_k(5,1)=-sum((1-PSF1(:)./u(:)).*du1_dsigy(:));
    B_k(6,1)=-sum((1-PSF1(:)./u(:)).*du1_dbg(:));
    
    %Calculate Hessian - setting second derivatives to zeros
    lbda=0;
    A_kl(1,1)=sum(du1_dIo(:).*du1_dIo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(2,1)=sum(du1_dIo(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,1)=sum(du1_dIo(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,1)=sum(du1_dIo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(5,1)=sum(du1_dIo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(6,1)=sum(du1_dIo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,2)=sum(du1_dxo(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,2)=sum(du1_dxo(:).*du1_dxo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(3,2)=sum(du1_dxo(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,2)=sum(du1_dxo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(5,2)=sum(du1_dxo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(6,2)=sum(du1_dxo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,3)=sum(du1_dyo(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,3)=sum(du1_dyo(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,3)=sum(du1_dyo(:).*du1_dyo(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(4,3)=sum(du1_dyo(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(5,3)=sum(du1_dyo(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(6,3)=sum(du1_dyo(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,4)=sum(du1_dsigx(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,4)=sum(du1_dsigx(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,4)=sum(du1_dsigx(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,4)=sum(du1_dsigx(:).*du1_dsigx(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(5,4)=sum(du1_dsigx(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(6,4)=sum(du1_dsigx(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,5)=sum(du1_dsigy(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,5)=sum(du1_dsigy(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,5)=sum(du1_dsigy(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,5)=sum(du1_dsigy(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(5,5)=sum(du1_dsigy(:).*du1_dsigy(:).*PSF1(:)./u(:).^2)*(1+lbda);
    A_kl(6,5)=sum(du1_dsigy(:).*du1_dbg(:).*PSF1(:)./u(:).^2);
    
    A_kl(1,6)=sum(du1_dbg(:).*du1_dIo(:).*PSF1(:)./u(:).^2);
    A_kl(2,6)=sum(du1_dbg(:).*du1_dxo(:).*PSF1(:)./u(:).^2);
    A_kl(3,6)=sum(du1_dbg(:).*du1_dyo(:).*PSF1(:)./u(:).^2);
    A_kl(4,6)=sum(du1_dbg(:).*du1_dsigx(:).*PSF1(:)./u(:).^2);
    A_kl(5,6)=sum(du1_dbg(:).*du1_dsigy(:).*PSF1(:)./u(:).^2);
    A_kl(6,6)=sum(du1_dbg(:).*du1_dbg(:).*PSF1(:)./u(:).^2)*(1+lbda);
    
    delta(:,i)=pinv(A_kl)*B_k;
    
    Io_test=Io+delta(1,i);
    xo_test=xo+delta(2,i);
    yo_test=yo+delta(3,i);
    sigx_test=sigx+delta(4,i);
    sigy_test=sigy+delta(5,i);
    bg_test=bg+delta(6,i);
    
    %     Test whether update gives a more favorable solution as measured by a greater likelihood
    passed=0;
    t=1;
    while (~passed && t<=10 &&~convergeflag)
        %Recalculate the objective function at the proposed parameter location
        Ex_test=(1./2).*(-erf((-0.5 + X - xo_test)./(sqrt(2).*sigx_test)) + erf((0.5 + X - xo_test)./(sqrt(2).*sigx_test)));
        Ey_test=(1./2).*(-erf((-0.5 + Y - yo_test)./(sqrt(2).*sigy_test)) + erf((0.5 + Y - yo_test)./(sqrt(2).*sigy_test)));
        u_test=sqrt(bg_test^2)+sqrt(Io_test^2).*Ex_test.*Ey_test;
        chi2MLE_test=2*sum(u_test(:)-PSF1(:))-2.*sum(PSF1(:).*log(u_test(:)./PSF1(:)));
        passed=chi2MLE_test<=(chi2MLE+funTol);%Compare with the previous value
        %If the value doesn't decrease (or doesn't stay within funTol of the previous value - suggesting convergence), try a smaller step
        if (~passed && t<10)
            lbda=2^t;
            t=t+1;
            A_kl(1,1)=sum(du1_dIo(:).*du1_dIo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(2,2)=sum(du1_dxo(:).*du1_dxo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(3,3)=sum(du1_dyo(:).*du1_dyo(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(4,4)=sum(du1_dsigx(:).*du1_dsigx(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(5,5)=sum(du1_dsigy(:).*du1_dsigy(:).*PSF1(:)./u(:).^2)*(1+lbda);
            A_kl(6,6)=sum(du1_dbg(:).*du1_dbg(:).*PSF1(:)./u(:).^2)*(1+lbda);
            delta(:,i)=pinv(A_kl)*B_k;
            Io_test=Io+delta(1,i);
            xo_test=xo+delta(2,i);
            yo_test=yo+delta(3,i);
            sigx_test=sigx+delta(4,i);
            sigy_test=sigy+delta(5,i);
            bg_test=bg+delta(6,i);
        else if t>=10
                warning('The smallest step size has been reached, but algorithm has not converged')
                pn=[0,0,0,0,0,0,0];
                pval = -1;
                return
            else %Update parameter values
                Io=Io_test;
                xo=xo_test;
                yo=yo_test;
                sigx=sigx_test;
                sigy=sigy_test;
                bg=bg_test;
                Ex=Ex_test;
                Ey=Ey_test;
                u=u_test;
                convergeflag=((chi2MLE-chi2MLE_test)/chi2MLE)<funTol;
            end
        end
    end
    chi2MLE=chi2MLE_test;
end
pn=[sqrt(Io^2),xo,yo,sigx,sigy,sqrt(bg^2)];
%Calculate Fisher information matrix
FSHR(1,1)=sum(du1_dIo(:).*du1_dIo(:)./u(:));
FSHR(1,2)=sum(du1_dIo(:).*du1_dxo(:)./u(:));
FSHR(1,3)=sum(du1_dIo(:).*du1_dyo(:)./u(:));
FSHR(1,4)=sum(du1_dIo(:).*du1_dsigx(:)./u(:));
FSHR(1,5)=sum(du1_dIo(:).*du1_dsigy(:)./u(:));
FSHR(1,6)=sum(du1_dIo(:).*du1_dbg(:)./u(:));


FSHR(2,1)=FSHR(1,2);
FSHR(2,2)=sum(du1_dxo(:).*du1_dxo(:)./u(:));
FSHR(2,3)=sum(du1_dxo(:).*du1_dyo(:)./u(:));
FSHR(2,4)=sum(du1_dxo(:).*du1_dsigx(:)./u(:));
FSHR(2,5)=sum(du1_dxo(:).*du1_dsigy(:)./u(:));
FSHR(2,6)=sum(du1_dxo(:).*du1_dbg(:)./u(:));


FSHR(3,1)=FSHR(1,3);
FSHR(3,2)=FSHR(2,3);
FSHR(3,3)=sum(du1_dyo(:).*du1_dyo(:)./u(:));
FSHR(3,4)=sum(du1_dyo(:).*du1_dsigx(:)./u(:));
FSHR(3,5)=sum(du1_dyo(:).*du1_dsigy(:)./u(:));
FSHR(3,6)=sum(du1_dyo(:).*du1_dbg(:)./u(:));

FSHR(4,1)=FSHR(1,4);
FSHR(4,2)=FSHR(2,4);
FSHR(4,3)=FSHR(3,4);
FSHR(4,4)=sum(du1_dsigx(:).*du1_dsigx(:)./u(:));
FSHR(4,5)=sum(du1_dsigx(:).*du1_dsigy(:)./u(:));
FSHR(4,6)=sum(du1_dsigx(:).*du1_dbg(:)./u(:));

FSHR(5,1)=FSHR(1,5);
FSHR(5,2)=FSHR(2,5);
FSHR(5,3)=FSHR(3,5);
FSHR(5,4)=FSHR(4,5);
FSHR(5,5)=sum(du1_dsigy(:).*du1_dsigy(:)./u(:));
FSHR(5,6)=sum(du1_dsigy(:).*du1_dbg(:)./u(:));

FSHR(6,1)=FSHR(1,6);
FSHR(6,2)=FSHR(2,6);
FSHR(6,3)=FSHR(3,6);
FSHR(6,4)=FSHR(4,6);
FSHR(6,5)=FSHR(5,6);
FSHR(6,6)=sum(du1_dbg(:).*du1_dbg(:)./u(:));

%Calculate CRLB
FSHRinv=inv(FSHR);
CRLB=diag(FSHRinv);

rChi=sum((PSF1(:)-u(:)).^2./PSF1(:));
pval = 0;
%

%diagnostics provides graphical plots of the function behavior in the neighborhood of the final parameter values
if diagnostics==1
    figure
    %Check molecule fit compared to raw data
    Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
    Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
    u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey;
    subplot(1,2,1)
    hold on
    imagesc(u)
    scatter(pn(3),pn(2),'bo','filled')
    axis equal
    title('fit image')
    subplot(1,2,2)
    hold on
    imagesc(PSF1)
    scatter(pn(3),pn(2),'bo','filled')
    axis equal
    title('original image')
    Irange=[sqrt(Io^2)-1000:1:sqrt(Io^2)+1000]';
    figure
    subplot(3,2,1)
    hold on
    for ii=1:length(Irange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        u=sqrt(bg^2)+sqrt(Irange(ii)^2).*Ex.*Ey;
        chi2MLE_I(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(Irange,chi2MLE_I,'b')
    plot(sqrt(Io^2),chi2MLE,'bo')
    xlabel('Io')
    ylabel('Log-Likelihood')
    subplot(3,2,2)
    hold on
    xrange=[xo-2:.01:xo+2]';
    for ii=1:length(xrange)
        Ex=(1./2).*(-erf((-0.5 + X - xrange(ii))./(sqrt(2).*sigx)) + erf((0.5 + X - xrange(ii))./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey;
        chi2MLE_x(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(xrange,chi2MLE_x,'r')
    plot(xo,chi2MLE,'ro')
    xlabel('xo')
    ylabel('Log-Likelihood')
    subplot(3,2,3)
    hold on
    yrange=[yo-2:.01:yo+2]';
    for ii=1:length(yrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yrange(ii))./(sqrt(2).*sigy)) + erf((0.5 + Y - yrange(ii))./(sqrt(2).*sigy)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey;
        chi2MLE_y(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(yrange,chi2MLE_y,'g')
    plot(yo,chi2MLE,'go')
    xlabel('yo')
    ylabel('Log-Likelihood')
    subplot(3,2,4)
    hold on
    sigxrange=[sigx-.2:.01:sigx+.2]';
    for ii=1:length(sigxrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigxrange(ii))) + erf((0.5 + X - xo)./(sqrt(2).*sigxrange(ii))));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey;
        chi2MLE_sigx(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(sigxrange,chi2MLE_sigx,'k')
    plot(sigx,chi2MLE,'ko')
    xlabel('sigx')
    ylabel('Log-Likelihood')
    subplot(3,2,5)
    hold on
    sigyrange=[sigy-.2:.01:sigy+.2]';
    for ii=1:length(sigyrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigyrange(ii))) + erf((0.5 + Y - yo)./(sqrt(2).*sigyrange(ii))));
        u=sqrt(bg^2)+sqrt(Io^2).*Ex.*Ey;
        chi2MLE_sigy(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(sigyrange,chi2MLE_sigy,'k')
    plot(sigy,chi2MLE,'ko')
    xlabel('sigy')
    ylabel('Log-Likelihood')
    subplot(3,2,6)
    hold on
    bgrange=[sqrt(bg^2)-10:.01:sqrt(bg^2)+10]';
    for ii=1:length(bgrange)
        Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
        Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
        u=sqrt(bgrange(ii)^2)+sqrt(Io^2).*Ex.*Ey;
        chi2MLE_bg(ii)=2*sum(u(:)-PSF1(:))-2.*sum(PSF1(:).*log(u(:)./PSF1(:)));
    end
    plot(bgrange,chi2MLE_bg,'c')
    plot(sqrt(bg^2),chi2MLE,'co')
    xlabel('bg')
    ylabel('Log-Likelihood')
end
end