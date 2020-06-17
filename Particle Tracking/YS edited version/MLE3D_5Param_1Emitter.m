function [pn,r,delta,u]=MLE3D_5Param_1Emitter(PSF1,p_initial)
iters=60;
r=0;
% delta=0;
Io=p_initial(1);
xo=p_initial(2);
yo=p_initial(3);
zo=p_initial(4);
bg=p_initial(5);
sigx=1.25;
sigy=1.11;
sigz=1.07;

[dimx,dimy,dimz]=size(PSF1);
[X,Y,Z]=meshgrid([1:1:dimy],[1:1:dimx],[1:1:dimz]);

for i=1:iters
    
    Ex=(1./2).*(-erf((-0.5 + X - xo)./(sqrt(2).*sigx)) + erf((0.5 + X - xo)./(sqrt(2).*sigx)));
    Ey=(1./2).*(-erf((-0.5 + Y - yo)./(sqrt(2).*sigy)) + erf((0.5 + Y - yo)./(sqrt(2).*sigy)));
    Ez=(1./2).*(-erf((-0.5 + Z - zo)./(sqrt(2).*sigz)) + erf((0.5 + Z - zo)./(sqrt(2).*sigz)));
    u=bg+Io.*Ex.*Ey.*Ez;
    
%Calculate 1st derivatives
du1_dIo=Ex.*Ey.*Ez;
du1_dxo=(1./2).*Io.*(sqrt(2./pi)./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx) - sqrt(2./pi)./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx)).*Ey.*Ez;
du1_dyo=(1./2).*Io.*(sqrt(2./pi)./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy) - sqrt(2./pi)./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy)).*Ex.*Ez;
du1_dzo=(1./2).*Io.*(sqrt(2./pi)./(exp((-0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz) - sqrt(2./pi)./(exp((0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz)).*Ex.*Ey;
      
du1_dbg=ones(dimx,dimy,dimz);

%Calculate 2nd derivatives
du1_dIo2=zeros(dimx,dimy,dimz);
du1_dxo2=(1./2).*Io.*((sqrt(2./pi).*(-0.5 + X - xo))./(exp((-0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^3) - (sqrt(2./pi).*(0.5 + X - xo))./(exp((0.5 + X - xo).^2./(2.*sigx.^2)).*sigx.^3)).*Ey.*Ez;
du1_dyo2=(1./2).*Io.*((sqrt(2./pi).*(-0.5 + Y - yo))./(exp((-0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^3) - (sqrt(2./pi).*(0.5 + Y - yo))./(exp((0.5 + Y - yo).^2./(2.*sigy.^2)).*sigy.^3)).*Ex.*Ez;
du1_dzo2=(1./2).*Io.*((sqrt(2./pi).*(-0.5 + Z - zo))./(exp((-0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz.^3) - (sqrt(2./pi).*(0.5 + Z - zo))./(exp((0.5 + Z - zo).^2./(2.*sigz.^2)).*sigz.^3)).*Ex.*Ey;
du1_dbg2=zeros(dimx,dimy,dimz);

%Update parameters
Io=Io-.8*((sum(du1_dIo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dIo2(:).*(PSF1(:)./u(:)-1)-du1_dIo(:).^2.*(PSF1(:)./u(:).^2))));
xo=xo-.8*((sum(du1_dxo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dxo2(:).*(PSF1(:)./u(:)-1)-du1_dxo(:).^2.*(PSF1(:)./u(:).^2))));
yo=yo-.8*((sum(du1_dyo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dyo2(:).*(PSF1(:)./u(:)-1)-du1_dyo(:).^2.*(PSF1(:)./u(:).^2))));
zo=zo-.8*((sum(du1_dzo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dzo2(:).*(PSF1(:)./u(:)-1)-du1_dzo(:).^2.*(PSF1(:)./u(:).^2))));
bg=bg-.8*((sum(du1_dbg(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dbg2(:).*(PSF1(:)./u(:)-1)-du1_dbg(:).^2.*(PSF1(:)./u(:).^2))));
delta(i,:)=[((sum(du1_dIo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dIo2(:).*(PSF1(:)./u(:)-1)-du1_dIo(:).^2.*(PSF1(:)./u(:).^2))))/Io,((sum(du1_dxo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dxo2(:).*(PSF1(:)./u(:)-1)-du1_dxo(:).^2.*(PSF1(:)./u(:).^2))))/xo,((sum(du1_dyo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dyo2(:).*(PSF1(:)./u(:)-1)-du1_dyo(:).^2.*(PSF1(:)./u(:).^2))))/yo,((sum(du1_dzo(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dzo2(:).*(PSF1(:)./u(:)-1)-du1_dzo(:).^2.*(PSF1(:)./u(:).^2))))/zo,((sum(du1_dbg(:).*(PSF1(:)./u(:)-1)))/(sum(du1_dbg2(:).*(PSF1(:)./u(:)-1)-du1_dbg(:).^2.*(PSF1(:)./u(:).^2))))/bg];
end
pn=[Io,xo,yo,zo,bg];
r=sum(sqrt((PSF1(:)-u(:)).^2))/length(u(:));
end

