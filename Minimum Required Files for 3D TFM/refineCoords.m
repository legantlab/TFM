function [r]=refineCoords(xl,xh,yl,yh,zl,zh,x,y,z,m,tops,thresh,nthresh,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv)
%This function was part of feature3dMB, but has been split up to allow for
%parallel implementation.  See feature3DMB for discriptions of all input
% and output arguments.

% Setup some result arrays
nmax=length(xl);
xc = zeros(1,nmax);
yc = zeros(1,nmax);
zc = zeros(1,nmax);
rg = zeros(1,nmax);


% Calculate the radius of gyration^2
for i=1:nmax
    clear temp
    temp=a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i)) - thresh(i);
    rg(i) = sum(sum(sum(((temp > 0).*temp) .* rmask )))/m(i);
end
% Calculate peak centroids
for i=1:nmax
    clear temp
    temp=a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i)) - thresh(i);
	xc(i) = sum(sum(sum( ((temp >0).*temp) .* xmask )));
	yc(i) = sum(sum(sum( ((temp >0).*temp) .* ymask )));  %before we have yc and xc exchanged
	zc(i) = sum(sum(sum( ((temp >0).*temp) .* zmask )));
end
% Correct for the 'offset' of the centroid masks
xc = xc ./ m - ((double(extent(1))+1.0)/2.0);
yc = yc ./ m - ((double(extent(2))+1.0)/2.0);
zc = zc ./ m - ((double(extent(3))+1.0)/2.0);
% 	Update the positions and correct for the width of the 'border'
x = x + xc';
y = y + yc';
z = z + zc';

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    xcn=xc;
    ycn=yc;
    zcn=zc;
% do 20 iteration of fracshift; Can be more or less;
for j=1:15
  for i=1:nmax
    suba(:,:,:,i)=fracshift3dMB(double(a(fix(xl(i)):fix(xh(i)), fix(yl(i)):fix(yh(i)), fix(zl(i)):fix(zh(i)))), -xcn(i), -ycn(i),-zcn(i));
    m(i)=sum(sum(sum(suba(:,:,:,i).*mask)));
    rg(i) = sum(sum(sum(suba(:,:,:,i) .* rmask )))/m(i);
  end
  for i=1:nmax
    xc(i)=sum(sum(sum(suba(:,:,:,i).*xmask)));
    yc(i)=sum(sum(sum(suba(:,:,:,i).*ymask)));
    zc(i)=sum(sum(sum(suba(:,:,:,i).*zmask)));
  end
xc = xc./ m - ((double(extent(1))+1.0)/2.0);
yc = yc./ m - ((double(extent(2))+1.0)/2.0);
zc = zc./ m - ((double(extent(3))+1.0)/2.0);
xcn=xc+xcn;
ycn=yc+ycn;
zcn=zc+zcn;
x = x + xc' ;
y = y + yc' ;
z = z + zc' ;
end

x = x - fix(extt(1)/2);
y = y - fix(extt(2)/2);
z = z - fix(extt(3)/2);

if inputv(3)==1 
r=[x,y,z,m',rg',tops',nthresh'];
else
    r=[x,y,z,m',rg',tops'];
end