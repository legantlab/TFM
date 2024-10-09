function [r,residual]=feature3dMLE(image,image_unfilt, diameter,masksz, xyzmax,inputv, sep, masscut, threshold)
%This program was originally written to find local maxima in the image and then perform a
%gaussian mask estimator fitting algorithm.  Now it finds the local
%maxima (unchanged from original code) and then performs gaussian maximum likelihood estimation of the bead centroids. 

%this program is written by Yongxiang Gao and Maria Kilfoil based on IDL
%code written by John C. Crocker and David G. Grier.
%varargin should contain input for [ separation, masscut, threshold], and
%whether there is input for it or not indicated in the logical input for
%inputv. inputv contains 3 logical input, for example [0,1,0].
%Last revisement is done on June 15, 2005
% 		Feature3d
% PURPOSE:
% Finds and measures roughly spheroidal 'features' within 
% 		a 3d image. Works best with dilute or separate feature,
% 		i.e. non close-packed structures.
% CATEGORY:
% 		Image Processing
%  CALLING SEQUENCE:
% 		f = feature3d( image, diameter ,separation, masscut, threshold )
%  INPUTS:
% 		image:	(nx,ny,nz) array which presumably contains some
% 			features worth finding
% 		diameter: a parameter which should be a little greater than
% 			the diameter of the largest features in the image.
% 			May be a single float number if the image is 
% 			isotropic, a 3-vector otherwise.
% 		separation: an optional parameter which specifies the 
% 			minimum allowable separation between feature 
% 			centers. The default value is diameter-1.
% 		masscut: Setting this parameter saves runtime by reducing
% 			the runtime wasted on low mass 'noise' features.
% 		threshold: Set this parameter to a number less than 1 to
% 			threshold each particle image by 
% 			(peak height)*(threshold).  Reduces pixel biasing
% 			with an particle specific threshold.
%  OUTPUTS:
%   r(:,1) = x position (in microns)
%   r(:,2) = y position (in microns)
%   r(:,3) = z position (in microns)
%   r(:,4) = Intensity (integrated area) of the fitted gaussian
%   r(:,5) = residual (in summed pixel counts)
% 
%  SIDE EFFECTS:
% 		Displays the number of features found on the screen.
%  RESTRICTIONS:
% 		To work properly, the image must consist of bright, 
% 		smooth regions on a roughly zero-valued background. 
% 		To find dark features, the image should be 
% 		inverted and the background subtracted. If the image
% 		contains a large amount of high spatial frequency noise,
% 		performance will be improved by first filtering the image.
% 		'bpass3d' will remove high spatial frequency noise, and 
% 		subtract the image background and thus may provides a useful 
% 		complement to using this program. Individual features 
% 		should NOT overlap or touch.  Furthermore, the maximum
% 		value of the top of the feature must be in the top 30th
% 		percentile of brightness in the entire image.
% 		For images where the particles are close packed, the
% 		system of bpass3d/feature3d is not ideal, but will give
% 		rough coordinates.  We often find setting 'sep' to roughly
% 	diameter/2 seems helpful to avoid particle loss.
%  PROCEDURE:
% 		First, identify the positions of all the local maxima in
% 		the image ( defined in a circular neighborhood with radius
% 		equal to 'separation' ). Then separate subimages around each bead and perform MLE fitting of the centroids.
% 
%  *********	       READ THE FOLLOWING IMPORTANT CAVEAT!        **********
% 		'feature3d' is capable of finding image features with sub-pixel
% 		accuracy, but only if used correctly- that is, if the 
% 		background is subtracted off properly and the centroid mask 
% 		is larger than the feature, so that clipping does not occur.
% 		It is an EXCELLENT idea when working with new data to plot
% 		a histogram of the x-positions mod 1, that is, of the
% 		fractional part of x in pixels.  If the resulting histogram
% 		is flat, then you're ok, if its strongly peaked, then you're
% 		doing something wrong- but probably still getting 'nearest
% 		pixel' accuracy.
% 
% 		For a more quantitative treatment of sub-pixel position 
% 		resolution see: 
% 		J.C. Crocker and D.G. Grier, J. Colloid Interface Sci.
% 		*179*, 298 (1996).
% 
%  MODIFICATION HISTORY:
% 		This code is inspired by feature_stats2 written by
% 			David G. Grier, U of Chicago, 			 1992.
% 		Generalized version of feature.pro			 1998.
% 		Improved local maximum routine, from fcp3d		 1999.
% 		Added 'threshold' keyword to reduce pixel biasing.	 1999.
% 		
% 	This code feature3d.pro is copyright 1999, John C. Crocker and 
% 	David G. Grier.  It should be considered 'freeware'- and may be
% 	distributed freely in its original form when properly attributed.
% -
% 
% 	produce a 3d, anisotropic parabolic mask
% 	anisotropic masks are 'referenced' to the x-axis scale.
% 	ratios less than one squash the thing relative to 'x'
% 	using float ratios allows precise 'diameter' settings 
% 	in an odd-sized mask.

% For anisotropy, make diameter a 3-vector, otherwise one # is ok.
% Image should consist of smooth well separated peaks on a zero
% or near zero background, with diameter set somewhat larger
% than the diameter of the peak!
if inputv(1)==1; separation=sep; end
if inputv(2)==0 ; masscut = 0; else; masscut=masscut; end
if (inputv(3)==1)
    threshold=threshold;
    if threshold <= 0 || threshold >= 0.9
        error('Threshold value must be between 0.0 and 0.9!');
    end
end

% make extents be the smallest odd integers bigger than diameter
if length(diameter) == 1 ;diameter = zeros(1,3)+diameter; end
extent = fix(diameter) + 1;
extent = extent + mod((extent+1),2);
extt=extent;
[nx, ny, nz] = size( image );
if not (inputv(1)==1)
       sep = diameter-1; 
   else
       sep = double(separation);
end
if size(sep) == 1 ; sep = zeros(1,3)+sep; end
% Put a border around the image to prevent mask out-of-bounds
a = zeros( nx+extent(1), ny+extent(2), nz+extent(3),'single' );
for i = 1:nz  % do as a loop to reduce memory piggage
a(fix(extent(1)/2)+1:fix((extent(1)/2))+nx,fix(extent(2)/2)+1:fix((extent(2)/2))+ny, fix(extent(3)/2)+i)=image(:,:,i);
end
nx = nx + extent(1); 
ny = ny + extent(2);
nz = nz + extent(3); 
% 	Find the local maxima in the filtered image
loc = llmx3dMB(a,sep,fix((extent/2)));
if loc(1)== -1 
disp('No features found!');		
loc(1)=-1;
end
%  	Set up some stuff....
nmax=length(loc(:,1));
x = loc(:,1)+1;
y = loc(:,2)+1;
z = loc(:,3)+1;
xs=(double(extent(1))+1.0)/2.0;
ys=(double(extent(2))+1.0)/2.0;
zs=(double(extent(3))+1.0)/2.0;

%Leave some space near the edge to avoid out of border
id=find(x>masksz(1)-2 & x<xyzmax(1)+extent(1)-masksz(1)+2 & y>masksz(2)-2 & y<xyzmax(2)+extent(2)-masksz(2)+2 & z>masksz(3)-2 & z<xyzmax(3)+extent(3)-masksz(3));
x=x(id); y=y(id); z=z(id);
extent = fix(masksz) + 1;
extent = extent + mod((extent+1),2);
xl = x - fix(extent(1)/2); 
xh = xl + extent(1) -1;
yl = y - fix(extent(2)/2); 
yh = yl + extent(2) -1;
zl = z - fix(extent(3)/2); 
zh = zl + extent(3) -1;

rsq   = lrsqd3dMB(extent,[1,1],masksz(2)/masksz(1),masksz(3)/masksz(1));
mask  = rsq < ((masksz(1)/2))^2 +1;
shell = ((mask) == (rsq > ((masksz(1)/2 -1))^2));
nask  = sum(mask(:));
rmask = (rsq.*mask)+(1/6);  

imask = mod(make_arr( extent(1), extent(2), extent(3))-1,extent(1)) + 1;
xmask = mask .* imask;
imask = mod(make_arr( extent(2), extent(1), extent(3))-1,extent(2)) + 1;
ymask = mask .*permute(imask,[2,1,3]);
imask = mod(make_arr( extent(3), extent(2), extent(1))-1,extent(3)) + 1;
zmask = mask .* permute(imask,[3,2,1]);
nmax=length(x);
m  = zeros(1,nmax);
pd = zeros(1,nmax);
thresh = zeros(1,nmax);
nthresh = zeros(1,nmax); 
for i=1:nmax
    tops(i)=a(x(i),y(i),z(i));
end
if inputv(3)==1
    thresh = tops*threshold;
end
if inputv(3)==1
    for i=1:nmax
        bb=a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i));
        nthresh(i)=sum(sum(sum(a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i)).*mask>thresh(i))))/sum(sum(sum(a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i)).*mask>0)));
    end
end

% Estimate the mass	
for i=1:nmax
    clear temp temp2
    temp=a(xl(i):xh(i),yl(i):yh(i),zl(i):zh(i))-thresh(i);
    m(i) = sum(sum(sum(((temp> 0).*temp) .* mask ))); %mass of each features
end
%do a masscut, and prevent divide by zeroes in the centroid calc.
w = find( m > masscut); %only those features with a total mass higher than masscut are considered to be a real feature
nmax=length(w);
disp(['Found ',num2str(nmax),' features'])
if nmax==0 
disp('No features found!');
nmax=-1;
end
xl = xl(w);
xh = xh(w);
yl = yl(w);
yh = yh(w);
zl = zl(w);
zh = zh(w);
x = x(w);
y = y(w);
z = z(w);
m = m(w);
tops = tops(w);
thresh = thresh(w);
nthresh = nthresh(w);
disp(strcat(num2str(nmax,'%01.0f'),' features found.'));

% Setup some result arrays
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
 a_unfilt = zeros( nx, ny, nz,'single' );
for i = 1:nz-fix(extent(3))  % do as a loop to reduce memory piggage
a_unfilt(fix(extent(1)/2)+1:nx-(fix(extent(1)/2)+1),fix(extent(2)/2)+1:ny-(fix((extent(2)/2))+1), fix(extent(3)/2)+i)=image_unfilt(:,:,i);
end
 subImage=double(a_unfilt(fix(xl(1)):fix(xh(1)), fix(yl(1)):fix(yh(1)), fix(zl(1)):fix(zh(1))));
 [subx,suby,subz]=size(subImage);
 t=1; 
 p=1;
 Io=zeros(nmax,1);
 xc=zeros(nmax,1);
 yc=zeros(nmax,1);
 zc=zeros(nmax,1);
 residual=zeros(nmax,1);
 %Perform MLE fitting of the bead centroids
 for i=1:nmax
     subImage=double(a_unfilt(fix(xl(i)):fix(xh(i)), fix(yl(i)):fix(yh(i)), fix(zl(i)):fix(zh(i))));
      [pn,residual(i),delta,u]=MLE3D_5Param_1Emitter(subImage,[max(subImage(:)),subx/2,suby/2,subz/2,1]);
                Io(p)=pn(1);
                xc(p)=fix(xl(i))-xs+pn(3);
                yc(p)=fix(yl(i))-ys+pn(2);
                zc(p)=fix(zl(i))-zs+pn(4);
      p=p+1;
 end
 

r=[xc(1:p-1),yc(1:p-1),zc(1:p-1),Io(1:p-1),residual(1:p-1)];
