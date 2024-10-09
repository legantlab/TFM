function [r]=feature3dMB_parallel(image, diameter,masksz, xyzmax,numproc,jm,inputv, sep, masscut, threshold)
% r=rand(200,6);
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
% 		f(:,1):	this contains the x centroid positions, in pixels.
% 		f(:,2): this contains the y centroid positions, in pixels.
% 		f(:,3): this contains the z centroid positions, in pixels.
% 		f(:,4): this contains the integrated brightness.
% 		f(:,5): this contains the squared radius of gyration.
% 		f(:,6): this contains the peak height of the feature.
% 		f(:,7): this contains the fraction of voxels above threshold.
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
% 		equal to 'separation' ). Around each of these maxima, place a
% 		circular mask, of diameter 'diameter', and calculate the x,y,z
%     	centroids, the total of all the pixel values.
% 		If the restrictions above are adhered to, and the features
% 		are more than about 5 pixels across, the resulting x
% 		and y values will have errors of order 0.1 pixels for
% 		reasonably noise free images.
% 		If 'threshold' is set, then the image within the mask is
% 		thresholded to a value of the peak height*threshold.  This
% 		is useful when sphere images are connected by faint bridges
% 		which can cause pixel biasing.
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

% Ridiculous conditional to split up for parallel processing
if numproc==1
    [r]=refineCoords(xl,xh,yl,yh,zl,zh,x,y,z,m,tops,thresh,nthresh,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv);
end

if numproc==2
    splitIndexR=floor(nmax/2);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:end);
    xh_2 = xh(splitIndexR+1:end);
    yl_2 = yl(splitIndexR+1:end);
    yh_2 = yh(splitIndexR+1:end);
    zl_2 = zl(splitIndexR+1:end);
    zh_2 = zh(splitIndexR+1:end);
    x_2 = x(splitIndexR+1:end);
    y_2 = y(splitIndexR+1:end);
    z_2 = z(splitIndexR+1:end);
    m_2 = m(splitIndexR+1:end);
    tops_2 = tops(splitIndexR+1:end);
    thresh_2 = thresh(splitIndexR+1:end);
    nthresh_2 = nthresh(splitIndexR+1:end);
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    r=cat(1,results1{1},results2{1});
    destroy(jobs)
end

if numproc==3
    splitIndexR=floor(nmax/3);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:end);
    xh_3 = xh(2*splitIndexR+1:end);
    yl_3 = yl(2*splitIndexR+1:end);
    yh_3 = yh(2*splitIndexR+1:end);
    zl_3 = zl(2*splitIndexR+1:end);
    zh_3 = zh(2*splitIndexR+1:end);
    x_3 = x(2*splitIndexR+1:end);
    y_3 = y(2*splitIndexR+1:end);
    z_3 = z(2*splitIndexR+1:end);
    m_3 = m(2*splitIndexR+1:end);
    tops_3 = tops(2*splitIndexR+1:end);
    thresh_3 = thresh(2*splitIndexR+1:end);
    nthresh_3 = nthresh(2*splitIndexR+1:end);
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    r=cat(1,results1{1},results2{1},results3{1});
    destroy(jobs)
end

if numproc==4
    splitIndexR=floor(nmax/4);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:3*splitIndexR);
    xh_3 = xh(2*splitIndexR+1:3*splitIndexR);
    yl_3 = yl(2*splitIndexR+1:3*splitIndexR);
    yh_3 = yh(2*splitIndexR+1:3*splitIndexR);
    zl_3 = zl(2*splitIndexR+1:3*splitIndexR);
    zh_3 = zh(2*splitIndexR+1:3*splitIndexR);
    x_3 = x(2*splitIndexR+1:3*splitIndexR);
    y_3 = y(2*splitIndexR+1:3*splitIndexR);
    z_3 = z(2*splitIndexR+1:3*splitIndexR);
    m_3 = m(2*splitIndexR+1:3*splitIndexR);
    tops_3 = tops(2*splitIndexR+1:3*splitIndexR);
    thresh_3 = thresh(2*splitIndexR+1:3*splitIndexR);
    nthresh_3 = nthresh(2*splitIndexR+1:3*splitIndexR);
    
    xl_4 = xl(3*splitIndexR+1:end);
    xh_4 = xh(3*splitIndexR+1:end);
    yl_4 = yl(3*splitIndexR+1:end);
    yh_4 = yh(3*splitIndexR+1:end);
    zl_4 = zl(3*splitIndexR+1:end);
    zh_4 = zh(3*splitIndexR+1:end);
    x_4 = x(3*splitIndexR+1:end);
    y_4 = y(3*splitIndexR+1:end);
    z_4 = z(3*splitIndexR+1:end);
    m_4 = m(3*splitIndexR+1:end);
    tops_4 = tops(3*splitIndexR+1:end);
    thresh_4 = thresh(3*splitIndexR+1:end);
    nthresh_4 = nthresh(3*splitIndexR+1:end);
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    job4 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task4=createTask(job4, @refineCoords, 1, {xl_4,xh_4,yl_4,yh_4,zl_4,zh_4,x_4,y_4,z_4,m_4,tops_4,thresh_4,nthresh_4,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    submit(job4)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    waitForState(job4,'finished')
    results4 = getAllOutputArguments(job4);
    r=cat(1,results1{1},results2{1},results3{1},results4{1});
    destroy(jobs)
end

if numproc==5
    splitIndexR=floor(nmax/5);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:3*splitIndexR);
    xh_3 = xh(2*splitIndexR+1:3*splitIndexR);
    yl_3 = yl(2*splitIndexR+1:3*splitIndexR);
    yh_3 = yh(2*splitIndexR+1:3*splitIndexR);
    zl_3 = zl(2*splitIndexR+1:3*splitIndexR);
    zh_3 = zh(2*splitIndexR+1:3*splitIndexR);
    x_3 = x(2*splitIndexR+1:3*splitIndexR);
    y_3 = y(2*splitIndexR+1:3*splitIndexR);
    z_3 = z(2*splitIndexR+1:3*splitIndexR);
    m_3 = m(2*splitIndexR+1:3*splitIndexR);
    tops_3 = tops(2*splitIndexR+1:3*splitIndexR);
    thresh_3 = thresh(2*splitIndexR+1:3*splitIndexR);
    nthresh_3 = nthresh(2*splitIndexR+1:3*splitIndexR);
    
    xl_4 = xl(3*splitIndexR+1:4*splitIndexR);
    xh_4 = xh(3*splitIndexR+1:4*splitIndexR);
    yl_4 = yl(3*splitIndexR+1:4*splitIndexR);
    yh_4 = yh(3*splitIndexR+1:4*splitIndexR);
    zl_4 = zl(3*splitIndexR+1:4*splitIndexR);
    zh_4 = zh(3*splitIndexR+1:4*splitIndexR);
    x_4 = x(3*splitIndexR+1:4*splitIndexR);
    y_4 = y(3*splitIndexR+1:4*splitIndexR);
    z_4 = z(3*splitIndexR+1:4*splitIndexR);
    m_4 = m(3*splitIndexR+1:4*splitIndexR);
    tops_4 = tops(3*splitIndexR+1:4*splitIndexR);
    thresh_4 = thresh(3*splitIndexR+1:4*splitIndexR);
    nthresh_4 = nthresh(3*splitIndexR+1:4*splitIndexR);
    
    xl_5 = xl(4*splitIndexR+1:end);
    xh_5 = xh(4*splitIndexR+1:end);
    yl_5 = yl(4*splitIndexR+1:end);
    yh_5 = yh(4*splitIndexR+1:end);
    zl_5 = zl(4*splitIndexR+1:end);
    zh_5 = zh(4*splitIndexR+1:end);
    x_5 = x(4*splitIndexR+1:end);
    y_5 = y(4*splitIndexR+1:end);
    z_5 = z(4*splitIndexR+1:end);
    m_5 = m(4*splitIndexR+1:end);
    tops_5 = tops(4*splitIndexR+1:end);
    thresh_5 = thresh(4*splitIndexR+1:end);
    nthresh_5 = nthresh(4*splitIndexR+1:end);
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    job4 = createJob(jm, 'Configuration', 'WRL Config');
    job5 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task4=createTask(job4, @refineCoords, 1, {xl_4,xh_4,yl_4,yh_4,zl_4,zh_4,x_4,y_4,z_4,m_4,tops_4,thresh_4,nthresh_4,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task5=createTask(job5, @refineCoords, 1, {xl_5,xh_5,yl_5,yh_5,zl_5,zh_5,x_5,y_5,z_5,m_5,tops_5,thresh_5,nthresh_5,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    submit(job4)
    submit(job5)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    waitForState(job4,'finished')
    results4 = getAllOutputArguments(job4);
    waitForState(job5,'finished')
    results5 = getAllOutputArguments(job5);
    r=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1});
    destroy(jobs)
end

if numproc==6
    splitIndexR=floor(nmax/6);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:3*splitIndexR);
    xh_3 = xh(2*splitIndexR+1:3*splitIndexR);
    yl_3 = yl(2*splitIndexR+1:3*splitIndexR);
    yh_3 = yh(2*splitIndexR+1:3*splitIndexR);
    zl_3 = zl(2*splitIndexR+1:3*splitIndexR);
    zh_3 = zh(2*splitIndexR+1:3*splitIndexR);
    x_3 = x(2*splitIndexR+1:3*splitIndexR);
    y_3 = y(2*splitIndexR+1:3*splitIndexR);
    z_3 = z(2*splitIndexR+1:3*splitIndexR);
    m_3 = m(2*splitIndexR+1:3*splitIndexR);
    tops_3 = tops(2*splitIndexR+1:3*splitIndexR);
    thresh_3 = thresh(2*splitIndexR+1:3*splitIndexR);
    nthresh_3 = nthresh(2*splitIndexR+1:3*splitIndexR);
    
    xl_4 = xl(3*splitIndexR+1:4*splitIndexR);
    xh_4 = xh(3*splitIndexR+1:4*splitIndexR);
    yl_4 = yl(3*splitIndexR+1:4*splitIndexR);
    yh_4 = yh(3*splitIndexR+1:4*splitIndexR);
    zl_4 = zl(3*splitIndexR+1:4*splitIndexR);
    zh_4 = zh(3*splitIndexR+1:4*splitIndexR);
    x_4 = x(3*splitIndexR+1:4*splitIndexR);
    y_4 = y(3*splitIndexR+1:4*splitIndexR);
    z_4 = z(3*splitIndexR+1:4*splitIndexR);
    m_4 = m(3*splitIndexR+1:4*splitIndexR);
    tops_4 = tops(3*splitIndexR+1:4*splitIndexR);
    thresh_4 = thresh(3*splitIndexR+1:4*splitIndexR);
    nthresh_4 = nthresh(3*splitIndexR+1:4*splitIndexR);
    
    xl_5 = xl(4*splitIndexR+1:5*splitIndexR);
    xh_5 = xh(4*splitIndexR+1:5*splitIndexR);
    yl_5 = yl(4*splitIndexR+1:5*splitIndexR);
    yh_5 = yh(4*splitIndexR+1:5*splitIndexR);
    zl_5 = zl(4*splitIndexR+1:5*splitIndexR);
    zh_5 = zh(4*splitIndexR+1:5*splitIndexR);
    x_5 = x(4*splitIndexR+1:5*splitIndexR);
    y_5 = y(4*splitIndexR+1:5*splitIndexR);
    z_5 = z(4*splitIndexR+1:5*splitIndexR);
    m_5 = m(4*splitIndexR+1:5*splitIndexR);
    tops_5 = tops(4*splitIndexR+1:5*splitIndexR);
    thresh_5 = thresh(4*splitIndexR+1:5*splitIndexR);
    nthresh_5 = nthresh(4*splitIndexR+1:5*splitIndexR);
    
    xl_6 = xl(5*splitIndexR+1:end);
    xh_6 = xh(5*splitIndexR+1:end);
    yl_6 = yl(5*splitIndexR+1:end);
    yh_6 = yh(5*splitIndexR+1:end);
    zl_6 = zl(5*splitIndexR+1:end);
    zh_6 = zh(5*splitIndexR+1:end);
    x_6 = x(5*splitIndexR+1:end);
    y_6 = y(5*splitIndexR+1:end);
    z_6 = z(5*splitIndexR+1:end);
    m_6 = m(5*splitIndexR+1:end);
    tops_6 = tops(5*splitIndexR+1:end);
    thresh_6 = thresh(5*splitIndexR+1:end);
    nthresh_6 = nthresh(5*splitIndexR+1:end);
    
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    job4 = createJob(jm, 'Configuration', 'WRL Config');
    job5 = createJob(jm, 'Configuration', 'WRL Config');
    job6 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    % set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task4=createTask(job4, @refineCoords, 1, {xl_4,xh_4,yl_4,yh_4,zl_4,zh_4,x_4,y_4,z_4,m_4,tops_4,thresh_4,nthresh_4,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task5=createTask(job5, @refineCoords, 1, {xl_5,xh_5,yl_5,yh_5,zl_5,zh_5,x_5,y_5,z_5,m_5,tops_5,thresh_5,nthresh_5,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task6=createTask(job6, @refineCoords, 1, {xl_6,xh_6,yl_6,yh_6,zl_6,zh_6,x_6,y_6,z_6,m_6,tops_6,thresh_6,nthresh_6,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    submit(job4)
    submit(job5)
    submit(job6)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    waitForState(job4,'finished')
    results4 = getAllOutputArguments(job4);
    waitForState(job5,'finished')
    results5 = getAllOutputArguments(job5);
    waitForState(job6,'finished')
    results6 = getAllOutputArguments(job6);
    r=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1});
    destroy(jobs)
end

if numproc==7
    splitIndexR=floor(nmax/7);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:3*splitIndexR);
    xh_3 = xh(2*splitIndexR+1:3*splitIndexR);
    yl_3 = yl(2*splitIndexR+1:3*splitIndexR);
    yh_3 = yh(2*splitIndexR+1:3*splitIndexR);
    zl_3 = zl(2*splitIndexR+1:3*splitIndexR);
    zh_3 = zh(2*splitIndexR+1:3*splitIndexR);
    x_3 = x(2*splitIndexR+1:3*splitIndexR);
    y_3 = y(2*splitIndexR+1:3*splitIndexR);
    z_3 = z(2*splitIndexR+1:3*splitIndexR);
    m_3 = m(2*splitIndexR+1:3*splitIndexR);
    tops_3 = tops(2*splitIndexR+1:3*splitIndexR);
    thresh_3 = thresh(2*splitIndexR+1:3*splitIndexR);
    nthresh_3 = nthresh(2*splitIndexR+1:3*splitIndexR);
    
    xl_4 = xl(3*splitIndexR+1:4*splitIndexR);
    xh_4 = xh(3*splitIndexR+1:4*splitIndexR);
    yl_4 = yl(3*splitIndexR+1:4*splitIndexR);
    yh_4 = yh(3*splitIndexR+1:4*splitIndexR);
    zl_4 = zl(3*splitIndexR+1:4*splitIndexR);
    zh_4 = zh(3*splitIndexR+1:4*splitIndexR);
    x_4 = x(3*splitIndexR+1:4*splitIndexR);
    y_4 = y(3*splitIndexR+1:4*splitIndexR);
    z_4 = z(3*splitIndexR+1:4*splitIndexR);
    m_4 = m(3*splitIndexR+1:4*splitIndexR);
    tops_4 = tops(3*splitIndexR+1:4*splitIndexR);
    thresh_4 = thresh(3*splitIndexR+1:4*splitIndexR);
    nthresh_4 = nthresh(3*splitIndexR+1:4*splitIndexR);
    
    xl_5 = xl(4*splitIndexR+1:5*splitIndexR);
    xh_5 = xh(4*splitIndexR+1:5*splitIndexR);
    yl_5 = yl(4*splitIndexR+1:5*splitIndexR);
    yh_5 = yh(4*splitIndexR+1:5*splitIndexR);
    zl_5 = zl(4*splitIndexR+1:5*splitIndexR);
    zh_5 = zh(4*splitIndexR+1:5*splitIndexR);
    x_5 = x(4*splitIndexR+1:5*splitIndexR);
    y_5 = y(4*splitIndexR+1:5*splitIndexR);
    z_5 = z(4*splitIndexR+1:5*splitIndexR);
    m_5 = m(4*splitIndexR+1:5*splitIndexR);
    tops_5 = tops(4*splitIndexR+1:5*splitIndexR);
    thresh_5 = thresh(4*splitIndexR+1:5*splitIndexR);
    nthresh_5 = nthresh(4*splitIndexR+1:5*splitIndexR);
    
    xl_6 = xl(5*splitIndexR+1:6*splitIndexR);
    xh_6 = xh(5*splitIndexR+1:6*splitIndexR);
    yl_6 = yl(5*splitIndexR+1:6*splitIndexR);
    yh_6 = yh(5*splitIndexR+1:6*splitIndexR);
    zl_6 = zl(5*splitIndexR+1:6*splitIndexR);
    zh_6 = zh(5*splitIndexR+1:6*splitIndexR);
    x_6 = x(5*splitIndexR+1:6*splitIndexR);
    y_6 = y(5*splitIndexR+1:6*splitIndexR);
    z_6 = z(5*splitIndexR+1:6*splitIndexR);
    m_6 = m(5*splitIndexR+1:6*splitIndexR);
    tops_6 = tops(5*splitIndexR+1:6*splitIndexR);
    thresh_6 = thresh(5*splitIndexR+1:6*splitIndexR);
    nthresh_6 = nthresh(5*splitIndexR+1:6*splitIndexR);
    
    xl_7 = xl(6*splitIndexR+1:end);
    xh_7 = xh(6*splitIndexR+1:end);
    yl_7 = yl(6*splitIndexR+1:end);
    yh_7 = yh(6*splitIndexR+1:end);
    zl_7 = zl(6*splitIndexR+1:end);
    zh_7 = zh(6*splitIndexR+1:end);
    x_7 = x(6*splitIndexR+1:end);
    y_7 = y(6*splitIndexR+1:end);
    z_7 = z(6*splitIndexR+1:end);
    m_7 = m(6*splitIndexR+1:end);
    tops_7 = tops(6*splitIndexR+1:end);
    thresh_7 = thresh(6*splitIndexR+1:end);
    nthresh_7 = nthresh(6*splitIndexR+1:end);
    
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    job4 = createJob(jm, 'Configuration', 'WRL Config');
    job5 = createJob(jm, 'Configuration', 'WRL Config');
    job6 = createJob(jm, 'Configuration', 'WRL Config');
    job7 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    % set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task4=createTask(job4, @refineCoords, 1, {xl_4,xh_4,yl_4,yh_4,zl_4,zh_4,x_4,y_4,z_4,m_4,tops_4,thresh_4,nthresh_4,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task5=createTask(job5, @refineCoords, 1, {xl_5,xh_5,yl_5,yh_5,zl_5,zh_5,x_5,y_5,z_5,m_5,tops_5,thresh_5,nthresh_5,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task6=createTask(job6, @refineCoords, 1, {xl_6,xh_6,yl_6,yh_6,zl_6,zh_6,x_6,y_6,z_6,m_6,tops_6,thresh_6,nthresh_6,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task7=createTask(job7, @refineCoords, 1, {xl_7,xh_7,yl_7,yh_7,zl_7,zh_7,x_7,y_7,z_7,m_7,tops_7,thresh_7,nthresh_7,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    submit(job4)
    submit(job5)
    submit(job6)
    submit(job7)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    waitForState(job4,'finished')
    results4 = getAllOutputArguments(job4);
    waitForState(job5,'finished')
    results5 = getAllOutputArguments(job5);
    waitForState(job6,'finished')
    results6 = getAllOutputArguments(job6);
    waitForState(job7,'finished')
    results7 = getAllOutputArguments(job7);
    r=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1});
    destroy(jobs)
end

if numproc==8
    splitIndexR=floor(nmax/8);
    xl_1 = xl(1:splitIndexR);
    xh_1 = xh(1:splitIndexR);
    yl_1 = yl(1:splitIndexR);
    yh_1 = yh(1:splitIndexR);
    zl_1 = zl(1:splitIndexR);
    zh_1 = zh(1:splitIndexR);
    x_1 = x(1:splitIndexR);
    y_1 = y(1:splitIndexR);
    z_1 = z(1:splitIndexR);
    m_1 = m(1:splitIndexR);
    tops_1 = tops(1:splitIndexR);
    thresh_1 = thresh(1:splitIndexR);
    nthresh_1 = nthresh(1:splitIndexR);
    
    xl_2 = xl(splitIndexR+1:2*splitIndexR);
    xh_2 = xh(splitIndexR+1:2*splitIndexR);
    yl_2 = yl(splitIndexR+1:2*splitIndexR);
    yh_2 = yh(splitIndexR+1:2*splitIndexR);
    zl_2 = zl(splitIndexR+1:2*splitIndexR);
    zh_2 = zh(splitIndexR+1:2*splitIndexR);
    x_2 = x(splitIndexR+1:2*splitIndexR);
    y_2 = y(splitIndexR+1:2*splitIndexR);
    z_2 = z(splitIndexR+1:2*splitIndexR);
    m_2 = m(splitIndexR+1:2*splitIndexR);
    tops_2 = tops(splitIndexR+1:2*splitIndexR);
    thresh_2 = thresh(splitIndexR+1:2*splitIndexR);
    nthresh_2 = nthresh(splitIndexR+1:2*splitIndexR);
    
    xl_3 = xl(2*splitIndexR+1:3*splitIndexR);
    xh_3 = xh(2*splitIndexR+1:3*splitIndexR);
    yl_3 = yl(2*splitIndexR+1:3*splitIndexR);
    yh_3 = yh(2*splitIndexR+1:3*splitIndexR);
    zl_3 = zl(2*splitIndexR+1:3*splitIndexR);
    zh_3 = zh(2*splitIndexR+1:3*splitIndexR);
    x_3 = x(2*splitIndexR+1:3*splitIndexR);
    y_3 = y(2*splitIndexR+1:3*splitIndexR);
    z_3 = z(2*splitIndexR+1:3*splitIndexR);
    m_3 = m(2*splitIndexR+1:3*splitIndexR);
    tops_3 = tops(2*splitIndexR+1:3*splitIndexR);
    thresh_3 = thresh(2*splitIndexR+1:3*splitIndexR);
    nthresh_3 = nthresh(2*splitIndexR+1:3*splitIndexR);
    
    xl_4 = xl(3*splitIndexR+1:4*splitIndexR);
    xh_4 = xh(3*splitIndexR+1:4*splitIndexR);
    yl_4 = yl(3*splitIndexR+1:4*splitIndexR);
    yh_4 = yh(3*splitIndexR+1:4*splitIndexR);
    zl_4 = zl(3*splitIndexR+1:4*splitIndexR);
    zh_4 = zh(3*splitIndexR+1:4*splitIndexR);
    x_4 = x(3*splitIndexR+1:4*splitIndexR);
    y_4 = y(3*splitIndexR+1:4*splitIndexR);
    z_4 = z(3*splitIndexR+1:4*splitIndexR);
    m_4 = m(3*splitIndexR+1:4*splitIndexR);
    tops_4 = tops(3*splitIndexR+1:4*splitIndexR);
    thresh_4 = thresh(3*splitIndexR+1:4*splitIndexR);
    nthresh_4 = nthresh(3*splitIndexR+1:4*splitIndexR);
    
    xl_5 = xl(4*splitIndexR+1:5*splitIndexR);
    xh_5 = xh(4*splitIndexR+1:5*splitIndexR);
    yl_5 = yl(4*splitIndexR+1:5*splitIndexR);
    yh_5 = yh(4*splitIndexR+1:5*splitIndexR);
    zl_5 = zl(4*splitIndexR+1:5*splitIndexR);
    zh_5 = zh(4*splitIndexR+1:5*splitIndexR);
    x_5 = x(4*splitIndexR+1:5*splitIndexR);
    y_5 = y(4*splitIndexR+1:5*splitIndexR);
    z_5 = z(4*splitIndexR+1:5*splitIndexR);
    m_5 = m(4*splitIndexR+1:5*splitIndexR);
    tops_5 = tops(4*splitIndexR+1:5*splitIndexR);
    thresh_5 = thresh(4*splitIndexR+1:5*splitIndexR);
    nthresh_5 = nthresh(4*splitIndexR+1:5*splitIndexR);
    
    xl_6 = xl(5*splitIndexR+1:6*splitIndexR);
    xh_6 = xh(5*splitIndexR+1:6*splitIndexR);
    yl_6 = yl(5*splitIndexR+1:6*splitIndexR);
    yh_6 = yh(5*splitIndexR+1:6*splitIndexR);
    zl_6 = zl(5*splitIndexR+1:6*splitIndexR);
    zh_6 = zh(5*splitIndexR+1:6*splitIndexR);
    x_6 = x(5*splitIndexR+1:6*splitIndexR);
    y_6 = y(5*splitIndexR+1:6*splitIndexR);
    z_6 = z(5*splitIndexR+1:6*splitIndexR);
    m_6 = m(5*splitIndexR+1:6*splitIndexR);
    tops_6 = tops(5*splitIndexR+1:6*splitIndexR);
    thresh_6 = thresh(5*splitIndexR+1:6*splitIndexR);
    nthresh_6 = nthresh(5*splitIndexR+1:6*splitIndexR);
    
    xl_7 = xl(6*splitIndexR+1:7*splitIndexR);
    xh_7 = xh(6*splitIndexR+1:7*splitIndexR);
    yl_7 = yl(6*splitIndexR+1:7*splitIndexR);
    yh_7 = yh(6*splitIndexR+1:7*splitIndexR);
    zl_7 = zl(6*splitIndexR+1:7*splitIndexR);
    zh_7 = zh(6*splitIndexR+1:7*splitIndexR);
    x_7 = x(6*splitIndexR+1:7*splitIndexR);
    y_7 = y(6*splitIndexR+1:7*splitIndexR);
    z_7 = z(6*splitIndexR+1:7*splitIndexR);
    m_7 = m(6*splitIndexR+1:7*splitIndexR);
    tops_7 = tops(6*splitIndexR+1:7*splitIndexR);
    thresh_7 = thresh(6*splitIndexR+1:7*splitIndexR);
    nthresh_7 = nthresh(6*splitIndexR+1:7*splitIndexR);
    
    xl_8 = xl(7*splitIndexR+1:end);
    xh_8 = xh(7*splitIndexR+1:end);
    yl_8 = yl(7*splitIndexR+1:end);
    yh_8 = yh(7*splitIndexR+1:end);
    zl_8 = zl(7*splitIndexR+1:end);
    zh_8 = zh(7*splitIndexR+1:end);
    x_8 = x(7*splitIndexR+1:end);
    y_8 = y(7*splitIndexR+1:end);
    z_8 = z(7*splitIndexR+1:end);
    m_8 = m(7*splitIndexR+1:end);
    tops_8 = tops(7*splitIndexR+1:end);
    thresh_8 = thresh(7*splitIndexR+1:end);
    nthresh_8 = nthresh(7*splitIndexR+1:end);
    
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    job3 = createJob(jm, 'Configuration', 'WRL Config');
    job4 = createJob(jm, 'Configuration', 'WRL Config');
    job5 = createJob(jm, 'Configuration', 'WRL Config');
    job6 = createJob(jm, 'Configuration', 'WRL Config');
    job7 = createJob(jm, 'Configuration', 'WRL Config');
    job8 = createJob(jm, 'Configuration', 'WRL Config');
    
    jobs=get(jm,'jobs');
    % set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});
    
    task1=createTask(job1, @refineCoords, 1, {xl_1,xh_1,yl_1,yh_1,zl_1,zh_1,x_1,y_1,z_1,m_1,tops_1,thresh_1,nthresh_1,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task2=createTask(job2, @refineCoords, 1, {xl_2,xh_2,yl_2,yh_2,zl_2,zh_2,x_2,y_2,z_2,m_2,tops_2,thresh_2,nthresh_2,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task3=createTask(job3, @refineCoords, 1, {xl_3,xh_3,yl_3,yh_3,zl_3,zh_3,x_3,y_3,z_3,m_3,tops_3,thresh_3,nthresh_3,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task4=createTask(job4, @refineCoords, 1, {xl_4,xh_4,yl_4,yh_4,zl_4,zh_4,x_4,y_4,z_4,m_4,tops_4,thresh_4,nthresh_4,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task5=createTask(job5, @refineCoords, 1, {xl_5,xh_5,yl_5,yh_5,zl_5,zh_5,x_5,y_5,z_5,m_5,tops_5,thresh_5,nthresh_5,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task6=createTask(job6, @refineCoords, 1, {xl_6,xh_6,yl_6,yh_6,zl_6,zh_6,x_6,y_6,z_6,m_6,tops_6,thresh_6,nthresh_6,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task7=createTask(job7, @refineCoords, 1, {xl_7,xh_7,yl_7,yh_7,zl_7,zh_7,x_7,y_7,z_7,m_7,tops_7,thresh_7,nthresh_7,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    task8=createTask(job8, @refineCoords, 1, {xl_8,xh_8,yl_8,yh_8,zl_8,zh_8,x_8,y_8,z_8,m_8,tops_8,thresh_8,nthresh_8,a,rmask,xmask,ymask,zmask,extent,mask,extt,inputv});
    
    submit(job1)
    submit(job2)
    submit(job3)
    submit(job4)
    submit(job5)
    submit(job6)
    submit(job7)
    submit(job8)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    waitForState(job4,'finished')
    results4 = getAllOutputArguments(job4);
    waitForState(job5,'finished')
    results5 = getAllOutputArguments(job5);
    waitForState(job6,'finished')
    results6 = getAllOutputArguments(job6);
    waitForState(job7,'finished')
    results7 = getAllOutputArguments(job7);
    waitForState(job8,'finished')
    results8 = getAllOutputArguments(job8);
    r=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1},results8{1});
    destroy(jobs)
end

meanMass=mean(r(:,4));
stdMass=std(r(:,4));
r=r(abs(r(:,4)-meanMass)<(2*stdMass),:);