function cord = gen_cord(ysize, xsize)
[X,Y] = meshgrid(1:xsize,1:ysize);
cord = zeros(ysize,xsize,2);
cord(:,:,1) = X;
cord(:,:,2) = Y;
end