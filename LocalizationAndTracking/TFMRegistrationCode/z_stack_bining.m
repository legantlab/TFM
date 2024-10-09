z_width = 2;
z_min = floor(min(beadpos_sort(:,3)));
z_max = ceil(max(beadpos_sort(:,3)));
z_nbin = ceil((z_max - z_min)/z_width);
z_maxprime = z_min + z_width*z_nbin;
z_hist = zeros(1,z_nbin+1);
z_hist_bin = linspace(z_min,z_maxprime,z_nbin+1);
beadpos_index = zeros(length(beadpos_sort),3);
for i = 1:length(beadpos_sort)
    index = ceil((beadpos_sort(i,3) - z_min)/z_width);
    beadpos_index(i,:) = [beadpos_sort(i,1),beadpos_sort(i,2),index*z_width+z_min];
    z_hist(index) = z_hist(index)+1;
end



plot(z_hist_bin,z_hist);

xy_scatter = cell(1,z_nbin+1);
for i = 1:z_nbin+1
    xy_scatter{1,i} = beadpos_index(beadpos_index(:,3) == z_hist_bin(i),:);
end

