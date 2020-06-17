Bead_coord = displacedPtsMat(1:3:end,1:3);
u_x=u(1:3:end);
u_y=u(2:3:end);
u_z=u(3:3:end);
n_bead = length(Bead_coord);
max_dist = 30;
min_dist = 3;
nbin = 10;
interval = (max_dist-min_dist)/(nbin-1);
dist_vec = linspace(min_dist,max_dist,nbin);
dist_count = zeros(nbin,1);
dist_disp = zeros(nbin,1);
for i = 1:n_bead
    coord = Bead_coord(i,1:3);
    dist = sqrt((pointforce_coord(1)-coord(1))^2+(pointforce_coord(2)-coord(2))^2+(pointforce_coord(3)-coord(3))^2);
    if dist < max_dist
        index = 1+ceil((dist - min_dist)/interval);
        dist_count(index) = dist_count(index)+1;
        dist_disp(index) = dist_disp(index)+sqrt(u_x(i)^2+u_y(i)^2+u_z(i)^2);
    end
end
dist_disp = dist_disp./dist_count;

figure
plot(dist_disp,dist_vec)
xlabel('Distance from load (um)')
ylabel('Averaged Displacement (um)')
title('Displacement vs Distance')


figure
plot(dist_count,dist_vec)
xlabel('Distance from load (um)')
ylabel('Number of beads)')
title('Bead Number vs Distance')
