pointforce_rad = 3;
pointforce_load = 0.03;
pointforce_direction_arr = [2 3];      % direction of point force, 1 for x, 2 for y and 3 for z
pointforce_coord_arr = [0 9 16; 0 10 18; 0 11 20; 0 12 22; 0 13 24];
%pointforce_coord_arr = [0 0 15; 0 2 15; 0 4 15; 0 6 15; 0 8 15];
%pointforce_coord_arr = [0 16 25; 0 18 25; 0 20 25; 0 28 25; 0 36 25; 0 48 25];

sum_name = '_channel_side';

Bead_coord = displacedPtsMat(1:3:end,1:3);
n_bead = length(Bead_coord);
max_dist = 30;
min_dist = 3;
nbin = 10;
interval = (max_dist-min_dist)/(nbin-1);

today = datestr(now, 'mm-dd');

for pointforce_direction = pointforce_direction_arr

    close all;
    if pointforce_direction == 1
        potinforce_direction_str = 'X';
    elseif pointforce_direction == 2
        pointforce_direction_str = 'Y';
    elseif pointforce_direction == 3
        pointforce_direction_str = 'Z';
    end
    
    f1 = figure;
    f2 = figure;
        
    for pointforce_coord = pointforce_coord_arr'      

        dist_vec = linspace(min_dist,max_dist,nbin);
        dist_count = zeros(nbin,1);
        dist_disp = zeros(nbin,1);        
        
        runName = ['04-20' sum_name num2str(pointforce_coord(2)) num2str(pointforce_coord(3)) '_R' num2str(pointforce_rad) '_L' num2str(abs(pointforce_load)*100, '%03d') '_' pointforce_direction_str '_uninoise'];
        BeadDispSaveAsFile = ['../COMSOL/Force_recon_simulation_multiload/BeadDisp_' runName '.mat'];
        load(BeadDispSaveAsFile)
        
        u_x=u(1:3:end);
        u_y=u(2:3:end);
        u_z=u(3:3:end);
        
        for i = 1:n_bead
            coord = Bead_coord(i,1:3);
            dist = sqrt((pointforce_coord(1)-coord(1))^2+(pointforce_coord(2)-coord(2))^2+(pointforce_coord(3)-coord(3))^2);
            if dist < max_dist
                index = 1+max(ceil((dist - min_dist)/interval),0);
                dist_count(index) = dist_count(index)+1;
                dist_disp(index) = dist_disp(index)+sqrt(u_x(i)^2+u_y(i)^2+u_z(i)^2);
            end
        end
        dist_disp = dist_disp./dist_count;
        dist_quapara = dist_disp.*sqrt(dist_count);
        
        figure(f1);
        hold on
        plot(dist_vec,dist_disp)
        
        figure(f2);
        hold on
        plot(dist_vec,dist_quapara)
        
    end

    figure(f1);
    legend(num2str(pointforce_coord_arr));
    xlabel('Distance from load (um)')
    ylabel('Averaged Displacement (um)')
    title('Displacement vs Distance')
    saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' today sum_name '_DispvsDist' pointforce_direction_str '.png'])
    
    
    figure(f2);
    legend(num2str(pointforce_coord_arr));
    xlabel('Distance from load (um)')
    ylabel('Number of beads)')
    title('Recon Quality vs Distance')
    saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' today sum_name '_QuaparavsDist' pointforce_direction_str '.png'])
    
end

