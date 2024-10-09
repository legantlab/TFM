x_min = 400;
x_max = 2300;
y_min = 500;
y_max = 600;
z_min = 160;
z_max = 180;
point_num = 10000;
add_point = zeros(point_num,3);
add_point(:,1) = x_min+(x_max-x_min)*rand(point_num,1);
add_point(:,2) = y_min+(y_max-y_min)*rand(point_num,1);
add_point(:,3) = z_min+(z_max-z_min)*rand(point_num,1);

beadpos_final = cat(1,beadpos,add_point);
idx = randperm(size(beadpos_final,1),20000);
beadpos_rand = beadpos_final(idx',:);

k_alt = boundary(beadpos_rand,1);

figure
scatter3(beadpos_rand(1:2:end,1),beadpos_rand(1:2:end,2),beadpos_rand(1:2:end,3),'.')
hold on
trisurf(k_alt,beadpos_rand(:,1),beadpos_rand(:,2),beadpos_rand(:,3),'Facecolor','red','FaceAlpha',0.1)

