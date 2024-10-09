cutoff_edge = 550;
%beadpos_sort_affine = sortrows(beadpos,2);

beadpos_left = beadpos_sort_affine(beadpos_sort_affine(:,2)<cutoff_edge,:);
beadpos_right = beadpos_sort_affine(beadpos_sort_affine(:,2)>cutoff_edge,:);

k_left = boundary(beadpos_left);
k_right = boundary(beadpos_right);
k_right = k_right+size(beadpos_left,1);
k_combine = cat(1,k_left,k_right);

figure
scatter3(beadpos(1:2:end,1),beadpos(1:2:end,2),beadpos(1:2:end,3),'.')
hold on

%trisurf(k_left,beadpos_left(:,1),beadpos_left(:,2),beadpos_left(:,3),'Facecolor','red','FaceAlpha',0.1)
%trisurf(k_right,beadpos_right(:,1),beadpos_right(:,2),beadpos_right(:,3),'Facecolor','red','FaceAlpha',0.1)

trisurf(k_combine,beadpos_sort_affine(:,1),beadpos_sort_affine(:,2),beadpos_sort_affine(:,3),'Facecolor','red','FaceAlpha',0.1)
