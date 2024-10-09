%script to compare displacement fields produced by different reference
%states. 
close all
%Start by comparing positions of beads in time 2 and 4

figure
scatter3(matches{4}(:,1), matches{4}(:,2),matches{4}(:,3))
hold on 
scatter3(matches{4}(:,4), matches{4}(:,5),matches{4}(:,6))
legend('DMSO localizations','Ideal Localizations')

%Compare quiver plots of DMSO reference with cell near channel (3)
scale = 2;
figure
quiver3(matches{1}(:,1), matches{1}(:,2),matches{1}(:,3),displacements{1}(:,1)*scale, ...
   displacements{1}(:,2)*scale,displacements{1}(:,3)*scale,0,'b')
hold on 
quiver3(matches{3}(:,1), matches{3}(:,2),matches{3}(:,3),displacements{3}(:,1)*scale, ...
   displacements{3}(:,2)*scale,displacements{3}(:,3)*scale,0,'r')
legend('DMSO reference', 'Cell')


%Compare histograms
figure
histogram(vecnorm(displacements{3}-mean(displacements{3}),2,2))
hold on 
histogram(vecnorm(displacements{2}-mean(displacements{2}),2,2))
legend('Cell Displacements', 'DMSO displacements')