%Calculate the drift correction fit based on all of the displacements and
%centroids from all timepoints. The idea is that the swelling will be
%common to all time points since it's unique to the reference frame.
%Won't correct for temp swelling between timepoints

allCents=cell2mat(matches(1));
allDisps=cell2mat(displacements(1));
[numBeads,~]=size(allCents);
CENTROIDSR=[[1:numBeads]',allCents(:,4:6)];

%Use one out of every 10 points to save time
polymodelx = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,1),9);
polymodely = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,2),5);
polymodelz = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,3),3);


dispsC={};
for ii=1:2
    index=ii;
[numBeads,~]=size(matches{index}(:,1));
CENTROIDSR=[[1:numBeads]',matches{index}(:,4:6)];
dispsC{index}(:,1)=polyvaln(polymodelx,CENTROIDSR(:,2:4));
dispsC{index}(:,2)=polyvaln(polymodely,CENTROIDSR(:,2:4));
dispsC{index}(:,3)=polyvaln(polymodelz,CENTROIDSR(:,2:4));
displacementsR_dedrift{index}=[displacements{index}(:,1)-dispsC{index}(:,1), ...
    displacements{index}(:,2)-dispsC{index}(:,2),displacements{index}(:,3)-dispsC{index}(:,3)];

if ii == 2
h = figure;
filename = 'testAnimated.gif';
scale=0.75;
quiverC3D(CENTROIDSR(:,2),CENTROIDSR(:,3),CENTROIDSR(:,4),displacementsR_dedrift{index}(:,1)*scale,displacementsR_dedrift{index}(:,2)*scale,displacementsR_dedrift{index}(:,3)*scale,0)
axis equal
set(gca,'Color','k')
% This command will rotate the plot by 180 degree to match the cell image
set(gca,'xdir','reverse','ydir','reverse')
% set(gcf, 'Position',  [100, 100, 1867, 1867])
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if ii == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end
%       close gcf
end
end


