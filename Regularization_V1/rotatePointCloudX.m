function [Final_Spots_rotated]=rotatePointCloudX(cloud,theta)
R=[1,0,0;0,cosd(theta),-sind(theta);0,sind(theta),cosd(theta)];
Final_Spots_rotated=[R*cloud']';
end