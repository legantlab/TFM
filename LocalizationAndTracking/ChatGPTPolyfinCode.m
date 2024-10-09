% Load the displacement field data
disp_field = displacements{1};
centroids = matches{1};
% Extract the x, y, and z components of the displacement field


% Define the order of the polynomial fit
poly_order_x = 5;
poly_order_y = 5;
poly_order_z = 3;

% Fit a polynomial to the x, y, and z components of the displacement field
p_x = polyfitn(centroids(:,1:3), disp_field(:,1), poly_order_x);
p_y = polyfitn(centroids(:,1:3), disp_field(:,2), poly_order_y);
p_z = polyfitn(centroids(:,1:3), disp_field(:,3), poly_order_z);

% Use the fitted polynomials to remove noise from the displacement field
disp_x_smooth = polyvaln(p_x, centroids(:,1:3));
disp_y_smooth = polyvaln(p_y, centroids(:,1:3));
disp_z_smooth = polyvaln(p_z, centroids(:,1:3));

% Remove structured noise from the displacement field
disp_x_smooth = smoothdata(disp_x_smooth, 'movmean', 5);
disp_y_smooth = smoothdata(disp_y_smooth, 'movmean', 5);
disp_z_smooth = smoothdata(disp_z_smooth, 'movmean', 5);

% Plot the original and smoothed displacement fields
figure;
subplot(1,2,1);
quiver3(centroids(:,1), centroids(:,2), centroids(:,3),disp_field(:,1), disp_field(:,2), disp_field(:,3),0);
title('Original Displacement Field');

subplot(1,2,2);
quiver3(centroids(:,1), centroids(:,2), centroids(:,3), disp_x_smooth, disp_y_smooth, disp_z_smooth,0);
title('Smoothed Displacement Field');
