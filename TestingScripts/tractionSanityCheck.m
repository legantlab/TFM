%Sanity check script for regularization and masking of traction maps

testDispField = correctedbeadDisps{10}(:,4:6);

temp=testDispField';
%temp(:,~coordinateMask) = 0; 
y=temp(:);

%Compute L-curve
fMinPseudo = tikhonov(U2,s2,V2,y,0,'tikh');
[reg_corner,rho,eta,reg_param] = l_curve(U2,s2,y,'tikh');

%Loop through some test regularization values and save resulting traction
%field. 

force_vector={};
reg_corner_compare = [0.01, 0.1, 1, 2, 3, 4, 5,6,7,8,9,10,15,20,25,30];%reg_corner*2;

for kk=1:length(reg_corner_compare)
    kk
    %temp=correctedbeadDisps{10}(:,4:6)';
    temp=testDispField';
    %temp(:,~coordinateMask) = 0;
    y=temp(:);
    reg_corner_timelapse=reg_corner_compare(kk); %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U2,s2,V2,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[elemCents2D,x_lambda(1:3:end),...
        x_lambda(2:3:end),x_lambda(3:3:end)];

    %Plot traction field
%     
%     scale = 1;
%     figure
%     %subplot(4,2,kk)
%     quiver(elemCents2D(:,1),elemCents2D(:,2),force_vector{1}(:,4)*scale,force_vector{1}(:,5)*scale,1)
%     %hold on 
%     axis equal
%     title(['Reg = ', num2str(reg_corner_timelapse)])

%     figure
%     quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{1}(:,4)*scale,force_vector{1}(:,5)*scale,force_vector{1}(:,6)*scale,0)
%     axis equal
end

%% Now compute residual error at each element. 
close all
meanResiduals = zeros(length(force_vector),1);
for j = 1:length(force_vector)
    %close all
    testForces = force_vector{j}*1;
    
%     figure
%     quiver3(testForces(:,1),testForces(:,2),testForces(:,3),testForces(:,4),testForces(:,5),...
%         testForces(:,6),1,'r')
%     hold on 
%     quiver3(testDispField(:,1),testDispField(:,2),testDispField(:,3),...
%         testDispField(:,4),testDispField(:,5),...
%         testDispField(:,6),1,'b')
%     title('Interpolated Displacements and Computed Tractions')
%     legend('Computed Tractions','Displacements')
%     axis equal
    
    %% Compute predicted displacements
    x_lam = zeros(length(testForces(:,1))*3,1);
    x_lam(1:3:end) = testForces(:,4);
    x_lam(2:3:end) = testForces(:,5);
    x_lam(3:3:end) = testForces(:,6);
    
    pred_disps = beadu*x_lam;
    
    
    %Plot
%     figure
%     quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),pred_disps(1:3:end), ...
%         pred_disps(2:3:end), pred_disps(3:3:end),0,'r')
%     hold on 
%     quiver3(testDispField(:,1),testDispField(:,2),testDispField(:,3),...
%         testDispField(:,4),testDispField(:,5),...
%         testDispField(:,6),0,'b')
%     %quiver3(elemCents2D(:,1), elemCents2D(:,2),elemCents2D(:,3),beadDisps{1}(:,4), ...
%     %    beadDisps{1}(:,5), beadDisps{1}(:,6),0,'b')
%     legend('Predicted Displacments','Interpolated Displacements')
%     axis equal tight
    
    %Plot residuals
    residuals = [testDispField(:,1) - pred_disps(1:3:end), testDispField(:,2) - pred_disps(2:3:end), ...
        testDispField(:,3) - pred_disps(3:3:end)];
    numBins = 20;
%     figure
%     subplot(1,3,1)
%     hist(vecnorm(residuals,1,2),numBins)
%     title('Residual error between computed and observed displacements')
%     allDisps = testDispField(:,1:3);
%     ylim([0,1800])
%     xlabel('RE (microns)')
%     ylabel('Count')
%     
%     subplot(1,3,2)
%     hist(vecnorm(allDisps,1,2),numBins)
%     title('Histogram of Raw Displacements')
%     xlabel('RE (microns)')
%     ylabel('Count')
%     ylim([0,1800])
%     
%     subplot(1,3,3)
%     corallDisps = testDispField(:,1:3);
%     hist(vecnorm(corallDisps,1,2),numBins)
%     title('Histogram of Corrected Displacements')
%     ylim([0,1800])
%     xlabel('RE (microns)')
%     ylabel('Count')
%    figure
%    histogram(vecnorm(residuals,1,2),25)
%    hold on 
%    histogram(vecnorm(testDispField,1,2),25)
%    ylabel('Counts')
%    xlabel('Displacement (um)')
%    title(['Histograms of interoplated displacements and residuals Reg = ',num2str(reg_corner_compare(j))])




%     %Plot 3Dresidual errors onto elements figure
%     quiver3(elemCents2D(:,1),
%     elemCents2D(:,2),elemCents2D(:,3),residuals(:,1), ...
%         residuals(:,2), residuals(:,3),0,'r')
%     hold on quiver3(elemCents2D(:,1),
%     elemCents2D(:,2),elemCents2D(:,3),beadDisps{1}(:,4), ...
%          beadDisps{1}(:,5), beadDisps{1}(:,6),0,'b')
%     quiver3(elemCents2D(:,1),
%     elemCents2D(:,2),elemCents2D(:,3),pred_disps(1:3:end), ...
%          pred_disps(2:3:end), pred_disps(3:3:end),0,'g')
%     title(['Residual error between computed and observed displacements ',
%     'RegValue = ', ...
%         num2str(reg_corner_compare(j))])
%     legend('Residual Error','Computed Displacements','Predicted
%     Displacements') axis equal

    %2D quiver Plots
    figure
    quiver(elemCents2D(:,1), elemCents2D(:,2),residuals(:,1), ...
        residuals(:,2),0,'r')
    hold on 
    quiver(elemCents2D(:,1), elemCents2D(:,2),testDispField(:,1), ...
         testDispField(:,2),0,'b')
    quiver(elemCents2D(:,1), elemCents2D(:,2),pred_disps(1:3:end), ...
         pred_disps(2:3:end),0,'g')
    quiver(testForces(:,1), testForces(:,2), testForces(:,4), testForces(:,5),1)
    title(['Residual error between computed and observed displacements ', 'RegValue = ', ...
        num2str(reg_corner_compare(j))])
    legend('Residual Error','Computed Displacements','Predicted Displacements','Traction Forces')
    axis equal
    %compute distribution fits of residuals
    pd = fitdist(vecnorm(residuals,1,2),'Normal');

    disp('Mean Normalized Residual Error')
    disp(mean(vecnorm(residuals,1,2)))
    meanResiduals(j) = mean(vecnorm(residuals,1,2));
end
testfit = fit(reg_corner_compare',meanResiduals,'spline');

figure
plot(reg_corner_compare,meanResiduals')
hold on 
%plot(logspace(-1,1.5,1000),testfit(logspace(-1,1.5,1000)))
xlabel('Regularization parameter')
ylabel('Mean Residual')
yline(0.3)
yline(0.4)
