%Troubleshooting function to plot forces at a few ranges of regularization
%values to assess results. 
reg_values = [10,5,2,1,0.5,0.1,0.01,0.001];
force_vector={};
for kk=1:length(reg_values)
    kk
    temp=beadDisps{1}(:,4:6)';
    temp(:,3) = temp(:,3);
    y=temp(:);
    reg_corner_timelapse=reg_values(kk); %L-curve picks out the wrong corner right now
    [x_lambda,rho,eta] = tikhonov2(U,s,V,y,reg_corner_timelapse);
    %Removing last row of x_Lamda corresponds to swelling I think. 
    force_vector{kk}=[elemCents2D,x_lambda(1:3:end-1),...
        x_lambda(2:3:end-1),x_lambda(3:3:end-1)];
    figure
    quiver3(elemCents2D(:,1),elemCents2D(:,2),elemCents2D(:,3),force_vector{kk}(:,4),...
        force_vector{kk}(:,5),force_vector{kk}(:,6),1)
    title(['Regularization lamda = ', num2str(reg_values(kk))])
end



