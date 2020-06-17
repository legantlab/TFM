%% Define conditions used to create inverse problem
%Choose COMSOL data files to create G and u
greensMatLoadFile = 'Z:\Patel\TractionForceCode_underDevelopment\COMSOL\Force_recon_simulation\GreensMat_block_channel_3_30_20.mat';   %load the the greens mat already run through (invert mat of Comsol file)
greensMatDataFile = 'Z:\Patel\TractionForceCode_underDevelopment\COMSOL\nodaltest_append_channel.csv';    %output from data query
greensMatSaveAsFile = 'Z:\Patel\TractionForceCode_underDevelopment\COMSOL\Force_recon_simulation\GreensMat_block_channel_3_30_20.mat'; %place to save the inver greens mat
COMSOLGeometryMeshFile = 'Z:\Patel\TractionForceCode_underDevelopment\COMSOL\meshchannel_3_29_20.mphtxt'; %by default, COMSOL saves this to '../data/meshData.mphtxt'

%particleTrackingDataFile = 'Y:\Regan\2019_10_10(TFM)\data_FLAT_20kPa_Analysis\cell10\cell10_tracking.mat'; %generated from the particle tracking matlab script


%% setup parameters
%pointforce_rad_arr = [10 5 3 2 1];             % radius of point force, in um
pointforce_rad = 3;
pointforce_load_arr = [1 0.5 0.3 0.1 0.05 0.03 0.01];           % load of point force, as a ratio of gel modulus
%pointforce_load_arr = -[1 0.01];
%pointforce_direction_arr = [2 3];           % direction of point force, 1 for x, 2 for y and 3 for z
pointforce_direction_arr = [2 3]; 
%pointforce_coord = [0 0 15];        %center of applied force
pointforce_coord_arr = [0 9 16; 0 10 18; 0 11 20; 0 12 22; 0 13 24];
%pointforce_coord_arr = [0 0 15; 0 2 15; 0 4 15; 0 6 15; 0 8 15];
%pointforce_coord_arr = [0 14 25; 0 16 25; 0 18 25; 0 20 25; 0 28 25; 0 36 25; 0 44 25;0 48 25];
channel_top_z = 25;           %top of channel z value
channel_bot_z = 15;           %bottom of channel z value
N_test = 10;

today = datestr(now, 'mm-dd');
sumName = [today '_channel_side_force_recon_variation_uninoise'];

force_var_vec = zeros(length(pointforce_direction_arr),length(pointforce_load_arr),size(pointforce_coord_arr,1));
force_arb_var_vec = zeros(length(pointforce_direction_arr),length(pointforce_load_arr),size(pointforce_coord_arr,1));  
direction_index = 0;
for pointforce_direction = pointforce_direction_arr
    direction_index = direction_index+1;
    if pointforce_direction == 1
        potinforce_direction_str = 'X';
    elseif pointforce_direction == 2
        pointforce_direction_str = 'Y';
    elseif pointforce_direction == 3
        pointforce_direction_str = 'Z';
    end

    load_index = 0;
  
    for pointforce_load = pointforce_load_arr
        position_index = 0;

        load_index = load_index+1;
            %for pointforce_rad = pointforce_rad_arr
        for pointforce_coord = pointforce_coord_arr'
            %%
            %Choose whether Greens Matrix should be loaded from another file or
            %    generated from a COMSOL csv file
            %enter -1 if you already have the data loaded and will not clear vars
            %    beforehand
            position_index = position_index+1;
            clc; close all;
            today = datestr(now, 'mm-dd');
            loadFromFile = -1;             % 1 is to load, will use greensMatLoadFil; 0 if need to make greensMat; -1 if already in workspace

            %Regularization method
            method = 'tikh';

            %Choose where displacement vector is taken from
            %  -1 = u is already defined
            %   0 = load in .mat file;
            %   1 = linear combination of Green's matrix columns (see next constant)
            dispVectorSource = 1;

            %if dispVectorSource == 1, give the linear combination of Green's matrix
            %   vectors that you want to test


            %greensColumns = [1];%round(linspace(1,1944,5));
            %weights = [1];%(-1).^round(linspace(1,1944,5)).*round(linspace(1,1944,5))./1944;

            %Standard deviation is defined in this script relative to the maximum power
            %   in the Green's matrix (p). Use noiseStdFac to define the noise standard
            %   deviation as 10^(p-noiseStdFac). Define as NaN for no noise
            noiseStdFac = NaN;

            %The name that figures at the end of the script are saved as
            runName = [today '_channel_side' num2str(pointforce_coord(2)) num2str(pointforce_coord(3)) '_R' num2str(pointforce_rad) '_L' num2str(abs(pointforce_load)*100, '%03d') '_' pointforce_direction_str '_uninoise'];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %START
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Setup

            % Load Green's matrix

            if loadFromFile == 1
            disp('Loading an already computed Greens Matrix')
            load(greensMatLoadFile)
            %     load('nonSingular_block_flipped_RPM_10_16_19.mat');
            disp('done')
            elseif loadFromFile == 0
            disp('Reading COMSOL mesh data')
            % Find centroids of where forces are applied (for plotting purposes)
            [meshPointsCoords, meshPointsIndex, meshSurfaceIndex, desiredBoundaryInd] = parseMeshData(COMSOLGeometryMeshFile);

            %the mesh centroids will be used to visually approximate where forces were
            %   applied
            meshCentroid = [];
            for i = 1:size(meshSurfaceIndex)
                %look for a meshElement that lives on the surface you applied forces to
                if find(desiredBoundaryInd == meshSurfaceIndex(i))% & sum(i==[1,6,12,17,23,28,34,39,45,50,56,61,67,72,78,83,89,94,100,105,111,116,122,127,133,138,144,149,155,160])
                     currMesh = meshPointsIndex(i,:)+1;
                    currCoord1 = meshPointsCoords(currMesh(1),:);
                    currCoord2 = meshPointsCoords(currMesh(2),:);
                    currCoord3 = meshPointsCoords(currMesh(3),:);
                    meshCentroid = [meshCentroid; (currCoord1+currCoord2+currCoord3)/3];
                end
            end
            disp('Reading COMSOL displacement data')
            [G, displacedPtsMat] = GreensMatParamSweep(greensMatDataFile);    
            disp(['Greens Matrix compiled. It is a ' num2str(size(G,1)) ' x ' num2str(size(G,2)) ' matrix.'])
            %%
            disp('Computing SVD')
            G = single(G);
            [U,s,V] = csvd(G);
            disp(['condition number = ' num2str(abs(max(s)/min(s)))]);
            nonSingular = s>(max(size(G))*eps(norm(G)));
            G = double(G);
            U = double(U);
            V = double(V);
            s = double(s);

            disp('saving Greens Matrix')
            save(greensMatSaveAsFile, 'G', 'U', 's', 'V', 'meshCentroid', 'displacedPtsMat','-v7.3')

            end


            %% Set up displacement vector u

%                 if dispVectorSource ~= -1
%                 % Load displacement vector
%                 if dispVectorSource == 0
%                     disp('loading vector from particle tracking file')
%                     load(BeadDispSaveAsFile);
% 
% 
%                 elseif dispVectorSource == 1
                    greensColumns=zeros(size(G,2),1);
                    disp('generating displacements vector from linear combo of Greens Matrix columns')
                    if abs(pointforce_coord(3)-channel_top_z) < 0.1   % force applied on top of the channel
                        mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2)<pointforce_rad) & abs(meshCentroid(:,3)-channel_top_z)<0.1);
                        %mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2)<pointforce_rad) & abs(meshCentroid(:,3)-channel_top_z)<0.1 & meshCentroid(:,2)<=pointforce_coord(2));        
                    elseif abs(pointforce_coord(3)-channel_bot_z)<0.1   % force applied on bottom of the channel
                        mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2)<pointforce_rad) & abs(meshCentroid(:,3)-channel_bot_z)<0.1);
                        %mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2)<pointforce_rad) & abs(meshCentroid(:,3)-channel_bot_z)<0.1 & meshCentroid(:,2)<=pointforce_coord(2));            
                    else     % force applied on side wall
                        mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2+(meshCentroid(:,3)-pointforce_coord(3)).^2)<pointforce_rad)&(meshCentroid(:,3)>channel_bot_z)&(meshCentroid(:,3)<channel_top_z));
                        %mesh_index = double((sqrt((meshCentroid(:,1)-pointforce_coord(1)).^2+(meshCentroid(:,2)-pointforce_coord(2)).^2+(meshCentroid(:,3)-pointforce_coord(3)).^2)<pointforce_rad)&(meshCentroid(:,3)>channel_bot_z)&(meshCentroid(:,3)<=pointforce_coord(3)));
                    end
                %        mesh_index = double(sqrt(meshCentroid(:,1).^2+meshCentroid(:,2).^2)<pointforce_rad);
                    greensColumns(pointforce_direction:3:end) = pointforce_load*mesh_index;
                    u = G*greensColumns;
                %       u = sum(repmat(weights,size(G,1),1).*G(:,greensColumns),2);
                %       uAllTimes = {u};
                
                force_var_mat = zeros(3*sum(mesh_index),N_test);
                force_arb_var_mat = zeros(3*sum(mesh_index),N_test);
                mesh_index_select = false(3*length(mesh_index),1);
                mesh_index_select(1:3:end) = logical(mesh_index~=0);
                mesh_index_select(2:3:end) = logical(mesh_index~=0);
                mesh_index_select(3:3:end) = logical(mesh_index~=0);
                
                for i = 1:N_test
                    close all;
                    % from flat gel displacement measurements, sigma x and y =0.0166 um, sigma z=0.09 um 
                    u(1:3:end) = u(1:3:end)+0.0166*randn(size(G,1)/3,1);
                    u(2:3:end) = u(2:3:end)+0.0166*randn(size(G,1)/3,1);
                    u(3:3:end) = u(3:3:end)+0.0166*randn(size(G,1)/3,1);       
                    BeadDispSaveAsFile = ['../COMSOL/Force_recon_simulation_multiload/BeadDisp_' runName '.mat'];
                    save(BeadDispSaveAsFile, 'u');
%                 end
%                 end



                %% Start Regularization for Each Time Point

                figure
                [reg_corner,rho,eta,reg_param] = l_curve(U,s,u,'tikh');
                %     saveas(gcf, ['../plots/' runName 'Lcurve_temp.png'])


                %% Warning/Error Throwing
                if size(G,1) < size(G,2)
                warning('system given by G is underdefined')
                end
                %     if sum(nonSingular) ~= size(G,2)
                %         warning('Greens matrix is not full rank')
                %     end

                if size(meshCentroid,1) ~= size(G,2)/3
                error('The given mesh does not match the given Greens matrix data')
                elseif size(u,1) ~= size(G,1)
                error('displacement vector dimension does not agree with Greens matrix')
                end


                %% Calculate Regularization Parameter
                %these functions and its subroutines are from "regtools" package by PC Hansen
                disp('Calculating Regularization Parameter')

                %Pseudo-inverse solution (i.e. solution without regularization - for
                %  comparison purposes)
                fMinPseudo = tikhonov(U,s,V,u,0,'tikh');

                %calculate l-curve and find regularization parameter
                %     figure
                %[reg_corner,rho,eta,reg_param] = l_curve(U,s,u,'tikh',fMinPseudo);
                [reg_corner,rho,eta,reg_param] = l_curve_ysedit(U,s,u,'tikh',fMinPseudo);
                saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' runName 'Lcurve_temp.png'])

                %% Find force profile that minimizes error
                disp('Error Optimization')

                %calculate solution based on regularization parameter found above

                if reg_corner < 1e-4 || reg_corner > 10
                    reg_corner = 5;
                end

                fMinTikh = tikhonov(U,s,V,u,reg_corner,'tikh');
                tractions{i}=fMinTikh;
                absErrorTikh = norm(fMinTikh-fMinPseudo);
                disp(['|Tikh - psuedoInv| = ' num2str(absErrorTikh)]);

                %     %calculate solution based on an arbitrary regularization parameter (for
                %     %  comparison purposes)
                regArb = 3;
                fMinTikhArb = tikhonov(U,s,V,u,regArb,'tikh');
                absErrorTikhArb = norm(fMinTikhArb-fMinPseudo);
                disp(['|Tikh - psuedoInv| = ' num2str(absErrorTikhArb)]);
                force_var_mat(:,i) = fMinTikh(mesh_index_select);
                force_arb_var_mat(:,i) = fMinTikhArb(mesh_index_select);
            end
            
            force_var = mean(std(force_var_mat,0,2))/abs(pointforce_load);
            force_arb_var = mean(std(force_arb_var_mat,0,2))/abs(pointforce_load);
            
            force_var_vec(direction_index,load_index,position_index)=force_var;
            force_arb_var_vec(direction_index,load_index,position_index)=force_arb_var;
%             %write solution to a file readable by Amira
%             text = buildAmiraMesh(fMinTikh, meshCentroid);
%             fid = fopen(['../COMSOL/Force_recon_simulation/AmiraTetraVectors_10_22_19_num_cell10' num2str(i) '.am'],'wt');
%             fprintf(fid, text);
%             fclose(fid);

            save(['../COMSOL/Force_recon_simulation_multiload/Varvec_' sumName '.mat'], 'force_var_vec');
            save(['../COMSOL/Force_recon_simulation_multiload/Vararbvec_' sumName '.mat'], 'force_arb_var_vec');
            
            %% Generate line profile (x=0, z = 100, along y direction)
            if pointforce_coord(3) < channel_top_z && pointforce_coord(3) > channel_bot_z
            mesh_index_lineprofile = abs(meshCentroid(:,2)-pointforce_coord(2))<3 & abs(meshCentroid(:,3)-pointforce_coord(3))<1;
            else
            mesh_index_lineprofile = abs(meshCentroid(:,2)-pointforce_coord(2))<3 & abs(meshCentroid(:,3)-pointforce_coord(3))<0.1;
            end
            meshCentroid_lineprof = meshCentroid(mesh_index_lineprofile,1);
            force_lineprof_x = abs(fMinTikh(1:3:end));
            force_lineprof_x = 100*force_lineprof_x(mesh_index_lineprofile);
            force_lineprof_y = abs(fMinTikh(2:3:end));
            force_lineprof_y = 100*force_lineprof_y(mesh_index_lineprofile);
            force_lineprof_z = abs(fMinTikh(3:3:end));
            force_lineprof_z = 100*force_lineprof_z(mesh_index_lineprofile);
            
            force_lineprof_arb_x = abs(fMinTikhArb(1:3:end));
            force_lineprof_arb_x = 100*force_lineprof_arb_x(mesh_index_lineprofile);
            force_lineprof_arb_y = abs(fMinTikhArb(2:3:end));
            force_lineprof_arb_y = 100*force_lineprof_arb_y(mesh_index_lineprofile);
            force_lineprof_arb_z = abs(fMinTikhArb(3:3:end));
            force_lineprof_arb_z = 100*force_lineprof_arb_z(mesh_index_lineprofile);            

            force_lineprof_theory = 100*abs(pointforce_load)*double(abs(meshCentroid_lineprof)<pointforce_rad);

            %% Plot
            %%{
            disp('plotting')
            figure
            d = 1;
            quiver3(displacedPtsMat(1:3:end,1),displacedPtsMat(1:3:end,2),displacedPtsMat(1:3:end,3),u(1:3:end),u(2:3:end),u(3:3:end),0)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title({'Displacements Used to Solve Inverse Problem'})
            %view([90 0]) %use this if you want a view of the yz plane
            % axis([domain_x domain_y domain_z])
            saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' runName 'displacements_temp.png'])


            figure
            subplot(1,2,1)
            hold on
            scatter(meshCentroid_lineprof,force_lineprof_theory,10,'b','filled')
            scatter(meshCentroid_lineprof,force_lineprof_z,10,'r','filled')
            scatter(meshCentroid_lineprof,force_lineprof_y,10,'c','filled')
            scatter(meshCentroid_lineprof,force_lineprof_x,10,'g','filled')
            xlabel('y')
            ylabel('%Emod')
            title({['Force line profile along y Nmesh=' num2str(abs(sum(mesh_index)))]})
            legend('sT','rTz','rTy','rTx')
            legend('boxoff')

            subplot(1,2,2)
            hold on
            scatter(meshCentroid_lineprof,force_lineprof_theory,10,'b','filled')
            scatter(meshCentroid_lineprof,force_lineprof_arb_z,10,'r','filled')
            scatter(meshCentroid_lineprof,force_lineprof_arb_y,10,'c','filled')
            scatter(meshCentroid_lineprof,force_lineprof_arb_x,10,'g','filled')
            xlabel('y')
            ylabel('%Emod')
            title({'Force line profile along y,fixed lambda'})
            legend('sT','rTz','rTy','rTx')
            legend('boxoff')            
            
            saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' runName 'forcelineprof.png'])

%             scale=0.5;
%             figure
%             subplot(1,3,1)
%             quiver3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),scale*d*greensColumns(1:3:end),scale*d*greensColumns(2:3:end),scale*d*greensColumns(3:3:end),0)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             title('Force Profile Applied')
%             %view([90 0]) %use this if you want a view of the yz plane
%             % axis([domain_x domain_y [0 10]])
% 
%             subplot(1,3,2)
%             quiver3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),scale*d*fMinTikh(1:3:end),scale*d*fMinTikh(2:3:end),scale*d*fMinTikh(3:3:end),0)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             title({['Force Profile Computed via L2 Regularization; |f_{min}-f_{pinv}| = ' num2str(absErrorTikh)]; ['Tikhonov regularization: \lambda = ' num2str(reg_corner)]})
%             %view([90 0]) %use this if you want a view of the yz plane
%             % axis([domain_x domain_y [0 10]])
% 
%             subplot(1,3,3)
%             quiver3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),scale*d*fMinTikhArb(1:3:end),scale*d*fMinTikhArb(2:3:end),scale*d*fMinTikhArb(3:3:end),0)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             title({['Force Profile Computed via L2 Regularization; |f_{min}-f_{pinv}| = ' num2str(absErrorTikhArb)]; ['Tikhonov regularization (arbitrary value): \lambda = ' num2str(regArb)]})
%             %view([90 0]) %use this if you want a view of the yz plane
%             % axis([domain_x domain_y [0 10]])
% 
%             saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' runName 'inverseSolution_temp.png'])


            figure
            subplot(2,2,1)
            scatter3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),10,abs(greensColumns(pointforce_direction:3:end)),'filled')
            colormap(jet)
            caxis([0 max(abs(greensColumns(pointforce_direction:3:end)))]);
            colorbar
            xlabel('x')
            ylabel('y')
            title('Force Profile Applied')
            %view([90 0]) %use this if you want a view of the yz plane
            % axis([domain_x domain_y [0 10]])

            subplot(2,2,2)
            scatter3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),10,abs(fMinTikh(pointforce_direction:3:end)),'filled')
            colormap(jet)
            caxis([0 max(abs(greensColumns(pointforce_direction:3:end)))]);
            colorbar
            xlabel('x')
            ylabel('y')
            title({['Force Profile Computed via L2 Regularization; |f_{min}-f_{pinv}| = ' num2str(absErrorTikh)]; ['Tikhonov regularization: \lambda = ' num2str(reg_corner)]})
            %view([90 0]) %use this if you want a view of the yz plane
            % axis([domain_x domain_y [0 10]])

            subplot(2,2,3)
            scatter3(meshCentroid(:,1),meshCentroid(:,2),meshCentroid(:,3),10,abs(fMinTikhArb(pointforce_direction:3:end)),'filled')
            colormap(jet)
            caxis([0 max(abs(greensColumns(pointforce_direction:3:end)))]);
            colorbar
            xlabel('x')
            ylabel('y')
            title({['Force Profile Computed via L2 Regularization; |f_{min}-f_{pinv}| = ' num2str(absErrorTikhArb)]; ['Tikhonov regularization (arbitrary value): \lambda = ' num2str(regArb)]})
            %view([90 0]) %use this if you want a view of the yz plane
            % axis([domain_x domain_y [0 10]])

            saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' runName 'Solutinoheatmap_temp.fig'])
            %}
        end
    end
end


figure
hold on
for i = 1:size(pointforce_coord_arr,1)
    plot(abs(pointforce_load_arr)*100,force_var_vec(1,:,i))
end
legend(num2str(pointforce_coord_arr));
xlabel('force load')
ylabel('scaled variation')
title('force variation in Y')
set(gca,'Xscale','log')
saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' sumName '_Y.png'])


figure
hold on
for i = 1:size(pointforce_coord_arr,1)
    plot(abs(pointforce_load_arr)*100,force_var_vec(2,:,i))
end
legend(num2str(pointforce_coord_arr));
xlabel('force load')
ylabel('scaled variation')    
title('force variation in Z')
set(gca,'Xscale','log')
saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' sumName '_Z.png'])

figure
hold on
for i = 1:size(pointforce_coord_arr,1)
    plot(abs(pointforce_load_arr)*100,force_arb_var_vec(1,:,i))
end
legend(num2str(pointforce_coord_arr));
xlabel('force load')
ylabel('scaled variation')
title('force variation in Y,fixed lambda')
set(gca,'Xscale','log')
saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' sumName '_fixlambda_Y.png'])


figure
hold on
for i = 1:size(pointforce_coord_arr,1)
    plot(abs(pointforce_load_arr)*100,force_arb_var_vec(2,:,i))
end
legend(num2str(pointforce_coord_arr));
xlabel('force load')
ylabel('scaled variation')    
title('force variation in Z,fixed lambda')
set(gca,'Xscale','log')
saveas(gcf, ['../COMSOL/Force_recon_simulation_multiload/' sumName '_fixlambda_Z.png'])