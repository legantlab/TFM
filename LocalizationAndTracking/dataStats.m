function [] = dataStats(cellMask,dispPoints, tracPoints, trans, tractions, displacements, time)
%dataStats Computes distances to boundary statistics 
%   Takes input TFM data and cellular mask at a given time point and
%   outputs histograms and scatter plots of magnitude data and magnitude
%   versus distance from the cell boundary. Saves the resulting plots in
%   the plots folder (subdirectory). 

            masks = cellMask;
%Retrieve list of points and their relative
            %displacements/tractions magnitudes. For each time point, the
            %points remain the same, but the magnitudes of the forces
            %change. Therefore we can generalize to the mesh/pts coords. We
            %also need to translate the coordinate axis to line up the
            %image and the mesh as before. 
            
            %dispPoints = dispPoints;%dispPtsCondensed;
            dispPoints(:,1) = dispPoints(:,1)+ trans(1);
            dispPoints(:,2) = dispPoints(:,2)+ trans(2);
            tracPoints = tracPoints(:,1:2);%meshCentroid(:,1:2); Only take xy dimension
            tracPoints(:,1) = tracPoints(:,1) + trans(1);
            tracPoints(:,2) = tracPoints(:,2) + trans(2);
            %Create a cell array for each time, containing the distances at
            %each row for a given time point.
            distances = cell(1,2);
            
            for j = 1:numel(distances(:,1))
                %Populate with for loops
                curMask = masks;
                distDisp = zeros(numel(dispPoints(:,1)),1);
                for k = 1:numel(dispPoints(:,1))
                   %Compute distance between mesh points and boundarys
                   temp = inf;
                   curPoint = dispPoints(k,1:3);
                   for l = 1:numel(curMask(:,1))
                       x = curMask(l,1) - curPoint(1);
                       y = curMask(l,2) - curPoint(2);
                       %z = 100 - curPoint(3);
                       %r = sqrt(x^2 + y^2 + z^2);
                       r = sqrt(x^2 + y^2);
                       if r < temp
                           temp = r;
                       end
                   end
                                      
                   distDisp(k,1) = temp;
                end
                
                distTrac = zeros(numel(tracPoints(:,1)),1);
                for m = 1:numel(tracPoints(:,1))
                    temp = inf;
                    curPoint = tracPoints(m,1:2);
                    for n = 1:numel(curMask(:,1))
                       x = curMask(n,1) - curPoint(1);
                       y = curMask(n,2) - curPoint(2);
                       r = sqrt(x^2 + y^2);
                       if r < temp
                           temp = r;
                       end
                    end
                    
                    distTrac(m,1) = temp;
                end
                
                distances{j,1} = distDisp;
                distances{j,2} = distTrac;
            
            end

            %Compute traction magnitudes at each point at each time
            magnitudes = cell(numel(masks),2);
            %Loop through displacments/tractions and compute magnitudes
            for o = 1:numel(magnitudes(:,1))
               curTractions = tractions;
               curDisplacements = displacements;
               %Get magnitudes from every three values
               xT = curTractions(1:3:end);
               yT = curTractions(2:3:end);
               zT = curTractions(3:3:end);
               tempTract = sqrt(xT.^2 + yT.^2 + zT.^2);
               
               xD = curDisplacements(1:3:end);
               yD = curDisplacements(2:3:end);
               zD = curDisplacements(3:3:end);
               tempDisp = sqrt(xD.^2 + yD.^2 + zD.^2);
               
               magnitudes{o,2} = tempTract;
               magnitudes{o,1} = tempDisp;
            end
           
            %Setup arrays for plots
            distDisps = [distances{1,1}, magnitudes{1,1}];
            distTracts = [distances{1,2}, magnitudes{1,2}];
            
            %Histograms of data
            figure
            histogram(magnitudes{1,1})
            title('Displacement Histogram')
            xlabel('Displacements')
            ylabel('Counts')
            saveas(gcf, ['./plots/', 'displacements_hist_', 't', ...
                num2str(time), '.png'])
            
            figure
            histogram(magnitudes{1,2})
            title('Traction Histogram')
            xlabel('Traction Forces')
            ylabel('Counts')
            saveas(gcf, ['./plots/', 'tractions_hist_', 't', ...
                num2str(time), '.png'])
            
            %Distance Plots
            figure
            scatter(distDisps(:,1), distDisps(:,2))
            title('Bead Displacments versus distance from Cell')
            xlabel('Distance from Cell Boundary')
            ylabel('Bead Displacement')
            saveas(gcf, ['./plots/', 'Displacement_dist_', 't', ...
                num2str(time), '.png'])
            
            figure
            scatter(distTracts(:,1), distTracts(:,2))
            title('Recoverd Traction Magnitude versus distance from Cell')
            xlabel('Distance from Cell Boundary')
            ylabel('Traction Force')
            saveas(gcf, ['./plots/', 'tractions_dist_', 't', ...
                num2str(time), '.png'])
end

