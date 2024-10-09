function [matches] = polyfitThermalExpansion(matches,displacements,xOrd, yOrd, zOrd, graph)
%polyfitThermalExpansion Fits displacements to polyfits of arbitrary order
%to reduce displacements due to swelling. *Deprecated* in favor of swelling
%correction via optistruct thermal expansion simulation.
%{ 
  Uses polyfitn to fits centroids and displacement data to higher order
  polynomials specified in function to reduce displacements due to swelling
  effect on substrates. This in hydrogels is often thermal or osmolarity
  based or in PDMS substrates can be a function of surface energy with
  addition of SDS to remove cells from substrate. 
        *Deprecated* in favor of swelling
        correction via optistruct thermal expansion simulation.
Inputs:
- matches: cell array containing matched bead positions from localization
and tracking code
- displacements: cell array containing displacements from localization and
tracking code
-xOrd/yOrd,zOrd: polyfitn order for x/y/z 
- graph: Flag to output graphs of the displacements pre and post correction
Outputs:
- matches: Cell array containing centroid positions for
each time point in TFM movie now corrected with polyfitn swelling
estimation.
%}
    
 
    allCents=cell2mat(matches(1:end-1));
    allDisps=cell2mat(displacements(1:end-1));
    [numBeads,~]=size(allCents);
    CENTROIDSR=[[1:numBeads]',allCents(:,4:6)];
    
    
    %Use one out of every 10 points to save time
    polymodelx = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,1),xOrd);
    polymodely = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,2),yOrd);
    polymodelz = polyfitn(CENTROIDSR(1:1:end,2:4),allDisps(1:1:end,3),zOrd);
    
    dispsC={};
    for ii=1:length(matches)
        
        [numBeads,~]=size(matches{ii}(:,1));
        CENTROIDSR=[[1:numBeads]',matches{ii}(:,4:6)];
        dispsC{ii}(:,1)=polyvaln(polymodelx,CENTROIDSR(:,2:4));
        dispsC{ii}(:,2)=polyvaln(polymodely,CENTROIDSR(:,2:4));
        dispsC{ii}(:,3)=polyvaln(polymodelz,CENTROIDSR(:,2:4));
        displacementsR_dedrift{ii}=[displacements{ii}(:,1)-dispsC{ii}(:,1),...
            displacements{ii}(:,2)-dispsC{ii}(:,2),displacements{ii}(:,3)-dispsC{ii}(:,3)];
    end
    
    %% Plot dedrift and polymodels

    if graph == 1
        scale = 1;
        [numBeads,~]=size(matches{1}(:,1));
        CENTROIDSR=[[1:numBeads]',matches{1}(:,4:6)];
        % figure
        % subplot(1,3,1)
        % testSpace = 1:100;
        % xmodel = polyvaln(polymodelx,testSpace);
        % ymodel = polyvaln(polymodely,testSpace);
        % zmodel = polyvaln(polymodelz,testSpace);
        % plot(testSpace,xmodel)
        % title('X Model')
        % subplot(1,3,2)
        % plot(testSpace,ymodel)
        % title('Y Model')
        % subplot(1,3,3)
        % plot(testSpace,zmodel)
        % title('Z Model')
        
        figure
        subplot(1,2,1)
        quiver3(CENTROIDSR(:,2),CENTROIDSR(:,3),CENTROIDSR(:,4),...
            displacementsR_dedrift{1}(:,1)*scale,displacementsR_dedrift{1}(:,2)*scale,...
            displacementsR_dedrift{1}(:,3)*scale,0)
        title('DeDrifted')
        subplot(1,2,2)
        quiver3(matches{1}(:,4),matches{1}(:,5),matches{1}(:,6),displacements{1}(:,1),...
            displacements{1}(:,2),displacements{1}(:,3),0)
        title('Original')
    end
    %% Replace matches displacements
    %If the displacements look better, then we can set the original matrices to
    %the dedrifted ones and carry on. --MH 2023/02/14
    for j = 1:length(matches)
        curDriftCorrectedDisp = displacementsR_dedrift{j};
        matches{j}(:,1) = matches{j}(:,1) + curDriftCorrectedDisp(:,1);
        matches{j}(:,2) = matches{j}(:,2) + curDriftCorrectedDisp(:,2);
        matches{j}(:,3) = matches{j}(:,3) + curDriftCorrectedDisp(:,3);
    
    end

end