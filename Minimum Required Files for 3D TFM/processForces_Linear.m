function [] = processForces_Linear(filestem,dataNameFEM,dataNameEXP,distFromCellBead,distFromCellSVD,multiple_time_pts,varargin)

if multiple_time_pts==1
    time_pt=varargin{1};
    num_time_pts=varargin{2};
    disp(['Will process ', num2str(num_time_pts-time_pt+1),' time points'])
else
    time_pt=1;
    num_time_pts=1;
    disp('Processing single time point')
end
for timept=time_pt:num_time_pts
 disp('Loading required datasets')
 eval(['load ''',filestem,'nodeList_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,'elemList_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,'voidElemConn_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,'repElems_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,'repNodes_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,'elemCheck_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),'.txt''']);
 eval(['load ''',filestem,dataNameEXP,num2str(timept,'%02.0f'),'.mat''']);
 eval(['load ''',filestem,dataNameFEM,num2str(timept,'%02.0f'),'.txt''']);
 eval([dataNameEXP,'=',dataNameEXP,num2str(timept,'%02.0f'),'(',dataNameEXP,num2str(timept,'%02.0f'),'(:,8)<=',num2str(distFromCellBead),',:);']);
 eval([dataNameEXP,'=',dataNameEXP,'(elemCheck_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),'==1,:);']);
 eval(['[beadDisplacements_',num2str(timept,'%02.0f'),']=[[1:1:length(',dataNameEXP,')]'',',dataNameEXP,'(:,1),',dataNameEXP,'(:,2),',dataNameEXP,'(:,3),[1:1:length(',dataNameEXP,')]'',',dataNameEXP,'(:,1)+',dataNameEXP,'(:,4),',dataNameEXP,'(:,2)+',dataNameEXP,'(:,5),',dataNameEXP,'(:,3)+',dataNameEXP,'(:,6)];']);
 eval(['[elemAreas_',num2str(timept,'%02.0f'),']=calcElemAreas(nodeList_m2r_',num2str(timept,'%02.0f'),',voidElemConn_m2r_',num2str(timept,'%02.0f'),');']);
 eval(['[elemNormals_',num2str(timept,'%02.0f'),']=calcElemNormals(nodeList_m2r_',num2str(timept,'%02.0f'),',voidElemConn_m2r_',num2str(timept,'%02.0f'),');']);
 eval(['[elemCents_',num2str(timept,'%02.0f'),']=elemCentFind(nodeList_m2r_',num2str(timept,'%02.0f'),',voidElemConn_m2r_',num2str(timept,'%02.0f'),');']);
 eval(['centMass_',num2str(timept,'%02.0f'),'=[sum(elemCents_',num2str(timept,'%02.0f'),'(:,1).*elemAreas_',num2str(timept,'%02.0f'),')/sum(elemAreas_',num2str(timept,'%02.0f'),'),sum(elemCents_',num2str(timept,'%02.0f'),'(:,2).*elemAreas_',num2str(timept,'%02.0f'),')/sum(elemAreas_',num2str(timept,'%02.0f'),'),sum(elemCents_',num2str(timept,'%02.0f'),'(:,3).*elemAreas_',num2str(timept,'%02.0f'),')/sum(elemAreas_',num2str(timept,'%02.0f'),')];']);
 
 disp('Generating Green''s matrix')
 eval(['[beadu_',num2str(distFromCellBead),'umFromCell_',num2str(timept,'%02.0f'),']=beadDispCalcLinear([',dataNameEXP,'(:,1),',dataNameEXP,'(:,2),',dataNameEXP,'(:,3)],repElems_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),',nodeList_m2r_',num2str(timept,'%02.0f'),',repNodes_',num2str(distFromCellBead),'umFromCell_m2r_',num2str(timept,'%02.0f'),',',dataNameFEM,num2str(timept,'%02.0f'),',elemList_m2r_',num2str(timept,'%02.0f'),');']);
 eval(['clear ',dataNameFEM])
 disp('Saving Green''s matrix')
 eval(['save ''',filestem,'beadu_',num2str(distFromCellBead),'umFromCell_',num2str(timept,'%02.0f'),'.mat'' ''beadu_',num2str(distFromCellBead),'umFromCell_',num2str(timept,'%02.0f'),''';']);
 
%  eval(['load ''',filestem2,'beadu_',num2str(distFromCellBead),'umFromCell_',num2str(timept,'%02.0f'),'.mat''']);
 eval(['[dist_expand]=distVectorExpand(',dataNameEXP,');']);
 if distFromCellBead~=distFromCellSVD
    eval(['beadu_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f'),'=beadu_',num2str(timept,'%02.0f'),'_',num2str(distFromCellBead),'(dist_expand<=',num2str(distFromCellSVD),',:);']);
    eval(['clear beadu_',num2str(distFromCellBead),'umFromCell_',num2str(timept,'%02.0f'),''])
 end
 
 disp('Calculating SVD')
 eval(['[U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',V_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),']=csvd(beadu_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f'),');']);

 disp('Saving SVD data')
 eval(['clear beadu_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f')])
 eval(['save ''',filestem,'U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),'.mat'' ''U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),''';']);
 eval(['save ''',filestem,'s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),'.mat'' ''s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),''';']);
 eval(['save ''',filestem,'V_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),'.mat'' ''V_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),''';']);
 
 disp('Calculating L-curve and computing tractions')
 eval(['[beadDisp_',num2str(timept,'%02.0f'),'] = genYlsq_3D(beadDisplacements_',num2str(timept,'%02.0f'),');']);
 eval(['beadDisp_',num2str(timept,'%02.0f'),' = beadDisp_',num2str(timept,'%02.0f'),'(dist_expand<=',num2str(distFromCellSVD),');'])
 eval(['[reg_corner_',num2str(timept,'%02.0f'),',rho_',num2str(timept,'%02.0f'),',eta_',num2str(timept,'%02.0f'),',reg_param_',num2str(timept,'%02.0f'),'] = l_curve(U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',beadDisp_',num2str(timept,'%02.0f'),',''Tikh'');']);
 eval(['[x_lambda_',num2str(timept,'%02.0f'),',rho_',num2str(timept,'%02.0f'),',eta_',num2str(timept,'%02.0f'),'] = tikhonov(U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',V_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),',beadDisp_',num2str(timept,'%02.0f'),',reg_corner_',num2str(timept,'%02.0f'),');']); 
 eval(['[F_',num2str(timept,'%02.0f'),'] = genFcomp(x_lambda_',num2str(timept,'%02.0f'),');']);
 eval(['[sigN_',num2str(timept,'%02.0f'),',tau_',num2str(timept,'%02.0f'),',Mag_',num2str(timept,'%02.0f'),'] = tractionComp(F_',num2str(timept,'%02.0f'),',elemNormals_',num2str(timept,'%02.0f'),');']);
 eval(['[theta_',num2str(timept,'%02.0f'),']=thetaCalc(F_',num2str(timept,'%02.0f'),',elemAreas_',num2str(timept,'%02.0f'),',elemCents_',num2str(timept,'%02.0f'),');']);
 eval(['[Tractions_',num2str(timept,'%02.0f'),']=[F_',num2str(timept,'%02.0f'),',sigN_',num2str(timept,'%02.0f'),',tau_',num2str(timept,'%02.0f'),',theta_',num2str(timept,'%02.0f'),'];']);
 eval(['beadDispR=U_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),'*diag(s_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),')*V_',num2str(distFromCellSVD),'_',num2str(timept,'%02.0f'),'''*x_lambda_',num2str(timept,'%02.0f'),';']);
 eval(['[bR]=genFcomp(beadDispR);']);
 eval(['[bM]=genFcomp(beadDisp_',num2str(timept,'%02.0f'),');']);
 eval(['bc=',dataNameEXP,'(',dataNameEXP,'(:,8)<=',num2str(distFromCellSVD),',1:3);']);
 eval(['beadDispMR_',num2str(timept,'%02.0f'),'=[bc,bM,bR,bM(:,1)-bR(:,1),bM(:,2)-bR(:,2),bM(:,3)-bR(:,3),sqrt((bM(:,1)-bR(:,1)).^2+(bM(:,2)-bR(:,2)).^2+(bM(:,3)-bR(:,3)).^2)];']);
 
 disp('Saving traction data')
 eval(['save ''',filestem,'Tractions_',num2str(timept,'%02.0f'),'.mat'' ''Tractions_',num2str(timept,'%02.0f'),''';']);
 eval(['save ''',filestem,'beadDisp_',num2str(timept,'%02.0f'),'.mat'' ''beadDisp_',num2str(timept,'%02.0f'),''';']);
  eval(['save ''',filestem,'beadDispMR_',num2str(timept,'%02.0f'),'.mat'' ''beadDispMR_',num2str(timept,'%02.0f'),''';']);
 eval(['[cellNodes_reNum_',num2str(timept,'%02.0f'),',elemConn_reNum_',num2str(timept,'%02.0f'),']=tecPlotRenumber(nodeList_m2r_',num2str(timept,'%02.0f'),',voidElemConn_m2r_',num2str(timept,'%02.0f'),');']);
 eval(['save ''',filestem,'cellNodes_reNum_',num2str(timept,'%02.0f'),'.mat'' ''cellNodes_reNum_',num2str(timept,'%02.0f'),''';']);
 eval(['save ''',filestem,'elemConn_reNum_',num2str(timept,'%02.0f'),'.mat'' ''elemConn_reNum_',num2str(timept,'%02.0f'),''';']);
 eval(['clear U_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f')])
 eval(['clear s_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f')])
 eval(['clear V_',num2str(distFromCellSVD),'umFromCell_',num2str(timept,'%02.0f')])
end

 eval(['fid = fopen([filestem,''tractions_TecplotFormat.dat''], ''w'')']);
 fprintf(fid,['TITLE =  "Computed cell Forces"\r\nVARIABLES = "X", "Y", "Z", "Fx", "Fy", "Fz", "M", "sigN","tau","theta"']);
 for timept=time_pt:num_time_pts
 eval(['[numCellNodes,dimCellNodes]=size(cellNodes_reNum_',num2str(timept,'%02.0f'),');']);
 eval(['[numCellElems,dimCellElems]=size(elemConn_reNum_',num2str(timept,'%02.0f'),');']);
 eval(['[numTractions,dimTractions]=size(Tractions_',num2str(timept,'%02.0f'),');'])
 fprintf(fid,['\r\nZONE T="Computed Cell Forces at t=',num2str(timept),'", NODES=',num2str(numCellNodes),', Elements=',num2str(numCellElems),', DATAPACKING=BLOCK, SolutionTime=',num2str(timept),' VARLOCATION=([4-10]=CellCentered), ZONETYPE=FETRIANGLE\r\n']);
 
 for k=2:dimCellNodes
     for j=1:numCellNodes
         eval(['fprintf(fid,mat2str(cellNodes_reNum_',num2str(timept,'%02.0f'),'(j,k)));']);
         fprintf(fid,' ');
         if mod(j,50)==0
             fprintf(fid,'\r\n')
         end
     end
     fprintf(fid,'\r\n')
 end
 
 fprintf(fid,'\r\n\r\n')
 for k=1:dimTractions
     for j=1:numTractions
         eval(['fprintf(fid,mat2str(Tractions_',num2str(timept,'%02.0f'),'(j,k)));']);
         fprintf(fid,' ');
         if mod(j,50)==0
             fprintf(fid,'\r\n')
         end
     end
     fprintf(fid,'\r\n')
 end
 
 fprintf(fid,'\r\n')
 for j=1:numCellElems
     for k=1:dimCellElems
         eval(['fprintf(fid,mat2str(elemConn_reNum_',num2str(timept,'%02.0f'),'(j,k)));']);
         fprintf(fid,'   ')
         if k==3
             fprintf(fid,'\r\n');
         end
     end
 end
 end
 fclose(fid)
 clear all