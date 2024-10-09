function [outputargs] = processBeadDefs_match_to_ref_2010_09_09(filestem,inputfile2D,relaxed_imageNameGFP,relaxed_imageNameRFP,stressed_imageNameGFP,stressed_imageNameRFP,dimensions,pxlscaling,multiple_time_pts,varargin)
%This function loads image stacks from the stressed and relaxed
%configurations, locates bead centroids and computes bead trajectories

%Input argument definitions
% filestem = a string indicating the directory where the image stackes are located
% inputfile2D = a string indicating the name of the inputfile containing the 2D node coordinates and element connectivity
% relaxed_imageNameGFP = the name of the relaxed image stack on the GFP channel
% stresed_imageNameGFP = the name of the stressed image stack on the GFP channel
% relaxed_imageNameGFP = the name of the relaxed image stack on the RFP channel
% stresed_imageNameGFP = the name of the stressed image stack on the RFP channel
% dimensions = a 1 by 3 vector indicating the [xdimensions, y dimensions and z dimensions] of the image stackes in pixels
% pxlscaling = a 1 by 3 vector indicating the [xscale,yscale,zscale] of the voxels in microns
% multiple_time_pts = a boolean indicating whether or not to analyze multiple time points (0 is single time point, 1 is multiple timepoints)
% varargin
%     1) index of the first time point to analyze
%     2) number of additional time points to analyze


jm = findresource('jobmanager'); %Initialize a jobmanager to manage parallel processing
numproc=6; %number of processors to use
fieldCutOff=0; %a value indicating how many microns from the cell to cutoff near and far field data for drift/swelling correction
outputargs=1; %initialize the outputargments

if multiple_time_pts==1
    time_pt=varargin{1};
    num_time_pts=varargin{2};
    disp(['Will process ', num2str(num_time_pts-time_pt+1),' time points'])
else
    time_pt=1;
    num_time_pts=1;
    disp('Processing single time point')
end

disp('Load relaxed image and locate bead centroids')
% eval(strcat('load ''',filestem,relaxed_imageNameRFP,'_relaxedCS.mat'''));
% eval(strcat('load ''',filestem,relaxed_imageNameGFP,'_relaxedCS.mat'''));
[relaxed_dataGFP]=loadData_Crocker(filestem,relaxed_imageNameGFP,dimensions);
[relaxed_BPdataGFP]=bpass3dMB(relaxed_dataGFP, [1,1,1], [3,3,3],[0,0]);
clear relaxed_dataGFP
[relaxedSpotsGFP]=feature3dMB_parallel(relaxed_BPdataGFP, 3,[3,3,3],[dimensions],numproc,jm,[1,0,0],1.5);
relaxedSpotsGFP=double(relaxedSpotsGFP);
clear relaxed_BPdataGFP

[relaxed_dataRFP]=loadData_Crocker(filestem,relaxed_imageNameRFP,dimensions);
[relaxed_BPdataRFP]=bpass3dMB(relaxed_dataRFP, [1,1,1], [3,3,3],[0,0]);
clear relaxed_dataRFP
[relaxedSpotsRFP]=feature3dMB_parallel(relaxed_BPdataRFP, 3,[3,3,3],[dimensions],numproc,jm,[1,0,0],1.5);
relaxedSpotsRFP=double(relaxedSpotsRFP);
clear relaxed_BPdataRFP
% 
% 
% 
relaxedCS_GFP=[[1:1:length(relaxedSpotsGFP)]',relaxedSpotsGFP(:,1)*pxlscaling(1)-pxlscaling(1),relaxedSpotsGFP(:,2)*pxlscaling(2)-pxlscaling(2),relaxedSpotsGFP(:,3)*pxlscaling(3)-pxlscaling(3)];
relaxedCS_RFP=[[1:1:length(relaxedSpotsRFP)]',relaxedSpotsRFP(:,1)*pxlscaling(1)-pxlscaling(1),relaxedSpotsRFP(:,2)*pxlscaling(2)-pxlscaling(2),relaxedSpotsRFP(:,3)*pxlscaling(3)-pxlscaling(3)];
relaxedCS=[relaxedCS_GFP;relaxedCS_RFP];
eval(['save ''',filestem,relaxed_imageNameGFP,'_relaxedCS.mat'' ''relaxedCS_GFP'';']);
eval(['save ''',filestem,relaxed_imageNameRFP,'_relaxedCS.mat'' ''relaxedCS_RFP'';']);
eval(['save ''',filestem,'relaxedCS.mat'' ''relaxedCS'';']);

for k=time_pt:num_time_pts
    disp(strcat('Load stressed image and locate bead centroids for time point_ ',num2str(k,'%02.0f'),''))
%     eval(strcat('load ''',filestem,stressed_imageNameRFP,'stressedCS_',num2str(k,'%02.0f'),'.mat'''));
%     eval(strcat('load ''',filestem,stressed_imageNameGFP,'stressedCS_',num2str(k,'%02.0f'),'.mat'''));
    eval(strcat('[stressed_dataGFP_',num2str(k,'%02.0f'),']=loadData_Crocker(''',filestem,''',''',stressed_imageNameGFP,num2str(k,'%02.0f'),''',[',num2str(dimensions),']);'));
    eval(strcat('[stressed_BPdataGFP_',num2str(k,'%02.0f'),']=bpass3dMB(stressed_dataGFP_',num2str(k,'%02.0f'),', [1,1,1], [3,3,3],[0,0]);'));
    eval(strcat('clear stressed_dataGFP_',num2str(k,'%02.0f')));
    eval(strcat('[stressedSpotsGFP_',num2str(k,'%02.0f'),']=feature3dMB_parallel(stressed_BPdataGFP_',num2str(k,'%02.0f'),', 3,[3,3,3], [',num2str(dimensions),'],',num2str(numproc),',jm,[1,0,0],1.5);'));
    eval(strcat('stressedSpotsGFP_',num2str(k,'%02.0f'),'=double(stressedSpotsGFP_',num2str(k,'%02.0f'),');'));
    eval(strcat('clear stressed_BPdataGFP_',num2str(k,'%02.0f')));
    
    eval(strcat('[stressed_dataRFP_',num2str(k,'%02.0f'),']=loadData_Crocker(''',filestem,''',''',stressed_imageNameRFP,num2str(k,'%02.0f'),''',[',num2str(dimensions),']);'));
    eval(strcat('[stressed_BPdataRFP_',num2str(k,'%02.0f'),']=bpass3dMB(stressed_dataRFP_',num2str(k,'%02.0f'),', [1,1,1], [3,3,3],[0,0]);'));
    eval(strcat('clear stressed_dataRFP_',num2str(k,'%02.0f')));
    eval(strcat('[stressedSpotsRFP_',num2str(k,'%02.0f'),']=feature3dMB_parallel(stressed_BPdataRFP_',num2str(k,'%02.0f'),', 3,[3,3,3], [',num2str(dimensions),'],',num2str(numproc),',jm,[1,0,0],1.5);'));
    eval(strcat('stressedSpotsRFP_',num2str(k,'%02.0f'),'=double(stressedSpotsRFP_',num2str(k,'%02.0f'),');'));
    eval(strcat('clear stressed_BPdataRFP_',num2str(k,'%02.0f')));
    eval(strcat('stressedCS_GFP_',num2str(k,'%02.0f'),'=[[1:1:length(stressedSpotsGFP_',num2str(k,'%02.0f'),')]'',stressedSpotsGFP_',num2str(k,'%02.0f'),'(:,1)*pxlscaling(1)-pxlscaling(1),stressedSpotsGFP_',num2str(k,'%02.0f'),'(:,2)*pxlscaling(2)-pxlscaling(2),stressedSpotsGFP_',num2str(k,'%02.0f'),'(:,3)*pxlscaling(3)-pxlscaling(3)];'));
    eval(strcat('stressedCS_RFP_',num2str(k,'%02.0f'),'=[[1:1:length(stressedSpotsRFP_',num2str(k,'%02.0f'),')]'',stressedSpotsRFP_',num2str(k,'%02.0f'),'(:,1)*pxlscaling(1)-pxlscaling(1),stressedSpotsRFP_',num2str(k,'%02.0f'),'(:,2)*pxlscaling(2)-pxlscaling(2),stressedSpotsRFP_',num2str(k,'%02.0f'),'(:,3)*pxlscaling(3)-pxlscaling(3)];'));
    eval(strcat('stressedCS_',num2str(k,'%02.0f'),'=[stressedCS_GFP_',num2str(k,'%02.0f'),';stressedCS_RFP_',num2str(k,'%02.0f'),'];'));
    eval(['save ''',filestem,stressed_imageNameGFP,'stressedCS_',num2str(k,'%02.0f'),'.mat'' ''stressedCS_GFP_',num2str(k,'%02.0f'),''';']);
    eval(['save ''',filestem,stressed_imageNameRFP,'stressedCS_',num2str(k,'%02.0f'),'.mat'' ''stressedCS_RFP_',num2str(k,'%02.0f'),''';']);
    eval(['save ''',filestem,'stressedCS_',num2str(k,'%02.0f'),'.mat'' ''stressedCS_',num2str(k,'%02.0f'),''';']);
    
    %Load finite element datasets and compute element centroids (files must
    %be in filestem directory and be named accordingly)
    eval(['filetext = fileread(''',filestem,inputfile2D,num2str(k,'%02.0f'),'.inp'');'])
    [data,triElemStart] = regexp(filetext, 'ELEMENT,TYPE=S3', 'match');
    [node_index] = regexp(filetext(1:triElemStart), '\s+\d+[\x2C\s]+[\x2D\dE\.]+[\x2C\s]+[\x2D\dE\.]+[\x2C\s]+[\x2D\dE\.]+[\r\n]?', 'match');
    eval(['nodeList_',num2str(k,'%02.0f'),'=str2num(cell2mat(node_index));']);
    [triElem_index] = regexp(filetext(triElemStart:end), '\s+\d+[\x2C\s]+[\x2D\dE\.]+[\x2C\s]+[\x2D\dE\.]+[\x2C\s]+[\x2D\dE\.]+[\r\n]?', 'match');
    eval(['voidElemConn_',num2str(k,'%02.0f'),'=str2num(cell2mat(triElem_index));']);
%     eval(['load ''',filestem,'nodeList_',num2str(k,'%02.0f'),'.txt''']);
%     eval(['load ''',filestem,'voidElemConn_',num2str(k,'%02.0f'),'.txt''']);
    eval(['[elemCents_',num2str(k,'%02.0f'),']=elemCentFind(nodeList_',num2str(k,'%02.0f'),', voidElemConn_',num2str(k,'%02.0f'),');'])
    
    
    disp('Seperate bead centroids into near and far field for noise correction')
    eval(['[relaxed_distGFP_',num2str(k,'%02.0f'),']=distanceCalcDisp(relaxedCS_GFP,elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
    eval(['[relaxed_distRFP_',num2str(k,'%02.0f'),']=distanceCalcDisp(relaxedCS_RFP,elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
    eval(['[stressed_distGFP_',num2str(k,'%02.0f'),']=distanceCalcDisp(stressedCS_GFP_',num2str(k,'%02.0f'),',elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
    eval(['[stressed_distRFP_',num2str(k,'%02.0f'),']=distanceCalcDisp(stressedCS_RFP_',num2str(k,'%02.0f'),',elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
    
    disp('Perform initial matching and drift correction')
    eval(['[matchedCoordinatesGFP_t',num2str(k,'%02.0f'),'_farfield]=displacementCalc_recursive_rough_2010_09_09(stressedCS_GFP_',num2str(k,'%02.0f'),'(stressed_distGFP_',num2str(k,'%02.0f'),'>fieldCutOff,:),relaxedCS_GFP(relaxed_distGFP_',num2str(k,'%02.0f'),'>(fieldCutOff-10),:),150,3,5,2,numproc,jm);'])
    eval(['[matchedCoordinatesRFP_t',num2str(k,'%02.0f'),'_farfield]=displacementCalc_recursive_rough_2010_09_09(stressedCS_RFP_',num2str(k,'%02.0f'),'(stressed_distRFP_',num2str(k,'%02.0f'),'>fieldCutOff,:),relaxedCS_RFP(relaxed_distRFP_',num2str(k,'%02.0f'),'>(fieldCutOff-10),:),150,3,5,2,numproc,jm);'])
    eval(['[matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield]=[matchedCoordinatesGFP_t',num2str(k,'%02.0f'),'_farfield;matchedCoordinatesRFP_t',num2str(k,'%02.0f'),'_farfield];']);
    eval(['[matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield_filt]=smoothFilt_2010_11_29(matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield,50,2,numproc,jm);'])
    eval(['[pre_matched_dist_',num2str(k,'%02.0f'),']=distanceCalcDisp(matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield_filt(:,1:4),elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
%     eval(['[matchedCoordinates_filt_t',num2str(k,'%02.0f'),'_farfield,misID]=movingAvgFilt(pre_matched_dist_',num2str(k,'%02.0f'),',matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield,64,6);'])
    eval(['[y_t',num2str(k,'%02.0f'),',stressedCS_Q_t',num2str(k,'%02.0f'),'] = quadDef_match_to_ref(matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield_filt,stressedCS_',num2str(k,'%02.0f'),');'])
    eval(['[stressedCS_GFP_Q_t',num2str(k,'%02.0f'),']=applyQuadDef(y_t',num2str(k,'%02.0f'),',stressedCS_GFP_',num2str(k,'%02.0f'),');']);
    eval(['[stressedCS_RFP_Q_t',num2str(k,'%02.0f'),']=applyQuadDef(y_t',num2str(k,'%02.0f'),',stressedCS_RFP_',num2str(k,'%02.0f'),');']);
    eval(strcat('stressedCS_Q_t',num2str(k,'%02.0f'),'=[stressedCS_GFP_Q_t',num2str(k,'%02.0f'),';stressedCS_RFP_Q_t',num2str(k,'%02.0f'),'];'));

    
    disp('Perform final matching')
    eval(['[matchedCoordinatesGFP_Q_t',num2str(k,'%02.0f'),']=displacementCalc_recursive_2010_09_09(stressedCS_GFP_Q_t',num2str(k,'%02.0f'),',relaxedCS_GFP,50,3,5,2,numproc,jm);'])
    eval(['[matchedCoordinatesRFP_Q_t',num2str(k,'%02.0f'),']=displacementCalc_recursive_2010_09_09(stressedCS_RFP_Q_t',num2str(k,'%02.0f'),',relaxedCS_RFP,50,3,5,2,numproc,jm);'])
    eval(['[matchedCoordinates_Q_t',num2str(k,'%02.0f'),']=[matchedCoordinatesGFP_Q_t',num2str(k,'%02.0f'),';matchedCoordinatesRFP_Q_t',num2str(k,'%02.0f'),'];']);
    eval(['[matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt]=smoothFilt_2010_11_29(matchedCoordinates_Q_t',num2str(k,'%02.0f'),',20,2,numproc,jm);'])
    
    disp('Perform final drift correction')
    eval(['[matched_dist_',num2str(k,'%02.0f'),']=distanceCalcDisp(matchedCoordinates_Q_t',num2str(k,'%02.0f'),'(:,1:4),elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
%     eval(['[matchedCoordinates_filt_Q_t',num2str(k,'%02.0f'),',misID]=movingAvgFilt(matched_dist_',num2str(k,'%02.0f'),',matchedCoordinates_Q_t',num2str(k,'%02.0f'),',64,6);']);
    eval(['[matched_dist_filt_',num2str(k,'%02.0f'),']=distanceCalcDisp(matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt(:,1:4),elemCents_',num2str(k,'%02.0f'),',numproc,jm);']);
    eval(['[y2_t',num2str(k,'%02.0f'),',shiftedMatches_t',num2str(k,'%02.0f'),'] = quadDef_match_to_ref(matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt(matched_dist_filt_',num2str(k,'%02.0f'),'>fieldCutOff,:),matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt(:,1:4));'])
    eval(['matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'=[shiftedMatches_t',num2str(k,'%02.0f'),',matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt(:,5:8)];'])
    eval(['[nodeList_Q_t',num2str(k,'%02.0f'),']=applyQuadDef(y_t',num2str(k,'%02.0f'),',nodeList_',num2str(k,'%02.0f'),');']);
    eval(['[nodeList_QQ_t',num2str(k,'%02.0f'),']=applyQuadDef(y2_t',num2str(k,'%02.0f'),',nodeList_Q_t',num2str(k,'%02.0f'),');']);
    eval(['[elemCents_QQ_t',num2str(k,'%02.0f'),']=elemCentFind(nodeList_QQ_t',num2str(k,'%02.0f'),', voidElemConn_',num2str(k,'%02.0f'),');'])
    eval(['[matched_dist_FINAL_',num2str(k,'%02.0f'),']=distanceCalcDisp(matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,1:4),elemCents_QQ_t',num2str(k,'%02.0f'),',numproc,jm);']);

    
    disp('Saving results to file')
    eval(['[props_driftCorr_t',num2str(k,'%02.0f'),',PrincS_driftCorr_t',num2str(k,'%02.0f'),']=driftStrain(matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,2:4),y_t',num2str(k,'%02.0f'),',y2_t',num2str(k,'%02.0f'),');'])
    eval(['save ''',filestem,'PrincS_driftCorr_t',num2str(k,'%02.0f'),'.mat'' ''PrincS_driftCorr_t',num2str(k,'%02.0f'),''';']);
    eval(['save ''',filestem,'quadCorr_t',num2str(k,'%02.0f'),'.mat'' ''y_t',num2str(k,'%02.0f'),''';']);
    eval(['save ''',filestem,'quadCorr2_t',num2str(k,'%02.0f'),'.mat'' ''y2_t',num2str(k,'%02.0f'),''';']);
    eval(['save ''',filestem,'matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield.mat'' ''matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield'';'])
    eval(['save ''',filestem,'matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield_filt.mat'' ''matchedCoordinates_t',num2str(k,'%02.0f'),'_farfield_filt'';'])
    eval(['save ''',filestem,'matchedCoordinates_Q_t',num2str(k,'%02.0f'),'.mat'' ''matchedCoordinates_Q_t',num2str(k,'%02.0f'),''';'])
    eval(['save ''',filestem,'matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt.mat'' ''matchedCoordinates_Q_t',num2str(k,'%02.0f'),'_filt'';'])
    eval(['save ''',filestem,'matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'.mat'' ''matchedCoordinates_QQ_t',num2str(k,'%02.0f'),''';'])
    eval(['vectorOutputsFinal_t',num2str(k,'%02.0f'),'=[matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,2),matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,3),matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,4),matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,6)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,2),matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,7)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,3),matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,8)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,4),sqrt((matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,6)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,2)).^2+(matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,7)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,3)).^2+(matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,8)-matchedCoordinates_QQ_t',num2str(k,'%02.0f'),'(:,4)).^2),matched_dist_FINAL_',num2str(k,'%02.0f'),'];'])
    eval(['save ''',filestem,'vectorOutputsFinal_t',num2str(k,'%02.0f'),'.mat'' ''vectorOutputsFinal_t',num2str(k,'%02.0f'),''';'])

    eval(['[props_cellStrain_',num2str(k,'%02.0f'),',PrincS_cellStrain_',num2str(k,'%02.0f'),',nodeLocs_cellStrain_',num2str(k,'%02.0f'),']=cellStrain(vectorOutputsFinal_t',num2str(k,'%02.0f'),',elemCents_QQ_t',num2str(k,'%02.0f'),',jm);']);
    eval(['save ''',filestem,'props_cellStrain_',num2str(k,'%02.0f'),'.mat'' ''props_cellStrain_',num2str(k,'%02.0f'),''';'])
    eval(['save ''',filestem,'PrincS_cellStrain_',num2str(k,'%02.0f'),'.mat'' ''PrincS_cellStrain_',num2str(k,'%02.0f'),''';'])
    eval(['save ''',filestem,'nodeLocs_cellStrain_',num2str(k,'%02.0f'),'.mat'' ''nodeLocs_cellStrain_',num2str(k,'%02.0f'),''';'])
    eval(['fid = fopen([filestem,inputfile2D,','''',num2str(k,'%02.0f'),''',''_m2r.inp''], ''w'')']);
    fprintf(fid,'**Adjusted cell node coordinates for 2D mesh\r\n**\r\n*NODE');
    eval(['[numCellNodes,dim]=size(nodeList_QQ_t',num2str(k,'%02.0f'),');']);
    for j=1:numCellNodes
        fprintf(fid,'\r\n');
        for p=1:dim
    eval(['fprintf(fid,mat2str(nodeList_QQ_t',num2str(k,'%02.0f'),'(j,p)));']);
    if p~=4
    fprintf(fid,',     ');
    end
        end
    end
    fprintf(fid,'\r\n*ELEMENT,TYPE=S3,ELSET=1');
    fprintf(fid,cell2mat(triElem_index));
    fclose(fid);
end
    
    
    eval(['fid = fopen([filestem,''beadDisplacements_TecplotFormat.dat''], ''w'');']);
for i=time_pt:num_time_pts
eval(['[numBeads,dimVectors]=size(vectorOutputsFinal_t',num2str(i,'%02.0f'),');']);
    fprintf(fid,['TITLE =  "Computed bead displacements for t=',num2str(i,'%02.0f'),'"\r\nVARIABLES = "X", "Y", "Z", "dx", "dy", "dz", "mag", "dist"\r\n']);
    fprintf(fid,['ZONE I=',num2str(numBeads),', SolutionTime=',num2str(i,'%02.0f'),' DATAPACKING=POINT\r\n']);
 for j=1:numBeads
     for k=1:dimVectors
         eval(['fprintf(fid,mat2str(vectorOutputsFinal_t',num2str(i,'%02.0f'),'(j,k)));']);
         fprintf(fid,'  ');
     end
     fprintf(fid,'\r\n');
 end
 fprintf(fid,'\r\n');
end
fclose(fid);
% eval(['dlmwrite(''',filestem,'nodeList_QQ_t',num2str(k,'%02.0f'),'.txt'', nodeList_QQ_t',num2str(k,'%02.0f'),', ''delimiter'', '','', ''newline'',''pc'');'])
end

