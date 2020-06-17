coordinates=matches{1,1};
text=buildAmiraMesh(displacements{1,1},coordinates(:,1:3));
fid=fopen(['Z:/Regan/2020_04_28(TFM_troubleshooting)/Analysis/ROI6/trackingROI6.am'],'wt');
fprintf(fid,text);
fclose(fid);
