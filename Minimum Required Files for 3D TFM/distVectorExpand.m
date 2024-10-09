function [dist_expand]=distVectorExpand(vectorOutputs)
 t=1;
 dist_expand=zeros(length(vectorOutputs)*3,1);
 for i=1:length(vectorOutputs)
     dist_expand(t:t+2,1)=vectorOutputs(i,8);
     t=t+3;
 end