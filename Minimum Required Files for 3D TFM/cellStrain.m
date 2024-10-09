function [props,PrincS,nodeLocs]=cellStrain(vecOutputs,elemCents,jm)

xmin=min(vecOutputs(:,1));
xmax=max(vecOutputs(:,1));
ymin=min(vecOutputs(:,2));
ymax=max(vecOutputs(:,2));
zmin=min(vecOutputs(:,3));
zmax=max(vecOutputs(:,3));

[Xi,Yi,Zi]=meshgrid(xmin:4:xmax,ymin:4:ymax,zmin:4:zmax);
[CXi,CYi,CZi]=meshgrid(xmin+2:4:xmax-2,ymin+2:4:ymax-2,zmin+2:4:zmax-2);
t=1;
[xdim,ydim,zdim]=size(Xi);
nodeLocs=zeros(xdim*ydim*zdim,4);
for i=1:xdim
    for j=1:ydim
        for k=1:zdim
            nodeLocs(t,1)=t;
            nodeLocs(t,2)=Xi(i,j,k);
            nodeLocs(t,3)=Yi(i,j,k);
            nodeLocs(t,4)=Zi(i,j,k);
            t=t+1;
        end
    end
end

ui=griddata3(vecOutputs(:,1),vecOutputs(:,2),vecOutputs(:,3),vecOutputs(:,4),Xi,Yi,Zi);
vi=griddata3(vecOutputs(:,1),vecOutputs(:,2),vecOutputs(:,3),vecOutputs(:,5),Xi,Yi,Zi);
wi=griddata3(vecOutputs(:,1),vecOutputs(:,2),vecOutputs(:,3),vecOutputs(:,6),Xi,Yi,Zi);


[props,PrincS]=mechCalcCubeNoFilt(Xi,Yi,Zi,ui,vi,wi,CXi,CYi,CZi);
[dist]=distanceCalcDisp([[1:1:length(PrincS)]',PrincS(:,1:3)],elemCents,6,jm);
PrincS(:,5)=dist;