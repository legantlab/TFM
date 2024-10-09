function [props_driftCorr,PrincS_driftCorr]=driftStrain(beadCoords,varargin)
xmin=min(beadCoords(:,1));
xmax=max(beadCoords(:,1));
ymin=min(beadCoords(:,2));
ymax=max(beadCoords(:,2));
zmin=min(beadCoords(:,3));
zmax=max(beadCoords(:,3));

[Xi,Yi,Zi]=meshgrid(xmin:8:xmax,ymin:8:ymax,zmin:8:zmax);
[CXi,CYi,CZi]=meshgrid(xmin+4:8:xmax-4,ymin+4:8:ymax-4,zmin+4:8:zmax-4);
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

shiftLocs=nodeLocs;
for i=1:length(varargin)
[shiftLocs]=applyQuadDef(varargin{i},shiftLocs);
end
ui=zeros(xdim,ydim,zdim);
vi=zeros(xdim,ydim,zdim);
wi=zeros(xdim,ydim,zdim);
t=1;
for i=1:xdim
    for j=1:ydim
        for k=1:zdim
            ui(i,j,k)=nodeLocs(t,2)-shiftLocs(t,2);
            vi(i,j,k)=nodeLocs(t,3)-shiftLocs(t,3);
            wi(i,j,k)=nodeLocs(t,4)-shiftLocs(t,4);
            t=t+1;
        end
    end
end


[props_driftCorr,PrincS_driftCorr]=mechCalcCubeNoFilt(Xi,Yi,Zi,ui,vi,wi,CXi,CYi,CZi);