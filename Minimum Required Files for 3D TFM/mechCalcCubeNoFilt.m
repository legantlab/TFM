function [props,PrincS]=mechCalcCubeNoFilt(xi,yi,zi,ui,vi,wi,cx,cy,cz)
% nodeReject=nodeLocs(elemCheck<1,:);
[xdim,ydim,zdim]=size(xi);
% PrincS=zeros(length(nodeReject),4);
props={};
t=1;
for i=1:xdim-1
    i
    for j=1:ydim-1
        for k=1:zdim-1
            %%Node 1
            N(:,1)=[xi(i,j,k),yi(i,j,k),zi(i,j,k)]';
            U(:,1)=[ui(i,j,k),vi(i,j,k),wi(i,j,k)]';
            %%Node 2
            N(:,4)=[xi(i+1,j,k),yi(i+1,j,k),zi(i+1,j,k)]';
            U(:,4)=[ui(i+1,j,k),vi(i+1,j,k),wi(i+1,j,k)]';
            %%Node 3
            N(:,3)=[xi(i+1,j+1,k),yi(i+1,j+1,k),zi(i+1,j+1,k)]';
            U(:,3)=[ui(i+1,j+1,k),vi(i+1,j+1,k),wi(i+1,j+1,k)]';
            %%Node 4
            N(:,2)=[xi(i,j+1,k),yi(i,j+1,k),zi(i,j+1,k)]';
            U(:,2)=[ui(i,j+1,k),vi(i,j+1,k),wi(i,j+1,k)]';
            %%Node 5
            N(:,5)=[xi(i,j,k+1),yi(i,j,k+1),zi(i,j,k+1)]';
            U(:,5)=[ui(i,j,k+1),vi(i,j,k+1),wi(i,j,k+1)]';
            %%Node 6
            N(:,6)=[xi(i,j+1,k+1),yi(i,j+1,k+1),zi(i,j+1,k+1)]';
            U(:,6)=[ui(i,j+1,k+1),vi(i,j+1,k+1),wi(i,j+1,k+1)]';
            %%Node 7
            N(:,7)=[xi(i+1,j+1,k+1),yi(i+1,j+1,k+1),zi(i+1,j+1,k+1)]';
            U(:,7)=[ui(i+1,j+1,k+1),vi(i+1,j+1,k+1),wi(i+1,j+1,k+1)]';
            %%Node 8
            N(:,8)=[xi(i+1,j,k+1),yi(i+1,j,k+1),zi(i+1,j,k+1)]';
            U(:,8)=[ui(i+1,j,k+1),vi(i+1,j,k+1),wi(i+1,j,k+1)]';
            
            epsilon = strain(0,0,0,abs(N(1,1)-N(1,2)),abs(N(2,1)-N(2,4)),abs(N(3,1)-N(3,5)),U');
               
            PrincS(t,1)=cx(i,j,k);
            PrincS(t,2)=cy(i,j,k);
            PrincS(t,3)=cz(i,j,k);
            props(t,1)={[cx(i,j,k),cy(i,j,k),cz(i,j,k)]};
            
%             if  (~inVoid(N(:,1),N(:,2),N(:,3),N(:,4),N(:,5),N(:,6),N(:,7),N(:,8),nodeReject) && ~max(max(isnan(epsilon))))
                props(t,2)={epsilon};
                PrincS(t,4)=principalStrain(epsilon);
                t=t+1;
%             else
%                 props(t,2)={-1};
%                 PrincS(t,4)=-1;
%                 t=t+1;
%             end
            %             N=[xi(i),yi(i),zi(i);xi(i),yi(i+1),zi(i);xi(i+1),yi(i+1),zi(i);xi(i+1),yi(i),zi(i);...
            %                 xi(i),yi(i),zi(i+1);xi(i),yi(i+1),zi(i+1);xi(i+1),yi(i+1)
            %                 ,zi(i+1);xi(i+1),yi(i),zi(i+1)]';
        end
    end
end
end

    function [found] = inVoid(N1,N2,N3,N4,N5,N6,N7,N8,nodeReject)
        f=0;
        for i=1:length(nodeReject)
            f=f+(sum((nodeReject(i,:)==N1'))==3);
            f=f+(sum((nodeReject(i,:)==N2'))==3);
            f=f+(sum((nodeReject(i,:)==N3'))==3);
            f=f+(sum((nodeReject(i,:)==N4'))==3);
            f=f+(sum((nodeReject(i,:)==N5'))==3);
            f=f+(sum((nodeReject(i,:)==N6'))==3);
            f=f+(sum((nodeReject(i,:)==N7'))==3);
            f=f+(sum((nodeReject(i,:)==N8'))==3);
        end
        found=(f>0);
    end