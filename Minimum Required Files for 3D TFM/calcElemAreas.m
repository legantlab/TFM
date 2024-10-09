function [elemAreas]=calcElemAreas(nodeList,voidElemConn)
    elemAreas=zeros(length(voidElemConn),1);
    for i=1:length(voidElemConn)
        N(:,1)=nodeList(voidElemConn(i,2),2:4)';
        N(:,2)=nodeList(voidElemConn(i,3),2:4)';
        N(:,3)=nodeList(voidElemConn(i,4),2:4)';
        V1=N(:,2)-N(:,1);
        V2=N(:,1)-N(:,3);
        elemAreas(i)=.5*norm(cross(V1,V2));
%         elemAreas(i)=.5*sqrt((N(1,1)*N(2,2)-N(1,1)*N(3,2)-N(1,2)*N(2,1)+N(1,2)*N(3,1)+N(2,1)*N(3,2)-N(3,1)*N(2,2)).^2+...
%             (N(1,2)*N(2,3)-N(1,2)*N(3,3)-N(1,3)*N(2,2)+N(1,3)*N(3,2)+N(2,2)*N(3,3)-N(3,2)*N(2,3)).^2+...
%             (N(1,3)*N(2,1)-N(1,3)*N(3,1)-N(1,1)*N(2,3)+N(1,1)*N(3,3)+N(2,3)*N(3,1)-N(3,3)*N(2,1)).^2);
    end