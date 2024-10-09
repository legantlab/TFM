function [elemNormals]=calcElemNormals(nodeList,voidElemConn)
    elemNormals=zeros(length(voidElemConn),3);
    for i=1:length(voidElemConn)
        N(:,1)=nodeList(voidElemConn(i,2),2:4)';
        N(:,2)=nodeList(voidElemConn(i,3),2:4)';
        N(:,3)=nodeList(voidElemConn(i,4),2:4)';
        n=(cross((N(:,2)-N(:,1))',(N(:,3)-N(:,1))'));
        elemNormals(i,:)=n/norm(n);
    end