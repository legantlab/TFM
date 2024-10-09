function [cellNodes_reNum,elemConn_reNum]=tecPlotRenumber(nodeList,voidElemConn)

    Nodes=voidElemConn(:,2:4);
    cellNodes=unique(Nodes(:));
    cellNodeIndex=zeros(max(cellNodes),1);
    for i=1:length(cellNodes)
        cellNodeIndex(cellNodes(i))=i;
    end
    cellNodes_Full=nodeList(cellNodes,:);
    
    for i=1:length(cellNodes)
        cellNodes_reNum(i,1)=i;
        cellNodes_reNum(i,2:4)=cellNodes_Full(i,2:4);
    end
    
    for i=1:length(voidElemConn)
        elemConn_reNum(i,1)=cellNodes_reNum(cellNodeIndex(voidElemConn(i,2)));
        elemConn_reNum(i,2)=cellNodes_reNum(cellNodeIndex(voidElemConn(i,3)));
        elemConn_reNum(i,3)=cellNodes_reNum(cellNodeIndex(voidElemConn(i,4)));
    end