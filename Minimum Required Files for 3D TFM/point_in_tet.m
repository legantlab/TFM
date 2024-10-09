function [repElems,beadCoords,elemCheck]=point_in_tet(pointList,tetList,nodeList)
[tetCents]=elemCentFind3D(nodeList,tetList);
[rPoints,cPoints]=size(pointList);
[rElems,cElems]=size(tetList);
elemCheck=zeros(rPoints,1);
t=1;
for i=1:rPoints
    dist=sum([(tetCents(:,1)-pointList(i,1)).^2,(tetCents(:,2)-pointList(i,2)).^2,(tetCents(:,3)-pointList(i,3)).^2],2).^0.5;
    [val,ind]=sort(dist);
    flag=0;
    k=1;
    while (~flag && k<=rElems)
%         k
        V1=nodeList(tetList(ind(k),2),2:4);
        V2=nodeList(tetList(ind(k),3),2:4);
        V3=nodeList(tetList(ind(k),4),2:4);
        V4=nodeList(tetList(ind(k),5),2:4);

        test1=det([V1,1;V2,1;V3,1;V4,1]);
        test2=det([pointList(i,:),1;V2,1;V3,1;V4,1]);
        if sign(test1)==sign(test2)
            test3=det([V1,1;pointList(i,:),1;V3,1;V4,1]);
            if sign(test3)==sign(test1)
                test4=det([V1,1;V2,1;pointList(i,:),1;V4,1]);
                if sign(test4)==sign(test1)
                    test5=det([V1,1;V2,1;V3,1;pointList(i,:),1]);
                    if sign(test5)==sign(test1) flag=1; end
                end
            end
        end
        k=k+1;
    end
    if flag
        repElems(t)=ind(k-1);
        beadCoords(t,:)=pointList(i,:);
        elemCheck(i)=1;
        t=t+1;
%         scatter3(nodeList(tetList(ind(k-1),2),2),nodeList(tetList(ind(k-1),2),3),nodeList(tetList(ind(k-1),2),4),'b.')
%         hold on
%         scatter3(nodeList(tetList(ind(k-1),3),2),nodeList(tetList(ind(k-1),3),3),nodeList(tetList(ind(k-1),3),4),'b.')
%         scatter3(nodeList(tetList(ind(k-1),4),2),nodeList(tetList(ind(k-1),4),3),nodeList(tetList(ind(k-1),4),4),'b.')
%         scatter3(nodeList(tetList(ind(k-1),5),2),nodeList(tetList(ind(k-1),5),3),nodeList(tetList(ind(k-1),5),4),'b.')
%         scatter3(pointList(i,1),pointList(i,2),pointList(i,3),'ro')
    end
end