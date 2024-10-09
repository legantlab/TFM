function [elem,face,index]=voidFaceFind(elemList,voidElemConn)
t=1;
for i=1:length(voidElemConn)
    a1(i)={find(elemList(:,2)==voidElemConn(i,2))};
    a2(i)={find(elemList(:,3)==voidElemConn(i,2))};
    a3(i)={find(elemList(:,4)==voidElemConn(i,2))};
    a4(i)={find(elemList(:,5)==voidElemConn(i,2))};
    
    b1(i)={find(elemList(:,2)==voidElemConn(i,3))};
    b2(i)={find(elemList(:,3)==voidElemConn(i,3))};
    b3(i)={find(elemList(:,4)==voidElemConn(i,3))};
    b4(i)={find(elemList(:,5)==voidElemConn(i,3))};
    
    c1(i)={find(elemList(:,2)==voidElemConn(i,4))};
    c2(i)={find(elemList(:,3)==voidElemConn(i,4))};
    c3(i)={find(elemList(:,4)==voidElemConn(i,4))};
    c4(i)={find(elemList(:,5)==voidElemConn(i,4))};
    i
end
a1=cell2mat(a1');
a2=cell2mat(a2');
a3=cell2mat(a3');
a4=cell2mat(a4');

b1=cell2mat(b1');
b2=cell2mat(b2');
b3=cell2mat(b3');
b4=cell2mat(b4');

c1=cell2mat(c1');
c2=cell2mat(c2');
c3=cell2mat(c3');
c4=cell2mat(c4');

node1=unique([a1;b1;c1]);
node2=unique([a2;b2;c2]);
node3=unique([a3;b3;c3]);
node4=unique([a4;b4;c4]);

[c12]=intersect(node1,node2);
[c123]=intersect(c12,node3);
[c124]=intersect(c12,node4);
[c13]=intersect(node1,node3);
[c134]=intersect(c13,node4);
[c23]=intersect(node2,node3);
[c234]=intersect(c23,node4);

potentialElems=unique([c123;c124;c134;c234]);%
% potentialElems=unique([node1;node2;node3;node4]);
len=length(voidElemConn);
for j=1:length(voidElemConn)
    for i=1:length(potentialElems)
    if(length(intersect(voidElemConn(j,2:4),elemList(potentialElems(i),2:5)))==3)
        elem(t,1)=potentialElems(i);
        [int,ai,bi]=intersect(voidElemConn(j,2:4),elemList(potentialElems(i),2:5));
        if (length(intersect(bi,[1,2,3]))==3)
            face(t,1)=1;
        end
        if (length(intersect(bi,[1,2,4]))==3)
            face(t,1)=2;
        end
        if (length(intersect(bi,[2,3,4]))==3)
            face(t,1)=3;
        end
        if (length(intersect(bi,[1,3,4]))==3)
            face(t,1)=4;
        end
        index(t,1)=j;
        t=t+1;
    end
    end
    j/len
end
