function [points]=dispInterpGrid_90um(vectorOutputsFinal)

x1=[1:1.8:91];
y1=[1:1.8:91];
z1=[0.5:2:4.5];
x2=[1.9:1.8:90.1];
y2=[1.9:1.8:90.1];
z2=[1.5:2:3.5];
[X1,Y1,Z1]=meshgrid(x1,y1,z1);
[X2,Y2,Z2]=meshgrid(x2,y2,z2);

u1=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,4),X1,Y1,Z1);
u1NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,4),X1,Y1,Z1,'nearest');
u1(isnan(u1))=u1NN(isnan(u1));
u2=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,4),X2,Y2,Z2);
u2NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,4),X2,Y2,Z2,'nearest');
u2(isnan(u2))=u2NN(isnan(u2));

v1=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,5),X1,Y1,Z1);
v1NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,5),X1,Y1,Z1,'nearest');
v1(isnan(v1))=v1NN(isnan(v1));
v2=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,5),X2,Y2,Z2);
v2NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,5),X2,Y2,Z2,'nearest');
v2(isnan(v2))=v2NN(isnan(v2));


w1=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,6),X1,Y1,Z1);
w1NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,6),X1,Y1,Z1,'nearest');
w1(isnan(w1))=w1NN(isnan(w1));
w2=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,6),X2,Y2,Z2);
w2NN=griddata(vectorOutputsFinal(:,1),vectorOutputsFinal(:,2),vectorOutputsFinal(:,3),vectorOutputsFinal(:,6),X2,Y2,Z2,'nearest');
w2(isnan(w2))=w2NN(isnan(w2));

t=1;
for i=1:length(x1)
for j=1:length(y1)
for k=1:length(z1)
points(t,:)=[X1(i,j,k),Y1(i,j,k),Z1(i,j,k),u1(i,j,k),v1(i,j,k),w1(i,j,k)];
t=t+1;
end
end
end
for i=1:length(x2)
for j=1:length(y2)
for k=1:length(z2)
points(t,:)=[X2(i,j,k),Y2(i,j,k),Z2(i,j,k),u2(i,j,k),v2(i,j,k),w2(i,j,k)];
t=t+1;
end
end
end

points(:,7)=sqrt(points(:,4).^2+points(:,5).^2+points(:,6).^2);
points(:,8)=points(:,3);