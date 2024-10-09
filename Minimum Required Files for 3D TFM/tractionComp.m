function [sigN,tau,Mag] = tractionComp(F,elemNormals)
sigN=zeros(length(F),1);
tau=zeros(length(F),1);
Mag=zeros(length(F),1);
t=zeros(1,3);
for i=1:length(F)
    sigN(i)=dot(F(i,1:3),elemNormals(i,:));
    t=F(i,1:3)-sigN(i)*elemNormals(i,:);
    tau(i)=sqrt(t(1).^2+t(2).^2+t(3).^2);
    Mag(i)=sqrt(F(i,1).^2+F(i,2).^2+F(i,3).^2);
end