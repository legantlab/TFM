function [F] = genFcomp(b)

j=1;
k=2;
m=3;
for i=1:length(b)/3;
    Fx(i)=-b(j);
    Fy(i)=-b(k);
    Fz(i)=-b(m);
    j=j+3;
    k=k+3;
    m=m+3;
end
Fx=Fx';
Fy=Fy';
Fz=Fz';
Mag=sqrt(Fx.^2+Fy.^2+Fz.^2);
F=[Fx,Fy,Fz,Mag];
