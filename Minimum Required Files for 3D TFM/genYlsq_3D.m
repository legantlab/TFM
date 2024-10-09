function [y] = genYlsq_3D(beadDisplacements)

j=1;
y=zeros(3*length(beadDisplacements),1);
for i=1:length(beadDisplacements)
    y(j)=beadDisplacements(i,6)-beadDisplacements(i,2);
    y(j+1)=beadDisplacements(i,7)-beadDisplacements(i,3);
    y(j+2)=beadDisplacements(i,8)-beadDisplacements(i,4);
    j=j+3;
end


