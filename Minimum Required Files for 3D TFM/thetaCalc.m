function [theta]=thetaCalc(F,elemAreas,elemCents)
    centMass=[sum(elemCents(:,1).*elemAreas)/sum(elemAreas),sum(elemCents(:,2).*elemAreas)/sum(elemAreas),sum(elemCents(:,3).*elemAreas)/sum(elemAreas)];
    theta=zeros(length(F),1);
    for i=1:length(F)
        theta(i)=acos(dot(F(i,1:3),elemCents(i,:)-centMass)/(norm(F(i,1:3),'fro')*norm(elemCents(i,:)-centMass,'fro')))*180/pi;
    end