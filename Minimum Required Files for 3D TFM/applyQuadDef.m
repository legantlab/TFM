function [CENTROIDSR]=applyQuadDef(y,CENTROIDS)
CS=CENTROIDS;
CSR=[CS(:,2:4)];
TQ=[y(1),y(2),y(3);y(4),y(5),y(6);y(7),y(8),y(9)]; %Quadratic Terms of single variayle
TQM=[y(10),y(11),y(12);y(13),y(14),y(15);y(16),y(17),y(18)]; %Mixed Quadratic Terms
TA=[y(19),y(20),y(21),y(22);y(23),y(24),y(25),y(26);y(27),y(28),y(29),y(30);,0,0,0,1]; %Affine Terms
affCorr=TA*[CSR';ones(1,length(CS))];

CSRA=TQ*CSR'.^2+TQM*[CSR(:,1)'.*CSR(:,2)';CSR(:,2)'.*CSR(:,3)';CSR(:,1)'.*CSR(:,3)']+affCorr(1:3,:);
% CENTROIDSR=mat2cell([CS(:,1),CSRA(1:3,:)'],ones(length(CS),1),4);
CENTROIDSR=[CS(:,1),CSRA(1:3,:)'];
end