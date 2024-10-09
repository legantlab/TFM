for ii=1:1:1110
    T=zeros(1110,1);
    T(ii)=1;
    temp=G(:,ii);
    quiver(displacedPtsMat(1:3:end,1),displacedPtsMat(1:3:end,2),temp(1:3:end),temp(2:3:end),0,'b')
    hold on
    scatter(displacedPtsMat(:,1),displacedPtsMat(:,2))
    hold off
    axis([-50,50,-100,100])
end