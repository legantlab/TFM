function [maxeps]=principalStrain(epsilon)
%     eps=[epsilon(1),epsilon(4)/2,epsilon(6)/2;epsilon(4)/2,epsilon(2),eps
%     ilon(5)/2;epsilon(6)/2,epsilon(5)/2,epsilon(3)];
if sum(isnan(epsilon))==0
    princ_eps=eig(epsilon);
    maxeps=max(abs(princ_eps));
else maxeps=-1;
end
end