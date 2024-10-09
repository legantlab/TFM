function [residual] = computeResNorm(residual,curVec,curDefVec, numNeighbors, numPerms, PV)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    for j=1:numNeighbors
        tempJ = zeros(1,numPerms);
        
        for k=1:numPerms
            curPV = PV(k,:);
            %residual(k,j)=norm(curVec-curDefVec((j-1)*3*numVectorsDef+curPV)); %Residual needs to be resliced in some other manner to get parfor to work
            tempJ(k) = norm(curVec-curDefVec((j-1)*3*numVectorsDef+curPV));
        end
        residual(:,j) = tempJ;
    end

end