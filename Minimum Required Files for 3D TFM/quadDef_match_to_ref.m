function [y,CENTROIDSR] = quadDef_match_to_ref(matchedCoordinates,CENTROIDS)
%This function 

%Input argument definitions

%Output argument definitions

%history = history of y, each column corresponding to one iteration of the
%solution process

% MC=cell2mat(matchedCoordinates);
MC=matchedCoordinates;
CENTSPre=[MC(:,6:8)];
CENTSPost=[MC(:,2:4)];

b0=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0];

options = optimset('Display', 'iter','GradObj','off','DerivativeCheck', 'off','UseParallel','Never', 'LargeScale', 'off','MaxFunEvals',100000,'TolFun', 1E-12,'TolX',1E-12,'TolCon',1E-12); %Define solver options
history.x = b0; % Set up shared variables with OUTFUN
% tic %Start solution timer
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
[y,resnorm,history,exitflag] = fmincon(@objective, b0,A,b,Aeq,beq,lb,ub,@nonlcon,options); %Execute solution process

% CS=cell2mat(CENTROIDS);
CS=CENTROIDS;
CSR=[CS(:,2:4)];
TQ=[y(1),y(2),y(3);y(4),y(5),y(6);y(7),y(8),y(9)]; %Quadratic Terms of single variable
TQM=[y(10),y(11),y(12);y(13),y(14),y(15);y(16),y(17),y(18)]; %Mixed Quadratic Terms
TA=[y(19),y(20),y(21),y(22);y(23),y(24),y(25),y(26);y(27),y(28),y(29),y(30);,0,0,0,1]; %Affine Terms
affCorr=TA*[CSR';ones(1,length(CS))];

CSRA=TQ*CSR'.^2+TQM*[CSR(:,1)'.*CSR(:,2)';CSR(:,2)'.*CSR(:,3)';CSR(:,1)'.*CSR(:,3)']+affCorr(1:3,:);
% CENTROIDSR=mat2cell([CS(:,1),CSRA(1:3,:)'],ones(length(CS),1),4);
CENTROIDSR=[CS(:,1),CSRA(1:3,:)'];
    
    function [err] = objective(b)
        %Calculate the residual to be minized
        TQ=[b(1),b(2),b(3);b(4),b(5),b(6);b(7),b(8),b(9)]; %Quadratic Terms of single variable
        TQM=[b(10),b(11),b(12);b(13),b(14),b(15);b(16),b(17),b(18)]; %Mixed Quadratic Terms
        TA=[b(19),b(20),b(21),b(22);b(23),b(24),b(25),b(26);b(27),b(28),b(29),b(30);,0,0,0,1]; %Affine Terms
        affCorr=TA*[CENTSPost';ones(1,length(MC))];
        
        adjustCoords=TQ*CENTSPost'.^2+TQM*[CENTSPost(:,1)'.*CENTSPost(:,2)';CENTSPost(:,2)'.*CENTSPost(:,3)';CENTSPost(:,1)'.*CENTSPost(:,3)']+affCorr(1:3,:);
        
        err = norm(CENTSPre'-adjustCoords,'fro');


% function stop = outfun(x,optimValues,state)
% %Output function to store the solution history after each iteration
% stop = false;
%    switch state
%        case 'iter'
%            % Concatenate current point with history
%            history.x = [history.x, x];
%        otherwise
%    end
% end
    end

    function [c,ceq]=nonlcon(b)
        c=[];
        ceq=[];
    end
% toc %Stop solution timer
end