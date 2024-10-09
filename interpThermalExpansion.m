function [newDisps] = interpThermalExpansion(beadDisps, refFrame)
%interpThermalExpansion Interpolates swelling profile to each timeframe and
%removes swelling from matches. 
%   Assumes interpolated matches to the same grid, then computes average
%   displacement profile from each time point, then removes average from
%   all time points. Optional graphing of results. If refFrame is nonzero,
%   then assume there is a frame with no cell that can be used instead of
%   averaging. 

%Compute average displacement profile
if refFrame ~= 0
    %Remove average from each disp
    newDisps = cell(length(beadDisps),1);
    avgDisp = beadDisps{refFrame}(:,4:6);
    for i = 1:length(newDisps)
        curDisps = beadDisps{i};
        curDisps(:,4:6) = curDisps(:,4:6) - avgDisp;
        newDisps{i} = curDisps;
    end


else
    allDispsX = zeros(length(beadDisps{1}(:,1)),length(beadDisps));
    allDispsY = zeros(length(beadDisps{1}(:,1)),length(beadDisps));
    allDispsZ = zeros(length(beadDisps{1}(:,1)),length(beadDisps));
    avgDisp = zeros(length(beadDisps{1}(:,1)),3);
    
    for i = 1:length(beadDisps)
        allDispsX(:,i) = beadDisps{i}(:,4);
        allDispsY(:,i) = beadDisps{i}(:,5);
        allDispsZ(:,i) = beadDisps{i}(:,6);
    end
    
    avgDisp(:,1) = mean(allDispsX,2);
    avgDisp(:,2) = mean(allDispsY,2);
    avgDisp(:,3) = mean(allDispsZ,2);
    
    %Remove average from each disp
    newDisps = cell(length(beadDisps),1);
    for i = 1:length(newDisps)
        curDisps = beadDisps{i};
        curDisps(:,4:6) = curDisps(:,4:6) - avgDisp;
        newDisps{i} = curDisps;
    end
end
%% Plotting
cents = beadDisps{1}(:,1:3);

figure
quiver3(cents(:,1),cents(:,2),cents(:,3), avgDisp(:,1), avgDisp(:,2), avgDisp(:,3))
axis equal
title('Average Displacement Profile')

figure
quiver3(cents(:,1),cents(:,2),cents(:,3), newDisps{1}(:,4), newDisps{1}(:,5), newDisps{1}(:,6),1)
axis equal
title('Example Swelling Removed Displacments')


end