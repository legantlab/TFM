function [yvals] = plotBoxCharts(xVals,yVals,grid, swarmFlag,options)
%plotBoxCharts Plots a bar chart with xVal labels using a matrix of yVals
%for each xVal. 
%   Input:
%       xVals: 1xn matrix containing labels for x axis box plot
%       yVals: nxm matrix containing n rows of data for each xVal and m
%       columns for replicates of data for each point. 
%       swarmFlag: If = 1, plot scatter points on top of bars. 
%       grid: Values for which to compute average/box chart bins
%       options: Options for plotting
%   Output:
%       Figure: Boxplot figure with swarmchart overlay
hold on 
roundDiff = abs(grid(1) - grid(2));
%Disassemble the x and yVal cells
gridBins = cell(length(grid),1);
yvals = zeros(length(grid),length(xVals));
for i = 1:length(xVals)
    bins = cell(length(grid),1);
    curX = xVals{i};
    curY = yVals{i};

    curX = round(curX/roundDiff) *roundDiff; %Rounds xvalues to grid. 
    [~,binID] = ismember(curX, grid);
    for j = 1:length(binID)
        curVal = binID(j);
        if curVal ~= 0
            bins{binID(j)} = [bins{binID(j)}, curY(j)];
        end
    end
    
    %now just loop through the bins and get average measurement to add to
    %gridbins. 
    for k = 1:length(bins)
        gridBins{k} = [gridBins{k}, mean(bins{k})];
    end
end

%Loop through each grid and plot associated yVals 
    for i = 1:length(grid)
        curX = grid(i) + options.offset;
        curY = gridBins{i};
        yvals(i,:) = curY;
        boxchart(curX*ones(length(curY),1), curY,"BoxFaceColor",options.boxColor,'MarkerStyle','none','BoxWidth',2)
        if swarmFlag == 1
            swarmchart(curX, curY,'filled',options.swarmColor,'MarkerFaceAlpha',options.swarmOpacity,'XJitterWidth',.5,'XJitter','rand')

        end
    end
end