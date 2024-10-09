function [averages] = avgBins(data,positions, sampleGrid)
%avgBins Averages data within bins specified by sampleGrid input
%   Detailed explanation goes here
averages = zeros(length(sampleGrid),1);
for i = 1:length(averages)
    if i == 1 
        averages(i) = mean(data(positions<sampleGrid(i)),'omitnan');
    elseif i == length(averages)
        averages(i) = mean(data(positions>sampleGrid(length(averages))),'omitnan');
    else
        averages(i) = mean(data(positions>sampleGrid(i-1) & positions<sampleGrid(i)),'omitnan');
    end

end

end