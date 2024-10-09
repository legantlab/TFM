function [nodeLocation, nodeVectorFull] = GreensMatrixBuilderTemp(greensMatFile)

today = datestr(now, 'mm-dd');
rawDispData = csvread(greensMatFile,9,0);
nodeLocation = rawDispData(:,1:3);
nodeVectorFull = rawDispData(:,4:end);

end


