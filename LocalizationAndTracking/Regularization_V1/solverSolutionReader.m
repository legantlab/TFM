%% Green's Function Generation Script (Depreciated)
%{
Author:         Keshav Patel
PI:             Wesley Legant
Date Created:   07/30/2018
Description:    This file takes data from COMSOL and generates a Green's 
matrix (i.e. the "G" in G*f=u, where f is a vector of forces and u is a 
vector of displacements). Within the COMSOL script, a study is done for each 
application of a single force of magnitude 1 applied in either the x, y, or 
z direction on the surface of each mesh element. The displacement in each 
cardinal direction at points chosen in COMSOL (either random or specified)
is then saved in a .csv file

The output of this script is the variable X, a grid built by the location
of the data samples from COMSOL, and V, a cell of grid values for the grid
locations.

Date Modified:  9/8/2019
Modification:   Removed from main code
%}

function [X, V] = solverSolutionReader(greensMatFile)
today = datestr(now, 'mm-dd');

%open file for reading
dispf = fopen(greensMatFile,'r');

%% Skip First Lines
%the column headers contain the component of the displacement reported and
%   the point the displacement is measured
headers = {};

%first lines give info on data table - will skip through most of it
%%{
while true 
    line = fgetl(dispf); %gets the next line of the csv
    lineCell = strsplit(line,',');
    
    %In the current csv setup, "% ele" is the first entry in the headers
    %   row; after that it is all data
    if strcmp(lineCell{1},'% ele') || strcmp(lineCell{1},'% coord')
        headers = strsplit(line,'"');
        line = fgetl(dispf);
        break
    end
end
%}
%% Read Header for Displacement Locations 
dispLocations = zeros(length(2:2:(length(headers)-1)),4);

%headers get parsed weirdly; the information we want in in every other
%   element of headers. 
%Data is grouped such that all x components are in
%   the 1st third, y components are in the 2nd third, ...; right now we
%   only care about the locations, so we are just moving through the first third
%%{
for i = 2:2:(length(headers)-1)
    currHead = headers{i};
    
    %scans each header for the relevant information and saves it to a cell
    num = textscan(currHead, 'Displacement field, %s component (µm), Point: (%f, %f, %f)');
    
    %displacement measurements
    dispLocations(i/2,1) = num{end-2};
    dispLocations(i/2,2) = num{end-1};
    dispLocations(i/2,3) = num{end};
    
    %component of displacement vector
    if strcmp(num{1},'Y')
        dispLocations(i/2,4) = 1;
    elseif strcmp(num{1},'Z')
        dispLocations(i/2,4) = 2;
    end
    
end

%displacedPtsNum = size(dispLocations,1)/3;
%}
fclose(dispf); %we have run through all the data at this point
%% Read in Displacement Data
%now we get to the actual meat of the data; again the columns are grouped 
%   by displacement component, which we will fix later. Rows are organized 
%   the way we want (go figure). Note we are also saving the row headers here
rawDispData = csvread(greensMatFile,5,2);
displacedPtsNum = size(rawDispData,2)/3;

if sum(sum(isnan(rawDispData)))
    error('NaNs found in data')
end

%% Organize G How We Want

[X.x, X.y, X.z] = meshgrid(sort(unique(dispLocations(:,1)')), sort(unique(dispLocations(:,2)')), sort(unique(dispLocations(:,3)')));
V = cell(size(rawDispData,1),1);

for i = 1:size(rawDispData,1)
    v.x = reshape(rawDispData(i,1:displacedPtsNum),size(X.x));
    v.y = reshape(rawDispData(i,displacedPtsNum+1:2*displacedPtsNum),size(X.x));
    v.z = reshape(rawDispData(i,2*displacedPtsNum+1:3*displacedPtsNum),size(X.x));
    V{i} = v;
end


end
