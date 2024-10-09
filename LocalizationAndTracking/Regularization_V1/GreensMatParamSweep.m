%% Green's Function Generation Script
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

The output of this script is the variable G, an n x m matrix operator, and
displacedPtsMat, an n x 4 matrix, where columns 1, 2, and 3 give the x, y,
and z components of the displacement point, and column 4 is an indicator
for the direction the point was displaced (0 = x, 1 = y, 2 = z).

Date Modified:  11/07/2018
Modification:   Cleaned and commented
%}

function [G, displacedPtsMat] = GreensMatParamSweep(greensMatDataFile)
today = datestr(now, 'mm-dd');

%open file for reading
dispf = fopen(greensMatDataFile,'r');

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
rawDispData = csvread(greensMatDataFile,5,2);
displacedPtsNum = size(rawDispData,2)/3;

%Check for NaNs in the COMSOL data. This can happen if the cutpoints lie
%outside the geometry or if a part of the study failed to run. In the
%former case, this may be intentional, like a cut point in a channel.
%Clicking "OK" in the dialog box will remove any columns with NaN in them.
% % if sum(sum(isnan(rawDispData)))
% %     f = uifigure;
% %     s = uiconfirm(f, ['warning: NaNs found in COMSOL data. Click "OK" ' ...
% %         'if you were expecting this (i.e. some cut points in COMSOL are in a channel)'],...
% %         'NaN Warning', 'Icon', 'warning');
% %     if strcmp(s, 'Cancel')
% %         error('NaNs found in COMSOL data')
% %     end
% %     % Delete f when done so we don't have an annoying box. 
% %     f.delete
% % end


%modified by Yu to take away grids that are not in the materials (e.g.
%grids in channels) 

%Is this still a problem, not sure I understand what
%this is doing really? Grids in the channel? Yu suggests that having the
%nan value freaks out Matlab and getting rid of them doesn't lose you
%information as the spatial information is contained in the mesh centroids
%and displacement coordinates. for now go ahead and leave them in. 

%This is a logical that checks what columns have nan values
raw_nan = ~any(~isnan(rawDispData), 1);
%this sets up an empty vector to? this actually removes the nans
rawDispData(:,raw_nan)=[];
displacedPtsNum = size(rawDispData,2)/3;
%And this gets rid of those nan values...might not be what we want? Not
%sure. 
dispLocations(raw_nan,:)=[];

%What is this actually doing? 
%rawDispData = rawDispData(:,1-sum(isnan(rawDispData)));

%% Organize G How We Want
%G should be grouped by measurement locations, and then each triple gives
%   the x, y, then z components of the measurement
order = repmat(1:displacedPtsNum,3,1);
order = order + repmat([0; displacedPtsNum; 2*displacedPtsNum],1,displacedPtsNum);
order = order(:);

%These are the final outputs for the function, not sure if leaving out the
%nan values is messing us up?
G = rawDispData(:,order)';
displacedPtsMat = dispLocations(order,:);


end
