%% Arbitrary Force Application Output
%{
Author:         Keshav Patel
PI:             Wesley Legant
Date Created:   10/10/2018
Description:    This file takes mesh data from COMSOL and returns the
locations, indices, and surfaces that define the trianglar mesh.

meshPointsCoord = the coordinates of all the vertices used to make meshes.
    The rows define the index of the vertex
meshPointsIndex = the vertex indices used to make a triangular mesh. The
    rows define the index of the mesh element
meshSurfaceIndex = the index of the surface the mesh lives on. The rows
    define the index of the mesh element

Date Modified:  11/07/2018
Modification:   Cleaned and commented
%}

function [meshPointsCoords, meshPointsIndex, meshSurfaceIndex, desiredBoundaryInd] = parseMeshData(filename)
%open file for reading
dispf = fopen(filename,'r');

%% Get to Coordinates of Mesh Vertices
while true
    line = fgetl(dispf); %gets the next line of the csv
    %The data after "# Mesh point coordinates" are the coordinates
    %if strcmp(line,'# Mesh point coordinates')     modified by Yu, using
    %current mesh file setting
    if strcmp(line,'# Mesh vertex coordinates')        
        line = fgetl(dispf); 
        break
    end
end

%% Save Vertex Coordinates
meshPointsCoords = [];

% data tables are separated by an empty line
while size(line,1)~=0
    %data is organized as a 3-tuple of coordinates
    lineDoubles = sscanf(line,'%f');
    meshPointsCoords = [meshPointsCoords; lineDoubles'];
    line = fgetl(dispf); %gets the next line of the csv
end

%% Get to Vertex Indices that Define Mesh Elements
while true
    line = fgetl(dispf); %gets the next line of the csv
    lineDoubles = sscanf(line,'%f');
    %This is the only place where data is organized in triples (after
    %   coordinates
    if length(lineDoubles) == 3
        break
    end
end

%% Save Vertex Indices
meshPointsIndex = [];

% data tables are separated by an empty line
while size(line,1)~=0
    size(line);
    lineDoubles = sscanf(line,'%f');
    meshPointsIndex = [meshPointsIndex; lineDoubles'];
    line = fgetl(dispf); %gets the next line of the csv
end

%% Get to Index of Surfaces Each Mesh Lives On
while true
    line = fgetl(dispf); %gets the next line of the csv
    lineDoubles = sscanf(line,'%f');
    %The data after "# Geometric entity indices" are the indices
    if strcmp(line,'# Geometric entity indices')
        line = fgetl(dispf);
        break
    end
end

%% Save Surface Indices
meshSurfaceIndex = [];
while size(line,1)~=0
    size(line);
    lineDoubles = sscanf(line,'%f');
    meshSurfaceIndex = [meshSurfaceIndex; lineDoubles];
    line = fgetl(dispf); %gets the next line of the csv
end

%% Get to Explicit Selection Indices
while true
    line = fgetl(dispf); %gets the next line of the csv
    lineDoubles = sscanf(line,'%f');
    %The data after "# Entities" are the indices
    if strcmp(line,'# Entities')
        line = fgetl(dispf);
        break
    end
end

%% Save Explicit Selection Indices
desiredBoundaryInd = [];
while size(line,1)~=0 & line ~= -1
    lineDoubles = sscanf(line,'%f');
    desiredBoundaryInd = [desiredBoundaryInd lineDoubles];
    line = fgetl(dispf); %gets the next line of the csv
end
fclose(dispf);

end
