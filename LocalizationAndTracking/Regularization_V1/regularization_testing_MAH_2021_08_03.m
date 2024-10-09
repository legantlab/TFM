%script to test regularization functions and step through code. 

clc; close all; clear;
today = datestr(now, 'mm-dd');


%We can step back, we don't need the GUIs at the moment so lets focus on
%just loading in things and see that they make sense/timing so that we can
%speed up the code.

%% Define conditions used to create inverse problem

%greensMatLoadFile = inputs.greensLoad;%'Z:\Regan\2019_10_10(TFM)\data_FLAT_20kPa_Analysis\GreensMat_block_10_22_19.mat';
%greensMatDataFile = inputs.nodalSoln;%'Z:\Regan\2019_10_10(TFM)\data_FLAT_20kPa_Analysis\nodalSolutions_10_22_19.csv';
%greensMatSaveAsFile = inputs.greensSave;%'Z:\Regan\2019_10_10(TFM)\data_FLAT_20kPa_Analysis\GreensMat_block_10_22_19.mat';
COMSOLGeometryMeshFile = 'X:\Max\2021_08_03_MatlabTestFiles\meshDatameshData_8_2_21.mphtxt'; %by default, COMSOL saves this to '../data/meshData.mphtxt'
%particleTrackingDataFile = inputs.loadPT;%Z:\Regan\2019_10_10(TFM)\data_FLAT_20kPa_Analysis\cell3\cell3_tracking.mat';



%Regularization method
method = 'tikh';

%% Setup

disp('Reading COMSOL mesh data')
% Find centroids of where forces are applied (for plotting purposes)
[meshPointsCoords, meshPointsIndex, meshSurfaceIndex, desiredBoundaryInd] = parseMeshData(COMSOLGeometryMeshFile);
    

% Lets see if we can visualize the geometry
%scatter3(meshPointsCoords(:,1), meshPointsCoords(:,2), meshPointsCoords(:,3)) 

%the mesh centroids will be used to visually approximate where forces were
%   applied
meshCentroid = [];

