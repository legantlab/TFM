function [] = firstTimeSetup()
%firstTimeSetup Creates .mat file with variables if not present
%   Creates .mat file containint preset user variables if not present in
%   initial folder. 

%Data variables
    isCZI = 0;
    postStem = 'test';
    postTreatmentFile = 'test'; 
    postTreatmentChannel = NaN;
    preStem = 'test';
    preTreatmentFile = {'test'};
    preTreatmentFileName = 'test'; 
    preTreatmentChannel = {NaN};
    analysisStem = 'test';

    dataFileName = 'test';
    
    %Voxel Sizes
    size_x = 1;
    size_y = 1;
    size_z = 1;
    
    %Bead parameters
    beadthres = 0.5;
    beadmin = 3;
    beadmax = 1000; 
    beadwin = 5;
    beadint = 4500;
    beadmult = 1.5;
   
    
    %TPT parameters
    tptknnFD = 20;
    tptknnFM = 6;
    tptfmThres = 2;
    tptoutThres = 6;

    
    save('userInputs.mat','isCZI','postStem','postTreatmentFile','postTreatmentChannel', ...
    'preStem','preTreatmentFile','preTreatmentFileName','preTreatmentChannel','analysisStem', ...
    'dataFileName', 'size_x', 'size_y', 'size_z','beadthres','beadmin','beadmax', ...
    'beadwin', 'beadmult', 'tptknnFD','tptknnFM','tptfmThres','tptoutThres','beadint', ...
    'tptoutThres');
end

