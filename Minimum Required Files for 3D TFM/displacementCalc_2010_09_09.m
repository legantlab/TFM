function [matchedCoordinates]=displacementCalc_temp_2010_09_09(reference,deformed,numNeighbors,numVectorsRef,numVectorsDef,tol,numProcessors,jm)
%This function matches beads between reference and deformed
%conditions and outputs the matched centroids.

%Input argument definitions
% reference = a matrix containing the centroids for the beads in the
% reference configuration of the form [bead index, xcoord, ycoord, zcoord].  Each row
% corresponds to one bead.  "bead index" is a numerical ordering to keep
% track of each bead (based on the order in which the centroids were
% found).

% deformed = a matrix containing the centroids for the beads in the deformed
% configuration of the form [bead index, xcoord, ycoord, zcoord].  Each row
% corresponds to one bead. Note, the length of reference and deformed need
% not be the same.

%numNeighbors = the number of neighbors to scan as a potential match.
%Usually between 20 and 100.  Should be larger if drift and swelling are
%present between the datasets.

%numVectorsRef = the number of feature vectors to generate for the
%reference dataset.  A good estimate is 2 or 3.

%numVectorsDef = the number of feature vectors to generate for the deformed
%dataset.  Should be 1 or 2 larger than numVectorsRef.

%tol = a cutoff value indicating how much more closely the feature vectors of the ideal match need to
%be when compaired to all other potential matches between the reference and deformed configurations in order to
%consider it a correct match. Typically a cutoff between 1.25 and 2 is
%appropriate.

%numProcessors = number of processors to use

%jm = the jobmanager for parallel processing

%Output arguement definitions
%matchedCoordinates = an array containing the matched beads in the reference and deformed configurations
%where each column is [bead_indexref, xcoordref ycoordref, zcoordref, bead_indexdef,
%xcoorddef ycoorddef, zcoorddef]


%Initialize variables and convert to cell arrays for nearest neighbor
%matching
matchedCoordinates={};
reference=mat2cell(cat(2,[1:1:length(reference)]',reference(:,2:4)),ones(length(reference),1),4);
deformed=mat2cell(cat(2,[1:1:length(deformed)]',deformed(:,2:4)),ones(length(deformed),1),4);

% Ridiculous conditional to split up for parallel processing
if numProcessors==1
neighborRef_Ref=nearestNeighborNew(reference,reference,numVectorsRef,1);
neighborDef_Def=nearestNeighborNew(deformed,deformed,numVectorsDef,1);
neighborRef_Def=nearestNeighborNew(reference,deformed,numNeighbors,1);
end

if numProcessors==2
splitIndexR=floor(length(reference)/2);
splitIndexD=floor(length(deformed)/2);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});

submit(job1)
submit(job2)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
neighborRef_Ref=cat(1,results1{1},results2{1});
neighborDef_Def=cat(1,results1{2},results2{2});
neighborRef_Def=cat(1,results1{3},results2{3});
destroy(jobs)
end

if numProcessors==3
splitIndexR=floor(length(reference)/3);
splitIndexD=floor(length(deformed)/3);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3});
destroy(jobs)
end

if numProcessors==4
splitIndexR=floor(length(reference)/4);
splitIndexD=floor(length(deformed)/4);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:3*splitIndexD);
deformed_4=deformed(3*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});
task4_2=createTask(job4, @nearestNeighborNew, 1, {deformed_4,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});
task4_3=createTask(job4, @nearestNeighborNew, 1, {reference_4,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)
submit(job4)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1},results4{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2},results4{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3},results4{3});
destroy(jobs)
end

if numProcessors==5
splitIndexR=floor(length(reference)/5);
splitIndexD=floor(length(deformed)/5);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:3*splitIndexD);
deformed_4=deformed(3*splitIndexD+1:4*splitIndexD);
deformed_5=deformed(4*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numVectorsRef,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});
task4_2=createTask(job4, @nearestNeighborNew, 1, {deformed_4,deformed,numVectorsDef,1});
task5_2=createTask(job5, @nearestNeighborNew, 1, {deformed_5,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});
task4_3=createTask(job4, @nearestNeighborNew, 1, {reference_4,deformed,numNeighbors,0});
task5_3=createTask(job5, @nearestNeighborNew, 1, {reference_5,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2},results4{2},results5{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3},results4{3},results5{3});
destroy(jobs)
end

if numProcessors==6
splitIndexR=floor(length(reference)/6);
splitIndexD=floor(length(deformed)/6);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:3*splitIndexD);
deformed_4=deformed(3*splitIndexD+1:4*splitIndexD);
deformed_5=deformed(4*splitIndexD+1:5*splitIndexD);
deformed_6=deformed(5*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numVectorsRef,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numVectorsRef,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});
task4_2=createTask(job4, @nearestNeighborNew, 1, {deformed_4,deformed,numVectorsDef,1});
task5_2=createTask(job5, @nearestNeighborNew, 1, {deformed_5,deformed,numVectorsDef,1});
task6_2=createTask(job6, @nearestNeighborNew, 1, {deformed_6,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});
task4_3=createTask(job4, @nearestNeighborNew, 1, {reference_4,deformed,numNeighbors,0});
task5_3=createTask(job5, @nearestNeighborNew, 1, {reference_5,deformed,numNeighbors,0});
task6_3=createTask(job6, @nearestNeighborNew, 1, {reference_6,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2},results4{2},results5{2},results6{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3},results4{3},results5{3},results6{3});
destroy(jobs)
end

if numProcessors==7
splitIndexR=floor(length(reference)/7);
splitIndexD=floor(length(deformed)/7);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:6*splitIndexR);
reference_7=reference(6*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:3*splitIndexD);
deformed_4=deformed(3*splitIndexD+1:4*splitIndexD);
deformed_5=deformed(4*splitIndexD+1:5*splitIndexD);
deformed_6=deformed(5*splitIndexD+1:6*splitIndexD);
deformed_7=deformed(6*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numVectorsRef,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numVectorsRef,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numVectorsRef,1});
task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});
task4_2=createTask(job4, @nearestNeighborNew, 1, {deformed_4,deformed,numVectorsDef,1});
task5_2=createTask(job5, @nearestNeighborNew, 1, {deformed_5,deformed,numVectorsDef,1});
task6_2=createTask(job6, @nearestNeighborNew, 1, {deformed_6,deformed,numVectorsDef,1});
task7_2=createTask(job7, @nearestNeighborNew, 1, {deformed_7,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});
task4_3=createTask(job4, @nearestNeighborNew, 1, {reference_4,deformed,numNeighbors,0});
task5_3=createTask(job5, @nearestNeighborNew, 1, {reference_5,deformed,numNeighbors,0});
task6_3=createTask(job6, @nearestNeighborNew, 1, {reference_6,deformed,numNeighbors,0});
task7_3=createTask(job7, @nearestNeighborNew, 1, {reference_7,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)
submit(job7)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);
waitForState(job7,'finished')
results7 = getAllOutputArguments(job7);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2},results4{2},results5{2},results6{2},results7{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3},results4{3},results5{3},results6{3},results7{3});
destroy(jobs)
end

if numProcessors==8
splitIndexR=floor(length(reference)/8);
splitIndexD=floor(length(deformed)/8);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:6*splitIndexR);
reference_7=reference(6*splitIndexR+1:7*splitIndexR);
reference_8=reference(7*splitIndexR+1:end);

deformed_1=deformed(1:splitIndexD);
deformed_2=deformed(splitIndexD+1:2*splitIndexD);
deformed_3=deformed(2*splitIndexD+1:3*splitIndexD);
deformed_4=deformed(3*splitIndexD+1:4*splitIndexD);
deformed_5=deformed(4*splitIndexD+1:5*splitIndexD);
deformed_6=deformed(5*splitIndexD+1:6*splitIndexD);
deformed_7=deformed(6*splitIndexD+1:7*splitIndexD);
deformed_8=deformed(7*splitIndexD+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
job8 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numVectorsRef,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numVectorsRef,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numVectorsRef,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numVectorsRef,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numVectorsRef,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numVectorsRef,1});
task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,reference,numVectorsRef,1});
task8_1=createTask(job8, @nearestNeighborNew, 1, {reference_8,reference,numVectorsRef,1});

task1_2=createTask(job1, @nearestNeighborNew, 1, {deformed_1,deformed,numVectorsDef,1});
task2_2=createTask(job2, @nearestNeighborNew, 1, {deformed_2,deformed,numVectorsDef,1});
task3_2=createTask(job3, @nearestNeighborNew, 1, {deformed_3,deformed,numVectorsDef,1});
task4_2=createTask(job4, @nearestNeighborNew, 1, {deformed_4,deformed,numVectorsDef,1});
task5_2=createTask(job5, @nearestNeighborNew, 1, {deformed_5,deformed,numVectorsDef,1});
task6_2=createTask(job6, @nearestNeighborNew, 1, {deformed_6,deformed,numVectorsDef,1});
task7_2=createTask(job7, @nearestNeighborNew, 1, {deformed_7,deformed,numVectorsDef,1});
task8_2=createTask(job8, @nearestNeighborNew, 1, {deformed_8,deformed,numVectorsDef,1});

task1_3=createTask(job1, @nearestNeighborNew, 1, {reference_1,deformed,numNeighbors,0});
task2_3=createTask(job2, @nearestNeighborNew, 1, {reference_2,deformed,numNeighbors,0});
task3_3=createTask(job3, @nearestNeighborNew, 1, {reference_3,deformed,numNeighbors,0});
task4_3=createTask(job4, @nearestNeighborNew, 1, {reference_4,deformed,numNeighbors,0});
task5_3=createTask(job5, @nearestNeighborNew, 1, {reference_5,deformed,numNeighbors,0});
task6_3=createTask(job6, @nearestNeighborNew, 1, {reference_6,deformed,numNeighbors,0});
task7_3=createTask(job7, @nearestNeighborNew, 1, {reference_7,deformed,numNeighbors,0});
task8_3=createTask(job8, @nearestNeighborNew, 1, {reference_8,deformed,numNeighbors,0});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)
submit(job7)
submit(job8)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);
waitForState(job7,'finished')
results7 = getAllOutputArguments(job7);
waitForState(job8,'finished')
results8 = getAllOutputArguments(job8);
neighborRef_Ref=cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1},results8{1});
neighborDef_Def=cat(1,results1{2},results2{2},results3{2},results4{2},results5{2},results6{2},results7{2},results8{2});
neighborRef_Def=cat(1,results1{3},results2{3},results3{3},results4{3},results5{3},results6{3},results7{3},results8{3});
destroy(jobs)
end

% t=1; %A counter to keep track of how many beads are matched

%For each bead in the reference configuration, re-sort the feature vectors of each of the
%neighbors in the deformed configuration to give the best match

%Convert cell arrays back into matrices
reference=cell2mat(reference);
deformed=cell2mat(deformed);
neighborRef_Ref=cell2mat(neighborRef_Ref);
neighborDef_Def=cell2mat(neighborDef_Def);
neighborRef_Def=cell2mat(neighborRef_Def);

%Initialize some variables
refVectors=zeros(length(reference),3*numVectorsRef); %A n by 3*nVR array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, ...,u3NVR_x, u3NVR_y, u3NVR_z].  Each row is a bead (up to n reference beads) and u's are Cartesian components of each feature vector.
matchIndex=zeros(length(reference),1); %A vector to store the indicies of the matched beads

%Generate feature vectors for each bead of the reference configuration
t=1;
for i=4:4:4*numVectorsRef
    refVectors(:,t:t+2)=[neighborRef_Ref(:,i+2:i+4)-neighborRef_Ref(:,2:4)]; 
    t=t+3;
end

%Generate feature vectors for each potential match in the deformed configuration
defVectors=zeros(length(reference),4*numNeighbors);%A n by 3*nVD*numNeighbors array of the form [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, ...,u3NVD*numNeighbors_x, u3NVD*numNeighbors_y, u3NVD*numNeighbors_z].
                                                   %Each row is vectors for the potential matches to a given bead in the reference configuration.
t=1;
for k=4:4:4*numNeighbors+1;
for i=4:4:4*numVectorsDef
    defVectors(:,t:t+2)=[neighborDef_Def(neighborRef_Def(:,k+1),i+2:i+4)-neighborDef_Def(neighborRef_Def(:,k+1),2:4)];
    t=t+3;
end
end

% Ridiculous conditional to split up for parallel processing
if numProcessors==1
index=matchFeatureVector_2010_09_09(refVectors,defVectors,numVectorsRef,numVectorsDef,numNeighbors,tol);
end

if numProcessors==2
splitIndexR=floor(length(refVectors)/2);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);

index=[results1{1};results2{1}];
destroy(jobs)
end

if numProcessors==3
splitIndexR=floor(length(refVectors)/3);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);

index=[results1{1};results2{1};results3{1}];
destroy(jobs)
end

if numProcessors==4
splitIndexR=floor(length(refVectors)/4);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:3*splitIndexR,:);
refVectors_4=refVectors(3*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:3*splitIndexR,:);
defVectors_4=defVectors(3*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});
task4_1=createTask(job4, @matchFeatureVector_2010_09_09, 1, {refVectors_4,defVectors_4,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)
submit(job4)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);

index=[results1{1};results2{1};results3{1};results4{1}];
destroy(jobs)
end

if numProcessors==5
splitIndexR=floor(length(refVectors)/5);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:3*splitIndexR,:);
refVectors_4=refVectors(3*splitIndexR+1:4*splitIndexR,:);
refVectors_5=refVectors(4*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:3*splitIndexR,:);
defVectors_4=defVectors(3*splitIndexR+1:4*splitIndexR,:);
defVectors_5=defVectors(4*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});
task4_1=createTask(job4, @matchFeatureVector_2010_09_09, 1, {refVectors_4,defVectors_4,numVectorsRef,numVectorsDef,numNeighbors,tol});
task5_1=createTask(job5, @matchFeatureVector_2010_09_09, 1, {refVectors_5,defVectors_5,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);


index=[results1{1};results2{1};results3{1};results4{1};results5{1}];
destroy(jobs)
end


if numProcessors==6
splitIndexR=floor(length(refVectors)/6);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:3*splitIndexR,:);
refVectors_4=refVectors(3*splitIndexR+1:4*splitIndexR,:);
refVectors_5=refVectors(4*splitIndexR+1:5*splitIndexR,:);
refVectors_6=refVectors(5*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:3*splitIndexR,:);
defVectors_4=defVectors(3*splitIndexR+1:4*splitIndexR,:);
defVectors_5=defVectors(4*splitIndexR+1:5*splitIndexR,:);
defVectors_6=defVectors(5*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});
task4_1=createTask(job4, @matchFeatureVector_2010_09_09, 1, {refVectors_4,defVectors_4,numVectorsRef,numVectorsDef,numNeighbors,tol});
task5_1=createTask(job5, @matchFeatureVector_2010_09_09, 1, {refVectors_5,defVectors_5,numVectorsRef,numVectorsDef,numNeighbors,tol});
task6_1=createTask(job6, @matchFeatureVector_2010_09_09, 1, {refVectors_6,defVectors_6,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);


index=[results1{1};results2{1};results3{1};results4{1};results5{1};results6{1}];
destroy(jobs)
end

if numProcessors==7
splitIndexR=floor(length(refVectors)/7);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:3*splitIndexR,:);
refVectors_4=refVectors(3*splitIndexR+1:4*splitIndexR,:);
refVectors_5=refVectors(4*splitIndexR+1:5*splitIndexR,:);
refVectors_6=refVectors(5*splitIndexR+1:6*splitIndexR,:);
refVectors_7=refVectors(6*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:3*splitIndexR,:);
defVectors_4=defVectors(3*splitIndexR+1:4*splitIndexR,:);
defVectors_5=defVectors(4*splitIndexR+1:5*splitIndexR,:);
defVectors_6=defVectors(5*splitIndexR+1:6*splitIndexR,:);
defVectors_7=defVectors(6*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});
task4_1=createTask(job4, @matchFeatureVector_2010_09_09, 1, {refVectors_4,defVectors_4,numVectorsRef,numVectorsDef,numNeighbors,tol});
task5_1=createTask(job5, @matchFeatureVector_2010_09_09, 1, {refVectors_5,defVectors_5,numVectorsRef,numVectorsDef,numNeighbors,tol});
task6_1=createTask(job6, @matchFeatureVector_2010_09_09, 1, {refVectors_6,defVectors_6,numVectorsRef,numVectorsDef,numNeighbors,tol});
task7_1=createTask(job7, @matchFeatureVector_2010_09_09, 1, {refVectors_7,defVectors_7,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)
submit(job7)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);
waitForState(job7,'finished')
results7 = getAllOutputArguments(job7);

index=[results1{1};results2{1};results3{1};results4{1};results5{1};results6{1};results7{1}];
destroy(jobs)
end

if numProcessors==8
splitIndexR=floor(length(refVectors)/8);
refVectors_1=refVectors(1:splitIndexR,:);
refVectors_2=refVectors(splitIndexR+1:2*splitIndexR,:);
refVectors_3=refVectors(2*splitIndexR+1:3*splitIndexR,:);
refVectors_4=refVectors(3*splitIndexR+1:4*splitIndexR,:);
refVectors_5=refVectors(4*splitIndexR+1:5*splitIndexR,:);
refVectors_6=refVectors(5*splitIndexR+1:6*splitIndexR,:);
refVectors_7=refVectors(6*splitIndexR+1:7*splitIndexR,:);
refVectors_8=refVectors(7*splitIndexR+1:end,:);

defVectors_1=defVectors(1:splitIndexR,:);
defVectors_2=defVectors(splitIndexR+1:2*splitIndexR,:);
defVectors_3=defVectors(2*splitIndexR+1:3*splitIndexR,:);
defVectors_4=defVectors(3*splitIndexR+1:4*splitIndexR,:);
defVectors_5=defVectors(4*splitIndexR+1:5*splitIndexR,:);
defVectors_6=defVectors(5*splitIndexR+1:6*splitIndexR,:);
defVectors_7=defVectors(6*splitIndexR+1:7*splitIndexR,:);
defVectors_8=defVectors(7*splitIndexR+1:end,:);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
job8 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');

task1_1=createTask(job1, @matchFeatureVector_2010_09_09, 1, {refVectors_1,defVectors_1,numVectorsRef,numVectorsDef,numNeighbors,tol});
task2_1=createTask(job2, @matchFeatureVector_2010_09_09, 1, {refVectors_2,defVectors_2,numVectorsRef,numVectorsDef,numNeighbors,tol});
task3_1=createTask(job3, @matchFeatureVector_2010_09_09, 1, {refVectors_3,defVectors_3,numVectorsRef,numVectorsDef,numNeighbors,tol});
task4_1=createTask(job4, @matchFeatureVector_2010_09_09, 1, {refVectors_4,defVectors_4,numVectorsRef,numVectorsDef,numNeighbors,tol});
task5_1=createTask(job5, @matchFeatureVector_2010_09_09, 1, {refVectors_5,defVectors_5,numVectorsRef,numVectorsDef,numNeighbors,tol});
task6_1=createTask(job6, @matchFeatureVector_2010_09_09, 1, {refVectors_6,defVectors_6,numVectorsRef,numVectorsDef,numNeighbors,tol});
task7_1=createTask(job7, @matchFeatureVector_2010_09_09, 1, {refVectors_7,defVectors_7,numVectorsRef,numVectorsDef,numNeighbors,tol});
task8_1=createTask(job8, @matchFeatureVector_2010_09_09, 1, {refVectors_8,defVectors_8,numVectorsRef,numVectorsDef,numNeighbors,tol});

submit(job1)
submit(job2)
submit(job3)
submit(job4)
submit(job5)
submit(job6)
submit(job7)
submit(job8)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
waitForState(job3,'finished')
results3 = getAllOutputArguments(job3);
waitForState(job4,'finished')
results4 = getAllOutputArguments(job4);
waitForState(job5,'finished')
results5 = getAllOutputArguments(job5);
waitForState(job6,'finished')
results6 = getAllOutputArguments(job6);
waitForState(job7,'finished')
results7 = getAllOutputArguments(job7);
waitForState(job8,'finished')
results8 = getAllOutputArguments(job8);

index=[results1{1};results2{1};results3{1};results4{1};results5{1};results6{1};results7{1};results8{1}];
destroy(jobs)
end

%Some comment here!
matchIndex=zeros(length(refVectors),1);
for i=1:length(matchIndex)
    if index(i)>0
        matchIndex(i)=neighborRef_Def(i,index(i));
    end
end
matchedCoordinates=[reference(matchIndex>0,1:4),deformed(matchIndex(matchIndex>0),1:4)];
end




