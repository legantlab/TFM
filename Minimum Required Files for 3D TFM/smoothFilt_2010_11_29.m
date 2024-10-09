function [matchedCoordinates_filt]=smoothFilt_2010_11_29(matchedCoordinates,numNeighbors,tol,numProcessors,jm)
% matchedCoordinates={zeros(200,4),zeros(200,4)};
%This function finds matches beads between reference and deformed
%conditions and outputs the vectors cooresponding to each bead displacement

%Input argument definitions
% reference = a cell array containing the centroids for the beads in the
% reference configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

% deformed = a cell array containing the centroids for the beads in the deformed
% configuration of the form {bead index, xcoord, ycoord, zcoord}.  Each row
% corresponds to one bead.

%Note, the length of reference and deformed need not be the same

% numNeighbors=50; %The number of neighbors to scan for initial feature vector generation
% numVectorsRef=3; %The number of feature vectors to match between reference and deformed configurations
% numVectorsDef=5; %The number of feature vectors to search in the deformed configuration

reference=mat2cell(cat(2,[1:1:length(matchedCoordinates)]',matchedCoordinates(:,2:4)),ones(length(matchedCoordinates),1),4);
neighborRef_Ref = nearestNeighborNew(reference,reference,numNeighbors,1);

if numProcessors==2
splitIndexR=floor(length(reference)/2);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});

submit(job1)
submit(job2)

waitForState(job1,'finished')
results1 = getAllOutputArguments(job1);
waitForState(job2,'finished')
results2 = getAllOutputArguments(job2);
neighborRef_Ref=cat(1,results1{1},results2{1});
destroy(jobs)
end


if numProcessors==3
splitIndexR=floor(length(reference)/3);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});

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
destroy(jobs)
end

if numProcessors==4
splitIndexR=floor(length(reference)/4);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numNeighbors,1});

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
destroy(jobs)
end

if numProcessors==5
splitIndexR=floor(length(reference)/5);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numNeighbors,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numNeighbors,1});

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
destroy(jobs)
end

if numProcessors==6
splitIndexR=floor(length(reference)/6);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numNeighbors,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numNeighbors,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numNeighbors,1});

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
destroy(jobs)
end

if numProcessors==7
splitIndexR=floor(length(reference)/7);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:6*splitIndexR);
reference_7=reference(6*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numNeighbors,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numNeighbors,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numNeighbors,1});
task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,reference,numNeighbors,1});

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
destroy(jobs)
end

if numProcessors==8
splitIndexR=floor(length(reference)/8);
reference_1=reference(1:splitIndexR);
reference_2=reference(splitIndexR+1:2*splitIndexR);
reference_3=reference(2*splitIndexR+1:3*splitIndexR);
reference_4=reference(3*splitIndexR+1:4*splitIndexR);
reference_5=reference(4*splitIndexR+1:5*splitIndexR);
reference_6=reference(5*splitIndexR+1:6*splitIndexR);
reference_7=reference(6*splitIndexR+1:7*splitIndexR);
reference_8=reference(7*splitIndexR+1:end);

job1 = createJob(jm, 'Configuration', 'WRL Config');
job2 = createJob(jm, 'Configuration', 'WRL Config');
job3 = createJob(jm, 'Configuration', 'WRL Config');
job4 = createJob(jm, 'Configuration', 'WRL Config');
job5 = createJob(jm, 'Configuration', 'WRL Config');
job6 = createJob(jm, 'Configuration', 'WRL Config');
job7 = createJob(jm, 'Configuration', 'WRL Config');
job8 = createJob(jm, 'Configuration', 'WRL Config');
jobs=get(jm,'jobs');
% set(jobs,'PathDependencies',{'O:\WRL Temp\Bead displacement scripts'});

task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,reference,numNeighbors,1});
task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,reference,numNeighbors,1});
task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,reference,numNeighbors,1});
task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,reference,numNeighbors,1});
task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,reference,numNeighbors,1});
task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,reference,numNeighbors,1});
task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,reference,numNeighbors,1});
task8_1=createTask(job8, @nearestNeighborNew, 1, {reference_8,reference,numNeighbors,1});

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
destroy(jobs)
end
% 
% 
% tic
% 'mapping reference configuration'
% neighborRef_Ref=nearestNeighborNew(reference,reference,numNeighbors,1); %A mapping of each bead and its neighbors in the reference configuration
% toc
% 'mapping deformed configuration'
% tic
% neighborDef_Def=nearestNeighborNew(deformed,deformed,numVectorsDef,1); %A mapping of each bead and its neighbors in the deformed configuration
% toc
% 'mapping reference to deformed configuration'
% tic
% neighborRef_Def=nearestNeighborNew(reference,deformed,numNeighbors,0); %A mapping of each bead in the reference configuration and its neighbors in the deformed configuration
% toc
t=1; %A counter to keep track of how many beads are matched

%For each bead in the reference configuration, re-sort the feature vectors of each of the
%neighbors in the deformed configuration to give the best match

t=1;
for i=1:length(reference)
    refbead=reference{i}; %index and centroid of the reference bead in the material
    refVector=[matchedCoordinates(refbead(1,1),6:8)-matchedCoordinates(refbead(1,1),2:4)];
    for j=2:numNeighbors+1
        neighbor_bead=neighborRef_Ref{i,j}; %index and centroid of the jth nearest neighbor to the reference bead in the deformed configuration
        neighborVector(j-1,:)=[matchedCoordinates(neighbor_bead(1,1),6:8)-matchedCoordinates(neighbor_bead(1,1),2:4)];
    end
    dist=sqrt(sum((bsxfun(@minus,neighborVector,refVector)).^2,2));
%     wAvg=sum(bsxfun(@times,1./dist,neighborVector))./sum(1./dist);
    if mean(neighborVector(:,1)-refVector(1))<(tol*mean(std(neighborVector(:,1))))&&mean(neighborVector(:,2)-refVector(2))<(tol*mean(std(neighborVector(:,2))))&&mean(neighborVector(:,3)-refVector(3))<(tol*mean(std(neighborVector(:,3))))
%     if mean(abs(mean(neighborVector)-refVector))<(tol*mean(std(neighborVector)))
        matchedCoordinates_filt(t,:)=matchedCoordinates(refbead(1,1),:);
        t=t+1;
    end
end



