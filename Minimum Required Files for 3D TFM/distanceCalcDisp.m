function [dist]=distanceCalcDisp(referenceDisp,elemCents,numProcessors,jm)
%This function computes the distance between the bead coordinates
%specificed in referenceDisp and the surface of the cell specified by
%elemCents

%Input argument definitions
% referenceDisp = a NX4 matrix containing the centroids for the beads in the
% following format [index,xcoord,ycoord,zcoord] where N is the number of
% beads

% elemCents = a Mx3 matrix containing the centroids of each facet on the cell
% surface in the following format [xcoord,ycoord,zcoord] where M is the
% number of facets on the cell surface

% numProcessors = an integer indicating the number of processors to use

% jm = a pointer for the current job manager (if using parallel
% implementation)

numNeighbors=1; %The number of nearest neighbors (should be set to one)
%Convert matrices to cell arrays for nearest neighbor matching
reference=mat2cell(referenceDisp(:,1:4),ones(length(referenceDisp),1),4);
elemCents=mat2cell([[1:1:length(elemCents)]',elemCents],ones(length(elemCents),1),4);

%Ridiculous conditional to split up for parallel processing and identify
%nearest elemCent to each bead
if numProcessors==1
    [result]=nearestNeighborNew(reference,elemCents,numNeighbors,0);
    distance=cell2mat(result);
end

if numProcessors==2
    splitIndexR=floor(length(reference)/2);
    reference_1=reference(1:splitIndexR);
    reference_2=reference(splitIndexR+1:end);
    
    job1 = createJob(jm, 'Configuration', 'WRL Config');
    job2 = createJob(jm, 'Configuration', 'WRL Config');
    jobs=get(jm,'jobs');
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    
    submit(job1)
    submit(job2)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    distance=cell2mat(cat(1,results1{1},results2{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    
    submit(job1)
    submit(job2)
    submit(job3)
    
    waitForState(job1,'finished')
    results1 = getAllOutputArguments(job1);
    waitForState(job2,'finished')
    results2 = getAllOutputArguments(job2);
    waitForState(job3,'finished')
    results3 = getAllOutputArguments(job3);
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,elemCents,numNeighbors,0});
    
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
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1},results4{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,elemCents,numNeighbors,0});
    task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,elemCents,numNeighbors,0});
    
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
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1},results4{1},results5{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,elemCents,numNeighbors,0});
    task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,elemCents,numNeighbors,0});
    task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,elemCents,numNeighbors,0});
    
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
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,elemCents,numNeighbors,0});
    task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,elemCents,numNeighbors,0});
    task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,elemCents,numNeighbors,0});
    task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,elemCents,numNeighbors,0});
    
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
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1}));
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
    
    task1_1=createTask(job1, @nearestNeighborNew, 1, {reference_1,elemCents,numNeighbors,0});
    task2_1=createTask(job2, @nearestNeighborNew, 1, {reference_2,elemCents,numNeighbors,0});
    task3_1=createTask(job3, @nearestNeighborNew, 1, {reference_3,elemCents,numNeighbors,0});
    task4_1=createTask(job4, @nearestNeighborNew, 1, {reference_4,elemCents,numNeighbors,0});
    task5_1=createTask(job5, @nearestNeighborNew, 1, {reference_5,elemCents,numNeighbors,0});
    task6_1=createTask(job6, @nearestNeighborNew, 1, {reference_6,elemCents,numNeighbors,0});
    task7_1=createTask(job7, @nearestNeighborNew, 1, {reference_7,elemCents,numNeighbors,0});
    task8_1=createTask(job8, @nearestNeighborNew, 1, {reference_8,elemCents,numNeighbors,0});
    
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
    distance=cell2mat(cat(1,results1{1},results2{1},results3{1},results4{1},results5{1},results6{1},results7{1},results8{1}));
    destroy(jobs)
end

%Compute the distance between each bead and the nearest element centroid of
%the cell (ie shortest distance between the bead and the cell surface)
dist=sqrt((distance(:,6)-distance(:,2)).^2+(distance(:,7)-distance(:,3)).^2+(distance(:,8)-distance(:,4)).^2);

end


