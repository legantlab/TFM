%Script to plot node locations of beadElems. 
testBeadElems = beadElems;

testNodeCoords = zeros(length(testBeadElems),12);

for i = 1:length(testNodeCoords)
    cur3DElem = testBeadElems(i);
    testNodes = newElem3DConn(cur3DElem, 2:5);

    %each node has 3 coords which we can get from newNodeCoords
    testNodeCoords(i,1:3) = newNodeCoords(testNodes(1),2:4);
    testNodeCoords(i,4:6) = newNodeCoords(testNodes(2),2:4);
    testNodeCoords(i,7:9) = newNodeCoords(testNodes(3),2:4);
    testNodeCoords(i,10:12) = newNodeCoords(testNodes(4),2:4);
end


 scatter3(testNodeCoords(:,1),testNodeCoords(:,2),testNodeCoords(:,3))
 hold on 
 scatter3(testNodeCoords(:,4),testNodeCoords(:,5),testNodeCoords(:,6))
 scatter3(testNodeCoords(:,7),testNodeCoords(:,8),testNodeCoords(:,9))
 scatter3(testNodeCoords(:,10),testNodeCoords(:,11),testNodeCoords(:,12))


 %% Sort beadElems by surface faces
cleanedBeadElems = [];
for j = 1:length(testBeadElems)
    curElem = testBeadElems(j);
    curElemNodes = newElem3DConn(curElem,2:5);
    if sum(ismember(curElemNodes,appliedNodeList(:,1))) == 4
        cleanedBeadElems = [cleanedBeadElems, [curElem; j]];
        %Also replacecleanedDisplacements
        
    end

end

cleanedBeadElemsIds = cleanedBeadElems';
newCleanBeadDisps = cell(length(beadDisps),1);
%Also reclean the beadDisps
for k = 1:length(beadDispsCleaned)
    curDisps = beadDispsCleaned{1};
    curDisps = curDisps(cleanedBeadElemsIds(:,2),:);
    newCleanBeadDisps{k} = curDisps;

end
