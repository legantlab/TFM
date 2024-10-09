function [ meshPointMatrix, elementMatrix ] = plotComsolMesh( meshName )
i = 1;
identifier = {'# Mesh vertex coordinates', '# Elements'};
fid = fopen(strcat(meshName, '.mphtxt'));
%type (strcat(meshName, '.mphtxt'));
while ~feof(fid)
    readLine = cellstr(fgets(fid));
    if strcmp(readLine, identifier(1))
        while(~strcmp(readLine, cellstr('')))
            readLine = cellstr(fgets(fid));
            meshPoint(i, :) = readLine;
            i = i+1;
        end
        meshPoint(end, :) = [];
        for i = 1:size(meshPoint, 1)
            meshPointMatrix(i, :) = str2num(cell2mat(meshPoint(i)));
        end
        clear 'meshPoint';
        i = 1;
       
    end
    if strcmp(readLine, identifier(2))
        while(~strcmp(readLine, cellstr('')))
            readLine = cellstr(fgets(fid));
            element(i, :) = readLine;
            i = i+1;
        end
        element(end, :) = [];
        elementMatrix = str2num(cell2mat(element(1)));
        for i = 1:size(element, 1)
            x = str2num(cell2mat(element(i)));
            elementMatrix = [elementMatrix; x];
        end
        %clear 'element';
        i = 1;
    end
end
fclose(fid);
figure(1), clf, hold on, axis equal,
% for i = 1:size(elementMatrix, 1)
%     plot([meshPointMatrix(elementMatrix(i, 1)+1,1),meshPointMatrix(elementMatrix(i, 2)+1,1)],[meshPointMatrix(elementMatrix(i, 1)+1,2),meshPointMatrix(elementMatrix(i, 2)+1,2)], '-k');
%     plot([meshPointMatrix(elementMatrix(i, 2)+1,1),meshPointMatrix(elementMatrix(i, 3)+1,1)],[meshPointMatrix(elementMatrix(i, 2)+1,2),meshPointMatrix(elementMatrix(i, 3)+1,2)], '-k');
%     plot([meshPointMatrix(elementMatrix(i, 3)+1,1),meshPointMatrix(elementMatrix(i, 1)+1,1)],[meshPointMatrix(elementMatrix(i, 3)+1,2),meshPointMatrix(elementMatrix(i, 1)+1,2)], '-k');
% end
end 