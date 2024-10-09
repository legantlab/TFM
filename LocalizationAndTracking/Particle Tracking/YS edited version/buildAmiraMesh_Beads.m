function text = buildAmiraMesh_Beads(coordinates)

text = ['# AmiraMesh ASCII 1.0 \n\ndefine Nodes ' num2str(size(coordinates,1)+1) ...
    '\ndefine Tetrahedra ' num2str(size(coordinates,1)-2) ...
    '\nNodes { float[3] Coordinates } = @1\nTetrahedra { int[4] Nodes } = @4'];

coordString = '0 0 -10';
tetraString = '';
% vectorString = '0 0 0';
for i = 1:size(coordinates,1)
    coordString = [coordString '\n' num2str(coordinates(i,1)) ' ' num2str(coordinates(i,2)) ...
        ' ' num2str(coordinates(i,3))];
%     vectorString = [vectorString '\n' num2str(vectors(3*i-2)) ' ' num2str(vectors(3*i-1)) ...
%         ' ' num2str(vectors(3*i))];
    if i > 2
        tetraString = [tetraString '\n1 ' num2str(i-2) ' ' num2str(i-1) ' ' num2str(i)];
    end
end

text = [text '@1\n' coordString '\n\n@4' tetraString];

end

