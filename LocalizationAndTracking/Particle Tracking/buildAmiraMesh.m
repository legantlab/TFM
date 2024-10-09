function text = buildAmiraMesh(vectors,coordinates)

text = ['# AmiraMesh ASCII 1.0 \n\ndefine Nodes ' num2str(size(coordinates,1)+1) ...
    '\ndefine Tetrahedra ' num2str(size(coordinates,1)-2) ...
    '\nNodes { float[3] Coordinates } = @1\nTetrahedra { int[4] Nodes } = @4' ...
    '\nNodes { float[3] Data } = @8\nField { float[3] Example } = Linear(@8)\n\n'];

coordString = '0 0 0';
tetraString = '';
vectorString = '0 0 0';
for i = 1:size(coordinates,1)
    coordString = [coordString '\n' num2str(coordinates(i,1)) ' ' num2str(coordinates(i,2)) ...
        ' ' num2str(coordinates(i,3))];
    vectorString = [vectorString '\n' num2str(vectors(i,1)) ' ' num2str(vectors(i,2)) ...
        ' ' num2str(vectors(i,3))];
    if i > 2
        tetraString = [tetraString '\n1 ' num2str(i-2) ' ' num2str(i-1) ' ' num2str(i)];
    end
end

text = [text '@1\n' coordString '\n\n@4' tetraString '\n\n@8\n' vectorString];

end

