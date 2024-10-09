function text = buildAmiraMeshNew(vectors,coordinates)

text = ['# AmiraMesh 3D ASCII 3.0 \n \nnTriangles ' num2str(size(coordinates,1)) ...
    '\n\nParameters {\n     ContentType "SurfaceField",\n     Encoding "OnTriangles" \n\n}' ...
    '\nTriangleData { float values } @1\n\n' '# Data section follows'];

coordString = '';
vecNorm = vecnorm(vectors,2,2);
for i = 1:size(coordinates,1)
    coordString = [coordString '\n' num2str(vecNorm(i))];
end
% for i = 1:size(coordinates,1)
%     coordString = [coordString '\n' num2str(coordinates(i,1)) ' ' num2str(coordinates(i,2)) ...
%         ' ' num2str(coordinates(i,3))];
%     vectorString = [vectorString '\n' num2str(vectors(3*i-2)) ' ' num2str(vectors(3*i-1)) ...
%         ' ' num2str(vectors(3*i))];
%     if i > 2
%         tetraString = [tetraString '\n ' num2str(i-2) ' ' num2str(i-1) ' ' num2str(i)];
%     end
% end

text = [text '\n@1' coordString];

end