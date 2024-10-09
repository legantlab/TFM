function gridVector = interpolateToGridTemp(nodeLocation, nodeVector, gridLocation)

vector_x = interp3(nodeLocation.x,nodeLocation.y,nodeLocation.z,V,gridLocation(:,1),gridLocation(:,2),gridLocation(:,3));
vector_y = interp3(nodeLocation.x,nodeLocation.y,nodeLocation.z,nodeVector(:,2),gridLocation(:,1),gridLocation(:,2),gridLocation(:,3));
vector_z = interp3(nodeLocation.x,nodeLocation.y,nodeLocation.z,nodeVector(:,3),gridLocation(:,1),gridLocation(:,2),gridLocation(:,3));

gridVector = zeros(3*size(vector_x,1),1);

end
