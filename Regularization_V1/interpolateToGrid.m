function gridVector = interpolateToGrid(nodeLocation, nodeVector, gridLocation)

gridVector = zeros(size(gridLocation));
for i = 1:size(gridLocation,1)
    thisGridLocation = gridLocation(i,:);
    nodeLocTemp = nodeLocation;
    
    [d,closestNodeIndex1] = min(sum(sqrt((nodeLocTemp - thisGridLocation).^2),2));
    if d < 1e-6
        gridVector(i,:) = nodeVector(closestNodeIndex1,:);
        continue;
    end
    
    nodeLocTemp(closestNodeIndex1,:) = 1000000;
    [d,closestNodeIndex2] = min(sum(sqrt((nodeLocTemp - thisGridLocation).^2),2));
    
    nodeLocTemp(closestNodeIndex2,:) = 1000000;
    [d,closestNodeIndex3] = min(sum(sqrt((nodeLocTemp - thisGridLocation).^2),2));
    
    nodeLocTemp(closestNodeIndex3,:) = 1000000;
    [d,closestNodeIndex4] = min(sum(sqrt((nodeLocTemp - thisGridLocation).^2),2));
    
    closestLocation1 = nodeLocation(closestNodeIndex1,:);
    closestLocation2 = nodeLocation(closestNodeIndex2,:);
    closestLocation3 = nodeLocation(closestNodeIndex3,:);
    closestLocation4 = nodeLocation(closestNodeIndex4,:);
    
    closestVector1 = nodeVector(closestNodeIndex1,:);
    closestVector2 = nodeVector(closestNodeIndex2,:);
    closestVector3 = nodeVector(closestNodeIndex3,:);
    closestVector4 = nodeVector(closestNodeIndex4,:);
    
    A = [closestLocation1 1; closestLocation2 1; closestLocation3 1; closestLocation4 1];
    bx = [closestVector1(1); closestVector2(1); closestVector3(1); closestVector4(1)];
    coeffx = A\bx;
    interpolator_x = sum(coeffx' .* [thisGridLocation 1]);
    
    by = [closestVector1(2); closestVector2(2); closestVector3(2); closestVector4(2)];
    coeffy = A\by;
    interpolator_y = sum(coeffy' .* [thisGridLocation 1]);
    
    bz = [closestVector1(3); closestVector2(3); closestVector3(3); closestVector4(3)];
    coeffz = A\bz;
    interpolator_z = sum(coeffz' .* [thisGridLocation 1]);
    
    gridVector(i,:) = [interpolator_x interpolator_y interpolator_z];
    %gridVector(i,:) = closestVector1 + (closestVector2 - closestVector1).*norm(thisGridLocation-closestLocation1)/norm(closestLocation1-closestLocation2);
    %gridVector(i,:) = closestVector1 + (closestVector2 - closestVector1).*(thisGridLocation-closestLocation1)./(closestLocation2-closestLocation1);
end

end

function mag = minimize(x,thisGridLocation,nodeLocation)
    i = round(x);
    mag = sum(sqrt((nodeLocation(i) - thisGridLocation).^2),2);
end