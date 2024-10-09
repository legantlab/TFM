%TestScript for removing thermal exansion of beads. 

testMatches = matches;
testDisplacements = displacements;

[new_matches] = polyfitThermalExpansion(testMatches,testDisplacements,12 ,12, 12, 1);
new_Disps = cell(length(new_matches),1);
for i = 1:length(new_matches)
    curMatches = new_matches{i};
    new_Disps{i} = curMatches(:,4:6) - curMatches(:,1:3);
end

%% Compute Swelling using a regression approach 3D
regX = MultiPolyRegress(testMatches{22}(:,1:3),testDisplacements{22}(:,1),7,'figure','range');
%now remove the fit from the displacement data
testExpressionX = regX.PolynomialExpression(testMatches{22}(:,1), testMatches{22}(:,2), testMatches{22}(:,3));
testSwellingRemoved = testDisplacements{22}(:,1) - testExpressionX;
cellTestSwellingRemovedExp = regX.PolynomialExpression(testMatches{1}(:,1), testMatches{1}(:,2), testMatches{1}(:,3));
cellTestSwellingRemoved = testDisplacements{1}(:,1) - cellTestSwellingRemovedExp;

% figure %Plot swelling correct
% testPointsX = abs(randn(1000,3)*200);
% testSwellingExp = regX.PolynomialExpression(testPointsX(:,1), testPointsX(:,2), testPointsX(:,3));
% scatter(testPointsX(:,1), testSwellingExp)

figure
subplot(1,4,1)
scatter(testMatches{22}(:,1), displacements{22}(:,1))
title('X displacements along X pre-correction')
subplot(1,4,2)
scatter(testMatches{22}(:,1), testSwellingRemoved)
title('X displacements along X post-correction')
subplot(1,4,3)
scatter(testMatches{1}(:,1), displacements{1}(:,1))
title('X displacements along X pre-correction - Cell')
subplot(1,4,4)
scatter(testMatches{1}(:,1), cellTestSwellingRemoved)
title('X displacements along X post-correction - Cell')

%% Compute swelling along 1D for X
RegX1 = MultiPolyRegress(testMatches{22}(:,1),testDisplacements{22}(:,1),9,'figure');
evalRegX1 = RegX1.PolynomialExpression(testMatches{22}(:,1));

figure
subplot(1,3,1)
scatter(testMatches{22}(:,1), testDisplacements{22}(:,1))
title('Original Swelling in X')
subplot(1,3,2)
scatter(testMatches{22}(:,1), evalRegX1)
title('Computed Swelling Correction in X')
subplot(1,3,3)
scatter(testMatches{22}(:,1), testDisplacements{22}(:,1) - evalRegX1)
title('Computed Swelling Correction in X - Corrected')

%% Swelling correction using a cubic spline
xFit = fit((testMatches{22}(:,1)), smooth(testDisplacements{22}(:,1),16),'smoothingspline');


testsplineXFit = xFit(testMatches{22}(:,1));
figure
subplot(1,3,1)
plot(xFit, testMatches{22}(:,1), smooth(testDisplacements{22}(:,1),5))
ylim([min(smooth(testDisplacements{22}(:,1))), max(smooth(testDisplacements{22}(:,1)))])
title('Raw X Data and Cubic Spline Fit')
subplot(1,3,2)
plot(testMatches{22}(:,1), testDisplacements{22}(:,1) - testsplineXFit)
ylim([min(smooth(testDisplacements{22}(:,1))), max(smooth(testDisplacements{22}(:,1)))])
title('Corrected Data')
subplot(1,3,3)
plot(testMatches{1}(:,1), testDisplacements{1}(:,1) - xFit(testMatches{1}(:,1)))
ylim([min(smooth(testDisplacements{22}(:,1))), max(smooth(testDisplacements{22}(:,1)))])
title('Corrected Data - Cell')

%% Swelling Correction Subtraction interpolation 1D
%Interpolate displacements onto another grid
interpSwellCell = interp1(testMatches{22}(:,1), testDisplacements{22}(:,1), testMatches{1}(:,1));
figure
subplot(1,3,1)
scatter(testMatches{22}(:,1), testDisplacements{22}(:,1))
title('Raw Displacments')
subplot(1,3,2)
scatter(testMatches{1}(:,1), interpSwellCell)
title('Interpolated Swelling on Cell Displacements')
subplot(1,3,3)
scatter(testMatches{1}(:,1), testDisplacements{1}(:,1) - interpSwellCell)
title('Subtracted Swelling Displacements')

%% Swelling Correction Subtraction Interpolation 1D
interpSwellCellX = scatteredInterpolant(testMatches{22}(:,1), testMatches{22}(:,2),...
    testMatches{22}(:,3),testDisplacements{22}(:,1));
testSwellingX = interpSwellCellX(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3));
    
figure
subplot(1,3,1)
quiver3(testMatches{22}(:,1), testMatches{22}(:,2),testMatches{22}(:,3), ...
    testDisplacements{22}(:,1), testDisplacements{22}(:,2),testDisplacements{22}(:,3))
title('Swelling Displacements')
subplot(1,3,2)
quiver3(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3), ...
    testDisplacements{1}(:,1), testDisplacements{1}(:,2),testDisplacements{1}(:,3))
title('Cell Time Frame 1 Displacements')
subplot(1,3,3)
quiver3(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3), ...
    testDisplacements{1}(:,1) - testSwellingX, testDisplacements{1}(:,2),testDisplacements{1}(:,3))
title('Cell Time Frame 1 Displacements - X Corrected')

%% %% Swelling Correction Subtraction Interpolation 3D
interpSwellCellX = scatteredInterpolant(testMatches{22}(:,1), testMatches{22}(:,2),...
    testMatches{22}(:,3),testDisplacements{22}(:,1));
testSwellingX = interpSwellCellX(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3));

interpSwellCellY = scatteredInterpolant(testMatches{22}(:,1), testMatches{22}(:,2),...
    testMatches{22}(:,3),testDisplacements{22}(:,2));
testSwellingY = interpSwellCellY(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3));

interpSwellCellZ = scatteredInterpolant(testMatches{22}(:,1), testMatches{22}(:,2),...
    testMatches{22}(:,3),testDisplacements{22}(:,3));
testSwellingZ = interpSwellCellZ(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3));

figure
subplot(1,3,1)
quiver3(testMatches{22}(:,1), testMatches{22}(:,2),testMatches{22}(:,3), ...
    testDisplacements{22}(:,1), testDisplacements{22}(:,2),testDisplacements{22}(:,3))
title('Swelling Displacements')
subplot(1,3,2)
quiver3(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3), ...
    testDisplacements{1}(:,1), testDisplacements{1}(:,2),testDisplacements{1}(:,3))
title('Cell Time Frame 1 Displacements')
subplot(1,3,3)
quiver3(testMatches{1}(:,1), testMatches{1}(:,2),testMatches{1}(:,3), ...
    testDisplacements{1}(:,1) - testSwellingX, testDisplacements{1}(:,2) - testSwellingY, ...
    testDisplacements{1}(:,3) - testSwellingZ)
title('Cell Time Frame 1 Displacements - XYZ Corrected')

%% Analysis of which displacements/directions correspond to swelling must dominately. 
figure
subplot(1,3,1)
scatter(testMatches{1}(:,1),displacements{1}(:,1))
title('X displacements along X')
subplot(1,3,2)
scatter(testMatches{1}(:,1),displacements{1}(:,2))
title('Y displacements along X')
subplot(1,3,3)
scatter(testMatches{1}(:,1),displacements{1}(:,3))
title('Z displacements along X')

figure
subplot(1,3,1)
scatter(testMatches{1}(:,2),displacements{1}(:,1))
title('X displacements along Y')
subplot(1,3,2)
scatter(testMatches{1}(:,2),displacements{1}(:,2))
title('Y displacements along Y')
subplot(1,3,3)
scatter(testMatches{1}(:,2),displacements{1}(:,3))
title('Z displacements along Y')

figure
subplot(1,3,1)
scatter(testMatches{1}(:,3),displacements{1}(:,1))
title('X displacements along Z')
subplot(1,3,2)
scatter(testMatches{1}(:,3),displacements{1}(:,2))
title('Y displacements along Z')
subplot(1,3,3)
scatter(testMatches{1}(:,3),displacements{1}(:,3))
title('Z displacements along Z')



%%
%Computes difference in displacement values after expansion fitting of a
%specific time
testTime = 22;
xChange = mean(testDisplacements{testTime}(:,1) - new_Disps{testTime}(:,1)); 
disp(['Change in X Disp = ', num2str(xChange), ' um'])
%% Plot
testTime = 22;
figure
subplot(1,2,1)
quiver3(testMatches{testTime}(:,1), testMatches{testTime}(:,2),testMatches{testTime}(:,3), ...
    testDisplacements{testTime}(:,1), testDisplacements{testTime}(:,2), testDisplacements{testTime}(:,3),1)

hold on
subplot(1,2,2)
quiver3(new_matches{testTime}(:,1), new_matches{testTime}(:,2),new_matches{testTime}(:,3), ...
    new_Disps{testTime}(:,1), new_Disps{testTime}(:,2), new_Disps{testTime}(:,3),1)