%Script to work out Wes' old tracking code
% ind = track{t-1}{1} ~= 0;
% map = track{t-1}{1}(ind);
% curMatchedBeads= x{t}{1}(map,:);
% curUnmatchedBeads = setxor(curMatchedBeads,x{t}{1},'rows');
% curUnmatchedRef = x{1}{1}(~ind,:);
% 
% Idx = knnsearch(curUnmatchedRef,curUnmatchedBeads);
% matchedRefBeads = curUnmatchedRef(Idx,:);
% distNewMatches = matchedRefBeads - curUnmatchedBeads;

%% Load a data set
close all
load('ExampleDisplacementField3.mat')
%t =3 ;
refBeads = [(1:length(x{1}{1}))',x{1}{1}];
strainBeads = [(1:length(x{3}{1}))',x{3}{1}];

figure
scatter(refBeads(:,2), refBeads(:,3), 'o','r')
hold on 
scatter(strainBeads(:,2), strainBeads(:,3), 'o','g')
axis equal

%% Do matching
tic
firstPass = displacementCalc_recursive_Max_2024_02_19(refBeads,strainBeads,10,2,3,2);

distNewMatches = firstPass(:,6:8) - firstPass(:,2:4);

% figure
% quiver(firstPass(:,6), firstPass(:,7),distNewMatches(:,1), distNewMatches(:,2),0)
% hold on 
% scatter(firstPass(:,2), firstPass(:,3), 'r','.')
% scatter(firstPass(:,6), firstPass(:,7), 'g','.')
% axis equal

%Implement a scheme where we compute the drift and apply that to the raw
%displacements then retrack...
xDrift = median(distNewMatches(:,1));
yDrift = median(distNewMatches(:,2));
zDrift = median(distNewMatches(:,3));
refBeadsCor = refBeads;
%refBeadsCor(:,2) = refBeadsCor(:,2) + xDrift;
%refBeadsCor(:,3) = refBeadsCor(:,3) + yDrift;
%refBeadsCor(:,4) = refBeadsCor(:,4) + zDrift;

strainBeadsCor = strainBeads;
strainBeadsCor(:,2) = strainBeadsCor(:,2) - xDrift;
strainBeadsCor(:,3) = strainBeadsCor(:,3) - yDrift;
strainBeadsCor(:,4) = strainBeadsCor(:,4) - zDrift;

% figure
% scatter(refBeadsCor(:,2), refBeadsCor(:,3), 'r','.')
% hold on 
% scatter(strainBeadsCor(:,2), strainBeadsCor(:,3), 'g','.')
% axis equal

dcMatches = displacementCalc_recursive_Max_2024_02_19(refBeadsCor,strainBeadsCor,8,3,5,1.5);

figure
scatter(dcMatches(:,2), dcMatches(:,3),'.','r')
hold on
scatter(dcMatches(:,6), dcMatches(:,7),'.','g')
dcDisps = dcMatches(:,6:8) - dcMatches(:,2:4);
quiver(dcMatches(:,2), dcMatches(:,3),dcDisps(:,1), dcDisps(:,2),0)
axis equal

%% Run smoothing function

matchedCoordinates_filt=smoothFilt_2024_02_20(dcMatches,20,1.5); 
filt_Disps = matchedCoordinates_filt(:,6:8) - matchedCoordinates_filt(:,2:4);
disp(['Filtered Beads = ',num2str(length(dcMatches) - length(matchedCoordinates_filt))])
matchedCoordinates_filt2=smoothFilt_2024_02_20(matchedCoordinates_filt,20,2); 

filt_Disps2 = matchedCoordinates_filt2(:,6:8) - matchedCoordinates_filt2(:,2:4);
disp(['Filtered Beads = ',num2str(length(matchedCoordinates_filt) - length(matchedCoordinates_filt2))])

figure
 scatter(refBeadsCor(:,2), refBeadsCor(:,3), 'r','o')
hold on 
scatter(strainBeadsCor(:,2), strainBeadsCor(:,3), 'g','o')
scatter(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3), 'r','.')
hold on
scatter(matchedCoordinates_filt2(:,6), matchedCoordinates_filt2(:,7), 'g','.')
quiver(matchedCoordinates_filt2(:,2), matchedCoordinates_filt2(:,3),filt_Disps2(:,1), filt_Disps2(:,2),0)
axis equal

%% 

%Load raw image for comparison
% I = loadtiff('S:\Max\2024-01-25\tiffs\F07\SubSample\Time_0002.ome.tiff');
% sizeI = size(I);
% figure
% xlim = [-1*sizeI(2)/2, sizeI(2)/2];
% ylim = [-1*sizeI(1)/2, sizeI(1)/2];
% imagesc(xlim, ylim, max(I,[],3),clim);
% imagesc(max(I,[],3))
% imagesc()
% figure
% quiver(matchedCoordinates_filt(:,2), matchedCoordinates_filt(:,3),filt_Disps2(:,1), filt_Disps2(:,2),0)
% hold on 
% scatter(matchedCoordinates_filt(:,2), matchedCoordinates_filt(:,3), 'r','.')
% scatter(matchedCoordinates_filt(:,6), matchedCoordinates_filt(:,7), 'g','.')
% scatter(refBeadsCor(:,2), refBeadsCor(:,3), 'o','r')
% scatter(strainBeadsCor(:,2), strainBeadsCor(:,3), 'o','g')
% legend('Matched Displacements','Matched Reference','Matched Strained','Total Reference','Total Strained')
% axis equal

toc