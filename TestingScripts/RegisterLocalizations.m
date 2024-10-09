%Register localizations to each other

swellingControlData = matches{2}(:,1:3);
s5Data = s5Matches{23}(:,1:3);

figure
title('Original Localizations')
scatter3(swellingControlData(:,1)*size_x, swellingControlData(:,2)*size_y, swellingControlData(:,3)*size_z,'r')
hold on 
scatter3(s5Data(:,1) - 1.5, s5Data(:,2), s5Data(:,3),'b')

