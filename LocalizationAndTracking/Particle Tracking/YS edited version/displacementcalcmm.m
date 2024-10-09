%calculate displacements, reference and deformed coordinates from Franck T-PT algorithm
referencecoord = cell(length(track{1,1}{1,1}),1);
deformedcoord = cell(length(track{1,1}{1,1}),1);
for i = 1:length(track{1,1}{1,1})-1
    if track{1,1}{1,1}(i) > 0;
        referencecoord{i} = x{1,1}{1,1}(i,:);
        deformedcoord{i} = x{1,2}{1,1}(track{1,1}{1,1}(i),:);
    else track{1,1}{1,1}(i) == 0;
        i = i+1;
    end
end

%remove empty cells
refcoord_empty_rem = referencecoord;
defcoord_empty_rem = deformedcoord;
emptyCells_r = cellfun('isempty', refcoord_empty_rem);
emptyCells_d = cellfun('isempty', defcoord_empty_rem);
refcoord_empty_rem(all(emptyCells_r,2),:) = [];
defcoord_empty_rem(all(emptyCells_d,2),:) = [];

%convert cells to matrix
refcoord = cell2mat(refcoord_empty_rem);
defcoord = cell2mat(defcoord_empty_rem);

%subtract out means of x, y, z and obtain displacements
refcoord_n = refcoord - mean(refcoord);
defcoord_n = defcoord - mean(defcoord);
disp = refcoord_n - defcoord_n;

% reformat ie main
disp2 = disp';
disp2 = disp2(:);
refcoord_n2 = refcoord_n';
refcoord_n2 = refcoord_n2(:);
defcoord_n2 = defcoord_n';
defcoord_n2 = defcoord_n2(:);

%plot
disp('plotting')
figure(2)
%quiver3(refcoord_n2(1:3:end),refcoord_n2(2:3:end),refcoord_n2(3:3:end),disp2(1:3:end),disp2(2:3:end),disp2(3:3:end))
quiver(refcoord_n2(1:3:end),refcoord_n2(2:3:end),disp2(1:3:end),disp2(2:3:end))
xlabel('x')
ylabel('y')
zlabel('z')
title('Displacements Used to Solve Inverse Problem')
axis equal

%csvwrite('I_0.csv',refcoord);
%save('refcoord','refcoord');
%save('refcoord_n','refcoord_n');
csvwrite('coords_um.csv',refcoord_n);
save('disp2','disp2');
save('refcoord_n2','refcoord_n2');
%csvwrite('refcoord_n',refcoord_n);