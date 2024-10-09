%Shape Function Test

x1 = rand(1);
y1 = rand(1);
z1 = rand(1);

x2 = rand(1);
y2 = rand(1);
z2 = rand(1);

x3 = rand(1);
y3 = rand(1);
z3 = rand(1);

ons = [1 1 1];
area = 0.5*sqrt(det([[x1,x2,x3];[y1,y2,y3];ons])^2 + det([[y1,y2,y3];[z1,z2,z3];ons])^2 + det([[z1,z2,z3];[x1,x2,x3];ons])^2);


scatter3(x1,y1,z1)
hold on
scatter3(x2,y2,z2)
scatter3(x3,y3,z3)

x = (x1 + x2 + x3)/3;
y = (y1 + y2 + y3)/3;
z = (z1 + z2 + z3)/3;

scatter3(x,y,z)

b = (y3-y + (((y1-y3)*x) - (y1-y3)*x3)/(x1-x3))/(((y1-y3)*(x2-x3))/(x1-x3) - (y2-y3));
a = (x-x3 - (x2-x3)*b)/(x1-x3);

xtest = (a*(x1) + b*(x2) + (1-a-b)*(x3));
ytest = (a*(y1) + b*(y2) + (1-a-b)*(y3));
ztest = (a*(z1) + b*(z2) + (1-a-b)*(z3));


        %updatedSub = zeros(length(connectMatrix),1);

rowind = 1;
test4 = [];
tic
           for j = 1:length(elemCents2D) 
               cursub = nodeDisps(numNodes*(i-1)+1:numNodes*i,1:5);
            curElemNodes = TR2.ConnectivityList(j,:);
%             N1 = TR2.Points(curElemNodes(1),:);
%             N2 = TR2.Points(curElemNodes(2),:);
%             N3 = TR2.Points(curElemNodes(3),:);
% 
%             elemCent = elemCents2D(j,:);
%                 x1 = N1(1);
%                 y1 = N1(2);
%                 z1 = N1(3);
%                 
%                 x2 = N2(1);
%                 y2 = N2(2);
%                 z2 = N2(3);
%                 
%                 x3 = N3(1);
%                 y3 = N3(2);
%                 z3 = N3(3);
%                 x = elemCent(1);
%                 y = elemCent(2);
%             b = (y3-y + (((y1-y3)*x) - (y1-y3)*x3)/(x1-x3))/(((y1-y3)*(x2-x3))/(x1-x3) - (y2-y3));
%             a = (x-x3 - (x2-x3)*b)/(x1-x3);

            N1Disp = curSub(curElemNodes(1) == curSub(:,1),3:5);
            N2Disp = curSub(curElemNodes(2) == curSub(:,1),3:5);
            N3Disp = curSub(curElemNodes(3) == curSub(:,1),3:5);

          
            xDisp = ((1/3)*(N1Disp(1)) + (1/3)*(N2Disp(1)) + (1/3)*(N3Disp(1)));
            yDisp = ((1/3)*(N1Disp(2)) + (1/3)*(N2Disp(2)) + (1/3)*(N3Disp(2)));
            zDisp = ((1/3)*(N1Disp(3)) + (1/3)*(N2Disp(3)) + (1/3)*(N3Disp(3)));

            

            test4 = [test4,[xDisp, yDisp, zDisp]];


           end
toc