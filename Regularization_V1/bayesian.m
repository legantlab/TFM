function [fMin, alpha, beta] = bayesian(G,u,reg_corner,fTikh,displacedPtsMat,UTilde,sTilde,VTilde,MTilde,w)

m = length(u);
uTilde = u - sum(u/3/m);
%Mbar = 1/3/m * sum(G);
%w = sqrt(1/(3*m-1) *sum((G-Mbar).^2));
%MTilde = (G-Mbar)./w;
%[UTilde,sTilde,VTilde] = csvd(MTilde);

% calculate beta
ufar = uTilde(sqrt(displacedPtsMat(:,2).^2+displacedPtsMat(:,3).^2) > .2);
beta = var(ufar);

Ef = @(f) (f'*f)/2;
Eu = @(f) (MTilde*f-uTilde)'*(MTilde*f-uTilde)/2;
K = @(f,a,b) a*Ef(f) + b*Eu(f);
%f_xi = @(f,a,b,i) (-K(f,a,b)+K(f+1e-11*(i==(1:size(MTilde,2))'),a,b))/1e-11;
%aij = @(f,a,b,i,j) (-f_xi(f,a,b,i)+f_xi(f+1e-11*(j==(1:size(MTilde,2))'),a,b,i))/1e-11;
aij = @(f,a,b,i,j) (K(f+1e-11*(i==(1:size(MTilde,2)) | j==(1:size(MTilde,2)))',a,b) - ...
                    K(f+1e-11*(j==(1:size(MTilde,2)))',a,b) - K(f+1e-11*(i==(1:size(MTilde,2)))',a,b) + ... 
                    K(f,a,b))/(1e-11)^2;

likelihood = @(f,a,b,A) -a*Ef(f) - b*Eu(f) - .5*log(det(A)) + size(MTilde,2)/3*log(a) ...
    + size(MTilde,1)/3*log(b) - size(MTilde,1)/3*log(2*pi);

alpha = reg_corner*beta;
fMP = tikhonov(UTilde,sTilde,VTilde,uTilde,alpha/beta,'Tikh');
A = zeros(size(MTilde,2));
for i = 1:size(MTilde,2)
    for j = 1:i
        entry = aij(fMP,alpha,beta,i,j);
        A(i,j) = entry;
        A(j,i) = entry;
    end
end
det(A)
Ef(fMP)
Eu(fMP)
asdf
currLikelihood = likelihood(fMP,alpha,beta,A)
fMin = fMP

aNav = alpha/10^(.001);
reverse = 0;
for i = 1:10
    [alpha beta]
    alpha/beta
    fMP = tikhonov(UTilde,sTilde,VTilde,uTilde,aNav/beta,'Tikh');
    A = zeros(size(MTilde,2));
    for i = 1:size(MTilde,2)
        for j = 1:i
            entry = aij(fMP,aNav,beta,i,j);
            A(i,j) = entry;
            A(j,i) = entry;
        end
    end
    l = likelihood(fMP,aNav,beta,A)
    if l > currLikelihood
        currLikelihood = l;
        alpha = aNav;
        fMin = fMP;
        reverse = 0;
        if aNav > alpha
            aNav = alpha/10^(.001);
        else 
            aNav = alpha*10^(.001);
        end
    else
        if reverse %the maximum is at our current location
            fMin = fMP;
            break;
        else
            aNav = alpha^2/aNav;
            reverse = 1;
        end
    end
end
fMin = fMin./w;
%{
for alpha = logspace(-12,1,10)
    alpha
    for b = logspace(-12,1,10)
        b
        fMP = tikhonov(U,s,V,u,alpha/b,'Tikh');
        A = zeros(size(G,2));
        for i = 1:size(G,2)
            for j = 1:i
                entry = aij(fMP,alpha,b,i,j);
                A(i,j) = entry;
                A(j,i) = entry;
            end
        end
        l = likelihood(fMP,alpha,b,A);
        if l > maxL
            maxL = l;
            alpha = alpha;
            beta = b;
            fMin = fMP;
        end
    end
end
%}


end