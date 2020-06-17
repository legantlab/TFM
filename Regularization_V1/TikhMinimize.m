function [fMin, reg_new] = TikhMinimize(u,U,s,V,reg_corner,pseudoInv,fTikh)
fMin=fTikh;
error = norm(fTikh-pseudoInv)*norm(fTikh);
reg_new = reg_corner;
reg_nav = reg_corner/10^(.05);
reverse = 0;
for i = 1:0
    [i reg_nav]
    fMP = tikhonov(U,s,V,u,reg_nav,'Tikh');
    errorNav = norm(fMP - pseudoInv)*norm(fMP);
    if errorNav < error
        error = errorNav;
        fMin = fMP;
        reverse = 0;
        if reg_nav < reg_new
            reg_new = reg_nav;
            reg_nav = reg_new/10^(.05);
        else 
            reg_new = reg_nav;
            reg_nav = reg_new*10^(.05);
        end
    else
        if reverse %the maximum is at our current location
            fMin = fMP;
            break;
        else
            reg_nav = reg_new^2/reg_nav;
            reverse = 1;
        end
    end
end
end