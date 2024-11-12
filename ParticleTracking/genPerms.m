function [ind]=genPerms(searchSpace,refSpace)
%This function scans all combinations of feature vectors in the
%deformed configuration to find the best matching (least error over all
%vectors) set to those in the reference configuration.

%Note: There is a better way to do this (using recursion?), but this function will work for
%search spaces up to 10 feature vectors

%Input arguments
%searchSpace = # of vectors to search in the deformed configuration
%refSpace = # of vectors to match in the reference configuration

%Output arguments
%ind = a matrix containing all possible permutations of indices for
%matching reference and deformed feature vectors
%   Author: Max Hockenberry
%   Last Update: 10/23/2024
counter=1; 
if refSpace==1 %match one feature vector
    for m=1:searchSpace
        ind(counter,1)=[m];
        counter=counter+1;
    end
end

if refSpace==2 %match two feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                ind(counter,1:2)=[m,n];
                counter=counter+1;
            end
        end
    end
end

if refSpace==3 %match three feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        ind(counter,1:3)=[m,n,o];
                        counter=counter+1;
                    end
                end
            end
        end
    end
end

if refSpace==4 %match four feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=n)
                        for p=1:searchSpace
                            if (p~=o)
                                ind(counter,1:4)=[m,n,o,p];
                                counter=counter+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

if refSpace==5 %match five feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        ind(counter,1:5)=[m,n,o,p,q];
                                        counter=counter+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
            
if refSpace==6 %match six feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        for r=1:searchSpace
                                            if (r~=m && r~=n && r~=o && r~=p && r~=q)
                                                ind(counter,1:6)=[m,n,o,p,q,r];
                                                counter=counter+1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%                             

if refSpace==7 %match seven  feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        for r=1:searchSpace
                                            if (r~=m && r~=n && r~=o && r~=p && r~=q)
                                                for s=1:searchSpace
                                                    if (s~=m && s~=n && s~=o && s~=p && s~=q && s~=r)
                                                        ind(counter,1:7)=[m,n,o,p,q,r,s];
                                                        counter=counter+1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if refSpace==8 %match eight feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        for r=1:searchSpace
                                            if (r~=m && r~=n && r~=o && r~=p && r~=q)
                                                for s=1:searchSpace
                                                    if (s~=m && s~=n && s~=o && s~=p && s~=q && s~=r)
                                                        for t=1:searchSpace
                                                            if (t~=m && t~=n && t~=o && t~=p && t~=q && t~=r && t~=s)
                                                                ind(counter,1:8)=[m,n,o,p,q,r,s,t];
                                                                counter=counter+1;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if refSpace==9 %match nine feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        for r=1:searchSpace
                                            if (r~=m && r~=n && r~=o && r~=p && r~=q)
                                                for s=1:searchSpace
                                                    if (s~=m && s~=n && s~=o && s~=p && s~=q && s~=r)
                                                        for t=1:searchSpace
                                                            if (t~=m && t~=n && t~=o && t~=p && t~=q && t~=r && t~=s)
                                                                for u=1:searchSpace
                                                                    if (u~=m && u~=n && u~=o && u~=p && u~=q && u~=r && u~=s && u~=t)
                                                                        ind(counter,1:9)=[m,n,o,p,q,r,s,t,u];
                                                                        counter=counter+1;
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if refSpace==10 %match ten feature vectors
    for m=1:searchSpace
        for n=1:searchSpace
            if (n~=m)
                for o=1:searchSpace
                    if (o~=m && o~=n)
                        for p=1:searchSpace
                            if (p~=m && p~=n && p~=o)
                                for q=1:searchSpace
                                    if (q~=m && q~=n && q~=o && q~=p)
                                        for r=1:searchSpace
                                            if (r~=m && r~=n && r~=o && r~=p && r~=q)
                                                for s=1:searchSpace
                                                    if (s~=m && s~=n && s~=o && s~=p && s~=q && s~=r)
                                                        for t=1:searchSpace
                                                            if (t~=m && t~=n && t~=o && t~=p && t~=q && t~=r && t~=s)
                                                                for u=1:searchSpace
                                                                    if (u~=m && u~=n && u~=o && u~=p && u~=q && u~=r && u~=s && u~=t)
                                                                        for v=1:searchSpace
                                                                            if (v~=m && v~=n && v~=o && v~=p && v~=q && v~=r && v~=s && v~=t && v~=u)
                                                                                ind(counter,1:10)=[m,n,o,p,q,r,s,t,u,v];
                                                                                counter=counter+1;
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end