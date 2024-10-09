function [repNodes]=reportNodeFind(repElems,elemList)
repElemsSort=unique(sort(repElems));
repNodes=elemList(repElemsSort,:);
repNodes=repNodes(:,2:end);
repNodes=unique(repNodes(:));
end