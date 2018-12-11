%% add all dependencies of a file to git

function [fList, pList] = add2git(fname)

[fList, pList] = matlab.codetools.requiredFilesAndProducts(fname)

for i = 1:length(fList)
    sysstr = ['git add ',fList{i}];  
    system(sysstr);
end

end