function varIndex= findIndexByLabel(MeasuredHeadLable,Variable_Symbol_Table)
[m,n] = size(MeasuredHeadLable);
varIndex =[];
for i = 1:m
    varIndexTemp =[];
    for  j = 1:n
        var_index = find(strcmp(Variable_Symbol_Table,MeasuredHeadLable{i,j}))
        varIndexTemp = [varIndexTemp var_index];
    end
     varIndex = [varIndex;varIndexTemp];
end
end