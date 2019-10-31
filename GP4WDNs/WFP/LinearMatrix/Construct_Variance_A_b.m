function [AnalysisMatrix,B] = Construct_Variance_A_b(NumberofX,A,demand_MC)
% another method for small matrix: equationsToMatrix
ASparse = sparse(A');
[columns, lines,values] = find(ASparse);
columnsCell = 3;
CellIndex = 1;
CountIndex = 2;
ValueIndex = 3;

[m,~] = size (lines);
CellMatrix = cell(NumberofX,columnsCell);
for i = 1:NumberofX
    range=int32(find(lines==i));
    colind = columns(range,1);
    CellMatrix{i,CellIndex} = colind;
    [count,~]= size(colind);
    CellMatrix{i,CountIndex} = count;
    CellMatrix{i,ValueIndex} = values(range,1);
end



[m,n_d]=size(demand_MC);
PesudoDemand = zeros(NumberofX,n_d);
for i = 1:m
    PesudoDemand(i,:) = demand_MC(i,:);
end


nSparse = double((NumberofX+1)*NumberofX*0.5);
AnalysisMatrix = sparse(nSparse,nSparse);
B = sparse(nSparse,1);

for i = 1:NumberofX
    for j = i:NumberofX
        index =  (NumberofX + NumberofX-(i-2))*(i-1)/2+1+(j-i);
        expMatrix_i = CellMatrix{i,CellIndex};
        expMatrix_j = CellMatrix{j,CellIndex};
        CountMatrix_i = CellMatrix{i,CountIndex};
        CountMatrix_j = CellMatrix{j,CountIndex};
        ValueMatrix_i = CellMatrix{i,ValueIndex};
        ValueMatrix_j = CellMatrix{j,ValueIndex};
        
        TempMatrix = sparse(CountMatrix_i*CountMatrix_j,nSparse);
        tempindex = 1;
        for m = 1:CountMatrix_i
            for n = 1:CountMatrix_j
                vars2Value = [expMatrix_i(m) expMatrix_j(n)];
                minCoord = min(vars2Value);
                maxCoord = max(vars2Value);
                ind = (NumberofX + NumberofX-(minCoord-2))*(minCoord-1)/2+1+(maxCoord-minCoord);
                AnalysisMatrix(index,ind) = (AnalysisMatrix(index,ind)) + ValueMatrix_i(m)*ValueMatrix_j(n);
            end
        end
        cov_value = cov(PesudoDemand(i,:),PesudoDemand(j,:));
        B(index,1)= cov_value(1,2);
    end
end



end