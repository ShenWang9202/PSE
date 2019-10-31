function [AnalysisMatrix,B] = Construct_Variance_A_b_over_Pipe(NumberofX,A,PesudoDemand)
% another method for small matrix: equationsToMatrix

%check 
[NumberofEq,NumberofX_A] = size(A);
if(NumberofX_A ~= NumberofX)
    disp('something is wrong, stopping')
    return
end


ASparse = sparse(A');
[columns, lines,values] = find(ASparse);
columnsCell = 3;

CellIndex = 1;
CountIndex = 2;
ValueIndex = 3;

[m,~] = size (lines);
CellMatrix = cell(NumberofEq,columnsCell);
for i = 1:NumberofEq
    range=int32(find(lines==i));
    colind = columns(range,1);
    CellMatrix{i,CellIndex} = colind;
    [count,~]= size(colind);
    CellMatrix{i,CountIndex} = count;
    CellMatrix{i,ValueIndex} = values(range,1);
end




% Construct sparse matrix
mSparse = double((NumberofEq+1)*NumberofEq*0.5);
nSparse = double((NumberofX+1)*NumberofX*0.5);
AnalysisMatrix = sparse(mSparse,nSparse);
B = sparse(mSparse,1);


for i = 1:NumberofEq
    for j = i:NumberofEq
        EquationIndex =  (NumberofEq + NumberofEq-(i-2))*(i-1)/2+1+(j-i);
        %W(EquationIndex,EquationIndex) = w(i)*w(j);
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
                AnalysisMatrix(EquationIndex,ind) = (AnalysisMatrix(EquationIndex,ind)) + ValueMatrix_i(m)*ValueMatrix_j(n);
            end
        end
        cov_value = cov(PesudoDemand(i,:),PesudoDemand(j,:));
        B(EquationIndex,1)= cov_value(1,2);
    end
end



end