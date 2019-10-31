function [AnalysisMatrix,B] = Construct_Variance_A_b_Pipe(NumberofX,A,demand_MC,coefficient_MC,TankMeasurement_MC,K_c,ForConstructA)
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

% index for each element.
JunctionCount = ForConstructA.JunctionCount;
PipeCount = ForConstructA.PipeCount;
PumpCount = ForConstructA.PumpCount;
TankCount = ForConstructA.TankCount;
ReservoirCount = ForConstructA.ReservoirCount;

% calculate the index, the order must be the same as the order in
% Construct_H_Q_A_b_Pipe.m which is A = [DemandA;PipeA;PumpA;TankA;ReservoirA];
JunctionHeadIndex = 1:JunctionCount;
BaseCount4Next = JunctionCount;
PipeFlowIndex = (BaseCount4Next + 1):(BaseCount4Next + PipeCount);

BaseCount4Next = BaseCount4Next + PipeCount;
PumpFlowIndex = (BaseCount4Next + 1):(BaseCount4Next + PumpCount);

BaseCount4Next = BaseCount4Next + PumpCount;
TankHeadIndex = (BaseCount4Next + 1):(BaseCount4Next + TankCount);

BaseCount4Next = BaseCount4Next + TankCount;
ReservoirHeadIndex = (BaseCount4Next + 1):(BaseCount4Next + ReservoirCount);


[m,n_d]=size(demand_MC);
PesudoDemand = zeros(NumberofX,n_d);
PesudoDemand(JunctionHeadIndex,:) = demand_MC;
coefficient_MC = coefficient_MC';
PesudoDemand(PipeFlowIndex,:) = K_c'.*coefficient_MC;

% no need to plus the elevation of tanks, because we are calculating
% covariance, the constant valvue doesn't make any sense
PesudoDemand(TankHeadIndex,:) = TankMeasurement_MC;

% assume the elecation of reservoirs has no uncertainty 
%PesudoDemand(ReservoirHeadIndex,:) = ReservoirMeasurement_MC;


nSparse = double((NumberofX+1)*NumberofX*0.5);
AnalysisMatrix = sparse(nSparse,nSparse);

cov_valu = cov(PesudoDemand');
B = sparse(cov_valu(triu(true(size(cov_valu)))));

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
    end
end



end