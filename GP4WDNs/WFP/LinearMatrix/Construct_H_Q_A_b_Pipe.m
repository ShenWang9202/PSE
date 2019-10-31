function [A,A1,b] = Construct_H_Q_A_b_Pipe(MassEnergyMatrixStruct,ForConstructA,Demand_known,K_pipe,K_c,K_pump)

%% Construct A and b matrix
%A = [DemandA;PipeA;PumpA;TankA;ReservoirA];

MassMatrix = MassEnergyMatrixStruct.MassMatrix;
EnergyPipeMatrix = MassEnergyMatrixStruct.EnergyPipeMatrix;
EnergyPumpMatrix = MassEnergyMatrixStruct.EnergyPumpMatrix;
EnergyValveMatrix = MassEnergyMatrixStruct.EnergyValveMatrix;

% Construct the A matrix for demand
Variance_Demand_A = MassMatrix;
DemandA_part1 = zeros(ForConstructA.JunctionCount,ForConstructA.h_end_index);
DemandA = [DemandA_part1 Variance_Demand_A];

% Construct the A matrix for pipes
% pipe part
FlowCount=ForConstructA.PipeCount+ForConstructA.PumpCount+ForConstructA.ValveCount;
PipeA_part2 =  zeros(ForConstructA.PipeCount,FlowCount);
PipeA1_part2 =  zeros(ForConstructA.PipeCount,FlowCount);
for i=1:ForConstructA.PipeCount
    PipeA_part2(i,i) = -K_pipe(i);
    PipeA1_part2(i,i) = -K_pipe(i)/1.852;
end
PipeA = [EnergyPipeMatrix PipeA_part2];
PipeA1 = [EnergyPipeMatrix PipeA1_part2];

% Construct the A matrix for pumps

PumpA_part2 =  zeros(ForConstructA.PumpCount,FlowCount);
for i=1:ForConstructA.PumpCount
    PumpA_part2(i,ForConstructA.PipeCount+i) = -K_pump(i);
end
PumpA = [EnergyPumpMatrix PumpA_part2];

% Construct the A matrix for valves TODO

% fixed tank part

TankHeadIndex = ForConstructA.TankHeadIndex;
[~,m] = size(TankHeadIndex);
[~,n] = size(DemandA);
TankA = zeros(m,n);
for i=1:m
    TankA(i,TankHeadIndex(i)) = 1;
end

% fixed reservoir part

ReservoirHeadIndex = ForConstructA.ReservoirHeadIndex;
[~,m] = size(ReservoirHeadIndex);
[~,n] = size(DemandA);
ReservoirA = zeros(m,n);
for i=1:m
    ReservoirA(i,ReservoirHeadIndex(i)) = 1;
end

% construct A done

A = [DemandA;PipeA;PumpA;TankA;ReservoirA];
A1 = [DemandA;PipeA1;PumpA;TankA;ReservoirA];



%% construnct b

% actually, we don't need b at all, no need to construct it if no
% validation of our linearizion  is performed


% demand part 
b = [];
b = [b;Demand_known';];

VarianceofLinearization = zeros(ForConstructA.NumofX-ForConstructA.JunctionCount,1);

b = [b;VarianceofLinearization];

end