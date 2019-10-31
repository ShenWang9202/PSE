function [A,b] = constructA_b(ForConstructA,ForConstructb,Demand,X0,base)

%% construct A
% demand part
DemandA_part1 = zeros(ForConstructA.JunctionCount,ForConstructA.h_end_index);
DemandA = [DemandA_part1 ForConstructA.MassMatrix];
% pipe part
FlowCount=ForConstructA.PipeCount+ForConstructA.PumpCount+ForConstructA.ValveCount;
PipeA_part2 =  zeros(ForConstructA.PipeCount,FlowCount);
for i=1:ForConstructA.PipeCount
    PipeA_part2(i,i) = -1;
end
PipeA = [ForConstructA.EnergyPipeMatrix PipeA_part2];
% pump part
% pump part
% getting initial values

PumpA_part2 = [];
PumpA=[];
C_1M =[];
PumpEquation = ForConstructb.PumpEquation;
if(~isempty(PumpEquation))
    h0_vector = PumpEquation(:,1);
    r_vector = PumpEquation(:,2);
    w_vector = PumpEquation(:,3);
    
    s = X0(ForConstructb.s_start_index:ForConstructb.s_end_index);
    q_pump = X0(ForConstructb.q_pump_start_index:ForConstructb.q_pump_end_index);
    
    C_1M = -h0_vector.*(s.^2);
    C_2M = -r_vector.* q_pump.^(w_vector-1).* s.^(2-w_vector);
    
    
    PumpA_part2 =  zeros(ForConstructA.PumpCount,FlowCount);
    for i=1:ForConstructA.PumpCount
        PumpA_part2(i,ForConstructA.PipeCount+i) = -C_2M(i);
    end
    PumpA = [ForConstructA.EnergyPumpMatrix PumpA_part2];
end
% valve part
% to do


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



%% construnct b

% demand part
b = [];
b = [b;Demand'];

% pipe part
C_estimate_b = CEstimateLinear_b(ForConstructb,X0);
b = [b;C_estimate_b];
% pump part
b = [b;C_1M];

% fixed tank part


b = [b; X0(TankHeadIndex)];

% fixed reservoir part
b = [b; X0(ReservoirHeadIndex)];









