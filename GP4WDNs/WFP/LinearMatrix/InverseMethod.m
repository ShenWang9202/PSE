% Inverse method

% Assume that each demand follow a normal distribution function, and we
% would like to see the the distribution of the states.
clc;
clear;
close all;
TestCase = 4;
[inpname, acc] = findInp_acc(TestCase);
[d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution,MassEnergyMatrix4GP,MC] = PrepareA_Monte_Carlo(inpname,TestCase);

%% sufficient scenario.
Demand_known=InitialParameter.Demand_known;
NumberofX= ForConstructA.NumofX;
HeadIndex = IndexInVar.JunctionHeadIndex;
FlowIndex = [IndexInVar.PumpFlowIndex IndexInVar.PipeFlowIndex IndexInVar.ValveFlowIndex];

deterministic = MC.determistic;
demand_MC = MC.demand_MC;
coefficient_MC = MC.coefficient_MC;
DemandIndex = MC.DemandIndex;
DemandVariance = MC.DemandVariance;
CoeffVariance= MC.CoeffVariance;
TankMeasurement_MC = MC.TankMeasurement_MC;
TankErrorVariance = MC.TankErrorVariance;
MCSolution = MC.MCSolution;

Variance = DemandVariance; % compatible with previous code when only considering demand uncertainty


ID_Index = Variable_Symbol_Table(:,1);
save('MonteCarloDemand.mat','DemandIndex','demand_MC','-v7');
save('MonteCarloData.mat','deterministic','HeadIndex','FlowIndex','ID_Index','MCSolution','DemandVariance','CoeffVariance','-v7');

PumpFlow = deterministic(IndexInVar.PumpFlowIndex);
PipeFlow = deterministic(IndexInVar.PipeFlowIndex);


%% Step 1, if all flows are positive, we can continue on, otherwise, update connection matrix at first.

NegativePipeIndex = find(PipeFlow<0);
MassEnergyMatrixStruct = UpdateConnectionMatrix(d,NegativePipeIndex);


%% Step 2 Get the solution of WFP, and linearizing around the solution

PipeFlow = abs(PipeFlow);
[K_pipe K_c] = LinearizePipe(ForConstructb,PipeFlow);
K_pump= [];
q = PumpFlow;
if(~isempty(IndexInVar.PumpEquation))
    r_vector = IndexInVar.PumpEquation(:,2);
    w_vector = IndexInVar.PumpEquation(:,3);
    K_pump = -r_vector.*w_vector.*(q.^(w_vector-1));
end

%% Step 3 Construct A and b

[A,A1,~] = Construct_H_Q_A_b_Pipe(MassEnergyMatrixStruct,ForConstructA,Demand_known,K_pipe,K_c,K_pump);

% index for each element.
JunctionCount = ForConstructA.JunctionCount;
PipeCount = ForConstructA.PipeCount;
PumpCount = ForConstructA.PumpCount;
TankCount = ForConstructA.TankCount;
ReservoirCount = ForConstructA.ReservoirCount;


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
variance_MC = zeros(NumberofX,1);
variance_MC(JunctionHeadIndex,:) = DemandVariance;
variance_MC(PipeFlowIndex,:) = CoeffVariance;
variance_MC(TankHeadIndex,:) = TankErrorVariance;

KC = ones(1,NumberofX);
KC(1,PipeFlowIndex) = K_c;
variance_MC = variance_MC';
A_inverse = A^(-1);

MC = MCSolution(1:NumberofX,:);
cov_MC = cov(MC');

Sigma_LM = [];
Sigma_MC = [];
for i = 1:NumberofX
    Sigma_LM = [Sigma_LM;sqrt(sum((A_inverse(i,:).*KC).^2.*variance_MC))];
    Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
end

Sigma_error = Sigma_MC - Sigma_LM
