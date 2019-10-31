% Extend Period using Inverse method

% Assume that each demand, measurement of tank, pipe coefficient follow a normal distribution function, and we
% would like to see the the distribution of the states.
clc;
clear;
close all;
TestCase = 5;
[inpname, acc] = findInp_EPS(TestCase);
[d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution,MassEnergyMatrix4GP,Inverse] = PrepareA_EPS(inpname,TestCase);

%% sufficient scenario.
Demand_known=InitialParameter.Demand_known;
NumberofX= ForConstructA.NumofX;
HeadIndex = IndexInVar.JunctionHeadIndex;
FlowIndex = [IndexInVar.PumpFlowIndex IndexInVar.PipeFlowIndex IndexInVar.ValveFlowIndex];

deterministic = Inverse.determistic;
%DemandIndex = Inverse.DemandIndex;
DemandVariance = Inverse.DemandVariance;
CoeffVariance= Inverse.CoeffVariance;
TankErrorVariance = Inverse.TankErrorVariance;

%% inverse
[~,time] = size(deterministic);
Variance_EPS = [];
for t = 1:time
    PumpFlow = deterministic(IndexInVar.PumpFlowIndex,t);
    PipeFlow = deterministic(IndexInVar.PipeFlowIndex,t);
    
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
    
    [A,A1,~] = Construct_H_Q_A_b_Pipe(MassEnergyMatrixStruct,ForConstructA,Demand_known(t,:),K_pipe,K_c,K_pump);
    
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
    
    % create variance
    variance_MC = zeros(NumberofX,1);
    variance_MC(JunctionHeadIndex,:) = DemandVariance(t,:);
    variance_MC(PipeFlowIndex,:) = CoeffVariance(t,:);
    variance_MC(TankHeadIndex,:) = TankErrorVariance(t,:);
    
    KC = ones(1,NumberofX);
    KC(1,PipeFlowIndex) = K_c;
    variance_MC = variance_MC';
    
    % find inverse
    A_inverse = A^(-1);
    
    Sigma_LM = [];
    Sigma_MC = [];
    for i = 1:NumberofX
        Sigma_LM = [Sigma_LM;sqrt(sum((A_inverse(i,:).*KC).^2.*variance_MC))];
    end
    Variance_EPS = [Variance_EPS Sigma_LM];
end
Variance_EPS

%% Plot each one to check
% ConfidenceLevel1 = 1.282;
% ConfidenceLevel2 = 1.96;
% for id = 1:17
% expectation = deterministic(id,1:24);
% stand_deviation = Variance_EPS(id,1:24);
% PlotConfidienceIntervalforid(expectation,stand_deviation,ConfidenceLevel1,ConfidenceLevel2,int2str(id));
% end
%PlotConfidienceInterval(expectation,stand_deviation,ConfidenceLevel)
%% Plot all selected head

ConfidenceLevel1 = 1.282;
ConfidenceLevel2 = 1.96;

SelectedHead = {'J3' 'J5' 'T8'};

expectation = deterministic(:,1:24);
stand_deviation = Variance_EPS(:,1:24);
PlotConfidienceInterval4SelectedHead(expectation,stand_deviation,ConfidenceLevel1,ConfidenceLevel2,SelectedHead,Variable_Symbol_Table);


%% Plot selected flow
SelectedHead = {'P2' 'P5' 'P8' 'PU9'};

expectation = deterministic(:,1:24);
stand_deviation = Variance_EPS(:,1:24);
PlotConfidienceInterval4SelectedFlow(expectation,stand_deviation,ConfidenceLevel1,ConfidenceLevel2,SelectedHead,Variable_Symbol_Table);

