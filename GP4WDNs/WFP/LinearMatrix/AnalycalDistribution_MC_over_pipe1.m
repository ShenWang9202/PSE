% Get the variances of the distribution of each demand(not zero)
load('MonteCarloData.mat')
load('MonteCarloDemand.mat')

CoeffVariance
DemandVariance
DemandIndex
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
%% remove 1.853 in K_pipe and verify the following
% q = PipeFlow;
% RealHeadloss = Headloss_pipe_R.* (q).^(1.852);
% LinearHeadloss = K_pipe.*q./1.852;
% b_pipe = RealHeadloss - LinearHeadloss;
% b = [Demand_known;];

A_save = A;
sigma = [];

Bcollect = [];
Wcollect = [];
WBcollect = [];
NewBcollect = [];

[JunctionMeasurement_MC, JunctionErrorVariance] = GenerateMCError(4, 1000, 0.01);

A = [A_save; 1 0 0 0 0];

A1 = A(1,:);
A2 = A(2:6,:);

w1 = 1;
w2 = diag([1 1 1 1 1]);


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



% Construct pesudo demands
[m,n_d]=size(demand_MC);
[NumberofEq,NumberofX_A] = size(A);
PesudoDemand = zeros(NumberofEq,n_d);
PesudoDemand(JunctionHeadIndex,:) = demand_MC;
coefficient_MC = coefficient_MC';
PesudoDemand(PipeFlowIndex,:) = K_c'.*coefficient_MC;

% no need to plus the elevation of tanks, because we are calculating
% covariance, the constant valvue doesn't make any sense
PesudoDemand(TankHeadIndex,:) = TankMeasurement_MC;
PesudoDemand(6,:) = JunctionMeasurement_MC;

b1 = PesudoDemand(1,:);
b2 = PesudoDemand(2:6,:);

A_NEW = A1'*w1*A1 + A2'*w2*A2;
b_NEW = A1'*w1*b1 + A2'*w2*b2;

x_new = A_NEW\b_NEW;
cov(x_new')

% 
% for w_i = 1e-2:-1e-5:1e-7
% %verify=A1*abs(deterministic(1:NumberofX));
% A = [A_save; 1 0 0 0 0];
% NumberofX  = 5;
% %w = [w_i;10e10;10e10;4;10e10;4];
% % weight of demand;pipe;pump;tank;reservoir;junction
% w = [w_i; %weight of demand
%     1; % pipe
%     1; % pump
%     1; % tank
%     1; % reservoir
%     1]; % junction
% %w = [w_i;1;1;1;1;1];
% %w = [1;10e10;10e10;4;4;4];
% %w = [1;10000000000000;10000000000000;4;4;0]; 
% 
% 
% 
% 
% % assume the elecation of reservoirs has no uncertainty 
% %PesudoDemand(ReservoirHeadIndex,:) = ReservoirMeasurement_MC;
% 
% A1 = A'*diag(w)*A;
% PesudoDemand1 = A'*diag(w)*PesudoDemand;
% 
% [AnalysisMatrix1,B] = Construct_Variance_A_b_over_Pipe(NumberofX,A1,PesudoDemand1);
% % cllectA = full(AnalysisMatrix1);
% % WB = full(W*B)
% % NewB = (AnalysisMatrix1'*W*B);
% % Bcollect  = [Bcollect full(B)];
% % Wcollect  = [Wcollect full(W)];
% % WBcollect  = [WBcollect full(WB)];
% % NewBcollect  = [NewBcollect full(NewB)];
% AnalysisSolution1 = AnalysisMatrix1\B;
% 
% Covar = zeros(NumberofX,NumberofX);
% ind = 1;
% for i = 1:NumberofX
%     for j = i:NumberofX
%          Covar(i,j) = AnalysisSolution1(ind);
%          if(i~=j)
%             Covar(j,i) = Covar(i,j);
%          end
%         ind = ind+1;
%     end
% end
% 
% MC = MCSolution(1:NumberofX,:);
% cov_MC = cov(MC');
% 
% Sigma_MC = [];
% Sigma_LM1 = [];
% 
% for i = 1:NumberofX
%     if(abs(Covar(i,i)) < 1e-3)
%         Covar(i,i) = 0;
%     end
%     Sigma_LM1 = [Sigma_LM1;sqrt(Covar(i,i))]; 
%     Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
% end
% 
% % Sigma_LM1
% %Sigma_LM1 - Sigma_MC
% sigma = [sigma;Sigma_LM1(4)]
% 
% end
