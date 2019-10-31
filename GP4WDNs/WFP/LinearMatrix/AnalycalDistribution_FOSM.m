% Get the variances of the distribution of each demand(not zero)
%  consider the demand uncertainty and pipe uncertainty
load('MonteCarloData.mat')
load('MonteCarloDemand.mat')

CoeffVariance
DemandVariance
DemandIndex
PumpFlow = deterministic(IndexInVar.PumpFlowIndex);
PipeFlow = deterministic(IndexInVar.PipeFlowIndex);


%% Step 1, if all flows are positive, we can continue on, otherwise, update connection matrix at first.

% We need to update the old A matrix according our WFP solution,
% Since we assume a direction before we solve the WFP, but for two random
% variable, their variables is the sum of the individual one regardless of
% addition or substraction.

% For example,  q12 + q23= d2, if we only list equation according to this,
% we can get sig23^2 + sig23^23 = sig2^2; If the solution is q12 = 100 and
% q23=-50,d2=50, then the previous one is wrong, and the correct one should be sig23^2 = sig23^23 +
% sig2^2, which means sig23^2 - sig23^23 = sig2^2.

% Now if we update the sign according to the solution at first, things can
% always be right. For example, if we know q23=-50, it equats to sig23^2 -
% new_sig23^23 = sig2^2 where new_sig23 = 50;

% Fine the negative index in PipeFlow vector
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
% 
% verify=A1*abs(deterministic(1:NumberofX));

% profile clear
% profile on



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
PesudoDemand(PipeFlowIndex,:) = K_c'.*coefficient_MC';

% no need to plus the elevation of tanks, because we are calculating
% covariance, the constant valvue doesn't make any sense
PesudoDemand(TankHeadIndex,:) = TankMeasurement_MC';

%AnalysisSolution = AnalysisMatrix\B;
AnalysisSolution = (A'*A)\(A'*PesudoDemand);

% % tic
% [AnalysisMatrix,B] = Construct_Variance_A_b_Pipe(NumberofX,A,demand_MC,coefficient_MC,TankMeasurement_MC,K_c,ForConstructA);
% % toc
Covar = cov(AnalysisSolution');

MC = MCSolution(1:NumberofX,:);
cov_MC = cov(MC');
% 
% 
Sigma_MC = [];
Sigma_LM = [];

for i = 1:NumberofX
    Sigma_LM = [Sigma_LM;sqrt(Covar(i,i))]; 
    Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
end
Sigma_error = Sigma_MC - Sigma_LM


%profile viewer