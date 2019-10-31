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
% tic
[AnalysisMatrix,B] = Construct_Variance_A_b_Pipe(NumberofX,A,demand_MC,coefficient_MC,TankMeasurement_MC,K_c,ForConstructA);
% toc



AnalysisSolution = AnalysisMatrix\B;

Covar = zeros(NumberofX,NumberofX);
ind = 1;
for i = 1:NumberofX
    for j = i:NumberofX
         Covar(i,j) = AnalysisSolution(ind);
         if(i~=j)
            Covar(j,i) = Covar(i,j);
         end
        ind = ind+1;
    end
end

MC = MCSolution(1:NumberofX,:);
cov_MC = cov(MC');


Sigma_MC = [];
Sigma_LM = [];

for i = 1:NumberofX
    if(abs(Covar(i,i)) < 1e-5)
        Covar(i,i) = 0;
    end
    Sigma_LM = [Sigma_LM;sqrt(Covar(i,i))]; 
    Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
end

Sigma_error = Sigma_MC - Sigma_LM

%profile viewer
