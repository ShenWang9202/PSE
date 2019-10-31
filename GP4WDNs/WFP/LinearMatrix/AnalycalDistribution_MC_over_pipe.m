% Get the variances of the distribution of each demand(not zero)
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

A_save = A;
sigma = [];

Bcollect = [];
Wcollect = [];
WBcollect = [];
NewBcollect = [];

[JunctionMeasurement_MC, JunctionErrorVariance] = GenerateMCError(4, 1000, 0.01);

for w_i = 1:-0.01:0
%verify=A1*abs(deterministic(1:NumberofX));
A = [A_save; 1 0 0 0 0];
NumberofX  = 5;
%w = [w_i;10e10;10e10;4;10e10;4];
% weight of demand;pipe;pump;tank;reservoir;junction
w = [w_i; %weight of demand
    1e2; % pipe
    1e2; % pump
    1e2; % tank
    1e2; % reservoir
    1e2]; % junction
%w = [w_i;1;1;1;1;1];
%w = [1;10e10;10e10;4;4;4];
%w = [1;10000000000000;10000000000000;4;4;0]; 


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
PesudoDemand(PipeFlowIndex,:) = 0*K_c'.*coefficient_MC;

% no need to plus the elevation of tanks, because we are calculating
% covariance, the constant valvue doesn't make any sense
PesudoDemand(TankHeadIndex,:) = TankMeasurement_MC;
PesudoDemand(6,:) = JunctionMeasurement_MC;

% assume the elecation of reservoirs has no uncertainty 
%PesudoDemand(ReservoirHeadIndex,:) = ReservoirMeasurement_MC;

A1 = A'*diag(w)*A;
PesudoDemand1 = A'*diag(w)*PesudoDemand;

[AnalysisMatrix1,B] = Construct_Variance_A_b_over_Pipe(NumberofX,A1,PesudoDemand1);

AnalysisSolution1 = AnalysisMatrix1\B;

Covar = zeros(NumberofX,NumberofX);
ind = 1;
for i = 1:NumberofX
    for j = i:NumberofX
         Covar(i,j) = AnalysisSolution1(ind);
         if(i~=j)
            Covar(j,i) = Covar(i,j);
         end
        ind = ind+1;
    end
end

MC = MCSolution(1:NumberofX,:);
cov_MC = cov(MC');

Sigma_MC = [];
Sigma_LM1 = [];

for i = 1:NumberofX
    if(abs(Covar(i,i)) < 1e-3)
        Covar(i,i) = 0;
    end
    Sigma_LM1 = [Sigma_LM1;sqrt(Covar(i,i))]; 
    Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
end

% Sigma_LM1
%Sigma_LM1 - Sigma_MC
sigma = [sigma;Sigma_LM1(4)]

end

fontsize = 35;
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create plot
plot(10:-0.1:0,sigma,'LineWidth',4)
view(axes1,[-180 -90]);

xlabel('Weight of $q_{12} - q_{23} =d_2$','FontSize',fontsize+3,'interpreter','latex');
ylabel('$\sigma_{q_{23}}$','FontSize',fontsize+2,'interpreter','latex');
% Create textbox
annotation(figure1,'textbox',...
    [0.210528880031072 0.767253968253973 0.330666089798455 0.111111111111111],...
    'String',{'Probabilistic'},...
    'LineWidth',2,...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize-5,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.214367651624162 0.427571428571434 0.33278025470332 0.111111111111111],...
    'String',{'Deterministic'},...
    'LineWidth',2,...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize-5,...
    'FitBoxToText','off');

% Create arrow
annotation(figure1,'arrow',[0.351219760052722 0.353166986564298],...
    [0.771619047619048 0.538888888888889],'LineWidth',2,'HeadWidth',fontsize-5,...
    'HeadLength',fontsize-5);

% Create arrow
annotation(figure1,'arrow',[0.305182341650671 0.800383877159308],...
    [0.308333333333334 0.305555555555556],'LineWidth',2,'HeadWidth',fontsize-5,...
    'HeadLength',fontsize-5);

% Create textbox
annotation(figure1,'textbox',...
    [0.463887805175021 0.322015873015881 0.177187050871043 0.111111111111111],...
    'String','Weight',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',fontsize-5,...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);

set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
print(figure1,'sigma','-depsc2','-r300');
