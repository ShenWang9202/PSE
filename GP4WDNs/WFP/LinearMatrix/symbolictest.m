% syms x y z a
% expr = [ x*y, x*z, a*x, y^2, y*z, a*y, z^2, a*z];
% result = (x+y+z)*(a+y+z);
% expand(result)
% [cxy, txy] = coeffs(result,expr)

% load('MonteCarloDemand.mat')
% load('MonteCarloData.mat')



flowsdf = MCSolution(9:17,:);
cov_MC= cov(flowsdf');

A = [];
b = [];
% demand part
MassMatrix = ForConstructA.MassMatrix;
A = [A;MassMatrix];
DetermisticSolution=Solution(:,demandpattern);
q_pipe = DetermisticSolution(9:16,1);
% LoopMatrix = ForConstructA.LoopMatrix;
LoopMatrix = [];
LoopMatrix = [LoopMatrix;
    0	-1	0	1	-1	0	0	0 0;
    0	0	-1	0	1	0	-1	0 0;
    1 0 0 1 0 0 0 1 -1;];

Headloss_pipe_R = REstimateLinear(ForConstructb,DetermisticSolution);

K_pipe = [];
[m,~] = size(q_pipe);
for i = 1:m
    K_pipe = [K_pipe 1.852*Headloss_pipe_R(i)*abs(q_pipe(i))^(0.852)];
end


PumpEquation = IndexInVar.PumpEquation;
h0 = PumpEquation(1);
r0 = PumpEquation(2);
nu = PumpEquation(3);


q_pump= DetermisticSolution(ForConstructb.q_pump_start_index:ForConstructb.q_pump_start_index)
K_pump = nu*r0*q_pump^(nu-1);%h0/((h0/(-r0))^(1/nu));
K_estimate = [K_pipe K_pump];
[m,~] = size(LoopMatrix);
for i = 1:m
    LoopMatrix(i,:) = LoopMatrix(i,:).*K_estimate;
end

A = [A;LoopMatrix];
disp('verify')
A*[q_pipe;q_pump]

NumberofVariable = ForConstructA.PipeCount + ForConstructA.PumpCount;
vv = sym('q',[1 NumberofVariable]);
combination = [];

for i = 1:NumberofVariable
    for j = i:NumberofVariable
        combination = [combination vv(i)*vv(j)] ;
    end
end

expMatrix = [];
for i = 1:NumberofVariable
    Ai = A(i,:);
    expMatrix=[expMatrix vv*Ai'];
end

[m,n_d]=size(demand_MC);
PesudoDemand = zeros(NumberofVariable,n_d);
for i = 1:m
    PesudoDemand(i,:) = demand_MC(i,:);
end
B = [];

expression = [];
for i = 1:NumberofVariable
    for j = i:NumberofVariable
        expression = [expression expMatrix(i)*expMatrix(j)] ;
        cov_value = cov(PesudoDemand(i,:),PesudoDemand(j,:));
        B = [B;cov_value(1,2)];
    end
end

AnalysisMatrix = zeros(45,45);

[~,n] = size(expression);
for i = 1:n
    result = expression(i);
    [cxy, txy] = coeffs(result,vv);
    [~,m] = size(cxy);
    for j = 1:m
        ind = find(combination == txy(j));
        AnalysisMatrix(i,ind) = cxy(j);
    end
end
% 
% A= [-1 1 0 0 0 0 0 0 1];
% e1 =  vv * A'
% result = e1*e1;
% expand(result)
% [cxy, txy] = coeffs(result,vv)
% 
% find(combination == txy(4))
AnalysisSolution = AnalysisMatrix\B;

Covar = zeros(9,9);
ind = 1;
for i = 1:NumberofVariable
    for j = i:NumberofVariable
         Covar(i,j) = AnalysisSolution(ind);
         if(i~=j)
            Covar(j,i) = Covar(i,j);
         end
        ind = ind+1;
    end
end


Sigma_MC = [];
Sigma_LM = [];

for i = 1:NumberofVariable
    Sigma_LM = [Sigma_LM;sqrt(Covar(i,i))]; 
    Sigma_MC = [Sigma_MC;sqrt(cov_MC(i,i))]; 
end








