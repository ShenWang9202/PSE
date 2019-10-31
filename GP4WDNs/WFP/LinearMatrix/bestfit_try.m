% https://www.mathworks.com/matlabcentral/fileexchange/40167-fitmethis

% x= gamrnd(5,0.5,1000,1);
% F= fitmethis(x);
% 
% mu = 685;
% sigma = 0.1*mu;
% MC_times = 1000;
% x = normrnd(mu,sigma,[1,MC_times]);
% F= fitmethis(x);
% 
% F=fitdist(x','Normal');


% assume they are normal
Variance_MonteCarlo = [];
mu_MonteCarlo = []
for i = IndexInVar.JunctionHeadIndex
    x = MCSolution(i,:);
    F=fitdist(x','Normal');
    variance = (F.sigma)^2;
    Variance_MonteCarlo = [Variance_MonteCarlo;variance];
    mu_MonteCarlo = [mu_MonteCarlo;F.mu];
end

for i = IndexInVar.ReservoirHeadIndex
    Variance_MonteCarlo = [Variance_MonteCarlo;0];
     mu_MonteCarlo = [mu_MonteCarlo;deterministic(i)];
end

for i = IndexInVar.TankHeadIndex
    Variance_MonteCarlo = [Variance_MonteCarlo;0];
     mu_MonteCarlo = [mu_MonteCarlo;deterministic(i)];
end

for i = IndexInVar.PipeFlowIndex
    x = MCSolution(i,:);
    F=fitdist(x','Normal');
    variance = (F.sigma)^2;
    Variance_MonteCarlo = [Variance_MonteCarlo;variance];
     mu_MonteCarlo = [mu_MonteCarlo;F.mu];
end

for i = IndexInVar.PumpFlowIndex
    x = MCSolution(i,:);
    F=fitdist(x','Normal');
    variance = (F.sigma)^2;
    Variance_MonteCarlo = [Variance_MonteCarlo;variance];
     mu_MonteCarlo = [mu_MonteCarlo;F.mu];
end

Variance_MonteCarlo_d = [];
mu_MonteCarlo_d=[];
for i = DemandIndex
    x = demand_MC(i,:);
    F=fitdist(x','Normal');
    variance = (F.sigma)^2;
    Variance_MonteCarlo_d = [Variance_MonteCarlo_d;variance];
     mu_MonteCarlo_d = [mu_MonteCarlo_d;F.mu];
end


% assume they are not normal
Variance_MonteCarlo = [];
mu_MonteCarlo = []
for i = IndexInVar.JunctionHeadIndex
    x = MCSolution(i,:);
    F=fitmethis(x')
end


for i = IndexInVar.PipeFlowIndex
    x = MCSolution(i,:);
    F=fitmethis(x')
end

for i = IndexInVar.PumpFlowIndex
    x = MCSolution(i,:);
    F=fitmethis(x');
end

Variance_MonteCarlo_d = [];
mu_MonteCarlo_d=[];
for i = DemandIndex
    x = demand_MC(i,:);
    F=fitdist(x','Normal');
    variance = (F.sigma)^2;
    Variance_MonteCarlo_d = [Variance_MonteCarlo_d;variance];
     mu_MonteCarlo_d = [mu_MonteCarlo_d;F.mu];
end








