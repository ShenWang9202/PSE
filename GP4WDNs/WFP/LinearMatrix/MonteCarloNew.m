% Assume that each demand follow a normal distribution function, and we
% would like to see the the distribution of the states.
clc;
clear;
close all;
TestCase = 5;% 4 1 17 5 23  
tic
[inpname, acc] = findInp_acc(TestCase);
[d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution,MassEnergyMatrix4GP,MC] = PrepareA_Monte_Carlo(inpname,TestCase);
toc
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

%% plot

%call PlotMonteCarlo to plot the result (python version).

% Matlab version
%plotMenteCarlo(Solution(:,demandpattern),MCSolution,IndexInVar,Variable_Symbol_Table);


%% Analycal Distribution

%symbolictest
% only consider demand
%AnalycalDistribution_MC

%consider both demand uncertainty and pipe parameter uncertainty
% tic
% AnalycalDistribution_MC_Pipe
% toc
tic
 AnalycalDistribution_FOSM
toc
% consider overdetermined (only demand)
% AnalycalDistribution_MC_over

% consider over-determined (all uncertainty) only test for 3-node example
% make the multiple_pipe_coefficient = 0 to make sure the impact of pipe
% uncertainty is removed

%AnalycalDistribution_MC_over_pipe

%% Version 2 Kxx Kbb methods.

tic
 AnalycalDistribution_Kxx
toc


%% plot error hist for Ctown
% [size_m,size_n] = size(Sigma_error);
% relativeError = zeros(size_m,size_n);
% for i = 1:size_m
%     abserror = abs(Sigma_error);
%     if(Sigma_MC(i) < 0.2)
%         relativeError(i) = 0;
%     else
%         relativeError(i) = abserror(i)/Sigma_MC(i)*100;
%     end
%     
% end
% PlotErrorDistributionforCTown_Variance(abs(relativeError))

