clc;
clear;
close all;
base = 1.001;
TestCase = 23;
[inpname, acc] = findInp_acc(TestCase);
[d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution] = PrepareA(inpname,TestCase);
%load('forctown')
BarSolution=[];
BarSolution=[BarSolution Solution(:,1)];

%% sufficient scenario.
Demand_known=InitialParameter.Demand_known;
[m,n] = size(Demand_known);
Error_All = cell(n,1);
Relative_Error_All = [];
IterateError_All = [];
XSolution = [];

M2FT = InitialParameter.M2FT;
LPS2GMP = InitialParameter.LPS2GMP;
NumberofX= ForConstructA.NumofX;

for i = 1:1
    %% initial conditions
    X0 = zeros(NumberofX,1);
    %     X0(JunctionHeadIndexInVar) = Solution(JunctionHeadIndexInVar,i);
    X0(IndexInVar.ReservoirHeadIndex) = Solution(IndexInVar.ReservoirHeadIndex,i)*M2FT;
    X0(IndexInVar.TankHeadIndex) = Solution(IndexInVar.TankHeadIndex,i)*M2FT;
    X0(IndexInVar.PumpSpeedIndex) = Solution(IndexInVar.PumpSpeedIndex,i);
    %     X0(IndexInVar.PumpFlowIndex) = Solution(IndexInVar.PumpFlowIndex,i)*LPS2GMP;
    %     X0(IndexInVar.PumpFlowIndex) = Solution(IndexInVar.PumpFlowIndex,i)*LPS2GMP;
    %     X0(IndexInVar.ValveFlowIndex) = Solution(IndexInVar.ValveFlowIndex,i)*LPS2GMP;
    X0(IndexInVar.PumpFlowIndex) = 100;%InitialParameter.PipeFLowAverage(IndexInVar.PumpFlowIndex)*LPS2GMP;
    X0(IndexInVar.PipeFlowIndex) = 0;%InitialParameter.PipeFLowAverage(IndexInVar.PipeFlowIndex)*LPS2GMP*rand(1,1);%(-1+rand(1,1)*(2));
    X0(IndexInVar.ValveFlowIndex) = InitialParameter.PipeFLowAverage(IndexInVar.ValveFlowIndex)*LPS2GMP;
    % make sure the flow of pump is greater than 0; can not be negative
    %X0(IndexInVar.PumpFlowIndex) = 1;% make sure the flow of pump is greater than 0; can not be negative
    X0
    Wsolution = [];
    Wsolution = [Wsolution;X0];
    C_estimate = [];
    Error =[];
    demand = Demand_known(:,i)';
    index = 1;
    IterateError = 10;
    % measure the head at 'J3' and 'J5'
    MeasuredHeadLable = { 'J2' };
    % and the measurements are
    MeasuredValue = [910];
    MeasuredIndex= findIndexByLabel(MeasuredHeadLable,Variable_Symbol_Table); 
    tic
   % while (IterateError > 0.1 & index < 200)
    while (IterateError > 0.01)
        DEMAND=demand;
        [A,b] = constructA_b(ForConstructA,ForConstructb,DEMAND,X0,base);
        %W = inv(A'*A)*A'*b;
        %W = inv(A)*b;
        W = A\b;
        W= [W;Solution(IndexInVar.PumpSpeedIndex,i)];
        Wsolution = [Wsolution W];
        X0 = W;
        if(mod(index,4)==0)
            tendency = Wsolution(:,index) - Wsolution(:,index-2);
            X0 = X0 + acc * tendency(:,1);
        end
        IterateError = norm(Wsolution(:,end)-Wsolution(:,end-1));
        Error = [Error;IterateError];
        index = index + 1;
    end
    % after solving sufficient measurement scenerio, now let's consider
    % over-determined scenerio.
    [A_over,b_over] = Measured(A,b,MeasuredIndex,MeasuredValue);
    w = [10e2 0 0 0 0 0;
            0  1 0 0 0 0;
            0  0 10e3 0 0 0;
            0  0 0 1 0 0;
            0  0 0 0 1 0;
            0  0 0 0 0 10e2];
    SE = (A_over'*w*A_over)\(A_over'*w*b_over);
end
toc
BarSolution=[BarSolution Wsolution(:,end)];
i=1;
Solution1(IndexInVar.JunctionHeadIndex) = Solution(IndexInVar.JunctionHeadIndex,i)*M2FT;
Solution1(IndexInVar.ReservoirHeadIndex) = Solution(IndexInVar.ReservoirHeadIndex,i)*M2FT;
Solution1(IndexInVar.TankHeadIndex) = Solution(IndexInVar.TankHeadIndex,i)*M2FT;
Solution1(IndexInVar.PumpSpeedIndex) = Solution(IndexInVar.PumpSpeedIndex,i);
Solution1(IndexInVar.PumpFlowIndex) = Solution(IndexInVar.PumpFlowIndex,i)*LPS2GMP;
Solution1(IndexInVar.PipeFlowIndex) =  Solution(IndexInVar.PipeFlowIndex,i)*LPS2GMP;
Solution1(IndexInVar.ValveFlowIndex) =  Solution(IndexInVar.ValveFlowIndex,i)*LPS2GMP;


ErrorWithEpanet = [];
[~,n] = size(Wsolution)
for i= 1:n
    ErrorWithEpanet = [ErrorWithEpanet;norm(Wsolution(:,i)-Solution1')];
end


figure;
plot(C_estimate','DisplayName','C_estimate')

% plot convergence and C_P
fontsize = 40;
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

pl1=plot(log10(ErrorWithEpanet),'LineWidth',5);
% Create xlabel
xlabel('Iteration','Interpreter','latex');

% Create ylabel
ylabel('Error','Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',40,'TickLabelInterpreter','latex','YTick',[-1 1 3 4]);
% Create legend
legend1 = legend([pl1],'$\log_{10}(||\boldmath \xi_{\mathrm{SE}}-\boldmath \xi_{\mathrm{EPANET}}||)$');
set(legend1,'Interpreter','latex','FontSize',50,'FontName','Helvetica Neue',...
    'Location','best');
set(gca,'FontSize',40,'TickLabelInterpreter','latex','YTick',[-1  4]);
hold off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 18 6])
print(figure1,'errorEPANET','-depsc2','-r300');

