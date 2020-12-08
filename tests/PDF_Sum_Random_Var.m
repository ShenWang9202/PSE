% please install CUPID toolbox before executing the code below.

% https://github.com/milleratotago/Cupid/commits/master/Cupid.mltbx

%% install cupid toolbox
toolboxFile = 'Cupid.mltbx';
agreeToLicense = true;
installedToolbox = matlab.addons.toolbox.installToolbox(toolboxFile,agreeToLicense);

 
%%
fontsize = 30;
%% Example 1
pd1 = makedist('Normal','mu',61.7,'sigma',4.48);
pd2 = makedist('Normal','mu',75,'sigma',7.5);

fig1 = figure
x = 40:1:80;
y = pdf(pd1,x);
plot(x,y,'LineWidth',2.5)
hold on

x = 55:1:95;
y = pdf(pd2,x);
plot(x,y,'LineWidth',2.5);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
legend1 = legend('PDF(X)','PDF(Y)');
set(legend1,'Interpreter','latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'example1_XY','-depsc2','-r300');
% MIX
mixtureDistributions(pd1,pd2)

%% Example 2
pd1 = makedist('Normal','mu',61.7,'sigma',4.48);
pd2 = makedist('Uniform','lower',50,'upper',100);

fig1 = figure
x = 40:1:80;
y = pdf(pd1,x);
plot(x,y,'LineWidth',2.5);
hold on

x = 40:1:110;
y = pdf(pd2,x);
plot(x,y,'LineWidth',2.5)

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
legend1 = legend('PDF(X)','PDF(Y)');
set(legend1,'Interpreter','latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'example2_XY','-depsc2','-r300');

% MIX
mixtureDistributions(pd1,pd2)


%% Example 3

data1 = normrnd(61.7,4.48,1000,1);
data2 = normrnd(75,7.5,1000,1);
% plot data set
fig1 = figure;
histogram(data1)
hold on
legend1 = legend('x dataset');
set(legend1,'Interpreter','latex','FontSize',fontsize);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'example3_X_data','-depsc2','-r300');

fig2 = figure;
histogram(data2)
hold on
legend2 = legend('y dataset');
set(legend2,'Interpreter','latex','FontSize',fontsize);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig2,'example3_Y_data','-depsc2','-r300');

% KDE
pd1 = fitdist(data1,'normal');
pd2 = fitdist(data2,'normal');

fig3 = figure
x = 40:1:80;
y = pdf(pd1,x);
plot(x,y,'LineWidth',2.5)
hold on

x = 55:1:95;
y = pdf(pd2,x);
plot(x,y,'LineWidth',2.5);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
legend1 = legend('PDF(X)','PDF(Y)');
set(legend1,'Interpreter','latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig3,'example3_XY','-depsc2','-r300');
% MIX
mixtureDistributions(pd1,pd2)

%% Example 4


data1 = normrnd(61.7,4.48,1000,1);
data2 = unifrnd(50,100,[1000,1]); 

% plot data set
fig1 = figure;
histogram(data1)
hold on
legend1 = legend('x dataset');
set(legend1,'Interpreter','latex','FontSize',fontsize);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'example4_X_data','-depsc2','-r300');

fig2 = figure;
histogram(data2)
hold on
legend2 = legend('y dataset');
set(legend2,'Interpreter','latex','FontSize',fontsize);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig2,'example4_Y_data','-depsc2','-r300');

% KDE
pd1 = fitdist(data1,'Kernel');
pd2 = fitdist(data2,'Kernel');


fig1 = figure
x = 40:1:80;
y = pdf(pd1,x);
plot(x,y,'LineWidth',2.5)
hold on

x = 40:1:110;
y = pdf(pd2,x);
plot(x,y,'LineWidth',2.5);

set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
legend1 = legend('PDF(X)','PDF(Y)');
set(legend1,'Interpreter','latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'example4_XY','-depsc2','-r300');

% MIX
mixtureDistributions(pd1,pd2)

%% The AN EXAMPLLE CODE to show the core algorithm of CUPID: 

x = [randn(30,1); 5+randn(30,1)];
[f,xi] = ksdensity(x,'function','cdf');
F=[spline(f,xi,0.05) spline(f,xi,0.995)]

% x data
x = [randn(100,1); 5+randn(100,1)];
% plot histogram of x data
fig1 = figure;
histogram(x,30)
hold on
legend1 = legend('x dataset');
set(legend1,'Interpreter','latex','FontSize',fontsize);
set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig1,'xdata','-depsc2','-r300');


% use KDE to estimate the pdf of x
pd = fitdist(x,'Kernel');
x_values = -5:1:10;
y = pdf(pd,x_values);

% plot the estimated pdf of x
fig2 = figure;
plot(x_values,y,'LineWidth',2.5)
hold on
legend2 = legend('PDF(X)');
set(legend2,'Interpreter','latex','FontSize',fontsize);
set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(fig2,'Xpdf','-depsc2','-r300');

% we can find the shortest confidence interval by solve an optimization
% problem

%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that a solution is found in 
%   the range LB <= X <= UB. Use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%   set UB(i) = Inf if X(i) is unbounded above.
t = fmincon(@(x)icdf(pd,x)-icdf(pd,x-0.95),0.95,[],[],[],[],[0.95],[1])   %0.95 is the confidence level you set
%t= 0.9744
[icdf(pd,t-0.95), icdf(pd,t)]  %the shortest confidence interval
set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);

