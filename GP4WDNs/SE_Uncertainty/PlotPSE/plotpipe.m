function plotpipe(Headloss_pipe_R)
fontsize = 30;
figure1 = figure;

% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure1);
hold(subplot1,'on');

q_pipe = 0:400;
mu = 200;
sigma = 20;
f = 1/(sqrt(2*pi)*sigma).* exp(-(q_pipe-mu).*(q_pipe-mu)/(2*sigma*sigma));

% Create plot
plot(q_pipe,f,'ZDataSource','','Parent',subplot1,'LineWidth',3,'DisplayName','PDF of $q$','Color',[0 0 0]);

% Create xlabel
xlabel({'$q$'},'Interpreter','latex');

% Create ylabel
ylabel({'Freqency'},'Interpreter','latex');

% Uncomment the following line to preserve the Y-limits of the axes
%ylim(subplot1,[0 0.005]);
box(subplot1,'on');
% Set the remaining axes properties
legend1=legend(subplot1,'show');
set(legend1,'Interpreter','latex','FontSize',fontsize-5,'Location','northoutside');
set(subplot1,'FontSize',fontsize,'XTick',[0    200]);
set(gca, 'TickLabelInterpreter', 'latex');
legend(subplot1,'show');




%%
q = 0:400;
h = Headloss_pipe_R(7).*q.^1.852;

% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure1);
hold(subplot2,'on');


q0 = 200;
k_pipe = 1.852*Headloss_pipe_R(7)*q0.^0.852;
nonlinear_headloss = Headloss_pipe_R(7)*q0^1.852;
b_pipe = nonlinear_headloss - k_pipe*q0;
h2 = k_pipe*q+b_pipe;

% Create multiple lines using matrix input to plot
plot1 = plot(q,h,'Parent',subplot2,'LineWidth',3);
hold(subplot2,'on');
plot2 = plot(q,h2,'Parent',subplot2,'LineWidth',3);
hold(subplot2,'on');
set(plot1,'DisplayName','Pipe head loss','Color',[0.87058824300766 0.490196079015732 0]);
set(plot2,'DisplayName','Linearized head loss','Color',[0 0.447058826684952 0.74117648601532]);
xlim([0 400]);
legend2 = legend(subplot2,'show');
set(legend2,'Interpreter','latex','FontSize',fontsize-5,'Location','northoutside');

% Create xlabel
xlabel({'$q$'},'Interpreter','latex');

% Create ylabel
ylabel('$\Delta h^{\mathrm{P}}$','Interpreter','latex');

box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'FontSize',fontsize,'TickLabelInterpreter','latex','XTick',...
    [0 200 400 ],'YTick',[0 100 200]);
set(gca, 'TickLabelInterpreter', 'latex');
legend(subplot2,'show');



%%
% Create subplot
%f
q_pipe = 0:400;
h_loss = Headloss_pipe_R(7).*q.^1.852;

newQ = ((h_loss)./Headloss_pipe_R(7)).^(1/1.852);
newQ2 = ((h_loss)./Headloss_pipe_R(7)).^(1/1.852-1);

mu = 200;
sigma = 20;
f = 1/(sqrt(2*pi)*sigma).* exp(-(newQ-mu).*(newQ-mu)/(2*sigma*sigma)).* (1/Headloss_pipe_R(7)).*newQ2 ;

%f2

newQ = (h_loss - b_pipe)./k_pipe;
f2 = 1/(sqrt(2*pi)*sigma).* exp(-(newQ-mu).*(newQ-mu)/(2*sigma*sigma)).* (1/k_pipe) ;

subplot3 = subplot(1,3,3,'Parent',figure1);
hold(subplot3,'on');

% Create plot
plot(h_loss,f,'DisplayName','PDF of $\Delta h^{\mathrm{P}}$','Parent',subplot3,'LineWidth',3,...
    'Color',[0.87058824300766 0.490196079015732 0]);
hold(subplot3,'on');
plot(h_loss,f2,'DisplayName','PDF of $\Delta h^{\mathrm{P}}$ (linear)','Parent',subplot3,'LineWidth',3,...
    'Color',[0 0.447058826684952 0.74117648601532]);

%f3
% q_pump = 700:1/6:900;
% h_increase = -(h0 - r.*q_pump.^nu);
% newQ = ((h_increase + h0)./r).^(1/nu);
% newQ2 = ((h_increase + h0)./r).^(1/nu-1);
% 
% mu = 800;
% sigma = 80;
% f3 = 1/(sqrt(2*pi)*sigma).* exp(-(newQ-mu).*(newQ-mu)/(2*sigma*sigma)).* (1/r*nu).*newQ2 ;
% 


% hold(subplot3,'on');
% plot(h_increase,f3,'Parent',subplot3,'LineWidth',3,...
%     'Color','r');


% Create xlabel
xlabel({'$\Delta h^{\mathrm{P}}$'},'Interpreter','latex');

% Create ylabel
ylabel({'Freqency'},'Interpreter','latex');

box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'FontSize',fontsize);
set(gca, 'TickLabelInterpreter', 'latex');
legend3=legend(subplot3,'show');
set(legend3,'Interpreter','latex','FontSize',fontsize-5,'Location','northoutside');


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 18 6])
print(figure1,'pipe','-depsc2','-r300');



end