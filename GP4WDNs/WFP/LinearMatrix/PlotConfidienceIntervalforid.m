function PlotConfidienceIntervalforid(Data,Data_sd,ConfidenceLevel1,ConfidenceLevel2,id)
% https://www.shuxuele.com/data/confidence-interval.html

% Get some random data 
% x  = linspace(0.3, pi-0.3, 10); 

% Data = sin(x) + randn(1, 10)/10; 
% Data_sd = 0.1+randn(1,10)/30; 

% 	ConfidenceLevel
%   80%	1.282
%   85%	1.440
%   90%	1.645
%   95%	1.960
%   99%	2.576
%   99.5%	2.807
%   99.9%	3.291

%ConfidenceLevel = 1.96;
% prepare it for the fill function
[~,n_Data] = size(Data);
x_ax = 1:n_Data;
X_plot = [x_ax, fliplr(x_ax)]; 
Y_plot1 = [Data-ConfidenceLevel1.*Data_sd, fliplr(Data+ConfidenceLevel1.*Data_sd)]; 
Y_plot2 = [Data-ConfidenceLevel2.*Data_sd, fliplr(Data+ConfidenceLevel2.*Data_sd)]; 

% plot a line + confidence bands 
% Create figure
figure1 = figure;
fontsize = 60;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

fill(X_plot, Y_plot2, 1,....
    'DisplayName','$80\%$ confidence interval',...
     'facecolor','blue', ... 
     'edgecolor','none', ... 
     'facealpha', 0.6); 
hold on 
% fill
fill(X_plot, Y_plot1, 1,....
    'DisplayName','$95\%$ confidence interval',...
     'facecolor','red', ... 
     'edgecolor','none', ... 
     'facealpha', 0.6); 
hold on 
 % plot
plot(x_ax, Data,'DisplayName','Prediction','LineWidth',1.5,'Color','green');


 % Create legend
%legend1 =legend(axes1,'show');
%set(legend1,'Orientation','horizontal','Location','north','FontSize',fontsize+3);
%set(legend1,'FontSize',fontsize-5,'Interpreter','latex');
% xlabel
xlabel('Time (h)','FontSize',fontsize+3,'interpreter','latex');
xlim([1 24])
set(axes1,'FontSize',fontsize,'TickLabelInterpreter','latex','XTick',...
    [4 8 12 16 20 24],'XTickLabel',{'4','8','12','16','20','24'});
% ylabel
ylabel('Demand','FontSize',fontsize+3,'interpreter','latex');
%set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 8])
['ConfidentLevel',id]
print(figure1,['ConfidentLevel',id],'-depsc2','-r300');