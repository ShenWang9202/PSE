function PlotConfidienceInterval4SelectedFlow(Data,Data_sd,ConfidenceLevel1,ConfidenceLevel2,SelectedFlow,Variable_Symbol_Table)
% https://www.shuxuele.com/data/confidence-interval.html

% 	ConfidenceLevel
%   80%	1.282
%   85%	1.440
%   90%	1.645
%   95%	1.960
%   99%	2.576
%   99.5%	2.807
%   99.9%	3.291

[m,n] = size(SelectedFlow);
varIndex = [];
% find the index of CaredElements
for  i = 1:n
    var_index = find(strcmp(Variable_Symbol_Table,SelectedFlow{i}))
    varIndex = [varIndex;var_index];
end

figure1 = figure;
subplot1 = subplot(n,1,1,'Parent',figure1);


FlowLable = {'','','','','','','','','','P23','P34','P45','P37','P46','P67','P78','PU12'}
% prepare it for the fill function
Datax = Data(varIndex(1),:);
Data_sd_x = Data_sd(varIndex(1),:);
[~,n_Data] = size(Datax);

%% plot first one
x_ax = 1:n_Data;
X_plot = [x_ax, fliplr(x_ax)];
Y_min1 = Datax-ConfidenceLevel1.*Data_sd_x;
Y_max1 = Datax+ConfidenceLevel1.*Data_sd_x;
Y_plot1 = [Y_min1, fliplr(Y_max1)];
Y_min2 = Datax-ConfidenceLevel2.*Data_sd_x;
Y_max2 = Datax+ConfidenceLevel2.*Data_sd_x;
Y_plot2 = [Y_min2, fliplr(Y_max2)];
fontsize = 60;
% Create axes
hold('on');

fill1 = fill(X_plot, Y_plot2, 1,....
    'DisplayName','$95\%$ confidence interval',...
    'facecolor','blue', ...
    'edgecolor','none', ...
    'facealpha', 0.6);
hold on
% fill
fill2 = fill(X_plot, Y_plot1, 1,....
    'DisplayName','$80\%$ confidence interval',...
    'facecolor','red', ...
    'edgecolor','none', ...
    'facealpha', 0.6);
hold on
% plot
plot1=plot(x_ax, Datax,'DisplayName','Expectation','LineWidth',2,'Color','green',...
        'MarkerFaceColor',[1 1 1],...
        'MarkerSize',10,...
        'Marker','diamond');
xlim([1 24])
ylim([min(Datax)-10 max(Datax)+10])
% lable
ylabel(FlowLable{varIndex(1)},'FontSize',fontsize+3,'interpreter','latex');
labley= floor(min(Datax)+ max(Datax))*0.5;
set(subplot1,'FontSize',fontsize,'TickLabelInterpreter','latex','YTick',...
            [labley],'YTickLabel',{string(labley)});
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% Create legend
legend1 =legend([plot1,fill1,fill2]);
set(legend1,'Orientation','horizontal','Location','northoutside');
set(legend1,'FontSize',fontsize-10,'Interpreter','latex');

for id = 2:n
    subplot1 = subplot(n,1,id,'Parent',figure1);
    % prepare it for the fill function
    Datax = Data(varIndex(id),:);
    Data_sd_x = Data_sd(varIndex(id),:);
    [~,n_Data] = size(Datax);
    
    %% plot
    x_ax = 1:n_Data;
    X_plot = [x_ax, fliplr(x_ax)];
    Y_min1 = Datax-ConfidenceLevel1.*Data_sd_x;
    Y_max1 = Datax+ConfidenceLevel1.*Data_sd_x;
    Y_plot1 = [Y_min1, fliplr(Y_max1)];
    Y_min2 = Datax-ConfidenceLevel2.*Data_sd_x;
    Y_max2 = Datax+ConfidenceLevel2.*Data_sd_x;
    Y_plot2 = [Y_min2, fliplr(Y_max2)];
    fontsize = 60;
    % Create axes
    hold('on');
    fill1 = fill(X_plot, Y_plot2, 1,....
        'DisplayName','$95\%$ confidence interval',...
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.6);
    hold on
    % fill
    fill2 = fill(X_plot, Y_plot1, 1,....
        'DisplayName','$80\%$ confidence interval',...
        'facecolor','red', ...
        'edgecolor','none', ...
        'facealpha', 0.6);
    hold on
    % plot
    plot1 = plot(x_ax, Datax,'DisplayName','Expectation','LineWidth',2,'Color','green',...
        'MarkerFaceColor',[1 1 1],...
        'MarkerSize',10,...
        'Marker','diamond');
    % lable
    ylabel(FlowLable{varIndex(id)},'FontSize',fontsize+3,'interpreter','latex');
    xlim([1 24])
    ylim([min(Y_min2) max(Y_max2)])
    set(subplot1,'FontSize',fontsize,'TickLabelInterpreter','latex');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    labley= floor((min(Y_min2)+ max(Y_max2))*0.5);
    set(subplot1,'FontSize',fontsize,'TickLabelInterpreter','latex','YTick',...
            [labley],'YTickLabel',{string(labley)});
    if(id==n)
        ylim([min(Y_min2)-2 max(Y_max2)])
        set(subplot1,'FontSize',fontsize,'TickLabelInterpreter','latex','XTick',...
            [4 8 12 16 20 24],'XTickLabel',{'4','8','12','16','20','24'});
        xlabel('Time (h)','FontSize',fontsize+3,'interpreter','latex');
    end
end
% xlabel
%set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 24 20])
print(figure1,['ConfidentLevel','Flow'],'-depsc2','-r300');