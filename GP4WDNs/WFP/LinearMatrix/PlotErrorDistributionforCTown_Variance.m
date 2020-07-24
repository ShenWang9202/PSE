function PlotErrorDistributionforCTown_Variance(error)
fontsize = 50;
figure2 = figure;
% Create axes
%axes2 = axes('Parent',figure2);
h1=histogram(error,'Normalization','probability');

xlabel('$\mathrm{RE(\%)}$','FontSize',fontsize+3,'interpreter','latex');
ylabel('$\mathrm{Frequency}$','FontSize',fontsize-10,'interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize+10);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
print(figure2,'Varianceerrordistribution','-depsc2','-r300');
end



