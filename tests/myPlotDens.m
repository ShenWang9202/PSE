% Update the PlotDens function in dGeneric.m (from CUPID toolbox) to plot the results


  function PlotDens(obj,varargin)  % NWJEFF: More work needed:
            % Extend to bounds if the tails are not too stretched.
            % Parameters to specify X range or CDF(X) range
            % Parameters to select out only PDF or CDF or Hazard
            % Nice numbers for PDF ticks
            fontsize = 30;
            [SkipFig, varargin] = ExtractNamei({'NoFig','SkipFig'},varargin);
            nvarargin = numel(varargin);
            nticks = 5;
            x = XsToPlot(obj);
            pdfy = PDF(obj,x);
            cdfy = CDF(obj,x);
%             if ~SkipFig
%                 figure;
%             end
            fig1 = figure;
            [AX,H1,H2] = plotyy(x,pdfy,x,cdfy);
            set(AX,'FontSize',fontsize)
            
            set(H1,'LineWidth',2.5);
            set(H2,'LineWidth',2.5);
            
            ylim(AX(2),[-0.02 1.02]);
            xlabel('X');
            title(obj.StringName,'interpreter', 'latex','FontSize',fontsize-10)
            %title('Convolution(Normal(61.7,4.62),Normal(75.3,7.73))','interpreter', 'latex','FontSize',fontsize-10)
            set(get(AX(1),'Ylabel'),'String','PDF(Z)');
            set(get(AX(2),'Ylabel'),'String','CDF(Z)');
            minpdfy = 0;
            maxpdfy = max(pdfy);
            set(AX(1),'YTick',minpdfy:(maxpdfy-minpdfy)/nticks:maxpdfy);
            set(AX(2),'YTick',0:0.2:1);
            if obj.DistType == 'd'
                set(H1,'linestyle','none','marker','o');
                set(H2,'linestyle','none','marker','*');
            end
            pos = get(AX,'Position');
            pos = pos{2};
            pos(1) = pos(1) + 0.2;
            pos(3) = pos(3) - 0.2;
            set(AX,'Position',pos)
            set(gca, 'TickLabelInterpreter', 'latex','FontSize',fontsize);
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
            print(fig1,'example4','-depsc2','-r300');
        end
