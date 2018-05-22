

%%
if SHOW_RESULTS == 1
    if mod(iR,MOD_VAL) == 0 || iR == 1
        clear ph
        f1 = figure('visible','off');box on, hold on,
        ccolors = get(gca,'colororder');
        plot([tB(1),tB(1)],[-25 65],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        plot([t(end),t(end)],[-25 65],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        ylim([-25 65])
        text(5,55,'Training', 'FontSize',12)
        %                     text(10+tA(1),230,'Validation', 'FontSize',12)
        text(3+tA(1),55,'Validation', 'FontSize',12)
        
        if eps~=0
            ph(4) = plot(tspan,x(:,1),'-','Color',0.7*ones(1,3),'LineWidth',1); %ccolors(1,:)+[0.15 0.3 0.25]
            plot(tspan,x(:,2),'-','Color',0.7*ones(1,3),'LineWidth',1); %ccolors(2,:)+[0.15 0.3 0.25]
            plot(tspan,x(:,3),'-','Color',0.7*ones(1,3),'LineWidth',1);
            ph(1) = plot([DataTrain.t;tA],[DataTrain.x(:,1);xA(:,1)],'-','Color',ccolors(1,:),'LineWidth',0.5);
            ph(2) = plot([DataTrain.t;tA],[DataTrain.x(:,2);xA(:,2)],'-','Color',ccolors(2,:),'LineWidth',0.5);
            ph(3) = plot([DataTrain.t;tA],[DataTrain.x(:,3);xA(:,3)],'-','Color',ccolors(3,:),'LineWidth',0.5);
        else
            ph(1) = plot([DataTrain.t;tA],[DataTrain.x(:,1);xA(:,1)],'-','Color',ccolors(1,:),'LineWidth',1);
            ph(2) = plot([DataTrain.t;tA],[DataTrain.x(:,2);xA(:,2)],'-','Color',ccolors(2,:),'LineWidth',1);
            ph(2) = plot([DataTrain.t;tA],[DataTrain.x(:,3);xA(:,3)],'-','Color',ccolors(3,:),'LineWidth',1);
            ph(3) = plot(t,x(:,1),'--','Color',[0 1 0],'LineWidth',1); % Training data
            plot(t,x(:,2),'--','Color',[0 1 0],'LineWidth',1);
            plot(t,x(:,3),'--','Color',[0 1 0],'LineWidth',1);
        end
        
        ph(5) = plot(tB,xB(:,1),'-.','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
        ph(6) = plot(tB,xB(:,2),'-.','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
        ph(7) = plot(tB,xB(:,3),'-.','Color',ccolors(3,:)-[0.1 0.2 0.09],'LineWidth',2);
        grid off
        xlim([0 tv(end)])
        xlabel('Time')
        ylabel('xi')
        set(gca,'LineWidth',1, 'FontSize',14)
        set(gcf,'Position',[100 100 300 200])
        set(gcf,'PaperPositionMode','auto')
        print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_iR',num2str(iR),'_noleg.eps']);
        
        if eps~=0
            lh = legend(ph([1,4,5]),'Truth','Training',ModelName,'Location','NorthWest');
            %                         lh.Position = [lh.Position(1)+0.13,lh.Position(2)-0.2,lh.Position(3:4)];
            lh.Position = [lh.Position(1)+0.02,lh.Position(2)-0.06,lh.Position(3:4)];
        else
            lh = legend(ph([1,4,5]),'Truth','Training',ModelName,'Location','NorthWest');
            lh.Position = [lh.Position(1)+0.13,lh.Position(2)-0.2,lh.Position(3:4)];
        end
        
        print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_iR',num2str(iR),'.eps']);
        close(f1);
    end
end