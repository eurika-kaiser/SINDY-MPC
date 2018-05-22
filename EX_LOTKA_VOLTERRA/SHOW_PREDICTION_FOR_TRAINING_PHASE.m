switch ModelName
    case 'SINDYc'
        xData = xSINDYc;
    case 'NARX'
        xData = xNARX;
    case 'DMDc'
        xData = xDMDc;
end
       

if SHOW_RESULTS == 1
if mod(iR,MOD_VAL) == 0 || (iR == 1 && Ntrain_vec(iN)<100) || (iR == 1 && mod(Ntrain_vec(iN),100)==0)
    clear ph
    f1 = figure('visible','off');box on, hold on
    ccolors = get(gca,'colororder');
    plot([tspan(Ntrain_vec(iN)) tspan(Ntrain_vec(iN))], [0 100],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
    if eps~=0
        ph(5) = plot(tspan,x(:,1),'-','Color',0.7*ones(1,3),'LineWidth',1.5); %ccolors(1,:)+[0.15 0.4 0.25]
        ph(6) = plot(tspan,x(:,2),'-','Color',0.7*ones(1,3),'LineWidth',1.5); %ccolors(2,:)+[0.15 0.4 0.25]
        ph(1) = plot(DataTrain.tspan,DataTrain.x(:,1),'-','Color',ccolors(1,:),'LineWidth',0.5);
        ph(2) = plot(DataTrain.tspan,DataTrain.x(:,2),'-','Color',ccolors(2,:),'LineWidth',0.5);
    else
        ph(1) = plot(DataTrain.tspan,DataTrain.x(:,1),'-','Color','k','LineWidth',0.5);
        ph(2) = plot(DataTrain.tspan,DataTrain.x(:,2),'-','Color','k','LineWidth',0.5);
%         ph(5) = plot(tspan,x(:,1),'-','Color','g','LineWidth',0.5); hold on %ccolors(1,:),ccolors(2,:)
%         ph(6) = plot(tspan,x(:,2),'-','Color','g','LineWidth',0.5);
    end
    ph(3) = plot(DataTrain.tspan,xData(:,1),'-.','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',1); %2
    ph(4) = plot(DataTrain.tspan,xData(:,2),'-.','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',1); %2
    xlim([0 tv(1)-dt]), ylim([0 100])
    xlabel('Time')
    ylabel('Population size')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_iR',num2str(iR),'_noleg.eps']);
    
    close(f1);
end
end

