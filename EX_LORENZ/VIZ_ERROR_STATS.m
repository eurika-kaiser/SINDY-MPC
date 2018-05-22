
errData = data2plot;

clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
symbols = {'o', 'd', 's'};
xx = 1:0.1:11;
yy = zeros(3,length(xx));
for iM = 2:3
    if iM == 2
        xm = zeros(1,N_ETA);
        for k = 1:N_ETA 
            TF = isnan(errData(1:Nr,k,iM));
            xm(k) = median(errData(TF==0,k,iM),1);
        end
    else
        xm = median(errData(1:Nr,:,iM),1);
    end
    yy(iM,:) = spline(1:11,xm,xx); 
end
figure; hold on
if shaded_region == 1
    fillh1 = fill([0.1 11.9 11.9 0.1], [10 10 ylimregion.*ones(1,2)],0.9*ones(1,3));
    fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5; 
end
ph(2)=plot(xx(1:end),yy(2,1:end),'-g','LineWidth',1);
ph(3)=plot(xx(1:end),yy(3,1:end),'-','Color',[0.7,0.7,1],'LineWidth',1);
ph(1) = [];
for iM = 2:3
    boxplot(errData(:,:,iM),'PlotStyle','compact','Colors',ccolors(iM,:), 'Symbol', symbols{iM}, ... 
        'Widths',0.1); hold on
    delete(findobj(gca,'Type','text'))
end
set(gca,'xtick',[1:2:11],'xticklabel',num2str(eta_vec(1:2:end)'),'Position',[0.2 0.22 0.75 0.73])
xlim([0 12]), ylim(yaxlim)


if exist('LOG_SCALE')
    if LOG_SCALE == 1
        set(gca,'YScale','log')
        ylim([10^-3 yaxlim(2)])
        set(gca,'ytick',[10.^[-3:1:1]])
    end
end

if exist('YSCALE_SET')
    if YSCALE_SET == 1
        ylim([yaxlim])
        if length([log10(yaxlim(1)):2:log10(yaxlim(2))])<3
            set(gca,'ytick',[10.^[log10(yaxlim(1)):1:log10(yaxlim(2))]])
        else
            set(gca,'ytick',[10.^[log10(yaxlim(1)):2:log10(yaxlim(2))]])
        end
    end
end

xt = xlabel('Eta'); xt.Position = [115 -20 -0.1];
ylabel(ytext)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'), 
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'_',PostName,'_noleg.eps']);

l1 = legend(ph,'SINDYc','NARX');
if exist('LegLoc')==1
    set(l1,'Location',LegLoc)
else
    set(l1,'Location','NorthWest')
    l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
end
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'_',PostName,'.eps']);
