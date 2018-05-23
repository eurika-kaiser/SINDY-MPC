if size(uv,2)~=size(u,2) && size(u,1)~=1
    u = u';
end

%% Show training and prediction
clear ph
h = figure;
ccolors = get(gca,'colororder');
subplot(2,1,1), box on, hold on
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
text(tA(1)-30,230,'Training', 'FontSize',12)
text(5+tA(1),230,'Validation', 'FontSize',12)
% text(5+tB(1),210,'turned on', 'FontSize',12)
plot([t;tA],[x(:,1);xA(:,1)],'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'--','Color','r','LineWidth',1.5);
if exist('xB_sparse')==1
    plot(tB,xB_sparse(:,1),'-.','Color','b','LineWidth',1.5);
end
grid on
if any(xB(:)<0)
    if min(xB(:,1))<0
        ylim([min(xB(:,1)) max([xA(:,1);xB(:,1)])])
    else
        ylim([0 max([xA(:,1);xB(:,1)])])
    end
else
    if max(xA(:,1))>250
        ylim([0 max(xA(:,1))])
    else
        ylim([0 250])
    end
end
ylabel('Prey, x_1','FontSize',13)
set(gca,'FontSize',13)
subplot(2,1,2), box on, hold on
plot([tA(1) tA(1)],[0 60],':')
ph(1) = plot([t;tA],[x(:,2);xA(:,2)],'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,2),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,2),'--','Color','r','LineWidth',1.5);
if exist('xB_sparse')==1
    plot(tB,xB_sparse(:,2),'-.','Color','b','LineWidth',1.5);
end
l1=legend(ph,'Training','Validation',ModelName);
set(l1,'Location','NorthWest')
grid on
if any(xB(:)<0)
    if min(xB(:,2))<0
        ylim([min(xB(:,2)) max([xA(:,2);xB(:,2)])])
    else
        ylim([0 max([xA(:,2);xB(:,2)])])
    end
else
    if max(xA(:,2))>60
        ylim([0 max(xA(:,2))])
    else
        ylim([0 60])
    end
end
ylabel('Predator, x_2','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)
set(gca,'FontSize',13)

set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2', '-loose','-cmyk', [figpath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Validation.eps']],'-r0');

%%

clear ph
figure,box on, hold on, 
ccolors = get(gca,'colororder');
if any(xB(:)<0)
    plot([tA(1),tA(1)],[min(xB(:)) max([xA(:);xB(:)])],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
    ylim([min(xB(:)) max([xA(:);xB(:)])])
    text(tA(1)-50,max([xA(:);xB(:)])-0.15*max([xA(:);xB(:)]),'Training', 'FontSize',12)
    text(10+tA(1),max([xA(:);xB(:)])-0.15*max([xA(:);xB(:)]),'Validation', 'FontSize',12)
else
    plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
    ylim([0 250])
    text(tA(1)-50,230,'Training', 'FontSize',12)
    text(10+tA(1),230,'Validation', 'FontSize',12)
end
plot(t,x(:,1),'Color',ccolors(1,:),'LineWidth',1);
plot(t,x(:,2),'Color',ccolors(2,:),'LineWidth',1);
ph(1) = plot([t;tA],[x(:,1);xA(:,1)],'Color',ccolors(1,:),'LineWidth',1);
ph(2) = plot([t;tA],[x(:,2);xA(:,2)],'Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(tB,xB(:,1),'-','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(4) = plot(tB,xB(:,2),'-','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
if exist('xB_sparse')==1
    ph(5) = plot(tB,xB_sparse(:,1),'-.','Color',ccolors(1,:)+[0 0.2 0.2],'LineWidth',2);
    ph(6) = plot(tB,xB_sparse(:,2),'-.','Color',ccolors(2,:)+[0.1 0.2 0.09],'LineWidth',2);
end
grid off
xlim([0 200])
xlabel('Time')
ylabel('Population size')
if exist('xB_sparse')==1
    lh = legend(ph([1,3,5]),'True','DMDc',ModelName);
else
    lh = legend(ph([1,3]),'True',ModelName);
end
lh.Position = [lh.Position(1)-0.06,lh.Position(2)-0.2,lh.Position(3:4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Validation_OneFig.eps']);

%% Actuation signal
clear ph
figure,box on, hold on, 
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(tv,uv,'-k','LineWidth',1);
grid off
ylim([min([u uv])+0.05*min([u uv]) max([u uv])+0.05*max([u uv])])
xlim([0 200])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Actuation_OneFig.eps']);