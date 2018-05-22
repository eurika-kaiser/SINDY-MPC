% Visualize validation results of system identification

%% Show training and prediction
clear ph
h = figure;
ccolors = get(gca,'colororder');
subplot(3,1,1), box on, hold on
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
text(tA(1)-30,230,'Training', 'FontSize',12)
text(5+tA(1),230,'Validation', 'FontSize',12)
plot([t;tA],[x(:,1);xA(:,1)],'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'--','Color','r','LineWidth',1.5);
grid on
ylim([-25 25])
ylabel('x_1','FontSize',13)
set(gca,'FontSize',13)
subplot(3,1,2), box on, hold on
plot([tA(1) tA(1)],[0 60],':')
plot([t;tA],[x(:,2);xA(:,2)],'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,2),'k','LineWidth',1.5);
plot(tB,xB(:,2),'--','Color','r','LineWidth',1.5);
grid on
ylim([-25 25])
ylabel('x_3','FontSize',13)
set(gca,'FontSize',13)

subplot(3,1,3), box on, hold on
plot([tA(1) tA(1)],[0 60],':')
ph(1) = plot([t;tA],[x(:,3);xA(:,3)],'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,3),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,3),'--','Color','r','LineWidth',1.5);
l1=legend(ph,'Training','Validation',ModelName);
set(l1,'Location','NorthWest')
grid on
ylim([0 50])
ylabel('x_3','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)
set(gca,'FontSize',13)

set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2', '-loose','-cmyk', [figpath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation.eps']],'-r0');

%%

clear ph
figure,box on, hold on, 
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[min(xB(:)) max([xA(:);xB(:)])],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
ylim([min([xA(:);xB(:)]) max([xA(:);xB(:)])])
plot(t,x(:,1),'Color',ccolors(1,:),'LineWidth',1);
plot(t,x(:,2),'Color',ccolors(2,:),'LineWidth',1);
ph(1) = plot([t;tA],[x(:,1);xA(:,1)],'Color',ccolors(1,:),'LineWidth',1);
ph(2) = plot([t;tA],[x(:,2);xA(:,2)],'Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot([t;tA],[x(:,3);xA(:,3)],'Color',ccolors(3,:),'LineWidth',1);
ph(4) = plot(tB,xB(:,1),'-.','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(5) = plot(tB,xB(:,2),'-.','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
ph(6) = plot(tB,xB(:,3),'-.','Color',ccolors(3,:)-[0.1 0.2 0.09],'LineWidth',2);
grid off
xlim([0 tB(end)])
xlabel('Time')
ylabel('xi')
t1 = text(tA(1)-5,40,'Training', 'FontSize',12);
t2 = text(1+tA(1),40,'Validation', 'FontSize',12);
lh = legend(ph([1,4]),'True',ModelName);
lh.Position = [lh.Position(1)-0.5,lh.Position(2)-0.15,lh.Position(3:4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_OneFig.eps']);

delete(lh); delete(t1), delete(t2);
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_OneFig_noleg.eps']);
%% Actuation signal
clear ph
figure,box on, hold on, 
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[min([u;uv]) max([u;uv])],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(tv,uv,'-k','LineWidth',1);
grid off
ylim([min([u;uv])+0.05*min([u;uv]) max([u;uv])+0.05*max([u;uv])])
xlim([0 tB(end)])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Actuation_OneFig.eps']);