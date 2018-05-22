clear ph lh

figure, hold on, box on
switch color_type 
    case 'axis'
        ccolors = get(gca,'colororder');
    case 'models'
        ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
    case 'error'
        err = max(abs(xTRUTH-xModel)./xTRUTH,[],2);
        ccolors = randn(size(err));
end
symbols = {'o', 'd', 's'};
ph(1) = plot3(xTRUTH(:,1),xTRUTH(:,2),xTRUTH(:,3),'-','Color',0.4*ones(1,3),'LineWidth',2); hold on
if strcmp(color_type,'error')
    color_line3(xModel(:,1),xModel(:,2),xModel(:,3),ccolors,'LineStyle','-');
else
    ph(2) = plot3(xModel(:,1),xModel(:,2),xModel(:,3),'--','Color',ccolors(iModel,:),'LineWidth',2);
end
view(30,20),
xlabel('x'), ylabel('y'), zlabel('z')
axis([-20 20 -50 50 0 50])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 400 400])
set(gcf,'PaperPositionMode','auto')

axis tight
axis off
print('-depsc2', '-opengl','-loose', '-cmyk', [figpath,filename,'_3D_noax.eps'])

if strcmp(color_type,'error')
    ph(2) = plot3([1000,1150],[1000,1150],[1000,1150],'-k','LineWidth',2); % dummy for legend
end
lh = legend(ph,'True',ModelName);
lh.Location = 'NorthOutside';
lh.Orientation = 'Horizontal';
lh.Box = 'off'; 
print('-depsc2', '-painters','-loose', '-cmyk', [figpath,filename,'_legend.eps']);