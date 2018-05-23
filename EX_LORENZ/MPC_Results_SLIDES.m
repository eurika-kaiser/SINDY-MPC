% Plots for paper
% Visualizations of MPC results

clear all; close all;

figpath = '../FIGURES/LORENZ_SLIDES/'; mkdir(figpath)
datapath = '../DATA2/LORENZ/';
SystemModel = 'LORENZ';
InputSignalTypeModel = 'sphs';
InputSignalType = 'sine2';
N =10;

% load(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'.mat']),'Results')
load(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'.mat']))

%% Show results
for iM = 1:3
    xplot = Results(iM).x';
    uplot = Results(iM).u';
    figure;
    color_line3(xplot(:,1),xplot(:,2),xplot(:,3),uplot,'LineWidth',2);
    colormap(parula)
    view(30,15)
    axis([-15 5 -25 5 5 45])
    % axis tight
    axis off
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelCollection{iM},'_',InputSignalType,'_PhasePlot_Control.eps']);
end
%%
iM = 2; ModelName = 'SINDYc';

figure; hold on
xplot  = [x];
plot3(xplot(:,1),xplot(:,2),xplot(:,3),'-','Color',0.7*ones(1,3),'LineWidth',2);
view(30,15)
axis([-20 20 -25 25 5 45])
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',num2str(iM),'_Training.eps']);

xplot  = [xC];
color_line3(xplot(:,1),xplot(:,2),xplot(:,3),log(sqrt(sum((xC-xA).^2,2)))-mean(sqrt(sum((xC-xA).^2,2))),'LineWidth',2);
colormap(jet)
view(30,15)
axis([-20 20 -25 25 5 45])
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',num2str(iM),'_Validation.eps']);


figure; hold on
xplot  = [x;xC];
plot3(xplot(:,1),xplot(:,2),xplot(:,3),'-','Color',0.7*ones(1,3),'LineWidth',2);
xplot = Results(2).x';
uplot = Results(2).u';
color_line3(xplot(:,1),xplot(:,2),xplot(:,3),uplot,'LineWidth',2);
colormap(parula)
view(30,15)
axis([-20 20 -25 25 5 45])
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',num2str(iM),'_Control.eps']);

%%
figure, box on
color_line3(10+(t(2)-t(1)).*[1:length(sqrt(sum((xC-xA).^2,2)))]',sqrt(sum((xC-xA).^2,2)),zeros(size(sqrt(sum((xC-xA).^2,2)))),log(sqrt(sum((xC-xA).^2,2)))-mean(sqrt(sum((xC-xA).^2,2))),'LineWidth',2)
colormap(jet)
xlabel('t'), ylabel('E')
axis tight
set(gca,'LineWidth',1, 'FontSize',14,'ytick',10.^[-3:2],'yscale','log')
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',num2str(iM),'_ValidationError.eps']);