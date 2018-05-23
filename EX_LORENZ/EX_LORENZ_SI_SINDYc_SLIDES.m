% LORENZ system
% Visualizations of results

clear all, close all, clc
figpath = '../FIGURES/LORENZ_SLIDES/'; mkdir(figpath)
addpath('../utils');


SystemModel = 'LORENZ';
ModelName = 'SINDYc';
%% Generate Data %{'sine2', 'chirp','prbs', 'sphs'}
InputSignalType = 'sphs'; %prbs; chirp; noise; sine2; sphs; mixed
ONLY_TRAINING_LENGTH = 1;
getTrainingData

%% SINDYc
u = u';

% Parameters
ModelName = 'SINDYc'; iModel = 2;
Nvar = 3;
polyorder = 4;
usesine = 0;
lambda = 0.1;
eps = 0;
xclean = x;

%Add noise
xclean = x;
eps = 0.05*std(x(:,1));
x = x + eps*randn(size(x));

%%
figure;
plot3(x(:,1),x(:,2),x(:,3),'-k','LineWidth',1);
color_line3(xclean(:,1),xclean(:,2),xclean(:,3),u,'LineWidth',2);
colormap(parula)
view(30,15)
axis tight
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_PhasePlot.eps']);

%%
polyorder = 2;
usesine = 0;
lambda = 0.5;
trainSINDYc

return
%%
ii = 1;
figure,
plot(t(1:size(xaug(:,ii),1)),xaug(:,ii),'-r','LineWidth',2);
hold on
plot(t(1:size(xaug(:,ii),1)),xclean(50:end-50-1,ii),'-k','LineWidth',1)
ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
xlim([1 4])
xlabel('t'), ylabel(['x_',num2str(ii)])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 200 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_eta',num2str(eps),'.eps']);

return
%%
ccol = {'b','r','g','k'};
tshort = t(1:size(dx(:,1),1));

for ii = 1:3
    
    xabs = max(abs([min(xaug(:,ii)), max(xaug(:,ii))]));
    figure;
    plot(tshort,xaug(:,ii),'-','Color',ccol{ii},'LineWidth',2)
    ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
    xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
    % daspect([1 30 1])
    xlabel('t'), ylabel(['x_',num2str(ii)])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 100])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_eta',num2str(eps),'.eps']);
    
    axis off
    axis tight
    % daspect([1 50 1])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 100])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_AxisOFF','_eta',num2str(eps),'.eps']);
    
    % dx
    xabs = max(abs([min(dx(:,ii)), max(dx(:,ii))]));
    
    figure;
    plot(tshort,dx(:,ii),'-','Color',ccol{ii},'LineWidth',2)
    ylim([min(dx(:,ii))-0.5*xabs  max(dx(:,ii))+0.5*xabs])
    xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
    % daspect([1 30 1])
    xlabel('t'), ylabel(['dx_',num2str(ii)])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 100])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_dx',num2str(ii),'_eta',num2str(eps),'.eps']);
    
    axis off
    axis tight
    % daspect([1 50 1])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 100])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_dx',num2str(ii),'_AxisOFF','_eta',num2str(eps),'.eps']);
end

ii = 4; % u
xabs = max(abs([min(xaug(:,ii)), max(xaug(:,ii))]));
figure;
plot(tshort,xaug(:,ii),'-','Color',ccol{ii},'LineWidth',2)
ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
% daspect([1 30 1])
xlabel('t'), ylabel(['u'])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_u','_eta',num2str(eps),'.eps']);

axis off
axis tight
% daspect([1 50 1])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_u','_AxisOFF','_eta',num2str(eps),'.eps']);

