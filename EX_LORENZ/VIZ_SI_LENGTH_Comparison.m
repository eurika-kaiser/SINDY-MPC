% LORENZ system
% Visualization for dependency on training length


clear all, close all, clc
SystemModel = 'LORENZ';

figpath = '../FIGURES/'; mkdir(figpath)
datapath = ['../DATA/EX_',SystemModel,'_Dependencies/'];
addpath('../utils');

%% Paramaters
NARXtraining = 'trainlm'; %'trainlm',trainbr
InputSignalType = 'sphs'; %sphs,sine2

%% Select case
% 1) Dependency on noise level
% Ntrain_vec = 1000; %1000;
% eta_vec = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
% Nr = 50;

% 2) Dependency on training length
Ntrain_vec = [5:15,20:5:95,100:100:1000];%,1500:500:3000];
eta_vec = 0;
Nr = 1;

%% Parameters
N_ETA = length(eta_vec);
N_LENGTHS = length(Ntrain_vec);
Nmodels = N_ETA*N_LENGTHS;
Nt = 1000;

LOG_SCALE = 1;
%% Load all models
%PostName = '_TrainLength';
ResultsALL(3) = struct('err', [],'PredHor_Ball', [], 'PredLength', [], 'RelErr', [], 'RelErrMax', [], 'Ntrain_vec', [], 'DataTrain', [], 'DataValid', [], 'Ttraining', []);

% ModelName = 'DMDc';
% datapath1 = [datapath,'DMDc/'];
% filename = fullfile(datapath1,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Nevol','_STATS.mat']);
% load(filename);
% ResultsALL(1) = Results;

ModelName = 'SINDYc';
datapath1 = [datapath,'SINDYc/'];
filename = fullfile(datapath1,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_Nevol','_STATS.mat']);
load(filename);
ResultsALL(2) = Results;

ModelName = 'NARX';
datapath1 = [datapath,'NARX/',NARXtraining,'/'];
filename = fullfile(datapath1,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_Nevol','_STATS.mat']);
load(filename);
ResultsALL(3) = Results;

%% Time series
clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 2;
ph(1) = plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.x(:,1),'-','Color',ccolors(1,:),'LineWidth',1);
ph(2) = plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.x(:,3),'-','Color',ccolors(3,:),'LineWidth',1);

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([-25 60])
xlabel('Time');
ylabel('xi')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_TS','_noleg.eps']);
l1 = legend(ph,'x1','x2','x3');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_TS','.eps']);


clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 2;
plot3(ResultsALL(iM).DataTrain.x(:,1), ResultsALL(iM).DataTrain.x(:,2), ResultsALL(iM).DataTrain.x(:,3),'-k','LineWidth',1);
axis equal, axis tight, %axis off
set(gca,'xtick', [], 'ytick', [], 'ztick', [])
view(40,30)
% xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 100])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 300])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_PhasePlot','.eps']);

clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 2;
plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.u,'-k','LineWidth',1);
%axis equal,% axis tight, %axis off
set(gca,'xtick', [], 'ytick', [])
% xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 100])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 50])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_Actuation','.eps']);

%% Mean squared error
clear ph
% Prediction over training/validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
    plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,2),'--','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];
xlim([0 max(ResultsALL(iM).DataTrain.t)]),
%ylim([0 12000])
xlabel('Training Time');
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE','_ZOOM.eps']);


% Prediction over training
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]),
ylim([0 4000])
xlabel('Training Time');
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','_ZOOM.eps']);

% Prediction over validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];

for iM = 2:3
    data2plot = ResultsALL(iM).err(:,2);
    if any(isnan(data2plot)==1)
        data2plot(find(isnan(ResultsALL(iM).err(:,2)))) = 10^6;
    end
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), data2plot,'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(2) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]),

% if log scale
if LOG_SCALE
    ylim([10^2 1*10^6])
    set(gca,'yscale','log')
    set(gca,'ytick',[10.^[2:2:7]])
end

xlabel('Training Time');
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','_ZOOM.eps']);

%% Relative error
clear ph
% Prediction over training/validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
    plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'--','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Time');
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','_ZOOM.eps']);

% Prediction over training
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];
xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Time');
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','_ZOOM.eps']);

% Prediction over validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    data2plot = ResultsALL(iM).RelErr(:,2);
    if any(isnan(ResultsALL(iM).err(:,2))==1)
        data2plot(find(isnan(ResultsALL(iM).err(:,2)))) = 10^6;
    end
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), data2plot,'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];
xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])


% if log scale
if LOG_SCALE
    ylim([10^0 10^6])
    set(gca,'yscale','log')
    set(gca,'ytick',[10.^[0:2:6]])
end

xlabel('Training Time');
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','.eps']);

delete(l1)
xlim([0 10]), %ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','_ZOOM.eps']);

%% Training time
clear ph
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).Ttraining,'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Length');
ylabel('Training Time [s]')

% if log scale
set(gca,'yscale','log')
ylim([10^-3 3*10^1])
set(gca,'ytick',[10^-3 10^-1 10^1])

set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_TrainTime','_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthWest')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_TrainTime','.eps']);

%% Prediction horizon
dt = 0.01;
iCase = 2;
if iCase == 1
    PostName1 = ['_train'];
elseif iCase == 2
    PostName1 = ['_valid'];
end
clear ph
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).PredLength(:,iCase)*dt,'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Length');
ylabel('Prediction horizon')


% if log scale
set(gca,'yscale','log')
ylim([10^-3 10^1])
set(gca,'ytick',[10^-3 10^-1 2.5 10^1])
set(gca,'xtick',[0.38 5 10])

set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_PredHorizon',PostName1,'_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','SouthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_PredHorizon',PostName1,'.eps']);

%%
clear ph
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 2:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).PredHor_Ball(:,iCase)*dt,'-','Color',ccolors(iM,:),'LineWidth',1);
end
ph(1) = [];

xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Length');
ylabel('Prediction horizon')

% if log scale
set(gca,'yscale','log')
ylim([10^-2 10^1])
set(gca,'ytick',[10.^[-2:1:1]])
set(gca,'xtick',[0 5 10])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_PredHorizonBall',PostName1,'_noleg.eps']);
% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','SouthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_TrainingLength_PredHorizonBall',PostName1,'.eps']);