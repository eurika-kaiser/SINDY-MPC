clear all, close all, clc
figpath = '../FIGURES/'; mkdir(figpath)
datapath = '../DATA/EX_LOTKA_Dependencies/';
addpath('../utils');

%% Paramaters
NARXtraining = 'trainlm'; %'trainlm',trainbr
InputSignalType = 'sphs'; %sphs,sine2

% Noise
% Ntrain_vec = 1000; %1000;
% eta_vec = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
% Nr = 50;

% Training length
Ntrain_vec = [5:15,20:5:95,100:100:1000];%,1500:500:3000];
eta_vec = 0;
Nr = 1;

N_ETA = length(eta_vec);
N_LENGTHS = length(Ntrain_vec);
Nmodels = N_ETA*N_LENGTHS;

Nt = 1000; % length(tv);

LOG_SCALE = 1;
%% Load all models
ResultsALL(3) = struct('err', [], 'RelErr', [], 'Ntrain_vec', [], 'DataTrain', [], 'DataValid', [], 'Ttraining', []);

ModelName = 'DMDc';
datapath1 = [datapath,'DMDc/'];
filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Nevol','_STATS.mat']);
load(filename);
ResultsALL(1) = Results;

ModelName = 'SINDYc';
datapath1 = [datapath,'SINDYc/'];
filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Nevol','_STATS.mat']);
load(filename);
ResultsALL(2) = Results;

ModelName = 'NARX';
datapath1 = [datapath,'NARX/',NARXtraining,'/'];
filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_Nevol','_STATS.mat']);
load(filename);
ResultsALL(3) = Results;

%% Time series
clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 1;
ph(1) = plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.x(:,1),'-','Color',ccolors(1,:),'LineWidth',1);
ph(2) = plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 100])
xlabel('Time'); 
ylabel('Population size')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_TS','_noleg.eps']);
l1 = legend(ph,'Predator','Prey');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_TS','.eps']);


clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 1;
plot(ResultsALL(iM).DataTrain.x(:,1), ResultsALL(iM).DataTrain.x(:,2),'-k','LineWidth',1);
axis equal, axis tight, %axis off
set(gca,'xtick', [], 'ytick', [])
% xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 100])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_PhasePlot','.eps']);

clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 1;
plot(ResultsALL(iM).DataTrain.t, ResultsALL(iM).DataTrain.u,'-k','LineWidth',1);
axis equal, axis tight, %axis off
set(gca,'xtick', [], 'ytick', [])
% xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 100])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 50])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_Actuation','.eps']);


%% Mean squared error
clear ph
% Prediction over training/validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
    plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,2),'--','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 12000])
xlabel('Training Time'); 
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE','_ZOOM.eps']);


% Prediction over training
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).err(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 12000])
xlabel('Training Time'); 
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_train','_ZOOM.eps']);

% Prediction over validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];

for iM = 1:3
    data2plot = ResultsALL(iM).err(:,2);
    if any(isnan(data2plot)==1)
        data2plot(find(isnan(ResultsALL(iM).err(:,2)))) = 10^6;
    end
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), data2plot,'-','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), 

% if log scale
if LOG_SCALE
    ylim([10^-1 5*10^4])
    set(gca,'yscale','log')
    set(gca,'ytick',[10.^[-1:2:4]])
end

xlabel('Training Time'); 
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2500])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_MSE_valid','_ZOOM.eps']);

%% Relative error
clear ph
% Prediction over training/validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
    plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'--','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 6])
xlabel('Training Time'); 
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr','_ZOOM.eps']);

% Prediction over training
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).RelErr(:,1),'-','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 6])
xlabel('Training Time'); 
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_train','_ZOOM.eps']);

% Prediction over validation
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    data2plot = ResultsALL(iM).RelErr(:,2);
    if any(isnan(ResultsALL(iM).err(:,2))==1)
        data2plot(find(isnan(ResultsALL(iM).err(:,2)))) = 10^6;
    end
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), data2plot,'-','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), ylim([0 6])


% if log scale
if LOG_SCALE
    ylim([5*10^-3 10^1])
    set(gca,'yscale','log')
    set(gca,'ytick',[10.^[-3:1:1]])
end

xlabel('Training Time'); 
ylabel('Avg. Rel. Error')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','.eps']);

delete(l1)
xlim([0 30]), ylim([0 2])
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_RelErr_valid','_ZOOM.eps']);

%% Training time
clear ph
figure, hold on, box on
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
for iM = 1:3
    ph(iM) = plot(ResultsALL(iM).DataTrain.t(ResultsALL(iM).Ntrain_vec), ResultsALL(iM).Ttraining,'-','Color',ccolors(iM,:),'LineWidth',1);
end

xlim([0 max(ResultsALL(iM).DataTrain.t)]), %ylim([0 6])
xlabel('Training Length'); 
ylabel('Training Time [s]')

% if log scale
set(gca,'yscale','log')
ylim([10^-3 3*10^0])
set(gca,'ytick',[10^-3 10^-2 10^-1 10^0])

set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_TrainTime','_noleg.eps']);
l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthWest')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_TrainingLength_TrainTime','.eps']);

