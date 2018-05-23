% LOTKA-VOLTERRA system
% System identification: NARX

clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

%% Generate Data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'sphs'; % prbs; chirp; noise; sine2; sphs; mixed
getTrainingData

%% NARX: Training
ModelName = 'NARX';
rng(2,'twister')
SUBSTRACT_MEAN = 0;

if SUBSTRACT_MEAN == 1
    xtrain = x' - repmat(xmean',[1 size(x',2)]);
    utrain = u; 
else
    xtrain = x';
    utrain = u;
end

% prepare training data
yt = con2seq(xtrain);
yi = con2seq(utrain);

% Neural network
stateDelays = 1;          % state delay vector
inputDelays = 1;          % input delay vector
hiddenSizes = [10];       % network structure (number of neurons per layer)

% Nonlinear autoregressive neural network
net = narxnet(inputDelays,stateDelays, hiddenSizes);

% Training parameters %nnstart
net.trainFcn = 'trainlm';%'trainbr'; %'trainlm'; trainscg
net.trainParam.min_grad = 1e-10;
net.trainParam.showCommandLine = 1;

% Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
[Us,Ui,Si,Ss] = preparets(net,yi,{},yt); %yt

% Train net with prepared training data in open-loop
tic
net = train(net,Us,Ss,Ui,Si);
toc
% view(net)

% Close loop for recursive prediction
netc = closeloop(net);

%% Prediction over training phase
% Prepare validation data / Get initial state from training data
[Us,Ui,Si,So] = preparets(netc,yi,{},yt); 

% Predict on validation data
predict = netc(Us,Ui,Si);
xNARX = cell2mat(predict)';

if SUBSTRACT_MEAN == 1
    xNARX = xNARX + repmat(xmean,[size(xNARX,1) 1]);
end
    

% Error
e = cell2mat(gsubtract(So,predict)); 
 
%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ph(1) = plot(tspan,x(:,1),'-','Color',ccolors(1,:),'LineWidth',1); hold on
ph(2) = plot(tspan,x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(tspan(2:end),xNARX(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(4) = plot(tspan(2:end),xNARX(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
xlim([0 100])
xlabel('Time')
ylabel('Population size')
legend(ph([1,3]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
if SUBSTRACT_MEAN == 1
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_SUBSTRACT_MEAN.eps']);
else
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.eps']);
end

%% Prediction
% prepare validation data
if SUBSTRACT_MEAN == 1
    xvalid = xv' - repmat(xmean',[1 size(xv',2)]);
    uvalid = uv; 
else
    xvalid = xv';
    uvalid = uv;
end

yt_valid = con2seq(xvalid);
yi_valid = con2seq(uvalid);
[Us,Ui,Si,So] = preparets(netc,yi_valid,{},yt_valid); 

% Reference
% tspan   = [100:dt:200];
xA      = xv;
tA      = tv;

% Predict on validation data
predict = netc(Us,Ui,Si);
xB = cell2mat(predict)';
tB = tA(2:end);

if SUBSTRACT_MEAN == 1
    xB = xB + repmat(xmean,[size(xB,1) 1]);
end
%% Show training and prediction
VIZ_SI_Validation

%% Save Data
Model.name = 'NARX';
Model.net = netc;
Model.xmean = xmean;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.dt = dt;
save(fullfile(datapath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')