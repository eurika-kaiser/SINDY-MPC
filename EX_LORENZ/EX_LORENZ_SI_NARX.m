% LORENZ system
% System identification: NARX

clear all, close all, clc
figpath = '../FIGURES/LORENZ/';
datapath = '../DATA/LORENZ/';
addpath('../utils');

SystemModel = 'LORENZ';

%% Generate Data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'sphs'; % prbs; chirp; noise; sine2; sphs; mixed
getTrainingData

%% NARX: Training
ModelName = 'NARX';
Nvar = 3;
rng(2,'twister')
SUBSTRACT_MEAN = 0;
xmean = mean(x)';

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
stateDelays = 1;        % state delay vector
inputDelays = 1;        % input delay vector
hiddenSizes = [10];     % network structure (number of neurons per layer)

% Nonlinear autoregressive neural network
net = narxnet(inputDelays,stateDelays, hiddenSizes);

% Training parameters %nnstart
net.trainFcn = 'trainlm';
net.trainParam.min_grad = 1e-10;
net.trainParam.showCommandLine = 1;

% Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
[Us,Ui,Si,Ss] = preparets(net,yi,{},yt); 

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
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];ccolors(2,:)-[0.1 0.2 0.09];ccolors(3,:)-[0.1 0.2 0.09]];
for i = 1:Nvar
    ph(i) = plot(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = plot(tspan(2:end),xNARX(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end
xlim([0 (length(tspan)-1)*dt]), ylim([-25 50])
xlabel('Time')
ylabel('xi')
legend(ph([1,4]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
if SUBSTRACT_MEAN == 1
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_SUBSTRACT_MEAN.eps']);
else
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);
end

%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_train'];
xModel = xNARX;
iModel = 3;
xTRUTH = x;
color_type = 'models';
VIZ_3D_MODELvsTRUTH

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
u = u'; uv = uv';
VIZ_SI_Validation

%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_valid'];
xModel = xB;
xTRUTH = xA;
iModel = 3;
color_type = 'models';
VIZ_3D_MODELvsTRUTH

%% Save Data
Model.name = 'NARX';
Model.net = netc;
Model.xmean = xmean;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')