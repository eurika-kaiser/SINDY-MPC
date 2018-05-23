% F8 Aircraft system
% System identification: NARX / NN

clear all, close all, clc
figpath = '../FIGURES/F8/'; mkdir(figpath)
datapath = '../DATA/F8/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'F8';
ModelName = 'NARX';
Nvar = 3;
%% Generate Data
ONLY_TRAINING_LENGTH = 1;
ENSEMBLE_DATA = 1;
InputSignalType = 'sine3'; % prbs; chirp; noise; sine2; sphs; mixed; sine3
getTrainingData

%% NARX: Training
rng(2,'twister') % for reproducibility
SUBSTRACT_MEAN = 0;


% Mean-correct data
if SUBSTRACT_MEAN==0
    xmean = zeros(1,Nvar); 
else
    xmean = mean(x)';
end
xtrain = x' - repmat(xmean',[1 size(x',2)]);
utrain = u; 

% Prepare training data for NN
yt = con2seq(xtrain);
yi = con2seq(utrain);

% Neural network
stateDelays = 1;          % state delay vector
inputDelays = 1;          % input delay vector
hiddenSizes = [15 15];    % network structure (number of neurons per layer)


% Nonlinear autoregressive neural network
net = narxnet(inputDelays,stateDelays, hiddenSizes);

% Training parameters 
net.trainFcn = 'trainlm';
net.trainParam.min_grad = 1e-10;
net.trainParam.showCommandLine = 1;

% Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
[Us,Ui,Si,Ss] = preparets(net,yi,{},yt); 

% Train net with prepared training data in open-loop
tic
net = train(net,Us,Ss,Ui,Si);
toc

% Close loop for recursive prediction
netc = closeloop(net);

%% Prediction over training phase
% Prepare validation data / Get initial state from training data
[Us,Ui,Si,So] = preparets(netc,yi,{},yt); 

% Predict on validation data
predict = netc(Us,Ui,Si);
xNARX = cell2mat(predict)';
xNARX = xNARX + repmat(xmean,[size(xNARX,1) 1]);

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
xlim([0 (length(tspan)-1)*dt]), ylim([-0.8 0.8])
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

%% Prediction
% Prepare validation data
xvalid = xv' - repmat(xmean',[1 size(xv',2)]);
uvalid = uv; 

yt_valid = con2seq(xvalid);
yi_valid = con2seq(uvalid);
[Us,Ui,Si,So] = preparets(netc,yi_valid,{},yt_valid); 

% Truth
xA      = xv;
tA      = tv;

% Predict on validation data
predict = netc(Us,Ui,Si);
xB = cell2mat(predict)';
tB = tA(2:end);
xB = xB + repmat(xmean,[size(xB,1) 1]);

%% Show training and prediction
VIZ_SI_Validation

%% Save Data
Model.name = 'NARX';
Model.net = netc;
Model.xmean = xmean;
% Model.umean = umean;
% Model.Ndelay = Ndelay;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')