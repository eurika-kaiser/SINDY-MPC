% LOTKA-VOLTERRA system
% System identification: NARX

clear all, close all, clc
figpath = '../../FIGURES/F8/'; mkdir(figpath)
datapath = '../../DATA/F8/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'F8';

%% Generate Data
ONLY_TRAINING_LENGTH = 1;
ENSEMBLE_DATA = 1;
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
% umean = mean(Hu(1,:));
% yt = con2seq(Hx(1:2,:)+repmat(xmean',[1 size(Hx,2)])); % Add mean again (was substracted)
% yi = con2seq(Hu(1,:)); %-umean


% prepare training data
yt = con2seq(xtrain);
yi = con2seq(utrain);

% Neural network
stateDelays = 1;     % state delay vector
inputDelays = 1;          % input delay vector
hiddenSizes = [15 15];                    % network structure (number of neurons per layer)


% Nonlinear autoregressive neural network
net = narxnet(inputDelays,stateDelays, hiddenSizes);
% net.layers{1}.transferFcn = 'logsig';
% net.layers{2}.transferFcn = 'radbas';
% net.layers{3}.transferFcn = 'purelin';

% Training parameters %nnstart
net.trainFcn = 'trainlm';%'trainbr'; %'trainlm'; trainscg
net.trainParam.min_grad = 1e-10;
net.trainParam.showCommandLine = 1;
% net.trainParam.epochs = 1000;
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
% net.performFcn = 'mse';  % Mean squared error

% Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
[Us,Ui,Si,Ss] = preparets(net,yi,{},yt); %yt

% Train net with prepared training data in open-loop
tic
net = train(net,Us,Ss,Ui,Si);
toc
% view(net)


% Plots
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotresponse(t,y)
%figure, ploterrcorr(e)
%figure, plotinerrcorr(x,e)

% Close loop for recursive prediction
netc = closeloop(net);

% % netc.trainParam.min_grad = 1e-10;
% netc.trainParam.showCommandLine = 1;
% netc.trainParam.epochs = 1000;
% 
% % Train net with prepared data in closed-loop
% [Us,Ui,Si,So] = preparets(netc,yi,{},yt); 
% netc=train(netc,Us,So,Ui);
% 
% % Predict on validation data
% predict = netc(Us,Ui,Si);
% 
% % Performance
% perfc = perform(netc,predict,So);
% view(netc);

% netc = removedelay(netc);
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

% %
% figure;
% plot(tspan(max(stateDelays):end-1),xNARX','-','LineWidth',1,'Color','k');%0.7*ones(1,3))
% grid on, hold on
 
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
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
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
% VIZ_3D_MODELvsTRUTH

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

%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_valid'];
xModel = xB;
xTRUTH = xA;
iModel = 3;
color_type = 'models';
% VIZ_3D_MODELvsTRUTH

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