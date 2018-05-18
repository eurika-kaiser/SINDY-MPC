% LOTKA-VOLTERRA system
% System identification: NARX

clear all, close all, clc
figpath = '../../FIGURES/F8/'; mkdir(figpath)
datapath = '../../DATA/F8/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'F8';

%% Generate Data
ENSEMBLE_DATA = 1;
InputSignalType = 'sine3'; % prbs; chirp; noise; sine2; sphs; mixed
getTrainingData

% sine2, Nic = 361;
% sine3, Nic = 250;

%% NARX: Training
ModelName = 'NARX';
Nvar = 3;
rng(2,'twister')

% Reshape & prepare data
M = length(tspan);
xtrain = cell(1,M);
utrain = cell(1,M);
for i = 1:M
    xtrain{i} = squeeze(x(i,:,:));
    utrain{i} = squeeze(u(i,:,:))';
end

% prepare training data
yt = xtrain;
yi = utrain;

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
net.trainParam.min_grad = 1e-6;
net.trainParam.goal = 1e-10;
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 300;
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
% iIC = find(ICinvalid==0,1,'first');
iIC = 1;
xpred = x(:,:,iIC)';
upred = u(:,1,iIC)';
yt = con2seq(xpred);
yi = con2seq(upred);
    
% Prepare validation data / Get initial state from training data
[Us,Ui,Si,So] = preparets(netc,yi,{},yt);

% Predict on validation data
predict = netc(Us,Ui,Si);
xNARX = cell2mat(predict)';

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
    ph(i) = plot(tspan,x(:,i,iIC),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = plot(tspan(2:end),xNARX(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end
% xlim([0 (length(tspan)-1)*dt]), ylim([-0.8 0.8])
xlabel('Time')
ylabel('xi')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,4]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);


%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_train'];
xModel = xNARX;
iModel = 3;
xTRUTH = x;
color_type = 'models';
% VIZ_3D_MODELvsTRUTH

%% Prediction
% prepare validation data
yt_valid = con2seq(xv(:,:,1)');
yi_valid = con2seq(uv(:,:,1)');
[Us,Ui,Si,So] = preparets(netc,yi_valid,{},yt_valid);

% Reference
% tspan   = [100:dt:200];
xA      = xv;
tA      = tv;

% Predict on validation data
predict = netc(Us,Ui,Si);
xB = cell2mat(predict)';
tB = tA(2:end);

%% Show training and prediction
% VIZ_SI_Validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];ccolors(2,:)-[0.1 0.2 0.09];ccolors(3,:)-[0.1 0.2 0.09]];
for i = 1:Nvar
    ph(i) = plot(tA,xA(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = plot(tA(2:end),xB(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end
% xlim([0 (length(tspan)-1)*dt]), ylim([-0.8 0.8])
xlabel('Time')
ylabel('xi')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,4]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation.eps']);
%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_valid'];
xModel = xB;
xTRUTH = xA;
iModel = 3;
color_type = 'models';
% VIZ_3D_MODELvsTRUTH

%% Save Data
if exist('xmean') == 0
    Model.xmean = zeros(Nvar,1);
    SUBSTRACT_MEAN = 0;
end
Model.name = 'NARX';
Model.net = netc;
% Model.umean = umean;
% Model.Ndelay = Ndelay;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')