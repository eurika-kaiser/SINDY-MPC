% HIV model
% System identification using NARX model
 
clear all, close all, clc
figpath = '../../FIGURES/HIV/';
datapath = '../../DATA/HIV/';
addpath('../utils');

SystemModel = 'HIV';
ModelName = 'NARX';

%% Generate Data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'prbs'; % prbs; chirp; noise; sine2; sphs; mixed
DATA_ENSEMBLE = 1;
if DATA_ENSEMBLE==1
    getTrainingData_Ensemble
else
    getTrainingData
end
return
%%
% figure, hold on
% for i = 1:16
% plot(xensemble(:,1,i))
% % plot(u(:,1,i))
% end
%%
% i=2
% figure,hold on
% plot(x(:,:,i),'-k','LineWidth',2)
% plot(1000*u(:,1,i),'--r','LineWidth',2)
% xforc = x;
% uforc = u;
%%
% figure,hold on
% subplot(5,1,1), hold on
% plot(squeeze(x(:,1,:)),'-k','LineWidth',2)
% plot(squeeze(xforc(:,1,:)),'--r','LineWidth',2)
% subplot(5,1,2), hold on
% plot(squeeze(x(:,2,:)),'-k','LineWidth',2)
% plot(squeeze(xforc(:,2,:)),'--r','LineWidth',2)
% subplot(5,1,3), hold on
% plot(squeeze(x(:,3,:)),'-k','LineWidth',2)
% plot(squeeze(xforc(:,3,:)),'--r','LineWidth',2)
% subplot(5,1,4), hold on
% plot(squeeze(x(:,4,:)),'-k','LineWidth',2)
% plot(squeeze(xforc(:,4,:)),'--r','LineWidth',2)
% subplot(5,1,5), hold on
% plot(squeeze(x(:,5,:)),'-k','LineWidth',2)
% plot(squeeze(xforc(:,5,:)),'--r','LineWidth',2)
x = abs(x);
%%
    %     xtrain = x;
    %     for i = 1:size(x,2)
    %         xtrain(:,i) = xtrain(:,i)-repmat(xmean(i),[size(x,1) 1]);
    %     end
Nvar = size(x,2);    
TRANSFORM_LOG = 1;
if TRANSFORM_LOG == 1
    xtrain = log(x);
else
    xtrain = x;
end
SUBSTRACT_MEAN = 0;
if SUBSTRACT_MEAN ==1
    xmean = xSTEADY1';
else
    xmean = zeros(Nvar,1);
end
NORMALIZE = 0;
xnorm = ones(Nvar,1);

Nvar = size(x,2);
if DATA_ENSEMBLE == 1
%     rng(10,'twister')
    rng(10,'twister')
    
    Nmodel = 5;
    % Reshape & prepare data
    Nshift = 0;
    M = size(xtrain,1)-Nshift; %(size(xtrain,1)-1)/2; %size(xtrain,1);
    yt = cell(1,M);
    yi = cell(1,M);
    for i = Nshift+1:size(xtrain,1)%M
        yt{i-Nshift} = squeeze(xtrain(i,1:Nmodel,:)) - repmat(xmean(1:Nmodel),[1 Nic]);
        yi{i-Nshift} = squeeze(u(i,:,:))';
    end
    
    % Train NN
    % Neural network
    stateDelays = [1];%[1,Ndelay];     % state delay vector
    inputDelays = [1];%[1,Ndelay];          % input delay vector
    hiddenSizes = [5];%[15 15];       % network structure (number of neurons per layer)
    
    % Nonlinear autoregressive neural network
    net = narxnet(inputDelays,stateDelays, hiddenSizes);
    
    % Training parameters %nnstart
    net.trainFcn = 'trainlm';%'trainb';%'trainlm';
    net.trainParam.min_grad = 1e-6;
    net.trainParam.goal = 1e-6;
    net.trainParam.showCommandLine = 1;
    net.trainParam.epochs = 100;
    for iL = 1:length(hiddenSizes)
        net.layers{iL}.transferFcn = 'purelin'; %tansig
        net.inputs{iL}.processFcns = {'mapminmax','mapstd'};
    end
    net.performParam.regularization = 0.5;
    net.outputs{length(hiddenSizes)+1}.processFcns = {'mapminmax','mapstd'};
%     net.layers{3}.transferFcn = 'purelin';
    % %     net.layers{1}.transferFcn = 'logsig';purelin
    %     net.inputs{2}.processParams{2} = struct('ymin',0,'ymax',1);
    
    %     netc = closeloop(net);
    % Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
    [Us,Ui,Si,Ss] = preparets(net,yi,{},yt); %yt
    %     [Us,Ui,Si,Ss] = preparets(netc,yi,{},yt);
    % yt: 1xM vector of 5x1 cells of states
    % yi: 1xM vector of 1x1 cells of inputs
    % Us: 2x(M-1) vector with 1:M-1 elements of (1)yi and (2)yt (SHIFTED INPUTS)
    % Ui: first column vector in Us, Ui{1} = Us{1,1}, Ui{2} = Us{2,1};
    % Si: Empty cell array 2x0
    % Ss: 1x(M-1) vector with 1:M(SHIFTED TARGETS)
    
    
    
    % Train net with prepared training data in open-loop
    tic
    net_old = net;
    net = train(net,Us,Ss,Ui,Si);
    
    toc
    
    % Close loop for recursive prediction
    netc = closeloop(net);
    netc2 = netc;
    %% Closed-loop training
%         % Close loop for recursive prediction
%         netc = closeloop(net);
%     
%         [Us,Ui,Si,Ss] = preparets(netc,yi,{},yt);
%         netc.trainParam.epochs = 500;
%         netc.trainParam.goal = 1e-6;
%         netc.trainParam.mu = 1e-3; %1e-3
%         netc.trainParam.mu_dec = 0.05; %0.1
%         netc.trainParam.mu_inc = 5; %10
%         netc2 = train(netc,Us,Ss,Ui,Si);
    
    %% Prediction over training phase
    iIC = 2;%Nic-1;
    xvalid = xtrain(:,1:Nmodel,iIC)';
    xvalid = xvalid - repmat(xmean(1:Nmodel),[1 size(xvalid,2)]); xvalid = xvalid./repmat(xnorm(1:Nmodel),[1 size(xvalid,2)]);
    uvalid = u(:,1,iIC)';
    xvalid = xvalid(:,Nshift+1:end);
    uvalid = uvalid(Nshift+1:end);
    tvalid = tspan(Nshift+1:end);
    
    yt = con2seq(xvalid);
    yi = con2seq(uvalid);
    %     [Us,Ui,Si,Ss] = preparets(net,yi,{},yt);
    
    % Prepare validation data / Get initial  tate from training data
    [Us,Ui,Si,So] = preparets(netc,yi,{},yt);
    
    %     Si{4} = yt{1};
    % Predict on validation data
    predict = netc2(Us,Ui,Si);
    xNARX = cell2mat(predict)';
    
    xNARX = xNARX + repmat(xmean(1:Nmodel)',[size(xNARX,1) 1]);
    xNARX = xNARX.*repmat(xnorm(1:Nmodel)',[size(xNARX,1) 1]);
    
    
    % Error
    e = cell2mat(gsubtract(So,predict));
    
    % %
    % figure;
    % plot(tspan(max(stateDelays):end-1),xNARX','-','LineWidth',1,'Color','k');%0.7*ones(1,3))
    % grid on, hold on
    xNARX = [xtrain(1:max(stateDelays),1:Nmodel,iIC);xNARX];
    
    
    % Show validation
    Nvar = 5;
    clear ph
    figure,box on,
    ccolors = get(gca,'colororder');
    ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
        ccolors(2,:)-[0.1 0.2 0.09];
        ccolors(3,:)-[0.1 0.2 0.09];
        ccolors(4,:)-[0.1 0.1 0.2];
        ccolors(5,:)-[0.1 0.2 0.09]];
    plot([tvalid(max(stateDelays)),tvalid(max(stateDelays))],[0 10],'-k','LineWidth',2),hold on,
    for i = 1:Nvar
        ph(i) = plot(tvalid,exp(xtrain(Nshift+1:end,i,iIC)),'-','Color',ccolors(i,:),'LineWidth',1); hold on
    end
    
    for i = 1:Nmodel
        ph(Nvar+i) = plot(tvalid(0+1:end),exp(xNARX(:,i)),'--','Color',ccolors_valid(i,:),'LineWidth',2);
    end
    ylim([0 10])
    xlabel('Time')
    ylabel('xi')
    legend(ph([1,6]),'True',ModelName)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    if SUBSTRACT_MEAN == 1
        print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_SUBSTRACT_MEAN.eps']);
    else
        print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);
    end
    
    return
    %     norm(x(2:end,:)-xNARX,'fro')
    
%         %% Closed-loop training
%         % Close loop for recursive prediction
%     
%         [Us,Ui,Si,Ss] = preparets(netc2,yi,{},yt);
%         netc.trainParam.epochs = 500;
%         netc.trainParam.goal = 1e-6;
%         netc.trainParam.mu = 1e-3; %1e-3
%         netc.trainParam.mu_dec = 0.05; %0.1
%         netc.trainParam.mu_inc = 5; %10
%         netc2 = train(netc2,Us,Ss,Ui,Si);
    
else
    uv = uv';
    
    %% NARX: Training
    rng(2,'twister')
    
    %% Prepare training data // 1 sample
    % prepare training data
    yt = con2seq([xtrain']); %;ones(1,size(x,1))
    yi = con2seq(u');
    
    %% Neural network
    stateDelays = [1];%[1,Ndelay];     % state delay vector
    inputDelays = [1];%[1,Ndelay];     % input delay vector
    hiddenSizes = [5];             % network structure (number of neurons per layer)
    
    % Nonlinear autoregressive neural network
    net = narxnet(inputDelays,stateDelays, hiddenSizes);
    
    % Training parameters %nnstart
    net.trainFcn = 'trainlm';%'trainbr'; %'trainlm'; trainscg
    net.trainParam.min_grad = 1e-10;
    net.trainParam.showCommandLine = 1;
    net.trainParam.goal = 1e-8;
    net.trainParam.epochs = 500;
    % net.divideParam.trainRatio = 70/100;
    % net.divideParam.valRatio = 15/100;
    % net.divideParam.testRatio = 15/100;
    % net.performFcn = 'mse';  % Mean squared error
    net.layers{1}.transferFcn = 'purelin'; %tansig,purelin,logsig
%     net.layers{2}.transferFcn = 'purelin';
%     net.layers{3}.transferFcn = 'purelin';
%     net.layers{4}.transferFcn = 'purelin';
    %     net.inputs{2}.processParams{2} = struct('ymin',-1,'ymax',1);
    
    %     netc = closeloop(net);
    
    % Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
    [Us,Ui,Si,Ss] = preparets(net,yi,{},yt);
    %     [Us,Ui,Si,Ss] = preparets(netc,yi,{},yt);
    
    % Train net with prepared training data in open-loop
    tic
    net = train(net,Us,Ss,Ui,Si);
    %     netc = train(netc,Us,Ss,Ui,Si);
    toc
    
    % Close loop for recursive prediction
    netc = closeloop(net);
    %% Train closed loop network
%         netc = closeloop(net);
%         netc.trainParam.min_grad = 1e-10;
%         netc.trainParam.showCommandLine = 1;
%         netc.trainParam.epochs = 500;
%         net.trainParam.goal = 1e-8;
%     
%         % % Train net with prepared data in closed-loop
%         [Us,Ui,Si,So] = preparets(netc,yi,{},yt);
%         netc=train(netc,Us,So,Ui);
%         %
%         % % Predict on validation data
%         % predict = netc(Us,Ui,Si);
%         %
%         % % Performance
%         % perfc = perform(netc,predict,So);
%         % view(netc);
%     
%         % netc = removedelay(netc);
    %% Prediction over training phase
    % Prepare validation data / Get initial state from training data
    %     xtrain = x;
    %     for i = 1:size(x,2)
    %         xtrain(:,i) = xtrain(:,i)-repmat(xmean(i),[size(x,1) 1]);
    %     end
    % prepare training data
    yt = con2seq([xtrain']); %;ones(1,size(x,1))
    yi = con2seq(u');
    
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
    
    % Show validation
    clear ph
    figure,box on,
    ccolors = get(gca,'colororder');
    ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
        ccolors(2,:)-[0.1 0.2 0.09];
        ccolors(3,:)-[0.1 0.2 0.09];
        ccolors(4,:)-[0.1 0.1 0.2];
        ccolors(5,:)-[0.1 0.2 0.09]];
    for i = 1:Nvar
        ph(i) = plot(tspan,xtrain(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
    end
    
    for i = 1:Nvar
        ph(Nvar+i) = plot(tspan(max(stateDelays)+1:end),xNARX(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
    end
    xlabel('Time')
    ylabel('xi')
    legend(ph([1,6]),'True',ModelName)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    if SUBSTRACT_MEAN == 1
        print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_SUBSTRACT_MEAN.eps']);
    else
        print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);
    end
    
end
return
%% Save Data
Model.name = 'NARX';
Model.net = netc;
Model.xmean = xmean;
Model.xnorm = xnorm;
% Model.umean = umean;
% Model.Ndelay = Ndelay;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.NORMALIZE = NORMALIZE;
Model.TRANSFORM_LOG  = TRANSFORM_LOG;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')