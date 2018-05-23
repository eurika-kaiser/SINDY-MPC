% HIV model
% System identification using neural network model
 
clear all, close all, clc
figpath = '../FIGURES/HIV/';
datapath = '../DATA/HIV/';
addpath('../utils');

SystemModel = 'HIV';
ModelName = 'NARX';

%% Generate Data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'prbs'; 
DATA_ENSEMBLE = 1;
getTrainingData

%% Parameters & data preparation
Nvar = size(x,2);    
TRANSFORM_LOG = 1;
if TRANSFORM_LOG == 1
    xtrain = log(abs(x));
else
    xtrain = x;
end
SUBSTRACT_MEAN = 0;
if SUBSTRACT_MEAN ==1
    xmean = xSTEADY1';
else
    xmean = zeros(Nvar,1);
end
NORMALIZE = 0; xnorm = ones(Nvar,1);

Nvar = size(x,2);

%% Run model identification
if DATA_ENSEMBLE == 1
    rng(2,'twister') % For reproducibility
    
    % Reshape & prepare data
    Nshift = 0;
    M = size(xtrain,1)-Nshift; 
    yt = cell(1,M);
    yi = cell(1,M);
    for i = Nshift+1:size(xtrain,1)%M
        yt{i-Nshift} = squeeze(xtrain(i,1:Nvar,:)) - repmat(xmean(1:Nvar),[1 Nic]);
        yi{i-Nshift} = squeeze(u(i,:,:))';
    end
    
    % Train NN
    % Neural network
    stateDelays = [1];     % state delay vector
    inputDelays = [1];     % input delay vector
    hiddenSizes = [5];     % network structure (number of neurons per layer)
    
    % Nonlinear autoregressive neural network, AR process not really used
    % as Ndelay = 1
    net = narxnet(inputDelays,stateDelays, hiddenSizes);
    
    % Training parameters 
    net.trainFcn            = 'trainlm';
    net.trainParam.min_grad = 1e-10;
    net.trainParam.goal     = 1e-8;
    net.trainParam.showCommandLine = 1;
    net.trainParam.epochs   = 100;
    for iL = 1:length(hiddenSizes)
        net.layers{iL}.transferFcn = 'purelin';  %'logsig';transig
    end
%     net.inputs{1}.processFcns = {'mapminmax','mapstd'};
%     net.inputs{2}.processFcns = {'mapminmax','mapstd'};
%     net.outputs{length(hiddenSizes)+1}.processFcns = {'mapminmax','mapstd'};
%     net.performParam.regularization = 0.5;

    % Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
    [Us,Ui,Si,Ss] = preparets(net,yi,{},yt); 
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

    %% Prediction over training phase
    netc = closeloop(net);
    iIC = 1; % validate on one trajectory from ensemble
    xvalid = xtrain(:,1:Nvar,iIC)';
    xvalid = xvalid - repmat(xmean(1:Nvar),[1 size(xvalid,2)]); xvalid = xvalid./repmat(xnorm(1:Nvar),[1 size(xvalid,2)]);
    uvalid = u(:,1,iIC)';
    xvalid = xvalid(:,Nshift+1:end);
    uvalid = uvalid(Nshift+1:end);
    tvalid = tspan(Nshift+1:end);
    
    yt = con2seq(xvalid);
    yi = con2seq(uvalid);
    
    % Prepare validation data / Get initial  tate from training data
    [Us,Ui,Si,So] = preparets(netc,yi,{},yt);
 
    % Predict on validation data
    predict = netc(Us,Ui,Si);
    xNARX = cell2mat(predict)';
    
    xNARX = xNARX + repmat(xmean(1:Nvar)',[size(xNARX,1) 1]);
    xNARX = xNARX.*repmat(xnorm(1:Nvar)',[size(xNARX,1) 1]);
    
    
    % Error
    e = cell2mat(gsubtract(So,predict));

    xNARX = [xtrain(1:max(stateDelays),1:Nvar,iIC);xNARX];
    
    
    %% Show prediction over training data
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
    
    for i = 1:Nvar
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
    
    
else
    uv = uv';
    
    rng(2,'twister') % for reproducibility
    
    %% Prepare training data 
    yt = con2seq([xtrain']); %;ones(1,size(x,1))
    yi = con2seq(u');
    
    %% Neural network
    stateDelays = [1];     % state delay vector
    inputDelays = [1];     % input delay vector
    hiddenSizes = [5];     % network structure (number of neurons per layer)
    
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
    
    % Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
    [Us,Ui,Si,Ss] = preparets(net,yi,{},yt);
    
    % Train net with prepared training data in open-loop
    tic
    net = train(net,Us,Ss,Ui,Si);
    toc
    
    % Close loop for recursive prediction
    netc = closeloop(net);
    
    %% Prediction over training phase
    % prepare training data
    yt = con2seq([xtrain']); 
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
    
    %% Show prediction over training stage
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

%% Save Data
Model.name = 'NARX';
Model.net = netc;
Model.xmean = xmean;
Model.xnorm = xnorm;
Model.stateDelays = stateDelays;
Model.inputDelays = inputDelays;
Model.hiddenSizes = hiddenSizes;
Model.SUBSTRACT_MEAN = SUBSTRACT_MEAN;
Model.NORMALIZE = NORMALIZE;
Model.TRANSFORM_LOG  = TRANSFORM_LOG;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')