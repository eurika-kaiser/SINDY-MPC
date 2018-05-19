% LOTKA-VOLTERRA system
% System identification: NARX

clear all, close all, clc
figpath = '../../FIGURES/HIV/';
datapath = '../../DATA/HIV/';
addpath('../utils');

SystemModel = 'HIV';
ModelName = 'NARX';
%% Generate Data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'sine2'; % prbs; chirp; noise; sine2; sphs; mixed
getTrainingData_Ensemble

% Data.x = x;
% 
% save('DATA_NARX_prbs.mat',)
% return
% xforc = x;
% uforc  = u;
% %%
% 
% figure, hold on, box on
% for i = 1:Nic
%     plot(x(:,4,i),'-k')
%     plot(xforc(:,4,i),'--r')
% end
%%
DATA_ENSEMBLE = 1;
Nvar = 5;
if DATA_ENSEMBLE == 1
    
    M = size(x,1);
    
    SUBSTRACT_MEAN = 0;
    NORMALIZE = 1;
    if NORMALIZE == 1
        xnorm = max(max(x,[],1),[],3)';
    else
        xnorm = ones(Nvar,1);
    end
    if  SUBSTRACT_MEAN == 1
        xmean = mean(mean(x,1),3)';
    else
        xmean = zeros(Nvar,1); % %zeros(Nvar,1); %mean(mean(x,3),1)'; %mean(xic);
    end
    
    
    % Reshape & prepare data
    xtrain = cell(1,M);
    utrain = cell(1,M);
    for i = 1:M
        xtrain{i} = squeeze(x(i,:,:));
        xtrain{i} = xtrain{i}./repmat(xnorm,[1 size(xtrain{i},2)]);
        utrain{i} = squeeze(u(i,1,:))';
    end
%     xtrain0 = x(1:end,:,:);
%     xtrain = xtrain - repmat(xmean,[1 size(xtrain,2)]);
%     xtrain = xtrain./repmat(xnorm,[1 size(xtrain,2)]);

    % Parameters
    rng(2,'twister')
        
    
    % prepare training data
%     yt = con2seq(xtrain);
%     yi = con2seq(utrain);
    yt = xtrain;
    yi = utrain;

    
    % Train NN
    % Neural network
    stateDelays = 1;%[1,Ndelay];     % state delay vector
    inputDelays = 1;%[1,Ndelay];          % input delay vector
    hiddenSizes = [10 10];%[15 15];       % network structure (number of neurons per layer)

    % Nonlinear autoregressive neural network
    net = narxnet(inputDelays,stateDelays, hiddenSizes);
    
    % Training parameters %nnstart
    net.trainFcn = 'trainlm';%'trainb';%'trainlm';
    net.trainParam.min_grad = 1e-8;
    net.trainParam.showCommandLine = 1;
    net.trainParam.epochs = 100;
%     net.layers{1}.transferFcn = 'logsig';

    
    % Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
    [Us,Ui,Si,Ss] = preparets(net,yi,{},yt); %yt
    % yt: 1xM vector of 5x1 cells of states
    % yi: 1xM vector of 1x1 cells of inputs
    % Us: 2x(M-1) vector with 1:M-1 elements of (1)yi and (2)yt (SHIFTED INPUTS)
    % Ui: first column vector in Us, Ui{1} = Us{1,1}, Ui{2} = Us{2,1}; 
    % Si: Empty cell array 2x0
    % Ss: 1x(M-1) vector with 1:M(SHIFTED TARGETS)
    
    % Train net with prepared training data in open-loop
    tic
    net = train(net,Us,Ss,Ui,Si);
    net_old = net;
%     for iIC = 2:Nic
%         % Update data stream
%         xtrain = x(:,:,iIC)'; xtrain = xtrain - repmat(xmean,[1 size(xtrain,2)]); xtrain = xtrain./repmat(xnorm,[1 size(xtrain,2)]);
%         utrain = u(:,1,iIC)';
%         yt = con2seq(xtrain);
%         yi = con2seq(utrain);
%         [Us,Ui,Si,Ss] = preparets(net,yi,{},yt);
%         % Re-train network
%         net_new = train(net_old,Us,Ss,Ui,Si);
%         net_old = net_new;
%     end
    toc
   
    % Close loop for recursive prediction
    netc = closeloop(net);
    
 
    %% Prediction over training phase
    iIC = 1;
    xpred = x(:,:,iIC)'; xpred = xpred - repmat(xmean,[1 size(xpred,2)]); xpred = xpred./repmat(xnorm,[1 size(xpred,2)]);
    upred = u(:,1,iIC)';
    yt = con2seq(xpred);
    yi = con2seq(upred);
%     [Us,Ui,Si,Ss] = preparets(net,yi,{},yt);

    % Prepare validation data / Get initial state from training data
    [Us,Ui,Si,So] = preparets(netc,yi,{},yt);
    
%     Si{4} = yt{1};
    % Predict on validation data
    predict = netc(Us,Ui,Si);
    xNARX = cell2mat(predict)';
    
    xNARX = xNARX + repmat(xmean',[size(xNARX,1) 1]);
    xNARX = xNARX.*repmat(xnorm',[size(xNARX,1) 1]);
    
    
    % Error
    e = cell2mat(gsubtract(So,predict));
    
    % %
    % figure;
    % plot(tspan(max(stateDelays):end-1),xNARX','-','LineWidth',1,'Color','k');%0.7*ones(1,3))
    % grid on, hold on
    xNARX = [x(1,:,Nic);xNARX];
    
    
    %% Show validation
    
    clear ph
    figure,box on,
    ccolors = get(gca,'colororder');
    ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
        ccolors(2,:)-[0.1 0.2 0.09];
        ccolors(3,:)-[0.1 0.2 0.09];
        ccolors(4,:)-[0.1 0.1 0.2];
        ccolors(5,:)-[0.1 0.2 0.09]];
    for i = 1:Nvar
        ph(i) = plot(tspan,x(:,i,Nic),'-','Color',ccolors(i,:),'LineWidth',1); hold on
    end
    
    for i = 1:Nvar
        ph(Nvar+i) = plot(tspan(1:end),xNARX(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
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
    
    
%     norm(x(2:end,:)-xNARX,'fro')
else
    u = u';
    uv = uv';
    %% NARX: Training
    ModelName = 'NARX';
    Nvar = 3;
    rng(2,'twister')
    SUBSTRACT_MEAN = 1;
    xmean = xref;%mean(x)';
    
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
    stateDelays = 1;%[1,Ndelay];%[1:50:501];%1:10:100;     % state delay vector
    inputDelays = 1;%[1,Ndelay];%[1:50:501];%1:10:100;          % input delay vector
    hiddenSizes = [10];%[15 15];%[10]; %[20 20];       % network structure (number of neurons per layer)
    %[15 15] does not work
    %[15] works
    %[10] better, after 1 period increasing phase difference, 47.10s
    %[10 10] wors, amplitude too small, phase shift
    %[20 20] very good, but takes long to train
    %[15 5] also works but not as good as [15]?
    %[15 15] bad
    %[15 20] good, but takes longer
    %[5 5 5] good, 86.97s
    
    % Nonlinear autoregressive neural network
    net = narxnet(inputDelays,stateDelays, hiddenSizes);
    
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
    ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
        ccolors(2,:)-[0.1 0.2 0.09];
        ccolors(3,:)-[0.1 0.2 0.09];
        ccolors(4,:)-[0.1 0.1 0.2];
        ccolors(5,:)-[0.1 0.2 0.09]];
    for i = 1:Nvar
        ph(i) = plot(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
    end
    
    for i = 1:Nvar
        ph(Nvar+i) = plot(tspan(2:end),xNARX(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
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
    % VIZ_SI_Validation
end

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
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')