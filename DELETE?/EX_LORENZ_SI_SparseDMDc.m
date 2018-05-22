% LORENZ system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../../FIGURES/LORENZ/';
datapath = '../../DATA/LORENZ/';
addpath('../utils');

SystemModel = 'LORENZ';
ModelName = 'sparseDMDc';
%% Generate Data
InputSignalType = 'sphs';%prbs; chirp; noise; sine2; sphs; mixed
Ndelay = 1;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

Nt = length(tspan)-1;

%% DMDc: B = unknown  and with time delay coordinates
ModelNumber = 1; % or 2
xrefs = [xref1,xref2];

Hu = getHankelMatrix_MV(u',1);
for i = 1:1%ModelNumber
    xmean{i} = xrefs(:,i)'; %mean(x);
    X   = x - repmat(xmean{i},[T 1]);
    Hx  = getHankelMatrix_MV(X,1);
    numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 3;
    r1 = size(Hx,1); r2 = size(Hx,1);
    [sysmodel_DMDc{i},U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);
end
%% Sparse DMDc
Nstates = size(Hx,1);
Ninputs = size(Hu,1);
epsilon = 1e-1;
lambda  = 0.1;

Y  = Hx(:,1:end-1)+randn(size(Hx(:,1:end-1)));
Yp = Hx(:,2:end)+randn(size(Hx(:,1:end-1)));
U  = Hu(1:end-1);
% using the l1 norm to promote sparsity
% cvx_begin %quiet
%     variable M( Nstates, Nstates+Ninputs);
%     minimize( norm(M(:),1));
%     subject to
%         norm(Xp - M(:,1:Nstates)*X - M(:,Nstates+1:Nstates+Ninputs)*U, 'fro') <= epsilon; %#ok<VUNUS>
% cvx_end
% A = M(:,1:Nstates);
% B = M(:,Nstates+1:Nstates+Ninputs);
% cvx_begin %quiet
%     variable M( Nstates, Nstates+Ninputs);
%     minimize( norm(M(:),1));
%     subject to
%         norm(Xp - M*[X;U], 'fro') <= epsilon; %#ok<VUNUS>
% cvx_end

M = size(Y,2);
% Gset = zeros(M,3);
for i = 1:Nstates
    [G{i}, FitInfo] = lasso([Y(:,:);U]',Yp(i,:)','Lambda',0.0000000); %,'Lambda',0.00000001
%     lassoPlot(G{i},FitInfo);
end

A = [G{1}(1:3)';G{2}(1:3)';G{3}(1:3)'];
B = [G{1}(4)';G{2}(4)';G{3}(4)'];
C = eye(Nstates,Nstates);
D = zeros(Nstates,1);

sysmodel_sparseDMDc{1} = ss(A,B,C,D,dt);

%% Prediction over training phase
for i = 1:ModelNumber
    [xDMDc{i},~] = lsim(sysmodel_DMDc{i},Hu',tspan(1:end-1),x(1,:)'-xmean{i}');
    xDMDc{i} = xDMDc{i} + repmat(xmean{i},[length(tspan)-1 1]);
    
    [xDMDc_sparse{i},~] = lsim(sysmodel_sparseDMDc{i},Hu',tspan(1:end-1),x(1,:)'-xmean{i}');
    xDMDc_sparse{i} = xDMDc_sparse{i} + repmat(xmean{i},[length(tspan)-1 1]);
end

%% Show validation
for i = 1:ModelNumber
    clear ph
    figure,box on,
    ccolors = get(gca,'colororder');
    ph(1) = plot(tspan,x(:,1),'-','Color','k','LineWidth',1); hold on
    ph(2) = plot(tspan,x(:,2),'-','Color','k','LineWidth',1);
    ph(3) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc{i}(:,1),'--','Color','r','LineWidth',2);
    ph(4) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc{i}(:,2),'--','Color','r','LineWidth',2);
    ph(5) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc_sparse{i}(:,1),':','Color','b','LineWidth',2);
    ph(6) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc_sparse{i}(:,2),':','Color','b','LineWidth',2);
    xlim([0 (length(tspan)-1)*dt]), ylim([-25 50])
    xlabel('Time')
    ylabel('Population size')
    % legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
    legend(ph([1,3,5]),'True','DMDc',ModelName)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
%     print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',i,'.eps']);
end

return
%% Prediction
% Reference
tspanV   = [10:dt:20];
xA      = xv;
tA      = tv;

% Model
for i = 1:ModelNumber
    if Ndelay == 1
        x0      = [x(end,1:3)];
        Hunew   = [u(end),uv(1:end)];
        [xBm{i},tB] = lsim(sysmodel_DMDc{i},Hunew,tspanV,[x0-[xmean{i}]]');
    elseif Ndelay > 1
        x0      = [x(end-Ndelay+1,1:3),x(end,1:3)];
        Hunew   = [ u(end-Ndelay+1:end),uv(1:end-Ndelay);
            u(end),uv(1:end-1)];
        [xBm,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
        xBm = xBm(:,4:6); xBm = xBm + repmat(xmean,[size(xBm,1) 1]);
    end
    
    xBm{i} = xBm{i} + repmat(xmean{i},[length(tB) 1]);
end
%% Show training and prediction
xB = xBm{1};
VIZ_SI_Validation

xB = xBm{2};
VIZ_SI_Validation

%% Save Data
Model.name = 'DelayDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.xmean = xmean;
Model.xrefs = xrefs;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')