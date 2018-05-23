% HIV system
% System identification: (Delay)DMDc


clear all, close all, clc
figpath = '../FIGURES/HIV/'; mkdir(figpath)
datapath = '../DATA/HIV/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'HIV';

% Model type selected based on DATA_ENSEMBLE and  Ndelay
%% Generate Data
InputSignalType = 'prbs';
Nvar = 5;
Ndelay = 1;  % Choose 1 for DMDc, and 10 for DelayDMDc in paper
DATA_ENSEMBLE = 0; % Choose 0 to reproduce results in paper

%%
if DATA_ENSEMBLE == 1 % ONLY FOR Ndelay = 1 implemented
    ModelName = 'DMDc';
    getTrainingData_Ensemble
    
    %% Reshape
    X = x(1:end-1,:,:);         % Array of snapshots
    Xp = x(2:end,:,:);          % time-shifted Array
    U = u(1:end-1,:,:);         % Array of inputs
    M = size(X,1);
    n = size(X,2);
    
    % Reorganize data from array to time-state matrix
    X_tmp = zeros(Nvar,Nic*M);
    Xp_tmp = zeros(Nvar,Nic*M);
    U_tmp = zeros(1,Nic*M);
    for iIC = 1:Nic
        for i = 1:Nvar
            X_tmp(i,(iIC-1)*M+1:iIC*M) = X(:,i,iIC);
            Xp_tmp(i,(iIC-1)*M+1:iIC*M) = Xp(:,i,iIC);
        end
        U_tmp(1,(iIC-1)*M+1:iIC*M) = U(:,1,iIC);
    end
    X = X_tmp; Xp = Xp_tmp; U = U_tmp;
    clear X_tmp Xp_tmp U_tmp
    M = size(X,2);
    
    %% Construct data matrices // Mean-correction
    xmean = zeros(1,5); %mean(X,2)';%xref; 
    X   = X - repmat(xmean',[1 M]);
    Xp  = Xp - repmat(xmean',[1 M]);
    
    r1 = size(X,1); r2 = size(Xp,1);
    [sysmodel_DMDc,Psi,Psi_p] = DelayDMDc_MV(zeros(5,1),U,r1,r2,dt,size(X,1),size(U,1),2,X,Xp);
    
    Nt = length(t)-1;
    
    %% Validation over training phase
    xDMDc_Valid = zeros(size(X));
    for iIC = 1:Nic
        [xDMDc,~] = lsim(sysmodel_DMDc,squeeze(u(1:Nt,1,iIC)),t(1:Nt),x0_ensemble(iIC,:));
        xDMDc = xDMDc + repmat(xmean,[Nt 1]);
        xDMDc_Valid(1:5,(iIC-1)*Nt+1:iIC*Nt) = xDMDc';
        disp(['PROGRESS: ',num2str(100*iIC/Nic),'%'])
    end
    %%
    figure, hold on, box on
    plot(X','-k','LineWidth',2)
    plot(xDMDc_Valid','--','LineWidth',2)
else
    
    getTrainingData
    %% DMDc: B = unknown  and with time delay coordinates
    Ndelay = 1; % 1 for DMDc, 10 for DelayDMDc in paper  
    if Ndelay == 1
        ModelName = 'DMDc';
    elseif Ndelay>1
        ModelName = 'DelayDMDc';
    end
    
    % Construct data matrices
    Hu = getHankelMatrix_MV(u,Ndelay);
    xmean = xref; %zeros(1,Nvar);%xref;%mean(x);%xref; %mean(x);
    X   = x - repmat(xmean,[T 1]);
    Hx  = getHankelMatrix_MV(X,Ndelay);
    numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 5;
    r1 = size(Hx,1); r2 = size(Hx,1);
    [sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);
    
    Nt = length(t)-Ndelay+1;
    %% Validation over training phase
    [xDMDc,~] = lsim(sysmodel_DMDc,Hu',tspan(1:Nt),Hx(:,1));
    xDMDc = xDMDc(:,end-Nvar+1:end);
    xDMDc = xDMDc + repmat(xmean,[Nt 1]);
    
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
        ph(i) = semilogy(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
    end
    for i = 1:Nvar
        ph(Nvar+i) = semilogy(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
    end
    xlabel('Time')
    ylabel('xi')
    legend(ph([1,6]),'True',ModelName)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);
    
end
%% Save Data & Model
Model.name = ModelName;
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.dt = dt;
Model.xmean = xmean;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')