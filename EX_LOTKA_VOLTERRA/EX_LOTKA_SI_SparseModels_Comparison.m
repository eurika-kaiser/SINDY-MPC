% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

ModelName = 'SparseModels';
%% Generate Data
InputSignalType = 'sphs';%prbs; chirp; noise; sine2; sphs; mixed
Ndelay = 1;%35;%35;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

DataTrain.x = x;
DataTrain.t = t;
DataTrain.tspan = tspan;
DataTrain.u = u;
DataTrain.xmean = xmean;
xstd = std(DataTrain.x(:,1));

%% Parameters
close all
rng(0,'twister')

eta_vec = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
N_ETA = length(eta_vec);

lambda_vec = [0,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1];
N_LAMBDA = length(lambda_vec);

model_collection = {'sDMDc','sEDMDc2', 'sEDMDc3','SINDYc'};
N_MODEL = length(model_collection);

N_REP = 30;

options_method.sparsify = 'ILSTH'; % 'LASSO'
options_method.usesine  = 0;

nVars = 2;
%% Initialization
TrainingError = inf*ones(N_ETA,N_LAMBDA,N_MODEL,N_REP);
ValidationError = inf*ones(N_ETA,N_LAMBDA,N_MODEL,N_REP);

PredictionResults_Training = zeros(1000,nVars,N_ETA,N_LAMBDA,N_MODEL,N_REP);
PredictionResults_Validation = zeros(1000,nVars,N_ETA,N_LAMBDA,N_MODEL,N_REP);

%% Model Identification

for iR = 1:N_REP
    counter = 0;
    
    tstart = tic;
    for iEta = 1:N_ETA
        
        % Add noise to data
        eps     = eta_vec(iEta)*xstd;
        x       = DataTrain.x + eps*randn(size(DataTrain.x));
        u       = DataTrain.u;
        M       = length(x);
        
        for iLambda = 1:N_LAMBDA
            lambda = lambda_vec(iLambda);
            
            %% Model Identification
            Hx  = getHankelMatrix_MV(x - repmat(DataTrain.xmean,[M 1]),1); % No time delay (=1)
            Hu  = getHankelMatrix_MV(u',1);
            
            % sparse DMDc
            options_method.order  = 1;
            yout = poolData(Hx',nVars,options_method.order,0);
            Y = yout(1:end-1,2:end)'; % remove constant
            Yp = yout(2:end,2:end)';
            U = Hu(1:end-1);
            sys_DMDc = SparseRegression(Y,Yp,U,dt, options_method, lambda);
            
            % sparse EDMDc
            options_method.order  = 2;
            yout = poolData(Hx',nVars,options_method.order,0);
            Y = yout(1:end-1,2:end)'; % remove constant
            Yp = yout(2:end,2:end)';
            U = Hu(1:end-1);
            sys_EDMDc2 = SparseRegression(Y,Yp,U,dt, options_method, lambda);
            
            % sparse EDMDc
            options_method.order  = 3;
            yout = poolData(Hx',nVars,options_method.order,0);
            Y = yout(1:end-1,2:end)'; % remove constant
            Yp = yout(2:end,2:end)';
            U = Hu(1:end-1);
            sys_EDMDc3 = SparseRegression(Y,Yp,U,dt, options_method, lambda);
            
            % SINDYc (not xmean-substracted)
            options_method.order  = 3;
            Xi = NonlinearSparseRegression(x,u',dt,options_method,lambda);
            
            %% Prediction on training data
            Nt = M-1;
            % sparse DMDc
            options_method.order  = 1;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xDMDc,~] = lsim(sys_DMDc,Hu,tspan(1:Nt),x0train);
            xDMDc = xDMDc(:,end-1:end);
            xDMDc = xDMDc + repmat(xmean,[Nt 1]);
            xPredTrain{1}.x = xDMDc;
            
            % sparse EDMDc
            options_method.order  = 2;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xEMD2,~] = lsim(sys_EDMDc2,Hu,tspan(1:Nt),x0train);
            xEMD2 = xEMD2(:,1:2);
            xEMD2 = xEMD2 + repmat(xmean,[Nt 1]);
            xPredTrain{2}.x = xEMD2;
            
            % sparse EDMDc
            options_method.order  = 3;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xEMD3,~] = lsim(sys_EDMDc3,Hu,tspan(1:Nt),x0train);
            xEMD3 = xEMD3(:,1:2);
            xEMD3 = xEMD3 + repmat(xmean,[Nt 1]);
            xPredTrain{3}.x = xEMD3;
            
            % SINDYc
            p.ahat = Xi(:,1:2);
            p.polyorder = options_method.order; p.usesine = options_method.usesine; p.dt = dt;
            [Ns,N] = size(Hx);
            xSINDYc = zeros(Ns,N); xSINDYc(:,1) = x(1,:)';
            for ct=1:N-1
                xSINDYc(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xSINDYc(:,ct),DataTrain.u(ct),dt,1,[],p);
            end
            xSINDYc = xSINDYc';
            xPredTrain{4}.x = xSINDYc;
            
            % Store results
            for iM = 1:N_MODEL
                PredictionResults_Training(:,:,N_ETA,N_LAMBDA,iM,N_REP) = xPredTrain{iM}.x;
            end
            
            %             clear ph
            %             figure,box on,
            %             ccolors = get(gca,'colororder');
            %             ph(1) = plot(tspan(1:Nt),x(1:Nt,1),'-','Color','k','LineWidth',1); hold on
            %             plot(tspan(1:Nt),x(1:Nt,2),'-','Color','k','LineWidth',1);
            %             for iM = 1:N_MODEL
            %                 for iVar = 1:nVars
            %                     ph(iM+1) = plot(tspan(1:Nt),xPredTrain{iM}.x(:,iVar),'--','Color',ccolors(iM,:),'LineWidth',1); hold on
            %                 end
            %             end
            %             legend(ph,model_collection)
            
            for iM = 1:N_MODEL
                TrainingError(iEta,iLambda,iM,iR) = norm(x(1:Nt,:)-xPredTrain{iM}.x,'fro');
            end
            
            %disp(['Training Error: ', num2str(squeeze(TrainingError(iEta,iLambda,:,iR))')])
            %             disp(['                DMDc: ', num2str(DMDc_err)])
            %             disp(['          sparseDMDc: ', num2str(sparseDMDc_err)])
            %% Prediction on validation data
            Nt = M-1;
            % sparse DMDc
            options_method.order  = 1;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xDMDc,~] = lsim(sys_DMDc,Hu,tspan(1:Nt),x0train);
            xDMDc = xDMDc(:,end-1:end);
            xDMDc = xDMDc + repmat(xmean,[Nt 1]);
            xPredTrain{1}.x = xDMDc;
            
            % sparse EDMDc
            options_method.order  = 2;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xEMD2,~] = lsim(sys_EDMDc2,Hu,tspan(1:Nt),x0train);
            xEMD2 = xEMD2(:,1:2);
            xEMD2 = xEMD2 + repmat(xmean,[Nt 1]);
            xPredTrain{2}.x = xEMD2;
            
            % sparse EDMDc
            options_method.order  = 3;
            x0train = poolData(Hx(:,1)',nVars,options_method.order,0)';
            x0train = x0train(2:end);
            [xEMD3,~] = lsim(sys_EDMDc3,Hu,tspan(1:Nt),x0train);
            xEMD3 = xEMD3(:,1:2);
            xEMD3 = xEMD3 + repmat(xmean,[Nt 1]);
            xPredTrain{3}.x = xEMD3;
            
            % SINDYc
            p.ahat = Xi(:,1:2);
            p.polyorder = options_method.order; p.usesine = options_method.usesine; p.dt = dt;
            [Ns,N] = size(Hx);
            xSINDYc = zeros(Ns,N); xSINDYc(:,1) = x(1,:)';
            for ct=1:N-1
                xSINDYc(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xSINDYc(:,ct),DataTrain.u(ct),dt,1,[],p);
            end
            xSINDYc = xSINDYc';
            xPredTrain{4}.x = xSINDYc;
            
            % Store results
            for iM = 1:N_MODEL
                PredictionResults_Validation(:,:,N_ETA,N_LAMBDA,iM,N_REP) = xPredTrain{iM}.x;
            end
            
            for iM = 1:N_MODEL
                ValidationError(iEta,iLambda,iM,iR) = norm(x(1:Nt,:)-xPredTrain{iM}.x,'fro');
            end
            
            
            
            counter = counter + 1;
        end
    end
    
    tstop = toc(tstart);
    disp(['STATUS (',num2str(iR),' of ', num2str(N_REP),'): ', num2str(100*counter/(N_ETA*N_LAMBDA)),'%', '[Time=',num2str(tstop),'s]'])
end

%% STATS
idxnan = isnan(TrainingError(:));
maxerr = max(TrainingError(:));
TrainingError(idxnan) = maxerr; % change later
ErrStats.mean = zeros(N_ETA,N_LAMBDA,N_MODEL);
ErrStats.median = zeros(N_ETA,N_LAMBDA,N_MODEL);
ErrStats.std = zeros(N_ETA,N_LAMBDA,N_MODEL);
ErrStats.min = zeros(N_ETA,N_LAMBDA,N_MODEL);
for iEta = 1:N_ETA
    for iLambda = 1:N_LAMBDA
        for iModel = 1:N_MODEL
            ErrStats.mean(iEta,iLambda,iModel) = mean(TrainingError(iEta,iLambda,iModel,:));
            ErrStats.median(iEta,iLambda,iModel) = median(TrainingError(iEta,iLambda,iModel,:));
            ErrStats.std(iEta,iLambda,iModel) = std(TrainingError(iEta,iLambda,iModel,:));
            ErrStats.min(iEta,iLambda,iModel) = min(TrainingError(iEta,iLambda,iModel,:));
        end
    end
end

%%
err_axis = [min(TrainingError(:)), 1e3];%max(TrainingError(:))
figure,
for iModel = 1:N_MODEL
    subplot(4,4,(iModel-1)*4+1)
    surf(lambda_vec,eta_vec,ErrStats.median(:,:,iModel)), shading interp
    xlabel('\lambda'), ylabel('\eta'), zlabel('median')
    %     zlim([0 1e6])
    caxis(err_axis);
    set(gca,'xscale','log')
    view(2)
    
    subplot(4,4,(iModel-1)*4+2)
    surf(lambda_vec,eta_vec,ErrStats.mean(:,:,iModel)) , shading interp
    xlabel('\lambda'), ylabel('\eta'), zlabel('mean')
    %     zlim([0 1e6])
    caxis(err_axis);
    set(gca,'xscale','log')
    view(2)
    
    subplot(4,4,(iModel-1)*4+3)
    surf(lambda_vec,eta_vec,ErrStats.std(:,:,iModel)) , shading interp
    xlabel('\lambda'), ylabel('\eta'), zlabel('std')
    %     zlim([0 1e6])
    caxis(err_axis);
    set(gca,'xscale','log')
    view(2)
    
    subplot(4,4,(iModel-1)*4+4)
    surf(lambda_vec,eta_vec,ErrStats.min(:,:,iModel)) , shading interp
    xlabel('\lambda'), ylabel('\eta'), zlabel('min')
    %     zlim([0 1e6])
    caxis(err_axis);
    set(gca,'xscale','log')
    view(2)
end
% for iEta = 1:N_ETA
%     for iLambda = 1:N_LAMBDA
%         plot3()
%     end
% end

%% TIme series stats
iEta = 1; iLambda = 1; iM = 4;
data = squeeze(PredictionResults_Training(:,:,iEta,iLambda,iM,N_REP));

xBmin = min(data(:,1:nVars,:),[],3);
xBmax = max(data(:,1:nVars,:),[],3);
clear ph
figure('visible','off'),box on, hold on,
ccolors = get(gca,'colororder');
plot([tB(1),tB(1)],[-25 65],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot([t(end),t(end)],[-25 65],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
ylim([-25 65])
t1 = text(5,55,'Training', 'FontSize',12);
t2 = text(5+tA(1),55,'Validation', 'FontSize',12);


X=[tB(2:end)',fliplr(tB(2:end)')];                %#create continuous x value array for plotting
Y=[xBmin(:,1)',flipud(xBmax(:,1))'];              %#create y values for out and then back
fillh1 = fill(X,Y,ccolors(1,:)-[0 0.2 0.2]);                  %#plot filled area
fillh1.EdgeColor = ccolors(1,:)-[0 0.2 0.2]; fillh1.FaceAlpha = 0.5;

X=[tB(2:end)',fliplr(tB(2:end)')];                %#create continuous x value array for plotting
Y=[xBmin(:,2)',flipud(xBmax(:,2))'];              %#create y values for out and then back
fillh = fill(X,Y,ccolors(2,:)-[0.1 0.2 0.09]);                  %#plot filled area
fillh.EdgeColor = ccolors(2,:)-[0.1 0.2 0.09]; fillh.FaceAlpha = 0.5;

X=[tB(2:end)',fliplr(tB(2:end)')];                %#create continuous x value array for plotting
Y=[xBmin(:,3)',flipud(xBmax(:,3))'];              %#create y values for out and then back
fillh = fill(X,Y,ccolors(3,:)-[0.1 0.2 0.09]);                  %#plot filled area
fillh.EdgeColor = ccolors(3,:)-[0.1 0.2 0.09]; fillh.FaceAlpha = 0.5;

%             plot(tB(2:end),xBmin(:,1),'-k','LineWidth',2)
%             plot(tB(2:end),xBmin(:,1),'--g','LineWidth',2)

if eps~=0
    ph(4) = plot(tspan,x(:,1),'-','Color',0.7*ones(1,3),'LineWidth',1); %ccolors(1,:)+[0.15 0.3 0.25]
    plot(tspan,x(:,2),'-','Color',0.7*ones(1,3),'LineWidth',1); %ccolors(2,:)+[0.15 0.3 0.25]
    plot(tspan,x(:,3),'-','Color',0.7*ones(1,3),'LineWidth',1);
    ph(1) = plot([DataTrain.t;tA],[DataTrain.x(:,1);xA(:,1)],'-','Color',ccolors(1,:),'LineWidth',1);
    ph(2) = plot([DataTrain.t;tA],[DataTrain.x(:,2);xA(:,2)],'-','Color',ccolors(2,:),'LineWidth',1);
    ph(3) = plot([DataTrain.t;tA],[DataTrain.x(:,3);xA(:,3)],'-','Color',ccolors(3,:),'LineWidth',1);
else
    ph(1) = plot([DataTrain.t;tA],[DataTrain.x(:,1);xA(:,1)],'-','Color',ccolors(1,:),'LineWidth',1);
    ph(2) = plot([DataTrain.t;tA],[DataTrain.x(:,2);xA(:,2)],'-','Color',ccolors(2,:),'LineWidth',1);
    ph(3) = plot([DataTrain.t;tA],[DataTrain.x(:,3);xA(:,3)],'-','Color',ccolors(3,:),'LineWidth',1);
    ph(4) = plot(t,x(:,1),'--','Color',[0 1 0],'LineWidth',1); % Training data
    plot(t,x(:,2),'--','Color',[0 1 0],'LineWidth',1);
    plot(t,x(:,3),'--','Color',[0 1 0],'LineWidth',1);
end

%             ph(4) = plot(tB,xB(:,1),'-.','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
%             ph(5) = plot(tB,xB(:,2),'-.','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
grid off
xlim([0 tv(end)])
xlabel('Time')
ylabel('xi')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto');

if eps~=0
    lh = legend([ph([1,4]),fillh1],'Truth','Training',ModelName,'Location','NorthWest');
    %             lh.Position = [lh.Position(1)+0.13,lh.Position(2)-0.2,lh.Position(3:4)];
    lh.Position = [lh.Position(1)+0.02,lh.Position(2)-0.06,lh.Position(3:4)];
else
    lh = legend(ph([1,4]),'Truth','Training',ModelName,'Location','NorthWest');
    %             lh.Position = [lh.Position(1)+0.13,lh.Position(2)-0.2,lh.Position(3:4)];
    lh.Position = [lh.Position(1)+0.02,lh.Position(2)-0.06,lh.Position(3:4)];
end

print('-depsc2', '-painters','-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.eps']);

delete(lh), delete(t1), delete(t2)
%             print('-depsc2', '-painters','-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Validation_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS_noleg.eps']);

return
%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ph(1) = plot(tspan,x(:,1),'-','Color','k','LineWidth',1); hold on
ph(2) = plot(tspan,x(:,2),'-','Color','k','LineWidth',1);
ph(3) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,1),'--','Color','r','LineWidth',2);
ph(4) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,2),'--','Color','r','LineWidth',2);
ph(5) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc_sparse(:,1),':','Color','b','LineWidth',2);
ph(6) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc_sparse(:,2),':','Color','b','LineWidth',2);
xlim([0 100]), ylim([0,120]);
xlabel('Time')
ylabel('Population size')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,3,5]),'True','DMDc',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.eps']);


DMDc_err = norm(x(1:end-1,:)-xDMDc,'fro');
sparseDMDc_err = norm(x(1:end-1,:)-xDMDc_sparse,'fro');
disp(['Training Error: '])
disp(['                DMDc: ', num2str(DMDc_err)])
disp(['          sparseDMDc: ', num2str(sparseDMDc_err)])


%% Prediction
% Reference
tspanV   = [100:dt:200];
xA      = xv;
tA      = tv;

% Model DMDc
x0      = [x(end,1:2)];
Hunew   = [u(end),uv(1:end)];
[xB,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
xB = xB + repmat(xmean,[length(tB) 1]);

% Model sparseDMDc / sparseEDMDc
if Nstates == 2
    [xB_sparse,tB] = lsim(sysmodel_sparseDMDc,Hunew,tspanV,[x0-[xmean]]');
    xB_sparse = xB_sparse(:,end-1:end);
    xB_sparse = xB_sparse + repmat(xmean,[length(tB) 1]);
elseif Nstates > 2
    x0dmdc = [x0(1)-xmean(1); x0(2)-xmean(2); (x0(1)-xmean(1)).^2; (x0(1)-xmean(1)).*(x0(2)-xmean(2)); (x0(2)-xmean(2)).^2];
    [xB_sparse,tB] = lsim(sysmodel_sparseDMDc,Hunew,tspanV,x0dmdc);
    xB_sparse = xB_sparse(:,1:2);
    xB_sparse = xB_sparse + repmat(xmean,[length(tB) 1]);
end


DMDc_err = norm(xA-xB(1:end-1,:),2);
sparseDMDc_err = norm(xA-xB_sparse(1:end-1,:),2);
disp(['Validation Error: '])
disp(['                DMDc: ', num2str(DMDc_err)])
disp(['          sparseDMDc: ', num2str(sparseDMDc_err)])
%% Show training and prediction
VIZ_SI_Validation

%% Save Data
Model.name = 'sparseDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.dt = dt;
save(fullfile(datapath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')