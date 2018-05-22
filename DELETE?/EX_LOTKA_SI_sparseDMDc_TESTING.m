% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../../FIGURES/';
datapath = '../../DATA/';
addpath('../utils');

ModelName = 'sparseDMDc';
%% Generate Data
InputSignalType = 'sine2';%prbs; chirp; noise; sine2; sphs; mixed
Ndelay = 1;%35;%35;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

eta = 0;
Hx = Hx + eta*randn(size(Hx));
n = size(Hx,1);
%% DMDc: B = unknown  and with time delay coordinates
numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 2;
r1 = size(Hx,1); r2 = size(Hx,1);
[sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);

%% Sparse DMDc
options_method.sparsify = 'LASSO'; % 'LASSO'ILSTH
options_method.subspace = 0;
options_method.meas_type = 'nonlinear'; % linear, nonlinear
options_method.order  = 4;
options_method.usesine = 0;
epsilon = 1e-1;
lambda  = .0001;%0.001;

% Reorganize measurements
U  = Hu(1:end-1);
switch options_method.meas_type
    case 'linear'
        % Linear / sparseDMDc / doesn't work well even with noise
        Y  = Hx(:,1:end-1);
        Yp = Hx(:,2:end);
    case 'nonlinear'
        % Nonlinear / sparseEDMDc / converges to DMDc with increasing noise (lambda>0 doesn't make it better, but worse)
        Y = poolData(Hx(:,1:end-1)',n,options_method.order,options_method.usesine)';
        Yp = poolData(Hx(:,2:end)',n,options_method.order,options_method.usesine)';
%         Y  = [Hx(1,1:end-1); Hx(2,1:end-1); Hx(1,1:end-1).^2; Hx(1,1:end-1).*Hx(2,1:end-1); Hx(2,1:end-1).^2];
%         Yp = [Hx(1,2:end); Hx(2,2:end); Hx(1,2:end).^2; Hx(1,2:end).*Hx(2,2:end); Hx(2,2:end).^2];
end

Nstates = size(Y,1);
Ninputs = size(U,1);

% Sparse regression 
switch options_method.sparsify
    case 'LASSO'
        % using the l1 norm to promote sparsity
        G = zeros(Nstates+Ninputs,Nstates);
        for i = 1:Nstates
            [G(:,i), FitInfo] = lasso([Y;U]',Yp(i,:)','Lambda',lambda); %,'Lambda',0.00000001
            %     lassoPlot(G{i},FitInfo);
        end
    case 'ILSTH'
        G = sparsifyDynamics([Y;U]',Yp',lambda,n);
end

G = G';
A = G(1:Nstates,1:Nstates)
B = G(1:Nstates,Nstates+1:Nstates+Ninputs)
C = eye(Nstates,Nstates);
D = zeros(Nstates,1);
sysmodel_sparseDMDc = ss(A,B,C,D,dt);

%% DEBUG
% % Nonlinear / EDMDc / identical to sparsifyDynamics with lambda=0
% Y  = [Hx(1,:); Hx(2,:); Hx(1,:).^2; Hx(1,:).*Hx(2,:); Hx(2,:).^2]; % 2
% % Y  = [Hx(1,:); Hx(2,:); Hx(1,:).^2; Hx(1,:).*Hx(2,:); Hx(1,:).^2]; % 3
% U  = Hu;
% r1 = size(Hx,1); r2 = size(Hx,1);
% [sysmodel_sparseDMDc,~,~] = DelayDMDc_MV(Y,U,size(Y,1),size(Y,1),dt,size(Y,1),size(U,1),2);

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

% DELETE?
% if options_method.subspace == 1
%     Gamma = U;
%     Omega = Y; %[Y; Gamma];
%     [m,n] = size(Omega);
%     [U,S,V] = svd(Omega,'econ');
%     [threshold, r1] = getSVDThreshold(diag(S),m,n);
% 	Ur = U(:,1:r1);Sr = S(1:r1,1:r1);Vr = V(:,1:r1);
%     Omega = Ur*Sr*Vr';
%     Y = Omega(1:Nstates,:); U = Gamma;
%     %U = Omega(Nstates+1:Nstates+Ninputs,:);
% 	%G = Yp*Vr*Sr^(-1)*Ur';
% end
%% Prediction over training phase
[xDMDc,~] = lsim(sysmodel_DMDc,Hu,tspan(1:Nt),Hx(:,1));
xDMDc = xDMDc(:,end-1:end);
xDMDc = xDMDc + repmat(xmean,[Nt 1]);

if Nstates == 2
    [xDMDc_sparse,~] = lsim(sysmodel_sparseDMDc,Hu,tspan(1:Nt),Hx(:,1));
    xDMDc_sparse = xDMDc_sparse(:,end-1:end);
    xDMDc_sparse = xDMDc_sparse + repmat(xmean,[Nt 1]);
elseif Nstates > 2
    %x0dmdc = [Hx(1,1); Hx(2,1); Hx(1,1).^2; Hx(1,1).*Hx(2,1); Hx(2,1).^2];
    x0dmdc = poolData(Hx(:,1)',n,options_method.order,options_method.usesine)';
    [xDMDc_sparse,~] = lsim(sysmodel_sparseDMDc,Hu,tspan(1:Nt),x0dmdc);
    xDMDc_sparse = xDMDc_sparse(:,2:3);
    xDMDc_sparse = xDMDc_sparse + repmat(xmean,[Nt 1]);
end
% Ns = 2;
% xDMDc = zeros(2*Ns,Nt);
% xDMDc(:,1) = Hx(:,1);
% for ct=1:Nt-1
%     % Obtain plant state at next prediction step.
%     xDMDc(:,ct+1) = sysmodel_DMDc.A*xDMDc(:,ct) + sysmodel_DMDc.B*Hu(:,ct); %[Hx(1:2,ct); xDMDc(3:4,ct)]
% end
% xDMDc = xDMDc(3:4,:);
% xDMDc = xDMDc + repmat(xmean',[1 Nt]);
% xDMDc = xDMDc';

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
% xlim([0 100]), ylim([0,120]);
xlabel('Time')
ylabel('Population size')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,3,5]),'True','DMDc',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.eps']);


DMDc_err = norm(x(1:end-1,:)-xDMDc,2);
sparseDMDc_err = norm(x(1:end-1,:)-xDMDc_sparse,2);
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
%     x0dmdc = [x0(1)-xmean(1); x0(2)-xmean(2); (x0(1)-xmean(1)).^2; (x0(1)-xmean(1)).*(x0(2)-xmean(2)); (x0(2)-xmean(2)).^2];
    x0dmdc = poolData([x0(1)-xmean(1), x0(2)-xmean(2)],n,options_method.order,options_method.usesine)';
    [xB_sparse,tB] = lsim(sysmodel_sparseDMDc,Hunew,tspanV,x0dmdc);
    xB_sparse = xB_sparse(:,2:3);
    xB_sparse = xB_sparse + repmat(xmean,[length(tB) 1]);
end


DMDc_err = norm(xA-xB(1:end-1,:),2);
sparseDMDc_err = norm(xA-xB_sparse(1:end-1,:),2);
disp(['Validation Error: '])
disp(['                DMDc: ', num2str(DMDc_err)])
disp(['          sparseDMDc: ', num2str(sparseDMDc_err)])
%% Show training and prediction
VIZ_SI_Validation

return
%% Save Data
Model.name = 'sparseDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.dt = dt;
save(fullfile(datapath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')