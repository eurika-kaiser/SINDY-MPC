% Automatic flight control of F8
% System identification: SINDYc

clear all, close all, clc
figpath = '../../FIGURES/HIV/'; mkdir(figpath)
datapath = '../../DATA/HIV/'; mkdir(datapath)
addpath('../utils');
% addpath('../../../../Documents/Academia/Sources/Code/psv-master')

% ReferenceModel = 'ZURAKOWSKI';
SystemModel = 'HIV';
ModelName = 'partialSINDYc'; %SINDYc,partialSINDYc

%% Generate Data %{'sine2', 'chirp','prbs', 'sphs'}
InputSignalType = 'prbs';%'sine2';%'unforced','sine2'; %prbs; noise; sine2; sphs;
ONLY_TRAINING_LENGTH = 1;
getTrainingData
% getTrainingData_Ensemble

%%
clear ph
figure,
semilogy(tspan,repmat(xSTEADY1,[length(tspan) 1]),'--k','LineWidth',1), hold on
% semilogy(tspan,repmat(xSTEADY2,[length(tspan) 1]),'--b','LineWidth',1)
for i = 1:5
    semilogy(tspan,squeeze(x(:,i,:)),'-','LineWidth',1);hold on
    %     ph(i)=semilogy(tspan,x(:,i,1),'-','LineWidth',1);hold on
end
% legend(ph,'x1','x2','x3','x4','x5','x6')

%% SINDYc
% Parameters
polyorder = 3;
usesine = 0;
Nvar = 5;
if strcmp(ModelName,'SINDYc')
    SelectVars = [1:Nvar];
    NSindy = Nvar;
elseif strcmp(ModelName,'partialSINDYc')
    SelectVars = [1,2];%[1,2,3];
    NSindy = length(SelectVars);
end
% lambda = 1e-1;
% lambda_vec = [1e1,1e0,1e-1,1e-1,1e-1]; % dt = 0.01
% lambda_vec = [1e1,3.1,3e0,1e-1,5e-1]; % dt = 1/12 % in paper
% lambda_vec = [1e1,30,3e0,1e-1]; % dt = 1/12, 1:3, was good
lambda_vec = [1e1,2.5e1]; % dt = 1/12, 1:2, bad
% lambda_vec = [5.e1,2e2,1e0,1e-1]; % dt = 1/12, 1,3,4, bad
% lambda_vec = [2.7e2,2e3,3e2,1e-1]; % dt = 1/12, 1,3,5, bad

% lambda2 = [1e-1.^[0.05:0.05:1],1:0.05:1e1,1e1:1:1e2];
% lambda2 = sort(unique(lambda2));
% xi2 = zeros(84,length(lambda2));
% for iL = 1:length(lambda2)
% lambda_vec = [1e1,lambda2(iL),3e0,1e-1,5e-1]; % dt = 1/12


%Add noise
%eps = 0.05*std(x(:,1));
%x = x + eps*randn(size(x));

trainSINDYc_HIV
% xi2(:,iL) = Xi(:,2);
% end
%%
% figure,hold all
% plot(Xi0(:,2),'-k','LineWidth',3)
% for i = 1:length(lambda2)
% plot(xi2(:,i))
% end

% err = zeros(size(xi2,1),1);
% for i = 1:length(lambda2)
%     err(i) = sum((xi2(:,i)-Xi0(:,2)).^2);
% end
%
% figure,semilogy(err)
% lambda_vec = [3e2,1e3, 1e0, 1e1, 1e0, 1e2];
% trainSINDYc_Ensemble

%% Truth
% KWON
% xi_truth.term{1} = {'1', 'x1', 'x1x4', 'x1x5', 'x1x4u'};
% xi_truth.coeff{1} = [lambda1, -d, -k1, -k2, k1];
% xi_truth.term{2} = {'x2', 'x1x4', 'x1x4u'};
% xi_truth.coeff{2} = [-delta, k1*(1-mu), -k1*(1-mu)];
% xi_truth.term{3} = {'x3', 'x1x4', 'x1x5', 'x1x4u'};
% xi_truth.coeff{3} = [-delta, k1*mu, k2, -k1*mu];
% xi_truth.term{4} = {'x2', 'x4'};
% xi_truth.coeff{4} = [alpha1, -c];
% xi_truth.term{5} = {'x3', 'x5'};
% xi_truth.coeff{5} = [alpha2, -c];

% ZURAKOWSKI
xi_truth.term{1} = {'1', 'x1', 'x1x2', 'x1x2u'};
xi_truth.coeff{1} = [lambda1, -d, -alpha1, eta*alpha1];
xi_truth.term{2} = {'x2', 'x1x2', 'x2x4', 'x2x5', 'x1x2u'};
xi_truth.coeff{2} = [-a, alpha1, -p1, -p2, -eta*alpha1];
xi_truth.term{3} = {'x3', 'x2x3', 'x1x2x3'};
xi_truth.coeff{3} = [-b2, -c2*q, c2];
xi_truth.term{4} = {'x4', 'x2x4'};
xi_truth.coeff{4} = [-b1, c1];
xi_truth.term{5} = {'x5', 'x2x3'};
xi_truth.coeff{5} = [-h, c2*q];

Theta0 = poolData([x,u],Nvar+1,polyorder,usesine);
% Construct true Xi
Xi0 = zeros(size(Theta0,2),Nvar+1);
for i = 1:Nvar
    for j = 1:length(xi_truth.term{i})
        idx = find(strcmp(yout,xi_truth.term{i}(j)));
        Xi0(idx,i) = xi_truth.coeff{i}(j);
    end
end

% Estimate Xi with known components
if NSindy == Nvar
    if exist('DX') == 0
        DX = dx;
    end
    Xi_est = zeros(size(Xi0));
    for ind = 1:Nvar
        idx = find(Xi0(:,ind)~=0);
        Xi_est(idx,ind) = Theta(:,idx)\DX(:,ind);
    end
    
    norm(Xi-Xi0,'fro')
end

return
%%
figure; box on
plot(t,squeeze(x(:,1,:)),'LineWidth',1.5)
ylabel('x1'),xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_x1.eps']);


figure; box on
plot(t,squeeze(x(:,2,:)),'LineWidth',1.5)
ylabel('x2'),xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_x2.eps']);


figure; box on
plot(t,squeeze(x(:,3,:)),'LineWidth',1.5)
ylabel('x3'),xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_x3.eps']);


figure; box on
plot(t,squeeze(x(:,4,:)),'LineWidth',1.5)
ylabel('x4')
xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_x4.eps']);

figure; box on
plot(t,squeeze(x(:,5,:)),'LineWidth',1.5)
ylabel('x5')
xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_x5.eps']);

figure; box on
plot(t,squeeze(u(:,1,1)),'LineWidth',1.5)%[1,3,4,15]
ylabel('u')
xlabel('Time / days')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
axis tight
ylim([0 0.8])
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_Training_u.eps']);
%%
figure,box on,hold on
cmap = jet(Nvar*2); cmap = cmap(1:2:end,:);
cmap_dark = cmap + [0.2,0.2,0.2; 0.2, 0.2,0; 0.2,0,0; 0.2,0,0.2; 0,0.2,0.2];
cmap_dark(1:5,:) = [0,0,1; 0.1,0.7,0.1; 0.8,0,1; 1,0,0; 0.8,0.7,0.3];

cmap_light = zeros(size(cmap));
cmap_light(1,:) = cmap(2,:);
cmap_light(2,:) = cmap(4,:);
cmap_light(5,:) = cmap(5,:);
cmap_light(3,:) = cmap_dark(3,:) + [0.1,0.1,0];
cmap_light(4,:) = cmap_dark(4,:) + [0,0.2,0.2];
for i = 1:NSindy
    plot(Xi0(:,i),'-','Color',cmap_dark(i,:),'LineWidth',2), hold on
    plot(Xi(:,i),'--','Color',cmap_light(i,:),'LineWidth',2)
end

%% Prediction  over training phase
if exist('Nic')==1
    forcing = @(x,t) [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.05*t)))+0.3];
    x0 = [lambda1/d, 0, 0, 0.3, 0];
    [tTRUE,xTRUE]=ode45(@(t,x) HIVsys_KWON(t,x,forcing(x,t)),tspan,x0,options);
else
    tTRUE = t;
    xTRUE = x;
end

options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,NSindy));
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tSINDYc,xSINDYc]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:NSindy),polyorder,usesine),tspan,x0(SelectVars),options);  % approximate
else
    p.ahat = Xi(:,1:NSindy);
    p.polyorder = polyorder;
    p.usesine = usesine;
    p.dt = dt;
    [N,Ns] = size(x);
    xSINDYc = zeros(NSindy,N); xSINDYc(:,1) = x0(SelectVars)';
    for ct=1:N-1
        xSINDYc(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xSINDYc(:,ct),u(ct),dt,1,[],p);
    end
    xSINDYc = xSINDYc';%(:,2:N+1)';
end

%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
    ccolors(2,:)-[0.1 0.2 0.09];
    ccolors(3,:)-[0.1 0.2 0.09];
    ccolors(4,:)-[0.1 0.1 0.2];
    ccolors(5,:)-[0.1 0.2 0.09]];
for i = 1:NSindy
    ph(i) = plot(tTRUE,xTRUE(:,SelectVars(i)),'-','Color',ccolors(SelectVars(i),:),'LineWidth',1); hold on
end
for i = 1:NSindy
    ph(NSindy+i) = plot(tSINDYc,xSINDYc(:,i),'--','Color',ccolors_valid(SelectVars(i),:),'LineWidth',2);
end
xlabel('Time')
ylabel('xi')
legend(ph([1,NSindy]),'True',ModelName,'Location','SouthEast')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);

%% Save Data
Model.name = 'SINDYc';
Model.polyorder = polyorder;
Model.usesine = usesine;
Model.Xi = Xi;
Model.dt = dt;
Model.SelectVars = SelectVars;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_',num2str(SelectVars')','.mat']),'Model')