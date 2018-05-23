% HIV system
% System identification: SINDYc

clear all, close all, clc
figpath = '../FIGURES/HIV/'; mkdir(figpath)
datapath = '../DATA/HIV/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'HIV';
ModelName = 'partialSINDYc'; %select: SINDYc, partialSINDYc

%% Generate Data 
InputSignalType = 'prbs';
ONLY_TRAINING_LENGTH = 1;
DATA_ENSEMBLE = 0;
getTrainingData

%% SINDYc
% Parameters
polyorder = 3;
usesine = 0;
Nvar = 5;
if strcmp(ModelName,'SINDYc')
    SelectVars = [1:Nvar];
    NSindy = Nvar;
    lambda_vec = [1e1,3.1,3e0,1e-1,5e-1]; % dt = 1/12 % in paper
elseif strcmp(ModelName,'partialSINDYc')
    SelectVars = [1,2,3];
    NSindy = length(SelectVars);
    lambda_vec = [1e1,30,3e0,1e-1]; % dt = 1/12, 1:3, in paper
end

trainSINDYc_HIV

%% Truth
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

%% Show training data
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

%% Show comparison of estimated coefficients with true coefficients
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

%% Prediction over training phase
tTRUE = t;
xTRUE = x;
options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,NSindy));
[tSINDYc,xSINDYc]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:NSindy),polyorder,usesine),tspan,x0(SelectVars),options);  % approximate

%% Show prediction over training phase
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
ylim([0 10])
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