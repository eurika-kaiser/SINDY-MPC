% LORENZ system

clear all, close all, clc

WORKING = 1;
SystemModel = 'LORENZ';

figpath = '../FIGURES/LORENZ/'; mkdir(figpath)
datapath = '../DATA/';
addpath('../utils');

%% Load Models
ModelCollection = {'DMDc','SINDYc', 'NARX'}; 
Nmodels = length(ModelCollection);
InputSignalType = 'sine2'; % prbs; chirp; noise; sine2; sphs; mixed
TrainAlg = 'trainbr'; %trainlm,trainbr
dep_trainlength = 1;
dep_noise = 0;

%% Select case
% 1) Dependency on training length, zero noise level
% Ntrain_vec = [5:15,20:5:95,100:100:1000];%,1500:500:3000];
% eta_vec = 0.0;
% Nr = 1;

% 2) Dependency on training length for fixed noise level
Ntrain_vec = [5:15,20:5:95,100:100:1000,1250,1500,2000,3000];%,1500:500:3000];
eta_vec = 0.05;
Nr = 1;

N_LENGTHS = length(Ntrain_vec);

%% Load data
ResultsALL(N_LENGTHS,1:Nmodels) = struct('x', [], 'u', [], 't', [], 'xref', [], 'J', [], 'eval_time', []);

% ModelName = 'DMDc';
% datapath1 = [datapath,'EX_',SystemModel,'_Dependencies/',ModelName,'/'];
% load(fullfile(datapath1,['EX_',SystemModel,'_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
% ResultsALL(:,1) = Results;

ModelName = 'SINDYc';
datapath1 = [datapath,'EX_',SystemModel,'_Dependencies/',ModelName,'/TL/'];
load(fullfile(datapath1,['EX_',SystemModel,'_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
ResultsALL(:,2) = Results;

ModelName = 'NARX';
datapath1 = [datapath,'EX_',SystemModel,'_Dependencies/',ModelName,'/',TrainAlg,'/TL/']; %%%% !!!!!!!!!! TL
load(fullfile(datapath1,['EX_',SystemModel,'_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
ResultsALL(:,3) = Results;


%% Show trajectories
cb_hrzntl = 1;

cmap = jet(N_LENGTHS);
for jModel = 2:Nmodels
    figure, hold on, box on
    for iM = 1:N_LENGTHS
        plot3(ResultsALL(iM,jModel).x(1,:),ResultsALL(iM,jModel).x(2,:),ResultsALL(iM,jModel).x(3,:),'-','Color',cmap(iM,:),'LineWidth',1)
    end
    plot3(ResultsALL(iM,jModel).xref(1,1),ResultsALL(iM,jModel).xref(2,1),ResultsALL(iM,jModel).xref(3,1),'ok','MarkerFaceColor','k')
    view(30,20)
    colormap(cmap),
    ch = colorbar;
    ch.Limits = [0 1];
    ch.Ticks = [0,0.25,0.5,0.75,1];
    ch.TickLabels = num2str( Ntrain_vec( [1, (N_LENGTHS-1)/4, (N_LENGTHS-1)/2+1, 3*(N_LENGTHS-1)/4+1, 4*(N_LENGTHS-1)/4+1] )' );
    if cb_hrzntl == 0
        ch.Position = [ch.Position(1),ch.Position(2)+0.14,ch.Position(3)-0.01,ch.Position(4)-0.1];
        set(gca,'Position',[0.08 0.2753 0.775 0.6497])
    elseif cb_hrzntl == 1
        ch.Location = 'southoutside';
        ch.Position = [0.2,0.1,0.75,ch.Position(4)];
        set(gca,'Position',[0.2 0.4 0.75 0.55])
    end
    xlabel('x1');
    ylabel('x2'); 
    zlabel('x3')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_State_',ModelCollection{jModel},'.eps']);
end

%% Show cost trajectory
cb_hrzntl = 1;

cmap = jet(N_LENGTHS);
for jModel = 2:Nmodels
    figure, hold on, box on
    for iM = 1:N_LENGTHS
        plot(ResultsALL(end,jModel).t,cumsum(ResultsALL(iM,jModel).J),'-','Color',cmap(iM,:),'LineWidth',1)
    end
    plot(ResultsALL(end,2).t,cumsum(ResultsALL(end,2).J),'--','Color','k','LineWidth',2)
    ylim([-10 10^6]), xlim([20 max(ResultsALL(end,jModel).t)])
    colormap(cmap),
    ch = colorbar;
    ch.Limits = [0 1];
    ch.Ticks = [0,0.25,0.5,0.75,1];
    ch.TickLabels = num2str( Ntrain_vec( [1, (N_LENGTHS-1)/4, (N_LENGTHS-1)/2+1, 3*(N_LENGTHS-1)/4+1, 4*(N_LENGTHS-1)/4+1] )' );
    if cb_hrzntl == 0
        ch.Position = [ch.Position(1),ch.Position(2)+0.14,ch.Position(3)-0.01,ch.Position(4)-0.1];
        set(gca,'Position',[0.08 0.2753 0.775 0.6497])
    elseif cb_hrzntl == 1
        ch.Location = 'southoutside';
        ch.Position = [0.2,0.1,0.75,ch.Position(4)];
        set(gca,'Position',[0.2 0.4 0.75 0.55])
    end
    ylim([10^3 10^6])
    set(gca,'yscale','log')
    xlabel('Time');
    ylabel('Cost')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost_',ModelCollection{jModel},'.eps']);
end


%% Show execution time
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
lstyles = {'-', '--', '-.'};
symbols = {'o', 'd', 's'};

Tx = zeros(length(N_LENGTHS),Nmodels);
for iM = 2:Nmodels
    for iL = 1:N_LENGTHS
        Tx(iL,iM) = ResultsALL(iL,iM).eval_time;
    end
end
figure, hold on, box on
for iM = 2:Nmodels
    plot(Ntrain_vec,Tx(:,iM)./60,symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:))
end
xlim([4 3001])
ylim([0.5*10^-5 10^2])
set(gca,'xtick',[10.^([0:1:3])])
set(gca,'ytick',[10.^([-5:2:2])])
set(gca,'XScale','log','YScale','log')
xlabel('Length of Training Data');
ylabel('Time [m]')
legend(ModelCollection{2:3},'Location','SouthEast')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_ExecTime','.eps']);
legend('off')
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_ExecTime','_noleg.eps']);

%% Flag unconverged results
for iM = 2:Nmodels
    for iL = 1:N_LENGTHS
       if ResultsALL(iL,iM).t(end) == 0
           % Flag not converged / ended early
           ResultsALL(iL,iM).J(end) = 10^5;
       else
           % Flag if far from end state at end of time
           %ResultsALL(iL,iM).J(end) = ResultsALL(iL,iM).J(end) + abs(ResultsALL(iL,iM).xref(:,end) - ResultsALL(iL,iM).x(:,end))'*Q*abs(ResultsALL(iL,iM).xref(:,end) - ResultsALL(iL,iM).x(:,end));
       end
    end
end

%% Show Performance over training lengths
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
lstyles = {'-', '--', '-.'};
symbols = {'o', 'd', 's'};

Nstart = 20;
xx = Ntrain_vec(Nstart):2:1000;
yModel = zeros(length(xx),Nmodels);
JendAll = zeros(length(N_LENGTHS),Nmodels);
Jend = cell(Nmodels,1);
xJend = cell(Nmodels,1);

for iM = 2:Nmodels
    count = 0;
    for iL = 1:N_LENGTHS,
        tmp = cumsum(ResultsALL(iL,iM).J);
        if tmp(end)<10^10
            count = count + 1;
            xJend{iM}(count) = Ntrain_vec(iL);
            Jend{iM}(count) = tmp(end);
        end
        JendAll(iL,iM) = tmp(end);
    end
    yModel(:,iM) = spline(xJend{iM}(Nstart:end),Jend{iM}(Nstart:end),xx);
end

clear ph
figure, hold on, box on

if eta_vec == 0.05
    fillh1 = fill([4.5*10^0 2.85*10^3 2.85*10^3 4.5*10^0], [2*10^4 2*10^4 6*10^4 6*10^4],0.9*ones(1,3)); % /w noise, not performing well above 
elseif  eta_vec == 0.0
    fillh1 = fill([4.5*10^0 2.85*10^3 2.85*10^3 4.5*10^0], [2*10^4 2*10^4 6*10^4 6*10^4],0.9*ones(1,3));
end
fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5; 

[~,idx_best] = min(JendAll(:,2));
plot(Ntrain_vec,JendAll(idx_best,2)*ones(size(Ntrain_vec)),'-','Color',ccolors(2,:),'LineWidth',1)

for iM = 2:Nmodels
    ph(iM) = plot(Ntrain_vec,JendAll(:,iM),symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:),'MarkerSize',4);
end
ph(1) = [];
if eta_vec == 0
%    iM = 1; plot([6.5,6.5],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 2; plot([13.5,13.5],[1*10^4 10^6],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 3; plot([100,100],[1*10^4 10^6],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
elseif eta_vec == 0.05
%    iM = 1; plot([7.5,7.5],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 2; plot([75,75],[1*10^4 10^6],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
%     iM = 3; plot([130,130],[1*10^4 10^6],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});    
end
xlim([4 3001])
ylim([1.8*10^4 5*10^5])
set(gca,'XScale','log','YScale','log', 'ytick', [10^4,10^5,10^7], 'xtick', [10^0,10^1,10^2,10^3]) %[10^3,5*10^3,10^4]
xlabel('Length of Training Data');
ylabel('Cost')
lh = legend(ph,ModelCollection{2:3},'Location','NorthEast');
lh.Position = [0.71 0.73 lh.Position(3)-0.05 lh.Position(4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','.eps']);

lh.Location = 'NorthOutside';
lh.Orientation = 'horizontal';
lh.Box = 'off';
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','_legout.eps']);
legend('off')
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','_noleg.eps']);

return
%% Show best case for all models
% eta = 0.05: 8    38    21 
% -->  Ntrain = 12        1250          65
% --> Jend = 9.0603e+03,  8.6549e+03, 2.0592e+04 
% Jend for best eta==0: 8(12), 38(1250), 21(65) 
% SINDY converged Ntrain = 200
% eta = 0.0: 12    25    37
% -->  Ntrain = 20          85        1000
% --> Jend = 9.0601e+03,  8.6595e+03, 8.9236e+03 
% SINDY converged Ntrain = 200
BestModelIDX = zeros(1,Nmodels);
Jend = zeros(N_LENGTHS,Nmodels);
for iM = 2:Nmodels
   for iL = 1:N_LENGTHS
        tmp = cumsum(ResultsALL(iL,iM).J);
        Jend(iL,iM) = tmp(end);
   end
   [~,BestModelIDX(iM)] = min(Jend(:,iM));
end

ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
lstyles = {'-', '--', '-.'};

figure,
hold on, box on
plot(ResultsALL(BestModelIDX(2),2).t,ResultsALL(1,2).xref(1,:),'-k')
for iM = 2:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,ResultsALL(BestModelIDX(iM),iM).x(1,:),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
ylim([-15 5])
xlabel('Time'), ylabel('x1')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_x1','.eps']);

figure,
hold on, box on
plot(ResultsALL(BestModelIDX(2),2).t,ResultsALL(1,2).xref(2,:),'-k')
for iM = 2:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,ResultsALL(BestModelIDX(iM),iM).x(2,:),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
ylim([-25 5])
xlabel('Time'), ylabel('x2')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_x2','.eps']);

figure,
hold on, box on
plot(ResultsALL(BestModelIDX(2),2).t,ResultsALL(1,2).xref(3,:),'-k')
for iM = 2:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,ResultsALL(BestModelIDX(iM),iM).x(3,:),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
ylim([0 45])
xlabel('Time'), ylabel('x3')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_x3','.eps']);

figure,
hold on, box on
for iM = 2:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,cumsum(ResultsALL(BestModelIDX(iM),iM).J),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
if eta_vec == 0
    ylim([7*10^3 9.5*10^3])
elseif eta_vec == 0.05
    ylim([1*10^4 3*10^4])
end
xlabel('Time'), ylabel('Cost')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_J','.eps']);


%% Legend // Dummy plot
figure, hold on, box on
for iM = 2:Nmodels
    plot(Ntrain_vec,ones(size(Ntrain_vec)),'-','Marker',symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:))
end

axis off
lh = legend(ModelCollection{2:3},'Location','NorthEast');
lh.Location = 'NorthOutside';
lh.Orientation = 'horizontal';
lh.Box = 'off';
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Legend','.eps']);


return
%% TESTING

% [~,idx_best_DMDc] = min(JendAll(:,1));
% % plot(Ntrain_vec,JendAll(idx_best,2)*ones(size(Ntrain_vec)),'-','Color',ccolors(2,:),'LineWidth',1)
% iM = 1;
% iL = idx_best_DMDc; 

% 
iM = 2;
iL = find(Ntrain_vec == 8); % DMDc

% iM = 1;
% iL = find(Ntrain_vec == 10); % DMDc



% iM = 3;
% iL = find(Ntrain_vec == 3000); % NARX: 85, 100, 14

figure, hold on
plot3(ResultsALL(end,2).x(1,:),ResultsALL(end,2).x(2,:),ResultsALL(end,2).x(3,:),'-k','LineWidth',1)
plot3(ResultsALL(iL,iM).x(1,:),ResultsALL(iL,iM).x(2,:),ResultsALL(iL,iM).x(3,:),'--r','LineWidth',1)

figure, hold on
plot(ResultsALL(end,2).t,ResultsALL(end,2).x(2,:),'-k','LineWidth',1)
plot(ResultsALL(end,iM).t,ResultsALL(iL,iM).x(2,:),'--r','LineWidth',1)
plot(ResultsALL(end,2).t,ResultsALL(end,2).x(1,:),'-k','LineWidth',1)
plot(ResultsALL(end,iM).t,ResultsALL(iL,iM).x(1,:),'--r','LineWidth',1)

figure, hold on
plot(ResultsALL(end,2).t,cumsum(ResultsALL(end,2).J),'-k','LineWidth',1)
plot(ResultsALL(end,iM).t,cumsum(ResultsALL(iL,iM).J),'--r','LineWidth',1)

figure, hold on
plot(ResultsALL(end,2).t,ResultsALL(end,2).u,'-k','LineWidth',1)
plot(ResultsALL(end,iM).t,ResultsALL(iL,iM).u,'--r','LineWidth',1)

