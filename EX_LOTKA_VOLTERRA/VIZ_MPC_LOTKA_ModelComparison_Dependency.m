% LOTKA-VOLTERRA system

clear all, close all, clc
figpath = '../FIGURES/'; mkdir(figpath)
datapath = '../DATA/';
addpath('../utils');

%% Load Models
ModelCollection = {'DMDc','SINDYc', 'NARX'}; 
Nmodels = length(ModelCollection);
InputSignalType = 'sphs'; % prbs; chirp; noise; sine2; sphs; mixed
TrainAlg = 'trainbr'; %trainlm,trainbr
dep_trainlength = 1;
dep_noise = 0;

% Ntrain_vec = [5:15,20:5:95,100:100:1000];%,1500:500:3000];
% eta_vec = 0.0;
% Nr = 1;

Ntrain_vec = [5:15,20:5:95,100:100:1000,1250,1500,2000,3000];%,1500:500:3000];
eta_vec = 0.05;
Nr = 1;

N_LENGTHS = length(Ntrain_vec);

%% Load data
ResultsALL(N_LENGTHS,1:Nmodels) = struct('x', [], 'u', [], 't', [], 'xref', [], 'J', [], 'eval_time', []);

ModelName = 'DMDc';
datapath1 = [datapath,'EX_LOTKA_Dependencies/',ModelName,'/'];
load(fullfile(datapath1,['EX_LOTKA_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
ResultsALL(:,1) = Results;

ModelName = 'SINDYc';
datapath1 = [datapath,'EX_LOTKA_Dependencies/',ModelName,'/'];
load(fullfile(datapath1,['EX_LOTKA_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
ResultsALL(:,2) = Results;

ModelName = 'NARX';
datapath1 = [datapath,'EX_LOTKA_Dependencies/',ModelName,'/',TrainAlg,'/'];
load(fullfile(datapath1,['EX_LOTKA_MPC_',ModelName,'_',InputSignalType,'_TrainLength','.mat']))
ResultsALL(:,3) = Results;

%% Correct cost
Q = [1,1]; R = 0.5; Ru = 0.5;
for iM = 1:Nmodels
    for iL = 1:N_LENGTHS
        ResultsALL(iL,iM).J = evalObjectiveFCN(ResultsALL(iL,iM).u,ResultsALL(iL,iM).x,ResultsALL(iL,iM).xref,diag(Q),R,Ru);
    end
end

%% Show trajectories
cb_hrzntl = 1;

cmap = jet(N_LENGTHS);
for jModel = 1:Nmodels
    figure, hold on, box on
    for iM = 1:N_LENGTHS
        plot(ResultsALL(iM,jModel).x(1,:),ResultsALL(iM,jModel).x(2,:),'-','Color',cmap(iM,:),'LineWidth',1)
    end
    plot(ResultsALL(iM,jModel).xref(1,1),ResultsALL(iM,jModel).xref(2,1),'ok','MarkerFaceColor','k')
    %     axis equal
    ylim([-20 100]), xlim([0 200])
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
    ylabel('x2')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_State_',ModelCollection{jModel},'.eps']);
end


%% Show cost trajectory
cb_hrzntl = 1;

cmap = jet(N_LENGTHS);
for jModel = 1:Nmodels
    figure, hold on, box on
    for iM = 1:N_LENGTHS
        plot(ResultsALL(end,jModel).t,cumsum(ResultsALL(iM,jModel).J),'-','Color',cmap(iM,:),'LineWidth',1)
    end
    plot(ResultsALL(end,2).t,cumsum(ResultsALL(end,2).J),'--','Color','k','LineWidth',2)
    ylim([-10 10^6]), xlim([200 max(ResultsALL(end,jModel).t)])
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
    print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost_',ModelCollection{jModel},'.eps']);
end


%% Show execution time
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
lstyles = {'-', '--', '-.'};
symbols = {'o', 'd', 's'};

Tx = zeros(length(N_LENGTHS),Nmodels);
for iM = 1:Nmodels
    for iL = 1:N_LENGTHS
        Tx(iL,iM) = ResultsALL(iL,iM).eval_time;
    end
end
figure, hold on, box on
for iM = 1:Nmodels
    plot(Ntrain_vec,Tx(:,iM)./60,symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:))
end
xlim([4 1001])
ylim([0.5*10^-5 10^2])
set(gca,'ytick',[10.^([-5:2:2])])
set(gca,'XScale','log','YScale','log')
xlabel('Length of Training Data');
ylabel('Time [m]')
legend(ModelCollection,'Location','SouthEast')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_ExecTime','.eps']);
legend('off')
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_ExecTime','_noleg.eps']);

%% Add penalization
Q = 1*eye(2); %20
for iM = 1:Nmodels
    for iL = 1:N_LENGTHS
       if ResultsALL(iL,iM).t(end) == 0
           % Flag not converged / ended early
           ResultsALL(iL,iM).J(end) = 10^6;
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

for iM = 1:Nmodels
%     Jend{iM} = zeros(length(N_LENGTHS),Nmodels);
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
    fillh1 = fill([4.5*10^0 2.85*10^3 2.85*10^3 4.5*10^0], [5.5*10^3 5.5*10^3 1.5*10^4 1.5*10^4],0.9*ones(1,3)); % /w noise, not performing well above 
elseif  eta_vec == 0.0
    fillh1 = fill([4.5*10^0 2.85*10^3 2.85*10^3 4.5*10^0], [5.5*10^3 5.5*10^3 2.5*10^4 2.5*10^4],0.9*ones(1,3));
end
fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5; 

[~,idx_best] = min(JendAll(:,2));
plot(Ntrain_vec,JendAll(idx_best,2)*ones(size(Ntrain_vec)),'-','Color',ccolors(2,:),'LineWidth',1)

for iM = 1:Nmodels
%     if iM <= Nmodels
%         plot(xx,yModel(:,iM),'-','Color',ccolors(iM,:), 'LineWidth',1, 'LineStyle',lstyles{iM}), hold on,
%     end
    ph(iM) = plot(Ntrain_vec,JendAll(:,iM),symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:),'MarkerSize',4);
end
if eta_vec == 0
    iM = 1; plot([6.5,6.5],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 2; plot([13.5,13.5],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 3; plot([100,100],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
elseif eta_vec == 0.05
    iM = 1; plot([7.5,7.5],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 2; plot([70,70],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});
    iM = 3; plot([130,130],[5*10^3 10^7],'Color',ccolors(iM,:),'LineStyle',lstyles{iM});    
end
xlim([4 3001])
% ylim([-100 5000])
% ylim([5*10^2 1.5*10^4])
% ylim([1*10^2 10^5])
ylim([5*10^3 10^7])
% set(gca,'YScale','log')
set(gca,'XScale','log','YScale','log', 'ytick', [10^3,10^4,10^5,10^6,10^7], 'xtick', [10^0,10^1,10^2,10^3]) %[10^3,5*10^3,10^4]
xlabel('Length of Training Data');
ylabel('Cost')
lh = legend(ph,ModelCollection,'Location','NorthEast');
lh.Position = [0.71 0.73 lh.Position(3)-0.05 lh.Position(4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','.eps']);

lh.Location = 'NorthOutside';
lh.Orientation = 'horizontal';
lh.Box = 'off';
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','_legout.eps']);
legend('off')
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Cost','_noleg.eps']);

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
for iM = 1:Nmodels
   for iL = 1:N_LENGTHS
        tmp = cumsum(ResultsALL(iL,iM).J);
        Jend(iL,iM) = tmp(end);
   end
   [~,BestModelIDX(iM)] = min(Jend(:,iM));
end

% [val,idx1] = min(abs(Jend(:,1)-9.0601e+03))
% Ntrain_vec(idx1)
% [val,idx2] = min(abs(Jend(:,2)-8.6595e+03))
% Ntrain_vec(idx2)
% [val,idx3] = min(abs(Jend(:,3)-8.9236e+03))
% Ntrain_vec(idx3)

ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
lstyles = {'-', '--', '-.'};

figure,
hold on, box on
plot(ResultsALL(BestModelIDX(1),1).t,ResultsALL(1,1).xref(1,:),'-k')
for iM = 1:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,ResultsALL(BestModelIDX(iM),iM).x(1,:),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
ylim([65 110])
xlabel('Time'), ylabel('Prey')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_x1','.eps']);

figure,
hold on, box on
plot(ResultsALL(BestModelIDX(1),1).t,ResultsALL(1,1).xref(2,:),'-k')
for iM = 1:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,ResultsALL(BestModelIDX(iM),iM).x(2,:),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
ylim([7 27])
xlabel('Time'), ylabel('Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_x2','.eps']);

figure,
hold on, box on
for iM = 1:Nmodels
    plot(ResultsALL(BestModelIDX(iM),iM).t,cumsum(ResultsALL(BestModelIDX(iM),iM).J),'Color',ccolors(iM,:),'LineStyle',lstyles{iM},'LineWidth',2)
end
if eta_vec == 0
    ylim([7*10^3 9.5*10^3])
elseif eta_vec == 0.05
    ylim([7*10^3 1*10^4])
end
set(gca,'yscale','linear','ytick',[7*10^3,8*10^3,9*10^3],'yticklabels',{'70','80','90'})
xlabel('Time'), ylabel('Cost')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_BestModels_J','.eps']);

%% Show State Performance over training lengths
Q = 1*eye(2); %20
Jstate = zeros(length(ResultsALL(iL,iM).J),N_LENGTHS,Nmodels);
for iM = 1:Nmodels
    for iL = 1:N_LENGTHS
        for it = 1:length(ResultsALL(iL,iM).xref(1,:))
            Jstate(it,iL,iM) = abs(ResultsALL(iL,iM).xref(:,it) - ResultsALL(iL,iM).x(:,it))'*Q*abs(ResultsALL(iL,iM).xref(:,it) - ResultsALL(iL,iM).x(:,it));
        end
    end
end


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

for iM = 1:Nmodels
%     Jend{iM} = zeros(length(N_LENGTHS),Nmodels);
    count = 0;
    for iL = 1:N_LENGTHS,
        tmp = cumsum(Jstate(:,iL,iM));
        if tmp(end)<10^10
            count = count + 1;
            xJend{iM}(count) = Ntrain_vec(iL);
            Jend{iM}(count) = tmp(end);
        end
        JendAll(iL,iM) = tmp(end);
    end
    yModel(:,iM) = spline(xJend{iM}(Nstart:end),Jend{iM}(Nstart:end),xx);
end

figure, hold on, box on
for iM = 1:Nmodels
%     if iM < Nmodels
%         plot(xx,yModel(:,iM),'-','Color',ccolors(iM,:), 'LineWidth',1, 'LineStyle',lstyles{iM}), hold on,
%     end
    plot(Ntrain_vec,JendAll(:,iM),symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:),'MarkerSize',4)
end
xlim([4 1001])
ylim([5*10^3 5*10^6])
% set(gca,'YScale','log')
set(gca,'XScale','log','YScale','log', 'ytick', [10^4,10^5,10^6], 'xtick', [10^0,10^1,10^2,10^3])
xlabel('Length of Training Data');
ylabel('Cost')
lh = legend(ModelCollection,'Location','NorthEast');
lh.Position = [0.71 0.73 lh.Position(3)-0.05 lh.Position(4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_CostState','.eps']);

lh.Location = 'NorthOutside';
lh.Orientation = 'horizontal';
lh.Box = 'off';
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_CostState','_legout.eps']);
legend('off')
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_CostState','_noleg.eps']);


return
%% Legend // Dummy plot
figure, hold on, box on
for iM = 1:Nmodels
    plot(Ntrain_vec,ones(size(Ntrain_vec)),'-','Marker',symbols{iM},'Color',ccolors(iM,:),'MarkerFaceColor',ccolors(iM,:))
end

% xlabel('Length of Training Data');
% ylabel('Cost')
axis off
lh = legend(ModelCollection,'Location','NorthEast');
lh.Location = 'NorthOutside';
lh.Orientation = 'horizontal';
lh.Box = 'off';
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_CTRLPERF_',InputSignalType,'_Eta',sprintf('%03g',100*eta_vec),'_TrainingLength_Legend','.eps']);


return
%% TESTING

% [~,idx_best_DMDc] = min(JendAll(:,1));
% % plot(Ntrain_vec,JendAll(idx_best,2)*ones(size(Ntrain_vec)),'-','Color',ccolors(2,:),'LineWidth',1)
% iM = 1;
% iL = idx_best_DMDc; 

% 
iM = 3;
iL = find(Ntrain_vec == 85); % DMDc

% iM = 1;
% iL = find(Ntrain_vec == 10); % DMDc



% iM = 3;
% iL = find(Ntrain_vec == 3000); % NARX: 85, 100, 14

figure, hold on
plot(ResultsALL(end,2).x(1,:),ResultsALL(end,2).x(2,:),'-k','LineWidth',1)
plot(ResultsALL(iL,iM).x(1,:),ResultsALL(iL,iM).x(2,:),'--r','LineWidth',1)

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

%% TESTING
figure,hold on
for iM = 1:Nmodels
    Jend = cumsum(ResultsALL(iM,1).J);
    semilogx(Ntrain_vec(iM),Jend(end),'or')
    Jend = cumsum(ResultsALL(iM,2).J);
    semilogx(Ntrain_vec(iM),Jend(end),'og')
    Jend = cumsum(ResultsALL(iM,3).J);
    semilogx(Ntrain_vec(iM),Jend(end),'ob')
end

%%
figure,hold on
for iM = 1:Nmodels
    Jend = cumsum(ResultsALL(iM,1).J);
    plot3(Ntrain_vec(iM).*ones(size(Jend)),ResultsALL(iM,1).t,Jend,'-r')
end

%%
Jend = zeros(Nmodels,1);
for iM = 1:Nmodels
    tmp = cumsum(ResultsALL(iM,1).J); Jend(iM) = tmp(end);
end

figure,plot(Jend)

figure,plot(ResultsALL(3,1).x(1,:))
figure,plot(ResultsALL(3,1).u)

