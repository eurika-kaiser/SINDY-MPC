clear all, close all, clc
figpath = '../FIGURES/'; mkdir(figpath)
datapath = '../DATA/EX_LOTKA_Dependencies/';
addpath('../utils');

%% Paramaters
NARXtraining = 'trainbr'; %'trainlm',trainbr
InputSignalType = 'sphs'; %sphs,sine2
Ntrain_vec = 3000; %1000;
N_LENGTHS = length(Ntrain_vec);

eta_vec = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
N_ETA = length(eta_vec);

Nmodels = N_ETA*N_LENGTHS;

Nr = 50;
Nt = 1000;
%% Load all models
ResultsALL(1:Nmodels,3) = struct('err', zeros(Nr,1), 'errM', zeros(Nr,1), 'xA', zeros(Nt,2), 'xB', zeros(Nt,2,Nr),'Ttraining',zeros(Nr,1));


count = 0;
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        count = count + 1;
        
        ModelName = 'DMDc';
        datapath1 = [datapath,'DMDc/'];
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
        load(filename);
        ResultsALL(count,1) = Results;
        
        ModelName = 'SINDYc';
        datapath1 = [datapath,'SINDYc/'];
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
        load(filename);
        ResultsALL(count,2) = Results;
        
        ModelName = 'NARX';
        datapath1 = [datapath,'NARX/',NARXtraining,'/'];
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
        load(filename);
        ResultsALL(count,3) = Results;
    end
end



%%
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
symbols = {'o', 'd', 's'};
xx = 1:0.1:11;%0:0.01:0.6;
yy = zeros(3,length(xx));
err = zeros(Nr,N_ETA,3);
for iM = 1:3
    for i = 1:N_ETA, err(1:Nr,i,iM) = ResultsALL(i,iM).err; end
    if iM == 2
        xm = zeros(1,N_ETA);
        for k = 1:N_ETA
            TF = isnan(err(1:Nr,k,iM));
            xm(k) = median(err(TF==0,k,iM),1);
        end
    else
        xm = median(err(1:Nr,:,iM),1);
    end
    yy(iM,:) = spline(1:11,xm,xx);
end
figure; hold on
data = err(1:Nr,eta_vec==0.4,2);
ylimregion = 300;
fillh1 = fill([0.1 11.9 11.9 0.1], [10 10 ylimregion.*ones(1,2)],0.9*ones(1,3));
fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5;
ph(1)=plot(xx,yy(1,:),'-r','LineWidth',1); hold on,
ph(2)=plot(xx(1:end),yy(2,1:end),'-g','LineWidth',1);
ph(3)=plot(xx(1:end),yy(3,1:end),'-','Color',[0.7,0.7,1],'LineWidth',1);
for iM = 1:3
    boxplot(err(:,:,iM),'PlotStyle','compact','Colors',ccolors(iM,:), 'Symbol', symbols{iM}, ... % eta_vec
        'Widths',0.1); hold on
    delete(findobj(gca,'Type','text'))
end
set(gca,'xtick',[1:2:11],'xticklabel',num2str(eta_vec(1:2:end)'),'Position',[0.2 0.22 0.75 0.73])
xlim([0 12]), ylim([0 1300])
xt = xlabel('Eta'); xt.Position = [115 -20 -0.1];
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'_noleg.eps']);

l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthWest')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'.eps']);


%%
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
symbols = {'o', 'd', 's'};
xx = 1:0.1:11;%0:0.01:0.6;
yy = zeros(3,length(xx));
err = zeros(Nr,N_ETA,3);
for iM = 1:3
    for i = 1:N_ETA, err(1:Nr,i,iM) = ResultsALL(i,iM).errM; end
    if iM == 2
        xm = zeros(1,N_ETA);
        for k = 1:N_ETA
            TF = isnan(err(1:Nr,k,iM));
            xm(k) = median(err(TF==0,k,iM),1);
        end
    else
        xm = median(err(1:Nr,:,iM),1);
    end
    yy(iM,:) = spline(1:11,xm,xx); %eta_vec(:)
end
figure; hold on
fillh1 = fill([0.1 11.9 11.9 0.1], [10 10 ylimregion.*ones(1,2)],0.9*ones(1,3));
fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5;
ph(1)=plot(xx,yy(1,:),'-r','LineWidth',1); hold on,
ph(2)=plot(xx(1:end),yy(2,1:end),'-g','LineWidth',1);
ph(3)=plot(xx(1:end),yy(3,1:end),'-','Color',[0.7,0.7,1],'LineWidth',1);
% plot([0,12],max(err(1:Nr,eta_vec==0.3,2)).*ones(1,2),'-k')
for iM = 1:3
    boxplot(err(:,:,iM),'PlotStyle','compact','Colors',ccolors(iM,:), 'Symbol', symbols{iM}, ... %, 'Labels', num2str(eta_vec'), eta_vec
        'Widths',0.1); hold on
    delete(findobj(gca,'Type','text'))
end
set(gca,'xtick',[1:2:11],'xticklabel',num2str(eta_vec(1:2:end)'),'Position',[0.2 0.22 0.75 0.73])
xlim([0 12]), ylim([0 1300])
xt = xlabel('Eta'); xt.Position = [115 -20 -0.1];
ylabel('MSE')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noiseM_',NARXtraining,'_noleg.eps']);

l1 = legend(ph,'DMDc','SINDYc','NARX');
set(l1,'Location','NorthWest')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noiseM_',NARXtraining,'.eps']);

%% Normalized error
clear LegLoc
LOG_SCALE = 0;
data2plot = zeros(Nr,N_ETA,3);
% for iN = 1:N_LENGTHS
for iNoise = 1:N_ETA
    for iR = 1:Nr
        for iM = 1:3
            data2plot(iR,iNoise,iM) = mean( abs(sum( (ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR))./ResultsALL(iNoise,iM).xA ,2)) );
        end
    end
end
% end

PostName = 'RelErr';
ytext = 'Avg. Rel. Error';
yaxlim = [0 1];
shaded_region = 0;
VIZ_ERROR_STATS

%% Training time
LOG_SCALE = 1;
data2plot = zeros(Nr,N_ETA,3);
% for iN = 1:N_LENGTHS
for iNoise = 1:N_ETA
    for iR = 1:Nr
        for iM = 1:3
            data2plot(iR,iNoise,iM) = ResultsALL(iNoise,iM).Ttraining(iR);
        end
    end
end
% end

PostName = 'TrainTime';
ytext = 'Training Time [s]';
yaxlim = [0 max(data2plot(:))];
shaded_region = 0;
LegLoc = 'NorthEast';
VIZ_ERROR_STATS
