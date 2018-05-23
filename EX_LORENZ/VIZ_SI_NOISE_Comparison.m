% LORENZ system
% Visualization for dependency on noise level

clear all, close all, clc

SystemModel = 'LORENZ';

figpath = '../FIGURES/'; mkdir(figpath)
datapath = ['../DATA/EX_',SystemModel,'_Dependencies/'];
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
Nt = 1000; % length(tv);
Nvar = 3;
%% Load all models
ResultsALL(1:Nmodels,3) = struct('err', zeros(Nr,1), 'errM', zeros(Nr,1), 'xA', zeros(Nt,Nvar), 'xB', zeros(Nt,Nvar,Nr),'Ttraining',zeros(Nr,1));

count = 0;
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        count = count + 1;
        
        
        %         ModelName = 'DMDc';
        %         datapath1 = [datapath,'DMDc/'];
        %         filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
        %         load(filename);
        %         ResultsALL(count,1) = Results;
        ResultsALL(count,1) = struct('err', zeros(Nr,1), 'errM', zeros(Nr,1), 'xA', zeros(Nt,Nvar), 'xB', zeros(Nt,Nvar,Nr),'Ttraining',zeros(Nr,1));
        
        ModelName = 'SINDYc';
        datapath1 = [datapath,'SINDYc/'];
        filename = fullfile(datapath1,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
        load(filename);
        ResultsALL(count,2) = Results;
        
        
        ModelName = 'NARX';
        datapath1 = [datapath,'NARX/',NARXtraining,'/'];
        filename = fullfile(datapath1,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'_STATS.mat']);
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
    for i = 1:N_ETA, err(1:Nr,i,iM) = ResultsALL(i,iM).err(1:Nr); end
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
data = err(1:Nr,eta_vec==0.4,2);% sine 0.4, sphs 0.45: 278.4062
% ylimregion = median(data(isnan(data)==0))
ylimregion = 300;
fillh1 = fill([0.1 11.9 11.9 0.1], [10 10 ylimregion.*ones(1,2)],0.9*ones(1,3));
fillh1.EdgeColor = 0.9*ones(1,3); fillh1.FaceAlpha = 0.5;
%ph(1)=plot(xx,yy(1,:),'-r','LineWidth',1); hold on,
ph(2)=plot(xx(1:end),yy(2,1:end),'-g','LineWidth',1);
ph(3)=plot(xx(1:end),yy(3,1:end),'-','Color',[0.7,0.7,1],'LineWidth',1);
ph(1) = [];
% plot([0,12],max(err(1:Nr,eta_vec==0.3,2)).*ones(1,2),'-k')
for iM = 2:3
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
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'_noleg.eps']);

% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthWest')
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise_',NARXtraining,'.eps']);


%%
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
symbols = {'o', 'd', 's'};
xx = 1:0.1:11;%0:0.01:0.6;
yy = zeros(3,length(xx));
err = zeros(Nr,N_ETA,3);
for iM = 1:3
    for i = 1:N_ETA, err(1:Nr,i,iM) = ResultsALL(i,iM).errM(1:Nr); end
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
% ph(1)=plot(xx,yy(1,:),'-r','LineWidth',1); hold on,
ph(2)=plot(xx(1:end),yy(2,1:end),'-g','LineWidth',1);
ph(3)=plot(xx(1:end),yy(3,1:end),'-','Color',[0.7,0.7,1],'LineWidth',1);
ph(1) = [];
% plot([0,12],max(err(1:Nr,eta_vec==0.3,2)).*ones(1,2),'-k')
for iM = 2:3
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
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noiseM_',NARXtraining,'_noleg.eps']);

% l1 = legend(ph,'DMDc','SINDYc','NARX');
l1 = legend(ph,'SINDYc','NARX');
set(l1,'Location','NorthWest')
% l1.Position = [l1.Position(1)-0.06 l1.Position(2)-0.07 l1.Position(3)-0.01 l1.Position(4)];
% l1.Position = [l1.Position(1)+0.08 l1.Position(2)-0.08 l1.Position(3)-0.01 l1.Position(4)];
l1.Position = [l1.Position(1)-0.01 l1.Position(2)+0.02 l1.Position(3)-0.01 l1.Position(4)];
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noiseM_',NARXtraining,'.eps']);

%% Normalized error
clear LegLoc
LOG_SCALE = 0;
data2plot = zeros(Nr,N_ETA,3);
% for iN = 1:N_LENGTHS
for iNoise = 1:N_ETA
    for iR = 1:Nr
        for iM = 1:3
            data2plot(iR,iNoise,iM) = mean( sum( abs(ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR))./abs(ResultsALL(iNoise,iM).xA) ,2)  );
            %data2plot(iR,iNoise,iM) = mean( abs(sum( (ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR))./ResultsALL(iNoise,iM).xA ,2)) );
        end
    end
end
% end

PostName = 'RelErr';
ytext = 'Avg. Rel. Error';
yaxlim = [0 14];
shaded_region = 0;
VIZ_ERROR_STATS


%%
iNoise = 1; iM = 3;
figure, hold on
plot(ResultsALL(iNoise,iM).xA,'-k')
for i = 1:Nr
    plot(ResultsALL(iNoise,iM).xB(:,:,iR),'--r');%,'Color',[0.7,0.7,0.7])
end
%%
ModelSelection = {'DMDc','SINDYc', 'NARX'};
dt = 0.01;
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
tspan = [1:Nt].*dt+30;
iM = 3;
data2plot = zeros(Nt,N_ETA,3);

for iM = 2:3
    PostName = ['_TS_',ModelSelection{iM}];
    for iNoise = 1:N_ETA
        f1 = figure('visible','off');hold on, box on
        xshift = [-35,0,+20];
        for iVar = 1:Nvar
            %data2plot(iT,iNoise,iM) = (ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR));
            %boxplot(squeeze(ResultsALL(iNoise,iM).xB(1:10:end,1,:))')
            Xstats = squeeze(ResultsALL(iNoise,iM).xB(1:1:end,iVar,:))';
            Ystats = prctile(Xstats,[25 50 75],1);
            
            X=[tspan(1:end),fliplr(tspan(1:end))];                %#create continuous x value array for plotting
            Y=[Ystats(1,:)+xshift(iVar),fliplr(Ystats(3,:))+xshift(iVar)];              %#create y values for out and then back
            fillh1 = fill(X,Y,ccolors(iM,:));                  %#plot filled area
            fillh1.EdgeColor = ccolors(iM,:);
            fillh1.FaceAlpha = 0.3;
            fillh1.EdgeAlpha = 0.3;
            
            plot(tspan,Ystats(2,:)+xshift(iVar),'-','Color',ccolors(iM,:),'LineWidth',1)
            %                 plot(tspan,Ystats(1,:),'--','Color',ccolors(iM,:))
            %                 plot(tspan,Ystats(3,:),'--','Color',ccolors(iM,:))
            plot(tspan,ResultsALL(iNoise,iM).xA(1:1:end,iVar)+xshift(iVar),'-k','LineWidth',1)%,'Color',[0.7,0.7,0.7])
        end
        ylim([-60 65])
        ylabel('xi'), xlabel('Time')
        set(gca,'LineWidth',1, 'FontSize',14)
        set(gcf,'Position',[100 100 300 200])
        set(gcf,'PaperPositionMode','auto'),
        drawnow
        print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise',sprintf('%03g',100*eta_vec(iNoise)),'_',NARXtraining,'_',PostName,'_noleg.eps']);
        close(f1);
        
    end
end

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

clear YSCALE_SET
PostName = 'TrainTime';
ytext = 'Training Time [s]';
yaxlim = [0 max(data2plot(:))];
shaded_region = 0;
LegLoc = 'NorthEast';
VIZ_ERROR_STATS

% Zoom
% PostName = 'TrainTimeZOOM';
% yaxlim = [0 0.04];
% VIZ_ERROR_STATS
%% Prediction horizon
dt = 0.01;
clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
ERR_RADIUS = 0.1;
PredHor_Ball = zeros(Nr,N_ETA,Nmodels);
PredHor_Max = zeros(Nr,N_ETA,Nmodels);
%Err_SMAPE = zeros(Nr,N_ETA,Nmodels); 
Err_Ratio = zeros(Nr,N_ETA,Nmodels); 
for iM = 2:3
    for iNoise = 1:N_ETA
        for iR = 1:Nr
            tmp =  (ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR))./ResultsALL(iNoise,iM).xA ;
            TF = abs(tmp)>0.1;
            idx = [];
            for iVar = 1:Nvar
                idx = [idx ; find(TF(:,iVar)==1,1,'first')];
            end
            PredHor_Max(iR,iNoise,iM) = max(idx);
            
%             tmp2 = sqrt(sum(tmp.^2,2));
            tmp2 = sqrt(sum(abs((ResultsALL(iNoise,iM).xA - ResultsALL(iNoise,iM).xB(:,:,iR))).^2,2));
            idx = find(tmp2>3,1,'first');
            PredHor_Ball(iR,iNoise,iM) = idx-1;
            
            %Err_SMAPE(iR,iNoise,iM) =  sum(sum(abs(sum(ResultsALL(iNoise,iM).xB(:,:,iR)-ResultsALL(iNoise,iM).xA)./sum(ResultsALL(iNoise,iM).xA + ResultsALL(iNoise,iM).xB(:,:,iR)))./Nt))./Nvar;
            Err_Ratio(iR,iNoise,iM) =  sum(sum( abs((ResultsALL(iNoise,iM).xB(:,:,iR)./ResultsALL(iNoise,iM).xA))  )./Nt)./Nvar;
            
        end
    end
end

%%
data2plot = (PredHor_Ball).*dt;
PostName = 'PredHoriz_Ball';
ytext = 'Pred. Horizon';
yaxlim = [10^-3 10^1];
%yaxlim = [0 10];
LOG_SCALE = 1;
YSCALE_SET = 1;
shaded_region = 0;
LegLoc = 'NorthEast';
VIZ_ERROR_STATS

return
data2plot = PredHor_Max;
PostName = 'PredHoriz_Max';
ytext = 'Pred. Horizon';
yaxlim = [10^-1 10^3];
LOG_SCALE = 1;
YSCALE_SET = 1;
shaded_region = 0;
LegLoc = 'NorthEast';
VIZ_ERROR_STATS

data2plot = Err_Ratio;
PostName = 'Err_Ratio';
ytext = 'Error';
yaxlim = [10^-1 10^1];
LOG_SCALE = 1;
YSCALE_SET = 1;
shaded_region = 0;
LegLoc = 'NorthEast';
VIZ_ERROR_STATS



