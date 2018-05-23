clear all, close all, clc
figpath = '../FIGURES/'; mkdir(figpath)
datapath = '../DATA/EX_LOTKA_Dependencies/'; mkdir(datapath)
addpath('../utils');

%% Paramaters
InputSignalType = 'sphs';
Ntrain_vec = 3000; %1000;
N_LENGTHS = length(Ntrain_vec);

eta_vec = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
N_ETA = length(eta_vec);

Nmodels = N_ETA*N_LENGTHS;

%% Load all models
ModelsNARX(1:Nmodels) = struct('name',[],'net',[],'stateDelays',[],'inputDelays',[],'hiddenSizes',[], ...
    'dt',[], 'Ttraining', [],'Err', [], 'ErrM', []);
ModelsDMDc(1:Nmodels) = struct('name',[],'sys',[],'Ndelay',[],'dt',[],'Err', [], 'ErrM', []);
ModelsSINDYc(1:Nmodels) = struct('name',[],'polyorder',[],'usesine',[],'Xi',[],'dt',[],'N',[], 'Err', [], 'ErrM', []);


ModelName = 'SINDYc';
datapath1 = [datapath,ModelName,'/'];
count = 0;
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        count = count + 1;
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'.mat']);
        load(filename);
        ModelsSINDYc(count) = Model;
    end
end

ModelName = 'NARX';
datapath1 = [datapath,ModelName,'/'];
count = 0;
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        count = count + 1;
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'.mat']);
        ModelsNARX(count) = Model;
    end
end


ModelName = 'DMDc';
datapath1 = [datapath,ModelName,'/'];
count = 0;
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        count = count + 1;
        filename = fullfile(datapath1,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_Eta',sprintf('%03g',100*eta_vec(iNoise)),'.mat']);
        ModelsDMDc(count) = Model;
    end
end
return
%%
clear ph
xx = 0:0.01:0.6;
tmp = zeros(length(N_ETA),1);
for i = 1:N_ETA, tmp(i) = ModelsDMDc(i).Err; end
yDMDc = spline(eta_vec(:),tmp,xx);
tmp = zeros(length(N_ETA),1);
for i = 1:N_ETA, tmp(i) = ModelsSINDYc(i).Err; end
ySINDYc = spline(eta_vec(:),tmp,xx(1:31));
% tmp = zeros(length(N_ETA),1);
% for i = 1:N_ETA, tmp(i) = ModelsNARX(i).Err; end
% yNARX = spline(eta_vec(:),tmp,xx);
figure, hold on, box on
plot(xx,yDMDc,'-r','LineWidth',1), hold on,
plot(xx(1:31),ySINDYc,'-g','LineWidth',1)
% plot(xx,yNARX,'-','Color',[0.7,0.7,1],'LineWidth',1)
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        ph(1) = plot(eta_vec(iNoise),ModelsDMDc(iNoise).Err,'or','MarkerFaceColor','r');
        ph(2) = plot(eta_vec(iNoise),ModelsSINDYc(iNoise).Err,'dg','MarkerFaceColor','g');
        ph(3) = plot(eta_vec(iNoise),ModelsNARX(iNoise).Err,'s','Color',[0.7,0.7,1],'MarkerFaceColor',[0.7,0.7,1]);
        %         plot(eta_vec(iNoise),ModelsNARX(iNoise).ErrM,'xr')
    end
end
xlabel('Eta')
ylabel('MSE')
l1 = legend(ph,'DMDc','SINDYc','NARX');
% set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.06 l1.Position(2)-0.07 l1.Position(3)-0.01 l1.Position(4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noise','.eps']);

%%
clear ph
xx = 0:0.01:0.6;
tmp = zeros(length(N_ETA),1);
for i = 1:N_ETA, tmp(i) = ModelsDMDc(i).ErrM; end
yDMDc = spline(eta_vec(:),tmp,xx);
tmp = zeros(length(N_ETA),1);
for i = 1:N_ETA, tmp(i) = ModelsSINDYc(i).ErrM; end
ySINDYc = spline(eta_vec(:),tmp,xx(1:31));
% tmp = zeros(length(N_ETA),1);
% for i = 1:N_ETA, tmp(i) = ModelsNARX(i).ErrM; end
% yNARX = spline(eta_vec(:),tmp,xx);
figure, hold on, box on
plot(xx,yDMDc,'-r','LineWidth',1), hold on,
plot(xx(1:31),ySINDYc,'-g','LineWidth',1)
% plot(xx,yNARX,'-','Color',[0.7,0.7,1],'LineWidth',1)
for iN = 1:N_LENGTHS
    for iNoise = 1:N_ETA
        ph(1) = plot(eta_vec(iNoise),ModelsDMDc(iNoise).ErrM,'or','MarkerFaceColor','r');
        ph(2) = plot(eta_vec(iNoise),ModelsSINDYc(iNoise).ErrM,'dg','MarkerFaceColor','g');
        ph(3) = plot(eta_vec(iNoise),ModelsNARX(iNoise).ErrM,'s','Color',[0.7,0.7,1],'MarkerFaceColor',[0.7,0.7,1]);
        %         plot(eta_vec(iNoise),ModelsNARX(iNoise).ErrM,'xr')
    end
end
xlabel('Eta')
ylabel('MSE')
l1 = legend(ph,'DMDc','SINDYc','NARX');
% set(l1,'Location','NorthEast')
l1.Position = [l1.Position(1)-0.06 l1.Position(2)-0.07 l1.Position(3)-0.01 l1.Position(4)];
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_PREDPERF_',InputSignalType,'_N',sprintf('%04g',Ntrain_vec(iN)),'_noiseM','.eps']);
