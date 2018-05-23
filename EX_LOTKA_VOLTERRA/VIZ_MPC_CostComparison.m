% LOTKA-VOLTERRA system
% Visualize cost function for different cases

clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

%% Load Data
% Nvec = [1,3,5,7,10,15,20,25,30,35,40,45,50];
Nvec = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,20,25,30,35,40,45,50];
Nn = length(Nvec);
Data(1:Nn) = struct('t',[],'x',[],'u',[],'J',[]);

for i = 1:Nn
    load(fullfile(datapath,['EX_LOTKA_MPC_SINDYc_N',num2str(Nvec(i)),'.mat']))
    Data(i) = Results;
end

Ton = 30;
tspan = [100 200];
clear Results

%% Correct cost
a = .5;
b = .025;
d = .005;
g = .5;
Q = [1,0]; R = 0.5; Ru = 0.5;
xref = repmat([g/d;a/b],[1 size(Data(iM).x,2)]) ;
for iM = 1:Nn
    Data(iM).J = evalObjectiveFCN(Data(iM).u,Data(iM).x,xref,diag(Q),R,Ru); %Data(iM).xref
    idx = find(Data(iM).u ~= 0,1,'first');
    Data(iM).J(1:idx) = 0;
end


%% Generate Data
ymax = max([max(cumsum(Data(1).J)) max(cumsum(Data(2).J))])+500;
ccolors = jet(2*Nn); ccolors = ccolors(1:2:end,:);

figure;hold on, box on,
plot([Ton+tspan(1),Ton+tspan(1)],[0.01 ymax],':','Color',[0.7,0.7,0.7],'LineWidth',1)
t1 = text(51+tspan(1),4.5*10^4,'Control', 'FontSize',12);
t2 = text(51+tspan(1),4*10^4,'turned on', 'FontSize',12);
for i = 1:Nn
    ph(i) = plot(Data(i).t+tspan(1),cumsum(Data(i).J),'-','Color',ccolors(i,:),'LineWidth',1.5);
end
xlabel('Time')
ylabel('Cost')
legend(ph(1:2:end),num2str(Nvec(1:2:end)'),'Location','eastoutside','Orientation','vertical')
axis tight
xlim([100 200]), ylim([0 5*10^4])
set(gca,'yscale','linear','xscale','linear','xtick',[100,150,200])
set(gca,'xtick',[50,100,150,200])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 400 245])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', [figpath,'EX_LOTKA_MPC-SINDYc_CostComparison_All_v2.eps']);

delete(t1); delete(t2);
xlim([100 200]), ylim([10^3 10^6])
t1 = text(31+tspan(1),6*10^5,'Control', 'FontSize',12);
t2 = text(31+tspan(1),4*10^5,'turned on', 'FontSize',12);
set(gca,'yscale','log','xscale','linear','xtick',[100,150,200])
print('-depsc2', [figpath,'EX_LOTKA_MPC-SINDYc_CostComparison_All_v2_log.eps']);
