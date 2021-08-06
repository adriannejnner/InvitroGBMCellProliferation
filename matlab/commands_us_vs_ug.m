load('dense_us_sims.mat')
dense_us_sims_infiltration = [infiltration_us_0;infiltration_us_000005;...
    infiltration_us_00001;infiltration_us_0001;infiltration_us_0005;infiltration_us_0009;...
    infiltration_us_001;infiltration_us_0011;infiltration_us_0015;infiltration_us_002;...
    infiltration_us_005;infiltration_us_01;infiltration_us_1];

dense_us_sims_infiltration_notstroma = [infiltration_us_notstroma_0;infiltration_us_notstroma_000005;...
    infiltration_us_notstroma_00001;infiltration_us_notstroma_0001;infiltration_us_notstroma_0005;...
    infiltration_us_notstroma_0009;infiltration_us_notstroma_001;infiltration_us_notstroma_0011;...
    infiltration_us_notstroma_0015;infiltration_us_notstroma_002;infiltration_us_notstroma_005;...
    infiltration_us_notstroma_01;infiltration_us_notstroma_1];

dense_us_sims_livecells = [livecells_us_0;livecells_us_000005;livecells_us_00001;livecells_us_0001;...
    livecells_us_0005;livecells_us_0009;livecells_us_001;livecells_us_0011;livecells_us_0015;...
    livecells_us_002;livecells_us_005;livecells_us_01;livecells_us_1];
    
    
dense_us_sims_deadcells = [deadcells_us_0;deadcells_us_000005;deadcells_us_00001;deadcells_us_0001;...
    deadcells_us_0005;deadcells_us_0009;deadcells_us_001;deadcells_us_0011;deadcells_us_0015;...
    deadcells_us_002;deadcells_us_005;deadcells_us_01;deadcells_us_1];

%%
load('sparse_us_sims.mat')
sparse_us_sims_infiltration = [infiltration_us_0;infiltration_us_000005;...
    infiltration_us_00001;infiltration_us_0001;infiltration_us_0005;infiltration_us_0009;...
    infiltration_us_001;infiltration_us_0011;infiltration_us_0015;infiltration_us_002;...
    infiltration_us_005;infiltration_us_01;infiltration_us_1];

sparse_us_sims_infiltration_notstroma = [infiltration_us_notstroma_0;infiltration_us_notstroma_000005;...
    infiltration_us_notstroma_00001;infiltration_us_notstroma_0001;infiltration_us_notstroma_0005;...
    infiltration_us_notstroma_0009;infiltration_us_notstroma_001;infiltration_us_notstroma_0011;...
    infiltration_us_notstroma_0015;infiltration_us_notstroma_002;infiltration_us_notstroma_005;...
    infiltration_us_notstroma_01;infiltration_us_notstroma_1];

sparse_us_sims_livecells = [livecells_us_0;livecells_us_000005;livecells_us_00001;livecells_us_0001;...
    livecells_us_0005;livecells_us_0009;livecells_us_001;livecells_us_0011;livecells_us_0015;...
    livecells_us_002;livecells_us_005;livecells_us_01;livecells_us_1];
    
    
sparse_us_sims_deadcells = [deadcells_us_0;deadcells_us_000005;deadcells_us_00001;deadcells_us_0001;...
    deadcells_us_0005;deadcells_us_0009;deadcells_us_001;deadcells_us_0011;deadcells_us_0015;...
    deadcells_us_002;deadcells_us_005;deadcells_us_01;deadcells_us_1];


%%
load('dense_ug_sims.mat')
dense_ug_sims_infiltration = [infiltration_ug_000001;infiltration_ug_000002;...
    infiltration_ug_00002;infiltration_ug_0001;infiltration_ug_00018;infiltration_ug_0004;...
    infiltration_ug_0008;infiltration_ug_002;infiltration_ug_02];

dense_ug_sims_infiltration_notstroma = [infiltration_ug_notstroma_000001;infiltration_ug_notstroma_000002;...
    infiltration_ug_notstroma_00002;infiltration_ug_notstroma_0001;infiltration_ug_notstroma_00018;...
    infiltration_ug_notstroma_0004;infiltration_ug_notstroma_0008;infiltration_ug_notstroma_002;infiltration_ug_notstroma_02];

dense_ug_sims_livecells = [livecells_ug_000001;livecells_ug_000002;livecells_ug_00002;livecells_ug_0001;...
    livecells_ug_00018;livecells_ug_0004;livecells_ug_0008;livecells_ug_002;livecells_ug_02];
    
    
dense_ug_sims_deadcells = [deadcells_ug_000001;deadcells_ug_000002;deadcells_ug_00002;deadcells_ug_0001;...
    deadcells_ug_00018;deadcells_ug_0004;deadcells_ug_0008;deadcells_ug_002;deadcells_ug_02];


%%
load('sparse_ug_sims.mat')
sparse_ug_sims_infiltration = [infiltration_ug_000001;infiltration_ug_000002;...
    infiltration_ug_00002;infiltration_ug_0001;infiltration_ug_00018;infiltration_ug_0004;...
    infiltration_ug_0008;infiltration_ug_002;infiltration_ug_02];

sparse_ug_sims_infiltration_notstroma = [infiltration_ug_notstroma_000001;infiltration_ug_notstroma_000002;...
    infiltration_ug_notstroma_00002;infiltration_ug_notstroma_0001;infiltration_ug_notstroma_00018;...
    infiltration_ug_notstroma_0004;infiltration_ug_notstroma_0008;infiltration_ug_notstroma_002;infiltration_ug_notstroma_02];

sparse_ug_sims_livecells = [livecells_ug_000001;livecells_ug_000002;livecells_ug_00002;livecells_ug_0001;...
    livecells_ug_00018;livecells_ug_0004;livecells_ug_0008;livecells_ug_002;livecells_ug_02];
    
    
sparse_ug_sims_deadcells = [deadcells_ug_000001;deadcells_ug_000002;deadcells_ug_00002;deadcells_ug_0001;...
    deadcells_ug_00018;deadcells_ug_0004;deadcells_ug_0008;deadcells_ug_002;deadcells_ug_02];


save('all_data.mat')


%%
load('all_data.mat')


%% plot end tumour fragment size for all simulations
xgrid_us = [0 0.00005 0.0001 0.001 0.005 0.009 0.01 0.011 0.015 0.02 0.05 0.1 1];

xgrid_ug = [0.00001 0.00002 0.0002 0.001 0.0018 0.004 0.008 0.02 0.2];

%sparse
figure
hold on 
yyaxis left
plot(xgrid_us,sparse_us_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (sparse)')
yyaxis right
plot(xgrid_us,dense_us_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (dense)')
xlabel('us')
set(gca,'FontSize',15)
set(gca,'xscale','log')

figure
hold on 
yyaxis left
plot(xgrid_ug,sparse_ug_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (sparse)')
yyaxis right
plot(xgrid_ug,dense_ug_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (dense)')
xlabel('ug')
set(gca,'FontSize',15)
set(gca,'xscale','log')


%% plot infiltration data for dense fragments ug vs us

figure
hold on 
bar(dense_ug_sims_infiltration_notstroma([3:3:27],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - ug - dense')
legend('0.00001','0.00002','0.0002','0.001','0.0018','0.004','0.008','0.02','0.2')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)


figure
hold on 
bar(dense_us_sims_infiltration_notstroma([3:3:39],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - us - dense')
legend('0','0.00005','0.0001','0.001','0.005','0.009','0.01','0.011','0.015','0.02','0.05','0.1','1')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)


figure
hold on 
bar(sparse_ug_sims_infiltration_notstroma([3:3:27],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - ug - sparse')
legend('0.00001','0.00002','0.0002','0.001','0.0018','0.004','0.008','0.02','0.2')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)


figure
hold on 
bar(sparse_us_sims_infiltration_notstroma([3:3:39],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - us - sparse')
legend('0','0.00005','0.0001','0.001','0.005','0.009','0.01','0.011','0.015','0.02','0.05','0.1','1')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)




%%
load('all_data.mat')
load('sparse_us_sims.mat')

figure
hold on 
bar([sparse_us_sims_infiltration_notstroma([36,39],:);infiltration_us_notstroma_125(3,:);infiltration_us_notstroma_15(3,:);infiltration_us_notstroma_2(3,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - us - sparse')
legend('0.1','1','1.25','1.5','2')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)
set(gca,'yscale','linear')

axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
bar([sparse_us_sims_infiltration_notstroma([36,39],2:5);infiltration_us_notstroma_125(3,2:5);infiltration_us_notstroma_15(3,2:5);infiltration_us_notstroma_2(3,2:5)]')
axis tight
ylim([0 5])
set(gca,'Xtick',linspace(1,4,4),'Xticklabel',{'50','100','150','200'})
set(gca,'FontSize',13)


load('dense_us_sims.mat')
figure
hold on 
bar([dense_us_sims_infiltration_notstroma([36,39],:);infiltration_us_notstroma_125(3,:);infiltration_us_notstroma_15(3,:);infiltration_us_notstroma_2(3,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - us - dense')
legend('0.1','1','1.25','1.5','2')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)


xgrid_us = [0.1 1 1.25 1.5 2];

figure
hold on 
yyaxis left

load('sparse_us_sims.mat')
plot(xgrid_us,[sparse_us_sims_livecells(12:13,end)',livecells_us_125(end),livecells_us_15(end),livecells_us_2(end)],'o:','LineWidth',2) 
ylabel('Number of live cells (sparse)')

yyaxis right
load('dense_us_sims.mat')
plot(xgrid_us,[sparse_us_sims_livecells(12:13,end)',livecells_us_125(end),livecells_us_15(end),livecells_us_2(end)],'o:','LineWidth',2) 
ylabel('Number of live cells (dense)')
xlabel('us')
set(gca,'FontSize',15)

