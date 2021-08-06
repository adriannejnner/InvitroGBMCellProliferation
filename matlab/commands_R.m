load('dense_R_sims.mat')
dense_R_sims_infiltration = [infiltration_R_05;infiltration_R_15;...
    infiltration_R_4;infiltration_R_10];

dense_R_sims_infiltration_notstroma = [infiltration_R_notstroma_05;infiltration_R_notstroma_15;...
    infiltration_R_notstroma_4;infiltration_R_notstroma_10];

dense_R_sims_livecells = [livecells_R_05;livecells_R_15;livecells_R_4;livecells_R_10];
    
    
dense_R_sims_deadcells = [deadcells_R_05;deadcells_R_15;deadcells_R_4;deadcells_R_10];

%%
load('sparse_R_sims.mat')

sparse_R_sims_infiltration = [infiltration_R_05;infiltration_R_15;...
    infiltration_R_4;infiltration_R_10];

sparse_R_sims_infiltration_notstroma = [infiltration_R_notstroma_05;infiltration_R_notstroma_15;...
    infiltration_R_notstroma_4;infiltration_R_notstroma_10];

sparse_R_sims_livecells = [livecells_R_05;livecells_R_15;livecells_R_4;livecells_R_10];
    
    
sparse_R_sims_deadcells = [deadcells_R_05;deadcells_R_15;deadcells_R_4;deadcells_R_10];

save('all_data_R.mat')


%%
load('all_data_R.mat')


%% plot end tumour fragment size for all simulations
xgrid_R = [0.5 1.5 4 10];

%sparse
figure
hold on 
yyaxis left
plot(xgrid_R,sparse_R_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (sparse)')
yyaxis right
plot(xgrid_R,dense_R_sims_livecells(:,end),'o:','LineWidth',2) 
ylabel('Number of live cells (dense)')
xlabel('R')
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'Xtick',[0.5 1.5 4 10])


%% plot infiltration data for dense fragments ug vs us

figure
hold on 
bar(dense_R_sims_infiltration_notstroma([3:3:12],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - R - dense')
legend('0.5', '1.5','4','10')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)


figure
hold on 
bar(sparse_R_sims_infiltration_notstroma([3:3:12],:)')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours - R - sparse')
legend('0.5', '1.5','4','10')
xlabel('\mu m from periphery')
xlim([1 14])
set(gca,'FontSize',15)





