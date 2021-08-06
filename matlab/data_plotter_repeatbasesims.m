Dense1 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun1\Dense_OV_sim1.mat');
Dense2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun2\Dense_OV_sim2.mat');
Dense3 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun3\Dense_OV_sim3.mat');
Dense4 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun4\Dense_OV_sim4.mat');
Dense5 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun5\Dense_OV_sim5.mat');
Dense6 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun6\Dense_OV_sim6.mat');
Dense7 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun7\Dense_OV_sim7.mat');
Dense8 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun8\Dense_OV_sim8.mat');
Dense9 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Dense_runs\Base sims\denserun9\Dense_OV_sim9.mat');

Sparse1 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun1\Sparse_OV_sim1.mat');
Sparse2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun2\Sparse_OV_sim2.mat');
Sparse3 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun3\Sparse_OV_sim3.mat');
Sparse4 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun4\Sparse_OV_sim4.mat');
Sparse5 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun5\Sparse_OV_sim5.mat');
Sparse6 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun6\Sparse_OV_sim6.mat');
Sparse7 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun7\Sparse_OV_sim7.mat');
Sparse8 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun8\Sparse_OV_sim8.mat');
Sparse9 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Sparse_runs\sparserun9\Sparse_OV_sim9.mat');

Patchy1 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun1\Patchy_sim1.mat');
Patchy2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun2\Patchy_sim2.mat');
Patchy3 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun3\Patchy_sim3.mat');
Patchy4 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun4\Patchy_sim4.mat');
Patchy5 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun5\Patchy_sim5.mat');
Patchy6 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\Patchy_runs\patchyrun6\Patchy_sim6.mat');

%% Dense

cmap = [165,0,38
215,48,39
244,109,67
253,174,97
254,224,144
224,243,248
171,217,233
116,173,209
69,117,180
49,54,149]/255;

figure
hold on 
plot(linspace(0,72,72)/12,Dense1.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.dead,'Color',cmap(3,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Dead cells')
set(gca,'FontSize',18)
saveas(gcf,'Dense_DeadCells_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);

Dense_Virus_Mat = [Dense1.dead;Dense2.dead;Dense3.dead;Dense4.dead;Dense5.dead;Dense6.dead;Dense7.dead;Dense8.dead;Dense9.dead];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Dead GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Dense_DeadCells_shade.png')

figure
hold on
plot(linspace(0,72,72)/12,Dense1.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.infected,'Color',cmap(2,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
set(gca,'FontSize',18)
saveas(gcf,'Dense_InfectedCells_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(2,:);
options.color_line = cmap(2,:);

Dense_Virus_Mat = [Dense1.infected;Dense2.infected;Dense3.infected;Dense4.infected;Dense5.infected;Dense6.infected;Dense7.infected;Dense8.infected;Dense9.infected];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Dense_InfectedCells_shade.png')


figure
hold on 
plot(linspace(0,72,72)/12,Dense1.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
saveas(gcf,'Dense_UninfectedCells_indiv.png')

fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(1,:);
options.color_line = cmap(1,:);

Dense_Virus_Mat = [Dense1.uninfected_live;Dense2.uninfected_live;Dense3.uninfected_live;Dense4.uninfected_live;Dense5.uninfected_live;Dense6.uninfected_live;Dense7.uninfected_live;Dense8.uninfected_live;Dense9.uninfected_live];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Dense_UninfectedCells_shade.png')

figure
hold on 
plot(linspace(0,72,72)/12,Dense1.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
saveas(gcf,'Dense_Virus_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Dense_Virus_Mat = [Dense1.extracellular_virus;Dense2.extracellular_virus;Dense3.extracellular_virus;Dense4.extracellular_virus;Dense5.extracellular_virus;Dense6.extracellular_virus;Dense7.extracellular_virus;Dense8.extracellular_virus;Dense9.extracellular_virus];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
box off
saveas(gcf,'Dense_Virus_shade.png')


% all together
fig1 = figure;
yyaxis left
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Dense_Infected_Mat = [Dense1.infected;Dense2.infected;Dense3.infected;Dense4.infected;Dense5.infected;Dense6.infected;Dense7.infected;Dense8.infected;Dense9.infected];
Dense_Dead_Mat = [Dense1.dead;Dense2.dead;Dense3.dead;Dense4.dead;Dense5.dead;Dense6.dead;Dense7.dead;Dense8.dead;Dense9.dead];
Dense_Uninfected_Mat = [Dense1.uninfected_live;Dense2.uninfected_live;Dense3.uninfected_live;Dense4.uninfected_live;Dense5.uninfected_live;Dense6.uninfected_live;Dense7.uninfected_live;Dense8.uninfected_live;Dense9.uninfected_live];
Dense_Virus_Mat = [Dense1.extracellular_virus;Dense2.extracellular_virus;Dense3.extracellular_virus;Dense4.extracellular_virus;Dense5.extracellular_virus;Dense6.extracellular_virus;Dense7.extracellular_virus;Dense8.extracellular_virus;Dense9.extracellular_virus];

options.color_area = cmap(1,:);
options.color_line = cmap(1,:);
plot_areaerrorbar(Dense_Uninfected_Mat,options)
hold on
%plot_areaerrorbar(Dense_Infected_Mat,options)
hold on
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);
plot_areaerrorbar(Dense_Dead_Mat,options)
hold on 
ylabel('Cells')
yyaxis right
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Dense')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Dense_shade_alltogether.png')


figure
hold on 
yyaxis left

plot(linspace(0,72,72)/12,Dense1.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)

plot(linspace(0,72,72)/12,Dense1.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.dead,':','Color',cmap(3,:),'LineWidth',1.5)

ylabel('Cells')
yyaxis right
plot(linspace(0,72,72)/12,Dense1.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense2.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense3.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense4.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense5.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense6.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense7.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense8.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Dense9.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Dense')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Dense_indiv_alltogether.png')


%% Sparse

cmap = [165,0,38
215,48,39
244,109,67
253,174,97
254,224,144
224,243,248
171,217,233
116,173,209
69,117,180
49,54,149]/255;

figure
hold on 
plot(linspace(0,72,72)/12,Sparse1.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.dead,'Color',cmap(3,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Dead cells')
set(gca,'FontSize',18)
saveas(gcf,'Sparse_DeadCells_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);

Dense_Virus_Mat = [Sparse1.dead;Sparse2.dead;Sparse3.dead;Sparse4.dead;Sparse5.dead;Sparse6.dead;Sparse7.dead;Sparse8.dead;Sparse9.dead];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Dead cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Sparse_DeadCells_shade.png')

figure
hold on
plot(linspace(0,72,72)/12,Sparse1.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.infected,'Color',cmap(2,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
saveas(gcf,'Sparse_InfectedCells_indiv.png')

fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);

Dense_Virus_Mat = [Sparse1.infected;Sparse2.infected;Sparse3.infected;Sparse4.infected;Sparse5.infected;Sparse6.infected;Sparse7.infected;Sparse8.infected;Sparse9.infected];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Sparse_InfectedCells_shade.png')


figure
hold on 
plot(linspace(0,72,72)/12,Sparse1.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
saveas(gcf,'Sparse_UninfectedCells_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(1,:);
options.color_line = cmap(1,:);

Dense_Virus_Mat = [Sparse1.uninfected_live;Sparse2.uninfected_live;Sparse3.uninfected_live;Sparse4.uninfected_live;Sparse5.uninfected_live;Sparse6.uninfected_live;Sparse7.uninfected_live;Sparse8.uninfected_live;Sparse9.uninfected_live];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Sparse_UninfectedCells_shade.png')

figure
hold on 
plot(linspace(0,72,72)/12,Sparse1.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
saveas(gcf,'Sparse_Virus_indiv.png')



fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Dense_Virus_Mat = [Sparse1.extracellular_virus;Sparse2.extracellular_virus;Sparse3.extracellular_virus;Sparse4.extracellular_virus;Sparse5.extracellular_virus;Sparse6.extracellular_virus;Sparse7.extracellular_virus;Sparse8.extracellular_virus;Sparse9.extracellular_virus];
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
box off
saveas(gcf,'Sparse_Virus_shade.png')



% all together
fig1 = figure;
yyaxis left
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Dense_Infected_Mat = [Sparse1.infected;Sparse2.infected;Sparse3.infected;Sparse4.infected;Sparse5.infected;Sparse6.infected;Sparse7.infected;Sparse8.infected;Sparse9.infected];
Dense_Dead_Mat = [Sparse1.dead;Sparse2.dead;Sparse3.dead;Sparse4.dead;Sparse5.dead;Sparse6.dead;Sparse7.dead;Sparse8.dead;Sparse9.dead];
Dense_Uninfected_Mat = [Sparse1.uninfected_live;Sparse2.uninfected_live;Sparse3.uninfected_live;Sparse4.uninfected_live;Sparse5.uninfected_live;Sparse6.uninfected_live;Sparse7.uninfected_live;Sparse8.uninfected_live;Sparse9.uninfected_live];
Dense_Virus_Mat = [Sparse1.extracellular_virus;Sparse2.extracellular_virus;Sparse3.extracellular_virus;Sparse4.extracellular_virus;Sparse5.extracellular_virus;Sparse6.extracellular_virus;Sparse7.extracellular_virus;Sparse8.extracellular_virus;Sparse9.extracellular_virus];

options.color_area = cmap(1,:);
options.color_line = cmap(1,:);
plot_areaerrorbar(Dense_Uninfected_Mat,options)
hold on
%plot_areaerrorbar(Dense_Infected_Mat,options)
hold on
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);
plot_areaerrorbar(Dense_Dead_Mat,options)
hold on 
ylabel('Cells')
yyaxis right
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Sparse')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Sparse_shade_alltogether.png')


figure
hold on 
yyaxis left

plot(linspace(0,72,72)/12,Sparse1.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)

plot(linspace(0,72,72)/12,Sparse1.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.dead,':','Color',cmap(3,:),'LineWidth',1.5)

ylabel('Cells')
yyaxis right
plot(linspace(0,72,72)/12,Sparse1.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse2.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse3.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse4.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse5.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse6.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse7.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse8.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Sparse9.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Sparse')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Sparse_indiv_alltogether.png')

%% Patchy

cmap = [165,0,38
215,48,39
244,109,67
253,174,97
254,224,144
224,243,248
171,217,233
116,173,209
69,117,180
49,54,149]/255;

figure
hold on 
plot(linspace(0,72,72)/12,Patchy1.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.dead,'Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.dead,'Color',cmap(3,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse7.dead,'Color',cmap(3,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse8.dead,'Color',cmap(3,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse9.dead,'Color',cmap(3,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Dead cells')
set(gca,'FontSize',18)
saveas(gcf,'Patchy_DeadCells_indiv.png')

fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);

Patchy_Virus_Mat = [Patchy1.dead;Patchy2.dead;Patchy3.dead;Patchy4.dead;Patchy5.dead;Patchy6.dead];
plot_areaerrorbar(Patchy_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Dead GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Patchy_DeadCells_shade.png')

figure
hold on
plot(linspace(0,72,72)/12,Patchy1.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.infected,'Color',cmap(2,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.infected,'Color',cmap(2,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse7.infected,'Color',cmap(2,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse8.infected,'Color',cmap(2,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Sparse9.infected,'Color',cmap(2,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
saveas(gcf,'Patchy_InfectedCells_indiv.png')

fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(2,:);
options.color_line = cmap(2,:);

Patchy_Virus_Mat = [Patchy1.infected;Patchy2.infected;Patchy3.infected;Patchy4.infected;Patchy5.infected;Patchy6.infected];
plot_areaerrorbar(Patchy_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Infected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Patchy_InfectedCells_shade.png')

figure
hold on 
plot(linspace(0,72,72)/12,Patchy1.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy7.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy8.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy9.uninfected_live,'Color',cmap(1,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
saveas(gcf,'Patchy_UninfectedCells_indiv.png')


fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(1,:);
options.color_line = cmap(1,:);

Patchy_Virus_Mat = [Patchy1.uninfected_live;Patchy2.uninfected_live;Patchy3.uninfected_live;Patchy4.uninfected_live;Patchy5.uninfected_live;Patchy6.uninfected_live];
plot_areaerrorbar(Patchy_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Cells')
title('Uninfected GBM cells')
set(gca,'FontSize',18)
box off
saveas(gcf,'Patchy_UninfectedCells_shade.png')


figure
hold on 
plot(linspace(0,72,72)/12,Patchy1.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy7.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy8.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
%plot(linspace(0,72,72)/12,Patchy9.extracellular_virus,'Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
saveas(gcf,'Patchy_Virus_indiv.png')

fig1 = figure;
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Patchy_Virus_Mat = [Patchy1.extracellular_virus;Patchy2.extracellular_virus;Patchy3.extracellular_virus;Patchy4.extracellular_virus;Patchy5.extracellular_virus;Patchy6.extracellular_virus];
plot_areaerrorbar(Patchy_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Oncolytic virus')
set(gca,'FontSize',18)
box off
saveas(gcf,'Patchy_Virus_shade.png')

% all together
fig1 = figure;
yyaxis left
options.handle = fig1;
options.alpha = 0.5;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,72,72)/12;
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);

Dense_Infected_Mat = [Patchy1.infected;Patchy2.infected;Patchy3.infected;Patchy4.infected;Patchy5.infected;Patchy6.infected];
Dense_Dead_Mat = [Patchy1.dead;Patchy2.dead;Patchy3.dead;Patchy4.dead;Patchy5.dead;Patchy6.dead];
Dense_Uninfected_Mat = [Patchy1.uninfected_live;Patchy2.uninfected_live;Patchy3.uninfected_live;Patchy4.uninfected_live;Patchy5.uninfected_live;Patchy6.uninfected_live];
Dense_Virus_Mat = [Patchy1.extracellular_virus;Patchy2.extracellular_virus;Patchy3.extracellular_virus;Patchy4.extracellular_virus;Patchy5.extracellular_virus;Patchy6.extracellular_virus];

options.color_area = cmap(1,:);
options.color_line = cmap(1,:);
plot_areaerrorbar(Dense_Uninfected_Mat,options)
hold on
%plot_areaerrorbar(Dense_Infected_Mat,options)
hold on
options.color_area = cmap(3,:);
options.color_line = cmap(3,:);
plot_areaerrorbar(Dense_Dead_Mat,options)
hold on 
ylabel('Cells')
yyaxis right
options.color_area = cmap(4,:);
options.color_line = cmap(4,:);
plot_areaerrorbar(Dense_Virus_Mat,options)
xlabel('Time (days)')
ylabel('Virions')
title('Patchy')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Patchy_shade_alltogether.png')

figure
hold on 
yyaxis left

plot(linspace(0,72,72)/12,Patchy1.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.uninfected_live,':','Color',cmap(1,:),'LineWidth',1.5)

plot(linspace(0,72,72)/12,Patchy1.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.dead,':','Color',cmap(3,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.dead,':','Color',cmap(3,:),'LineWidth',1.5)

ylabel('Cells')
yyaxis right
plot(linspace(0,72,72)/12,Patchy1.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy2.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy3.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy4.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy5.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
plot(linspace(0,72,72)/12,Patchy6.extracellular_virus,':','Color',cmap(4,:),'LineWidth',1.5)
xlabel('Time (days)')
ylabel('Virions')
title('Patchy')
set(gca,'FontSize',18)
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gcf,'Patchy_indiv_alltogether.png')
