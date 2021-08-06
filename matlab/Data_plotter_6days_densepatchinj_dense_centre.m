Dense_Centre = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Centre injection\Dense_ 6 days\Dense_Centre_6day.mat');
Dense_Periphery = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\us_1_6 days\Dense_periphery_us_1.mat');
Patchy_dense_injections = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\patchy\patch_injections_dense_6days\Patchy_dense_injections_6days.mat');

colmap = [102,194,165
252,141,98
141,160,203
231,138,195
166,216,84]/255;

figure
hold on 
plot(Dense_Centre.live_cells/Dense_Centre.live_cells(1)*100,'Color',colmap(1,:),'LineWidth',2)
plot(Dense_Periphery.live_cells/Dense_Periphery.live_cells(1)*100,'Color',colmap(2,:),'LineWidth',2)
plot(Patchy_dense_injections.live_cells/Patchy_dense_injections.live_cells(1)*100,'Color',colmap(3,:),'LineWidth',2)
legend('1. Dense centre','2. Dense periphery','3. Centre of patches')
xlabel('Time (days)')
ylabel('% GBM cells (uninfected + infected)')
set(gca,'FontSize',19)
set(gca,'XTick',[0 24 48 72],'XTickLabels',{'0','2','4','6'})  
xlim([0 72])


figure
hold on 
plot(Dense_Periphery.live_cells,'LineWidth',2)
plot(Dense_Periphery.infected_cells,'LineWidth',2)
plot(Dense_Periphery.dead_cells,'LineWidth',2)
legend('Uninfected GBM cells','Infected GBM cells','Dead cells')
set(gca,'FontSize',19)
set(gca,'XTick',[0 24 48 72],'XTickLabels',{'0','2','4','6'})  
xlim([0 72])
ylabel('Cell count')
xlabel('Time (days)')

figure
hold on 
yyaxis left
plot(Dense_Periphery.virion,'LineWidth',2)
ylabel('Virions')
yyaxis right
plot(Dense_Periphery.chemokine,'LineWidth',2)
ylabel('Chemokine')
legend('Total virions','Total chemokine')
set(gca,'FontSize',19)
set(gca,'XTick',[0 24 48 72],'XTickLabels',{'0','2','4','6'})  
xlim([0 72])
xlabel('Time (days)')
ax = gca;
ax.YAxis(1).Color = [0.49 0.18 0.56];
ax.YAxis(2).Color = [0.47 0.67 0.19];

