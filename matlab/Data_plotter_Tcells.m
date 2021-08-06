control_case = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\Immune control with tumour specific antigen\contol_immune_tumour_antigen.mat');
base_case = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\base case\dense_homogeneous_us_0point01.mat');
dense_immune_double = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\double immune\dense_immune_double.mat');
dense_immune_double_killers = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\double immune killers\dense_immune_double_killers.mat');
dense_immune_hundredfold_killers = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\hundred fold immune killers\dense_immune_hundredfold_killers.mat');%_v2
dense_immune_hundredfold = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\hundredfold immune\dense_immune_hundredfold.mat');
dense_immune_tenfold = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\tenfold immune\dense_immune_tenfold.mat');
dense_immune_tenfold_killers = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\tenfold immune killers\dense_immune_tenfold_killers.mat');
base_case_killers = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\base case killers\base_case_killers.mat');
%base_case_killers = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM_OV project\PhysiCell-GBM-OV project\matlab\dense\immune\base case killers 2\base_case_killers_2.mat');

colmap = [255,255,204
199,233,180
127,205,187
65,182,196
29,145,192
34,94,168
12,44,132]/255;

dense_initial_cell_number = base_case.live_cells(1);

figure
hold on 
plot(base_case.live_cells/dense_initial_cell_number*100,'Color',colmap(1,:),'LineWidth',1)
plot(dense_immune_double.live_cells/dense_initial_cell_number*100,'Color',colmap(2,:),'LineWidth',1)
plot(dense_immune_tenfold.live_cells/dense_initial_cell_number*100,'Color',colmap(3,:),'LineWidth',1)
plot(dense_immune_hundredfold.live_cells/dense_initial_cell_number*100,'Color',colmap(4,:),'LineWidth',1)
plot(base_case_killers.live_cells/dense_initial_cell_number*100,'Color',colmap(1,:),'LineWidth',1)
plot(dense_immune_double_killers.live_cells/dense_initial_cell_number*100,'Color',colmap(5,:),'LineWidth',1)
plot(dense_immune_tenfold_killers.live_cells/dense_initial_cell_number*100,'Color',colmap(6,:),'LineWidth',1)
plot(dense_immune_hundredfold_killers.live_cells/dense_initial_cell_number*100,'Color',colmap(7,:),'LineWidth',1)
legend('u_s = 0','u_s = 0.00005','u_s = 0.0001','u_s = 0.001','u_s = 0.01','u_s = 0.1','u_s = 1')

figure
hold on 
plot(base_case.infected_cells,'Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double.infected_cells,'Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold.infected_cells,'Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold.infected_cells,'Color',[1,0.64,0.45],'LineWidth',2)
plot(base_case_killers.infected_cells,'--','Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double_killers.infected_cells,'--','Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold_killers.infected_cells,'--','Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold_killers.infected_cells,'--','Color',[1,0.64,0.45],'LineWidth',2)
xlabel('Time (hours)')
ylabel('No. of infected cells')
set(gca,'FontSize',19)
set(gca,'XTick',[0 12 24 36],'XTickLabels',{'0','24','48','72'})  
xlim([0 36])

figure
hold on 
b = bar([base_case.live_cells(end)/dense_initial_cell_number*100,... 
	dense_immune_double.live_cells(end)/dense_initial_cell_number*100,... 
	dense_immune_tenfold.live_cells(end)/dense_initial_cell_number*100,...
	dense_immune_hundredfold.live_cells(end)/dense_initial_cell_number*100;...
    base_case_killers.live_cells(end)/dense_initial_cell_number*100,...
	dense_immune_double_killers.live_cells(end)/dense_initial_cell_number*100,...
	dense_immune_tenfold_killers.live_cells(end)/dense_initial_cell_number*100,...
	dense_immune_hundredfold_killers.live_cells(end)/dense_initial_cell_number*100],'FaceColor','flat')
ylim([20 45])
set(gca,'FontSize',19)
ylabel('% Fragment remaining')
legend('Base case','Double','Ten fold','Hundred fold')
set(gca,'XTick',[1 2.4],'XTickLabels',{'Virus-specific','GBM & virus-specific'})
xtickangle(35)
b(1).CData = [0.31,0.54,0.68];
b(2).CData = [0.26,0.4,0.55];
b(3).CData = [0.93,0.46,0.39];
b(4).CData = [1,0.64,0.45];

%%
figure
subplot(1,2,1)
hold on 
plot(base_case.TH_number,'Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double.TH_number,'Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold.TH_number,'Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold.TH_number,'Color',[1,0.64,0.45],'LineWidth',2)
plot(dense_immune_double_killers.TH_number,'--','Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold_killers.TH_number,'--','Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold_killers.TH_number,'--','Color',[1,0.64,0.45],'LineWidth',2)
xlim([0 36])
ylabel('No. of THs')
xlabel('Time (hours)')
set(gca,'FontSize',19)
set(gca,'yscale','log')
subplot(1,2,2)
hold on 
plot(base_case.CTL_number,'Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double.CTL_number,'Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold.CTL_number,'Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold.CTL_number,'Color',[1,0.64,0.45],'LineWidth',2)
plot(dense_immune_double_killers.CTL_number,'--','Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold_killers.CTL_number,'--','Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold_killers.CTL_number,'--','Color',[1,0.64,0.45],'LineWidth',2)
xlim([0 36])
ylabel('No. of CTLs')
xlabel('Time (hours)')
set(gca,'FontSize',19)
legend('Base case','Double','Ten fold','Hundred fold','Double (GBM antigen)','Ten fold (GBM antigen)','Hundredful (GBM antigen)')
set(gca,'yscale','log')


figure
subplot(1,2,1)
hold on 
plot(base_case.TH_number./base_case.TH_number(1)*100,'Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double.TH_number./dense_immune_double.TH_number(1)*100,'Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold.TH_number./dense_immune_tenfold.TH_number(1)*100,'Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold.TH_number./dense_immune_hundredfold.TH_number(1)*100,'Color',[1,0.64,0.45],'LineWidth',2)
plot(dense_immune_double_killers.TH_number./dense_immune_double_killers.TH_number(1)*100,'--','Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold_killers.TH_number./dense_immune_tenfold_killers.TH_number(1)*100,'--','Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold_killers.TH_number./dense_immune_hundredfold_killers.TH_number(1)*100,'--','Color',[1,0.64,0.45],'LineWidth',2)
xlim([0 36])
ylabel('% initial THs')
xlabel('Time (hours)')
set(gca,'FontSize',19)
%set(gca,'yscale','log')
subplot(1,2,2)
hold on 
plot(base_case.CTL_number./base_case.CTL_number(1)*100,'Color',[0.31,0.54,0.68],'LineWidth',2)
plot(dense_immune_double.CTL_number/dense_immune_double.CTL_number(1)*100,'Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold.CTL_number/dense_immune_tenfold.CTL_number(1)*100,'Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold.CTL_number/dense_immune_hundredfold.CTL_number(1)*100,'Color',[1,0.64,0.45],'LineWidth',2)
plot(dense_immune_double_killers.CTL_number/dense_immune_double_killers.CTL_number(1)*100,'--','Color',[0.26,0.4,0.55],'LineWidth',2)
plot(dense_immune_tenfold_killers.CTL_number/dense_immune_tenfold_killers.CTL_number(1)*100,'--','Color',[0.93,0.46,0.39],'LineWidth',2)
plot(dense_immune_hundredfold_killers.CTL_number/dense_immune_hundredfold_killers.CTL_number(1)*100,'--','Color',[1,0.64,0.45],'LineWidth',2)
xlim([0 36])
ylabel('% initial CTLs')
xlabel('Time (hours)')
set(gca,'FontSize',19)
legend('Base case','Double','Ten fold','Hundred fold','Double (GBM antigen)','Ten fold (GBM antigen)','Hundredful (GBM antigen)')
%set(gca,'yscale','log')
