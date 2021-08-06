base_case_dense = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\us_0point01\dense_homogeneous_us_0point01.mat');
base_case_sparse = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\us_0point01\sparse_homogeneous_us_0point01.mat');
base_case_patch_40_60 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\patchy\proportions\40percent\patch_40_60.mat');
base_case_sparse = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\us_0point01\sparse_homogeneous_us_0point01.mat');
centre_injection_dense = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Centre injection\Dense\Centre_injection_dense.mat');
centre_injection_sparse = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Centre injection\Sparse\Centre_injection_sparse.mat');
GBM_cells_only = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\GBM cells only\GBM_cells_only.mat');
patchy_injections_sparse = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\patchy\patch_injections_sparse\patchy_injections_sparse.mat');
patchy_injections_dense = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\patchy\patch_injections_dense\patchy_injections_dense.mat');

colmap = [102,194,165
252,141,98
141,160,203
231,138,195
166,216,84]/255;

dense_initial_cell_number = base_case_dense.live_cells(1);
sparse_initial_cell_number = base_case_sparse.live_cells(1);

figure
hold on 
plot(base_case_dense.live_cells/dense_initial_cell_number*100,'Color',colmap(1,:),'LineWidth',2)
plot(base_case_sparse.live_cells/sparse_initial_cell_number*100,'Color',colmap(2,:),'LineWidth',2)
plot(base_case_patch_40_60.live_cells/base_case_patch_40_60.live_cells(1)*100,'Color',colmap(3,:),'LineWidth',2)
plot(patchy_injections_sparse.live_cells/patchy_injections_sparse.live_cells(1)*100,'Color',colmap(4,:),'LineWidth',2)
plot(patchy_injections_dense.live_cells/patchy_injections_dense.live_cells(1)*100,'Color',colmap(5,:),'LineWidth',2)
legend('Fully dense (periphery)','Fully sparse (periphery)','Patchy 40:60 (periphery)','Patchy 40:60 (sparse patch centre)','Patchy 40:60 (dense patch centre)')
xlabel('Time (hours)')
ylabel('% GBM cells (uninfected + infected)')
set(gca,'FontSize',19)
set(gca,'XTick',[0 12 24 36],'XTickLabels',{'0','24','48','72'})  
xlim([0 36])


figure
hold on 
plot(base_case_dense.live_cells/dense_initial_cell_number*100,'Color',colmap(1,:),'LineWidth',2)
plot(base_case_sparse.live_cells/sparse_initial_cell_number*100,'Color',colmap(2,:),'LineWidth',2)
plot(centre_injection_dense.live_cells/dense_initial_cell_number*100,'Color',colmap(3,:),'LineWidth',2)
plot(centre_injection_sparse.live_cells/sparse_initial_cell_number*100,'Color',colmap(4,:),'LineWidth',2)
plot(GBM_cells_only.live_cells/GBM_cells_only.live_cells(1)*100,'Color',colmap(5,:),'LineWidth',2)
legend('Base case (dense)','Base case (sparse)','Centre injection (dense)','Centre injection (sparse)','GBM cells only')
xlabel('Time (hours)')
ylabel('% GBM cells (uninfected + infected)')
set(gca,'FontSize',19)
set(gca,'XTick',[0 12 24 36],'XTickLabels',{'0','24','48','72'})  
xlim([0 36])


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

