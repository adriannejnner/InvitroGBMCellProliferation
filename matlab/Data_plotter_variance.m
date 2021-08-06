stroma_uptake_rate1_1 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run1\stroma_uptake_rate1_1.mat');
stroma_uptake_rate1_2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run2\stroma_uptake_rate1_2.mat');
stroma_uptake_rate1_3 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run3\stroma_uptake_rate1_3.mat');
stroma_uptake_rate1_4 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run4\stroma_uptake_rate1_4.mat');
stroma_uptake_rate1_5 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run5\stroma_uptake_rate1_5.mat');
stroma_uptake_rate1_6 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\MC\StromaUptake1run6\stroma_uptake_rate1_6.mat');

initial_cell_number = stroma_uptake_rate1_1.live_cells(1);

rep_mat_live_cells = [stroma_uptake_rate1_1.live_cells;...
    stroma_uptake_rate1_2.live_cells;... 
    stroma_uptake_rate1_3.live_cells;... 
    stroma_uptake_rate1_4.live_cells;... 
    stroma_uptake_rate1_5.live_cells;... 
    stroma_uptake_rate1_6.live_cells]./initial_cell_number*100;

mean_vec_live_cells = mean(rep_mat_live_cells);
std_vec_live_cells = std(rep_mat_live_cells);

fig1 = figure;
options.handle = fig1;
options.alpha = 0.1;
options.line_width = 1.5;
options.error = 'std';
options.x_axis = linspace(0,37,36);

options.color_area = 'y';
options.color_line = 'y';
plot_areaerrorbar(rep_mat_live_cells,options)

rep_mat_infected_cells = [stroma_uptake_rate1_1.infected_cells;...
    stroma_uptake_rate1_2.infected_cells;... 
    stroma_uptake_rate1_3.infected_cells;... 
    stroma_uptake_rate1_4.infected_cells;... 
    stroma_uptake_rate1_5.infected_cells;... 
    stroma_uptake_rate1_6.infected_cells];

mean_vec_infected_cells = mean(rep_mat_infected_cells);
std_vec_infected_cells = std(rep_mat_infected_cells);

colmap = [255,255,204
199,233,180
127,205,187
65,182,196
29,145,192
34,94,168
12,44,132]/255;


figure
hold on 
plot(stroma_uptake_rate1_1.live_cells/initial_cell_number*100,'Color',colmap(1,:),'LineWidth',1)
plot(stroma_uptake_rate1_2.live_cells/initial_cell_number*100,'Color',colmap(2,:),'LineWidth',1)
meanplot(stroma_uptake_rate1_4.live_cells/initial_cell_number*100,'Color',colmap(4,:),'LineWidth',1)
plot(stroma_uptake_rate1_5.live_cells/initial_cell_number*100,'Color',colmap(5,:),'LineWidth',1)
plot(stroma_uptake_rate1_6.live_cells/initial_cell_number*100,'Color',colmap(6,:),'LineWidth',1)
legend('u_s = 0','u_s = 0.00005','u_s = 0.0001','u_s = 0.001','u_s = 0.01','u_s = 0.1','u_s = 1')

figure
hold on 
b = bar([stroma_uptake_rate1_1.live_cells(end)/initial_cell_number*100;... 
	stroma_uptake_rate1_2.live_cells(end)/initial_cell_number*100;... 
	stroma_uptake_rate1_3.live_cells(end)/initial_cell_number*100;...
	stroma_uptake_rate1_4.live_cells(end)/initial_cell_number*100;...
	stroma_uptake_rate1_5.live_cells(end)/initial_cell_number*100;...
	stroma_uptake_rate1_6.live_cells(end)/initial_cell_number*100],'FaceColor','flat')
set(gca,'FontSize',19)
ylabel('% Fragment remaining')
legend('Base case','Double','Ten fold','Hundred fold')
set(gca,'XTick',[1 2.4],'XTickLabels',{'Virus-specific','GBM & virus-specific'})
xtickangle(35)

