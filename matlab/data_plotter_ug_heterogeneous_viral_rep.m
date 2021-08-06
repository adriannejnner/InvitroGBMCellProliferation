Dense_ug_0point00002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Dense_ug_0point00002\Dense_ug_0point00002.mat');
Dense_ug_0point0002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Dense_ug_0point0002\Dense_ug_0point0002.mat');
Dense_ug_0point002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\us\Dense_us_0point01_rep1\Dense_us_0point01_rep1.mat');
Dense_ug_0point02 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Dense_ug_0point02\Dense_ug_0point02.mat');
Dense_ug_0point2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Dense_ug_0point2\Dense_ug_0point2.mat');

Sparse_ug_0point00002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Sparse_ug_0point00002\Sparse_ug_0point00002.mat');
Sparse_ug_0point0002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Sparse_ug_0point0002\Sparse_ug_0point0002.mat');
Sparse_ug_0point002 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\us\Sparse_us_0point01_rep1\Sparse_us_0point01_rep1.mat');
Sparse_ug_0point02 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Sparse_ug_0point02\Sparse_ug_0point02.mat');
Sparse_ug_0point2 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\Heterogeneous_rep_sims\ug\Sparse_ug_0point2\Sparse_ug_0point2.mat');

dense_initial_cell_number = Dense_ug_0point00002.live_cells(1);
sparse_initial_cell_number = Sparse_ug_0point00002.live_cells(1);


figure
hold on 
plot(linspace(0,72,36),Dense_ug_0point00002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point0002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point02.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point2.dead_cells(1:36),'LineWidth',2)
legend('u_g = 0.00002','u_g = 0.0002','u_g = 0.002','u_g = 0.02','u_g = 0.2')
xlabel('Time (hours)')
ylabel('Dead cells')
set(gca,'FontSize',18)
title('Dense')



figure
hold on 
plot(linspace(0,72,36),Dense_ug_0point00002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point0002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point02.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Dense_ug_0point2.live_cells(1:36),'LineWidth',2)
legend('u_g = 0.00002','u_g = 0.0002','u_g = 0.002','u_g = 0.02','u_g = 0.2')
xlabel('Time (hours)')
ylabel('GBM cells (uninfected+infected)')
set(gca,'FontSize',18)
title('Dense')


figure
hold on 
plot(linspace(0,72,36),Sparse_ug_0point00002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point0002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point002.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point02.dead_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point2.dead_cells(1:36),'LineWidth',2)
legend('u_g = 0.00002','u_g = 0.0002','u_g = 0.002','u_g = 0.02','u_g = 0.2')
xlabel('Time (hours)')
ylabel('Dead cells')
set(gca,'FontSize',18)
title('Sparse')



figure
hold on 
plot(linspace(0,72,36),Sparse_ug_0point00002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point0002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point002.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point02.live_cells(1:36),'LineWidth',2)
plot(linspace(0,72,36),Sparse_ug_0point2.live_cells(1:36),'LineWidth',2)
legend('u_g = 0.00002','u_g = 0.0002','u_g = 0.002','u_g = 0.02','u_g = 0.2')
xlabel('Time (hours)')
ylabel('GBM cells (uninfected+infected)')
set(gca,'FontSize',18)
title('Sparse')


figure
hold on 
yyaxis left
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Dense_ug_0point00002.live_cells(36),...
    Dense_ug_0point0002.live_cells(36),...
    Dense_ug_0point002.live_cells(36),...
    Dense_ug_0point02.live_cells(36),...
    Dense_ug_0point2.live_cells(36)]/dense_initial_cell_number*100,':','LineWidth',2)
ylabel('% Initial GBM cells ')
yyaxis right
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Sparse_ug_0point00002.live_cells(36),...
    Sparse_ug_0point0002.live_cells(36),...
    Sparse_ug_0point002.live_cells(36),...
    Sparse_ug_0point02.live_cells(36),...
    Sparse_ug_0point2.live_cells(36)]/sparse_initial_cell_number*100,':','LineWidth',2)
set(gca,'xscale','log')
xlabel('Stroma uptake rate, u_s')
ylabel('% Initial GBM cells ')
set(gca,'FontSize',18)
legend('Dense','Sparse')


figure
hold on 
yyaxis left
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Dense_ug_0point00002.live_cells(36),...
    Dense_ug_0point0002.live_cells(36),...
    Dense_ug_0point002.live_cells(36),...
    Dense_ug_0point02.live_cells(36),...
    Dense_ug_0point2.live_cells(36)],':','LineWidth',2)
ylabel('Number GBM cells at 72 hrs ')
yyaxis right
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Sparse_ug_0point00002.live_cells(36),...
    Sparse_ug_0point0002.live_cells(36),...
    Sparse_ug_0point002.live_cells(36),...
    Sparse_ug_0point02.live_cells(36),...
    Sparse_ug_0point2.live_cells(36)],':','LineWidth',2)
set(gca,'xscale','log')
xlabel('Stroma uptake rate, u_s')
ylabel('Number GBM cells at 72 hrs')
set(gca,'FontSize',18)
legend('Dense','Sparse')

figure
hold on 
yyaxis left
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Dense_ug_0point00002.uninfected_cells(36),...
    Dense_ug_0point0002.uninfected_cells(36),...
    Dense_ug_0point002.uninfected_cells(36),...
    Dense_ug_0point02.uninfected_cells(36),...
    Dense_ug_0point2.uninfected_cells(36)]./dense_initial_cell_number*100,':','LineWidth',2)
ylabel('% initial GBM cells uninfected')
yyaxis right
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Sparse_ug_0point00002.uninfected_cells(36),...
    Sparse_ug_0point0002.uninfected_cells(36),...
    Sparse_ug_0point002.uninfected_cells(36),...
    Sparse_ug_0point02.uninfected_cells(36),...
    Sparse_ug_0point2.uninfected_cells(36)]./sparse_initial_cell_number*100,':','LineWidth',2)
set(gca,'xscale','log')
xlabel('Stroma uptake rate, u_s')
ylabel('% initial GBM cells uninfected')
set(gca,'FontSize',18)
legend('Dense','Sparse')


figure
hold on 
yyaxis left
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Dense_ug_0point00002.infected_cells(36),...
    Dense_ug_0point0002.infected_cells(36),...
    Dense_ug_0point002.infected_cells(36),...
    Dense_ug_0point02.infected_cells(36),...
    Dense_ug_0point2.infected_cells(36)],':','LineWidth',2)
ylabel('No. infected cells (virus>1)')
yyaxis right
plot([0.00002 0.0002 0.002 0.02 0.2],[...
    Sparse_ug_0point00002.infected_cells(36),...
    Sparse_ug_0point0002.infected_cells(36),...
    Sparse_ug_0point002.infected_cells(36),...
    Sparse_ug_0point02.infected_cells(36),...
    Sparse_ug_0point2.infected_cells(36)],':','LineWidth',2)
set(gca,'xscale','log')
xlabel('Stroma uptake rate, u_s')
ylabel('No. infected cells (virus>1)')
set(gca,'FontSize',18)
legend('Dense','Sparse')
