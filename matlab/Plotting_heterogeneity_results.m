
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

tcount = 13;
K = [A2 num2str(tcount-1,'%d') B];
MCDS = read_MultiCellDS_xml(K);
% cancer cell = 2, stroma cell = 4

live_cells_numDr00001(1) = length(MCDS.discrete_cells.live_cells); % number live cells
infected_cells_numDr00001(1) = length(intersect(find(MCDS.discrete_cells.custom.intracellular_virus_amount >1),MCDS.discrete_cells.live_cells)); %infected cell number
dead_cells_numDr00001(1) = length(MCDS.discrete_cells.dead_cells); % number dead cells

tcount = 25;
K = [A2 num2str(tcount-1,'%d') B];
MCDS = read_MultiCellDS_xml(K);
% cancer cell = 2, stroma cell = 4

live_cells_numDr00001(2) = length(MCDS.discrete_cells.live_cells); % number live cells
infected_cells_numDr00001(2) = length(intersect(find(MCDS.discrete_cells.custom.intracellular_virus_amount >1),MCDS.discrete_cells.live_cells)); %infected cell number
dead_cells_numDr00001(2) = length(MCDS.discrete_cells.dead_cells); % number dead cells

tcount = 37;   
K = [A2 num2str(tcount-1,'%d') B];
MCDS = read_MultiCellDS_xml(K);
% cancer cell = 2, stroma cell = 4

live_cells_numDr00001(3) = length(MCDS.discrete_cells.live_cells); % number live cells
infected_cells_numDr00001(3) = length(intersect(find(MCDS.discrete_cells.custom.intracellular_virus_amount >1),MCDS.discrete_cells.live_cells)); %infected cell number
dead_cells_numDr00001(3) = length(MCDS.discrete_cells.dead_cells); % number dead cells

save('sparse_us_00001.mat','live_cells_numDr00001','infected_cells_numDr00001','dead_cells_numDr00001')

%%
load('sparse_us_001.mat')
load('sparse_us_0001.mat')
load('sparse_us_00001.mat')

figure
bar([live_cells_numDr001;live_cells_numDr0001;live_cells_numDr00001]);  

figure
bar([infected_cells_numDr001;infected_cells_numDr0001;infected_cells_numDr00001]); 

figure
bar([dead_cells_numDr001;dead_cells_numDr0001;dead_cells_numDr00001]);  

figure
bar([live_cells_numDr001(end) live_cells_numDr0001(end) live_cells_numDr00001(end);infected_cells_numDr001(end) infected_cells_numDr0001(end) infected_cells_numDr00001(end); dead_cells_numDr001(end) dead_cells_numDr0001(end) dead_cells_numDr00001(end)]);
set(gca,'XtickLabels',{'Live cells','Infected cells','Dead cells'})
ylabel('Number of cells')
set(gca,'FontSize',15)

figure
bar([live_cells_numDr001(1) live_cells_numDr0001(1) live_cells_numDr00001(1);infected_cells_numDr001(1) infected_cells_numDr0001(1) infected_cells_numDr00001(1); dead_cells_numDr001(1) dead_cells_numDr0001(1) dead_cells_numDr00001(1)]);
set(gca,'XtickLabels',{'Live cells','Infected cells','Dead cells'})
ylabel('Number of cells')
set(gca,'FontSize',15)

figure
hold on 
plot(ind0,'LineWidth',2)
xlabel('Hours')
ylabel('Live cells')
set(gca,'FontSize',15,'xtick',[0 12 24 36],'xticklabels',{'0','24','48','72'})
xlim([0 36])

figure
hold on 
plot(dead,'LineWidth',2)
xlabel('Hours')
ylabel('Dead cells')
set(gca,'FontSize',15,'xtick',[0 12 24 36],'xticklabels',{'0','24','48','72'})
xlim([0 36])

%% dosage
clear all
 load('tri_D0_increaseby10.mat')
 load('tri_2inj_0and1day.mat')
 load('tri_3inj_500and1000.mat')
 load('tri_3inj_500and1000.mat')

 timetotal =96;
 figure
 hold on
 plot([1:timetotal]*60,tri_D0_increaseby10,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_2inj_0and1day,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_3inj_500and1000,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_3inj_500and1000,'o-','LineWidth',1,'MarkerSize',3)


 
plot(tri_Dv_10)