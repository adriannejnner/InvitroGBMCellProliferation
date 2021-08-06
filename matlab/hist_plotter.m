tcount = 1;

%heterogeneity_R_0point5
%heteorogeneity_R_1point5
%heterogeneity_R_10

A = 'output0000000';
B = '.xml';
K = [A num2str(tcount-1,'%d') B];


MCDS = read_MultiCellDS_xml(K);


GBMcells = find(MCDS.discrete_cells.metadata.type == 2);
STROMAcells = find(MCDS.discrete_cells.metadata.type == 4);

figure
[n,xout] = hist(MCDS.discrete_cells.custom.special_virus_replication_rate(GBMcells),100);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
xlabel('Virus replication rate \nu')
ylabel('Frequency')
set(gca,'FontSize',15)
title('GBM')

figure
[n,xout] = hist(MCDS.discrete_cells.custom.special_virus_uptakerate(GBMcells),100);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
xlabel('GBM cell virus uptake rate u_g')
ylabel('Frequency')
set(gca,'FontSize',15)
title('GBM')

figure
[n,xout] = hist(MCDS.discrete_cells.custom.special_virus_uptakerate(STROMAcells),100);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
xlabel('Stroma cell virus uptake rate u_s')
ylabel('Frequency')
set(gca,'FontSize',15)
title('Stroma')

