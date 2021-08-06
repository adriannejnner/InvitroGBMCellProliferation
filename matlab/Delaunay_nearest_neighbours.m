%% d
points_at_edge_less_less = [];

tcount = 37;
A2 = 'output000000';
B = '.xml';
K = [A2 num2str(tcount-1,'%d') B];
MCDS = read_MultiCellDS_xml(K);
    
P = MCDS.discrete_cells.custom.intracellular_virus_amount;
locs_infected = find( P > 0.5); %finds location of cells with more than 1 virus inside
locs_GBM = find(MCDS.discrete_cells.metadata.type==2); %finds location
locs_GBM_alive = intersect(MCDS.discrete_cells.live_cells,locs_GBM);
locs_infected_GBM_alive = intersect(locs_infected,locs_GBM_alive);
locs_infected_GBM_alive_and_dead = intersect(locs_infected,locs_GBM);



P_sub = MCDS.discrete_cells.custom.intracellular_virus_amount(MCDS.discrete_cells.live_cells);
locs_infected_sub = find( P > 1); %finds location of cells with more than 1 virus inside
locs_GBM_alive_sub = find(MCDS.discrete_cells.metadata.type(MCDS.discrete_cells.live_cells)==2); %finds location
locs_infected_GBM_alive_sub = intersect(locs_infected_sub,locs_GBM_alive_sub);


%create delaunay triangulation
DT = delaunay(MCDS.discrete_cells.state.position(:,1:2));

%for each live GBM cell determine whether it has a nearest neighbour that is
%uninfected
for ii = 1:length(locs_infected_GBM_alive) %iterate through the infected GBM cell locations

    [row_neighbours col_neighbours] = find(DT==locs_infected_GBM_alive(ii));%find location of infected cell index in DT

    for ncount = 1:length(row_neighbours) %check each infected cell neighbour to see if it is uninfected
        neighbour_index = DT(row_neighbours(ncount),col_neighbours(ncount));
        if MCDS.discrete_cells.metadata.type(neighbour_index) == 2 && MCDS.discrete_cells.custom.intracellular_virus_amount(neighbour_index)<=10 %&& isempty(find(MCDS.discrete_cells.live_cells==neighbour_index))==0
            disp('uninfected cell neighbour')
            points_at_edge_less_less = [points_at_edge_less_less;MCDS.discrete_cells.state.position(locs_infected_GBM_alive(ii),1:2)];
            ncount = length(row_neighbours);
        end
    end
    
end

points_at_edge_less_less = unique(points_at_edge_less_less,'rows');
 
save('sparse_heterogeneous_rep_infiltrating_edge_sig_0point01.mat','points_at_edge_less_less')

%%
dense_homogeneous = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\us_0point01\dense_homogeneous.mat');
dense_heterogeneous_rep_infiltrating_edge_sig_0point0001 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\heterogeneous_rep_0point0001\dense_heterogeneous_rep_infiltrating_edge_sig_0point0001.mat');
dense_heterogeneous_rep_infiltrating_edge_sig_0point001 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\heterogeneous_rep_0point001\dense_heterogeneous_rep_infiltrating_edge_sig_0point001.mat');
dense_heterogeneous_rep_infiltrating_edge_sig_0point01 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\dense\heterogeneous_rep_0point01\dense_heterogeneous_rep_infiltrating_edge_sig_0point01.mat');

distance2 = sqrt(dense_heterogeneous_rep_infiltrating_edge_sig_0point0001.points_at_edge_less_less(:,1).^2+dense_heterogeneous_rep_infiltrating_edge_sig_0point0001.points_at_edge_less_less(:,2).^2);
distance3 = sqrt(dense_heterogeneous_rep_infiltrating_edge_sig_0point001.points_at_edge_less_less(:,1).^2+dense_heterogeneous_rep_infiltrating_edge_sig_0point001.points_at_edge_less_less(:,2).^2);
distance4 = sqrt(dense_heterogeneous_rep_infiltrating_edge_sig_0point01.points_at_edge_less_less(:,1).^2+dense_heterogeneous_rep_infiltrating_edge_sig_0point01.points_at_edge_less_less(:,2).^2);
distance1 = sqrt(dense_homogeneous.points_at_edge_less_less(:,1).^2+dense_homogeneous.points_at_edge_less_less(:,2).^2);

figure
hold on
histogram(distance2,'FaceAlpha',0.5)
histogram(distance3,'FaceAlpha',0.5)
histogram(distance4,'FaceAlpha',0.5)
histogram(distance1,'FaceAlpha',0.5)
xlabel('\mu from periphery')
ylabel('Frequency')
legend('\sigma = 0.0001','\sigma = 0.001','\sigma = 0.01','Homogenous')
set(gca,'FontSize',18')

%%
sparse_homogeneous = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\us_0point01\sparse_homogeneous.mat');
sparse_heterogeneous_rep_infiltrating_edge_sig_0point0001 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\heterogeneity_rep_0poin0001\sparse_heterogeneous_rep_infiltrating_edge_sig_0point0001.mat');
sparse_heterogeneous_rep_infiltrating_edge_sig_0point001 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\heterogeneity_rep_0poin001\sparse_heterogeneous_rep_infiltrating_edge_sig_0point001.mat');
sparse_heterogeneous_rep_infiltrating_edge_sig_0point01 = load('C:\Users\adria\PhysiCell_V.1.6.0 GBM-OV project\PhysiCell-GBM-OV project\matlab\sparse\heterogeneity_rep_0poin01\sparse_heterogeneous_rep_infiltrating_edge_sig_0point01.mat');

distance2 = sqrt(sparse_heterogeneous_rep_infiltrating_edge_sig_0point0001.points_at_edge_less_less(:,1).^2+sparse_heterogeneous_rep_infiltrating_edge_sig_0point0001.points_at_edge_less_less(:,2).^2);
distance3 = sqrt(sparse_heterogeneous_rep_infiltrating_edge_sig_0point001.points_at_edge_less_less(:,1).^2+sparse_heterogeneous_rep_infiltrating_edge_sig_0point001.points_at_edge_less_less(:,2).^2);
distance4 = sqrt(sparse_heterogeneous_rep_infiltrating_edge_sig_0point01.points_at_edge_less_less(:,1).^2+sparse_heterogeneous_rep_infiltrating_edge_sig_0point01.points_at_edge_less_less(:,2).^2);
distance1 = sqrt(sparse_homogeneous.points_at_edge_less_less(:,1).^2+sparse_homogeneous.points_at_edge_less_less(:,2).^2);

figure
hold on
histogram(distance2,'FaceAlpha',0.5)
histogram(distance3,'FaceAlpha',0.5)
histogram(distance4,'FaceAlpha',0.5)
histogram(distance1,'FaceAlpha',0.5)
xlabel('\mu from periphery')
ylabel('Frequency')
legend('\sigma = 0.0001','\sigma = 0.001','\sigma = 0.01','Homogenous')
set(gca,'FontSize',18')

