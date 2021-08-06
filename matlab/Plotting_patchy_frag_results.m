
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

for tcount = 2:37
    if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    elseif tcount<1001
        K = [A3 num2str(tcount-1,'%d') B];
    else
        K = [A4 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    P = MCDS.discrete_cells.custom.intracellular_virus_amount;
    % cancer cell = 2, stroma cell = 4

     locs_infected = find( P > 1); %finds location of cells with more than 1 virus inside
     number_infected = length(locs_infected); %gives the total number of cells with more than 1 virus inside

     position_infected = MCDS.discrete_cells.state.position(locs_infected,1:2);

     if isempty(position_infected)==1
        distance_to_edge(tcount) = 0; 
     else
          Radius = 1270;
         distance_to_edge(tcount) = max(Radius-sqrt(position_infected(:,1).^2+position_infected(:,2).^2));
     end
 
     if tcount == 37
            distance_to_center = sqrt(position_infected(:,1).^2+position_infected(:,2).^2);
            band1_loc = find(distance_to_center>Radius-50);
            band2_loc = find(distance_to_center>Radius-2*50 & distance_to_center<=Radius-50);
            band3_loc = find(distance_to_center>Radius-3*50 & distance_to_center<=Radius-2*50);
            band4_loc = find(distance_to_center>Radius-4*50 & distance_to_center<=Radius-3*50);
            band5_loc = find(distance_to_center>Radius-5*50 & distance_to_center<=Radius-4*50);
            band6_loc = find(distance_to_center>Radius-6*50 & distance_to_center<=Radius-5*50);
            band7_loc = find(distance_to_center>Radius-7*50 & distance_to_center<=Radius-6*50);
            band8_loc = find(distance_to_center>Radius-8*50 & distance_to_center<=Radius-7*50);
            band9_loc = find(distance_to_center>Radius-9*50 & distance_to_center<=Radius-8*50);
            band10_loc = find(distance_to_center>Radius-10*50 & distance_to_center<=Radius-9*50);
            band11_loc = find(distance_to_center>Radius-11*50 & distance_to_center<=Radius-10*50);
            band12_loc = find(distance_to_center>Radius-12*50 & distance_to_center<=Radius-11*50);
            band13_loc = find(distance_to_center>Radius-13*50 & distance_to_center<=Radius-12*50);
            band14_loc = find(distance_to_center>Radius-14*50 & distance_to_center<=Radius-13*50);
            band15_loc = find(distance_to_center>Radius-15*50 & distance_to_center<=Radius-14*50);
            band16_loc = find(distance_to_center>Radius-16*50 & distance_to_center<=Radius-15*50);

            infiltration_hete_uptake_sig_0point01 = [sum(P(band1_loc)), sum(P(band2_loc)), sum(P(band3_loc)),...
                sum(P(band4_loc)), sum(P(band5_loc)), sum(P(band6_loc)), sum(P(band7_loc)),...
                sum(P(band8_loc)), sum(P(band9_loc)), sum(P(band10_loc)), sum(P(band11_loc)),...
                sum(P(band12_loc)), sum(P(band13_loc)), sum(P(band14_loc)), sum(P(band15_loc)), sum(P(band16_loc))];
            
            locs_stroma = find(MCDS.discrete_cells.metadata.type==4);
            common = intersect(locs_infected,locs_stroma);
            locs_GBM = find(MCDS.discrete_cells.metadata.type==2);
            interr = intersect(locs_infected,locs_GBM);
            position_infected_notstroma =  MCDS.discrete_cells.state.position(interr,1:2);% MCDS.discrete_cells.state.position(setxor(locs_infected,common),1:2);
            distance_to_center = sqrt(position_infected_notstroma(:,1).^2+position_infected_notstroma(:,2).^2);
            band1_loc = find(distance_to_center>Radius-50);
            band2_loc = find(distance_to_center>Radius-2*50 & distance_to_center<=Radius-50);
            band3_loc = find(distance_to_center>Radius-3*50 & distance_to_center<=Radius-2*50);
            band4_loc = find(distance_to_center>Radius-4*50 & distance_to_center<=Radius-3*50);
            band5_loc = find(distance_to_center>Radius-5*50 & distance_to_center<=Radius-4*50);
            band6_loc = find(distance_to_center>Radius-6*50 & distance_to_center<=Radius-5*50);
            band7_loc = find(distance_to_center>Radius-7*50 & distance_to_center<=Radius-6*50);
            band8_loc = find(distance_to_center>Radius-8*50 & distance_to_center<=Radius-7*50);
            band9_loc = find(distance_to_center>Radius-9*50 & distance_to_center<=Radius-8*50);
            band10_loc = find(distance_to_center>Radius-10*50 & distance_to_center<=Radius-9*50);
            band11_loc = find(distance_to_center>Radius-11*50 & distance_to_center<=Radius-10*50);
            band12_loc = find(distance_to_center>Radius-12*50 & distance_to_center<=Radius-11*50);
            band13_loc = find(distance_to_center>Radius-13*50 & distance_to_center<=Radius-12*50);
            band14_loc = find(distance_to_center>Radius-14*50 & distance_to_center<=Radius-13*50);
            band15_loc = find(distance_to_center>Radius-15*50 & distance_to_center<=Radius-14*50);
            band16_loc = find(distance_to_center>Radius-16*50 & distance_to_center<=Radius-15*50);

            infiltration_notstroma_hete_uptake_sig_0point01 = [sum(P(interr(band1_loc))), sum(P(interr(band2_loc))), sum(P(interr(band3_loc))),...
                sum(P(interr(band4_loc))), sum(P(interr(band5_loc))), sum(P(interr(band6_loc))), sum(P(interr(band7_loc))),...
                sum(P(interr(band8_loc))), sum(P(interr(band9_loc))), sum(P(interr(band10_loc))), sum(P(interr(band11_loc))),...
                sum(P(interr(band12_loc))), sum(P(interr(band13_loc))), sum(P(interr(band14_loc))), sum(P(interr(band15_loc))), sum(P(interr(band16_loc)))];
     end
    
    live_cells_hete_uptake_sig_0point01(tcount-1) = length(MCDS.discrete_cells.live_cells); % number live cells
    infected_cells_hete_uptake_sig_0point01(tcount-1) = length(intersect(find(MCDS.discrete_cells.custom.intracellular_virus_amount >1),MCDS.discrete_cells.live_cells)); %infected cell number
    dead_cells_hete_uptake_sig_0point01(tcount-1) = length(MCDS.discrete_cells.dead_cells); % number dead cells

end

save('hete_uptake_sig_0point01','infiltration_hete_uptake_sig_0point01','infiltration_notstroma_hete_uptake_sig_0point01',...
    'live_cells_hete_uptake_sig_0point01','infected_cells_hete_uptake_sig_0point01','dead_cells_hete_uptake_sig_0point01')

%%
load('hete_rep_sig_0point0001.mat')
load('hete_rep_sig_0point01.mat')
load('hete_ug_sig_0point0001.mat')
load('hete_ug_sig_0point01.mat')
load('hete_uptake_sig_0point0001.mat')
load('hete_uptake_sig_0point01.mat')
load('hete_us_0point01_sig_0point0001.mat')
load('hete_us_0point01_sig_0point01.mat')
load('hete_us_0point1_sig_0point0001.mat')
load('hete_us_0point1_sig_0point01.mat')
load('homo_us_0point01.mat')
load('homo_us_0point1.mat')
load('homo_us_1.mat')

figure
hold on 
plot(live_cells_hete_rep_sig_0point0001./17451*100,'LineWidth',2)
plot(live_cells_hete_rep_sig_0point01./17451*100,'LineWidth',2)
plot(live_cells_hete_ug_sig_0point0001./17451*100,'LineWidth',2)
plot(live_cells_hete_ug_sig_0point01./17451*100,'LineWidth',2)
plot(live_cells_hete_uptake_sig_0point0001./17451*100,'LineWidth',2)
plot(live_cells_hete_uptake_sig_0point01./17451*100,'LineWidth',2)
legend('Hete \nu, \sigma = 0.0001','Hete \nu, \sigma = 0.01','Hete u_g, \sigma = 0.0001','Hete u_g \sigma = 0.01',...
    'Hete u_s & u_g, \sigma = 0.0001','Hete u_s & u_g, \sigma = 0.01')
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'XTick',linspace(1,37,5),'XTickLabels',{'0','18','36','54','72'})
set(gca,'FontSize',16)

figure 
hold on
plot(live_cells_hete_us_0point01_sig_0point0001./17451*100,'LineWidth',2)
plot(live_cells_hete_us_0point01_sig_0point01./17451*100,'LineWidth',2)
plot(live_cells_hete_us_0point1_sig_0point0001./17451*100,'LineWidth',2)
plot(live_cells_hete_us_0point1_sig_0point01./17451*100,'LineWidth',2)
legend('Hete u_s = 0.01 \sigma = 0.0001','Hete u_s = 0.01 \sigma = 0.01',...
    'Hete u_s = 0.1 \sigma = 0.0001','Hete u_s = 0.1 \sigma = 0.01')
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'XTick',linspace(1,37,5),'XTickLabels',{'0','18','36','54','72'})
set(gca,'FontSize',16)

figure
hold on
plot(live_cells_homo_us_0point01./17451*100,'LineWidth',2)
plot(live_cells_homo_us_0point1./17451*100,'LineWidth',2)
plot(live_cells_homo_us_1./17451*100,'LineWidth',2)
legend('Homo u_s = 0.01','Homo u_s = 0.1','Homo u_s = 1')
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'XTick',linspace(1,37,5),'XTickLabels',{'0','18','36','54','72'})
set(gca,'FontSize',16)

figure
hold on
bar([infiltration_notstroma_hete_rep_sig_0point0001;...
   infiltration_notstroma_hete_rep_sig_0point01]')
legend('Hete \nu, \sigma = 0.0001','Hete \nu, \sigma = 0.01')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)

figure
hold on
bar([infiltration_notstroma_hete_ug_sig_0point0001;...
   infiltration_notstroma_hete_ug_sig_0point01;...
   infiltration_notstroma_hete_uptake_sig_0point0001;...
   infiltration_notstroma_hete_uptake_sig_0point01]')
legend('Hete u_g, \sigma = 0.0001','Hete u_g, \sigma = 0.01',...
    'Hete u_s & u_g, \sigma = 0.0001', 'Hete u_s & u_g, \sigma = 0.01')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)

figure
hold on
bar([infiltration_notstroma_hete_us_0point01_sig_0point0001;...
   infiltration_notstroma_hete_us_0point01_sig_0point01;...
   infiltration_notstroma_hete_us_0point1_sig_0point0001;...
   infiltration_notstroma_hete_us_0point1_sig_0point01]')
legend('Hete u_s = 0.01 \sigma = 0.0001','Hete u_s = 0.01 \sigma = 0.01',...
    'Hete u_s = 0.1 \sigma = 0.0001','Hete u_s = 0.1 \sigma = 0.01')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)


figure
hold on
bar([infiltration_notstroma_homo_us_0point01;...
   infiltration_notstroma_homo_us_0point1;...
   infiltration_notstroma_homo_us_1]')
legend('Homo u_s = 0.01','Homo u_s = 0.1','Homo u_s = 1')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)
