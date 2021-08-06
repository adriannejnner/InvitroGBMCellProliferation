timetotal =36;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';
 
 for tcount = 2:timetotal+1
   
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    
 MCDS = read_MultiCellDS_xml(K); 
 P = MCDS.discrete_cells.custom.intracellular_virus_amount;
 
 locs_infected = find( P > 1); %finds location of cells with more than 1 virus inside
 number_infected = length(locs_infected); %gives the total number of cells with more than 1 virus inside
 
 position_infected = MCDS.discrete_cells.state.position(locs_infected,1:2);
 
 if isempty(position_infected)==1
    distance_to_edge(tcount) = 0; 
 else
      Radius = 1270;
     distance_to_edge(tcount) = max(Radius-sqrt(position_infected(:,1).^2+position_infected(:,2).^2));
 end
 
     if tcount == 12 || tcount == 24 || tcount == 36
         for ii = 1:length(locs_infected)
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

            infiltration(tcount/(12),:) = [sum(P(band1_loc)), sum(P(band2_loc)), sum(P(band3_loc)),...
                sum(P(band4_loc)), sum(P(band5_loc)), sum(P(band6_loc)), sum(P(band7_loc)),...
                sum(P(band8_loc)), sum(P(band9_loc)), sum(P(band10_loc)), sum(P(band11_loc)),...
                sum(P(band12_loc)), sum(P(band13_loc)), sum(P(band14_loc)), sum(P(band15_loc)), sum(P(band16_loc))];
            
            locs_stroma = find(MCDS.discrete_cells.metadata.type==4);
            common = intersect(locs_infected,locs_stroma);
            position_infected_notstroma =  MCDS.discrete_cells.state.position(setxor(locs_infected,common),1:2);
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

            infiltration_notstroma(tcount/(12),:) = [sum(P(band1_loc)), sum(P(band2_loc)), sum(P(band3_loc)),...
                sum(P(band4_loc)), sum(P(band5_loc)), sum(P(band6_loc)), sum(P(band7_loc)),...
                sum(P(band8_loc)), sum(P(band9_loc)), sum(P(band10_loc)), sum(P(band11_loc)),...
                sum(P(band12_loc)), sum(P(band13_loc)), sum(P(band14_loc)), sum(P(band15_loc)), sum(P(band16_loc))];
         end
     end
     
    ind1 = find( MCDS.discrete_cells.metadata.type == 4); %vein cells
    ind0(tcount-1) = length(MCDS.discrete_cells.live_cells)-length(ind1);
    dead(tcount-1) = length(MCDS.discrete_cells.dead_cells);
     
 end
infiltration_us_125 = infiltration;
infiltration_us_notstroma_125 = infiltration_notstroma;

livecells_us_125 = ind0;
deadcells_us_125 = dead;

%%

save('dense_us_sims.mat','infiltration_us_125','infiltration_us_notstroma_125','livecells_us_125','deadcells_us_125','-append')
save('dense_us_sims.mat','infiltration_us_15','infiltration_us_notstroma_15','livecells_us_15','deadcells_us_15','-append')
save('dense_us_sims.mat','infiltration_us_2','infiltration_us_notstroma_2','livecells_us_2','deadcells_us_2','-append')

%%

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

%%
figure
hold on 
bar([infiltration_1(1,:);infiltration_2(1,:);infiltration_3(1,:);infiltration_4(1,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('24 hours')
figure
hold on 
bar([infiltration_1(2,:);infiltration_2(2,:);infiltration_3(2,:);infiltration_4(2,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('48 hours')
figure
hold on 
bar([infiltration_1(3,:);infiltration_2(3,:);infiltration_3(3,:);infiltration_4(3,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours')


figure
hold on 
bar([infiltration_1notstroma(1,:);infiltration_2notstroma(1,:);infiltration_3notstroma(1,:);infiltration_4notstroma(1,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('24 hours')
figure
hold on 
bar([infiltration_1notstroma(2,:);infiltration_2notstroma(2,:);infiltration_3notstroma(2,:);infiltration_4notstroma(2,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('48 hours')
figure
hold on 
bar([infiltration_1notstroma(3,:);infiltration_2notstroma(3,:);infiltration_3notstroma(3,:);infiltration_4notstroma(3,:)]')
ylabel('Total intracellular virions')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
title('72 hours')

%%
