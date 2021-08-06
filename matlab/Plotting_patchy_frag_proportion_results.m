
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

            infiltration_prop_70_30 = [sum(P(band1_loc)), sum(P(band2_loc)), sum(P(band3_loc)),...
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

            infiltration_notstroma_prop_70_30 = [sum(P(interr(band1_loc))), sum(P(interr(band2_loc))), sum(P(interr(band3_loc))),...
                sum(P(interr(band4_loc))), sum(P(interr(band5_loc))), sum(P(interr(band6_loc))), sum(P(interr(band7_loc))),...
                sum(P(interr(band8_loc))), sum(P(interr(band9_loc))), sum(P(interr(band10_loc))), sum(P(interr(band11_loc))),...
                sum(P(interr(band12_loc))), sum(P(interr(band13_loc))), sum(P(interr(band14_loc))), sum(P(interr(band15_loc))), sum(P(interr(band16_loc)))];
     end
    
    live_cells_prop_70_30(tcount-1) = length(MCDS.discrete_cells.live_cells); % number live cells
    infected_cells_prop_70_30(tcount-1) = length(intersect(find(MCDS.discrete_cells.custom.intracellular_virus_amount >1),MCDS.discrete_cells.live_cells)); %infected cell number
    dead_cells_prop_70_30(tcount-1) = length(MCDS.discrete_cells.dead_cells); % number dead cells

end

save('prop_70_30','infiltration_prop_70_30','infiltration_notstroma_prop_70_30',...
    'live_cells_prop_70_30','infected_cells_prop_70_30','dead_cells_prop_70_30')

%%
load('prop_0_100.mat')
load('prop_10_90.mat')
load('prop_20_80.mat')
load('prop_30_70.mat')
load('prop_40_60.mat')
load('prop_50_50.mat')
load('prop_60_40.mat')
load('prop_70_30.mat')
load('prop_80_20.mat')
load('prop_90_10.mat')
load('prop_100_0.mat')

figure
hold on 
plot(live_cells_prop_0_100/(live_cells_prop_0_100(1))*100,'LineWidth',2)
plot(live_cells_prop_50_50/(live_cells_prop_50_50(1))*100,'LineWidth',2)
plot(live_cells_prop_100_0/(live_cells_prop_100_0(1))*100,'LineWidth',2)
legend('0:100','50:50','100:0')
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'XTick',linspace(1,37,5),'XTickLabels',{'0','18','36','54','72'})
set(gca,'FontSize',16)

figure
hold on 
bar([live_cells_prop_0_100(end)/(live_cells_prop_0_100(1))*100;...
    live_cells_prop_10_90(end)/(live_cells_prop_10_90(1))*100;...
    live_cells_prop_20_80(end)/(live_cells_prop_20_80(1))*100;...
    live_cells_prop_30_70(end)/(live_cells_prop_30_70(1))*100;...
    live_cells_prop_40_60(end)/(live_cells_prop_40_60(1))*100;...
    live_cells_prop_50_50(end)/(live_cells_prop_50_50(1))*100;...
    live_cells_prop_60_40(end)/(live_cells_prop_60_40(1))*100;...
    live_cells_prop_70_30(end)/(live_cells_prop_70_30(1))*100;...
    live_cells_prop_80_20(end)/(live_cells_prop_80_20(1))*100;...
    live_cells_prop_90_10(end)/(live_cells_prop_90_10(1))*100;...
    live_cells_prop_100_0(end)/(live_cells_prop_100_0(1))*100])
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'FontSize',16)
set(gca,'XTick',linspace(1,11,11),'XTickLabels',{'0:100','10:90','20:80','30:70','40:60','50:50','60:40','70:30','80:20','90:10','100:0'})


figure
hold on 
plot([live_cells_prop_0_100(end);...
    live_cells_prop_10_90(end);...
    live_cells_prop_20_80(end);...
    live_cells_prop_30_70(end);...
    live_cells_prop_40_60(end);...
    live_cells_prop_50_50(end);...
    live_cells_prop_60_40(end);...
    live_cells_prop_70_30(end);...
    live_cells_prop_80_20(end);...
    live_cells_prop_90_10(end);...
    live_cells_prop_100_0(end)],'o-','LineWidth',2)
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'FontSize',16)
set(gca,'XTick',linspace(1,11,11),'XTickLabels',{'0:100','10:90','20:80','30:70','40:60','50:50','60:40','70:30','80:20','90:10','100:0'})



figure
hold on 
plot([live_cells_prop_0_100(end)/(live_cells_prop_0_100(1))*100;...
    live_cells_prop_10_90(end)/(live_cells_prop_10_90(1))*100;...
    live_cells_prop_20_80(end)/(live_cells_prop_20_80(1))*100;...
    live_cells_prop_30_70(end)/(live_cells_prop_30_70(1))*100;...
    live_cells_prop_40_60(end)/(live_cells_prop_40_60(1))*100;...
    live_cells_prop_50_50(end)/(live_cells_prop_50_50(1))*100;...
    live_cells_prop_60_40(end)/(live_cells_prop_60_40(1))*100;...
    live_cells_prop_70_30(end)/(live_cells_prop_70_30(1))*100;...
    live_cells_prop_80_20(end)/(live_cells_prop_80_20(1))*100;...
    live_cells_prop_90_10(end)/(live_cells_prop_90_10(1))*100;...
    live_cells_prop_100_0(end)/(live_cells_prop_100_0(1))*100],'o-','LineWidth',2)
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'FontSize',16)
set(gca,'XTick',linspace(1,11,11),'XTickLabels',{'0:100','10:90','20:80','30:70','40:60','50:50','60:40','70:30','80:20','90:10','100:0'})

figure
hold on 
bar([live_cells_prop_0_100(end);...
    live_cells_prop_10_90(end);...
    live_cells_prop_20_80(end);...
    live_cells_prop_30_70(end);...
    live_cells_prop_40_60(end);...
    live_cells_prop_50_50(end);...
    live_cells_prop_60_40(end);...
    live_cells_prop_70_30(end);...
    live_cells_prop_80_20(end);...
    live_cells_prop_90_10(end);...
    live_cells_prop_100_0(end)])
ylabel('% initial tumour remaining')
xlabel('Time (hours)')
set(gca,'FontSize',16)
set(gca,'XTick',linspace(1,11,11),'XTickLabels',{'0:100','10:90','20:80','30:70','40:60','50:50','60:40','70:30','80:20','90:10','100:0'})

figure
hold on
bar([infiltration_notstroma_prop_0_100;...
   infiltration_notstroma_prop_10_90;...
   infiltration_notstroma_prop_20_80;...
   infiltration_notstroma_prop_30_70;...
   infiltration_notstroma_prop_40_60;...
   infiltration_notstroma_prop_50_50]')
legend('0:100','10:90','20:80','30:70','40:60','50:50')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)


figure
hold on
bar([infiltration_notstroma_prop_60_40;...
   infiltration_notstroma_prop_70_30;...
   infiltration_notstroma_prop_80_20;...
   infiltration_notstroma_prop_90_10;...
   infiltration_notstroma_prop_100_0]')
legend('60:40','70:30','80:20','90:10','100:0')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)


figure
hold on
bar([infiltration_notstroma_prop_0_100;...
   infiltration_notstroma_prop_50_50;...
   infiltration_notstroma_prop_100_0]')
legend('0:100','50:50','100:0')
ylabel('Total intracellular virions')
xlabel('\mu from periphery')
set(gca,'Xtick',linspace(1,16,16),'Xticklabel',{'0','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750'})
set(gca,'FontSize',16)

%%


figure
hold on 
b = bar([live_cells_prop_0_100(end)/(live_cells_prop_0_100(1))*100;...
    live_cells_prop_10_90(end)/(live_cells_prop_10_90(1))*100;...
    live_cells_prop_20_80(end)/(live_cells_prop_20_80(1))*100;...
    live_cells_prop_30_70(end)/(live_cells_prop_30_70(1))*100;...
    live_cells_prop_40_60(end)/(live_cells_prop_40_60(1))*100;...
    live_cells_prop_50_50(end)/(live_cells_prop_50_50(1))*100;...
    live_cells_prop_60_40(end)/(live_cells_prop_60_40(1))*100;...
    live_cells_prop_70_30(end)/(live_cells_prop_70_30(1))*100;...
    live_cells_prop_80_20(end)/(live_cells_prop_80_20(1))*100;...
    live_cells_prop_90_10(end)/(live_cells_prop_90_10(1))*100;...
    live_cells_prop_100_0(end)/(live_cells_prop_100_0(1))*100],'FaceColor','flat','EdgeColor',[0.5 0.5 0.5])
ylabel('% Fragment remaining')
xlabel('Time (hours)')
set(gca,'FontSize',16)
set(gca,'XTick',linspace(1,11,11),'XTickLabels',{'0:100','10:90','20:80','30:70','40:60','50:50','60:40','70:30','80:20','90:10','100:0'})

ylim([60 100])

initial_GBM_prop_0_100 = 5602;
initial_GBM_prop_10_90 = 6311;
initial_GBM_prop_20_80 = 7084;
initial_GBM_prop_30_70 = 7927;
initial_GBM_prop_40_60 = 8836;
initial_GBM_prop_50_50 = 9850;
initial_GBM_prop_60_40 = 9825;
initial_GBM_prop_70_30 = 10518;
initial_GBM_prop_80_20 = 11221;
initial_GBM_prop_90_10 = 11929;
initial_GBM_prop_100_0 = 12609;

minGBMcells = min([initial_GBM_prop_0_100(1),initial_GBM_prop_10_90(1),initial_GBM_prop_20_80(1),initial_GBM_prop_30_70,initial_GBM_prop_40_60(1),initial_GBM_prop_50_50(1),initial_GBM_prop_60_40(1),initial_GBM_prop_70_30(1),initial_GBM_prop_80_20(1),initial_GBM_prop_90_10(1),initial_GBM_prop_100_0(1)]);
maxGBMcells = max([initial_GBM_prop_0_100(1),initial_GBM_prop_10_90(1),initial_GBM_prop_20_80(1),initial_GBM_prop_30_70,initial_GBM_prop_40_60(1),initial_GBM_prop_50_50(1),initial_GBM_prop_60_40(1),initial_GBM_prop_70_30(1),initial_GBM_prop_80_20(1),initial_GBM_prop_90_10(1),initial_GBM_prop_100_0(1)]);

colvec = linspace(minGBMcells,maxGBMcells+1,10);

colmap = [255,247,251
236,226,240
208,209,230
166,189,219
103,169,207
54,144,192
2,129,138
1,108,89
1,70,54]/255;


b.CData(1,:) = colmap(find(colvec>initial_GBM_prop_0_100,1)-1,:);
b.CData(2,:) = colmap(find(colvec>initial_GBM_prop_10_90,1)-1,:)
b.CData(3,:) = colmap(find(colvec>initial_GBM_prop_20_80,1)-1,:)
b.CData(4,:) = colmap(find(colvec>initial_GBM_prop_30_70,1)-1,:)
b.CData(5,:) = colmap(find(colvec>initial_GBM_prop_40_60,1)-1,:)
b.CData(6,:) = colmap(find(colvec>initial_GBM_prop_50_50,1)-1,:)
b.CData(7,:) = colmap(find(colvec>initial_GBM_prop_60_40,1)-1,:)
b.CData(8,:) = colmap(find(colvec>initial_GBM_prop_70_30,1)-1,:)
b.CData(9,:) = colmap(find(colvec>initial_GBM_prop_80_20,1)-1,:)
b.CData(10,:) = colmap(find(colvec>initial_GBM_prop_90_10,1)-1,:)
b.CData(11,:) = colmap(find(colvec>initial_GBM_prop_100_0,1)-1,:)

colormap(colmap)
c = colorbar
set(c,'Ticks',linspace(0,1,6),'TickLabels',{'0.56\times 10^4','0.4\times 10^4','0.84\times 10^4','0.98\times 10^4','1.1\times 10^4','1.3\times 10^4'})

