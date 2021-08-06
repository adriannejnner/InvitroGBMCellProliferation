timetotal =168;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';
 
 tcount = 1;
 K = [A num2str(tcount-1,'%d') B];
 MCDS = read_MultiCellDS_xml(K); 
 P = MCDS.discrete_cells.state.position;
 locs_TH = find( MCDS.discrete_cells.metadata.type == 1); %vein cells
 Initial_TH = length(locs_TH);
 locs_CTL = find( MCDS.discrete_cells.metadata.type == 3); %vein cells
 Initial_CTL = length(locs_CTL);
 locs_cancer = find( MCDS.discrete_cells.metadata.type == 2); %vein cells
 Initial_cancer = length(locs_cancer);
 
 Initial_dead = length(MCDS.discrete_cells.dead_cells);
 

 for tcount = 3:timetotal
   
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    
 MCDS = read_MultiCellDS_xml(K); 
 P = MCDS.discrete_cells.state.position;
 locs_TH = find( MCDS.discrete_cells.metadata.type == 1); %vein cells
 fold_TH(tcount-2) = length(locs_TH)/Initial_TH;
 locs_CTL = find( MCDS.discrete_cells.metadata.type == 3); %vein cells
 fold_CTL(tcount-2)  = length(locs_CTL)/Initial_CTL;
 locs_cancer = find( MCDS.discrete_cells.metadata.type == 2); %vein cells
 fold_cancer(tcount-2)  = length(locs_cancer)/Initial_cancer;
 
 fold_dead(tcount-2)  = length(MCDS.discrete_cells.dead_cells);%/Initial_dead;

 deadcells(tcount-2) = length(MCDS.discrete_cells.dead_cells);
    
 end
 
figure
hold on
plot(fold_cancer*Initial_cancer,'LineWidth',1)

x = [0 1 2 3];
SE = [0 0.75*1e6 0.4*1e6 0.6*1e6];
y = [1.2*1e6, 2.4*1e6, 3.48*1e6, 4.28*1e6];

errorbar(x*24*2,y/1000,SE/1000,'o','LineWidth',2)
set(gca,'FontSize',15)
xlabel('Time (hours)')
ylabel('Fold increase')
legend('Cancer','CTL','TH','Dead cells')
set(gca,'Xtick',linspace(0,144,4),'Xticklabel',{'0','24','48','72'})



figure
hold on
yyaxis left
plot(fold_CTL*Initial_CTL,'LineWidth',1)
ylabel('No. of CTLs')
yyaxis right
plot(fold_TH*Initial_TH,'LineWidth',1)
ylabel('No. of TH cells')
set(gca,'FontSize',15)
legend('CTL','TH')
set(gca,'Xtick',linspace(0,144,4),'Xticklabel',{'0','24','48','72'})
xlabel('Time (hours)')

figure
hold on 
plot(deadcells)
ylabel('No. of dead cells')
set(gca,'Xtick',linspace(0,144,4),'Xticklabel',{'0','24','48','72'})
xlabel('Time (hours)')
set(gca,'FontSize',15)

STOP
%%
%sim1_v0_less.mat
sim2_v0_less.mat
sim3_v0_less.mat
%sim3_v0_more.mat % don't think this is actually more...
sim3_v0_more_linear.mat
sim1_v0_more.mat
sim2_v0_more.mat
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


%% 
 clear all
 load('tri_Dv_200per.mat')
 load('tri_Dv_third.mat')
 load('tri_beta_10times.mat')
 load('tri_cT_01times.mat')
 load('tri_TRAIL_st01_stau500.mat')
 load('tri_cT_2.mat')
 load('tri_Dv_10.mat')

 
 timetotal =96;
 figure
 hold on
 plot([1:timetotal]*60,noofTcellssT01stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_Dv_10,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_Dv_200per,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_Dv__third,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_beta_10times,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_cT_01times,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,tri_cT_2,'o-','LineWidth',1,'MarkerSize',3)
 legend('Original simulation','D_V \times 10','D_V \times 3','D_V\times 1/3','c_I \times 10','c_R \times 0.1','c_R \times 2')
  set(gca,'FontSize',13)
 xlabel('Time (minutes)')
 ylabel('No. of tumour cells')
%% TUMOUR ST STAU PERT IMAGE TRI
 clear all
 
 load('tri_TRAIL_st00001_stau500.mat')
 load('tri_TRAIL_st001_stau50.mat')
 load('tri_TRAIL_st001_stau500.mat')
 load('tri_TRAIL_st01_stau500.mat')
 load('tri_TRAIL_st01_stau10.mat')
 load('tri_TRAIL_st01_stau50.mat')
 load('tri_TRAIL_st01_stau100.mat')
 
 timetotal =96;
 figure
 hold on
 plot([1:timetotal]*60,noofTcellssT00001stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT001stau50,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT001stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau500,'o-','LineWidth',1,'MarkerSize',3)
% plot([1:timetotal]*60,noofTcellssT01stau10,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau50,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau100,'o-','LineWidth',1,'MarkerSize',3)
 set(gca,'FontSize',13)
 xlabel('Time (minutes)')
 ylabel('No. of tumour cells')
 legend( 's_T = 0.0001, s_\tau = 500',...
     's_T = 0.01, s_\tau = 50',...
     's_T = 0.01, s_\tau = 500',...
     's_T = 0.1, s_\tau = 500',... %'s_T = 0.1, s_\tau = 10',...
     's_T = 0.1, s_\tau = 50',...
     's_T = 0.1, s_\tau = 100')


 %% TUMOUR ST STAU PERT IMAGE CIRC
 clear all
 
 load('circ_TRAIL_st00001_stau500.mat')
 load('circ_TRAIL_st00001_stau100.mat') 
 load('circ_TRAIL_st001_stau50.mat') 
 load('circ_TRAIL_st01_stau10.mat') 
 load('circ_TRAIL_st01_stau50.mat') 
 load('circ_TRAIL_st01_stau500.mat') 
 load('circ_TRAIL_st001_stau500.mat') 
 
 timetotal =96;
 figure
 hold on
% plot([1:timetotal]*60,noofTcellssT00001stau100,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT00001stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT001stau50,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT001stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau50,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellssT01stau10,'o-','LineWidth',1,'MarkerSize',3)
 set(gca,'FontSize',13)
 xlabel('Time (minutes)')
 ylabel('No. of tumour cells')
 %'s_T = 0.0001, s_\tau = 100',...
 legend( 's_T = 0.0001, s_\tau = 500',...
     's_T = 0.01, s_\tau = 50',...
     's_T = 0.01, s_\tau = 500',...
     's_T = 0.1, s_\tau = 500',...
     's_T = 0.1, s_\tau = 50',...
     's_T = 0.1, s_\tau = 10')
 
 %% VIRUS DOSAGE IMAGE
 
 load('sT001stau500')
 load('VIRUSday012')
 
 figure
 hold on
 plot([1:timetotal]*60,noofTcellssT001stau500,'o-','LineWidth',1,'MarkerSize',3)
 plot([1:timetotal]*60,noofTcellsVIRUSday012,'o-','LineWidth',1,'MarkerSize',3)
 set(gca,'FontSize',13)
 xlabel('Time (minutes)')
 ylabel('No. of tumour cells')
 
 legend('Single injection day 0', 'Three injections day 0, 1 and 2')
 
 
 