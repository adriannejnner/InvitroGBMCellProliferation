timetotal = 37;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

 for tcount = 1:timetotal
    clf
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    Xmesh = MCDS.mesh.X(:,:,k);
    Ymesh = MCDS.mesh.Y(:,:,k);
    Datamesh = MCDS.continuum_variables(2).data(:,:,k);

    Dens1(tcount) = sum(sum(Datamesh))*abs(Xmesh(1,1)-Xmesh(1,2))*abs(Ymesh(1,1)-Ymesh(2,1));
    ind1 = find( MCDS.discrete_cells.metadata.type == 4);
    ind0(tcount) = length(MCDS.discrete_cells.live_cells)-length(ind1);
    dead(tcount) = length(MCDS.discrete_cells.dead_cells);
    cellviab1(tcount) = ind0(tcount)/(ind0(tcount)+dead(tcount));
    
    
 end
 
 live_cells = ind0;
save('contol_immune_tumour_antigen.mat','live_cells')


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
% 
% figure
% hold on 
% plot(cellviab1,'LineWidth',2)
% xlabel('Hours')
% ylabel('Cell Viability %')
% set(gca,'FontSize',15,'xtick',[0 24 48 72])
% xlim([0 72])
% 
% figure
% hold on 
% plot(Dens1,'LineWidth',2)
% xlabel('Hours')
% ylabel('Secreted stTRAIL (\mu g/mL) %')
% set(gca,'FontSize',15,'xtick',[0 24 48 72])
% xlim([0 72])

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