
timetotal = 72;
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
    uninfected_live(tcount) = length(MCDS.discrete_cells.live_cells)-length(ind1);
    ind0(tcount) = length(MCDS.discrete_cells.live_cells)-length(ind1);
    dead(tcount) = length(MCDS.discrete_cells.dead_cells);
    cellviab1(tcount) = ind0(tcount)/(ind0(tcount)+dead(tcount));
    
    if isempty(intersect(MCDS.discrete_cells.live_cells,find(MCDS.discrete_cells.custom.intracellular_virus_amount>0)))==1
        infected(tcount) = 0;
    else
        infected(tcount) = length(intersect(MCDS.discrete_cells.live_cells,find(MCDS.discrete_cells.custom.intracellular_virus_amount>0)));
    end
    
    extracellular_virus(tcount) = sum(sum(MCDS.continuum_variables(2).data(:,:,k)))*20*20*20;
    
end
 
save('Dense_OV_sim1.mat','uninfected_live','dead','infected','extracellular_virus')

figure
hold on 
yyaxis left
plot(uninfected_live)
plot(dead)
plot(infected)
yyaxis right
plot(extracellular_virus)