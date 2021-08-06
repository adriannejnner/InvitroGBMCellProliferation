timetotal = 37;
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
locs_GBM = find(MCDS.discrete_cells.metadata.type==2); 
locs_GBM_alive = intersect(MCDS.discrete_cells.live_cells,locs_GBM);
locs_infected_GBM_alive = intersect(locs_infected,locs_GBM_alive);
locs_infected_GBM_alive_and_dead = intersect(locs_infected,locs_GBM);
    
locs_CTL = find(MCDS.discrete_cells.metadata.type==3); 
locs_TH = find(MCDS.discrete_cells.metadata.type==1); 

locs_uninfected = find( P <= 1); %finds location of cells with less than 1 virus inside
locs_GBM_alive_uninfected = intersect(locs_GBM_alive,locs_uninfected);

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
            position_infected_notstroma =  MCDS.discrete_cells.state.position(locs_infected_GBM_alive_and_dead,1:2);
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

            infiltration(tcount/(12),:) = [sum(P(locs_infected_GBM_alive_and_dead(band1_loc))), sum(P(locs_infected_GBM_alive_and_dead(band2_loc))), sum(P(locs_infected_GBM_alive_and_dead(band3_loc))),...
                sum(P(locs_infected_GBM_alive_and_dead(band4_loc))), sum(P(locs_infected_GBM_alive_and_dead(band5_loc))), sum(P(locs_infected_GBM_alive_and_dead(band6_loc))), sum(P(locs_infected_GBM_alive_and_dead(band7_loc))),...
                sum(P(locs_infected_GBM_alive_and_dead(band8_loc))), sum(P(locs_infected_GBM_alive_and_dead(band9_loc))), sum(P(locs_infected_GBM_alive_and_dead(band10_loc))), sum(P(locs_infected_GBM_alive_and_dead(band11_loc))),...
                sum(P(locs_infected_GBM_alive_and_dead(band12_loc))), sum(P(locs_infected_GBM_alive_and_dead(band13_loc))), sum(P(locs_infected_GBM_alive_and_dead(band14_loc))), sum(P(locs_infected_GBM_alive_and_dead(band15_loc))), sum(P(locs_infected_GBM_alive_and_dead(band16_loc)))];
         end
     end
     
    live_cells(tcount-1) = length(locs_GBM_alive); % number live cells that are also GBM cells (this includes infected cells)
    infected_cells(tcount-1) = length(locs_infected_GBM_alive); % number of cells that are alive and infected >1 virus and GBM cells
    dead_cells(tcount-1) = length(MCDS.discrete_cells.dead_cells); % number of dead cells
    uninfected_cells(tcount-1) = length(locs_GBM_alive_uninfected); 
    
    CTL_number(tcount-1) = length(locs_CTL);
    TH_number(tcount-1) = length(locs_TH);
    
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    deltax = abs(MCDS.mesh.X(1,1,k)-MCDS.mesh.X(1,2,k));
    deltay = abs(MCDS.mesh.Y(1,1,k)-MCDS.mesh.Y(2,1,k));
    virion(tcount-1) = sum(sum(MCDS.continuum_variables(2).data(:,:,k)))*deltax*deltay*20;%virion
    chemokine(tcount-1) = sum(sum(MCDS.continuum_variables(4).data(:,:,k)))*deltax*deltay*20;%assembled virion
   
 end

save('dense_immune_tenfold_killers_v2.mat','infiltration','live_cells','infected_cells','dead_cells','virion','chemokine','uninfected_cells','CTL_number','TH_number')





