A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

virus_index = 2;
chemokine_index = 4;

tcount = 13;

K = [A2 num2str(tcount-1,'%d') B];
MCDS = read_MultiCellDS_xml(K);
k = find( MCDS.mesh.Z_coordinates == 0 ); 

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(virus_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 

title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(virus_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19)

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(chemokine_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 


title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(chemokine_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19)


locs_CTL = find(MCDS.discrete_cells.metadata.type==3); 
CTL_positions = MCDS.discrete_cells.state.position(locs_CTL,1:2);

locs_TH = find(MCDS.discrete_cells.metadata.type==1); 
TH_positions = MCDS.discrete_cells.state.position(locs_TH,1:2);

distance_to_centre = sqrt(MCDS.discrete_cells.state.position(:,1).^2+MCDS.discrete_cells.state.position(:,2).^2);
overall_tumour_radius = max(distance_to_centre);
figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(CTL_positions(:,1),CTL_positions(:,2),50,'filled','MarkerEdgeColor',[0 0 139]/255,'MarkerFaceColor',[127 255 212]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500]) 
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')

figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(TH_positions(:,1),TH_positions(:,2),'filled','MarkerEdgeColor',[139 0 0]/255,'MarkerFaceColor',[255,165,0]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500])
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')

%%
    
tcount = 25;
K = [A2 num2str(tcount-1,'%d') B];

MCDS = read_MultiCellDS_xml(K);
k = find( MCDS.mesh.Z_coordinates == 0 ); 

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(virus_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 


title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(virus_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19)  

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(chemokine_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 


title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(chemokine_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19) 


locs_CTL = find(MCDS.discrete_cells.metadata.type==3); 
CTL_positions = MCDS.discrete_cells.state.position(locs_CTL,1:2);

locs_TH = find(MCDS.discrete_cells.metadata.type==1); 
TH_positions = MCDS.discrete_cells.state.position(locs_TH,1:2);

distance_to_centre = sqrt(MCDS.discrete_cells.state.position(:,1).^2+MCDS.discrete_cells.state.position(:,2).^2);
overall_tumour_radius = max(distance_to_centre);
figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(CTL_positions(:,1),CTL_positions(:,2),50,'filled','MarkerEdgeColor',[0 0 139]/255,'MarkerFaceColor',[127 255 212]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500]) 
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')

figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(TH_positions(:,1),TH_positions(:,2),'filled','MarkerEdgeColor',[139 0 0]/255,'MarkerFaceColor',[255,165,0]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500])
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')
%%

tcount = 36;
K = [A2 num2str(tcount-1,'%d') B];

MCDS = read_MultiCellDS_xml(K);
k = find( MCDS.mesh.Z_coordinates == 0 ); 

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(virus_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 


title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(virus_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19)

figure
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(chemokine_index).data(:,:,k) , 20 ) ;
axis image;
colorbar;
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 


title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(chemokine_index).name , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units) );  
set(gca,'FontSize',19)

locs_CTL = find(MCDS.discrete_cells.metadata.type==3); 
CTL_positions = MCDS.discrete_cells.state.position(locs_CTL,1:2);

locs_TH = find(MCDS.discrete_cells.metadata.type==1); 
TH_positions = MCDS.discrete_cells.state.position(locs_TH,1:2);

distance_to_centre = sqrt(MCDS.discrete_cells.state.position(:,1).^2+MCDS.discrete_cells.state.position(:,2).^2);
overall_tumour_radius = max(distance_to_centre);
figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(CTL_positions(:,1),CTL_positions(:,2),50,'filled','MarkerEdgeColor',[0 0 139]/255,'MarkerFaceColor',[127 255 212]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500]) 
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')

figure
hold on 
fill(overall_tumour_radius.*(cos(0:0.01:2*pi)),overall_tumour_radius.*(sin(0:0.01:2*pi)),[0.95 0.95 0.95],'EdgeColor','none')
scatter(TH_positions(:,1),TH_positions(:,2),'filled','MarkerEdgeColor',[139 0 0]/255,'MarkerFaceColor',[255,165,0]/255,'LineWidth',1.5)
ylim([-1500 1500])
xlim([-1500 1500])
pbaspect([1 1 1])
set(gca,'FontSize',19)
ylabel('y (micron)')
xlabel('x (micron)')
