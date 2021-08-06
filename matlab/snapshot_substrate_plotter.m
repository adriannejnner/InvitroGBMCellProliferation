timetotal =13;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000'; 
A4 = 'output0000'; 
B = '.xml';

var = 2;

figure
hold on 

tcount = timetotal;
clf
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
k = find( MCDS.mesh.Z_coordinates == 0 ); 
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
MCDS.continuum_variables(var).data(:,:,k) , 20 ) ;
axis image;
colorbar; 
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 

title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(var).name , ...
 MCDS.continuum_variables(var).units , ...
 MCDS.metadata.current_time , ...
 MCDS.metadata.time_units, ...
 MCDS.mesh.Z_coordinates(k), ...
 MCDS.metadata.spatial_units ) ); 