% initial number of GBM cells in patchy tumours

A = 'output0000000';
B = '.xml';
tcount = 2;
K = [A num2str(tcount-1,'%d') B];

MCDS = read_MultiCellDS_xml(K);
locs_GBM = find(MCDS.discrete_cells.metadata.type==2); 
locs_GBM_alive = intersect(MCDS.discrete_cells.live_cells,locs_GBM);
length(locs_GBM_alive)