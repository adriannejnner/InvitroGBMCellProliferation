# params_run.py:
#
# This Python script provides simple parameter exploration functionality. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default 
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder, then
# runs the simulation (optionally, in the background) which writes simulation results into the new folder.
# 
# Author: Randy Heiland

import xml.etree.ElementTree as ET
from shutil import copyfile
import os
import sys
import glob
import subprocess
import numpy as np
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt
from random import seed
from random import random

# print(len(sys.argv))
exec_pgm = 'GBM_OV_immune_stroma'#sys.argv[1]

sequential_flag = 1  # if =1, do runs sequentially, i.e., not in background
#omp_num_threads 4
#initial_GBM_cells 12556
#initial_stroma_cells 2861 
#folder TESTINGPARAMRUN
#virus_burst_number 660

xml_file_in = 'config/PhysiCell_settings.xml'
xml_file_out = 'config/tmp.xml'
copyfile(xml_file_in, xml_file_out)
tree = ET.parse(xml_file_out)
xml_root = tree.getroot()
first_time = True
output_dirs = []

# load data points to compare to number of cells
data_vec = np.array([347, 618, 799, 861, 1007]) # simulating only 0.0069 of the dish
     
count = 0
while count < 50:

    #set up folder to store this run
    key = 'folder'
    val = 'Results_'+str(count)

    folder_name = val
    output_dirs.append(folder_name)
    if (not os.path.exists(folder_name)):
        print("--- parsed 'folder', makedir " + folder_name)
        os.makedirs(folder_name)
    xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

    try:
        xml_root.find('.//' + key).text = val
    except:
        print("--- Error: could not find ",key," in .xml\n")
        sys.exit(1)    
    
    # set parameter values to simulate
    key = 'GBM_cell_proliferation_rate'
    val = str(random()*0.000573549)
        
    try:
        xml_root.find('.//' + key).text = val
    except:
        print("--- Error: could not find ",key," in .xml\n")
        sys.exit(1)  
    
    key = 'run_it'
    val = 'dummy'
    
    # write the new config file with the new parameter value and simulate
    print('\n\n---> write config file (and start sim): ', xml_file_out)
    tree.write(xml_file_out)   # will create folder_name/config.xml
    log_file = folder_name + ".log"  
    if sequential_flag > 0:
        # cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
        cmd =  exec_pgm + " " + xml_file_out + " > " + log_file 
        print("Doing sequential run...")
        print("cmd = ",cmd)
        os.system(cmd)   # <------ Execute the simulation
    else:  # put (multiple) runs in the background
        with open(log_file,"w") as outf:
            subprocess.Popen([exec_pgm, xml_file_out],stdout=outf)
          
    # output number of cells       
    data_dir = 'Results_'+str(count)
    print('# data_dir = ',data_dir)
    os.chdir(data_dir)
    xml_files = glob.glob('output*.xml')
    os.chdir('..')
    xml_files.sort()
    
    ds_count = len(xml_files)
    print("# ----- ds_count = ",ds_count)
    mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]

    tval = np.linspace(0, mcds[-1].get_time(), ds_count)
    print('type(tval)= ',type(tval))
    print('tval= ',repr(tval))
    
    # finds cells that are cell type 2, i.e. GBM cell type, and still cycling (so not dead)
    yval4 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 2) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )
    print('Number GBM cells = ',repr(yval4))       
    
    # calculate vector of residuals
    residual = np.subtract(data_vec,yval4)
    print("Residual = ", residual)
              
    # update counter      
    count = count+1 

print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")
