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
import subprocess
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
    val = str(random()*0.00073549)
    
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
          
    # update counter      
    count = count+1 

print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")
