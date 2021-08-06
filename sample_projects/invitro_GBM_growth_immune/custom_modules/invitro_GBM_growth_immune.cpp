/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./invitro_GBM_growth_immune.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>


Cell_Definition TH_cell; 
Cell_Definition CTL_cell; 
Cell_Definition cancer_cell; 

void create_TH_cells( void )
{
	TH_cell = cell_defaults;
	
	TH_cell.name = "TH cell";
	TH_cell.type = 1;
	
	// proliferation 
	TH_cell.functions.cycle_model = Ki67_basic;
	TH_cell.phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 

	TH_cell.phenotype.secretion.uptake_rates[1] = 0.0;
	TH_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	TH_cell.phenotype.motility.migration_speed = 4;
	TH_cell.phenotype.geometry.radius = 3.6;
	
	TH_cell.phenotype.volume.total = 185.66;
	TH_cell.phenotype.volume.fluid_fraction = 0.75;
	TH_cell.phenotype.volume.fluid = TH_cell.phenotype.volume.fluid_fraction*TH_cell.phenotype.volume.total;
	TH_cell.phenotype.volume.solid = TH_cell.phenotype.volume.total-TH_cell.phenotype.volume.fluid;
	TH_cell.phenotype.volume.nuclear = 95.21;
	TH_cell.phenotype.volume.nuclear_solid = 23.8;
	TH_cell.phenotype.volume.nuclear_fluid = TH_cell.phenotype.volume.nuclear - TH_cell.phenotype.volume.nuclear_solid;
	TH_cell.phenotype.volume.cytoplasmic = TH_cell.phenotype.volume.total - TH_cell.phenotype.volume.nuclear;
	TH_cell.phenotype.volume.cytoplasmic_fluid = TH_cell.phenotype.volume.fluid_fraction*TH_cell.phenotype.volume.cytoplasmic;
	TH_cell.phenotype.volume.cytoplasmic_solid = TH_cell.phenotype.volume.cytoplasmic-TH_cell.phenotype.volume.cytoplasmic_fluid;
	TH_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.05;
	TH_cell.phenotype.volume.target_solid_cytoplasmic = TH_cell.phenotype.volume.cytoplasmic_solid;
	TH_cell.phenotype.volume.target_solid_nuclear = TH_cell.phenotype.volume.nuclear_solid;
	TH_cell.phenotype.volume.target_fluid_fraction = TH_cell.phenotype.volume.fluid_fraction;
	TH_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = TH_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	//TH_cell.phenotype.mechanics.set_relative_equilibrium_distance( 24/10.75 );
	
	TH_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.9026*1e-5; 	
	TH_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = 0.00143; 	
	
	return;
}
void create_CTL_cells( void )
{
	CTL_cell = cell_defaults;
	
	CTL_cell.name = "CTL cell";
	CTL_cell.type = 3;
	
	// turn off proliferation 
	CTL_cell.functions.cycle_model = Ki67_basic;
	CTL_cell.phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive); 

	CTL_cell.phenotype.secretion.uptake_rates[1] = 0.0;
	CTL_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	CTL_cell.phenotype.motility.migration_speed = 4;
	CTL_cell.phenotype.geometry.radius = 3.6;
	
	CTL_cell.phenotype.volume.total = 185.66;
	CTL_cell.phenotype.volume.fluid_fraction = 0.75;
	CTL_cell.phenotype.volume.fluid = CTL_cell.phenotype.volume.fluid_fraction*CTL_cell.phenotype.volume.total;
	CTL_cell.phenotype.volume.solid = CTL_cell.phenotype.volume.total-CTL_cell.phenotype.volume.fluid;
	CTL_cell.phenotype.volume.nuclear = 96.23;
	CTL_cell.phenotype.volume.nuclear_solid = 24.06;
	CTL_cell.phenotype.volume.nuclear_fluid = CTL_cell.phenotype.volume.nuclear - CTL_cell.phenotype.volume.nuclear_solid;
	CTL_cell.phenotype.volume.cytoplasmic = CTL_cell.phenotype.volume.total - CTL_cell.phenotype.volume.nuclear;
	CTL_cell.phenotype.volume.cytoplasmic_fluid = CTL_cell.phenotype.volume.fluid_fraction*CTL_cell.phenotype.volume.cytoplasmic;
	CTL_cell.phenotype.volume.cytoplasmic_solid = CTL_cell.phenotype.volume.cytoplasmic-CTL_cell.phenotype.volume.cytoplasmic_fluid;
	CTL_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.03;
	CTL_cell.phenotype.volume.target_solid_cytoplasmic = CTL_cell.phenotype.volume.cytoplasmic_solid;
	CTL_cell.phenotype.volume.target_solid_nuclear = CTL_cell.phenotype.volume.nuclear_solid;
	CTL_cell.phenotype.volume.target_fluid_fraction = CTL_cell.phenotype.volume.fluid_fraction;
	CTL_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = CTL_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.2206*1e-5; 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = 0.00143; 	
	
	return;
}
void create_cancer_cells( void )
{
	cancer_cell = cell_defaults;
	
	static int virus_index = microenvironment.find_density_index( "virus");
	cancer_cell.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cancer_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	
	cancer_cell.phenotype.geometry.radius = 10.75;
		
	cancer_cell.phenotype.motility.migration_speed = 0.0;//0.05;
	
	cancer_cell.phenotype.volume.total = 4/3*3.1416*(cancer_cell.phenotype.geometry.radius)*(cancer_cell.phenotype.geometry.radius)*(cancer_cell.phenotype.geometry.radius);//5203.7;
	cancer_cell.phenotype.volume.fluid_fraction = 0.75;
	cancer_cell.phenotype.volume.fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.total;
	cancer_cell.phenotype.volume.solid = cancer_cell.phenotype.volume.total-cancer_cell.phenotype.volume.fluid;
	cancer_cell.phenotype.volume.nuclear = 740;
	cancer_cell.phenotype.volume.nuclear_solid = 185;
	cancer_cell.phenotype.volume.nuclear_fluid = cancer_cell.phenotype.volume.nuclear - cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.cytoplasmic = cancer_cell.phenotype.volume.total - cancer_cell.phenotype.volume.nuclear;
	cancer_cell.phenotype.volume.cytoplasmic_fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.cytoplasmic;
	cancer_cell.phenotype.volume.cytoplasmic_solid = cancer_cell.phenotype.volume.cytoplasmic-cancer_cell.phenotype.volume.cytoplasmic_fluid;
	cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic/cancer_cell.phenotype.volume.nuclear;//6.0321;
	cancer_cell.phenotype.volume.target_solid_cytoplasmic = cancer_cell.phenotype.volume.cytoplasmic_solid;
	cancer_cell.phenotype.volume.target_solid_nuclear = cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.target_fluid_fraction = cancer_cell.phenotype.volume.fluid_fraction;
	cancer_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	cancer_cell.phenotype.volume.calcified_fraction = 0; 
	cancer_cell.phenotype.volume.calcification_rate = 0;
	
	cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.35*cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength;
	
	cancer_cell.name = "cancer cell";
	cancer_cell.type = 2; 
	
	return;
}

void create_cell_types( void )
{
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	//cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	double C0 = parameters.doubles("C0");
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;//0.00073549;//0.0000064812*(1-9740/512720);
	
 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 
	
	// set equilbrium distance mechanics
	//std::cout<<cell_defaults.phenotype.mechanics.set_absolute_equilibrium_distance<<std::endl;
	//cell_defaults.phenotype.mechanics.set_relative_equilibrium_distance( 73.46 );
	
	// reduce cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.0;//0.05;
	
	//std::cout<<cell_defaults.phenotype.motility.migration_speed<<std::endl;
	
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "bias_phenotype", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "replication_start_time", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "infection_time_start", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "infection_time_length", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "time_of_injection", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	
	
	// turn off secretion from these cells (oxygen and virus)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	static int virus_index = microenvironment.find_density_index( "virus");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = 10; 
	cell_defaults.phenotype.molecular.fraction_released_at_death[virus_index] = 1;
	static int wall_index = microenvironment.find_density_index( "wall");
	
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
	
	// update cell and phenotype based on virus dynamics only
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3; 
	
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3;//ECM_movement;//virus_dynamics; 
	cell_defaults.functions.update_phenotype = only_proliferation;
	//cell_defaults.functions.update_phenotype = infection_dynamics;
	cell_defaults.phenotype.motility.is_motile = true; 
	
	// change the max cell-cell adhesion distance 
	//cell_defaults.phenotype.mechanics.set_relative_maximum_adhesion_distance(parameters.doubles("max_relative_cell_adhesion_distance") );
	
	cell_defaults.name = "holder cell"; 
	cell_defaults.type = 0; 
	
	//create the immune and cancer cell types
	
	create_cancer_cells();
	create_CTL_cells();
	create_TH_cells();
		
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters
	
	// make sure not override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	
	initialize_microenvironment(); 	
	
	static int virus_index = microenvironment.find_density_index( "virus");
	static int oxygen_index = microenvironment.find_density_index( "oxygen");
	static int wall_index = microenvironment.find_density_index( "wall");
	
	microenvironment.diffusion_coefficients[virus_index] = 1.4766; 
	microenvironment.decay_rates[virus_index] = 0;
	
	microenvironment.diffusion_coefficients[oxygen_index] = 0; 
	microenvironment.decay_rates[oxygen_index] = 0;
	
	microenvironment.diffusion_coefficients[wall_index] = 0; 
	microenvironment.decay_rates[wall_index] = 0;
	
	std::vector<double> ve(3);
    ve[0] = 0;
	ve[1] = 10;
	ve[2] = 0;
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		if( ECMdense[1]<-678.33/2+5) //white
		{microenvironment(n)[virus_index] = 0;}//600;}	
		else
		{microenvironment(n)[virus_index] = 0;}
		
		microenvironment(n)[oxygen_index] = 0;//-2/678.33*(ECMdense[1]+678.33/2)+10;
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>400*400)//1270*1270)
		{	
			//get_default_microenvironment()->add_dirichlet_node(n,ve);
			microenvironment(n)[virus_index] = 0;//60;
			microenvironment(n)[wall_index] = 10;
			
		}
		else
		{   microenvironment(n)[wall_index] = 0;}
		
	}
	return; 
}

void setup_tissue_circle_immune( void )
{
	
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
		
	double Radius = 400;//1270;
	
	double start_of_boundary = -Radius;
	double end_of_boundary = Radius;
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	/*
		
	//UNPROLIFERATIVE TH CELLS
	x = -18.0725;
	y = 922.2948;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 954.1439;
	y = 56.7601;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
    x = -115.9499;
	y = -349.9367;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
    x = 773.5733;
	y = -613.7738;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
    x = 753.4041;
	y = -219.3846;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = -292.5135;
	y = -91.2428;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = 449.1283;
	y = 407.6992;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 452.4130;
	y = 616.9196;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  60.4108;
	y =  -1021.8;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 534.5535;
	y =  -882.0731;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  15.5367;
	y = 400.0095;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = 412.5616;
	y = -921.4520;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  64.9108;
	y = 1021.7;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 669.1877;
	y = -318.7248;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -483.3989;
	y = 759.6808;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = 134.3949;
	y = 347.2437;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = -47.5182;
	y = 611.5124;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -738.2632;
	y = -659.5292;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -937.1317;
	y = 156.6883;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -614.4497;
	y = 789.3966;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  437.2258;
	y = -793.7856;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = -170.1507;
	y = -132.5528;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = -885.0803;
	y = -277.3879;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  838.1983;
	y = -501.9633;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -260.7010;
	y =  604.0757;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  103.6960;
	y = -899.8739;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -23.0385;
	y = -922.1671;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -458.5337;
	y = 474.7504;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -765.6844;
	y = -326.6257;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  360.4622;
	y = 228.2000;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  771.4487;
	y = 319.7186;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = -195.5453;
	y = -5.9414;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x =  114.0970;
	y = -561.2747;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x =  203.4536;
	y = -89.5478;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  235.8782;
	y = 265.3679;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = -868.0050;
	y = -405.4009;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -811.6261;
	y = 187.1180;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	//PROLIFERATIVE TH CELLS
	x = 454.1934;
	y = -46.7240;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	
	x = -427.8415;
	y = 616.7429;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	*/
	//UNPROLIFERATIVE CTL CELLS
	x = 27.9074;
	y = 119.1812;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = 225.7304;
	y = -728.5426;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = -258.8621;
	y =  265.2323;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = -809.2029;
	y =  -60.6306;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 473.7135;
	y =  887.8817;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -366.2289;
	y = -332.2614;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -197.1932;
	y =  446.5596;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -393.1067;
	y =  -597.4458;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -264.8455;
	y =  -840.0492;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -62.2700;
	y = -690.4383;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 968.1414;
	y =  -363.0428;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 592.2715;
	y =   125.8000;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 231.3138;
	y =  809.8071;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 547.2984;
	y =  -545.0700;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 215.9077;
	y =  557.6151;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	*/
	x = 141.7866;
	y = -321.7799;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	/*
	x = -557.2198;
	y = -155.5629;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 809.8675;
	y =  128.2654; 
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 682.2515;
	y =  547.2600;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -717.2219;
	y =  432.8312;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 588.0764;
	y =  345.4055;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//std::cout<<pCell->phenotype.phase.index<<std::endl;
	
	// PROLIFERATIVE CTL cell
	x = 736.2712;
	y = -90.0571;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	*/
		
	double GBM_NO = 974;//1.2*1e3;//1.2*1e6;
	//GBM cells 	
	for( int i=0; i<GBM_NO; i++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		
		pCell = create_cell( cancer_cell );
				
		pCell->assign_position( x , y , 0.0 );
		
		
		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			
			double speed_var = UniformRandom();
			
			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
				pCell->phenotype.motility.migration_speed = speed_vec[k];
				k = 12;
				}
			}
		} 
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP
	
		double go_stop_var = UniformRandom();
		
		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}
			
		}
	}

	return;
	
}

void only_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
	//double max_pressure = 0.0001012*268*pCell->phenotype.mechanics.cell_cell_repulsion_strength;
	//double max_pressure = 0.0026;//0.0001012*168*pCell->phenotype.mechanics.cell_cell_repulsion_strength;
	
	double R = 21.5/2;
	double SA = 4*3.1416*R*R;//8.41271*8.41271;
	double s = 24.6523;
	double pressure = 6*(1-1/(2*R)*s)*(1-1/(2*R)*s);
	
	double pressure_scale = 0.027288820670331; 
	
	double max_pressure = pressure/pressure_scale;
	
	//tumour cell proliferation
	if(pCell->type ==2)
	{	
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		double U = all_cells->size();
		//{
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00073549*(1-pCell->state.simple_pressure*pCell->state.simple_pressure/max_pressure);
		
		if( pCell->state.simple_pressure*pCell->state.simple_pressure>max_pressure)
		{pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;}
	
		
	}
	
	cell_movement( pCell, phenotype, dt);
	//infection_dynamics( pCell, phenotype, dt );
		
	return;
	
}

void cell_movement( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	//cell movement
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
	
	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; // 1 = go, 2 = stop
	
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
	
	if( wall_amount>2 )
	{
		//std::cout<<"hit the boundary"<<std::endl;
		//pCell->phenotype.motility.is_motile = false;
		//pCell->is_movable = false;
		
		pCell->phenotype.motility.migration_speed = 0;
	}
	else if(pCell->type != 2)
	{
		pCell->phenotype.motility.migration_speed = 4;
	}
	else		
	{
		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
		{
			// assign new type (stop = 2, or go = 1)
			double new_type_rand = UniformRandom();
			if(new_type_rand<=0.5)// GO
			{
				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type
				
				double speed_var = UniformRandom();
			
				for( int k=0; k<12; )
				{
					if( speed_var> speed_cumul[k] )
					{k++;}
					else
					{
						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
						k = 12;
					}
				}
			} 
			else
			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type

			// assign persistence time - needs to be a real time!
			double go_stop_var = UniformRandom();
			for( int j=0; j<8; )
			{
				if( go_stop_var> go_times_cumul[j] )
				{j++;}
					else
				{
					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
					j = 8;
				}
			}
		}
			
	}
	return;
}

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
			
	//cell infection
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	
	double u = 0.0020276;//uptake rate
	double Vvoxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	double nstar = 10;//infection threshol
	double nu = 0.4886;//repliation rate
	double alp = 6600;//virus burst number
	
	double b = 1;
	
	if(pCell->type == 2)
	{
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		b = 	0.00073549;//pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index );
	}
	else
	{
		b = 	0.00143;
	}

			
		
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
	
	// update uptake rate
	pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(n+nstar/2);
	
		if( n > 1 && n < alp-300 ) // update amount inside due to replication 
		{
				
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n+dt*(nu*b*n); 
				if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<0 )
				{std::cout<<"NEGATIVE HELP1"<<u*p*Vvoxel+nu*n*(1-n/alp)<<n<<" "<<alp<<" "<<nstar<<" "<<dt<<std::endl;}
		}
		else if( n > alp-300)//cell dies
		{
				pCell->start_death( apoptosis_model_index );	
		}
	}
	/*else if( pCell->phenotype.death.dead==true && pCell->custom_data.variables[intracellular_virus_index].value>0)//dead cell secretes virus
	{
		//std::cout<<"Cell died"<<std::endl;	
		pCell->phenotype.secretion.secretion_rates[virus_signal_index] = 0;
		double virus_amount_in_cell = pCell->custom_data.variables[intracellular_virus_index].value;
		if(virus_amount_in_cell/Vvoxel+pCell->nearest_density_vector()[virus_signal_index]>pstar)
		{
			//std::cout<<TRAIL_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[TRAIL_signal_index]<<std::endl;
			double amount_addble = pstar-pCell->nearest_density_vector()[virus_signal_index];
			pCell->nearest_density_vector()[virus_signal_index] += amount_addble;
			pCell->custom_data.variables[intracellular_virus_index].value = virus_amount_in_cell-amount_addble*Vvoxel;
		}
		else
		{pCell->nearest_density_vector()[virus_signal_index] += virus_amount_in_cell/Vvoxel;}
	}*/
	
	return;
}

/*

if(internal_virus<nstar)
	{
			uptake_rate =u*p*Vvoxel;
	}
   else if(internal_virus>=nstar)//infected so not replicating
   {
	    pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;
	   //std::cout<<bias_phenotype<<std::endl;
		pCell->phenotype.molecular.internalized_total_substrates[v_index] = internal_virus+replication_rate*6;//virus replicates at linear rate
		if(pCell->phenotype.molecular.internalized_total_substrates[v_index]>1e2) 
		{pCell->phenotype.secretion.uptake_rates[v_index] = 0;}
   
		if(pCell->phenotype.molecular.internalized_total_substrates[v_index]>1e4)
		{pCell->start_death( apoptosis_model_index );}	   
	
		pCell->custom_data.variables[bias_phenotype_index].value = 0;
   }*/

void virus_dynamics_TRAIL_MODEL3( Cell* pCell, Phenotype& phenotype, double dt )
{
	// PROLIFERATION ^^
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
    //static int pressure_store_index =pCell->phenotype.death.find_death_model_index( "pressure_store" );	
	//pCell->custom_data.variables[pressure_store_index].value = pCell->state.simple_pressure;
		double U = all_cells->size();
	pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00064812*(1-U/5127);	
	
	
	
	// PROLIFERATION ^^
	
	
	// find microenvironment virus index
	//static int virus_signal_index	= microenvironment.find_density_index( "virus" ); 
	
	// find microenvironment TRAIL index
	static int TRAIL_signal_index	= microenvironment.find_density_index( "TRAIL" ); 
	
	// sample amount of virus in microenvironment
	//double virus_amount = pCell->nearest_density_vector()[virus_signal_index];
	// sample amount of TRAIL in microenvironment
	double TRAIL_amount = pCell->nearest_density_vector()[TRAIL_signal_index];
	
	// find intracellular virus index
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	// find intracellular virus index
	static double intracellular_TRAIL_index = pCell->custom_data.find_variable_index( "intracellular_TRAIL_amount" );
	
	// find replication start time index
	static double replication_start_index = pCell->custom_data.find_variable_index( "replication_start_time" );
	
	// find death model index 
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	
	// parameters 
	static double infection_density_capacity = parameters.doubles("infection_density_capacity");
	static double intrinsic_infection_rate = parameters.doubles("intrinsic_infection_rate"); //2; 
	static double viral_replication_rate = parameters.doubles("viral_replication_rate");	
	static double virus_replication_minimum = parameters.doubles("virus_replication_minimum");	
	static double TRAIL_secretion_rate = parameters.doubles("TRAIL_secretion_rate");
	static double TRAIL_generation_rate = parameters.doubles("TRAIL_generation_rate");
	static double TRAIL_killing_level = parameters.doubles("TRAIL_killing_level");
	static double M = parameters.doubles("M");
	
	double c_I = intrinsic_infection_rate;//rate of infection  0.1
	double K = infection_density_capacity;//capacity of cell
	double c_R = 0;//viral_replication_rate; // 0.1
	double c_T = TRAIL_generation_rate; // 0.1
	double n_Istar = virus_replication_minimum;//10
	//double rho_E = virus_amount;
	double s_T =  TRAIL_secretion_rate;
	int alpha = 1;
	double L = 3000;//K/2+1.3*K/2;
	double s_tau = parameters.doubles("TRAIL_secretion_time_commence");
	
	double V_voxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	
	double n_I = pCell->custom_data.variables[intracellular_virus_index].value;
	//double n_E = virus_amount*V_voxel;
	double T_I = pCell->custom_data.variables[intracellular_TRAIL_index].value;
	double replication_start_time = pCell->custom_data.variables[replication_start_index].value;
	double T_E = pCell->nearest_density_vector()[TRAIL_signal_index];
	double pstarT = phenotype.secretion.saturation_densities[TRAIL_signal_index]; 

	double n_T = infection_density_capacity;

	double n_I_zero = 0;
	if( n_I < n_T )
	{n_I_zero = n_I;}
	else
	{n_I_zero = n_T;}
	
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
		if( replication_start_time = 0)
		{
			pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
		}
		
		if( PhysiCell_globals.current_time>s_tau+replication_start_time )
		{
			pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I+dt*(c_T*n_I_zero-s_T*T_I*(pstarT-T_E));
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
			//std::cout<<"here 2 "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<" "<<s_T*T_I/V_voxel<<" "<<c_T*n_I_zero-s_T*T_I*(pstarT-T_E)<<std::endl;
		}
		else
		{
			pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I+dt*(c_T*n_I_zero);
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
		}
		
		if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
		{	
			std::cout<<"NEGATIVE HELP TRAIL7 "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<" "<<T_I<<" "<<c_T<<std::endl;	
			system("pause");
		}
	}
	else if( pCell->phenotype.death.dead==true && pCell->custom_data.variables[intracellular_TRAIL_index].value>0)
	{
		//std::cout<<"Cell died"<<std::endl;	
		pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
		double TRAIL_amount_in_cell = pCell->custom_data.variables[intracellular_TRAIL_index].value;
		if(TRAIL_amount_in_cell/V_voxel+pCell->nearest_density_vector()[TRAIL_signal_index]>pstarT)
		{
			//std::cout<<TRAIL_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[TRAIL_signal_index]<<std::endl;
			double amount_addble = pstarT-pCell->nearest_density_vector()[TRAIL_signal_index];
			pCell->nearest_density_vector()[TRAIL_signal_index] += amount_addble;
			pCell->custom_data.variables[intracellular_TRAIL_index].value = TRAIL_amount_in_cell-amount_addble*V_voxel;
		}
		else
		{pCell->nearest_density_vector()[TRAIL_signal_index] += TRAIL_amount_in_cell/V_voxel;}
	}
	
	/*	if( pCell->phenotype.death.dead==false && virus_amount>1e-2)
		{// cell becomes infected
	
			//std::cout<<n_I<<" "<<n_Istar<<std::endl;
			if( n_I > n_Istar && n_I< K )
			{
				
				//std::cout<<"NEGATIVE HELP0"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
				
				
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
				if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
				{
					std::cout<<"NEGATIVE HELP TRAIL"<<std::endl;
					system("pause");
				}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL7"<<T_I<<c_T<<std::endl;
						system("pause");
					}
				}
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP1"<<n_I+dt*(c_I*n_E*(1-n_I/K)+c_R)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<" "<<n_E<<" "<<dt<<" "<<c_I<<" "<<c_R<<std::endl;
					system("pause");
				}
			}
			else if( n_I > K)
			{
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL6"<<std::endl;
						system("pause");
					}
								}
								else
								{
									pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);//-s_T*T_I);
									pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;//s_T*T_I/V_voxel;
									
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL5"<<std::endl;
						system("pause");
					}
				}
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP2"<<" "<<pCell->custom_data.variables[intracellular_virus_index].value<<" "<<n_I+dt*(c_R)<<" "<<T_I+dt*(c_T-s_T*T_I)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
			}
			else
			{
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP3"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
				
			}
			}
		else if( n_I>0 && pCell->phenotype.death.dead==false )
		{// cell becomes infected
	
			
	
			if( n_I > n_Istar && n_I< K )
			{
				
				//std::cout<<"NEGATIVE HELP0"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
				
				
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL4"<<std::endl;
						system("pause");
					}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL3"<<T_I<<" "<<c_T<<" "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<std::endl;
						system("pause");
					}
				}
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP1"<<n_I+dt*(c_I*n_E*(1-n_I/K)+c_R)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<" "<<n_E<<" "<<dt<<" "<<c_I<<" "<<c_R<<std::endl;
					system("pause");
				}
			}
			else if( n_I > K)
			{
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL2"<<std::endl;
						system("pause");
					}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);//-s_T*T_I);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;//s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL1"<<std::endl;
						system("pause");
					}
				}
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP2"<<" "<<pCell->custom_data.variables[intracellular_virus_index].value<<" "<<n_I+dt*(c_R)<<" "<<T_I+dt*(c_T-s_T*T_I)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
			}
			else
			{
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP3"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
				
			}
			
			/*if( n_I>L/5 && n_I<L )
			{
				pCell->phenotype.death.rates[apoptosis_model_index] = n_I*n_I*n_I/(L/2*L/2*L/2+n_I*n_I*n_I);
				
			}
			else if(n_I > L)
			{
				//pCell->phenotype.death.rates[apoptosis_model_index] = 1000;
				pCell->start_death( apoptosis_model_index );				
				
			}
			else
			{
				pCell->phenotype.death.rates[apoptosis_model_index]=0;
			}*/
		//}
		/*else if( T_E>1e-5 && pCell->phenotype.death.dead==false)
		{
			double M = 0.0005;
			pCell->phenotype.death.rates[apoptosis_model_index] = T_E*T_E/(M*M+T_E*T_E);
			//if (pCell->phenotype.death.rates[apoptosis_model_index]>0.01)
			//std::cout<<"apop rate"<<pCell->phenotype.death.rates[apoptosis_model_index] <<std::endl;
			
		}*/
		/*else if( pCell->phenotype.death.dead==true)
		{
			//std::cout<<"Cell died"<<std::endl;
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->phenotype.secretion.secretion_rates[virus_signal_index] = 0;		
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
			if( pCell->custom_data.variables[intracellular_virus_index].value>0)
				{
					double virus_amount_in_cell = pCell->custom_data.variables[intracellular_virus_index].value;
					//std::cout<<"Cell dies: " << virus_amount_in_cell<<std::endl;
					pCell->phenotype.secretion.saturation_densities[virus_signal_index] = 5;
					pstarV = 5;
					double TRAIL_amount_in_cell = pCell->custom_data.variables[intracellular_TRAIL_index].value;
					if(TRAIL_amount_in_cell/V_voxel+pCell->nearest_density_vector()[TRAIL_signal_index]>pstarT)
					{
						std::cout<<TRAIL_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[TRAIL_signal_index]<<std::endl;
						double amount_addble = pstarT-pCell->nearest_density_vector()[TRAIL_signal_index];
						pCell->nearest_density_vector()[TRAIL_signal_index] += amount_addble;
						pCell->custom_data.variables[intracellular_TRAIL_index].value = TRAIL_amount_in_cell-amount_addble*V_voxel;
					}
					else
					{pCell->nearest_density_vector()[TRAIL_signal_index] += TRAIL_amount_in_cell/V_voxel;}
					
					if( virus_amount_in_cell/V_voxel+pCell->nearest_density_vector()[virus_signal_index]>pstarV )
					{
						//std::cout<<virus_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[virus_signal_index]<<std::endl;
						double amount_addble = pstarV-pCell->nearest_density_vector()[virus_signal_index];
						pCell->nearest_density_vector()[virus_signal_index] += 0;//amount_addble;
						pCell->custom_data.variables[intracellular_virus_index].value = virus_amount_in_cell-amount_addble*V_voxel;
					
					}
					else
					{pCell->nearest_density_vector()[virus_signal_index] += 0;//virus_amount_in_cell/V_voxel;
					pCell->custom_data.variables[intracellular_virus_index].value=0;
					//pCell->phenotype.secretion.saturation_densities[virus_signal_index] = 5;
					pCell->custom_data.variables[intracellular_virus_index].value=0;
					}
				}
		}*/
	
	/*int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	if (pCell->state.simple_pressure > 10)//30
	{
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0;
	}
	else 
	{
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0001428;
	}
	//std::cout<<pCell->state.simple_pressure <<std::endl;
	
	*/
	if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
	{
		std::cout<<"NEGATIVE HELP TRAIL"<<std::endl;
	}
	
		
	 if( T_E>0.01 && pCell->phenotype.death.dead==false)
	{
		double M = 0.8;
		pCell->phenotype.death.rates[apoptosis_model_index] = (T_E*T_E*T_E/(M*M*M+T_E*T_E*T_E))*0.05;
	}
		
		
	
	return;
}

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int v_index = microenvironment.find_density_index( "virus");
	static int infection_density_capacity = parameters.doubles("infection_density_capacity");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	double p_min = 1;
	double p_max = 6600;
	
	
	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	
	if( n_I > 6600)
	{
		std::cout<<n_I<<std::endl;	
	}	
	
	if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			
			/*if( bias_phenotype<0.9 && n_I>1)
			{
				
				int oncoprotein = (int) round( (1.0/(p_min-p_max)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_max) * 230.0 ); 
				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein,255); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
				output[0].assign( szTempString );
				output[1]="blue";
				output[2].assign( szTempString );
				output[3]="blue";
							
				return output;
			}
			else if( pCell->phenotype.molecular.internalized_total_substrates[v_index] <0)
			{
				std::cout<<"NEGATIVE"<<std::endl;
				std::vector< std::string > output( 4, "darkgrey" ); 
				return output;
			}
			else
			{*/
		if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{ 
					output[0] = "orange";
					output[1] = "orange";
					output[2] = "coral";
					output[3] = "coral";
					return output; 
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	

					output[0] = "darkred";
					output[1] = "darkred";
					output[2] = "firebrick";
					output[3] = "firebrick";
			}		
			//}
		}
		if( pCell->type == 3)
		{
			if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{   
			    output[0] = "aquamarine";
				output[1] = "lightsteelblue";
				output[2] = "lightskyblue";
				output[3] = "aquamarine";
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	
		
				output[0] = "darkslateblue";
				output[1] = "darkblue";
				output[2] = "darkblue";
				output[3] = "aquamarine";
			}
		}
		if( pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if( n_I>1)//n_I > 1 )
			{
				double p_min = 1;
				double p_max = 6600;//V_0;
				int oncoprotein1 = (int) round((210-139)*(1.0/(p_max-p_min)) * (n_I-1)); 
				int oncoprotein2 = (int) round((180-69)*(1.0/(p_max-p_min)) * (n_I-1)); 
				int oncoprotein3 = (int) round((140-19)*(1.0/(p_max-p_min)) * (n_I-1)); 
				
				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", 210-oncoprotein1, 180-oncoprotein2, 140-oncoprotein3); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
				
				output[0].assign( szTempString );
				output[1]="brown";
				output[2].assign( szTempString );
				output[3]="brown";
				
				return output;
			}
			else
			{
					output[0] = "orchid";//"rgb(255,230,230)";
					output[1] = "orchid";
					output[2] = "plum";//"rgb(255,230,230)";
					output[3] = "plum";

					return output; 
					 
			}
		}
	if( pCell->phenotype.death.dead == true )
	{ 
			output[0] = "rgb(255, 255, 224)";
			output[1] = "rgb(255, 255, 224)";
			output[2] = "rgb(255, 228, 181)";
			output[3] = "rgb(255, 228, 181)";
			
		
		return output; 
	}
	 
	return output;
}

