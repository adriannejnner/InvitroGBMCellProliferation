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

#include "./GBM_OV_immune_stroma_double_immune.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>


Cell_Definition TH_cell; 
Cell_Definition CTL_cell; 
Cell_Definition cancer_cell; 
Cell_Definition stroma_cell; 

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
	TH_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate"); 	
	TH_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("TH_quiescent_transistion_rate"); 

	// cell actions
	TH_cell.phenotype.secretion.uptake_rates[1] = 0.0;
	TH_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	TH_cell.phenotype.motility.migration_speed = parameters.ints("TH_migration_speed");//4
	
	// cell morphology
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
	
	// update function
	TH_cell.functions.update_phenotype = TH_functions;
	
	return;
}
void create_CTL_cells( void )
{
	CTL_cell = cell_defaults;
	
	CTL_cell.name = "CTL cell";
	CTL_cell.type = 3;
	
	// proliferation 
	CTL_cell.functions.cycle_model = Ki67_basic;
	CTL_cell.phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive); 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate");//7.2206*1e-5; 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("CTL_quiescent_transistion_rate");//0.00143; 	
	
	// cell actions
	CTL_cell.phenotype.secretion.uptake_rates[1] = 0.0;
	CTL_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	CTL_cell.phenotype.motility.migration_speed = parameters.ints("CTL_migration_speed");
	
	//cell morphology
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
	
	// cell update phenotype		
	CTL_cell.functions.update_phenotype = CTL_functions;
	
	
	return;
}void create_stroma_cells( void )
{
	stroma_cell = cell_defaults;
	
	stroma_cell.name = "stroma cell";
	stroma_cell.type = 4;
	
	// turn off proliferation 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live); 
	stroma_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0; 

	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	stroma_cell.phenotype.secretion.uptake_rates[virus_index] = parameters.doubles("stroma_virus_uptake_rate")*0.5;
	stroma_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	stroma_cell.phenotype.motility.migration_speed = 0;
	stroma_cell.phenotype.motility.is_motile = false;
	stroma_cell.phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;
	
	// cell morphology
	stroma_cell.phenotype.geometry.radius = parameters.doubles("stroma_radius");//7.5;
	stroma_cell.phenotype.volume.total = 1767;
	stroma_cell.phenotype.volume.fluid_fraction = 0.75;
	stroma_cell.phenotype.volume.fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.total;
	stroma_cell.phenotype.volume.solid = stroma_cell.phenotype.volume.total-stroma_cell.phenotype.volume.fluid;
	stroma_cell.phenotype.volume.nuclear = 500;
	stroma_cell.phenotype.volume.nuclear_solid = 125;
	stroma_cell.phenotype.volume.nuclear_fluid = stroma_cell.phenotype.volume.nuclear - stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.cytoplasmic = stroma_cell.phenotype.volume.total - stroma_cell.phenotype.volume.nuclear;
	stroma_cell.phenotype.volume.cytoplasmic_fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.cytoplasmic;
	stroma_cell.phenotype.volume.cytoplasmic_solid = stroma_cell.phenotype.volume.cytoplasmic-stroma_cell.phenotype.volume.cytoplasmic_fluid;
	stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 2.53;
	stroma_cell.phenotype.volume.target_solid_cytoplasmic = stroma_cell.phenotype.volume.cytoplasmic_solid;
	stroma_cell.phenotype.volume.target_solid_nuclear = stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.target_fluid_fraction = stroma_cell.phenotype.volume.fluid_fraction;
	stroma_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update phenotype
	stroma_cell.functions.update_phenotype = stroma_function;
	
	return;
}
void create_cancer_cells( void )
{
	
	cancer_cell = cell_defaults;
	
	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	cancer_cell.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cancer_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	cancer_cell.phenotype.motility.migration_speed = 0.0;
	
	// cell morphology
	cancer_cell.phenotype.geometry.radius = 10.75;
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
	
	//update phenoypt
	cancer_cell.functions.update_phenotype = cancer_cell_proliferation_infection_movement;
	
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
	
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;//0.00073549;//0.0000064812*(1-9740/512720);
	
 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 

	// reduce cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.0;//0.05;
		
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 ); // amount of intracellular virus
	cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 ); // how long cells will persist in move or stop phenotype
	cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "rep_rate", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached 
	cell_defaults.custom_data.add_variable( "special_virus_uptakerate" ,"min",0); // creating a distribution for the uptake rates
	cell_defaults.custom_data.add_variable( "special_virus_replication_rate","min",0.0); 
	
	Parameter<double> paramD;
	
	paramD = parameters.doubles[ "elastic_coefficient" ]; 
	cell_defaults.custom_data.add_variable( "elastic coefficient" , paramD.units, paramD.value ); 
	
	// turn off secretion from these cells (oxygen and virus)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	static int virus_index = microenvironment.find_density_index( "virus");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = parameters.doubles("virus_saturation_density");//;//10; 
	cell_defaults.phenotype.molecular.fraction_released_at_death[virus_index] = 1;//1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[0] = 1;//1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[2] = 1;//1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[3] = 1;//1;
	
	static int wall_index = microenvironment.find_density_index( "wall");
	
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
	
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	cell_defaults.phenotype.secretion.secretion_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[chemokine_index] = parameters.ints("chemokine_saturation_density");//5; 
	
	// update cell and phenotype based on virus dynamics only
	
	cell_defaults.functions.update_phenotype = cancer_cell_proliferation_infection_movement;
	cell_defaults.phenotype.motility.is_motile = true; 
		
	cell_defaults.name = "holder cell"; 
	cell_defaults.type = 0; 
	
	//create cell types
	create_cancer_cells();
	create_CTL_cells();
	create_TH_cells();
	create_stroma_cells();
		
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
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")+10)*(parameters.doubles("tumour_radius")+10))//1280*1280)//420*420)//
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 1;
			microenvironment(n)[chemokine_index] = 0;
			
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")-10)*(parameters.doubles("tumour_radius")-10))//420*420)//1260*1260)//1260*1260)400*400)//
		{
			
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[virus_index] = parameters.doubles("initial_virus_density");
			microenvironment(n)[wall_index] = 3.5;
			microenvironment(n)[chemokine_index] = 0;
		}
		else if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>250*250 )
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 5;
			microenvironment(n)[chemokine_index] = 0;
				
		}
		else
		{   
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 10;
			microenvironment(n)[chemokine_index] = 0;
		}
			
	}
	return; 
}

void setup_tissue_circle_immune( void )
{
	double Radius = parameters.doubles("tumour_radius");
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	// setting up distributions for movement and persistance of cells
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
				
	immune_cell_placement();
		
	double GBM_NO = parameters.ints("initial_GBM_cells");
	double stroma_NO = parameters.ints("initial_stroma_cells");
	std::default_random_engine generator;
/*
	std::default_random_engine generator;
	double nu_mean1 = parameters.doubles("GBM_virus_uptake_rate")*10;
	double nu_variance1 = 0.0001;
	double shape1 = nu_mean1*nu_mean1/nu_variance1;
	double scale1 = nu_variance1/nu_mean1;
	std::gamma_distribution<double> distribution1(shape1,scale1);
	
	double nu_mean2 = parameters.doubles("stroma_virus_uptake_rate")*10;
	double nu_variance2 = 0.0001;
	double shape2 = nu_mean2*nu_mean2/nu_variance2;
	double scale2 = nu_variance2/nu_mean2;
	std::gamma_distribution<double> distribution2(shape2,scale2);
*/

	double nu_mean = 10;
	double nu_variance = 0.0001;
	double shape = nu_mean*nu_mean/nu_variance;
	double scale = nu_variance/nu_mean;
	std::gamma_distribution<double> distribution(shape,scale);
		
	//GBM cells 	
	for( int i=0; i<GBM_NO; i++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		
		pCell = create_cell( cancer_cell );		
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index( "virus");
		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
		int	virus_uptake_rate_index = pCell->custom_data.find_variable_index( "special_virus_uptakerate"); 
		int	virus_replication_rate_index = pCell->custom_data.find_variable_index( "special_virus_replication_rate"); 
				
		pCell->custom_data.variables[virus_uptake_rate_index].value = parameters.doubles("GBM_virus_uptake_rate")*distribution(generator);//distribution1(generator);
		pCell->custom_data.variables[virus_replication_rate_index].value = parameters.doubles("virus_replication_rate");//distribution(generator);
		
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
	
	for( int l=0; l<stroma_NO; l++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		pCell = create_cell( stroma_cell );
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index("virus");
		//pCell->phenotype.secretion.uptake_rates[virus_signal_index] = NormalRandom( 0.2, 0.3);
		//if( pCell->phenotype.secretion.uptake_rates[virus_signal_index] < 0)
		//{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;}
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("stroma_virus_uptake_rate")*distribution(generator);//distribution2(generator);
		
		int	virus_uptake_rate_index = pCell->custom_data.find_variable_index( "special_virus_uptakerate"); 
		pCell->custom_data.variables[virus_uptake_rate_index].value = pCell->phenotype.secretion.uptake_rates[virus_signal_index];
		
	}
	
	return;
	
}

void cancer_cell_proliferation_infection_movement( Cell* pCell, Phenotype& phenotype, double dt )
{

	double R = 21.5/2;
	double SA = 4*3.1416*R*R;
	double s = parameters.doubles("maximum_cell_density");
	double pressure = 6*(1-1/(2*R)*s)*(1-1/(2*R)*s);
	double pressure_scale = 0.027288820670331; 
	double max_pressure = pressure/pressure_scale;
	
	//tumour cell proliferation
	if(pCell->type ==2)
	{	
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("GBM_cell_proliferation_rate");

		if( pCell->state.simple_pressure*pCell->state.simple_pressure>max_pressure) // if cell under too much pressure -> no proliferation
		{pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;}
	
	}
	
	cell_movement( pCell, phenotype, dt);
	infection_dynamics( pCell, phenotype, dt );
		
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
	
	if( wall_amount<2 & pCell->type == 2 )
	{
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

	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
		
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	// double u = parameters.doubles("GBM_virus_uptake_rate");//0.0020276;//uptake rate
	int index_specal_uptake = pCell->custom_data.find_variable_index( "special_virus_uptakerate" );
	double u = pCell->custom_data.variables[index_specal_uptake].value;
	double Vvoxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	double nstar = parameters.doubles("infection_threshold");//10;//infection threshol
	//double nu = parameters.doubles("virus_replication_rate");//10*0.4886;//repliation rate
	double alp = parameters.doubles("virus_burst_number");//1000;//virus burst number
	
	
	int	virus_replication_rate_index = pCell->custom_data.find_variable_index( "special_virus_replication_rate"); 
	double nu = pCell->custom_data.variables[virus_replication_rate_index].value;
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];		
		
	double pmax = parameters.doubles("pmax");//0.0125;
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
		if(p<pmax)
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(n/Vvoxel+nstar/Vvoxel);}
		else
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*pmax*pmax/(n/Vvoxel+nstar/Vvoxel)/p;}
		
		if( n > 1 && n <= alp) // update amount inside due to replication 
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n+dt*(nu*n); 
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		
		}
		else if( n > alp-1)
		{
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
			pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->start_death( apoptosis_model_index );	
		}
		else if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<0 )
		{std::cout<<"NEGATIVE INTRACELLULAR VIRUS"<<std::endl;}
	}
	if( pCell->phenotype.death.dead == true && pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index]>0)
	{
		virus_induced_lysis(pCell, phenotype, dt );
	}
	
	return;
}

void virus_induced_lysis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double pstar = pCell->phenotype.secretion.saturation_densities[virus_signal_index];
	double delta_V = parameters.doubles("viral_decay_from_burst_cell");//0.1466;
	double Vvoxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	
	if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]> 1 ) 
	{
			double amount_to_add = (n-n*exp(-delta_V*dt))/Vvoxel;
			if( amount_to_add > pstar-p )
			{
				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
			}
			else
			{
				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n*exp(-delta_V*dt);
			}	
			//std::cout<<amount_to_add<<" old p: "<<p<<" new p: "<<pCell->nearest_density_vector()[virus_signal_index]<<" old n: "<<n<<" new n: "<<pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<<std::endl;
			
	}

	return;
}

void stroma_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	if(pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]>1e5)
	{
		pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
		pCell->custom_data.variables[intracellular_virus_index].value = 0;
		//stroma_cell.phenotype.secretion.uptake_rates[virus_signal_index] = 0;
	}
	else
	{
		pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	}
			
	//std::cout<<pCell->phenotype.secretion.uptake_rates[virus_signal_index]<<" "<<pCell->custom_data.variables[intracellular_virus_index].value<<std::endl;
	return;
}


void TH_functions( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	std::vector<double> ae_ini(3);
    
	
	
	if( wall_amount<2 )// Make TH cells that have left the diameter of the tumour make cell turn around
	{
		pCell->phenotype.motility.migration_speed = 4;
		ae_ini = -1*pCell->position;
		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
	}
	else 
	{
		pCell->phenotype.motility.migration_speed = 4;
		pCell->phenotype.motility.migration_bias = 0;
	}
	
	// TH secretion of cytokines
	std::vector<Cell*> nearby = pCell->cells_in_my_container();
	
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	Cell* pC = NULL;
	bool stop = false;
	int i=0;
	
	double nstar = parameters.doubles("infection_threshold");

	while( !stop && i < nearby.size() )
	{
		pC = nearby[i];
		if( pC->phenotype.molecular.internalized_total_substrates[virus_signal_index] > nstar &&
			pC->phenotype.death.dead == false &&
			pC != pCell && pC->type != 4)
			{ stop = true; }
		
		i++;
	
		if( stop == false )
			{ pC = NULL; }
	}
	
	if( pC )
	{
		//std::cout << "infected cell found by TH" << std::endl;
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = parameters.doubles("chemokine_secretion_rate");//0.0417;//0.1;	
			
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate")*parameters.doubles("TH_prolif_increase_due_to_stimulus"); 	
		pCell->phenotype.motility.migration_speed = 0.1;
	
	}
	else
	{
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;//0.1;	
			
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate"); 	
		//pCell->phenotype.motility.is_motile = true;
	}
	
	// TH increased proliferation ?
	
	
	return;
}

void CTL_functions( Cell* pCell, Phenotype& phenotype, double dt )
{	

	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	if( pCell->type != 3 )
	{std::cout<<"wrong cell type"<<std::endl;}

	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment lifetime" ); 
	
	// is CTL docked to infected cell - increase prolif - time for docking
	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
		
		// attempt to kill my attached cell
		bool dettach_me = false; 
		
		if( immune_cell_attempt_apoptosis( pCell, pCell->state.neighbors[0], dt ) )
		{
			
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] ); 
			dettach_me = true; 
		}
			
		// if I dettach, resume motile behavior 
		if( dettach_me )
		{
			
			dettach_cells( pCell, pCell->state.neighbors[0] ); 
			phenotype.motility.is_motile = true; 
			CTL_movement( pCell, phenotype, dt);
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate")*parameters.doubles("CTL_prolif_increase_due_to_stimulus"); 	
			//std::cout<<"detached and should be proliferating"<<std::endl;
			
			
		}
		return; 
	}
	
	//has CTL encounted a new infected cell - start docking clock
	// if this returns non-NULL, we're now attached to a cell 
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	
	// is there a gradient for the CTL to follow - start random walk in the direction of the chemattractant
	phenotype.motility.is_motile = true; 
	CTL_movement( pCell, phenotype, dt);
	return;
}
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	static double max_attachment_distance = parameters.doubles("max_attachment_distance"); // 18.0; 
	
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double internal_virus = pTarget->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	
	double kill_time = parameters.doubles("time_to_kill_cell"); // how long the cell needs to attach for  the infected cell to be killed
	double nstar = parameters.doubles("infection_threshold"); // how long the cell needs to attach for  the infected cell to be killed
	// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX CHANGED SO THAT CTLS KILL ANY CELL
	//if( internal_virus > nstar && pTarget->phenotype.death.dead == false && pTarget->type!=4)
	if( pTarget->phenotype.death.dead == false && pTarget->type!=4)
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 
	
	
		attach_cells( pAttacker, pTarget ); 
		pAttacker->custom_data[attach_lifetime_i] = PhysiCell_globals.current_time + kill_time;
		
		return true; 
	}
	
	return false; 
}

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return; 
}

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

void CTL_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	
	std::vector<double> ae_ini(3);
    
	
	// TH movement
	if( wall_amount<2 )
	{
		pCell->phenotype.motility.migration_speed = 1;
		pCell->phenotype.motility.migration_bias = 1;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(wall_index);
		ae_ini = -1*pCell->position;
	
		pCell->phenotype.motility.migration_bias = 1;
		//std::cout<<pCell->nearest_gradient(wall_index)<<" "<<pCell->position<<" "<<ae_ini<<" "<<wall_amount<<std::endl;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
		
		return;
	}
	else if(chemokine_amount>1e-8)// sample chemotaxis gradient and random walk in that direction
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("CTL_min_speed")+(parameters.doubles("CTL_max_speed")-parameters.doubles("CTL_min_speed"))*(chemokine_amount/(parameters.doubles("Chemokine_EC50")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("CTL_chemokine_migration_bias");//0.85;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate");
		
		return;
	}
	else if(chemokine_amount>1e-3)
	{
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate")*parameters.doubles("CTL_prolif_increase_due_to_stimulus");//7.9026*1e-5*1e1; 	
		pCell->phenotype.motility.migration_speed = parameters.doubles("CTL_min_speed")+(parameters.doubles("CTL_max_speed")-parameters.doubles("CTL_min_speed"))*(chemokine_amount/(parameters.doubles("Chemokine_EC50")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("CTL_chemokine_migration_bias");//0.85;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
	}
	else
	{
		pCell->phenotype.motility.migration_bias = 0;
		pCell->phenotype.motility.migration_speed = 4;
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.9026*1e-5;
		return;
	}
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	
	// CTL kills cell if it has been attached for enough time
	if( pAttacker->custom_data[attach_lifetime_i] < PhysiCell_globals.current_time )
	{ 
		//std::cout << "cell killed" << std::endl; 
		return true; 
	}
	return false; 
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	static int virus_index = microenvironment.find_density_index( "virus" ); 
	
	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	pTarget->phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;
	return true; 
}

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int v_index = microenvironment.find_density_index( "virus");
	static int infection_density_capacity = parameters.doubles("infection_density_capacity");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	double p_min = 1;
	double p_max = parameters.doubles("virus_burst_number");
	
	
	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	if(n_I>7000 && pCell->type ==2)
	{std::cout<<" above 6600: "<<n_I<<std::endl;}
	
	if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			
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
				//double p_max = 6600;//V_0;
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
					output[0] = "rgb(104, 55, 99)";//"orchid";//"rgb(255,230,230)";
					output[1] = "rgb(104, 55, 99)";
					output[2] = "rgb(85, 50, 70)";//"plum";//"rgb(255,230,230)";
					output[3] = "rgb(85, 50, 70)";

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
	if( pCell->type == 4)
	{
		output[0] = "rgb(234, 172, 199)";//"rgb(255,230,230)";
		output[1] = "rgb(234, 172, 199)";
		output[2] = "rgb(243, 186, 211)";//"rgb(255,230,230)";
		output[3] = "rgb(243, 186, 211)";
	 
	}
	return output;
}

void immune_cell_placement(void)
{
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	//UNPROLIFERATIVE TH CELLS
	//
	double x = -886.1874 ;
	double y = -377.9556;
	
	Cell* pCell = NULL; 
	
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -92.9598;
	y = -408.5643;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 715.7926;
	y = 108.4208;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 305.4827;
	y = -100.0779;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 914.8197;
	y = 334.7098;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 818.4357;
	y = 102.3158;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -657.9876;
	y = 552.5539;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 396.205;
	y = 790.0036;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -197.4969;
	y = -709.7589;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -591.2807;
	y = -732.2582;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 49.2169;
	y = -452.7280;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -26.8695;
	y = 465.632;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -289.5373;
	y = 513.8335 ;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -315.2772;
	y = -772.8922;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -120.1230;
	y = -608.1192;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -312.1967;
	y = 974.7604;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 924.1434;
	y = 206.5514;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -329.2369;
	y = 386.8689;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 434.6186;
	y = 366.5498;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -109.1024;
	y = -992.4774;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 588.8850;
	y = 88.1708;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 293.6073;
	y = -769.3062;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -758.5694;
	y = 469.0668;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 711.9720;
	y = 522.8015;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -202.5984;
	y = 317.7974;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -550.3849;
	y = 31.0353;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 641.0909;
	y = -685.7513;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 208.1485;
	y = -182.9283;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -773.5824;
	y = 598.8495;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 589.6732;
	y = -64.6089;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -631.5368;
	y = -88.5340;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 358.2984;
	y = -894.3028;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -65.4904;
	y = -508.6133;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 422.3507;
	y = -808.1449;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  69.0509;
	y = -696.4078;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 188.2291;
	y = -34.1874;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -890.0763;
	y = 425.9101;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 749.9007;
	y = -614.8952;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -69.6065;
	y = -283.4312;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  375.0652;
	y = 20.8637;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 14.8366;
	y = -587.9312;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -338.3966;
	y = -646.7752;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -617.0420;
	y = -438.5693;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -808.1436;
	y = -480.9575;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -11.5636;
	y = -922.3531;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -239.5675;
	y = -578.7236;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -521.0613;
	y = 768.1352;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 107.7491;
	y = 797.8105;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -447.2199;
	y = -727.9198;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 105.9531;
	y = 441.1681;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 644.2658;
	y = -538.5330;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 908.9122;
	y = 26.160;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 701.5531;
	y = -4.7266;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 509.1125;
	y = 700.8723;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 22.525;
	y = 700.6533;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 840.2498;
	y = 552.8631;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 285.0359;
	y = -535.5775;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-758.6877;
	y = -605.8141;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =285.0359;
	y = -535.5775;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-758.6877;
	y =-605.8141;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =826.7548;
	y =-303.953;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =142.712;
	y =668.9715;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-566.4609;
	y =651.0637;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-136.700;
	y = 828.8107;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-310.4904;
	y =249.8461;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =249.6918;
	y =605.5657;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =1.0007e+03;
	y = 112.6398;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 223.6434;
	y = -885.2210;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 49.8545;
	y = 571.3492;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -100.4299;
	y = 233.0141;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 109.5735;
	y = -954.7865;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	pCell->phenotype.cycle.data.current_phase_index = 0;
	
	x = -622.7696;
	y = 428.8238;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -950.7404;
	y = 289.7073;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -746.3485;
	y =  -354.7596;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -58.8064;
	y =  -710.556;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
		
	///////////////////////////////////////////////////////////////////////////////////////////////
	//PROLIFERATIVE TH CELLS
	x = 18.3263;
	y = 202.6912;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	x = 85.522;
	y = 928.9369;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	x = 294.9299;
	y = 408.8944;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	x = 490.6082;
	y = -357.8950;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	///////////////////////////////////////////////////////////////////////////////////////////////
	//UNPROLIFERATIVE CTL CELLS
	x = -229.7546;
	y = -99.6055;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 357.8347;
	y = 634.1085;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	
	x = 563.1912;
	y = -209.8626;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -465.3017;
	y = -404.2420;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -739.3612;
	y = -195.8691;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -310.6741;
	y = 91.5128;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-836.8711;
	y = 5.789;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =170.3629;
	y =296.8521;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =502.7784;
	y =-503.5691;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =717.2875;
	y =-384.5512;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-293.0794;
	y =-299.6722;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =660.4160;
	y =663.8240;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-619.8710;
	y =-589.8832;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =426.8735;
	y =928.6956;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-167.0797;
	y =687.6530;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-1.0237e3;
	y =137.0964;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =502.0188;
	y =229.9418;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =329.8713;
	y =168.4815;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-263.1634;
	y =-450.6436;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =144.7437;
	y =101.6446;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =395.2043;
	y =-576.2011;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-430.9143;
	y =850.4354;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =76.5861;
	y =-245.4335;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-456.6180;
	y =-557.6040;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-55.1681;
	y =-23.9806;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =806.7518;
	y = 242.4464;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-644.9479;
	y =796.8378;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-157.3983;
	y =509.2217;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =662.5205;
	y =290.8544;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =583.3736;
	y =843.4895;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-107.6750;
	y = -839.4925;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-423.3672;
	y =579.6212;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-432.9000;
	y =-70.883;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-530.5991;
	y =224.9672;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =675.3453;
	y =394.918;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-321.6393;
	y =731.7346;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =503.9698;
	y =498.1822;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-749.3654;
	y =210.6520;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 194.5936;
	y =-670.9127;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 499.6782;
	y =-731.5829 ;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =337.0153;
	y =-298.361;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-479.2653;
	y =412.2822;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =-479.2653;
	y =412.2822;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// PROLIFERATIVE CTL cell
	x = -340.8754;
	y = -928.9261;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	
	x = 244.3090;
	y = -461.5348;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	//*/
	return;
}
