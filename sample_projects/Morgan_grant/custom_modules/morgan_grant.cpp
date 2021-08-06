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

#include "./morgan_grant.h"
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
	TH_cell.phenotype.motility.migration_speed = 0.5;
	//TH_cell.phenotype.geometry.radius = 7.5/2;
	
	TH_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.9026*1e-4; 	
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
	CTL_cell.phenotype.motility.migration_speed = 0.5;
	//CTL_cell.phenotype.geometry.radius = 9.7/2;
	
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.2206*1e-4; 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = 0.00143; 	
	
	return;
}
void create_cancer_cells( void )
{
	cancer_cell = cell_defaults;
	
	cancer_cell.phenotype.secretion.uptake_rates[1] = 0.0;
	cancer_cell.phenotype.secretion.secretion_rates[1] = 0.0;
	cancer_cell.phenotype.geometry.radius = 21.5/2;
	cancer_cell.name = "cancer cell";
	cancer_cell.type = 2; 
	
	return;
}

void create_cell_types( void )
{
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	double C0 = parameters.doubles("C0");
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0008648125;//0.0000064812*(1-9740/512720);
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 
	
	// reduce cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.05;
	
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "bias_phenotype", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "replication_start_time", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "infection_time_start", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "infection_time_length", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "time_of_injection", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "uval", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "nuval", "dimensionless", -1.0 );
	
	// turn off secretion from these cells (oxygen and virus)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	static int virus_index = microenvironment.find_density_index( "virus");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = 10; 
	
	static int wall_index = microenvironment.find_density_index( "wall");
	
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
	
	// update cell and phenotype based on virus dynamics only
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3; 
	
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3;//ECM_movement;//virus_dynamics; 
	cell_defaults.functions.update_phenotype = infection_dynamics;
	cell_defaults.phenotype.motility.is_motile = true; 
	
	//create the immune and cancer cell types
	create_TH_cells();
	create_CTL_cells();
	create_cancer_cells();
		
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
	
	microenvironment.diffusion_coefficients[virus_index] = 0.1766; 
	microenvironment.decay_rates[virus_index] = 0;
	
	microenvironment.diffusion_coefficients[oxygen_index] = 0; 
	microenvironment.decay_rates[oxygen_index] = 0;
	
	microenvironment.diffusion_coefficients[wall_index] = 0; 
	microenvironment.decay_rates[wall_index] = 0;
	
	std::vector<double> ve(3);
    ve[0] = 0;
	ve[1] = 60;
	ve[2] = 10;
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		if( ECMdense[1]<-678.33/2+5) //white
		{microenvironment(n)[virus_index] = 0;}//600;}	
		else
		{microenvironment(n)[virus_index] = 0;}
		
		microenvironment(n)[oxygen_index] = 0;//-2/678.33*(ECMdense[1]+678.33/2)+10;
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>1144.4*1144.4)
		{	
			microenvironment(n)[wall_index] = 2;
			get_default_microenvironment()->add_dirichlet_node(n,ve);
			microenvironment(n)[virus_index] = 30;
			
		}
		else
		{   microenvironment(n)[wall_index] = 0;}
		
	}
	return; 
}

	
void setup_tissue_square( void )
{
	double C0 = 2;//2*1e4;//parameters.doubles("C0");//5127;
		
	static int virus_index = microenvironment.find_density_index( "virus");
	double MOI = parameters.doubles("MOI");
	
	double start_of_boundary = -979;//6913/2;
	double end_of_boundary = 979;//6913/2;
	double interval_lengths = (end_of_boundary-start_of_boundary)/(sqrt(C0)+1);
	double start_of_xboundary = start_of_boundary+interval_lengths;
		
	double length_x = microenvironment.mesh.bounding_box[3] - 
		microenvironment.mesh.bounding_box[0]; 
		
	double length_y = microenvironment.mesh.bounding_box[4] - 
		microenvironment.mesh.bounding_box[1]; 
		
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	double xgrid = 0.0;
	double ygrid = 0.0;
		
	
	//placing vien cells and normal cells in a hexagonal grid
	for( int i=0; i<sqrt(C0); i++ )
	{	
		if( i % 2==0 )
		{x = start_of_xboundary+i*interval_lengths;}//18;
		else
		{x = start_of_xboundary+1/2*interval_lengths+i*interval_lengths;}
		
		for( int j=0; j<sqrt(C0); j++ )
		{
			y = start_of_boundary+j*interval_lengths;//18;
			if(UniformRandom() < 0.3679)
			{pCell = create_cell( cancer_cell );}
		    else if(UniformRandom()<0.36) 
			{pCell = create_cell( CTL_cell );}
		else
		{pCell = create_cell( TH_cell );}
			double x = microenvironment.mesh.bounding_box[0] + UniformRandom() * length_x; 
			double y = microenvironment.mesh.bounding_box[1] + UniformRandom() * length_y; 
		
			pCell->assign_position( x , y , 0.0 );
			pCell->phenotype.molecular.internalized_total_substrates[virus_index] = 0; 
			double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
		    pCell->custom_data.variables[bias_phenotype_index].value =  UniformRandom();
		}
		
	}

	return;
	
}

void setup_tissue_circle( void )
{
	double C0 = 9740;//2*1e4;//parameters.doubles("C0");//5127;
		
	double Radius = 1144.4;
	
	double start_of_boundary = -Radius;//6913/2;
	double end_of_boundary = Radius;//6913/2;
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	//placing vien cells and normal cells in a hexagonal grid
	for( int i=0; i<C0; i++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		
		if(UniformRandom() < 0.9)
			{pCell = create_cell( cancer_cell );}
		    else if(UniformRandom()<0.95) 
			{pCell = create_cell( CTL_cell );}
		else
		{pCell = create_cell( TH_cell );}
		
		pCell->assign_position( x , y , 0.0 );
		double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
		pCell->custom_data.variables[bias_phenotype_index].value = UniformRandom();
		
		
	}

	return;
	
}

void setup_tissue_circle_immune( void )
{
	double C0 = 9740;//2*1e4;//parameters.doubles("C0");//5127;
		
	double Radius = 1144.4;
	
	double start_of_boundary = -Radius;//6913/2;
	double end_of_boundary = Radius;//6913/2;
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
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
	
    x = -115.9499;
	y = -349.9367;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	x = -292.5135;
	y = -91.2428;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	x = 134.3949;
	y = 347.2437;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	x = -170.1507;
	y = -132.5528;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	x = -195.5453;
	y = -5.9414;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  114.0970;
	y = -561.2747;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	//UNPROLIFERATIVE CTL CELLS
	x = 27.9074;
	y = 119.1812;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 225.7304;
	y = -728.5426;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -258.8621;
	y =  265.2323;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	x = 141.7866;
	y = -321.7799;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
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
	
	double GBM_NO = C0 - 61;
	//placing vien cells and normal cells in a hexagonal grid
	for( int i=0; i<GBM_NO; i++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		
		pCell = create_cell( cancer_cell );
				
		pCell->assign_position( x , y , 0.0 );
		
		
	}

	return;
	
}


void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
	//tumour cell proliferation
	if(pCell->type ==2)
	{	
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		double U = all_cells->size();
		if( pCell->state.simple_pressure>0.0241*15)// ccr (from macklins paper on DCIS)
		{
			pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;
		}
		else
		{
			pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0008648125;//0.0000064812*(1-U/512720);
		//std::cout<<pCell->state.simple_pressure<<std::endl;
		}
	}
	
	//cell movement
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	if( wall_amount>9 )
	{
		//std::cout<<"hit the boundary"<<std::endl;
		pCell->phenotype.motility.is_motile = false;
		pCell->is_movable = false;
	}
		
	//cell infection
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	
	double n = pCell->custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	
	double ubase = 0.000615;//uptake rate
	double Vvoxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	double nstar = 25;//infection threshol
	double nubase = 0.09;//repliation rate
	double alp = 1700;//virus burst number
	double sv = 10;//secretion rate
	double pstar = 10;
	
	static double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
	double bias_phenotype = pCell->custom_data.variables[bias_phenotype_index].value;
		
	static double uindex  = pCell->custom_data.find_variable_index( "uval" );
	double u = pCell->custom_data.variables[uindex].value; 
	/*if( u<0 )
	{
		pCell->custom_data.variables[uindex].value = NormalRandom(0.0015,0.0005); 
		u = pCell->custom_data.variables[uindex].value;
		if(u<0)
		{
			pCell->custom_data.variables[uindex].value = 1e-5; 
			u = pCell->custom_data.variables[uindex].value;
		}
		
	}*/
	

	static double nuindex  = pCell->custom_data.find_variable_index( "nuval" );
	double nu = pCell->custom_data.variables[nuindex].value; 
	if( nu<0 )
	{
		pCell->custom_data.variables[uindex].value = NormalRandom(0.0002,0.01); 
		nu = pCell->custom_data.variables[uindex].value;
		if(nu<0)
		{
			pCell->custom_data.variables[uindex].value = 1e-6; 
			nu = pCell->custom_data.variables[uindex].value;
		}
		
	}
	
	if( u<0 )
	{
		pCell->custom_data.variables[uindex].value = NormalRandom(0.00015,0.05); 
		u = pCell->custom_data.variables[uindex].value;
		if(u<1e-5)
		{
			pCell->custom_data.variables[uindex].value = 1e-10; 
			u = pCell->custom_data.variables[uindex].value;
		}
		
	}
	
	
	pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u;
	if( pCell->phenotype.death.dead == false & bias_phenotype>0.01)// cell not dead	
	{	
	
		if( n < nstar )
		{
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u;
			pCell->custom_data.variables[intracellular_virus_index].value = n+dt*(u*p*Vvoxel+nu*n*(1-n/alp)); 
			if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
			{std::cout<<"NEGATIVE HELP1"<<u*p*Vvoxel+nu*n*(1-n/alp)<<n<<" "<<alp<<" "<<nstar<<" "<<dt<<std::endl;}
		}
		else if(n< alp-100) //no more uptake but replication
		{
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->custom_data.variables[intracellular_virus_index].value = n+dt*(nu*n*(1-n/alp)); 
		}			
		else if( n > alp-100) //cell dies
		{
			pCell->start_death( apoptosis_model_index );	
		}
		else
		{std::cout<<"something wrong "<<n<<std::endl;}
	}
	else if( pCell->phenotype.death.dead==true && pCell->custom_data.variables[intracellular_virus_index].value>0)//dead cell secretes virus
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
	}
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
	double p_max = 3000;
	
	static double  V_0 = parameters.doubles("no_of_viruses_in_initial_injection")/221; //CIRCLE - no of vein cells
	
	double n_I = pCell->custom_data.variables[intracellular_virus_index].value;
	
	if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			if( n_I>1)//n_I > 1 )
			{
				double p_min = 1;
				double p_max = 2000;//V_0;
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
		//if(pCell-> phenotype.cycle.data.current_phase_index==0)
			//{ 
					output[0] = "orange";
					output[1] = "orange";
					output[2] = "coral";
					output[3] = "coral";
					return output; 
			/*}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	

					output[0] = "darkred";
					output[1] = "darkred";
					output[2] = "firebrick";
					output[3] = "firebrick";
			}*/		
			//}
			}
		}
		if( pCell->type == 3 && pCell->phenotype.death.dead==false)
		{
			if( n_I>1)//n_I > 1 )
			{
				double p_min = 1;
				double p_max = 1600;//V_0;
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
			//if(pCell-> phenotype.cycle.data.current_phase_index==0)
			//{   
			    output[0] = "aquamarine";
				output[1] = "lightsteelblue";
				output[2] = "lightskyblue";
				output[3] = "aquamarine";
			/*}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	
		
				output[0] = "darkslateblue";
				output[1] = "darkblue";
				output[2] = "darkblue";
				output[3] = "aquamarine";
			}*/
			}
		}
		if( pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if( n_I>1)//n_I > 1 )
			{
				double p_min = 1;
				double p_max = 1600;//V_0;
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

