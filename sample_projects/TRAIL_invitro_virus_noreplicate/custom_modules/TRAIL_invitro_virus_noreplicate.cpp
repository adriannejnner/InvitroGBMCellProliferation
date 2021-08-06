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

#include "./TRAIL_invitro_virus_noreplicate.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>

void create_cell_types( void )
{
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	double C0 = parameters.doubles("C0")/100;
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00064812*(1-C0/6515.5);//0.00064812*(1-C0/5127);
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	
 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 
	
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "intracellular_TRAIL_amount", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "replication_start_time", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "infection_time_start", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "infection_time_length", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "time_of_injection", "dimensionless", 0.0 );
	
	// turn off secretion from these cells (oxygen and virus)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	cell_defaults.phenotype.secretion.secretion_rates[1] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[1] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[1] = 2600; 
	
	cell_defaults.phenotype.geometry.radius = 10.25;
	
	// update cell and phenotype based on virus dynamics only
	//cell_defaults.functions.update_phenotype = nothing_function;//virus_dynamics_TRAIL_MODEL3; 
	
	cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3;//ECM_movement;//virus_dynamics; 
	
	cell_defaults.phenotype.motility.is_motile = true; 
		
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
	
	return; 
}	
void setup_tissue_square( void )
{
		
	double C0 = parameters.doubles("C0")/100;//5127;

	double MOI = parameters.doubles("MOI");
	
	double start_of_boundary = -979.28;//6913/2;
	double end_of_boundary = 979.28;//6913/2;
	double interval_lengths = (end_of_boundary-start_of_boundary)/(sqrt(C0)+1);
	double start_of_xboundary = start_of_boundary+interval_lengths;
		
	
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	double xgrid = 0.0;
	double ygrid = 0.0;
		
	std::default_random_engine generator;
	std::poisson_distribution<int> distribution(MOI);
	
	
	//placing vien cells and normal cells in a hexagonal grid
	for( int i=0; i<sqrt(C0)-1; i++ )
	{	
		if( i % 2==0 )
		{x = start_of_xboundary+i*interval_lengths;}//18;
		else
		{x = start_of_xboundary+1/2*interval_lengths+i*interval_lengths;}
		
		for( int j=0; j<sqrt(C0)+1; j++ )
		{
			y = start_of_boundary+j*interval_lengths+0.5*interval_lengths;//18;
			pCell = create_cell(  );
			//std::cout<<number<<std::endl;
			int number = distribution(generator);
			
			if(number>125)
			{number = 125;}
			int intracellular_TRAIL_index = pCell->custom_data.find_variable_index( "intracellular_TRAIL_amount" );
			pCell->custom_data.variables[intracellular_TRAIL_index].value = number;
			int intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
			pCell->custom_data.variables[intracellular_virus_index].value = number;
		
			pCell->assign_position( x , y , 0.0 );
		}
		
	}

	
	
	
	

	return;
	
}
void nothing_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	return;
	
}
void virus_dynamics_TRAIL_MODEL3( Cell* pCell, Phenotype& phenotype, double dt )
{
	// PROLIFERATION ^^
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
    //static int pressure_store_index =pCell->phenotype.death.find_death_model_index( "pressure_store" );	
	//pCell->custom_data.variables[pressure_store_index].value = pCell->state.simple_pressure;
		double U = all_cells->size();
	pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00064812*(1-U/6515.5);	
	
	
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
	static double viral_replication_rate = parameters.doubles("viral_replication_rate");	
	static double TRAIL_secretion_rate = parameters.doubles("TRAIL_secretion_rate");
	static double TRAIL_generation_rate = parameters.doubles("TRAIL_generation_rate");
	static double TRAIL_killing_level = parameters.doubles("TRAIL_killing_level");
	static double M = parameters.doubles("M");
	
	double n_T = infection_density_capacity;//capacity of cell
	double c_T = TRAIL_generation_rate; // 0.1
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
	
	double n_I_zero = 0;
	if( n_I < n_T )
	{n_I_zero = n_I;}
	else
	{n_I_zero = n_T;}
	
	if( pCell->phenotype.death.dead == false && T_I > 0)// cell not dead	
	{			
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
	
	
	if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
	{
		std::cout<<"NEGATIVE HELP TRAIL"<<std::endl;
	}
	
		
	 if( T_E>0.0001 && pCell->phenotype.death.dead==false)
	{
		double M = 0.1;
		pCell->phenotype.death.rates[apoptosis_model_index] = (T_E*T_E/(M*M+T_E*T_E))*0.05;
	}
		
		
	
	return;
}

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	static int infection_density_capacity = parameters.doubles("infection_density_capacity");
	
	double p_min = 1;
	double p_max = 125;
	//static double  V_0 = parameters.doubles("no_of_viruses_in_initial_injection");
	//static double  V_0 = parameters.doubles("no_of_viruses_in_initial_injection")/1215; //TRIANGLE - no of vein cells
	static double  V_0 = parameters.doubles("no_of_viruses_in_initial_injection")/221; //CIRCLE - no of vein cells
	
	
	double n_I = pCell->custom_data.variables[intracellular_virus_index].value; // amount of virus
	
	if( pCell->phenotype.death.dead==false && pCell->type==0)
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->custom_data.variables[intracellular_virus_index].value-p_min) * 255.0 ); 
		char szTempString [128]; // ceates a character array that can store 128
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
		
		
		if(pCell->custom_data.variables[intracellular_virus_index].value>1)
		{
			
			int oncoprotein = (int) round( (1.0/(p_min-p_max)) * (pCell->custom_data.variables[intracellular_virus_index].value-p_max) * 230.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein,255); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			output[0].assign( szTempString );
			output[1].assign( szTempString );//="blue";
			output[2].assign( szTempString );
			output[3]="blue";
			
			/*int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->custom_data.variables[intracellular_virus_index].value-p_min) * 255.0 ); 
			char szTempString [128];
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
			output[0].assign( szTempString );
			output[1].assign( szTempString );
			output[2].assign( szTempString );
			output[3].assign( szTempString );
			//sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
			//output[2].assign( szTempString );*/

				
			return output;
		}
		else if( pCell->custom_data.variables[intracellular_virus_index].value <0)
		{
			std::cout<<"NEGATIVE"<<std::endl;
			std::vector< std::string > output( 4, "darkgrey" ); 
			return output;
		}
		else
		{
				output[0] = "rgb(255, 182, 193)";
				output[1] = "rgb(255, 182, 193)";
				output[2] = "rgb(255, 160, 122)";
				output[3] = "rgb(255, 160, 122)";
				return output; 
		}
	}
	
	if( pCell->type == 1 )
	{
		if( n_I > 0 )
		{
			double p_min = 1;
			double p_max = V_0;
			
			int oncoprotein = (int) round( (1.0/(p_min-p_max)) * (pCell->custom_data.variables[intracellular_virus_index].value-p_max) * 230.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", 255, oncoprotein, oncoprotein); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			output[0].assign( szTempString );
			output[1]="red";
			output[2].assign( szTempString );
			output[3]="red";
			
			return output;
		}
		else
		{
			    output[0] = "red";//"rgb(255,230,230)";
				output[1] = "red";
				output[2] = "red";//"rgb(255,230,230)";
				output[3] = "red";
				//output[0] = "saddlebrown";//rgb(255,230,230)";
				//output[1] = "red";
				//output[2] = "peru";//rgb(255,230,230)";
				//output[3] = "red";
				return output; 
				
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

