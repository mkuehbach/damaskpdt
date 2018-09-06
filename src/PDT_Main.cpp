/*
	Copyright Max-Planck-Institut f\"uer Eisenforschung, GmbH, D\"usseldorf
	Programmer: Markus K\"uhbach, 2017

	damaskpdt --- the DAMASK parallelized property distribution tracker --- is an MPI/OpenMP/SIMD-
	parallelized tool for datamining the spatial distribution of state variable values at interfaces
	which evolve during straining in 3d full-field crystal plasticity simulations.

	It supplements the python post-processing routines of the DAMASK crystal plasticity solver
	package <https://damask.mpie.de/>. The novelty is the rigorous parallel and locality-aware
	spatial-distributed storing of snapshot data to make the processing of terabyte-scale
	crystal plasticity simulation snapshot data efficient, scalable.

	This enables to study and quantify distributional trends, which other techniques and
	analysis workflows have so far not practically been able to process.

	The source code was developed by Markus K\"uhbach during his postdoc with the Max-Planck-Institut f\"ur
	Eisenforschung GmbH in D\"usseldorf for the purpose of quantifying the potential accumulation of dislocation
	density, stresses, and strain in front of interfaces. The model implements a distance based approach
	m.kuehbach at mpie.de

	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the project "Consistent physically-based modeling of dynamic recrystallization
	under hot working conditions" project number 315419526.

	This file is part of damaskpdt.

	damaskpdt is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	damask_pdt is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with damask_pdt.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "PDT_OptionParser.h"


//parameter handed over
#define SIMID				1
#define CONTROLFILE			2

void info() {
	//MK::provide definitions and other pieces of information relevant to the user
	cout << "Starting up DAMASK PDT v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << endl;
	cout << "Utilizing Cartesian global coordinate system, right-handed, x,y,z!" << endl;
	Performance::displaySettings();
}

//polite information
void reporting( const string what, const int rank, const unsigned int thr, const bool verboseall ) {
	if (verboseall == false) { //not all output only master
		if (rank != MASTER)
			return; //non master do not verbose
	}
	//master or all remain
	cout << "VERBOSE::" << rank << "," << thr << " " << what << endl;
}

//firing warnings but staying alive
void complaining( const string what, const int rank, const unsigned int thr ) {
	cout << "WARNING::" << rank << "," << thr << " " << what << endl;
}

//stopping in case of fatal errors
void stopping( const string what, const int rank, const unsigned int thr ) {
	cout << "ERROR::" << rank << "," << thr << " " << what << endl;
}


bool init_tasks_and_settings(  int pargc, char** pargv )
{
	//precision and format of console I/O
	cout << setprecision(18) << scientific;

	string croak = "";
	string s = "";

	Settings::SimID = stoul( pargv[SIMID] );
	try { Settings::readUserInput(pargv[CONTROLFILE]); }
	catch (exception& croak) {
		stopping( "Unable to parse control file!", MASTER, MASTER);
		stopping( croak.what(), MASTER, MASTER );
		return false;
	}
	if ( Settings::checkUserInput() == false ) {
		croak = "Control file settings failed the validity check!";
		stopping( croak, MASTER, MASTER );
		return false;
	}

	Settings::displaySettings();

	s = "Input is valid under SimulationID " + to_string(Settings::SimID);
	reporting(s, MASTER, MASTER, false);
	return true;
}


void AnalysisDebug( int nranks, int rank, int pargc, char** pargv )
{
	//MK::mind that we are in an MPI parallel environment identified as MPI_COMM_WORLD
	specOutIncr* worker = NULL;
	try { worker = new specOutIncr; }
	catch (bad_alloc &croak) {
		stopping( croak.what(), rank, 0 );
		return;
	}
	worker->set_myrank(rank);
	worker->set_nranks(nranks);
	worker->init_mpi_derivedtypes();

//read header and data structure of the spectralOut binary with all processes
	//##MK::can later be improved by reading by only one process and broadcasting...

	worker->specout_read_header();
	worker->specout_read_structure_homogenization();
	worker->specout_read_structure_crystallite();
	worker->specout_read_structure_constitutive();

	//MK::we have to take into account that not every increment has rawdata written to spectralOutfile and that frequency is load-case-dependent
	int localhealth = worker->map_increments2ranks(); //MK::if returning 1 rank remains healthy, else 0

	MPI_Barrier(MPI_COMM_WORLD);
	int globalhealth = 0;
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if ( globalhealth != worker->get_nranks() ) {
		stopping("Either work was not properly distributable or not all processes able to read the heavydata", worker->get_myrank(), 0);
	}

	//adjust analysis range if file contains less data
	if ( worker->head.get_lincr() < Settings::IncrementLast ) {
		Settings::IncrementLast = worker->head.get_lincr();
		string warn = "Narrowing analysis range to " + to_string(Settings::IncrementLast) + " because file does not contain more increments!";
		complaining( warn, worker->get_myrank(), 0);
	}

	//process in parallel different increments
	unsigned int incr = worker->head.sincr;
	unsigned int targetincr = Settings::IncrementFirst;
	for ( unsigned int loadcase = 0; loadcase < worker->head.loadcases; ++loadcase) {
		//unsigned int loadcase_dump_freq = worker->head.loadcasesmeta.at(loadcase).freq;
		for ( unsigned int li = 0; li < worker->head.loadcasesmeta.at(loadcase).nincr; ++li ) {
			if ( incr == targetincr ) {
				//for ( unsigned int incr = Settings::IncrementFirst; incr <= Settings::IncrementLast; incr = incr + Settings::IncrementOffset) {
				//still possible that no dump for this converged increment exists
				if ( worker->incr2rank.find(incr) != worker->incr2rank.end() ) {
					//so if an entry exists there should a worker assigned to process the dataset
					if ( worker->incr2rank.find(incr)->second == worker->get_myrank() ) { //I check whether there is a cherry to eat for me

						//parallelizes loading of simulation data and clarifying integration point grid
						worker->thisincrement = incr;
						worker->thiswrittenincrement = static_cast<unsigned int>(worker->incr2wincr.find(incr)->second);

						worker->map_meshelements2threads();

						worker->grid_initial_configuration();

						worker->specout_read_heavydata2();

						if ( worker->specout_check_heavydata2() == true ) {
							cout << "Heavy data are as expected!" << endl;
						}
						else {
							//##MK::implement error management
							cout << "Heavy data are faulty!" << endl;
						}

						worker->specout_auxiliarytasks();

						worker->pp3rve1withguard.owner = worker; //MK::equip specific instance of specOutIncr classes with address to calling object worker
						worker->pp3rve27.owner = worker;
						worker->grains.owner = worker;

						worker->analyze_addRVEAverages( loadcase, li, incr );
						//RVE averages are required for correct mapping of periodic images of integration points

						//##MK::DEBUG, currently go out after having data for the flow curve
						//i.e. either compute flow curve or more advanced point-wise stuff...

						worker->analyze_ipgrid_displacements();

						if(Settings::VisIPGridWithTextureID > 0) {
							worker->write_ipgrid_textureid();
						}

						//work through individual analysis tasks for this increment
						if ( Settings::AnalyzeAddStrainTensors == true ) {
							worker->analyze_addStrainTensors();
						}

						if ( Settings::AnalyzeAddCauchy == true ) {
							worker->analyze_addCauchy();
						}

						//vonMises practically of interest only when at least strain or cauchy computed
						//but its output is optional
						if (Settings::AnalyzeAddStrainTensors == true && Settings::AnalyzeAddCauchy == true) {
							if ( Settings::AnalyzeAddVonMises == true ) {
								worker->analyze_addVonMises();
							}
						}

						if ( Settings::AnalyzeStateVarSpatialDistr == true ||
								Settings::AnalyzeReconstructGrains == true) {
							worker->bounded_volume_hierarchy();

							if ( Settings::AnalyzeStateVarSpatialDistr == true ) {
								if ( Settings::StateVarThrshldModel == E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD )
									worker->analyze_svar_closestuip_disoriangle();
							}

							if ( Settings::AnalyzeReconstructGrains == true ||
									(Settings::AnalyzeStateVarSpatialDistr == true &&
											Settings::StateVarThrshldModel == E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION) ) {
								worker->analyze_identify_grains();
								//worker->debug_signed_distance_function();

								if (Settings::AnalyzeStateVarSpatialDistr == true &&
										Settings::StateVarThrshldModel == E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION) {

									worker->analyze_mesh_grains();
								}

								if (Settings::AnalyzeStateVarSpatialDistr == true &&
										Settings::StateVarThrshldModel == E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION)
									worker->analyze_svar_grainbased_sdf();

								//worker->analyze_reconstruct_grains(); //MK::deprecated and inconsistent as the similar in its approach implemented in addGrainID.py
							}
						}

						//##MK::add further analysis task implementations if desired
						//else flow curve data are at RVE-scale, these however, were already computed, so

						//release resources in this MPI process for reinit and reusage within next ranks increment to process
						worker->free_increment_heavydata();
						worker->free_increment_bvh_xyzm2();
						worker->free_increment_bvh_p3dm1();
						worker->free_increment_grains();
						worker->free_increment_sdf();

						worker->report_rveshapes();

						string mess = "Processed global increment " + to_string(incr);
						reporting( mess, worker->get_myrank(), 0, true );
					}
				}

				targetincr += Settings::IncrementOffset;
				if ( targetincr > Settings::IncrementLast )
					targetincr = Settings::IncrementLast; //thus enforcing a stop criterion, as incr keeps increasing...
			}

			incr++;
			//##MK::targetincr

		} //next increment of current loadcase
	} //next loadcase

	//MK::wait for all MPI processes to compute increment data before doing I/O of multiple increments!
	MPI_Barrier(MPI_COMM_WORLD);

	//perform results I/O
	worker->report_flowcurve2();

	//##MK::DEBUG
	worker->spit_profiling();

	//delete MPI working instance
	if ( worker != NULL ) {
		delete worker;
		worker = NULL;
	}
}

void Debugging()
 {
	 vector<quat> oris;
	 //MK::already passive to interpret quaternion components!
	 oris.push_back(quat(8.897931877474624107e-01, 1.007390959234611971e-01,-2.509506738574025508e-01,-3.676186568731902304e-01));
	 oris.push_back(quat(9.476387593060630055e-01,-1.934032399262330115e-01,2.537172454462384708e-01,-1.426632433652256192e-02));
	 oris.push_back(quat(9.087577692838431087e-01,3.213676556156696007e-01,1.424273189308490339e-01,2.249368922887307332e-01));
	 oris.push_back(quat(9.584333869365586622e-01,2.601558784926337187e-01,1.169568934324672349e-01,-6.741421932303652727e-03));
	 oris.push_back(quat(9.005609600548302174e-01,-2.545728269126454335e-01,-1.537348448602691986e-01,3.170934097369666471e-01));
	 oris.push_back(quat(9.310051492271317342e-01,-2.300408556699746754e-01,-5.543603878545275837e-02,-2.779162867465125863e-01));
	 oris.push_back(quat(9.100665056774119854e-01,-1.037114195853210685e-01,3.200364774697847525e-01,2.420734388164769713e-01));
	 oris.push_back(quat(9.038947261999651372e-01,1.235016564242449322e-01,3.066708332955571659e-01,2.714307735209229855e-01));
	 oris.push_back(quat(9.116240679404760128e-01,4.633192292250594771e-03,-3.802997584418089105e-01,-1.558595072812154758e-01));
	 oris.push_back(quat(9.029047338058159289e-01,1.365429353072289964e-01,-2.955059867668757834e-01,2.807049701620363824e-01));
	 for (unsigned int i = 0; i < 10; ++i) {
		 real_ori o1[4] = {oris[i].q0,oris[i].q1,oris[i].q2,oris[i].q3};
		 //active2passive( o1 );
		 for(unsigned int j = 0; j < 10; ++j) {
			 if ( i > j ) { //i != j
				 real_ori o2[4] = {oris[j].q0,oris[j].q1,oris[j].q2,oris[j].q3};
				 //active2passive( o2 );
				 real_ori disori = disorientation_angle_fcc( o1, o2);
				 //if ( disori <= DEGREE2RADIANT(15.0) )
					 cout << i << "\t\t" << j << "\t\t" << disori << "\t\t";
					 cout << o1[0] << ";" << o1[1] << ";" << o1[2] << ";" << o1[3];
					 cout << "---" << o2[0] << ";" << o2[1] << ";" << o2[2] << ";" << o2[3] << endl;
			 }
		 }
cout << endl;
	 }
 }


int main(int argc, char** argv)
{
	//Debugging();
	//return 0;

//SETUP
	info();
	if ( argc < 2 ) {
		stopping( "We need at least a simulation id <unsigned long> and a XML control file <*.xml>!", MASTER, MASTER);
		return 0;
	}
	init_tasks_and_settings( argc, argv );

//go MPI parallel with sufficient hybrid MPI/OpenMP functionality
	int nr, r, supportlevel_provided;
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	MPI_Init_thread( &argc, &argv, supportlevel_desired, &supportlevel_provided);
	if ( supportlevel_provided < supportlevel_desired ) {
		stopping( "Insufficient threading capabilities of the MPI library!", MASTER, MASTER);
		return 0;
	}
	else { //sufficiently supportive MPI loaded
		reporting("Initialized MPI_COMM_WORLD", MASTER, MASTER, false);
		MPI_Comm_size(MPI_COMM_WORLD, &nr);
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		reporting( "MPI firing up rank " + to_string(r), r, 0, true );
		reporting( "MPI ranks in total " + to_string(nr), r, 0, true );
	}

//START THE ANALYSIS
	switch ( Settings::AnalysisMode )
	{
		case E_DEBUG :
			AnalysisDebug( nr, r, argc, argv );
			break;
		default :
			reporting("Nothing to do for me :(!", r, 0, false);
	}

//exit parallel environment and return control to operating system
	reporting("Stopping MPI_COMM_WORLD", r, 0, false);
	MPI_Finalize();
	return 0;
}
