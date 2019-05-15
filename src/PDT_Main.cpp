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
#define SIMID						1
#define CONTROLFILE					2
#define DAMASKSPECTRALFILENAME		3

void info() {
	//MK::provide definitions and other pieces of information relevant to the user
	cout << "Starting up DAMASK PDT v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << "\n";
	if ( VERSION_BETASTAGE == 1 )
		cout << "DAMASK PDT betastage" << "\n";
	cout << "Utilizing Cartesian global coordinate system, right-handed, x,y,z!" << "\n";
	Performance::displaySettings();
}

//polite information
void reporting( const string what, const int rank, const unsigned int thr, const bool verboseall ) {
	if (verboseall == false) { //not all output only master
		if (rank != MASTER)
			return; //non master do not verbose
	}
	//master or all remain
	cout << "VERBOSE::" << rank << "," << thr << " " << what << "\n";
}

//firing warnings but staying alive
void complaining( const string what, const int rank, const unsigned int thr ) {
	cerr << "WARNING::" << rank << "," << thr << " " << what << "\n";
}

//stopping in case of fatal errors
void stopping( const string what, const int rank, const unsigned int thr ) {
	cerr << "ERROR::" << rank << "," << thr << " " << what << "\n";
}


bool init_tasks_and_settings(  int pargc, char** pargv )
{
	//precision and format of console I/O
	cout << setprecision(18) << scientific;

	string croak = "";
	string s = "";

	Settings::SimID = stoul( pargv[SIMID] );
	try {
		Settings::readUserInput( pargv[CONTROLFILE] );
	}
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

	//##MK::optional overwriting of Settings::SpectralOutFilename to use same CONTROLFILE for multiple spectral out files
	if ( pargc > 3 ) {
		Settings::SpectralOutFilename = pargv[DAMASKSPECTRALFILENAME];
	}

	s = "Input is valid under SimulationID " + to_string(Settings::SimID) + " input read from " + Settings::SpectralOutFilename;
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
	double tic = MPI_Wtime();

	worker->set_myrank(rank);
	worker->set_nranks(nranks);
	worker->init_mpi_derivedtypes();
	worker->parse_taskspecific_opts();

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot(); //tictoc.get_memoryconsumption();
	worker->tictoc.prof_elpsdtime_and_mem( "ParameterSetup", APT_UTL, APT_IS_SEQ, mm, tic, toc, 0);

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
		unsigned int li_last = (worker->head.loadcasesmeta.size()==1) ?
				worker->head.loadcasesmeta.at(loadcase).nincr + 1 : worker->head.loadcasesmeta.at(loadcase).nincr;

		for ( unsigned int li = 0; li < li_last; ++li ) {
			if ( incr == targetincr ) {
				//for ( unsigned int incr = Settings::IncrementFirst; incr <= Settings::IncrementLast; incr = incr + Settings::IncrementOffset) {
				//still possible that no dump for this converged increment exists
				if ( worker->incr2rank.find(incr) != worker->incr2rank.end() ) {
					//so if an entry exists there should a worker assigned to process the dataset
					if ( worker->incr2rank.find(incr)->second == worker->get_myrank() ) { //I check whether there is a cherry to eat for me

cout << "li/incr\t\t" << li << "\t\t" << incr << "\n";

						double incrtic = MPI_Wtime();

						//parallelizes loading of simulation data at process scale
						//clarifying integration point grid configuration
						worker->thisincrement = incr;
						worker->thiswrittenincrement = static_cast<unsigned int>(worker->incr2wincr.find(incr)->second);

						worker->map_meshelements2threads();

						worker->grid_initial_configuration();

						worker->specout_read_heavydata2();

						//##MK::check status
						bool status_of_heavydata = true; //####MK::worker->specout_check_heavydata2(); //MK::true if okay, false otherwise
						if ( status_of_heavydata == true ) {
							auto it = worker->incr2healthy.find(incr);
							if ( it != worker->incr2healthy.end() )
								it->second = true;
							cout << "Heavy data are as expected!" << "\n";

							worker->pp3rve1withguard.owner = worker; //MK::equip specific instance of specOutIncr classes with address to calling object worker
							worker->pp3rve27.owner = worker;
							worker->grains.owner = worker;

							//##MK::changed for Matthew, process first
							if ( Settings::AnalyzeAddStrainTensors == true ) {
								worker->analyze_addStrainTensors_mp();
								worker->analyze_addRVEAverages2( loadcase, li, incr );
							}

							cout << "StrainTensors analyzed" << "\n";

/*
							//##MK::computation of Cauchy either based on averaged F & P or better per mp then average
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
*/

							//##MK::Markus your Paper13 DAMASKPDT_VoroComposer_FINAL

							worker->analyze_ipgrid_displacements();

							if(Settings::VisIPGridWithTextureID > 0) {
								worker->write_ipgrid_textureid();
							}

							if ( Settings::AnalyzeStateVarSpatialDistr == true ) {
								worker->bounded_volume_hierarchy();

								//if ( Settings::StateVarThrshldModel == E_DISORIANGLE ) {
								//	worker->analyze_svar_closestuip_disoriangle();
								//}
								//else {
								if ( Settings::AnalyzeReconstructGrains == true) {
									//so either or method Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION
									//Settings::StateVarThrshldModel == E_VORONOITESSELLATION

									worker->analyze_identify_grains();

									worker->analyze_boxup_grains();

									if ( Settings::StateVarThrshldModel == E_DISORIANGLE ||
											Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION ||
												Settings::StateVarThrshldModel == E_VORONOITESSELLATION ) {

										worker->analyze_build_grains();

										switch(Settings::StateVarThrshldModel)
										{
											case E_DISORIANGLE:
												worker->analyze_svar_closestuip_disoriangle();
												break;
											case E_SIGNEDDISTANCEFUNCTION:
												worker->analyze_svar_grainbased_sdf();
												break;
											case E_VORONOITESSELLATION:
												worker->analyze_svar_grainbased_cgeom();
												break;
											default:
												break;
										}

										//if ( Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION ) {
										//	worker->analyze_svar_grainbased_sdf();
										//}
										//if ( Settings::StateVarThrshldModel == E_VORONOITESSELLATION ) {
										//	worker->analyze_svar_grainbased_cgeom();
										//}
									}
								}
								//}
							}

/*
							worker->analyze_addRVEAverages( loadcase, li, incr );
							//RVE averages are required for correct mapping of periodic images of ip

							//##MK::DEBUG, currently go out after having data for the flow curve
							//i.e. either compute flow curve or more advanced point-wise stuff...

							worker->analyze_ipgrid_displacements();

							if(Settings::VisIPGridWithTextureID > 0) {
								worker->write_ipgrid_textureid();
							}

							//work through individual analysis tasks for this increment

							//##MK::additional RVE averaged tensorial quantities?

							if ( Settings::AnalyzeStateVarSpatialDistr == true ||
									Settings::AnalyzeReconstructGrains == true) {
								worker->bounded_volume_hierarchy();


								if ( Settings::AnalyzeStateVarSpatialDistr == true ) {
									if ( Settings::StateVarThrshldModel == E_DISORIANGLE )
										worker->analyze_svar_closestuip_disoriangle();
								}

								if ( Settings::AnalyzeReconstructGrains == true ||
										(Settings::AnalyzeStateVarSpatialDistr == true &&
												Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION) ) {
									worker->analyze_identify_grains();
									//worker->debug_signed_distance_function();

									if (Settings::AnalyzeStateVarSpatialDistr == true &&
											Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION) {

										worker->analyze_boxup_grains();
									}

									if (Settings::AnalyzeStateVarSpatialDistr == true &&
											Settings::StateVarThrshldModel == E_SIGNEDDISTANCEFUNCTION)
										worker->analyze_svar_grainbased_sdf();

									//worker->analyze_reconstruct_grains(); //MK::deprecated and inconsistent as the similar in its approach implemented in addGrainID.py
								}
							}
*/
							//##MK::add further analysis task implementations if desired
							//else flow curve data are at RVE-scale, these however, were already computed, so

						} //done processing healthy heavydata
						else { //##MK::dataset is incomplete, likely because DAMASK_spectral convergence issues
							auto it = worker->incr2healthy.find(incr);
							if ( it != worker->incr2healthy.end() )
								it->second = false;

							cout << "Heavy data are faulty!" << "\n";
						}

						//release resources in this MPI process for reinit and reusage within next ranks increment to process
						double mtic = MPI_Wtime();

						worker->free_increment_heavydata();
//cout << "free heavydata" << "\n";
						worker->free_increment_bvh_xyzm2();
//cout << "free bvh xyzm2" << "\n";
						worker->free_increment_bvh_p3dm1();
//cout << "free bvh p3dm1" << "\n";
						worker->free_increment_grains();
//cout << "free grains" << "\n";
						worker->free_increment_sdf();
//cout << "free sdf" << "\n";
						//worker->report_rveshapes();

						double mtoc = MPI_Wtime();
						memsnapshot mm = worker->tictoc.get_memoryconsumption(); //document that memory consumption was reduced again
						worker->tictoc.prof_elpsdtime_and_mem( "ResetMemory", APT_UTL, APT_IS_SEQ, mm, mtic, mtoc, incr );

						double incrtoc = MPI_Wtime();
						string mess = "Processed global increment " + to_string(incr) + " took " + to_string(incrtoc-incrtic) + " seconds";
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
	//##MK::
	//worker->report_flowcurve2();

	//##MK::for Matthew Kasemer
	double utic = MPI_Wtime();

	worker->report_flowcurve3( "BasedOnFT", EPS_RIGHTCAUCHYGREEN_LN_FT );
	worker->report_flowcurve3( "BasedOnFP", EPS_RIGHTCAUCHYGREEN_LN_FP );

	double utoc = MPI_Wtime();
	mm = memsnapshot(); //worker->tictoc.get_memoryconsumption();
	worker->tictoc.prof_elpsdtime_and_mem( "ReportFlowcurves", APT_UTL, APT_IS_SEQ, mm, utic, utoc, 0);

	//##MK::DEBUG
	worker->tictoc.spit_profiling( Settings::SimID, worker->get_myrank() );

	//delete MPI working instance
	if ( worker != NULL ) {
		delete worker;
		worker = NULL;
	}
}


int main(int argc, char** argv)
{
//INTERNAL PROFILING
	//timeval time;
	//gettimeofday(&time, NULL);
	//double globaltic = time.tv_sec + time.tv.usec / 1000000.0;
	double globaltic = omp_get_wtime();

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

	//gettimeofday(&time, NULL);
	//double globaltoc = time.tv_sec + time.tv_usec / 1000000.0;
	double globaltoc = omp_get_wtime();
	cout << "DAMASKPDT took in total " << (globaltoc - globaltic) << " seconds" << endl;

	return 0;
}
