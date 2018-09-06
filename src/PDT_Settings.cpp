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


#include "PDT_Settings.h"

using namespace rapidxml;


E_ANALYSIS_MODE Settings::AnalysisMode = E_DEBUG;
//E_THRESHOLDING_MODE Settings::SpatialDistrThrshldModel = E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD;	//which is the decisive quantity against which one thresholds distance when quantifying the accumulation of state variable values
string Settings::SpectralOutFilename = "";
unsigned int Settings::IncrementFirst = 0;
unsigned int Settings::IncrementOffset = 1;
unsigned int Settings::IncrementLast = 0;

bool Settings::AnalyzeFlowCurve = false;
bool Settings::AnalyzeAddStrainTensors = false;
bool Settings::AnalyzeAddCauchy = false;
bool Settings::AnalyzeAddVonMises = false;
bool Settings::AnalyzeAddDisplacements = false;
bool Settings::AnalyzeStateVarSpatialDistr = false;
bool Settings::AnalyzeReconstructGrains = false;

unsigned int Settings::GrainReconModel = E_HIERARCHICAL_CLUSTERING;
unsigned int Settings::StateVarThrshldModel = E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD;
bool Settings::ThrshldAgainstDefConfig = false;

unsigned long Settings::SimID = 0;
real_xyz Settings::KernelRadius = 0.0;
real_xyz Settings::DBScanKernelRadius = 0.0; //eps
real_xyz Settings::PVTessKernelRadius = 0.0; //<2.0 min(ip/ip distance undeformed configuration)
real_xyz Settings::PVTessCubeVoxelEdgeLen = 0.0;
real_ori Settings::CriticalDisoriAngle = DEGREE2RADIANT(15.0);
real_ori Settings::CommDetectionDisoriAngleWght = 75.0; //S. Dancette et al., MSMSE 24 2016
real_xyz Settings::CommDetectionDistanceWght = 1.0;
//real_ori Settings::GrainReconLocalDisoriAngle = DEGREE2RADIANT(15.0);
real_xyz Settings::GrainReconInspectionRadius = 0.0;
long double Settings::GrainReconLouvainPrecision = LOUVAIN_DEFAULT_PRECISION;

//real_xyz Settings::PhysicalDomainSize = 1.0; //in meter

unsigned int Settings::PVTessGuardZoneEdgeLen = 3; 		//in multiples of generalized coordinate PVTessCubeVoxelEdgeLen

unsigned int Settings::VisIPGridWithTextureID = 0;
unsigned int Settings::VisIPGridWithGrainID = 0;
unsigned int Settings::VisIPGridWithPerImages = 0;
unsigned int Settings::VisHigherOrderNeighbor = 0;
unsigned int Settings::VisPeriodicReplicaIPs = 0;
unsigned int Settings::VisDBScanClusterIDs = 0;
unsigned int Settings::VisReconstructedGrain = 0;
unsigned int Settings::VisRVE27BVH = 0;
unsigned int Settings::VisGrainLocalVoxelGrid = 0;
unsigned int Settings::VisGrainFinalReconstruction = 0;

unsigned int Settings::VisGrainQuaternionClouds = 0;
double Settings::DebugDouble = 0.0;
unsigned int Settings::DebugUnsignedInt = 0;


void Settings::readUserInput(string filename) {

	//find the desired .xml file
	if ( 0 == filename.compare("") )
		filename = string("PDT.Input.xml");

	ifstream file( filename );
	if (file.fail()) {
		throw runtime_error(string("Unable to locate input file ") + filename);
	}

	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	std::string fnm(rootNode->name()); //convert const char* to string
	string key = "Parameters";
	if (0 != fnm.compare(key)) {
		throw runtime_error("Undefined parameters file!");
	}

//META DATA
	unsigned long mode = 0;
	if (0 != rootNode->first_node("AnalysisMode")) {
		mode = stoul( rootNode->first_node("AnalysisMode")->value());
	}
	if ( mode == E_DEBUG )
		AnalysisMode = E_DEBUG;

	if (0 != rootNode->first_node("SpectralOutFilename")) {
		SpectralOutFilename = rootNode->first_node("SpectralOutFilename")->value();
	}

	if (0 != rootNode->first_node("IncrementFirst")) {
		IncrementFirst = static_cast<unsigned int>(stoi(rootNode->first_node("IncrementFirst")->value()));
	}
	if (0 != rootNode->first_node("IncrementOffset")) {
		IncrementOffset = static_cast<unsigned int>(stoi(rootNode->first_node("IncrementOffset")->value()));
	}
	if (0 != rootNode->first_node("IncrementLast")) {
		IncrementLast = static_cast<unsigned int>(stoi(rootNode->first_node("IncrementLast")->value()));
	}

//DEFINE WHAT WE HAVE TO DO
	if (0 != rootNode->first_node("AnalyzeFlowCurve")) {
		if (  static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeFlowCurve")->value())) == 1 )
			AnalyzeFlowCurve = true;
	}
	if (0 != rootNode->first_node("AnalyzeAddStrainTensors")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeAddStrainTensors")->value())) == 1 )
			AnalyzeAddStrainTensors = true;
	}
	if (0 != rootNode->first_node("AnalyzeAddCauchy")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeAddCauchy")->value())) == 1 )
			AnalyzeAddCauchy = true;
	}
	if (0 != rootNode->first_node("AnalyzeAddVonMises")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeAddVonMises")->value())) == 1 )
			AnalyzeAddVonMises = true;
	}
	if (0 != rootNode->first_node("AnalyzeAddDisplacements")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeAddDisplacements")->value())) == 1 )
			AnalyzeAddDisplacements = true;
	}
	if (0 != rootNode->first_node("AnalyzeStateVarSpatialDistr")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeStateVarSpatialDistr")->value())) == 1 )
			AnalyzeStateVarSpatialDistr = true;
	}
	if (0 != rootNode->first_node("AnalyzeReconstructGrains")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("AnalyzeReconstructGrains")->value())) == 1 )
			AnalyzeReconstructGrains = true;
	}

	if (0 != rootNode->first_node("GrainReconModel")) {
		unsigned long what = static_cast<unsigned long>(stoul(rootNode->first_node("GrainReconModel")->value() ));
		if ( what == E_HIERARCHICAL_CLUSTERING )
			GrainReconModel = E_HIERARCHICAL_CLUSTERING;
		if ( what == E_INITIAL_TEXTUREID )
			GrainReconModel = E_INITIAL_TEXTUREID;
	}
	if (0 != rootNode->first_node("SpatialDistrThrshldModel")) {
		unsigned long what = static_cast<unsigned long>(stoul(rootNode->first_node("SpatialDistrThrshldModel")->value()));
		if ( what == E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD )
			StateVarThrshldModel = E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD;
		if ( what == E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION )
			StateVarThrshldModel = E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION;
	}

	if (0 != rootNode->first_node("ThrshldAgainstDefConfiguration")) {
		if ( static_cast<unsigned long>(stoul(rootNode->first_node("ThrshldAgainstDefConfiguration")->value())) == 1 )
			ThrshldAgainstDefConfig = true;
		else
			ThrshldAgainstDefConfig = false;
	}


//ANALYSES-SPECIFIC SETTINGS

	if (0 != rootNode->first_node("KernelRadius"))
#ifdef SINGLE_PRECISION
		KernelRadius = stof( rootNode->first_node("KernelRadius")->value() );
#else
		KernelRadius = stod( rootNode->first_node("KernelRadius")->value() );
#endif

	if (0 != rootNode->first_node("CriticalDisoriAngle"))
#ifdef SINGLE_PRECISION
		CriticalDisoriAngle = DEGREE2RADIANT(stof( rootNode->first_node("CriticalDisoriAngle")->value() ));
#else
		CriticalDisoriAngle = DEGREE2RADIANT(stod( rootNode->first_node("CriticalDisoriAngle")->value() ));
#endif

	if (0 != rootNode->first_node("CommDetectionDisoriAngleWeighting"))
#ifdef SINGLE_PRECISION
		CommDetectionDisoriAngleWght = stof( rootNode->first_node("CommDetectionDisoriAngleWeighting")->value() );
#else
		CommDetectionDisoriAngleWght = stod( rootNode->first_node("CommDetectionDisoriAngleWeighting")->value() );
#endif

		if (0 != rootNode->first_node("CommDetectionDistanceWeighting"))
	#ifdef SINGLE_PRECISION
			CommDetectionDistanceWght = stof( rootNode->first_node("CommDetectionDistanceWeighting")->value() );
	#else
			CommDetectionDistanceWght = stod( rootNode->first_node("CommDetectionDistanceWeighting")->value() );
	#endif

/*
	if (0 != rootNode->first_node("GrainReconLocalDisoriAngle"))
#ifdef SINGLE_PRECISION
		GrainReconLocalDisoriAngle = DEGREE2RADIANT(stof( rootNode->first_node("GrainReconLocalDisoriAngle")->value() ));
#else
		GrainReconLocalDisoriAngle = DEGREE2RADIANT(stod( rootNode->first_node("GrainReconLocalDisoriAngle")->value() ));
#endif
*/

	if (0 != rootNode->first_node("GrainReconInspectionRadius"))
#ifdef SINGLE_PRECISION
		GrainReconInspectionRadius = stof( rootNode->first_node("GrainReconInspectionRadius")->value() );
#else
		GrainReconInspectionRadius = stod( rootNode->first_node("GrainReconInspectionRadius")->value() );
#endif

	if (0 != rootNode->first_node("GrainReconLouvainPrecision"))
		GrainReconLouvainPrecision = static_cast<long double>(stod( rootNode->first_node("GrainReconLouvainPrecision")->value() ));

	if (0 != rootNode->first_node("DBScanKernelRadius"))
#ifdef SINGLE_PRECISION
		DBScanKernelRadius = stof( rootNode->first_node("DBScanKernelRadius")->value() );
#else
		DBScanKernelRadius = stod( rootNode->first_node("DBScanKernelRadius")->value() );
#endif

	if (0 != rootNode->first_node("PVTessKernelRadius"))
#ifdef SINGLE_PRECISION
		PVTessKernelRadius = stof( rootNode->first_node("PVTessKernelRadius")->value() );
#else
		PVTessKernelRadius = stod( rootNode->first_node("PVTessKernelRadius")->value() );
#endif

	if (0 != rootNode->first_node("PVTessCubeVoxelEdgeLen"))
#ifdef SINGLE_PRECISION
		PVTessCubeVoxelEdgeLen = stof( rootNode->first_node("PVTessCubeVoxelEdgeLen")->value() );
#else
		PVTessCubeVoxelEdgeLen = stod( rootNode->first_node("PVTessCubeVoxelEdgeLen")->value() );
#endif

	if (0 != rootNode->first_node("PVTessGuardZoneEdgeLen"))
		PVTessGuardZoneEdgeLen = static_cast<unsigned int>(stoi(rootNode->first_node("PVTessGuardZoneEdgeLen")->value() ));


//VISUALIZATION FILES WRITTEN OR NOT?
	if (0 != rootNode->first_node("VisIPGridWithTextureID"))
		VisIPGridWithTextureID = static_cast<unsigned int>(stoi(rootNode->first_node("VisIPGridWithTextureID")->value() ));
	if (0 != rootNode->first_node("VisIPGridWithGrainID"))
		VisIPGridWithGrainID = static_cast<unsigned int>(stoi(rootNode->first_node("VisIPGridWithGrainID")->value() ));
	if (0 != rootNode->first_node("VisIPGridWithPerImages"))
		VisIPGridWithPerImages = static_cast<unsigned int>(stoi(rootNode->first_node("VisIPGridWithPerImages")->value() ));
	if (0 != rootNode->first_node("VisHigherOrderNeighbor"))
		VisHigherOrderNeighbor = static_cast<unsigned int>(stoi(rootNode->first_node("VisHigherOrderNeighbor")->value() ));
	if (0 != rootNode->first_node("VisPeriodicReplicaIPs"))
		VisPeriodicReplicaIPs = static_cast<unsigned int>(stoi(rootNode->first_node("VisPeriodicReplicaIPs")->value() ));
	if (0 != rootNode->first_node("VisDBScanClusterIDs"))
		VisDBScanClusterIDs = static_cast<unsigned int>(stoi(rootNode->first_node("VisDBScanClusterIDs")->value() ));
	if (0 != rootNode->first_node("VisReconstructedGrain"))
		VisReconstructedGrain = static_cast<unsigned int>(stoi(rootNode->first_node("VisReconstructedGrain")->value() ));
	if (0 != rootNode->first_node("VisRVE27BVH"))
		VisRVE27BVH = static_cast<unsigned int>(stoi(rootNode->first_node("VisRVE27BVH")->value() ));
	if (0 != rootNode->first_node("VisGrainLocalVoxelGrid"))
		VisGrainLocalVoxelGrid = static_cast<unsigned int>(stoi(rootNode->first_node("VisGrainLocalVoxelGrid")->value() ));
	if (0 != rootNode->first_node("VisGrainFinalReconstruction"))
		VisGrainFinalReconstruction = static_cast<unsigned int>(stoi(rootNode->first_node("VisGrainFinalReconstruction")->value() ));


//##MK::DEBUGGING and DEVELOPMENT OF NEW FUNCTIONALITY
	if (0 != rootNode->first_node("VisGrainQuaternionClouds"))
		VisGrainQuaternionClouds = static_cast<unsigned int>(stoi(rootNode->first_node("VisGrainQuaternionClouds")->value() ));

	if (0 != rootNode->first_node("DebugDouble"))
		DebugDouble = stod( rootNode->first_node("DebugDouble")->value() );
	if (0 != rootNode->first_node("DebugUnsignedInt"))
		DebugUnsignedInt = stod( rootNode->first_node("DebugUnsignedInt")->value() );
}


bool Settings::checkUserInput()
{
	//MK::check user input for validity and good sense
	//default mode is E_DEBUG
	if ( Settings::IncrementOffset == 0 ) return false;

	if ( Settings::KernelRadius < EPSILON ) return false;
	if ( Settings::KernelRadius > 0.50 ) return false;
	if ( Settings::CriticalDisoriAngle <= DEGREE2RADIANT(0.1) || Settings::CriticalDisoriAngle >= DEGREE2RADIANT(20.0) ) return false;
	if ( Settings::CommDetectionDisoriAngleWght < EPSILON ) return false;
	//##MK::define reasonable values here
	if ( Settings::CommDetectionDistanceWght < EPSILON ) return false;
	//if ( Settings::GrainReconLocalDisoriAngle <= DEGREE2RADIANT(0.5) || Settings::GrainReconLocalDisoriAngle >= (20.0) ) return false;
	if ( Settings::GrainReconInspectionRadius < EPSILON ) return false;
	if ( Settings::GrainReconInspectionRadius >= Settings::KernelRadius ) return false;
	if ( Settings::GrainReconLouvainPrecision < EPSILON ) return false;
	if ( Settings::DBScanKernelRadius < EPSILON ) return false;
	if ( Settings::PVTessKernelRadius < EPSILON ) return false;
	if ( Settings::PVTessCubeVoxelEdgeLen < EPSILON ) return false;
	if ( Settings::PVTessCubeVoxelEdgeLen > Settings::PVTessKernelRadius) return false;
	if ( Settings::PVTessGuardZoneEdgeLen > 5) return false; //MK::implicitly takes care about potentially negative values as userinput because of bit wraparound
	//###MK::during candidate inspection distances will be inspected in order of their distance
 	//if ( Settings::PhysicalDomainSize <= 0.0 ) return false;

	//contains file type specification?
	string key (".");
	if (SpectralOutFilename.find(key) != string::npos) { //contains an ending specifier, check if that reads as 'spectralOut'
		if( SpectralOutFilename.compare(SpectralOutFilename.size()-11,11,"spectralOut") != 0 ) {//does not read as such
			return false;
		}
		string stripped = SpectralOutFilename.substr(0,SpectralOutFilename.size()-12);
		SpectralOutFilename = stripped;
	}
	//contains no file type specification leave it with this
	return true;
}

void Settings::displaySettings()
{
	//set precision of cout
#ifdef SINGLE_PRECISION
	cout << setprecision(9) << endl;
#else
	cout << setprecision(18) << endl;
#endif
	cout << Settings::SpectralOutFilename << endl;
	cout << "Analyzing increments in [" << Settings::IncrementFirst << ", ";
	cout << Settings::IncrementOffset << ", " << Settings::IncrementLast << "]" << endl;

	cout << "AnalysisID\t\t\t" << Settings::SimID << endl;

	cout << "KernelRadius\t\t\t" << Settings::KernelRadius << endl;
	cout << "CriticalDisoriAngle\t\t" << RADIANT2DEGREE(Settings::CriticalDisoriAngle) << endl;
	//cout << "GrainReconLocalDisoriAngle\t" << RADIANT2DEGREE(Settings::GrainReconLocalDisoriAngle) << endl;
	//cout << "PhysDomainSize\t\t\t" << Settings::PhysicalDomainSize << endl;

	cout << "The following tasks are planned:" << endl;
		cout << "\t-->addFlowcurve" << endl;
	if ( Settings::AnalyzeAddStrainTensors == true )
		cout << "\t-->addStrainTensors" << endl;
	if ( Settings::AnalyzeAddCauchy == true )
		cout << "\t-->addCauchy" << endl;
	if ( Settings::AnalyzeAddVonMises == true )
		cout << "\t-->addVonMises" << endl;
	if ( Settings::AnalyzeAddDisplacements == true)
		cout << "\t-->addDisplacements" << endl;
	if ( Settings::AnalyzeStateVarSpatialDistr == true ) {
		cout << "\t-->addAnalyzeStateVariableSpatialDistributions" << endl;
		if ( Settings::StateVarThrshldModel == E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD )
			cout << "\t-->-->Closest neighboring unique ip with disori angle larger threshold" << endl;
		if ( Settings::StateVarThrshldModel == E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION )
			cout << "\t-->-->Hierarchical clustering-based grain detection with signed distance function" << endl;
		if ( Settings::ThrshldAgainstDefConfig == true)
			cout << "\t-->-->Deformed configuration of integration point grid" << endl;
		else
			cout << "\t-->-->Initial configuration of integration point grid" << endl;
	}
	if ( Settings::AnalyzeReconstructGrains == true ) {
		cout << "\t-->Grain reconstruction" << endl;
		cout << "\t-->-->Modularity-based community detection" << endl;
		cout << "\t-->-->DisoriAngleWeightingFactor " << Settings::CommDetectionDisoriAngleWght << endl;
		//cout << "\t-->-->DistanceWeightingFactor " << Settings::CommDetectionDistanceWght << endl; //##MK::currently disabled

		//cout << "\t-->-->LocalDisoriAngle " << RADIANT2DEGREE(Settings::GrainReconLocalDisoriAngle) << endl;
		cout << "\t-->-->InspectionRadius " << Settings::GrainReconInspectionRadius << endl;
		cout << "\t-->-->LouvainPrecision " << Settings::GrainReconLouvainPrecision << endl;
		cout << "\t-->-->DBScanPerImageRadius " << Settings::DBScanKernelRadius << endl;

		cout << "\t-->Voxelization prior to signed-distance function computation" << endl;
		cout << "\t-->-->PVCubeVoxelEdgeLenght " << Settings::PVTessCubeVoxelEdgeLen << endl;
		cout << "\t-->-->PVTessKernelRadius " << Settings::PVTessKernelRadius << endl;
		cout << "\t-->-->PVTessGuardZone " << Settings::PVTessGuardZoneEdgeLen << " voxel" << endl;
		//##MK::add further parameter
	}

}


unsigned int RangeLimits::DummyID = std::numeric_limits<unsigned int>::max();
unsigned int RangeLimits::UIntMax = RangeLimits::DummyID - 1;
unsigned int RangeLimits::NXYZMax = 1024*1024*1024;
//hardly any DAMASK simulation in the near future is expected at all to become as large...


void RangeLimits::displaySettings()
{
#ifdef SINGLE_PRECISION
	cout << setprecision(9) << endl;
#else
	cout << setprecision(18) << endl;
#endif

	cout << "DummyID\t\t\t\t" << RangeLimits::DummyID << endl;
	cout << "UIntMax\t\t\t\t" << RangeLimits::UIntMax << endl;
	cout << "NXYZMax\t\t\t\t" << RangeLimits::NXYZMax << endl;
}

