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

#ifndef __PDT_SETTINGS_H__
#define __PDT_SETTINGS_H__


#include "PDT_Profiler.h"

//add thirdparty XML reader header library functionality by M. Kalicinski
#include "thirdparty/RapidXML/rapidxml.hpp"


enum E_ANALYSIS_MODE {
	E_DEBUG		//developer, comma to separate elements except the last E1, E2, E3
};

enum E_THRESHOLDING_MODE {
	E_CLOSESTNBORIP_DISORIANGLE_THRESHOLD,
	E_HIERARCHICALCOMMUNITY_SIGNEDDISTANCEFUNCTION //E_SLIPSYSTEM_COMPATIBILITY_BASED
};

enum E_GRAINRECON_MODE {
	E_HIERARCHICAL_CLUSTERING,
	E_INITIAL_TEXTUREID
};

class Settings {
public:

	static E_ANALYSIS_MODE AnalysisMode;					//the work to do
	//static E_THRESHOLDING_MODE SpatialDistrThrshldModel;	//which is the decisive quantity against which one thresholds distance when quantifying the accumulation of state variable values
	static string SpectralOutFilename;
	static unsigned int IncrementFirst;						//which strain increments to analyze
	static unsigned int IncrementOffset;
	static unsigned int IncrementLast;

	//##MK::several analysis modes internally require calls to the IntelMKL library
	static bool AnalyzeFlowCurve;
	static bool AnalyzeAddStrainTensors;
	static bool AnalyzeAddCauchy;
	static bool AnalyzeAddVonMises;						//von Mises true strain and stresses
	static bool AnalyzeAddDisplacements;				//project from initial coordinates to deformed configuration
	static bool AnalyzeStateVarSpatialDistr;			//characterize spatial distribution of state variable close to interfaces
	static bool AnalyzeReconstructGrains;

	static unsigned int GrainReconModel;
	static unsigned int StateVarThrshldModel;
	static bool ThrshldAgainstDefConfig;
	static unsigned long SimID;							//consistent identifier of the data analysis run
	static real_xyz KernelRadius;
	static real_ori CriticalDisoriAngle;
	static real_ori CommDetectionDisoriAngleWght;
	static real_xyz CommDetectionDistanceWght;
	static real_ori GrainReconLocalDisoriAngle;
	static real_xyz GrainReconInspectionRadius;
	static long double GrainReconLouvainPrecision;

	static real_xyz DBScanKernelRadius;
	static real_xyz PVTessKernelRadius;
	static real_xyz PVTessCubeVoxelEdgeLen;
	static unsigned int PVTessGuardZoneEdgeLen;

//DEBUG and VISUALIZATION OPTIONS
	static unsigned int VisIPGridWithTextureID;
	static unsigned int VisIPGridWithGrainID;
	static unsigned int VisIPGridWithPerImages;
	static unsigned int VisHigherOrderNeighbor;
	static unsigned int VisPeriodicReplicaIPs;
	static unsigned int VisDBScanClusterIDs;
	static unsigned int VisReconstructedGrain;
	static unsigned int VisRVE27BVH;
	static unsigned int VisGrainLocalVoxelGrid;
	static unsigned int VisGrainFinalReconstruction;

	static unsigned int VisGrainQuaternionClouds;
	static double DebugDouble;
	static unsigned int DebugUnsignedInt;

	static void readUserInput(string filename = "");
	static bool checkUserInput();
	static void displaySettings();
};


class RangeLimits {
public:
	static unsigned int DummyID;
	static unsigned int UIntMax;
	static unsigned int NXYZMax;

	static void displaySettings();
};


#endif
