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


#include "PDT_Performance.h"

size_t Performance::MPIIOStripeSizeMin = 1*1024*1024; //B
size_t Performance::MPIIOStripeSizeMax = 10*1024*1024; //B
size_t Performance::MPIIOStripeSize = Performance::MPIIOStripeSizeMax;

size_t Performance::CachePageSize = 4096; //B
size_t Performance::CacheIntegerWidth = Performance::CachePageSize / sizeof(int);
size_t Performance::CacheSingleWidth = Performance::CachePageSize / sizeof(float);
size_t Performance::CacheDoubleWidth = Performance::CachePageSize / sizeof(double);

size_t Performance::SIMDRegisterSize = 256 / 8; //bits/(bits/byte)
size_t Performance::SIMDIntegerWidth = Performance::SIMDRegisterSize / sizeof(int);
size_t Performance::SIMDSingleWidth = Performance::SIMDRegisterSize / sizeof(float);
size_t Performance::SIMDDoubleWidth = Performance::SIMDRegisterSize / sizeof(double);

size_t Performance::ThreadBufferAllocSize = 1024*1024; //sizeof MB equivalent

//##MK::add function to change globally at runtime for all processes as in the Settings class


void Performance::displaySettings()
{
#ifdef SINGLE_PRECISION
	cout << setprecision(9) << endl;
#else
	cout << setprecision(18) << endl;
#endif

	cout << "MPIIOStripeSize\t\t\t" << Performance::MPIIOStripeSize << " Bytes" << endl;
	cout << "CacheIntegerWidth\t\t" << Performance::CacheIntegerWidth << " elements" << endl;
	cout << "CacheSingleWidth\t\t" << Performance::CacheSingleWidth << " elements" << endl;
	cout << "CacheDoubleWidth\t\t" << Performance::CacheDoubleWidth << " elements" << endl;
	cout << "SIMDIntegerWidth\t\t" << Performance::SIMDIntegerWidth << " elements" << endl;
	cout << "SIMDSingleWidth\t\t\t" << Performance::SIMDSingleWidth << " elements" << endl;
	cout << "SIMDDoubleWidth\t\t\t" << Performance::SIMDDoubleWidth << " elements" << endl;
	cout << "SizeTMaxValue\t\t\t" << std::numeric_limits<size_t>::max() << " elements" << endl;
	cout << endl;
}
