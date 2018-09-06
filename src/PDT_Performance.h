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

#ifndef __PDT_PERFORMANCE_H__
#define __PDT_PERFORMANCE_H__

//redefine allocator classes, to allow single-source switch between default STL and boostSIMD allocator
//##MK::this is an intermediate solution and very likely requires improvement as an STL default allocated
//std::vector<T, std::allocator<T>> invokes new, which is agnostic of strict enough alignment for SIMD
//infact new allocates with alignment(sizeof(T)) while alignment to cache boundaries
//i.e. alignment(32) if not (64) may be necessary


/*typedef std::allocator<unsigned int> ui_allo;
typedef std::allocator<int> i_allo;
typedef std::allocator<real> r_allo;
typedef std::allocator<real_ori> r_ori_allo;
typedef std::allocator<real_rho> r_rho_allo;
typedef std::allocator<real_xyz> r_xyz_allo;
typedef std::allocator<real_m33> r_m33_allo;*/
//##MK::replace with boostSIMD allocator versions



#include "PDT_Settings.h"

#define INCREMENTAL_KERNEL_RADIUS_STEPPING		2

class Performance {
public:

	static size_t MPIIOStripeSizeMin;		//most efficient portion to read on RAID system by single process at a time in bytes
	static size_t MPIIOStripeSizeMax;
	static size_t MPIIOStripeSize;

	static size_t CachePageSize;			//target system page size, typically 4096B
	static size_t CacheIntegerWidth;		//how many int32 or uint32 values fitting in page
	static size_t CacheSingleWidth;			//how many single precision floating point
	static size_t CacheDoubleWidth;			//how many double precision

	static size_t SIMDRegisterSize;			//how many bytes in SIMD register, CPU architecture dependent 2017-th Xeons typically 256Bit i.e 256Bit/32Bit = 8, 512Bit/64Bit=8
	static size_t SIMDIntegerWidth;
	static size_t SIMDSingleWidth;			//how many single precision floating point in SIMD register
	static size_t SIMDDoubleWidth;			//how many double precision ...

	static size_t ThreadBufferAllocSize;

	static void displaySettings();
};

#endif
