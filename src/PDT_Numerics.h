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


#ifndef __PDT_NUMERICS_H__
#define __PDT_NUMERICS_H__

//single or double precision for elementary floating point values
//#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
	typedef float real;
	#define EPSILON (1.0e-6)
#else
	#define DOUBLE_PRECISION
	typedef double real;
	#define EPSILON (1.0e-12)
#endif

//user-derived adjustable precision for different core functionalities
typedef real voro_real;			//tessellation library interface
typedef size_t voro_uint;
typedef int voro_int;			//MK::int is demanded by Voro++ library!
typedef real real_ori;			//orientation math
typedef real real_rho;			//microstructure state variable values related to dislocations
typedef real real_xyz;			//positional data
typedef real real_m33;			//tensor data, deformation gradients, etc...

typedef real real_sdf;			//accuracy when computing signed-distance function values



//type range
#define UCHARMX							numeric_limits<unsigned char>::max()
#define UCHARMI							numeric_limits<unsigned char>::lowest()
#define INT16MX							numeric_limits<short>::max()
#define INT16MI							numeric_limits<short>::min()
#define INT32MX							numeric_limits<int>::max()
#define INT32MI							numeric_limits<int>::lowest()
#define UINT16MX						numeric_limits<unsigned short>::max()
#define UINT16MI						numeric_limits<unsigned short>::lowest()
#define UINT32MX						numeric_limits<unsigned int>::max()
#define UINT32MI						numeric_limits<unsigned int>::lowest()
#define UINT64MX						numeric_limits<size_t>::max()
#define UINT64MI						numeric_limits<size_t>::lowest()
#define UINT64XX						-1
#define F64MX							numeric_limits<double>::max()
#define F64MI							numeric_limits<double>::lowest()

//type range user types
#define VOROUINTMX						numeric_limits<voro_uint>::max()
#define VOROUINTMI						numeric_limits<voro_uint>::lowest()
#define VOROUINTXX						-1
#define VOROINTMX						numeric_limits<voro_int>::max()
#define VOROINTMI						numeric_limits<voro_int>::lowest()
#define VOROFMX							numeric_limits<voro_real>::max()
#define VOROFMI							numeric_limits<voro_real>::lowest()




//signum function
template <typename T> int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}


//SDF initialization
#define	FASTSWEEPING_INITIAL_VALUE		(-1.0)

//tessellation specific
#define VORO_COMPOSER_VERTEX_MAPPING_HASHVAL	2097152		//2^21 if spreading the bit range of size_t for x,y, z
#define VOROZERO 				static_cast<voro_real>(0.0)


//random
#define MT19937SEED						(-1)
#define MT19937WARMUP					(700000)



//flags
#define YES					(0xFF)
#define NO					(0x00)


#define LOUVAIN_DEFAULT_PRECISION 			(0.000001L)

#endif
