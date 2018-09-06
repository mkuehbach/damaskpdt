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

#ifndef __PDT_CRYSTALLOGRAPHY_H__
#define __PDT_CRYSTALLOGRAPHY_H__

//#include "PDT_IntelMKL.h"
#include "PDT_Datatypes.h"

struct slipsysdata_fcc
{
	real_rho b2;	//fcc primary slip system related data utilizing the Schmid-Boas notation
	real_rho b4;
	real_rho b5;
	real_rho c1;
	real_rho c3;
	real_rho c5;
	real_rho a2;
	real_rho a3;
	real_rho a6;
	real_rho d1;
	real_rho d4;
	real_rho d6;
    slipsysdata_fcc() :
    	b2(numeric_limits<real_rho>::max()),
		b4(numeric_limits<real_rho>::max()),
		b5(numeric_limits<real_rho>::max()),
		c1(numeric_limits<real_rho>::max()),
		c3(numeric_limits<real_rho>::max()),
		c5(numeric_limits<real_rho>::max()),
		a2(numeric_limits<real_rho>::max()),
		a3(numeric_limits<real_rho>::max()),
		a6(numeric_limits<real_rho>::max()),
		d1(numeric_limits<real_rho>::max()),
		d4(numeric_limits<real_rho>::max()),
		d6(numeric_limits<real_rho>::max()) {}
    slipsysdata_fcc(
    		const real_rho _b2,
			const real_rho _b4,
			const real_rho _b5,
			const real_rho _c1,
			const real_rho _c3,
			const real_rho _c5,
			const real_rho _a2,
			const real_rho _a3,
			const real_rho _a6,
			const real_rho _d1,
			const real_rho _d4,
			const real_rho _d6 ) :
				b2(_b2), b4(_b4), b5(_b5),
				c1(_c1), c3(_c3), c5(_c5),
				a2(_a2), a3(_a3), a6(_a6),
				d1(_d1), d4(_d4), d6(_d6) {}
	~slipsysdata_fcc() {}
};

std::ostream& operator << (std::ostream& in, slipsysdata_fcc const & val);

#endif
