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

#ifndef __PDT_TENSORMATH_H__
#define __PDT_TENSORMATH_H__

#include "PDT_DatatypesMPI.h"
//#include "PDT_IntelMKL.h"

bv3x3 leftmult( t3x3 const & defgrad, bv3x3 const & bvecs );

t3x3 transpose(t3x3 const &in);

t3x3 diag(t3x1 const &in);

real_m33 det(t3x3 const &in);

real_m33 dot(t3x1 const &in1, t3x1 const &in2);

void killnoise(t3x3 &m, real_m33 thrshld);

real_m33 trace(t3x3 const &in);

void eye(t3x3 &out);

t3x3 eye(void);

void zerot3x1(t3x1 &out);

t3x1 zerot3x1(void);

t3x3 failt3x3(void);

void zerot3x3(t3x3 &out);

t3x3 zerot3x3(void);

t3x3 add(t3x3 const &in1, t3x3 const &in2);

t3x3 substrct(t3x3 const &in1, t3x3 const &in2);

//void log( t3x3 &inout);

t3x1 log(t3x1 const &in);

real_m33 sum(t3x3 const &in);

t3x3 mult(real_m33 sc, t3x3 const &in);

t3x3 dyadic(t3x3 const &in1, t3x3 const &in2);

t3x3 dot(t3x3 const &in1, t3x3 const &in2);

#endif
