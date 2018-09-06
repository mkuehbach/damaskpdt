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

#ifndef __PDT_DTPMPI_H__
#define __PDT_DTPMPI_H__

//#include "PDT_Datatypes.h"
#include "PDT_Crystallography.h"
//container datastructures for MPI intra-process communication to reduce amount of MPI messages and file interaction


struct MPI_Tensor3x3_Double
{
	double a11;
	double a12;
	double a13;
	double a21;
	double a22;
	double a23;
	double a31;
	double a32;
	double a33;
	MPI_Tensor3x3_Double() :
		a11(numeric_limits<double>::max()), a12(numeric_limits<double>::max()), a13(numeric_limits<double>::max()),
		a21(numeric_limits<double>::max()), a22(numeric_limits<double>::max()), a23(numeric_limits<double>::max()),
		a31(numeric_limits<double>::max()), a32(numeric_limits<double>::max()), a33(numeric_limits<double>::max()) {}
	MPI_Tensor3x3_Double(t3x3 const & in) :
		a11(static_cast<double>(in.a11)), a12(static_cast<double>(in.a12)), a13(static_cast<double>(in.a13)),
		a21(static_cast<double>(in.a21)), a22(static_cast<double>(in.a22)), a23(static_cast<double>(in.a23)),
		a31(static_cast<double>(in.a31)), a32(static_cast<double>(in.a32)), a33(static_cast<double>(in.a33)) {}
	~MPI_Tensor3x3_Double() {}
};


//typedef struct
struct MPI_VonMises_Double
{
	double eps;
	double sigma;
	MPI_VonMises_Double() :
		eps(numeric_limits<double>::max()), sigma(numeric_limits<double>::max()) {}
	MPI_VonMises_Double( vMises const & in) :
		eps(static_cast<double>(in.vMisesEquivStrain)),
		sigma(static_cast<double>(in.vMisesEquivStress)) {}
	~MPI_VonMises_Double() {}
};


struct MPI_DisloSpatDistr_Double
{
	double d;

	double b2;	//fcc primary slip system related data utilizing the Schmid-Boas notation
	double b4;
	double b5;
	double c1;
	double c3;
	double c5;
	double a2;
	double a3;
	double a6;
	double d1;
	double d4;
	double d6;
    MPI_DisloSpatDistr_Double() :
    	d(numeric_limits<double>::max()),
    	b2(static_cast<double>(0.0)),
		b4(static_cast<double>(0.0)),
		b5(static_cast<double>(0.0)),
		c1(static_cast<double>(0.0)),
		c3(static_cast<double>(0.0)),
		c5(static_cast<double>(0.0)),
		a2(static_cast<double>(0.0)),
		a3(static_cast<double>(0.0)),
		a6(static_cast<double>(0.0)),
		d1(static_cast<double>(0.0)),
		d4(static_cast<double>(0.0)),
		d6(static_cast<double>(0.0)) {}
    MPI_DisloSpatDistr_Double( const real_xyz _d, const slipsysdata_fcc _sd ) :
    	d(static_cast<double>(_d)),
		b2(static_cast<double>(_sd.b2)), b4(static_cast<double>(_sd.b4)), b5(static_cast<double>(_sd.b5)),
		c1(static_cast<double>(_sd.c1)), c3(static_cast<double>(_sd.c3)), c5(static_cast<double>(_sd.c5)),
		a2(static_cast<double>(_sd.a2)), a3(static_cast<double>(_sd.a3)), a6(static_cast<double>(_sd.a6)),
		d1(static_cast<double>(_sd.d1)), d4(static_cast<double>(_sd.d4)), d6(static_cast<double>(_sd.d6)) {}
    MPI_DisloSpatDistr_Double( const real_xyz _d, vector<real_rho> const & in, size_t rs ) :
    	d(static_cast<double>(_d)),
		b2(static_cast<double>(in[rs+0])), b4(static_cast<double>(in[rs+1])), b5(static_cast<double>(in[rs+2])),
		c1(static_cast<double>(in[rs+3])), c3(static_cast<double>(in[rs+4])), c5(static_cast<double>(in[rs+5])),
		a2(static_cast<double>(in[rs+6])), a3(static_cast<double>(in[rs+7])), a6(static_cast<double>(in[rs+8])),
		d1(static_cast<double>(in[rs+9])), d4(static_cast<double>(in[rs+10])), d6(static_cast<double>(in[rs+11])) {}
};



#endif
