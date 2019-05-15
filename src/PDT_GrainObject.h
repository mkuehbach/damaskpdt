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

#ifndef __PDT_GRAINOBJECT_H__
#define __PDT_GRAINOBJECT_H__


//#include "PDT_Crystallography.h"
#include "PDT_IntelMKL.h"

struct lvwtedge
{
	double wt;				//weight of edge ij
	unsigned int src;		//ip i
	unsigned int dest;		//ip j
	float disori;			//disorientation angle
	float distance;			//distance float

	lvwtedge() : wt(0.0), src(numeric_limits<unsigned int>::max()), dest(numeric_limits<unsigned int>::max()),
			disori(0.f), distance(0.f) {}
	lvwtedge( const double _wt, const unsigned int _src, const unsigned int _dst,
			const float _dis, const float _dist ) :
		wt(_wt), src(_src), dest(_dst), disori(_dis), distance(_dist) {}
	~lvwtedge(){}
};


struct grain
{
	squat ori;
	unsigned int gid;
	unsigned int np;
    grain() :
    	ori(squat()), gid(numeric_limits<unsigned int>::max()), np(0) {}
    grain( const squat _ori, const unsigned int _gid) :
    	ori(_ori), gid(_gid), np(0) {}
    grain( const squat _ori, const unsigned int _gid, const unsigned int _np) :
    	ori(_ori), gid(_gid), np(_np) {}
    ~grain() {}
};

std::ostream& operator << (std::ostream& in, grain const & val);



#endif
