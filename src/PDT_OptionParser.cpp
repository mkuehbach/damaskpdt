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

void option::add_key(const string _k, const string _expl)
{
	keys.push_back(_k);
	explanation.push_back(_expl);
}


bool option::exists_key(const string _k)
{
	//check all options if they exist
	for (size_t i = 0; i < this->keys.size(); ++i) {
		if ( _k.compare( keys.at(i)) != 0 ) { //option does not exist
			continue;
		}
		else {
			return true;
		}
	}
	//all options checked but not returned, so key not existent
	return false;
}


string option::get_expl(const string _k)
{
	//check all options if they exist
	for (size_t i = 0; i < this->keys.size(); ++i) {
		if ( _k.compare( keys.at(i)) != 0 ) { //option does not exist
			continue;
		}
		else {
			return explanation.at(i);
		}
	}
	//all options checked but not returned, so key not existent
	string serror = _k + "---> such a key does not exist!";
	return serror;
}



optparser::optparser()
{
}


optparser::~optparser()
{
}


void optparser::define_options()
{
	string k = "--h";
	string e = "damaskpdt v" + to_string(VERSION_MAJOR) + "." + to_string(VERSION_MINOR) + "." + to_string(VERSION_REVISION);
	opts.push_back( option() );
	opts.back().add_key(k,e);
}
