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

#ifndef __PDT_VTKIO_H__
#define __PDT_VTKIO_H__


#include "PDT_OriMath.h"

bool vtk_grfinalrecon( vxlgrd const & ggrd, vxlgrd const & lgrd,
		vector<unsigned int> const & gid, vector<unsigned int> const & uipid,
			vector<real_sdf> const & sdfbuf1, vector<real_sdf> const & sdfbuf2,
				const unsigned int id, const unsigned int increment  );

bool vtk_vxlgrdm2( vxlgrd const & ggrd, vxlgrd const & lgrd,
		vector<unsigned int> const & m1, vector<real_sdf> const & m2,
			const unsigned int id, const unsigned int increment,
				const string what, const string descr, const string m1what, const string m2what );

bool vtk_p3d(	vector<p3d> const & pp3, const unsigned int id,
					const unsigned int increment, const string what );

bool vtk_p3dm1(	vector<p3dm1> const & pp3, const unsigned int id,
					const unsigned int increment, const string what );

bool vtk_bvh_p3dm1( vector<vector<p3dm2>*> const & pp3, const unsigned int increment,
		const string what, const string dsrc );


bool vtk_p3dm1( vector<p3d> const & pp3, vector<unsigned int> const & m1,
		const unsigned int id, const unsigned int increment, const string what );


bool vtk_p3dm3( 	vector<vector<p3d>*> const & xyz,
						vector<vector<unsigned int>*> const & m,
								vector<vector<unsigned char>*> const & img,
									const unsigned int increment);


bool vtk_nbhd3d( vector<p3d> const & p,
					vector<unsigned int> const & pid,
						vector<real_xyz> const & dist,
							vector<real_ori> const & theta,
								vector<unsigned int> const & textureid,
									const unsigned int increment);

bool vtk_p3dm3( vector<p3d> const & p,
					vector<unsigned int> const & textureid,
						const unsigned int increment);

bool vtk_gp3d( vector<p3d> const & p,
					vector<unsigned int> const & grainid,
						const unsigned int increment);

#endif
