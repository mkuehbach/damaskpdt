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


#ifndef __PDT_XDMF_H__
#define __PDT_XDMF_H__

#include "PDT_HDF5.h"

#define XDMF_HEADER_LINE1				"<?xml version=\"1.0\" ?>"
#define XDMF_HEADER_LINE2				"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"
#define XDMF_HEADER_LINE3				"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">"

#define WRAPPED_XDMF_SUCCESS				+1
#define WRAPPED_XDMF_IOFAILED				-1

/*
class xdmfHdl
{
//coordinating instance handling all (sequential) writing of XDMF metafiles detailing HDF5 additional metadata
//for visualization for Paraview or VisIt
public:
	xdmfHdl();
	~xdmfHdl();

	//file generation and closing
	int create_volrecon_file( const string xmlfn, const size_t nions, const string h5ref );
	int create_crystalloxyz_file( const string xmlfn, const size_t npoints, const string h5ref );
	int create_debugsynthesis_file( const string xmlfn, const size_t nions, const string h5ref );
	//int create_iondistance_file( const string xmlfn, const size_t nions, const string h5ref );
	int create_tipsurface_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const string h5ref );
	int create_voronoicell_vis_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref );
	int create_voronoicell_vol_file( const string xmlfn, const size_t ncells, const string h5ref );

	int create_voronoicell_debug_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref );


	int create_tapsim_recon_vs_synth_file( const string xmlfn, const size_t nions, const string h5ref );

private:
	ofstream xdmfout;
};

*/

#endif
