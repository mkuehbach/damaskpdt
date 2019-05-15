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



#include "PDT_XDMF.h"

/*

xdmfHdl::xdmfHdl()
{
}


xdmfHdl::~xdmfHdl()
{
}


int xdmfHdl::create_volrecon_file( const string xmlfn, const size_t nions, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << nions << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << 3*nions << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Iontype\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 1\" DataType=\"UInt\" Precision=\"1\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLRECON_IONTYPE_IDS << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		if ( Settings::IOIonTipSurfDists == true ) {
			xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"SqrdDist\">" << "\n";
			xdmfout << "	    <DataItem Dimensions=\"" << nions << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
			xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_SURFDISTSQR << "\n";
			xdmfout << "        </DataItem>" << "\n";
			xdmfout << "      </Attribute>" << "\n";
		}
 		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_crystalloxyz_file( const string xmlfn, const size_t npoints, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"crystxyz\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << npoints << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << 3*npoints << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_CRYSTALLO_MATPOINT_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << npoints << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_CRYSTALLO_MATPOINT_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"PointID\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << npoints << " 1\" DataType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_CRYSTALLO_MATPOINT_IDS << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
 		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";


		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_debugsynthesis_file( const string xmlfn, const size_t nions, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"debugsynth\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << nions << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << 3*nions << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_DEBUG_APTOIM_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_DEBUG_APTOIM_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"GrainID\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 1\" DataType=\"UInt\" Precision=\"2\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_DEBUG_APTOIM_GID << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
 		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_tipsurface_file( const string xmlfn, const size_t topo_nelements,
		const size_t topo_dims, const size_t geom_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"alphashape\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_SURFRECON_ASHAPE_HULL_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" <<  PARAPROBE_SURFRECON_ASHAPE_HULL_GEOM << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_voronoicell_vis_file( const string xmlfn, const size_t topo_nelements,
		const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoi\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType\"UInt\" Precision=\"8\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_TOPOLOGY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_GEOMETRY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}



int xdmfHdl::create_voronoicell_vol_file( const string xmlfn, const size_t ncells, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoicells\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Polyvertex\" NodesPerElement=\"" << ncells << "\">" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << " 3\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_CELLPOS << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute Name=\"Volume\" AttributeType=\"Scalar\" Center=\"Node\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_VOL << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "      <Attribute Name=\"ThreadID\" AttributeType=\"Scalar\" Center=\"Node\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << "\" NumberType=\"UChar\" Precision=\"1\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_THREADID << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_voronoicell_debug_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoi\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType=\"UInt\" Precision=\"8\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_TOPOLOGY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_GEOMETRY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
#ifndef VALIDZONE_IONS_ONLY
		xdmfout << "      <Attribute Name=\"ThreadID\" AttributeType=\"Scalar\" Center=\"Cell\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << attr_dims << " 1\" DataType=\"Int\" Precision=\"2\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_THREADIDATTR << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
#endif
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}



int xdmfHdl::create_tapsim_recon_vs_synth_file( const string xmlfn, const size_t nions, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"tapsim recon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << nions << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << 3*nions << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << TAPSIM_ION_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << TAPSIM_RECON_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Distance\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << TAPSIM_ONE_TO_ONE_DIST << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Vector\" Center=\"Node\" Name=\"DiffVec\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << TAPSIM_ONE_TO_ONE_DXYZ << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";


		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}

*/
