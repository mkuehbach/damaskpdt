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


#include "PDT_VTKIO.h"

bool vtk_grfinalrecon( vxlgrd const & ggrd, vxlgrd const & lgrd,
		vector<unsigned int> const & gid, vector<unsigned int> const & uipid,
			vector<real_sdf> const & sdfbuf1, vector<real_sdf> const & sdfbuf2,
				const unsigned int id, const unsigned int increment  )
{
	//writes center positions of cubic voxel specified by lgrd with piecewise mark data m1 and m2 per voxel
	if (gid.size() != lgrd.nxyz || uipid.size() != lgrd.nxyz || sdfbuf1.size() != lgrd.nxyz ||  sdfbuf2.size() != lgrd.nxyz ) {
		cerr << "Inconsistent length of mark data per voxel lgrd.nxyz " << lgrd.nxyz << " sizes() gid/uipid/sdf1/sdf2 ";
		cerr << gid.size() << ";" << uipid.size() << ";" << sdfbuf1.size() << ";" << sdfbuf2.size() << "\n";
		return false;
	}

	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + ".FinalReconGrainID" + to_string(id) + ".vtk";
	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {

		size_t nvertices = lgrd.nxyz;
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT Final grain reconstruction voxelated " << id << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";

		real_xyz h = lgrd.dcell;
		real_xyz hhalf = static_cast<real_xyz>(0.5) * h;
		vxl org = lgrd.origin_discr_lnk;
		p3d globalorigin = ggrd.origin_cntum;

		for(size_t z = 0; z < lgrd.nz; ++z) {
			//go back to origin in global coordinate system bottom face walk half cell length add integer full increments of cell length
			real_xyz cz = globalorigin.z + hhalf + (h * static_cast<real_xyz>(org.z + z));
			for(size_t y = 0; y < lgrd.ny; ++y) {
				real_xyz cy = globalorigin.y + hhalf + (h * static_cast<real_xyz>(org.y + y));
				for(size_t x = 0; x < lgrd.nx; ++x) {
					real_xyz cx = globalorigin.x + hhalf + (h * static_cast<real_xyz>(org.x + x));
					vtk << cx << " " << cy << " " << cz << "\n"; //MK::centers of gravity of voxels not location of axis perpendicular voxel walls!
				}
			}
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 4\n";
		vtk << "GrainID 1 " << nvertices << " double\n"; //MK::this is a hack, unsigned int is unknown type...
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << gid[i] << "\n";
		}
		vtk << "\n";

		vtk << "UniqueIPID 1 " << nvertices << " double\n"; //MK::this is a hack, unsigned int is unknown type...
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << uipid[i] << "\n";
		}
		vtk << "\n";

#ifdef SINGLE_PRECISION
		vtk << "SDFinit 1 " << nvertices << " float\n";
#else
		vtk << "SDFinit 1 " << nvertices << " double\n";
#endif
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << sdfbuf1[i] << "\n";
		}
		vtk << "\n";

#ifdef SINGLE_PRECISION
		vtk << "SDFfsm 1 " << nvertices << " float\n";
#else
		vtk << "SDFfsm 1 " << nvertices << " double\n";
#endif
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << sdfbuf2[i] << "\n";
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write voxel-related data of grain " << id << " to VTK file" << "\n";
		return false;
	}
}


bool vtk_vxlgrdm2( vxlgrd const & ggrd, vxlgrd const & lgrd,
		vector<unsigned int> const & m1, vector<real_sdf> const & m2,
			const unsigned int id, const unsigned int increment,
				const string what, const string descr, const string m1what, const string m2what )
{
	//writes center positions of cubic voxel specified by lgrd with piecewise mark data m1 and m2 per voxel
	if (m1.size() != lgrd.nxyz || m2.size() != lgrd.nxyz) {
		cerr << "Missing mark data per voxel lgrd.nxyz " << lgrd.nxyz << " m1.size() " << m1.size() << " m2.size() " << m2.size() << "\n";
		return false;
	}

	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + "." + what + to_string(id) + ".vtk";
	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {

		size_t nvertices = lgrd.nxyz;
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT " << descr << " " << id << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";

		real_xyz h = lgrd.dcell;
		real_xyz hhalf = static_cast<real_xyz>(0.5) * h;
		vxl org = lgrd.origin_discr_lnk;
		p3d globalorigin = ggrd.origin_cntum;

		for(size_t z = 0; z < lgrd.nz; ++z) {
			real_xyz cz = globalorigin.z + hhalf + (h * static_cast<real_xyz>(org.z + z));
			for(size_t y = 0; y < lgrd.ny; ++y) {
				real_xyz cy = globalorigin.y + hhalf + (h * static_cast<real_xyz>(org.y + y));
				for(size_t x = 0; x < lgrd.nx; ++x) {
					real_xyz cx = globalorigin.x + hhalf + (h * static_cast<real_xyz>(org.x + x));
					vtk << cx << " " << cy << " " << cz << "\n"; //MK::centers of gravity of voxels not location of axis perpendicular voxel walls!
				}
			}
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 2\n";
		vtk << m1what << " 1 " << nvertices << " float\n"; //MK::this is a hack, unsigned int is unknown type...
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << m1.at(i) << "\n";
		}
		vtk << "\n";

#ifdef SINGLE_PRECISION
		vtk << m2what << " 1 " << nvertices << " float\n";
#else
		vtk << m2what << " 1 " << nvertices << " double\n";
#endif
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << m2.at(i) << "\n";
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write voxel with content to VTK file" << "\n";
		return false;
	}
}


bool vtk_p3d(	vector<p3d> const & pp3, const unsigned int id,
					const unsigned int increment, const string what )
{
	//writes 3d positions to VTK file
	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + "." + what + ".vtk";

	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {
		size_t nvertices = pp3.size();
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT RVE27 Periodic and Unique IPs for grain " << id << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << pp3[i].x << " " << pp3[i].y << " " << pp3[i].z << "\n";
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "GrainID 1 " << nvertices << " int\n"; //double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << id << "\n";
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write integration point positions to VTK file" << "\n";
		return false;
	}
}


bool vtk_p3dm1(	vector<p3dm1> const & pp3, const unsigned int id,
					const unsigned int increment, const string what )
{
	//writes 3d positions to VTK file
	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + "." + what + ".vtk";

	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {
		size_t nvertices = pp3.size();
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT RVE27 Periodic and Unique IPs for grain " << id << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << pp3[i].x << " " << pp3[i].y << " " << pp3[i].z << "\n";
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 2\n";
		vtk << "UniqueIPID 1 " << nvertices << " int\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << pp3[i].m1 << "\n";
		}
		vtk << "\n";

		vtk << "GrainID 1 " << nvertices << " int\n"; //double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << id << "\n";
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write integration point positions to VTK file" << "\n";
		return false;
	}
}


bool vtk_p3dm1( vector<p3d> const & pp3, vector<unsigned int> const & m1,
		const unsigned int id, const unsigned int increment, const string what )
{
	//writes 3d positions and marks to VTK file
	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + "." + what + ".vtk";

	if ( pp3.size() != m1.size() ) {
		cerr << "Input datasets dissimilar in length" << "\n";
		return false;
	}

	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {
		size_t nvertices = pp3.size();
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT RVE27 DBScan Replica Clustering for grain " << id << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << pp3[i].x << " " << pp3[i].y << " " << pp3[i].z << "\n";
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "DBClusterID 1 " << nvertices << " int\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << m1[i] << "\n";
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write integration point positions and marks to VTK file" << "\n";
		return false;
	}
}


bool vtk_bvh_p3dm1( vector<vector<p3dm2>*> const & pp3, const unsigned int increment, const string what, const string dscr)
{
	//writes points and bins of a bounded volume hierarchy to VTK file
	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(increment) + "." + what + ".vtk";

	size_t nvertices = 0;
	for(size_t i = 0; i < pp3.size(); ++i) {
		if ( pp3.at(i) != NULL ) {
			nvertices += pp3.at(i)->size();
			continue;
		}
		//not continued ?
		cerr << "Input datasets pp3 bucket " << i << " at least is empty" << "\n";
		return false;
	}

	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT " << dscr << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";
		for(size_t b = 0; b < pp3.size(); ++b) {
			if( pp3.at(b) != NULL) {
				vector<p3dm2>* these = pp3.at(b);
				for(size_t i = 0; i < these->size(); ++i) {
					vtk << (*these)[i].x << " " << (*these)[i].y << " " << (*these)[i].z << "\n";
				}
			}
		}
		vtk << "\n";
		//add point identifier
		vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << 1 << " " << i << "\n";
		}
		//add marks as field data for coloring i.e. in Paraview

		//##MK::check formating of field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 2\n";
		vtk << "BVHBinID 1 " << nvertices << " int\n";
		for(size_t b = 0; b < pp3.size(); ++b) {
			if( pp3.at(b) != NULL) {
				vector<p3dm2>* these = pp3.at(b);
				for(size_t i = 0; i < these->size(); ++i) {
					vtk << b << "\n"; //(*these)[i] << "\n";
				}
			}
		}
		vtk << "\n";

		vtk << "UIPMark 1 " << nvertices << " int\n";
		for(size_t b = 0; b < pp3.size(); ++b) {
			if( pp3.at(b) != NULL) {
				vector<p3dm2>* these = pp3.at(b);
				for(size_t i = 0; i < these->size(); ++i) {
					vtk << these->at(i).m1 << "\n";
				}
			}
		}
		vtk << "\n";

		vtk.flush();
		vtk.close();
		return true;
	}
	else {
		cerr << "Unable to write RVE27 points and binIDs to VTK file" << "\n";
		return false;
	}
}


bool vtk_p3dm3( 	vector<vector<p3d>*> const & xyz,
							vector<vector<unsigned int>*> const & m,
								vector<vector<unsigned char>*> const & img,
									const unsigned int increment )
{
	//writes marks and 3d positions of integration points to VTK file

	//remember that the marked points are organized in buckets with each potentially dissimilar number of p3dm3 objects...
	//hence check first whether all spatial bin aka buckets are consistent for all labels
	if (m.size() == xyz.size() && img.size() == xyz.size() ) { //write consistent datasets only

		//MK::a key impracticality of the VTK file is that the number of geometrical entities has to be known a priori
		//as otherwise one needs to overwrite the nvertices length line after writing to file, which is inconvenient
		//count now how many p3dm3 objects in total contained in the bucket ensemble
		unsigned int np[3] = {0, 0, 0};
		size_t ns = 0;
		for ( size_t b = 0; b < m.size(); ++b ) {
			ns = xyz.at(b)->size();
			if ( ns == m.at(b)->size() && ns == img.at(b)->size() ) {
				np[0] += xyz.at(b)->size();		//spatial coordinates
				np[1] += m.at(b)->size();		//integration point ID
				np[2] += img.at(b)->size();		//periodicImage point ID
			}
			else {
				//##MK::DEBUG
				cerr << "Bin length " << b << " is inconsistent!" << "\n";
				return false;
			}
		} //check consistency for all bins

		//write to file
		string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) +
				".Incr." + to_string(increment) + ".IPGridBinning.vtk";

		ofstream vtk;
		vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
		if ( vtk.is_open() == true ) {
			size_t nvertices = np[0];
			//construct header and point coordinates
			vtk << "# vtk DataFile Version 2.0\n";
			vtk << "DAMASKPDT IntegrationPointGrid to BVH spatial bin assignment\n";
			vtk << "ASCII\n";
			vtk << "DATASET POLYDATA\n";
			vtk << "\n";
			vtk << "POINTS " << nvertices << " double\n";
			for (size_t b = 0; b < m.size(); ++b) {
				for ( size_t i = 0; i < xyz.at(b)->size(); ++i ) {
					vtk << xyz.at(b)->at(i).x << " " << xyz.at(b)->at(i).y << " " << xyz.at(b)->at(i).z << "\n";
				}
			}
			vtk << "\n";
			//add point identifier
			vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
			size_t j = 0;
			for (size_t b = 0; b < m.size(); ++b) {
				for ( size_t i = 0; i < m.at(b)->size(); ++i ) {
					vtk << 1 << " " << j << "\n";
					j++;
				}
			}
			//add marks as field data for coloring i.e. in Paraview

			//##MK::check formating of field data

			vtk << "POINT_DATA " << nvertices << "\n";
			vtk << "FIELD FieldData 3\n";
			vtk << "BVH_BINID 1 " << nvertices << " int\n"; //double\n";
			for (size_t b = 0; b < m.size(); ++b) {
				for ( size_t i = 0; i < m.at(b)->size(); ++i ) {
					vtk << b << "\n";
				}
			}
			vtk << "\n";

			vtk << "IP_GLOBALID 1 " << nvertices << " int\n"; //##MK::coming from uint...
			for (size_t b = 0; b < m.size(); ++b) {
				for ( size_t i = 0; i < m.at(b)->size(); ++i ) {
					vtk << static_cast<int>(m.at(b)->at(i)) << "\n";
				}
			}
			vtk << "\n";

			vtk << "BVH_MOOREPERIMG 1 " << nvertices << " int\n";
			for (size_t b = 0; b < m.size(); ++b) {
				for ( size_t i = 0; i < m.at(b)->size(); ++i ) {
					vtk << static_cast<int>(img.at(b)->at(i)) << "\n";
				}
			}

			//vtk << "\n";
			vtk.flush();
			vtk.close();
			return true;
		}
		else {
			cerr << "Unable to write integration point positions to VTK file" << "\n";
			return false;
		}
	}
	else {
		cerr << "Input data containers have dissimilar number of points or marks" << "\n";
		return false;
	}
}


bool vtk_nbhd3d( vector<p3d> const & p,
					vector<unsigned int> const & pid,
						vector<real_xyz> const & dist,
							vector<real_ori> const & theta,
								vector<unsigned int> const & textureid,
									const unsigned int increment)
{
	//writes a local higher-order neighborhood of integration points (ip) about p.at(0) and their metadata,
	//distance to the point, disorientation angle to central ip, texture index, converged DAMASK spectral increment

	unsigned int nvertices = p.size();
	if (p.size() == pid.size() && p.size() == dist.size()
			&& p.size() == theta.size() && p.size() == textureid.size() ) { //write consistent datasets only
		//MK::a key impracticality of the VTK file is that the number of geometrical entities has to be known a priori

		//write to file
		string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) +
				".Incr." + to_string(increment) + ".HOrderForIP." + to_string(pid.at(0)) + ".vtk";

		ofstream vtk;
		vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
		if ( vtk.is_open() == true ) {
			//construct header and point coordinates
			vtk << "# vtk DataFile Version 2.0\n";
			vtk << "DAMASKPDT HigherOrder neighborhood for integration point " << to_string(pid.at(0)) << "\n";
			vtk << "ASCII\n";
			vtk << "DATASET POLYDATA\n";
			vtk << "\n";

			vtk << "PointsNonPeriodicCoordinates " << nvertices << " double\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << p.at(i).x << " " << p.at(i).y << " " << p.at(i).z << "\n";
			}
			vtk << "\n";
			//add point identifier
			vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
			for ( size_t j = 0; j < nvertices; j++ ) {
				vtk << 1 << " " << j << "\n";
			}
			//add marks as field data for coloring i.e. in Paraview

			//##MK::check formating of field data

			vtk << "POINT_DATA " << nvertices << "\n";
			vtk << "FIELD FieldData 4\n";
			vtk << "IntegrationPointID 1 " << nvertices << " int\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << pid.at(i) << "\n";
			}
			vtk << "\n";

			vtk << "Distance 1 " << nvertices << " double\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << dist.at(i) << "\n";
			}
			vtk << "\n";

			vtk << "DisoriAngle 1 " << nvertices << " double\n";
			for (size_t i = 0; i < nvertices; ++i) {
				if ( isnan(theta.at(i)) == false )
					vtk << theta.at(i) << "\n";
				else
					vtk << 360.0 << "\n"; //MK::prevent corrupting VTK file data formats
			}
			vtk << "\n";

			vtk << "TextureID 1 " << nvertices << " int\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << textureid.at(i) << "\n";
			}
			vtk << "\n";

			//vtk << "\n";
			vtk.flush();
			vtk.close();
			return true;
		}
		else {
			cerr << "Unable to write higher order neighbor environment to VTK file" << "\n";
			return false;
		}
	}
	else {
		cerr << "Input data containers have dissimilar number of entries" << "\n";
		return false;
	}
}


bool vtk_p3dm3( vector<p3d> const & p,
					vector<unsigned int> const & textureid,
						const unsigned int increment)
{
	//writes integration point grid with one scalar value uint32 per point

	unsigned int nvertices = p.size();
	if (p.size() == textureid.size() ) { //write consistent datasets only
		//MK::a key impracticality of the VTK file is that the number of geometrical entities has to be known a priori

		//write to file
		string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) +
				".Incr." + to_string(increment) + ".IPGridTextureID.vtk";

		ofstream vtk;
		vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
		if ( vtk.is_open() == true ) {
			//construct header and point coordinates
			vtk << "# vtk DataFile Version 2.0\n";
			vtk << "DAMASKPDT Integration point grid with texture ID assignment\n";
			vtk << "ASCII\n";
			vtk << "DATASET POLYDATA\n";
			vtk << "\n";

			vtk << "PointsNonPeriodicCoordinates " << nvertices << " double\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << p.at(i).x << " " << p.at(i).y << " " << p.at(i).z << "\n";
			}
			vtk << "\n";
			//add point identifier
			vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
			for ( size_t j = 0; j < nvertices; j++ ) {
				vtk << 1 << " " << j << "\n";
			}
			//add marks as field data for coloring i.e. in Paraview

			//##MK::check formating of field data

			vtk << "POINT_DATA " << nvertices << "\n";
			vtk << "FIELD FieldData 1\n";
			vtk << "TextureID 1 " << nvertices << " int\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << textureid.at(i) << "\n";
			}
			vtk << "\n";

			//vtk << "\n";
			vtk.flush();
			vtk.close();
			return true;
		}
		else {
			cerr << "Unable to write integration point grid texture ID to VTK file" << "\n";
			return false;
		}
	}
	else {
		cerr << "Input data containers have dissimilar number of entries" << "\n";
		return false;
	}
}


bool vtk_gp3d( vector<p3d> const & p,
					vector<unsigned int> const & grainid,
						const unsigned int increment)
{
	//writes integration point grid with one scalar value uint32, which is the grain id per point

	unsigned int nvertices = p.size();
	if (p.size() == grainid.size() ) { //write consistent datasets only
		//MK::a key impracticality of the VTK file is that the number of geometrical entities has to be known a priori

		//write to file
		string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) +
				".Incr." + to_string(increment) + ".IPGridGrainID.vtk";

		ofstream vtk;
		vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
		if ( vtk.is_open() == true ) {
			//construct header and point coordinates
			vtk << "# vtk DataFile Version 2.0\n";
			vtk << "DAMASKPDT Integration point grid with grain ID assignment\n";
			vtk << "ASCII\n";
			vtk << "DATASET POLYDATA\n";
			vtk << "\n";

			vtk << "PointsNonPeriodicCoordinates " << nvertices << " double\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << p.at(i).x << " " << p.at(i).y << " " << p.at(i).z << "\n";
			}
			vtk << "\n";
			//add point identifier
			vtk << "VERTICES " << nvertices << " " << 2*nvertices << "\n";
			for ( size_t j = 0; j < nvertices; j++ ) {
				vtk << 1 << " " << j << "\n";
			}
			//add marks as field data for coloring i.e. in Paraview

			//##MK::check formating of field data

			vtk << "POINT_DATA " << nvertices << "\n";
			vtk << "FIELD FieldData 1\n";
			vtk << "GrainID 1 " << nvertices << " int\n";
			for (size_t i = 0; i < nvertices; ++i) {
				vtk << grainid.at(i) << "\n";
			}
			vtk << "\n";

			//vtk << "\n";
			vtk.flush();
			vtk.close();
			return true;
		}
		else {
			cerr << "Unable to write integration point grid grain ID to VTK file" << "\n";
			return false;
		}
	}
	else {
		cerr << "Input data containers have dissimilar number of entries" << "\n";
		return false;
	}
}
