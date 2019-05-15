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


#ifndef __PDT_VOROCOMPOSER_H__
#define __PDT_VOROCOMPOSER_H__

#include "PDT_VTKIO.h"

#include "thirdparty/VoroRycroft/voro++-0.4.6/src/voro++.hh"
using namespace voro;


struct tess_domain_info
{
	voro_aabb3d bounds;
	voro_p3i blocks;
	voro_real pdensity;
	int initialmem;						//how much point object space allocated initially per Voro++ spatial bucket?
	bool periodicx;
	bool periodicy;
	bool periodicz;
	tess_domain_info() : bounds(voro_aabb3d()), blocks(voro_p3i()), pdensity(static_cast<voro_real>(0.0)),
			initialmem(1), periodicx(false), periodicy(false), periodicz(false) {}
	tess_domain_info( const voro_aabb3d _bx, const voro_p3i _bpart, const voro_real _pdens, const int _imem,
			const bool _perx, const bool _pery, const bool _perz ) :
			bounds(_bx), blocks(_bpart), pdensity(_pdens), initialmem(_imem),
			periodicx(_perx), periodicy(_pery), periodicz(_perz) {}
};

ostream& operator<<(ostream& in, tess_domain_info const & val);


struct voro_facet
{
	voro_uint v_s;		//index start
	voro_uint v_e; 		//index one past end
	voro_uint owner;	//belonging to which cell MK::currently interior_v array position reference
	voro_p3d ounormal;	//outer unit normal
	voro_p3d centrum; 	//centrum a point inside the convex polygon contour of the facet MK::Voronoi cells!
	voro_p3d inplane_u; //normal vector in the plane
	voro_p3d inplane_v; //normal vector in the plane
	voro_uint v_s_polyproj; //index start
	voro_uint v_e_polyproj; //index one past end, inplane_u inplane_v reversed order and projected contour
	voro_facet() : v_s(VOROUINTXX), v_e(VOROUINTXX), owner(VOROUINTXX), ounormal(voro_p3d(1.0,0.0,0.0)), centrum(voro_p3d(0.0,0.0,0.0)),
			inplane_u(voro_p3d(0.0,0.0,0.0)), inplane_v(voro_p3d(0.0,0.0,0.0)), v_s_polyproj(VOROUINTXX), v_e_polyproj(VOROUINTXX) {}
	voro_facet( const voro_uint _ns, const voro_uint _ne, const voro_uint _nown, const voro_p3d _centr ) :
		v_s(_ns), v_e(_ne), owner(_nown), ounormal(voro_p3d(0.0,0.0,0.0)), centrum(_centr),
		inplane_u(voro_p3d(0.0,0.0,0.0)), inplane_v(voro_p3d(0.0,0.0,0.0)), v_s_polyproj(VOROUINTXX), v_e_polyproj(VOROUINTXX) {}
};


class microstructural_object
{
public:
	microstructural_object();
	~microstructural_object();

	bool construct_object( vector<voro_p3dm2> const & in, const voro_uint refpid );
	void identify_exterior_facets( const unsigned int clpid, const voro_uint refpid,
			vector<int> const & nbors, vector<int> const & fvert, vector<double> const & verts );
	void build_facet_bvh();
	void build_facet_normals_and_reverse_contours();
/*
	void identify_normaldistances();  	//incorrect geometrical shadows not considered and unstable numerics
	void identify_normaldistances2();	//shadows considered as well as edge and vertex checks but all against all
*/
	void identify_normaldistances3( const voro_real threshold ); 	//like 2 but aabb pruning, using prereversed and projected contour
	bool visualize_object_hull( const unsigned int simid, const unsigned int incr, const unsigned int gid );
	bool visualize_distances( const unsigned int simid, const unsigned int incr, const unsigned int gid );

	squat ori;
	tess_domain_info tess_info;
	vector<voro_uint> cid2pid;		//mapping local Voronoi cell ID to point object ID
	vector<voro_uint> cid2uip;		//mapping local Voronoi cell ID to unique material point ID
	vector<voro_p3d> interior_v;	//an interior support point of the complex at which there is a Voronoi cell
	vector<voro_real> interior_d;	//normal distance to exterior hull
	vector<unsigned int> interior_stats;	//how many facets to test
	vector<unsigned int> interior_uip;		//reference to a unique material point ID
	vector<voro_facet> exterior_f;	//information about exterior facets
	vector<voro_uint> exterior_n;	//vertex references how the individual facet polygons can be constructed from exterior_v
	vector<voro_p3d> exterior_v;	//vertices building the exterior facets

	vector<voro_p3d> exterior_v_aligned; //possible duplication of vertices in exterior_v but laid out linearly in cache according to exterior_f.at(i)->v_s, v_e

	vector<point_xy> clockwise_projv;	//clockwise sorted projected contour segments linear in cache according to exterior_f.at(i)->v_s_polyproj, v_e_polyproj

	//vector<objecthull> hull;
	Tree* bvh;						//an (axis-aligned bounding box) bounded volume hierarchy on the facet to reduce distancing querying costs
};


class voroComposer
{
public:
	voroComposer();
	voroComposer( const voro_uint nobjects );
	~voroComposer();

	vector<microstructural_object*> obj;
};


#endif
