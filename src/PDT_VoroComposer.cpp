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



#include "PDT_VoroComposer.h"


ostream& operator<<(ostream& in, tess_domain_info const & val)
{
	in << "\n";
	in << "Tessellation domain" << "\n";
	in << val.bounds;
	in << val.blocks;
	in << "PDensity " << val.pdensity << "\n";
	in << "PeriodicX " << val.periodicx << "\n";
	in << "PeriodicY " << val.periodicy << "\n";
	in << "PeriodicZ " << val.periodicz << "\n";
	in << "InitAlloc " << val.initialmem << "\n";
	in << "\n";
	return in;
}


microstructural_object::microstructural_object()
{
	ori = squat();
	tess_info = tess_domain_info();
	bvh = NULL;
}


microstructural_object::~microstructural_object()
{
	delete bvh;
}


bool microstructural_object::construct_object( vector<voro_p3dm2> const & in, const voro_uint refpid )
{
	if ( in.size() > VOROINTMX ) {
		#pragma omp critical
		{
			cerr << "Too many points in the input point cloud, object is not constructable!" << "\n";
		}
		return false;
	}

	ori = squat();
	tess_info = tess_domain_info();

	//get extremal bounds from the point cloud
	voro_aabb3d aabb = voro_aabb3d();
	for( auto it = in.begin(); it != in.end(); ++it ) { //bounds are changed infrequently
		aabb.xmi = ( it->x > aabb.xmi ) ? aabb.xmi : it->x;
		aabb.xmx = ( it->x < aabb.xmx ) ? aabb.xmx : it->x;
		aabb.ymi = ( it->y > aabb.ymi ) ? aabb.ymi : it->y;
		aabb.ymx = ( it->y < aabb.ymx ) ? aabb.ymx : it->y;
		aabb.zmi = ( it->z > aabb.zmi ) ? aabb.zmi : it->z;
		aabb.zmx = ( it->z < aabb.zmx ) ? aabb.zmx : it->z;
	}
	aabb.add_epsilon_guard( EPSILON );
	tess_info.bounds = aabb;
	tess_info.periodicx = false;
	tess_info.periodicy = false;
	tess_info.periodicz = false;
	tess_info.pdensity = 5.0; //based on empirical result by Rycroft, http://math.lbl.gov/voro++/examples/timing_test/
	tess_info.initialmem = 1;
	tess_info.blocks = aabb.blockpartitioning( in.size(), tess_info.pdensity );
	//too few blocks too many points to test against per block
	//too many blocks too much memory overhead and cache trashing (L1,L2 and not to forget "page" TLB cache)

	#pragma omp critical
	{
		cout << tess_info << "\n";
	}

	//build Voronoi tessellation of point cloud in the local tessellation domain aabb
	double mytic = omp_get_wtime();

	//Voro++ tessellation domain
	container con(	tess_info.bounds.xmi, tess_info.bounds.xmx,
					tess_info.bounds.ymi, tess_info.bounds.ymx,
					tess_info.bounds.zmi, tess_info.bounds.zmx,
					tess_info.blocks.x, tess_info.blocks.y, tess_info.blocks.z,
					tess_info.periodicx, tess_info.periodicy, tess_info.periodicz,
					tess_info.initialmem );

	//##MK::Voro++ use int, only when start at 0 here, we can use cl.pid() to reference the particleID via an index on cid2pid array
	cid2pid.reserve( in.size() );
	cid2uip.reserve( in.size() );
	int ni = static_cast<int>(in.size());
	for ( int i = 0; i < ni; i++ ) {
		//##MK::given that Voro++ uses a spatial bucket the cells are not necessarily processed in the same order than listed in input array in!
		con.put( i, in[i].x, in[i].y, in[i].z ); //Voro++ dictates to use int and double
		cid2pid.push_back( in[i].m1 ); //MK::pass m1 the grain ID, not m2 the unique material point ID !
		cid2uip.push_back( in[i].m2 ); //MK::pass m2
	}

	#pragma omp critical
	{
		cout << "Number of particles in con " << con.total_particles() << "\n";
	}

	c_loop_all cl(con);
	voronoicell_neighbor c;

	if (cl.start()) { //this is incremental computing of the tessellation, we holds at no point the entire tessellation in memory
		do {
			voro_uint current_pid = cid2pid.at(cl.pid()); //which particleID for this cell
			voro_uint current_uip = cid2uip.at(cl.pid()); //##MK::fuse together with cid2pid
			if ( current_pid != refpid ) { //we build the Voronoi complex only from those Voronoi cells of points with a specific particleID
				continue;
			}
			else { //a cell belonging to the complex, identify its facets and whether they support the inside of the complex
				//constitute an exterior facet
				if ( con.compute_cell(c, cl) == false ) {
					#pragma omp critical
					{
						cerr << "Tessellation of a cell for object " << refpid << " failed" << "\n";
					}
					return false;
				}
				else { //constructing this cell with neighborhood information was successfully

					vector<int> first_neigh;
					c.neighbors(first_neigh);

					//was the guard zone about the complex chosen large enough to have no cell within making domain wall contact?
					bool domainwall_contact = false;
					for( auto jt = first_neigh.begin(); jt != first_neigh.end(); ++jt ) {
						if ( *jt >= 0 ) { //Voro++ internally assigns negative neighbor IDs from -1, -6 to identify the domain walls
							continue;
						}
						else {
							domainwall_contact = true;
							break;
						} //no need to test further neighbors
					}

					if ( domainwall_contact == true ) {
						#pragma omp critical
						{
							cerr << "Tessellation of object " << refpid << " was successful but local domain wall contact" << "\n";
						}
						return false;
					}
					else {
						//cell definately supports the volume of the complex
						double x, y, z;
						cl.pos(x, y, z);

						interior_v.push_back( voro_p3d( static_cast<voro_real>(x),
								static_cast<voro_real>(y), static_cast<voro_real>(z) ) );

						interior_uip.push_back( current_uip );

						//identify whether the cell has contact with the complex exterior boundary
						//MK::for cells without domain wall contact this is the case only if at least one neighbor cell has
						//a particleID which is not the refpid
						bool hull_contact = false;
						for( auto jt = first_neigh.begin(); jt != first_neigh.end(); ++jt ) {
							//inside all *jt >= 0
							voro_uint current_pid = cid2pid.at(*jt); //##MK::cache miss? most cells are inside
							if ( current_pid == refpid ) { //the more likely case as most cells are inside for a large complex
								continue;
							}
							else {
								hull_contact = true;
								break;
							} //no need to test further there is at least one exterior facet
						}

						if ( hull_contact == true ) { //MK::store all exterior facets of this cell wrt to world Cartesian coordinate system
							vector<int> first_f_vert;
							c.face_vertices(first_f_vert);
							vector<double> first_v;
							c.vertices(x, y, z, first_v); //right hand rule wrt outer unit normal, anticlockwise

							//unsigned int clpid = static_cast<unsigned int>(cl.pid());
							unsigned int clpid = static_cast<unsigned int>(interior_v.size()-1);

							identify_exterior_facets( clpid, refpid, first_neigh, first_f_vert, first_v );
						}
					} //done with processing hull contribution
				} //done with cell construction
			} //done with processing the geometrical contribution to the exterior hull of the current cell to the target complex
		} while (cl.inc());
	}

	double mytoc = omp_get_wtime();

	#pragma omp critical
	{
		cout << "Processing the exterior hull completed" << "\n";
		cout << "Complex hull exterior facets " << exterior_f.size() << "\n";
		cout << "Complex hull exterior vertices " << exterior_v.size() << "\n";
		cout << "Complex interior support " << interior_v.size() << "\n";
		cout << "Computing the complex took " << (mytoc-mytic) << " seconds sequentially" << "\n";
	}
	return true;
}


void microstructural_object::identify_exterior_facets( const unsigned int clpid, const voro_uint refpid,
		vector<int> const & nbors, vector<int> const & fvert, vector<double> const & verts )
{
	//when reporting the cell facet vertices in the global coordinate system the Voro++ library reports
	//a list of all the unique 1d coordinate values not the coordinate value triplets!
	//for visualization of the cell geometry with at least to my knowledge all tools, though coordinate value triplets are required.

	//this is only a cell-local ID remapping to avoid duplicates but for a single Voronoi cell
	//however now we patch together facets from multiple Voronoi cells into the polygon hull of the microstructural_object
	//a Voronoi cell may contribute multiple facets to the complex
	voro_uint exterior_v_offset = exterior_v.size();

	map<voro_uint,voro_uint> old2new;
	vector<voro_uint> fvert_reindexed; //references to positions on exterior_v !

	//because Voro++ reports vertices as a set of the unique 1d coordinate values
	//keys of the maps are build from the old verts indices where int indices u,v,w are converted into implicit 1d key
	//values of the map are the new indices referring to p3d verts_triplets specific to only the current cell
	//for subsequent I/O these packets of p3d triplets are fused together in a large array using an implicit position offset on the 1d array
	voro_uint ii = 0;
	int j = 0; //does not need to be so large because only single Voronoi cell, however complex may be build from millions of facets...
	for( auto it = nbors.begin(); it != nbors.end(); ++it ) { //O(N) processing
		//only facets to neighboring cells with particleID different than refpid
		voro_uint current_pid = cid2pid.at(*it);
		if ( current_pid != refpid ) { //a facet bounding a voronoi cell about a particle with a different ID can be an exterior one

			int nvertices_thisfacet = fvert[j]; //##MK::should be voro_int
			if ( nvertices_thisfacet > 2 ) { //do not add degenerate facets

				voro_p3d centr = voro_p3d( 0.0, 0.0, 0.0 );

				voro_uint fp_past = exterior_n.size(); //append remapped vertex IDs from cell to cell complex

				for( int k = 0; k < nvertices_thisfacet; k++) { //generating 3d points with implicit indices 0, 1, 2, ....
					voro_uint l = static_cast<voro_uint>(3*fvert[j+k+1]);
					//hashing combinations of x,y,z
					voro_uint implicitkey = l + (l+1)*VORO_COMPOSER_VERTEX_MAPPING_HASHVAL +
							(l+2)*VORO_COMPOSER_VERTEX_MAPPING_HASHVAL*VORO_COMPOSER_VERTEX_MAPPING_HASHVAL; //##MK::robustness?
					auto jt = old2new.find( implicitkey );
					if ( jt != old2new.end() ) { //O(lg(N)) found integer triplet exists already
						//fvert_reindexed.push_back( exterior_v_offset + jt->second );
						exterior_n.push_back( exterior_v_offset + jt->second ); //reference already existent locally remapped vertex and map to global wrt to complex
	//cout << jt->second << "\n";
					}
					else { //integer triplet does not yet exist so create
						old2new.insert( make_pair( implicitkey, ii ) );
						//fvert_reindexed.push_back( exterior_v_offset + ii ); //mapping from locally ID remapped to complex globally ID remapped
						exterior_n.push_back( exterior_v_offset + ii ); //reference new locally remapped vertex and map to global
	//cout << ii << "\n";
						++ii;
						//add this unique locally remapped vertex to the complex
						exterior_v.push_back( voro_p3d( static_cast<voro_real>(verts[l]),
									static_cast<voro_real>(verts[l+1]), static_cast<voro_real>(verts[l+2]) )  );
					}

					//layout plssobily even duplicates in memory to avoid frequent cache misses when probing facet polygon contour distances
					exterior_v_aligned.push_back( voro_p3d( static_cast<voro_real>(verts[l]),
									static_cast<voro_real>(verts[l+1]), static_cast<voro_real>(verts[l+2]) )  );

					centr.x += verts[l];
					centr.y += verts[l+1];
					centr.z += verts[l+2];
				}
				voro_uint fp_now = exterior_n.size();

				//Voro++ stores facets in anticlockwise order but we need them in clockwise order for boost
				//##MK::reorder them already here

				//register a new non-singular at least triangular exterior facet
				centr.x /= static_cast<voro_real>(nvertices_thisfacet);
				centr.y /= static_cast<voro_real>(nvertices_thisfacet);
				centr.z /= static_cast<voro_real>(nvertices_thisfacet);

				exterior_f.push_back( voro_facet( fp_past, fp_now, static_cast<voro_uint>(clpid), centr ) );
				/*
				if ( (fp_now - fp_past) > 2 ) { //so new points were added, ie. non singular facet added
					exterior_f.push_back( voro_facet( fp_past, fp_now, static_cast<voro_uint>(clpid) ) );
				}
				else {
					exterior_f.push_back( voro_facet( fp_past, fp_now, static_cast<voro_uint>(clpid) ) );
					cerr << "Attempting to add a singular facet" << "\n";
				}*/
			}
			else {
				#pragma omp critical
				{
					cerr << "Attempting to add a singular facet" << "\n";
				}
			}
		} //done with current exterior facet, check if next facet is also an exterior
		//MK::in every case, especially also after discarding the facet advance facet vertex counter, all facets need processing because all can in principle contribute to hull
		j += fvert[j] + 1;
	} //done with all facets of the Voronoi cell


	/*
	//second step: pass heavy data to threadlocal buffer
	//##MK::could at some point be fused with first step to save loop overhead...
	ii = 0;
	j = 0;
	for( size_t i = 0; i < nbors.size(); ++i ) { //O(N)lg(N) identification of how many faces of certain type
		unsigned int nvertices_thisface = fvert[j];
		io_topo.push_back( 3 ); //XDMF keyword to visualize an n-polygon
		io_topo.push_back( nvertices_thisface ); //promotion uint23 to size_t no problem
		//for( int k = 0; k < fvert[j]; k++) { //it is essential to sweep through this loop nest in the same order than above for carrying over the prescribed order on fvert_reindexed
		for( unsigned int k = 0; k < nvertices_thisface; k++, ii++ ) {
			io_topo.push_back( fvert_reindexed.at(ii) ); //promotion of int to size_t not a problem if int >= 0
		}
		j += fvert[j] + 1;
	}

	for( auto it = vtriplets_unique.begin(); it != vtriplets_unique.end(); ++it ) {
		io_geom.push_back( it->x );
		io_geom.push_back( it->y );
		io_geom.push_back( it->z );
	}
	*/
}


void microstructural_object::build_facet_bvh()
{
	try {
		bvh = new Tree;
	}
	catch (bad_alloc &croak) {
		cerr << "Unable to allocate BVH tree for upon microstructural object constructor" << "\n";
		return;
	}

	for ( size_t ft = 0; ft < exterior_f.size(); ++ft ) {
		voro_aabb3d ftbox = voro_aabb3d();
		size_t vs = exterior_f[ft].v_s;
		size_t ve = exterior_f[ft].v_e;
		for ( size_t i = vs; i < ve; ++i ) {
			voro_p3d here = exterior_v_aligned[i];
			ftbox.possibly_enlarge_me( here );
		}

		bvh->insertParticle( static_cast<unsigned int>(ft),
				trpl(ftbox.xmi, ftbox.ymi, ftbox.zmi), trpl(ftbox.xmx, ftbox.ymx, ftbox.zmx) );
	}
}


void microstructural_object::build_facet_normals_and_reverse_contours()
{
	double mytic = omp_get_wtime();

	//compute outer unit normals
	for( auto ft = exterior_f.begin(); ft != exterior_f.end(); ++ft ) {
		size_t prev = ft->v_e - 1;
		size_t here = ft->v_s;
		size_t next = ft->v_s + 1;
		/*voro_p3d prevp3d = exterior_v.at(exterior_n.at(prev)); //##MK::
		voro_p3d herep3d = exterior_v.at(exterior_n.at(here));
		voro_p3d nextp3d = exterior_v.at(exterior_n.at(next));*/
		voro_p3d prevp3d = exterior_v_aligned[prev];
		voro_p3d herep3d = exterior_v_aligned[here];
		voro_p3d nextp3d = exterior_v_aligned[next];
		//take to direction vectors in the facet plane
		voro_p3d next_here = voro_p3d( nextp3d.x - herep3d.x, nextp3d.y - herep3d.y, nextp3d.z - herep3d.z );
		voro_p3d prev_here = voro_p3d( prevp3d.x - herep3d.x, prevp3d.y - herep3d.y, prevp3d.z - herep3d.z );
		//cross product is a normal to that plane if the facet polygon is not degenerated
		//the Voro++ has assured this
		voro_p3d normal = cross( next_here, prev_here );
		normal.normalize();
		//this can be either the inner or the outer unit normal now, so we have to enforce consistency next

		if ( isnan(normal.x) == true || isnan(normal.y) == true || isnan(normal.z) ) {
			#pragma omp critical
			{
				cerr << "ounormal is nan " << normal << "\n";
			}
		}
		//now take a test point on a known "side" of the facet, ie. the support of the Voronoi cell
		//Voronoi cel is convex
		voro_p3d cellcenter = interior_v.at(ft->owner);
		voro_p3d cellcenter_here = voro_p3d( cellcenter.x - herep3d.x, cellcenter.y - herep3d.y, cellcenter.z - herep3d.z );

		//##MK::do not normalize ! cellcenter_here.normalize();

		//projection on normal if positive on the same side, given that cellcenter is inside this will be the inner normal
		voro_real projection = dot( cellcenter_here, normal );

		//##MK::Voro++ returns facet contour polygon vertices in anticlockwise order wrt to outer unit normal
		if ( projection > static_cast<voro_real>(0.0) ) { //positive on the same side than is the cellcenter so for a Voronoi cell this is
			//candidate normal was the inner unit normal, so flip signs
			ft->ounormal = voro_p3d( -normal.x, -normal.y, -normal.z );
		}
		else { //normal is sought after outer unit normal
			ft->ounormal = voro_p3d( normal.x, normal.y, normal.z );
		}

/*
		//compute in plane normal vector u
		//take inplane_u projection vector using the largest projected world coordinate vector
		//x - dot(x,oun)*oun
		int largest = 0;
		voro_real val = fabs(ft->ounormal.x);
		if ( fabs(ft->ounormal.y) >= val ) {
			largest = 1;
			val = fabs(ft->ounormal.y);
		}
		if ( fabs(ft->ounormal.z) >= val ) {
			largest = 2;
			val = fabs(ft->ounormal.z);
		}

		voro_p3d wx = voro_p3d( 1.0, 0.0, 0.0 );
		if ( largest == 1 )
			wx = voro_p3d( 0.0, 1.0, 0.0 );
		if ( largest == 2 )
			wx = voro_p3d( 0.0, 0.0, 1.0 );
		voro_real d = dot( wx, ft->ounormal );

		ft->inplane_u = voro_p3d( 	wx.x - d*ft->ounormal.x,
									wx.y - d*ft->ounormal.y,
									wx.z - d*ft->ounormal.z   );
		ft->inplane_u.normalize();
*/

		ft->inplane_u = next_here;
		ft->inplane_u.normalize();

		voro_real dot_u_n = dot( ft->inplane_u, ft->ounormal );
		if ( dot_u_n > EPSILON ) {
			#pragma omp critical
			{
				cerr << "dot_u_n > EPS " << dot_u_n << "\n";
				cerr << "ounormal " << ft->ounormal << "\n";
				cerr << "inplane_u" << ft->inplane_u << "\n";
			}
		}

		if ( isnan(ft->inplane_u.x) == true || isnan(ft->inplane_u.y) == true || isnan(ft->inplane_u.z) ) {
			#pragma omp critical
			{
				cerr << "inplane u is nan " << ft->inplane_u << "\n";
				cerr << "ounormal " << ft->ounormal << "\n";
			}
/*
			cerr << "wx " << wx << "\n";
			cerr << "d " << d << "\n";
			cerr << "cellcenter_here " << cellcenter_here << "\n";
*/
		}


		//second in plane normal vector v
		ft->inplane_v = cross( ft->inplane_u, ft->ounormal );
		ft->inplane_v.normalize();

		if ( isnan(ft->inplane_v.x) == true || isnan(ft->inplane_v.y) == true || isnan(ft->inplane_v.z) ) {
			#pragma omp critical
			{
				cerr << "inplane_v is nan " << ft->inplane_v << "\n";
				cerr << "ounormal " << ft->ounormal << "\n";
			}
/*
	cerr << "wx " << wx << "\n";
	cerr << "d " << d << "\n";
	cerr << "cellcenter_here " << cellcenter_here << "\n";
*/
		}

		//hyperphobic
		if ( fabs(dot(ft->inplane_u, ft->ounormal)) < EPSILON && fabs(dot(ft->inplane_v, ft->ounormal)) < EPSILON )
			continue;
		else {
			#pragma omp critical
			{
				cerr << "Facet inplane u or inplane_v not perpendicular to ounormal " <<
						fabs(dot(ft->inplane_u, ft->ounormal)) << "\t\t" << fabs(dot(ft->inplane_v, ft->ounormal)) << "\n";
			}
		}

		//reverse the ordering of the individual facet contour polygon vertices
		//from anticlockwise wrt to ounormal to clockwise
		voro_p3d uu = ft->inplane_u;
		voro_p3d vv = ft->inplane_v;

		voro_uint fppast = clockwise_projv.size();
		//and project them to a Cartesian coordinate system in the individual facet plane for boost covered_by tested point in/on 2d polygon contour
		for( voro_uint i = ft->v_e-1; i > ft->v_s; i-- ) {
			voro_p3d thisone = exterior_v_aligned[i];
			clockwise_projv.push_back( point_xy(dot(thisone, uu), dot(thisone, vv)) );
		}
		//dont forget to close the contour by adding the first point
		voro_p3d thisone = exterior_v_aligned[ft->v_e-1];
		clockwise_projv.push_back( point_xy(dot(thisone,uu), dot(thisone,vv)) );

		voro_uint fpnow = clockwise_projv.size();

		ft->v_s_polyproj = fppast;
		ft->v_e_polyproj = fpnow;
	}

	double mytoc = omp_get_wtime();
	/*
	#pragma omp critical
	{
		cout << "Computing normals and sorting reversed projected contours took " << (mytoc-mytic) << " seconds sequentially" << "\n";
	}
	*/
}


/*
void microstructural_object::identify_normaldistances()
{
	double mytic = omp_get_wtime();

	for( auto pt = interior_v.begin(); pt != interior_v.end(); ++pt ) {
		voro_real shortest_sqr_normal_distance = VOROFMX;

		for( auto ft = exterior_f.begin(); ft != exterior_f.end(); ++ft ) {
			//shot ray from pt along ounormal of facet get if possible the intersection with the plane
			//l0 interior point currently l0
			//l ray directions ft->ounormal
			//t scaling such that (l0 + l*t - p0)n = 0
			//p0 point on the facet
			voro_p3d nou = ft->ounormal;
			voro_p3d cen = ft->centrum;
			voro_real denom = dot( nou, nou );
			if ( denom > EPSILON ) {  //##MK::given that nou is always a normal vector there is always a possible plane so no if necessary!
				voro_p3d l0 = *pt;
				voro_p3d p0 = exterior_v.at(exterior_n.at(ft->v_s)); //##MK::two likely cache misses could be reduced if for each facet storing the point
				//voro_p3d p0l0 = voro_p3d( p0.x - l0.x, p0.y - l0.y, p0.z - l0.z );
				voro_p3d l0p0 = voro_p3d( l0.x-p0.x, l0.y-p0.y, l0.z-p0.z );


				voro_real t = dot( l0p0, nou ) / denom;
				voro_p3d hitpoint = voro_p3d( pt->x + t*nou.x, pt->y + t*nou.y, pt->z + t*nou.z);

				//MK::Voronoi cells! are concex, complex of Voronoi cells is not, but facets of a Voronoi cell are convex polygons!
				bool hitpoint_inside = true; //attempt to reject
				for( size_t i = ft->v_s; i < ft->v_e; i++ ) {
					//MK::only triangular facets!
					//size_t prev = ( i > ft->v_s ) ? ft->v_s - 1 : ft->v_e - 1;
					size_t here = i;
					voro_p3d herep3d = exterior_v.at(exterior_n.at(here));
					voro_p3d nextp3d = voro_p3d();
					voro_p3d nexthere = voro_p3d();
					bool wraparound = ( (i+1) < (ft->v_e - 1) ) ? false : true;
					if ( wraparound == false ) {
						size_t next = i+1;
						nextp3d = exterior_v.at(exterior_n.at(next));
						nexthere = voro_p3d( nextp3d.x - herep3d.x,  nextp3d.y - herep3d.y,  nextp3d.z - herep3d.z );
					}
					else {
						size_t next = ft->v_e - 1;
						nextp3d = exterior_v.at(exterior_n.at(next));
						//MK::mind wraparound! change order to have nexthere still pointing in circulation direction!
						nexthere = voro_p3d( herep3d.x - nextp3d.x, herep3d.y - nextp3d.y, herep3d.z - nextp3d.z );
					}

					voro_p3d ft_contour_segmentnormal = cross( nou, nexthere );
					ft_contour_segmentnormal.normalize();

					//compute outer unit normal to facet contour segment normal
					voro_p3d cen_here = voro_p3d( cen.x - herep3d.x, cen.y - herep3d.y, cen.z - herep3d.z );
					//###MK::umlaufsinn about cen assuming here anticlockwise!
					voro_real sign = dot( ft_contour_segmentnormal, cen_here );
					if ( sign >= 0.0 ) { //##MK::inconsistent outer unit normal to the contour segment with respect to umlaufsinn and cen
						//flip sign
						ft_contour_segmentnormal.x = -ft_contour_segmentnormal.x;
						ft_contour_segmentnormal.y = -ft_contour_segmentnormal.y;
						ft_contour_segmentnormal.z = -ft_contour_segmentnormal.z;
					}
					//is hitshere "behind" infinite line in the facet plane to which such contour segment could be imaginarily extended
					//and hence indicating that the contour segment in fact circulate (anticlockwise) about the point hitshere?
					voro_p3d hitpoint_here = voro_p3d( hitpoint.x - herep3d.x, hitpoint.y - herep3d.y, hitpoint.z - herep3d.z );
					voro_real sign_decision = dot( hitpoint_here, ft_contour_segmentnormal );
					if ( sign_decision > 0.0 ) { //if == 0.0 then on the line if negative with respect to outer contour segment normal vector then behind so no reject
						hitpoint_inside = false;
						//assumption that hitpoint is inside facet contour polygon has to be rejected
						break;
					}
				}

				if ( hitpoint_inside == true ) {
					voro_real dist = SQR(hitpoint.x - pt->x)+SQR(hitpoint.y - pt->y)+SQR(hitpoint.z - pt->z);
					if ( dist > shortest_sqr_normal_distance ) { //likely most will not be the closest
						continue;
					} //test next facet
					else {
						shortest_sqr_normal_distance = dist;
//cout << shortest_sqr_normal_distance << "\n";
					}
				}
			}
			else { //should not be encountered
				#pragma omp critical
				{
					cerr << "Normal distance computation dot(nou,nou) <= EPSILON!" << "\n";
				}
				return;
			}

			//get the is the intersectionpoint inside the facet polygon contour
		}

		//##MK::take square root!
		interior_d.push_back( shortest_sqr_normal_distance );

//cout << "Reporting shortest ___" << interior_d.back() << "___\n";
	}

	double mytoc = omp_get_wtime();
	cout << "Computing the distances took " << (mytoc-mytic) << " seconds sequentially" << "\n";
}


void microstructural_object::identify_normaldistances2()
{
	double mytic = omp_get_wtime();

	vector<unsigned int> iotopo;
	vector<float> ioattr;
	vector<float> iogeom;
	unsigned int faultypid = 0;

	size_t iid = 0;
	for( auto pt = interior_v.begin(); pt != interior_v.end(); ++pt, iid++ ) {
		voro_real shortest_sqr_normal_distance = VOROFMX;
		//##MK::bounding box pre-checks to prune search space

//cout << "---->" << pt->x << ";" << pt->y << ";" << pt->z << "\n";
		for( auto ft = exterior_f.begin(); ft != exterior_f.end(); ++ft ) {
			//arbitrary point on the facet plane and outer unit normal
			voro_p3d r0 = ft->centrum;
			voro_p3d nou = ft->ounormal;
			//difference vector r0-pt
			voro_p3d v3 = voro_p3d( r0.x - pt->x, r0.y - pt->y, r0.z - pt->z );
			//normal projection of pt on the facet plane
			voro_real d = dot( nou, v3 );
			voro_p3d hit = voro_p3d( 	pt->x + nou.x*d,
										pt->y + nou.y*d,
										pt->z + nou.z*d );

			//project hit point + contour polygon of the facet vertices from 3d Cartesian to 2d Cartesian and point in polygon inclusion
			//using Boost library

			//Voro++ uses right-hand rule, ie we circulate anticlockwise but Boost by default circulates clockwise!
			//so we have to reverse order the points ##MK

			vector<point_xy> contour_pts;
			voro_p3d uu = ft->inplane_u;
			voro_p3d vv = ft->inplane_v;
//cout << ft->v_s << "\t\t" << ft->v_e << "\n";

//if ( iid == 154 ) {
//	cout << "Contour how it is " << "\n";
//	cout << "r0" << r0;
//	cout << "nou" << nou;
//	cout << "v3" << v3;
//	cout << "d " << d << "\n";
//	cout << "hit " << hit;
//	cout << "uu " << uu;
//	cout << "vv " << vv;
//}

			for( size_t i = ft->v_e-1; i > ft->v_s; i-- ) {
				//voro_p3d thisone = exterior_v.at(exterior_n.at(i));

//if ( iid == 154 ) {
//	cout << thisone.x << ";" << thisone.y << ";" << thisone.z << "\n";
//}


				voro_p3d thisone = exterior_v_aligned[i];
				contour_pts.push_back( point_xy(dot(thisone, uu), dot(thisone, vv)) );
			}

			//dont forget to close the contour by adding the first point

//			voro_p3d thisone = exterior_v.at(exterior_n.at(ft->v_e-1));
//if ( iid == 154 ) {
//	cout << thisone.x << ";" << thisone.y << ";" << thisone.z << "\n";
//}

			voro_p3d thisone = exterior_v_aligned[ft->v_e-1];
			contour_pts.push_back( point_xy(dot(thisone,uu), dot(thisone,vv)) );

//if ( iid == 154 ) {
//	for( auto dt = contour_pts.begin(); dt != contour_pts.end(); ++dt )
//		cout << dt->x() << "\t\t" << dt->y() << "\n";
//}

			//project the hit point on the 3d plane also in 2d the same way
			point_xy hitproj = point_xy( dot(hit,uu), dot(hit,vv) );

//if ( iid == 154 ) {
//	cout << "Hitproj" << "\n";
//	cout << hitproj.x() << "\t\t" << hitproj.y() << "\n";
//}

			//check if hitproj is in polygon
			boost::geometry::model::polygon<point_xy> contour_poly;
			boost::geometry::assign_points( contour_poly, contour_pts );

			//fran fran_;
			//bool inside = boost::geometry::within( hitproj, contour_poly, fran_);
			bool inside = boost::geometry::covered_by( hitproj, contour_poly );


//if ( iid == 154 ) {
//	if ( inside == true)
//		cout << "in" << "\n";
//	else
//		cout << "out" << "\n";
//}


			//https://www.boost.org/doc/libs/1_62_0/libs/geometry/doc/html/geometry/reference/algorithms/within/within_3_with_strategy.html
			//https://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/covered_by/covered_by_3_with_strategy.html
			if ( inside == false ) { //most likely case, not inside test the next facet candidate
				//continue;
				//MK::there are geometrical shadows which are not covered by this solution and hence we have to evaluate the distance to all edges and vertices as well
				vector<voro_p3d> edges;
				voro_p3d a, b;
				for( size_t i = ft->v_s; i < ft->v_e-1; i++ ) {
					a = exterior_v_aligned[i];
					b = exterior_v_aligned[i+1];
					edges.push_back( voro_p3d( b.x-a.x, b.y-a.y, b.z-a.z ) );
				}
				a = exterior_v_aligned[ft->v_e-1];
				b = exterior_v_aligned[ft->v_s];
				edges.push_back( voro_p3d( b.x-a.x, b.y-a.y, b.z-a.z ) );

				for( auto edg = edges.begin(); edg != edges.end(); ++edg ) {
					int largest = 0; voro_real ba_max = edg->x;
					if ( edg->y >= ba_max ) {
						largest = 1; ba_max = edg->y;
					}
					if ( edg->z >= ba_max ) {
						largest = 2; ba_max = edg->z;
					}

					voro_p3d ortho = voro_p3d(1.0, 1.0, 1.0);
					if ( largest == 0 )
						ortho.x = (-edg->y - edg->z ) / edg->x;
					else if ( largest == 1 )
						ortho.y = (-edg->x - edg->z ) / edg->y;
					else
						ortho.z = (-edg->x - edg->y ) / edg->z;

					ortho.normalize();

					voro_real M[3][3] = { {edg->x, -ortho.x, pt->x - a.x},
					                      {edg->y, -ortho.y, pt->y - a.y},
										  {edg->z, -ortho.z, pt->z - a.z} };

					to_reduced_row_echelon_form(M);

					if ( M[2][0] < EPSILON && M[2][1] < EPSILON && M[2][2] < EPSILON ) {
						//there is a solution get beta M[1][2]
						//orthogonal projection from p on the polygon contour segment shoots past the segment
						voro_real beta = M[1][2];
						//voro_p3d pproj = voro_p3d( p.x + beta*ortho.x, p.y + beta*ortho.y, p.z + beta*ortho.z );
						//voro_real d = SQR(pproj.x-p.x)+SQR(pproj.y-p.y)+SQR(pproj.z-p.z);
						voro_real dist = SQR(beta*ortho.x)+SQR(beta*ortho.y)+SQR(beta*ortho.z);
						if ( dist > shortest_sqr_normal_distance )
							continue;
						else
							shortest_sqr_normal_distance = dist;
					}
				} //test next contour edge

				//finally test all vertices
				for( size_t i = ft->v_s; i < ft->v_e; i++ ) {
					voro_p3d vertex = exterior_v_aligned[i];
					voro_real dist = SQR(vertex.x-pt->x)+SQR(vertex.y-pt->y)+SQR(vertex.z-pt->z);
					if ( dist > shortest_sqr_normal_distance )
						continue;
					else
						shortest_sqr_normal_distance = dist;
				}
			}
			else { //inside compute distance between hitpoint hit and original material point
				voro_real dist = SQR(hit.x - pt->x)+SQR(hit.y - pt->y)+SQR(hit.z - pt->z);
				if ( dist > shortest_sqr_normal_distance ) { //likely most will not be the closest
					continue;
				} //test next facet
				else {
					shortest_sqr_normal_distance = dist;
//cout << shortest_sqr_normal_distance << "\n";
				}
			}

//			if ( warmedup == false ) {
//				voro_real R = sqrt(shortest_sqr_normal_distance) + EPSILON;
//				hedgesAABB e_aabb( trpl(pt->x-R, pt->y-R, pt->z-R), trpl(pt->x+R, pt->y+R, pt->z+R) );
//				facet_candidates.clear();
//
//				facet_candidates = bvh->query( e_aabb );
//				//first shot facet_candidates.)
//			}

		} //probe the next facet

		//##MK::take square root!
		if ( shortest_sqr_normal_distance > 50.0 ) {
			iotopo.push_back( 1 );
			iotopo.push_back( 1 );
			iotopo.push_back( faultypid );
			iogeom.push_back( pt->x );
			iogeom.push_back( pt->y );
			iogeom.push_back( pt->z );
			ioattr.push_back( shortest_sqr_normal_distance );
			faultypid++;
cout << iid << "\n";
		}

		interior_d.push_back( shortest_sqr_normal_distance );

cout << "Reporting shortest ___" << interior_d.back() << "___\n";
	}


	//##MK::debug faulty cases
	ofstream xdmfout;
	string xmlfn = "FaultyCases.xdmf";
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << iotopo.size() / 3 << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iotopo.size() << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iotopo.begin(); it != iotopo.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iogeom.size() / 3 << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iogeom.begin(); it != iogeom.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Cell\" Name=\"SQRDist\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ioattr.size() << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ioattr.begin(); it != ioattr.end(); ++it )
				xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";


		xdmfout.flush();
		xdmfout.close();
	}
	else {
		cerr << "Unable to open " << xmlfn << "\n";
	}


	double mytoc = omp_get_wtime();
	cout << "Computing the distances took " << (mytoc-mytic) << " seconds sequentially" << "\n";
}
*/


//#define REPRUNING_FRACTION 0.90
void microstructural_object::identify_normaldistances3( const voro_real threshold )
{
	double mytic = omp_get_wtime();
	size_t globalstats = 0;
/*
	vector<unsigned int> iotopo;
	vector<float> ioattr;
	vector<float> iogeom;
	unsigned int faultypid = 0;
*/

	size_t iid = 0;
	for( auto pt = interior_v.begin(); pt != interior_v.end(); ++pt, iid++ ) {
		unsigned int stats = 0;
		voro_real shortest_sqr_normal_distance = VOROFMX;
		bool reprune = false;
		//##MK::bounding box pre-checks to prune search space

		vector<unsigned int> facet_candidates;
		unsigned int f = 0;
		unsigned int nf = static_cast<unsigned int>(exterior_f.size());
		facet_candidates.reserve( nf );
		for( unsigned int f = 0; f < nf; ++f )
			facet_candidates.push_back( f );

		while( f < nf ) {
		//for(   ; f < nf;   ) {

			voro_facet & ft = exterior_f[facet_candidates[f]];
			stats++;

			//arbitrary point on the facet plane and outer unit normal
			voro_p3d r0 = ft.centrum;
			voro_p3d nou = ft.ounormal;
			//difference vector r0-pt
			voro_p3d v3 = voro_p3d( r0.x - pt->x, r0.y - pt->y, r0.z - pt->z );
			//normal projection of pt on the facet plane
			voro_real d = dot( nou, v3 );
			voro_p3d hit = voro_p3d( 	pt->x + nou.x*d,
										pt->y + nou.y*d,
										pt->z + nou.z*d );

			//project hit point + contour polygon of the facet vertices from 3d Cartesian to 2d Cartesian and point in polygon inclusion
			//using Boost library

			//Voro++ uses right-hand rule, ie we circulate anticlockwise but Boost by default circulates clockwise!
			//so we have to reverse order the points ##MK

			voro_p3d uu = ft.inplane_u;
			voro_p3d vv = ft.inplane_v;

			//project the hit point on the 3d plane also in 2d the same way
			point_xy hitproj = point_xy( dot(hit,uu), dot(hit,vv) );

			vector<point_xy> contour_pts;
			for( voro_uint i = ft.v_s_polyproj; i < ft.v_e_polyproj; ++i ) {
				contour_pts.push_back( clockwise_projv[i] );
			}

			/*
			for( size_t i = ft.v_e-1; i > ft.v_s; i-- ) {
				voro_p3d thisone = exterior_v_aligned[i];
				contour_pts.push_back( point_xy(dot(thisone, uu), dot(thisone, vv)) );
			}
			//dont forget to close the contour by adding the first point
			voro_p3d thisone = exterior_v_aligned[ft.v_e-1];
			contour_pts.push_back( point_xy(dot(thisone,uu), dot(thisone,vv)) );
			*/

			//check if hitproj is in polygon
			boost::geometry::model::polygon<point_xy> contour_poly;
			boost::geometry::assign_points( contour_poly, contour_pts );

			//fran fran_;
			//bool inside = boost::geometry::within( hitproj, contour_poly, fran_);
			bool inside = boost::geometry::covered_by( hitproj, contour_poly );

			//https://www.boost.org/doc/libs/1_62_0/libs/geometry/doc/html/geometry/reference/algorithms/within/within_3_with_strategy.html
			//https://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/covered_by/covered_by_3_with_strategy.html
			if ( inside == false ) { //most likely case, not inside test the next facet candidate
				//continue;
				//MK::there are geometrical shadows which are not covered by this solution and hence we have to evaluate the distance to all edges and vertices as well
				vector<voro_p3d> edges;
				voro_p3d a, b;
				for( size_t i = ft.v_s; i < ft.v_e-1; i++ ) {
					a = exterior_v_aligned[i];
					b = exterior_v_aligned[i+1];
					edges.push_back( voro_p3d( b.x-a.x, b.y-a.y, b.z-a.z ) );
				}
				a = exterior_v_aligned[ft.v_e-1];
				b = exterior_v_aligned[ft.v_s];
				edges.push_back( voro_p3d( b.x-a.x, b.y-a.y, b.z-a.z ) );

				for( auto edg = edges.begin(); edg != edges.end(); ++edg ) {
					int largest = 0; voro_real ba_max = edg->x;
					if ( edg->y >= ba_max ) {
						largest = 1; ba_max = edg->y;
					}
					if ( edg->z >= ba_max ) {
						largest = 2; ba_max = edg->z;
					}

					voro_p3d ortho = voro_p3d(1.0, 1.0, 1.0);
					if ( largest == 0 )
						ortho.x = (-edg->y - edg->z ) / edg->x;
					else if ( largest == 1 )
						ortho.y = (-edg->x - edg->z ) / edg->y;
					else
						ortho.z = (-edg->x - edg->y ) / edg->z;

					ortho.normalize();

					voro_real M[3][3] = { {edg->x, -ortho.x, pt->x - a.x},
					                      {edg->y, -ortho.y, pt->y - a.y},
										  {edg->z, -ortho.z, pt->z - a.z} };

					to_reduced_row_echelon_form(M);

					if ( M[2][0] < EPSILON && M[2][1] < EPSILON && M[2][2] < EPSILON ) {
						//there is a solution get beta M[1][2]
						//orthogonal projection from p on the polygon contour segment shoots past the segment
						voro_real beta = M[1][2];
						//voro_p3d pproj = voro_p3d( p.x + beta*ortho.x, p.y + beta*ortho.y, p.z + beta*ortho.z );
						//voro_real d = SQR(pproj.x-p.x)+SQR(pproj.y-p.y)+SQR(pproj.z-p.z);
						voro_real dist = SQR(beta*ortho.x)+SQR(beta*ortho.y)+SQR(beta*ortho.z);
						if ( dist > shortest_sqr_normal_distance )
							continue;
						else {
							if ( reprune == false && ((dist / shortest_sqr_normal_distance) <= threshold ) ) { //REPRUNING_FRACTION) ) {
								reprune = true;
							}
							shortest_sqr_normal_distance = dist;
						}
					}
				} //test next contour edge

				//finally test all vertices
				for( size_t i = ft.v_s; i < ft.v_e; ++i ) {
					voro_p3d vertex = exterior_v_aligned[i];
					voro_real dist = SQR(vertex.x-pt->x)+SQR(vertex.y-pt->y)+SQR(vertex.z-pt->z);
					if ( dist > shortest_sqr_normal_distance ) {
						continue;
					}
					else {
						if ( reprune == false && ((dist / shortest_sqr_normal_distance) <= threshold ) ) { //REPRUNING_FRACTION) ) {
							reprune = true;
						}
						shortest_sqr_normal_distance = dist;
					}
				}
			}
			else { //inside compute distance between hitpoint hit and original material point
				voro_real dist = SQR(hit.x - pt->x)+SQR(hit.y - pt->y)+SQR(hit.z - pt->z);
				if ( dist > shortest_sqr_normal_distance ) { //likely most will not be the closest
					//f++;
					//continue;
				} //test next facet
				else {
					if ( reprune == false && ((dist / shortest_sqr_normal_distance) <= threshold ) ) { //REPRUNING_FRACTION) ) {
						reprune = true;
					}
					shortest_sqr_normal_distance = dist;
//cout << shortest_sqr_normal_distance << "\n";
				}
			}

			if ( reprune == false ) {
				f++;
			}
			else {
//cout << "Repruning " << shortest_sqr_normal_distance << "\n";
				facet_candidates.clear();
				voro_real R = sqrt(shortest_sqr_normal_distance) + EPSILON;
				hedgesAABB e_aabb( trpl(pt->x-R, pt->y-R, pt->z-R), trpl(pt->x+R, pt->y+R, pt->z+R) );
				facet_candidates = bvh->query( e_aabb );
				f = 0;
				nf = static_cast<unsigned int>(facet_candidates.size());
				reprune = false;
			}

		} //probe the next facet
/*
		//##MK::take square root!
		if ( shortest_sqr_normal_distance > 50.0 ) {
			iotopo.push_back( 1 );
			iotopo.push_back( 1 );
			iotopo.push_back( faultypid );
			iogeom.push_back( pt->x );
			iogeom.push_back( pt->y );
			iogeom.push_back( pt->z );
			ioattr.push_back( shortest_sqr_normal_distance );
			faultypid++;
cout << iid << "\n";
		}
*/
		interior_d.push_back( sqrt(shortest_sqr_normal_distance) );
		interior_stats.push_back( stats );

		globalstats += stats;

//cout << "Reporting shortest ___" << interior_d.back() << "___ " << stats << " facets tested" << "\n";
	}

/*
	//##MK::debug faulty cases
	ofstream xdmfout;
	string xmlfn = "FaultyCases.xdmf";
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << iotopo.size() / 3 << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iotopo.size() << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iotopo.begin(); it != iotopo.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iogeom.size() / 3 << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iogeom.begin(); it != iogeom.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Cell\" Name=\"SQRDist\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ioattr.size() << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ioattr.begin(); it != ioattr.end(); ++it )
				xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();
	}
	else {
		cerr << "Unable to open " << xmlfn << "\n";
	}
*/
	double mytoc = omp_get_wtime();
	/*#pragma omp critical
	{
		cout << "Total number of facets tested for distancing " << globalstats << "\n";
		cout << "Computing the distances took " << (mytoc-mytic) << " seconds sequentially" << "\n";
	}*/
}


bool microstructural_object::visualize_object_hull( const unsigned int simid, const unsigned int incr, const unsigned int gid )
{
	double mytic = omp_get_wtime();

	vector<unsigned int> iotopo;
	vector<unsigned int> ioattr;
	for( auto ft = exterior_f.begin(); ft != exterior_f.end(); ++ft ) {
		iotopo.push_back( 3 ); //XDMF topo key polygon

		size_t nvertices = ft->v_e - ft->v_s;
		iotopo.push_back( nvertices ); //XDMF polygon key
		for( size_t vv = 0; vv < nvertices; ++vv ) {
			iotopo.push_back( static_cast<unsigned int>(exterior_n[ft->v_s + vv]) );
		}

		ioattr.push_back( ft->owner );
	}

	vector<float> iogeom;
	for( auto vt = exterior_v.begin(); vt != exterior_v.end(); ++vt ) {
		iogeom.push_back( static_cast<float>(vt->x) );
		iogeom.push_back( static_cast<float>(vt->y) );
		iogeom.push_back( static_cast<float>(vt->z) );
	}

	//compute outer unit normals
	vector<float> ionormals;
	for( auto ft = exterior_f.begin(); ft != exterior_f.end(); ++ft ) {
		/*
		size_t prev = ft->v_e - 1;
		size_t here = ft->v_s;
		size_t next = ft->v_s + 1;
		voro_p3d prevp3d = exterior_v.at(exterior_n.at(prev));
		voro_p3d herep3d = exterior_v.at(exterior_n.at(here));
		voro_p3d nextp3d = exterior_v.at(exterior_n.at(next));
		voro_p3d next_here = voro_p3d( nextp3d.x - herep3d.x, nextp3d.y - herep3d.y, nextp3d.z - herep3d.z );
		voro_p3d prev_here = voro_p3d( prevp3d.x - herep3d.x, prevp3d.y - herep3d.y, prevp3d.z - herep3d.z );
		voro_p3d normal = cross( next_here, prev_here );
		normal.normalize();

		voro_p3d cellcenter = interior_v.at(ft->owner);
		voro_p3d cellcenter_here = voro_p3d( cellcenter.x - herep3d.x, cellcenter.y - herep3d.y, cellcenter.z - herep3d.z );
		cellcenter_here.normalize();

		//projection on normal if positive on the same side, given that cellcenter is inside this will be the inner normal
		voro_real projection = dot( cellcenter_here, normal );

		//##MK::umlaufsinn beachten ggf. nicht consistent fuer Voro++
		//##MK::Umlaufsinn for Voro++ naturally is anticlockwise
		//##MK::some Boost functions though need clockwise, so for this we have the contour reversal
		if ( projection > static_cast<voro_real>(0.0) ) { //positive on the same side than is the cellcenter so for a Voronoi cell this is
			//an inner unit normal vector to the facet, so flip the normal
			ft->ounormal = voro_p3d( -normal.x, -normal.y, -normal.z );
		}
		else { //normal is already an outer unit normal
			ft->ounormal = voro_p3d( normal.x, normal.y, normal.z );
		}
		 */

		//##MK::debugging inner unit normal
		ionormals.push_back( ft->ounormal.x );
		ionormals.push_back( ft->ounormal.y );
		ionormals.push_back( ft->ounormal.z );

		/*
		//compute in plane normal vector u
		//take inplane_u projection vector using the largest projected world coordinate vector
		//x - dot(x,oun)*oun
		int largest = 0;
		voro_real val = fabs(ft->ounormal.x);
		if ( fabs(ft->ounormal.y) >= val ) {
			largest = 1;
			val = fabs(ft->ounormal.y);
		}
		if ( fabs(ft->ounormal.z) >= val ) {
			largest = 2;
			val = fabs(ft->ounormal.z);
		}

		voro_p3d wx = voro_p3d( 1.0, 0.0, 0.0 );
		if ( largest == 1 )
			wx = voro_p3d( 0.0, 1.0, 0.0 );
		if ( largest == 2 )
			wx = voro_p3d( 0.0, 0.0, 1.0 );
		voro_real d = dot( wx, ft->ounormal );

		ft->inplane_u = voro_p3d( 	wx.x - d*ft->ounormal.x,
									wx.y - d*ft->ounormal.y,
									wx.z - d*ft->ounormal.z   );
		ft->inplane_u.normalize();

		//second in plane normal vector v
		ft->inplane_v = cross( ft->inplane_u, ft->ounormal );
		ft->inplane_v.normalize();

		//hyperphobic
		if ( fabs(dot(ft->inplane_u, ft->ounormal)) < EPSILON && fabs(dot(ft->inplane_v, ft->ounormal)) < EPSILON )
			continue;
		else
			cerr << "Facet inplane u or inplane_v not perpendicular to ounormal" << "\n";
		*/
	}

	ofstream xdmfout;
	string xmlfn = "DAMASKPDT.SimID." + to_string(simid) + ".Incr." + to_string(incr) + ".GrainID." + to_string(gid) + ".VoroFacetHull.xdmf";
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << exterior_f.size() << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iotopo.size() << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iotopo.begin(); it != iotopo.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << exterior_v.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iogeom.begin(); it != iogeom.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Cell\" Name=\"Clpid\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << exterior_f.size() << " 1\" DataType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ioattr.begin(); it != ioattr.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "      <Attribute AttributeType=\"Vector\" Center=\"Cell\" Name=\"Ounormal\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << exterior_f.size() << " 3\" DataType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ionormals.begin(); it != ionormals.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();
	}
	else {
		cerr << "Unable to open " << xmlfn << "\n";
		return false;
	}

	double mytoc = omp_get_wtime();
	cout << "Visualizing the complex to XDMF file took " << (mytoc-mytic) << " seconds sequentially" << "\n";
	return true;
}


bool microstructural_object::visualize_distances( const unsigned int simid, const unsigned int incr, const unsigned int gid )
{
	double mytic = omp_get_wtime();

	vector<unsigned int> iotopo = vector<unsigned int>( 3*interior_v.size(), 1 ); //implicitly set XDMF geometry type keyword and object keyword single point
	vector<float> iogeom;
	vector<float> ioattr1;
	vector<unsigned int> ioattr2;
	iogeom.reserve( 3*interior_v.size() );
	ioattr1.reserve( 1*interior_v.size() );
	ioattr2.reserve( 1*interior_v.size() );

	for( size_t i = 0; i < interior_v.size(); ++i ) {
		iotopo[3*i+2] = i;

		iogeom.push_back( interior_v[i].x );
		iogeom.push_back( interior_v[i].y );
		iogeom.push_back( interior_v[i].z );

		ioattr1.push_back( interior_d[i] );
		ioattr2.push_back( interior_stats[i] );
	}

	ofstream xdmfout;
	string xmlfn = "DAMASKPDT.SimID." + to_string(simid) + ".Incr." + to_string(incr) + ".GrainID." + to_string(gid) + ".Distancing.xdmf";
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << interior_v.size() << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iotopo.size() << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iotopo.begin(); it != iotopo.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << iogeom.size() / 3 << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = iogeom.begin(); it != iogeom.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";

		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Dist\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ioattr1.size() << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ioattr1.begin(); it != ioattr1.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"FacetChecks\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ioattr2.size() << " 1\" DataType=\"UInt\" Precision=\"4\" Format=\"XML\">" << "\n";
		for( auto it = ioattr2.begin(); it != ioattr2.end(); ++it )
			xdmfout << *it << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";

		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();
	}
	else {
		cerr << "Unable to open " << xmlfn << "\n";
		return false;
	}

	double mytoc = omp_get_wtime();
	cout << "Visualizing the complex to XDMF file took " << (mytoc-mytic) << " seconds sequentially" << "\n";
	return true;
}



voroComposer::voroComposer()
{
}


voroComposer::voroComposer( const voro_uint nobjects )
{
	for( voro_uint i = 0; i < nobjects; ++i ) {
		microstructural_object* msobj = NULL;
		try {
			msobj = new microstructural_object();
		}
		catch (bad_alloc &croak) {
			cerr << "Allocation of the " << i << "-th object failed" << "\n";
		}
		obj.push_back( msobj );
	}
}


voroComposer::~voroComposer()
{
	for( auto it = obj.begin(); it != obj.end(); ++it )
		delete *it;
}
