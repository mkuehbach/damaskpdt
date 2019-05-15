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


#ifndef __PDT_SPECOUTREADER_H__
#define __PDT_SPECOUTREADER_H__

//#include "PDT_GrainObject.h"
//#include "PDT_Crystallography.h"
//#include "PDT_IntelMKL.h"
//#include "PDT_VTKIO.h"
//#include "PDT_OriMath.h"
#include "LOUVAIN_Core.h"

//##MK::uninclude uper include lower
//#include "PDT_GrainHdl.h"

class mesh;
class homog;
class cryst;
class consti;
class memRegion;
class rveAverageResults;
class specOutHeader;
class grainpool;
class grGeomHdl;
class specOutIncr;



class taskspecifics
{
public:
	taskspecifics(){
		whichstraintensor = RIGHTCAUCHYGREEN; //the default
		whichstrainmodel = LNSTRAIN;
	};
	~taskspecifics(){};
	//bookkeeps details what should be done specifically for a task

	string get_strainopts(){
		string which = "";
		switch (whichstraintensor) {
			case RIGHTCAUCHYGREEN:
				which += "Right";
				break;
			case LEFTCAUCHYGREEN:
				which += "Left";
				break;
			default:
				break;
		}
		switch (whichstrainmodel) {
			case LNSTRAIN:
				which += "Ln";
				break;
			case BIOTSTRAIN:
				which += "Biot";
				break;
			case GREENSTRAIN:
				which += "Green";
				break;
			default:
				break;
		}
		return which;
	}

	STRAINTENSOR_TYPE whichstraintensor;
	STRAINVALUE_KIND whichstrainmodel;
};


class dlayoutnode
{
	//functionality to bundle pairs, triplets, and higher order groups of spectralOut data layout data
	//MK::one could utilize a std::map container for this which is fast for random access but keeps no element insertion history
	//on the contrary, std::vector or std::list have these implicit insertion history but allow to store T types only, therefore
	//the choice to make a small dummy class bundling the node pieces of information
public:
	dlayoutnode() {
		key = "";
		valstr = "";
		valui = 0;
	};
	dlayoutnode(const string _k, const string _v) : key(_k), valstr(_v) {
		try {
			valui = static_cast<size_t>(stoul(valstr));
		}
		catch (invalid_argument &croak) {
			valui = 0;
		}
	};
	~dlayoutnode() {};

	string key;
	string valstr;
	size_t valui;
};


class specmeta
{
	//construct maintaining the metadata to a spectralOut file such as its internal data layout and structure
	friend class specOutIncr;

public:
	specmeta();
	~specmeta();

	inline unsigned int constitutive_getvalui(const string keyword) {
		for(size_t kw = 0; kw < const_meta.size(); ++kw) {
			if (const_meta.at(kw).key.compare(keyword) != 0)
				continue; //most likely not the keyword
			else
				return const_meta.at(kw).valui; //0 if formally existent but originated from an invalid_argumen
		}
		return 0;
	}

	//storage order in these vectors corresponds with spectralOut file layout
	vector<dlayoutnode> homog_meta;
	vector<dlayoutnode> cryst_meta;
	vector<dlayoutnode> const_meta;

	//##MK::introducing a dictionary-based bookkeeping of metadata
};


class mesh
{
	friend class specOutIncr;
	//functionality instance taking care for the location of the mesh elements' (elem) integration points (ips)
	//their periodic replication, the displacement calculation as well as the mapping of mesh element integration point ids eipid
public:
	mesh();
	~mesh();

	memRegion* owner;

//private:
	bool init();
	void eipid2xyz0();

	vector<p3d> xyz0;				//coordinates of mesh integration point in initial configuration!
	vector<d3d> dxyz_avg;			//average displacement of integration point p between deformed to initial configuration
	vector<d3d> dxyz_flu;			//displacements which allow to map x0 to x in the deformed configuration via x0 + dx --> x

	bv3x3 perOffset;				//the displacement vector for reconstructing the periodic images of the ips in the deformed configuration
									//by applying the RVE-volume averaged deformation gradient F to the initial RVE domain dimensions LX, LY, LZ
};


class homog
{
	//functionality instance taking care of the homogenization method data
	friend class specOutIncr;
public:
	homog();
	~homog();

	memRegion* owner;

private:
	bool init();

	vector<unsigned int> gid;			//id of the grain assigned to this material point
};


class cryst
{
	//functionality instance taking care of the crystallite specific data
	friend class specOutIncr;
public:
	cryst();
	~cryst();

	memRegion* owner;

//##MK::better data encapsulation
//private:
	bool init();
	real_xyz compute_total_V();
	t3x3 compute_average_Fp();
	t3x3 compute_average_F();
	t3x3 compute_average_P();
	real_xyz compute_summate_equalweighting( vector<real_xyz> const & in );
	t3x3 compute_average_equalweighting( vector<t3x3> const & in );

	//possible content of a spectralOut or HDF5 file from DAMASK
	vector<unsigned int> PhaseID;
	vector<unsigned int> TextureID;
	vector<real_xyz> V;
	vector<squat> q;
	vector<t3x3> F;
	vector<t3x3> Fe;
	vector<t3x3> Fp;
	vector<t3x3> P;

	//crystallite-focused results from the postProcessing
	vector<t3x3> straintensor1; //##MK::changed for Matthew EPS_RIGHTCAUCHYGREEN_LN_FT
	vector<t3x3> straintensor2; //##MK::changed for Matthew EPS_RIGHTCAUCHYGREEN_LN_FP
	vector<t3x3> stresstensor;
	vector<vMises> scalars;

	//derived results
	vector<unsigned int> GrainID;
};


class consti
{
	friend class specOutIncr;
	//functionality instance taking care of constitutive model specific data
public:
	consti();
	~consti();

	memRegion* owner;

//##MK::better data encapsulation
//private:
	bool init();

	vector<real_rho> rho_e;
	vector<real_rho> rho_d;

	size_t rho_e_mult;						//multiplicity, for instance storing rho_e with 12 values per el ip set rho_e_mult = 12
	size_t rho_d_mult;
};


//MK::basic idea: three layer work partitioning
//all MPI processes read spectralOut header, but get assigned different strain increments,
//possible, because so far (10/2018) DAMASK spectral simulations with more than even
//256^3 elements in particular for constitutive models with explicit dislocation density tracking are practically not executable
//apart from this, a hypothetical 1600^3 DAMASK simulation with 225 data elements per mesh element would already,
//with the present I/O concept generate a 6.7TB snapshot file --- per increment!
//therefore we agree that MPI processes work on different strain increments,
//the MPI processes individually spawn (OpenMP) threads which work on spatial disjoint memory contiguous data
//of the process' increment
//the threads individually may utilize SIMD where practical at core runtime-decisive parts functionalities
//this design suggest to map a spectralOut file to a node of multiple multicore processors
//mapping each processor an MPI process, the underlying cores of the processors a thread, and the SIMD lanes of
//each core the SIMD commands, this allows distinct resource utilization closer to theoretical peak performance

class memRegion
{
	friend class specOutIncr;

	//functionality splitting the heavy data along the z directory into xy slabs of mesh elements
	//allows OpenMP thread-parallelized heavy data allocation in thread local memory via clear first touch and
	//lower accessing of metadata through specOutIncr handler
	//thus an MPI process operates via its proxy specOutIncr who reads metadata from specOutHeader in all threads
	//but processes threadlocal data in independent data containers, the memRegion class objects, boosted individually by SIMD
public:
	memRegion();
	~memRegion();

	void db_distr_homogenization2( const double* raw, const size_t eipfirst, const size_t eipn );
	void db_distr_crystallite2( const double* raw, const size_t eipfirst, const size_t eipn );
	void db_distr_constitutive2( const double* raw, const size_t eipfirst, const size_t eipn );


	specOutIncr* owner;

	size_t eipid_start;		//which elements we process first index
	size_t eipid_end;			//first after last index, i.e. [eip_start,eip_end[
	size_t eipid_n;

//private:
	mesh grid;
	homog homogenization;
	cryst crystallite;
	consti constitutive;
};


class rveAverageResults
{
	friend class specOutIncr;
public:
	//collective instance carrying RVE average values
	rveAverageResults();
	rveAverageResults(
			const real_xyz Vtot,
			t3x3 const & Fp,
			t3x3 const & F, t3x3 const & P,
			t3x3 const & eps, t3x3 const & cauchy,
			vMises const & vm, const unsigned int lcid,
			const unsigned int lincrid, const unsigned int gincrid);
	/*rveAverageResults(
				const real_xyz Vtot,
				MPI_Tensor3x3_Double const & Fp,
				MPI_Tensor3x3_Double const & F,
				MPI_Tensor3x3_Double const & P,
				MPI_Tensor3x3_Double const & eps,
				MPI_Tensor3x3_Double const & cauchy,
				vMises const & vm, const unsigned int lcid,
				const unsigned int lincrid, const unsigned int gincrid);*/

	~rveAverageResults();

	//RVE-averaged quantities
	real_xyz Vtotal;			//total volume of all grid points
	t3x3 Fpavgrve;				//plastic deformation gradient
	t3x3 Favgrve;				//deformation gradient
	t3x3 Pavgrve;				//piola stress
	t3x3 Strainavgrve;			//strain tensor
	t3x3 Cauchyavgrve;			//cauchy stress
	vMises Equivavgrve;			//equivalency stress and strain according to von Mises

	unsigned int loadcaseID;	//in which loadcase
	unsigned int localincrID;	//which increment of this loadcase
	unsigned int globalincrID;	//which increment given that there are potentially multiple loadcases

	//MK::utilize these quantities to get point on flowcurve, i.e. point in equivalent true strain-stress space
};


class rveAverageResults2
{
	friend class specOutIncr;
public:
	//collective instance carrying RVE average values
	rveAverageResults2();
	~rveAverageResults2();

	//dictionary of RVE-averaged quantities
	map<unsigned int, real_m33> RVEAvgScalar;
	map<unsigned int, t3x3> RVEAvgTensorial;
	whenTaken info;
};

class specOutHeader
{
	friend class specOutIncr;
public:
	specOutHeader();
	~specOutHeader();

	specOutIncr* owner;

	unsigned int get_nincr(const unsigned int lc){ return loadcasesmeta.at(lc).nincr; }
	unsigned int get_lincr() { return lincr; }
//private:
 	string load;						//loadFile
	string wdir;						//working dir
	string geom;						//geometryFile

	unsigned int N;						//RVE mesh element count along one dimension
	unsigned int Nip;					//number of integration points per mesh element
	unsigned int Ncp;					//N^3 for DAMASK_spectral

	unsigned int matpres; 				//material point results
	unsigned int loadcases;
	vector<struct lcasemeta> loadcasesmeta;
	//unsigned int freqs;				//one for each loadcase
	//double times;						//one for each loadcase
	//unsigned int logscales;			//one for each loadcase
	//unsigned int nincr;				//one for each loadcase

	unsigned int sincr;					//starting increment
	unsigned int lincr;					//id of last increment sincr+[loadcase 0 all increments]+[loadcase 1 all increments]+....+loadcase[N-1 all increments] //##MK::mind different loadcases potentially have different frequency stepping!

	unsigned int NX;					//quick access
	unsigned int NY;
	unsigned int NZ;
	unsigned int NXY;
	unsigned int NXYZ;

	v3x1 L;								//RVE physical extent along edges parallel coordinate axes, initial undeformed configuration!

	size_t DataElementsPerIncrement;	//total number of materialpoint_result elements (all ips, all el)
	size_t FirstByteAfterHeader;		//byte offset specifying where data parts in spectralOut start
};



class bvh_xyzm2
{
	friend class specOutIncr;

public:
	bvh_xyzm2();
	~bvh_xyzm2();

	void aabb3d_about_deformed();
	void aabb3d_add_guardzone();
	void spatialpartitioning_init();
	inline unsigned int binning_x( const real_xyz x );
	inline unsigned int binning_y( const real_xyz y );
	inline unsigned int binning_z( const real_xyz z );
	inline unsigned int binning_xyz( const p3d candidate );
	void spatialpartitioning_deformedrve_perform();
	void build_bvh_xyzm2();
	void find_higherorder_neighbors( const p3d p, vector<dist>& results, real_xyz r );
	void destroy_bvh_xyzm2();
	inline aabb3d get_owindeformed();

	specOutIncr* owner;

private:
	aabb3d owin_deformed;
	aabb3d owin_final;
	sqb spatialdistr;

	//MK::splitting into individual coordinate vectors for SIMD purposes
	vector<vector<p3d>*> pxyz;			//the outer vector wraps a thread-local storage to a vector of mesh element coordinates in the deformed configuration
	vector<vector<unsigned int>*> pid;	//the metadata of this gridpoint
										//MK::such carrying of metadata is practical because as we have periodic images of the elements to implement the kernel scanning,
										//these periodic images would requiring the storage of significant redundant physical element data, instead, we utilize that only
										//the original mesh elements have relevant metadata, thereby reducing the total amount of metadata, thus the memory footprint, thus
										//reducing the amount of memory operations to load metadata for processing points
	vector<vector<unsigned char>*> pimage;	//which periodic image is it?
	bool healthy;
};


class bvh_p3dm1
{
	friend class specOutIncr;

public:
	bvh_p3dm1();
	~bvh_p3dm1();

	void spatialpartitioning_init();
	inline unsigned int binning_x( const real_xyz x );
	inline unsigned int binning_y( const real_xyz y );
	inline unsigned int binning_z( const real_xyz z );
	inline unsigned int binning_xyz( const p3d candidate );

	void build_bvh_p3dm1();
	nbp3d find_nearest_neighbor( const p3d p, real_xyz r );
	void destroy_bvh_p3dm1();
	inline aabb3d get_owin();
	inline sqb get_spatialdistr();

	specOutIncr* owner;
private:
	aabb3d owin;
	sqb spatialdistr;

	//##MK::splitting into individual coordinate vectors for SIMD purposes, consider to split in point data and label data as during finding the nearest
	vector<vector<p3dm2>*> pp3_points;
	bool healthy;
};



class grainpool
{
	friend class specOutIncr;
public:
	grainpool();
	~grainpool();

	void build(vector<unsigned int> const & ip2gr);
	void disori_uppertriangle_matrix();
	void report();
	void destroy();

	specOutIncr* owner;

	vector<vector<unsigned int>*> ipsupport; //for each grain, which unique ips build it? necessary as well for orientation averaging of grains which as soon as deformation has been applied are composed of ips with different orientations requiring an averaging procedure
	vector<grain> polyxx;
};

//helper to report shape and base vector of RVE for plotting for scientific paper
class shape
{
public:
	shape();
	~shape();
	vector<hexahedron> sp;
	vector<string> nm;

	void add( const aabb3d & in, const string str );
	void add( const bv3x3 & in, const string str );
};


class cluster
{
public:
	cluster();
	~cluster();

	unsigned int lid;
	unsigned int cnt;
	aabb3d box;
	bool chosen;
};


class grGeomHdl
{
public:
	grGeomHdl();
	~grGeomHdl();

	bool debug_duplicate_check( const real_xyz x, const real_xyz y, const real_xyz z );
	void compute_all_replica();
	//void update_threadlocalbox( aabb3d & current_extremal );
	void build_sparse_bvh();
	void spatial_range_query_sbvh_noclear_nosort(const size_t thisip, vector<unsigned int> & results, real_xyz r );
	void dbscan_ips2grreplicates1( const real_xyz eps ); //##MK::seems there is a flaw here somewhere....
	void dbscan_ips2grreplicates2( const real_xyz eps );
	void dbscan_pick_representative();
	void build_localgrid();
	void voxelize_via_pvtessellation();
	void compute_sgndistfun_coarse();
	void compute_sgndistfun_fsm();

	void build_contour_from_vorotess();

	inline real_sdf m_inputDistance_getValueAt( const int row, const int col, const int dep );
	inline real_sdf m_outputDistance_getValueAt( const int row, const int col, const int dep );
	inline void m_outputDistance_setValueAt( const int row, const int col, const int dep, const real_sdf val );
	void debug_sdf_sphere_seed( const real_xyz R);
	void debug_sdf_sphere_guess();
	void debug_sdf_sphere_fsm();
	void debug_sdf_sphere_exact( const real_xyz R );
	void debug_sdf_sphere_report();

	specOutIncr* owner;

	unsigned int cgid;			//ID of grain which this Hdl object works on
	unsigned int nuip;			//how many unique ips were during the community detection identified for the grain?
	aabb3d localbox;			//an axis-aligned bounding box about the ips, local because only bounding all ips of that grain cgid

	vector<p3dm1> ips;			//all integration points (periodic/unique) the mark is the uniqueip eid (eipid) from which the point was
								//was generated by evaluating periodicity, supporting the grain within the simulated RVE
								//volume as well as its periodic replica in the general, i.e. deformed configuration
								//not ordered! as upon construction of the periodic ips the final size of the to be binned container may not be known

	map<unsigned int, sbvhrange> sbvh_locator;
	vector<p3d> sbvh_points;
	vector<unsigned int> sbvh_uipref;
	vector<unsigned char> sbvh_isrepresentative;		//true if ip is part of the representatively chosen dbscan cluster
	//##MK::split into two linear memory pieces because for spatial range queries additional integer mark not required, hence improving cache line reutilization
	//in fact the uip reference are needed later: when computing signed distance function related to refer to simulated state variable values at a unique ip


	vector<unsigned int> lbl;	//labels of ips in sbvh as ordered as in sbvh resulting from the DBScan algorithm
	vector<cluster> dbscanres;
	vector<p3d> thegrain;		//the representative one chosen

	aabb3d grainfence;			//axis-aligned bounding box about the ips sampling thegrain, i.e. the chosen representative of the periodic replicates
	vxlgrd localgrid;			//a local vxlgrd embedded in the global specOutIncr::thegrid, local implicit x+y*NX+z*NXY convention
	size_t nuvxl;				//how many voxel of localgrid with GrainIDField[i] == cgid

	vector<unsigned int> GrainIDField;	//##MK::DEBUG
	vector<unsigned int> UIPIDField;
	vector<unsigned char> IsRepresentativeField;		//##MK::build bitmap in future storing if ip in RVE27 that is closest to this voxel is included in the point cloud of the representative dbscan cluster chosen for this grain with id cgid

	vector<real_sdf> SgnDistField1;
	vector<real_sdf> SgnDistField2;

	microstructural_object contour;

	bool alloc_success;
	bool ips_success;
	bool sbvh_success;
	bool dbscan_success;
	bool dbpick_success;
	bool vxlinit_success;
	bool vxlfill_success;
	bool sdfinit_success;
	bool sdfsprd_success;

	bool healthy;
};




/*class analyzerAddStrainTensors
{
	friend class specOutIncr;

public:
	analyzerAddStrainTensors();
	~analyzerAddStrainTensors();

	specOutIncr* owner;
};*/



class specOutIncr
{
	//functionality instance per process to read DAMASK *.spectralOut binary files
	//each increment has a different mesh because strictly speaking we are interested
	//in the relative displacement of material points in the global sample coordinate system,
	//as we for instance care for a true distance to a boundary, a boundary which may become
	//tilted, displaced or otherwise changes location during deformation, or becomes generated in-situ

public:
	specOutIncr();
	~specOutIncr();

	void init_mpi_derivedtypes();
	void parse_taskspecific_opts();

	bool specout_read_header();
	bool specout_read_structure_homogenization();
	bool specout_read_structure_crystallite();
	bool specout_read_structure_constitutive();

	int map_increments2ranks();
	void map_meshelements2threads();
	void grid_initial_configuration();
	real_xyz rve_volume_total();
	t3x3 rve_volume_averaged_defpgradient();
	t3x3 rve_volume_averaged_defgradient();
	t3x3 rve_volume_averaged_piolastress();
	t3x3 rve_volume_averaged_strain();
	t3x3 rve_volume_averaged_truestrain( t3x3 const & Fav, t3x3 const & Pav );
	t3x3 rve_volume_averaged_cauchystress( t3x3 const & Fav, t3x3 const & Pav );
	real_xyz rve_summated_equalweighting( const unsigned int which );
	t3x3 rve_averaged_equalweighting( const unsigned int which );
	vMises rve_volume_averaged_scalars( t3x3 const & eps, t3x3 const & cau);
	void analyze_addRVEAverages( const unsigned int lc, const unsigned int li, const unsigned int glbincrid );
	void analyze_addRVEAverages2( const unsigned int lc, const unsigned int li, const unsigned int glbincrid );
	//void analyze_optionalRVEAverages( const unsigned int lc, const unsigned int li, const unsigned int glbincrid );

	void analyze_addStrainTensors_mp();
	void analyze_addCauchy();
	void analyze_addVonMises();

	void analyze_addDisplacements();
	void analyze_ignoreDisplacements();
	void analyze_ipgrid_displacements();

	//void specout_read_heavydata1();
	void specout_read_heavydata2();
	bool specout_check_heavydata2();
	void write_ipgrid_textureid();
	void bounded_volume_hierarchy();
	void analyze_svar_closestuip_disoriangle();
	void hierarchical_community_detection( vector<lvwtedge> const & edgs, vector<unsigned int> & uip2community );
	void analyze_identify_grains();
	void analyze_svar_grainbased_sdf();
	void analyze_svar_grainbased_cgeom();

	bool check_if_myworkpackage(unsigned int const i, unsigned int const tid, unsigned int const ntid);
	bool init_threadlocalmemory_sdf();
	void compute_perips_and_uips();
	void init_global_bvh_p3dm1();
	void build_grainlocal_sbvh();
	void discern_replica_via_dbscan();
	void pick_one_replica_from_dbscan();
	void write_dbscan_result();
	void init_global_voxelgrid_csys();
	void init_grainlocal_voxelgrids_csys();
	void extract_grainlocal_cuboidalregion_from_rve27();
	void write_grainlocal_vxlgrids();
	void fill_global_bvh_p3dm1();
	void voxelize_this_replica();
	void approximate_signed_distance_function();
	void spread_signed_distance_function();


	void analyze_boxup_grains();
	void analyze_build_grains();
	//void analyze_reconstruct_grains();

	inline void eid_write_gid( const size_t eid, const unsigned int gid );
	inline squat eid2quaternion( const size_t eid );
	inline p3d eid2p3d( const size_t eid );
	inline unsigned int eid2textureid( const size_t eid );
	inline unsigned int eid2gid( const size_t eid );
	dist closestBoundaryInKernelIfAny( const size_t central_eid, vector<dist> const & candidates ); //, vector<size_t>& suspicious );
	/*void write_histograms(vector<spatdist> const & buf);*/
	void write_searchefficiency(vector<unsigned int> const & ntested );
	/*
	void write_histograms2_ascii(vector<real_xyz> const & bufd, vector<slipsysdata_fcc> const & bufedge, vector<slipsysdata_fcc> const & bufdipo);
	void write_histograms2_binary( vector<MPI_DisloSpatDistr_Double> const & bufedge, vector<MPI_DisloSpatDistr_Double> const & bufdipo);
	*/
	void write_histograms2_binary_stress( vector<MPI_Tensor3x3SpatDistr_Double> const & bufstress );
	void write_histograms2_binary_ori( vector<MPI_GrainOriSpatDistr_Double> const & bufori );


	void compute_edge_weights( const size_t central_eid, vector<dist> const & candidates, vector<lvwtedge> & edgs );
	void write_edge_weights( vector<lvwtedge> const & edgs );
	/*
	unsigned int anyGrainAlreadyAssigned( const size_t central_eid, vector<dist> const & candidates,
			vector<unsigned int> const & p2g, vector<grain> & grpool );
	unsigned int anyGrainAlreadyExistent1( const size_t central_eid, vector<dist> const & candidates,
			vector<unsigned int> const & p2g, vector<grain> & grpool );
	*/
	//unsigned int anyGrainAlreadyExistent2( const size_t central_eid, vector<dist> const & candidates,
	//			vector<unsigned int> const & p2gtmp, vector<grain> const & oldgrpool, vector<grain> & newgrpool );

	void write_reconstructed_grains1( vector<grain> const & grpool );
	void write_identified_grains();
	void write_grainid_and_quaternions();
	void write_grainid_and_quaternions_mpiio();
	void report_rveshapes();

	void free_increment_heavydata();
	void free_increment_bvh_xyzm2();
	void free_increment_bvh_p3dm1();
	void free_increment_grains();
	void free_increment_sdf();


	//void report_flowcurve();
	void report_flowcurve2();
	void report_flowcurve3( const string fbase, const unsigned int which );


	void debug_signed_distance_function();

	specOutHeader head;				//pertains to multiple increments
	specmeta dlayout;

	vector<memRegion*> db;			//pertains to individual increments therefore requiring resetting after processing each increment
	bvh_xyzm2 pp3rve1withguard;
	bvh_p3dm1 pp3rve27;
	grainpool grains;
	vector<grGeomHdl*> sdf;
	vxlgrd thegrid;

	map<int, bool>incr2healthy;
	map<int, int> incr2rank;
	map<int, int> incr2wincr;

	//task-specific command line switch settings
	taskspecifics taskops;

	inline string get_prefix() {
		string prefix = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".Rank." + to_string(get_myrank()) + ".Incr." + to_string(Settings::IncrementFirst) + "." + to_string(Settings::IncrementOffset) + "." + to_string(Settings::IncrementLast);
		return prefix;
	}
	inline int get_myrank(){ return myrank; }
	inline int get_nranks(){ return nranks; }
	void set_myrank(const int r) { myrank = r; }
	void set_nranks(const int nr) { nranks = nr; }


	//results
	unsigned int thisincrement;
	unsigned int thiswrittenincrement;
	vector<rveAverageResults> avg;
	vector<rveAverageResults2> avg2;
	bv3x3 rveBaseIni;
	bv3x3 rveBaseDef;

	shape rveShape;

	//profiling
	//vector<plog> tictoc;
	profiler tictoc;
	//void spit_profiling();

private:
	//int which_db( const size_t ip, const size_t e);

	int myrank;
	int nranks;
	bool healthy;

	MPI_Datatype MPI_StatusInfo_Type;
	MPI_Datatype MPI_Tensor3x3_Double_Type;
	MPI_Datatype MPI_VonMises_Double_Type;
	MPI_Datatype MPI_DisloSpatDistr_Double_Type;
	MPI_Datatype MPI_Tensor3x3SpatDistr_Double_Type;
	MPI_Datatype MPI_GrainOriSpatDistr_Double_Type;
	MPI_Datatype MPI_GrainQuat_Double_Type;
};


#endif
