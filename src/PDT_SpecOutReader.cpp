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

#include "PDT_SpecOutReader.h"

inline bool SortDistByDistanceAsc(const dist &dist1 , const dist &dist2) {
	return dist1.d < dist2.d;
}


class diary
{
public:
	diary() { simid = 0; }
	diary(unsigned int _sid) : simid(_sid) {}
	~diary(){};

	unsigned int simid;
	vector<string> events;				//allowing to report what was done
};

ostream& operator << (ostream& in, diary const & val) {
	//append as in the postResults.py scripts which analysis tasks were conducted

	//##MK::implement a check which ticks internally what was really done...
	cout << "DAMASK PDT v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << endl;
	cout << "Utilizing Cartesian global coordinate system, right-handed, x,y,z!" << endl;
	in << "Task monitoring remains to be implemented #####" << endl;
	return in;
}


specmeta::specmeta()
{
}

specmeta::~specmeta()
{
}



mesh::mesh()
{
	owner = NULL;
	perOffset = bv3x3(); 	//identity unit cube
}


mesh::~mesh()
{
	//do not deallocate owner, only backreference
	//vector xyz0,dxyz_avg and dxyz_flu auto deconstruct
}


bool mesh::init()
{
	return true;
}


void mesh::eipid2xyz0()
{
	//what are the IDs to transform?
	size_t eipid_s = owner->eipid_start;
	size_t eipid_e = owner->eipid_end;
	size_t eipid_n = owner->eipid_n;

	//get memory to store their coordinates
	xyz0.reserve(eipid_n);

	//global mesh dimensions
/*	unsigned int dim[4] = {
			owner->owner->head.NX,
			owner->owner->head.NY,
			owner->owner->head.NZ,
			owner->owner->head.NXY }; //##MK::reduce size to speedup access

	real_xyz Ld[3] = { 	Settings::PhysicalDomainSize * owner->owner->head.L.x,
						Settings::PhysicalDomainSize * owner->owner->head.L.y,
						Settings::PhysicalDomainSize * owner->owner->head.L.z };*/

	size_t dim[4] = {
			static_cast<size_t>(owner->owner->head.NX),
			static_cast<size_t>(owner->owner->head.NY),
			static_cast<size_t>(owner->owner->head.NZ),
			static_cast<size_t>(owner->owner->head.NXY) };
	real_xyz scaler[3] = {
			static_cast<real_xyz>(1.0) / static_cast<real_xyz>(dim[0]),
			static_cast<real_xyz>(1.0) / static_cast<real_xyz>(dim[1]),
			static_cast<real_xyz>(1.0) / static_cast<real_xyz>(dim[2])  };

//cout << scaler[0] << ";" << scaler[1] << ";" << scaler[2] << endl;

	//perform mapping knowing in advance that indexing in DAMASK_spectral is x+y*NX+z*NX*NY
	//size_t Nnodes = static_cast<size_t>(dim[0]+1) * static_cast<size_t>(dim[1]+1) * static_cast<size_t>(dim[2]+1);
	//size_t Nelements = static_cast<size_t>(dim[0]) * static_cast<size_t>(dim[1]) * static_cast<size_t>(dim[2]);

/*
	size_t a = static_cast<size_t>(dim[0]+1); //how many nodes
	size_t b = static_cast<size_t>(dim[1]+1);
	size_t c = static_cast<size_t>(dim[2]+1);
	size_t ab = a*b;

	Ld[0] /= static_cast<real_xyz>(owner->owner->head.NX); //how many elements
	Ld[1] /= static_cast<real_xyz>(owner->owner->head.NY);
	Ld[2] /= static_cast<real_xyz>(owner->owner->head.NZ);
*/

	//real_xyz xx = 0.0;
	//real_xyz yy = 0.0;
	//real_xyz zz = 0.0;

	for ( size_t eipid = eipid_s; eipid < eipid_e; ++eipid ) {
		//##MK::mesh element integer cell coordinate
/*		size_t zz = eipid / (dim[3]);
		size_t rem = eipid - zz*dim[3];
		size_t yy = rem / dim[0];
		size_t xx = rem - yy*dim[0];
		xyz0.push_back( p3d(
				Ld[0]*static_cast<real_xyz>(xx),
				Ld[1]*static_cast<real_xyz>(yy),
				Ld[2]*static_cast<real_xyz>(zz)) );*/

		//##incorrect compared to DAMASK_spectral, implement in range 1/dim 0

		//##Debug for 117 type elements and non-nodal analyses only
		xyz0.push_back( p3d(
				0.5*scaler[0] + static_cast<real_xyz>( eipid         % dim[0]) * scaler[0],
				0.5*scaler[1] + static_cast<real_xyz>((eipid/dim[0]) % dim[1]) * scaler[1],
				0.5*scaler[2] + static_cast<real_xyz>((eipid/dim[3]) % dim[2]) * scaler[2] )
		);

/*
    	//##MK::DAMASK postResults node coordinate
		//xx = Ld[0] * static_cast<real_xyz>((eipid % a)); 		//[self.size[0] *       (n%a) / self.grid[0],
		//yy = Ld[1] * static_cast<real_xyz>(((eipid/a) % b)); 	// self.size[1] *   ((n/a)%b) / self.grid[1],
		//zz = Ld[2] * static_cast<real_xyz>(((eipid/ab) % c)); 	//self.size[2] * ((n/a/b)%c) / self.grid[2],

		xyz0.push_back( p3d(
				Ld[0] * static_cast<real_xyz>((eipid % a)),
				Ld[1] * static_cast<real_xyz>(((eipid/a) % b)),
				Ld[2] * static_cast<real_xyz>(((eipid/ab) % c))) ); //###-->this gives nodal coordinates not ip coordinates!
*/

//##MK::DEBUG
//p3d debug = xyz0.back();
//cout << debug.x << "\t\t\t" << debug.y << "\t\t\t" << debug.z << endl;
	}
}



homog::homog()
{
	owner = NULL;
}


homog::~homog()
{
	//do not deallocate owner, only backreference
	//vector gid auto deconstruct
}


bool homog::init()
{
/*	size_t n = owner->eip_n;
	try { eid = new unsigned int[n]; }
	catch (bad_alloc &exc) {
		stopping( "homog unable to allocate eid", 0, 0);
		return false;
	}
	eid_n = n;
	eid_nextfreeslot = 0;

	try { gid = new unsigned int[n]; }
	catch (bad_alloc &exc) {
		stopping( "homog unable to allocate gid", 0, 0);
		return false;
	}
	gid_n = n;
	gid_nextfreeslot = 0;*/

	return true;
}


cryst::cryst()
{
	owner = NULL;
}


cryst::~cryst()
{
	//do not deallocate owner, only backreference
	//vector PhaseID,TextureID,V,q,F,Fp auto deconstruct
}


bool cryst::init()
{
	return true;
}


real_xyz cryst::compute_total_V()
{
	//volume average values of deformation gradient tensor
	real_xyz tmp = 0.0;
	
	//##MK::OpenMP and SIMD

	size_t i = 0;
	size_t ni = V.size();
	for ( 	; i < ni; ++i ) { //summating
		tmp += V.at(i);
	}
	
	return tmp;
}

t3x3 cryst::compute_average_Fp()
{
	//volume average values of deformation gradient tensor
	t3x3 tmp = t3x3 (	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0) );

	//##MK::OpenMP and SIMD

	real_m33 w = static_cast<real_m33>(1.0); //##MK::strictly spreaking should be ip volume

	size_t i = 0;
	size_t ni = Fp.size();
	for ( 	; i < ni; ++i ) {
		tmp.add( Fp.at(i), w );
	}

	if ( ni > 0 ) {	//most likely averaging
		real_m33 scaler = static_cast<real_m33>(ni);
		tmp.div( scaler );
		return tmp;
	}
	//implicit else identity matrix
	return t3x3 (	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0) );
}


t3x3 cryst::compute_average_F()
{
	//volume average values of deformation gradient tensor
	t3x3 tmp = t3x3 (	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0) );

	//##MK::OpenMP and SIMD

	real_m33 w = static_cast<real_m33>(1.0); //##MK::strictly spreaking should be ip volume

	size_t i = 0;
	size_t ni = F.size();
	for ( 	; i < ni; ++i ) {
		tmp.add( F.at(i), w );
	}

	if ( ni > 0 ) {	//most likely averaging
		real_m33 scaler = static_cast<real_m33>(ni);
		tmp.div( scaler );
		return tmp;
	}
	//implicit else identity matrix
	return t3x3 (	static_cast<real_m33>(1.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(1.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(1.0) );
}


t3x3 cryst::compute_average_P()
{
	//volume average values of Piola Kirchhoff stress tensor
	//##MK::which Piola stress
	t3x3 tmp = t3x3 (	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
						static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0) );

	//##MK::OpenMP and SIMD
	real_m33 w = static_cast<real_m33>(1.0); //##MK::strictly spreaking should be ip volume

	size_t i = 0;
	size_t ni = P.size();
	for ( 	; i < ni; ++i ) { //summating
		tmp.add( P.at(i), w );
	}

	if ( ni > 0 ) {	//most likely averaging
		real_m33 scaler = static_cast<real_m33>(ni);
		tmp.div( scaler );
		return tmp;
	}
	//implicit else zero matrix
	return t3x3 (	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0),	static_cast<real_m33>(0.0) );
}


consti::consti()
{
	owner = NULL;

	rho_e_mult = 1;
	rho_d_mult = 1;
}


consti::~consti()
{
	//do not deallocate owner, only backreference
	//rho_e and rho_d auto destruct
}


bool consti::init()
{
	return true;
}


memRegion::memRegion()
{
	owner = NULL;

	eipid_start = 0;		//which elements we process first index
	eipid_end = 0;		//first after last index, i.e. [eip_start,eip_end[
	eipid_n = 0;
}


memRegion::~memRegion()
{
	//MK::owner is only a backreference
	//individual class objects have own destructors automatically executed
}


void memRegion::db_distr_homogenization( const double* raw, const size_t eipfirst, const size_t eipn )
{
}


/*
void memRegion::db_distr_crystallite( const double* raw, const size_t eipfirst, const size_t eipn )
{
	//raw, is a buffer of doubles, which by construction contains complete packages of homog+cryst+consti for en Elements <=> mesh element ips

	//efirst the first DAMASK mesh integration point for which raw buffers data
	//en the total number of complete raw buffer data packets

	if ( (eipfirst + eipn) < eipid_start ) { //(no index overlap, leftcase) raw contains data from ips with ids much smaller than what I should take care of
		return;
	}
	else if ( eipfirst >= eipid_end ) { //(no index overlap, rightcase) raw contains data from ip with ids much larger than what I should take care of
		return;
	} //for most threads these are the most likely entered cases, however we still perform the datapartitioning on the in parallel to allow the writing of the data in controllable threadlocal memory, which in particular for ccNUMA system is crucial
	else { //there is partial overlap, but I care only for the data from ip with ids within my id range to process i.e. (eip_start,eip_end(
		//we know that the first byte of raw begins always with data for an integration points not a fraction of data for one ip on two or more threads...
		//but that integration point may still to be stored in another memory region (partial index overlap) so first of all we
		//need to think in multiples of

		size_t dbls_per_eip = static_cast<size_t>(owner->head.matpres);

		//initially pointing to first byte of raw, offset is a byte offset
		size_t offset_start = 0 + ( (eipid_start - eipfirst) * dbls_per_eip);
		size_t offset_end = eipn * dbls_per_eip; //stopping after having read so many bytes, ie up to the end of the input buffer raw!
		//potentially however eip_next requires me to stop processing earlier as otherwise I would also process ids past eip_end(
		if ( (eipfirst + eipn) > eipid_end ) {
			//raw contains data past the range of damask integration point grid ids for which I have to take care of
			offset_end = (eipid_end - eipfirst) * dbls_per_eip;
			//#####MK???????????????????????????????????????????????????????????
		}

//cout << "cryst-->tid/dblspereip/offsetstart/offsetend/eipfirst/eipn = " << omp_get_thread_num() << ";" << dbls_per_eip << ";" << offset_start << ";" << offset_end << ";" << eipfirst << ";" << eipn << endl;

		//transfer values to database
		for ( size_t here = offset_start; here < offset_end; here = here + dbls_per_eip ) { //here is a byte offset, but we work in packages of head.matpres results
			//spectralOut files are read strictly forward iterating
			//thus we can pushback to db.arrays as integration points with lower ids that I take care of
			//then were already pushed back

			//###detailed data layout depends on case, evaluate owner->dlayout following the vectors forward iterating

			//###DEBUG
			size_t i = 0; //offset_start;
			i+=1; //skip n homo, knowing there is no homo
			i+=1; //1 ip per el
			i+=1; //not interested in Ncrystallite results

			//mind order of owner->dlayout. ...
			crystallite.PhaseID.push_back( static_cast<unsigned int>(raw[here+i]) ); i+=1;
//cout << crystallite.PhaseID.back() << endl;
			crystallite.TextureID.push_back( static_cast<unsigned int>(raw[here+i]) ); i+=1;
//cout << crystallite.TextureID.back() << endl;
			crystallite.V.push_back( static_cast<real_xyz>(raw[here+i]) ); i+=1;
//cout << crystallite.V.back() << endl;
			real_ori qtmp[4] = { 	static_cast<real_ori>(raw[here+i+0]),
									static_cast<real_ori>(raw[here+i+1]),
									static_cast<real_ori>(raw[here+i+2]),
									static_cast<real_ori>(raw[here+i+3])   };
			//##MK::at least even in DAMASK v2.0.1 internally quaternions specified actively interpreted rotations
			// will change with next release to become consistent with
			// D. Rowenhorst, A. D. Rollett, G. S. Rohrer, M. Groeber, M. Jackson, P. J. Konijnenberg, M. de Graef
			// Modelling and Simulation in Materials Science and Engineering 2015, Vol 23,
			// doi: 10.1088/0965-0393/23/8/083501
			// as such currently we have to reinterpret
			//
//cout << "before active quaternion \t\t" << qtmp[0] << ";" << qtmp[1] << ";" << qtmp[2] << ";" << qtmp[3] << endl;
/*
active2passive( qtmp );
//cout << "now passive quaternion \t\t" << qtmp[0] << ";" << qtmp[1] << ";" << qtmp[2] << ";" << qtmp[3] << endl;

			crystallite.q.push_back( quat( qtmp[0], qtmp[1], qtmp[2], qtmp[3] ) ); i+=4;

//cout << crystallite.q.back().q0 << ";" << crystallite.q.back().q1 << ";" << crystallite.q.back().q2 << ";" << crystallite.q.back().q3 << endl;
			i+=3; //skip euler
			i+=4; //skip rotation

			crystallite.F.push_back( t3x3(
					static_cast<real_m33>(raw[here+i+0]),
					static_cast<real_m33>(raw[here+i+1]),
					static_cast<real_m33>(raw[here+i+2]),
					static_cast<real_m33>(raw[here+i+3]),
					static_cast<real_m33>(raw[here+i+4]),
					static_cast<real_m33>(raw[here+i+5]),
					static_cast<real_m33>(raw[here+i+6]),
					static_cast<real_m33>(raw[here+i+7]),
					static_cast<real_m33>(raw[here+i+8])) ); i+=9;
//cout << crystallite.F.at(crystallite.F.size()-1) << endl;
			i+=9; //skip fe
			crystallite.Fp.push_back( t3x3(
					static_cast<real_m33>(raw[here+i+0]),
					static_cast<real_m33>(raw[here+i+1]),
					static_cast<real_m33>(raw[here+i+2]),
					static_cast<real_m33>(raw[here+i+3]),
					static_cast<real_m33>(raw[here+i+4]),
					static_cast<real_m33>(raw[here+i+5]),
					static_cast<real_m33>(raw[here+i+6]),
					static_cast<real_m33>(raw[here+i+7]),
					static_cast<real_m33>(raw[here+i+8])) ); i+=9;
			//nothing more we are interested in

			crystallite.P.push_back( t3x3(
					static_cast<real_m33>(raw[here+i+0]),
					static_cast<real_m33>(raw[here+i+1]),
					static_cast<real_m33>(raw[here+i+2]),
					static_cast<real_m33>(raw[here+i+3]),
					static_cast<real_m33>(raw[here+i+4]),
					static_cast<real_m33>(raw[here+i+5]),
					static_cast<real_m33>(raw[here+i+6]),
					static_cast<real_m33>(raw[here+i+7]),
					static_cast<real_m33>(raw[here+i+8])) ); i+=9;
//cout << crystallite.P.back() << endl;

		} //proceed with next element aka ip
	}
}
*/


void memRegion::db_distr_crystallite2( const double* raw, const size_t eipfirst, const size_t eipn )
{
	//raw, is a buffer of doubles, which by construction contains complete packages of homog+cryst+consti for en Elements <=> mesh element ips

	//efirst the first DAMASK mesh integration point for which raw buffers data
	//en the total number of complete raw buffer data packets

	if ( (eipfirst + eipn) < eipid_start ) { //(no index overlap, leftcase) raw contains data from ips with ids much smaller than what I should take care of
		return;
	}
	else if ( eipfirst >= eipid_end ) { //(no index overlap, rightcase) raw contains data from ip with ids much larger than what I should take care of
		return;
	} //for most threads these are the most likely entered cases, however we still perform the datapartitioning on the in parallel to allow the writing of the data in controllable threadlocal memory, which in particular for ccNUMA system is crucial
	else { //there is at least partial if not complete overlap,
		//but I care only for the data from ip with ids within my id range to process i.e. (eip_start,eip_end(
		//we know that the first byte of raw begins always with data for an integration point NOT a fraction of data for one ip on two or more threads...
		//but that integration point may still to be stored in another memory region (in case of, referred to partial index overlap) so first of all we
		//need to think in multiples of

		size_t dbls_per_eip = static_cast<size_t>(owner->head.matpres);

		//initially pointing to first byte of raw, offset is a byte offset
		size_t offset_start = 0; //assume first that raw contains data from ip ids completely inside (eipid_start,eipid_end(
		if (eipid_start >= eipfirst) { //if not initial portion of data on raw is for another memory
			offset_start = 0 + ( (eipid_start - eipfirst) * dbls_per_eip);
		}
		size_t offset_end = eipn * dbls_per_eip; //stopping after having read so many bytes, ie up to the end of the input buffer raw!
		//potentially however eip_next requires me to stop processing earlier as otherwise I would also process ids past eip_end(
		if ( (eipfirst + eipn) > eipid_end ) {
			//raw contains data past the range of damask integration point grid ids for which I have to take care of
			offset_end = (eipid_end - eipfirst) * dbls_per_eip;
			//#####MK???????????????????????????????????????????????????????????
		}

cout << "cryst-->tid/dblspereip/offsetstart/offsetend/eipfirst/eipn = " << omp_get_thread_num() << ";" << dbls_per_eip << ";" << offset_start << ";" << offset_end << ";" << eipfirst << ";" << eipn << endl;

		//transfer values to database
		for ( size_t here = offset_start; here < offset_end; here = here + dbls_per_eip ) { //here is a byte offset, but we work in packages of head.matpres results
			//spectralOut files are read strictly forward iterating
			//thus we can pushback to db.arrays as integration points with lower ids that I take care of
			//then were already pushed back

			//###detailed data layout depends on case, evaluate owner->dlayout following the vectors forward iterating

			//###DEBUG
			size_t i = 0; //offset_start;
			i+=1; //skip number of results for homogenization 
			i+=0; //skip homogenization knowing there is no homo ##MK::DEBUG
			i+=1; //skip how many grains per material point
			i+=1; //skip how many results for crystallite as we know the structure in more detail given the *.outputCrystallite file	
			
			//##MK::mind order of owner->dlayout. ...
			//##MK::mind order of owner->dlayout. ...
			//##MK::mind order of owner->dlayout. ...
			
			//MK::utilize that order in *.outputCrystallite labeling meta information file
			//encodes in strictly ascending order the internal structure of a single binary data block per eip
			unsigned int nkeys = owner->dlayout.crystallite.size();
			for ( unsigned int key = 0; key < nkeys; ++key) {
				
				string thiskey = owner->dlayout.crystallite.at(key).key;
				unsigned int len = owner->dlayout.crystallite.at(key).valui;
				
				if (thiskey.compare("phase") == 0) {
					crystallite.PhaseID.push_back( static_cast<unsigned int>(raw[here+i]) ); 
					i+=len; //1;
//cout << crystallite.PhaseID.back() << endl;
				}
				else if (thiskey.compare("texture") == 0) {
					crystallite.TextureID.push_back( static_cast<unsigned int>(raw[here+i]) ); 
					i+=len; //1;
//cout << crystallite.TextureID.back() << endl;
				}
				else if (thiskey.compare("volume") == 0) {
					crystallite.V.push_back( static_cast<real_xyz>(raw[here+i]) );
					i+=len; //1;
//cout << crystallite.V.back() << endl;				
				}
				else if (thiskey.compare("orientation") == 0) {
					real_ori qtmp[4] = { 	static_cast<real_ori>(raw[here+i+0]),
											static_cast<real_ori>(raw[here+i+1]),
											static_cast<real_ori>(raw[here+i+2]),
											static_cast<real_ori>(raw[here+i+3])   };
					//##MK::at least even in DAMASK v2.0.1 internally quaternions specified actively interpreted rotations
					/* will change with next release to become consistent with
					 * D. Rowenhorst, A. D. Rollett, G. S. Rohrer, M. Groeber, M. Jackson, P. J. Konijnenberg, M. de Graef
					 * Modelling and Simulation in Materials Science and Engineering 2015, Vol 23,
					 * doi: 10.1088/0965-0393/23/8/083501
					 * as such currently we have to reinterpret
					 */
//cout << "before active quaternion \t\t" << qtmp[0] << ";" << qtmp[1] << ";" << qtmp[2] << ";" << qtmp[3] << endl;
					active2passive( qtmp );
//cout << "now passive quaternion \t\t" << qtmp[0] << ";" << qtmp[1] << ";" << qtmp[2] << ";" << qtmp[3] << endl;

					crystallite.q.push_back( quat( qtmp[0], qtmp[1], qtmp[2], qtmp[3] ) ); 
//cout << "stored in crystallite \t\t" << crystallite.q.back().q0 << ";" << crystallite.q.back().q1 << ";" << crystallite.q.back().q2 << ";" << crystallite.q.back().q3 << endl;
					i+=len; //4;
				}
				else if (thiskey.compare("eulerangles") == 0) {
					i+=len; //3; //skip eulerangles
				}
				else if (thiskey.compare("grainrotation") == 0) {
					i+=len; //4; //skip rotation
				}
				else if (thiskey.compare("f") == 0) {
					crystallite.F.push_back( t3x3(
							static_cast<real_m33>(raw[here+i+0]),
							static_cast<real_m33>(raw[here+i+1]),
							static_cast<real_m33>(raw[here+i+2]),
							static_cast<real_m33>(raw[here+i+3]),
							static_cast<real_m33>(raw[here+i+4]),
							static_cast<real_m33>(raw[here+i+5]),
							static_cast<real_m33>(raw[here+i+6]),
							static_cast<real_m33>(raw[here+i+7]),
							static_cast<real_m33>(raw[here+i+8])) ); 
					i+=len; //9;
//cout << crystallite.F.back() << endl;					
				}
				else if (thiskey.compare("fe") == 0) {
					i+=len; //9; //skip fe					
				}
				else if (thiskey.compare("fp") == 0) {
					crystallite.Fp.push_back( t3x3(
						static_cast<real_m33>(raw[here+i+0]),
						static_cast<real_m33>(raw[here+i+1]),
						static_cast<real_m33>(raw[here+i+2]),
						static_cast<real_m33>(raw[here+i+3]),
						static_cast<real_m33>(raw[here+i+4]),
						static_cast<real_m33>(raw[here+i+5]),
						static_cast<real_m33>(raw[here+i+6]),
						static_cast<real_m33>(raw[here+i+7]),
						static_cast<real_m33>(raw[here+i+8])) );
					i+=len; //9;
				}
				else if (thiskey.compare("e") == 0) {
					i+=len; //9;
				}
				else if (thiskey.compare("ee") == 0) {
					i+=len; //9;
				}
				else if (thiskey.compare("p") == 0) {
					crystallite.P.push_back( t3x3(
						static_cast<real_m33>(raw[here+i+0]),
						static_cast<real_m33>(raw[here+i+1]),
						static_cast<real_m33>(raw[here+i+2]),
						static_cast<real_m33>(raw[here+i+3]),
						static_cast<real_m33>(raw[here+i+4]),
						static_cast<real_m33>(raw[here+i+5]),
						static_cast<real_m33>(raw[here+i+6]),
						static_cast<real_m33>(raw[here+i+7]),
						static_cast<real_m33>(raw[here+i+8])) );
					i+=len; //9;
//cout << crystallite.P.back() << endl;
				}
				else if (thiskey.compare("lp") == 0) {
					i+=len; //9;
				}
				else { //nothing more we are interested in
					i+=0;
				}				
			} //collect values for next keyword

			//add dummy value for derived quantity
			crystallite.GrainID.push_back( numeric_limits<unsigned int>::max() );
		} //proceed with next element aka ip
	}
}


/*
void memRegion::db_distr_constitutive( const double* raw, const size_t eipfirst, const size_t eipn )
{
	//see detailed explanation in db_distr_crystallite
	if ( (eipfirst + eipn) < eipid_start ) { return; }
	else if ( eipfirst >= eipid_end ) { return; }
	else {
		size_t dbls_per_eip = static_cast<size_t>(owner->head.matpres);

		//initially pointing to first byte of raw, offset is a byte offset
		size_t offset_start = 0 + ( (eipid_start - eipfirst) * dbls_per_eip);
		size_t offset_end = eipn * dbls_per_eip;
		if ( (eipfirst + eipn) > eipid_end ) {
			offset_end = (eipid_end - eipfirst) * dbls_per_eip;
		}

cout << "const-->tid/dblspereip/offsetstart/offsetend/eipfirst/eipn = " << omp_get_thread_num() << ";" << dbls_per_eip << ";" << offset_start << ";" << offset_end << ";" << eipfirst << ";" << eipn << endl;

		//set multiplier
		constitutive.rho_e_mult = owner->dlayout.constitutive_getvalui( "edge_density" );
		constitutive.rho_d_mult = owner->dlayout.constitutive_getvalui( "dipole_density" );

//cout << "emult/dmult " << constitutive.rho_e_mult << ";" << constitutive.rho_d_mult << endl;

		//transfer values to database
		for ( size_t here = offset_start; here < offset_end; here = here + dbls_per_eip ) {

			size_t i = 0; //offset_start;

			size_t homogres = static_cast<size_t>(raw[here+i]);
//cout << i << "\t\t" << homogres << endl;
			i+=homogres; //skip homogenization completely
			i+=1;
			i+=1; //1 ip per el
			size_t crystres = static_cast<size_t>(raw[here+i]);
			i+=crystres; //skip Ncrystallite results completelyted in Ncrystallite results
//cout << i << "\t\t" << crystres << endl;
			i+=1;

			//beginning of constitutive section
			i+=1; //skip how many results there are

			//keep skipping until found what desired
			for (unsigned int r = 0; r < owner->dlayout.constitutive.size(); ++r) { //order of dlayout has sequence as the corresponding data in the raw buffer
				if (owner->dlayout.constitutive.at(r).key.compare("edge_density") == 0) {
					for (unsigned int sys = 0; sys < constitutive.rho_e_mult; ++sys) {
						constitutive.rho_e.push_back(static_cast<real_rho>(raw[here+i+sys]));
//cout << constitutive.rho_e.back() << "__";
					}
//cout << endl << endl;
					i+= constitutive.rho_e_mult;
				}
				if (owner->dlayout.constitutive.at(r).key.compare("dipole_density") == 0) {
					for (unsigned int sys = 0; sys < constitutive.rho_d_mult; ++sys) {
						constitutive.rho_d.push_back(static_cast<real_rho>(raw[here+i+sys]));
//cout << constitutive.rho_d.back() << "__";
					}
//cout << endl;
					i+= constitutive.rho_d_mult;
				}
				//anything else we dont care
			}

//break;
		} //proceed with next element aka ip
	}
}
*/


void memRegion::db_distr_constitutive2( const double* raw, const size_t eipfirst, const size_t eipn )
{
	//see detailed explanation in db_distr_crystallite
	if ( (eipfirst + eipn) < eipid_start ) { return; }
	else if ( eipfirst >= eipid_end ) { return; }
	else {
		size_t dbls_per_eip = static_cast<size_t>(owner->head.matpres);

		//initially pointing to first byte of raw, offset is a byte offset
		size_t offset_start = 0;
		if (eipid_start >= eipfirst) {
			offset_start = 0 + ( (eipid_start - eipfirst) * dbls_per_eip);
		}
		size_t offset_end = eipn * dbls_per_eip;
		if ( (eipfirst + eipn) > eipid_end ) {
			offset_end = (eipid_end - eipfirst) * dbls_per_eip;
		}

cout << "const-->tid/dblspereip/offsetstart/offsetend/eipfirst/eipn = " << omp_get_thread_num() << ";" << dbls_per_eip << ";" << offset_start << ";" << offset_end << ";" << eipfirst << ";" << eipn << endl;

		//set multiplier
		constitutive.rho_e_mult = owner->dlayout.constitutive_getvalui( "edge_density" );
		constitutive.rho_d_mult = owner->dlayout.constitutive_getvalui( "dipole_density" );

cout << "emult/dmult " << constitutive.rho_e_mult << ";" << constitutive.rho_d_mult << endl;

		//transfer values to database
		for ( size_t here = offset_start; here < offset_end; here = here + dbls_per_eip ) {

			size_t i = 0; //offset_start;
			size_t homogres = static_cast<size_t>(raw[here+i]);
			i+=1; //skip number of results for homogenization
//cout << i << "\t\t" << homogres << endl;
			i+=homogres; //skip homogenization completely

			i+=1; //skip how many grains per material point
			size_t crystres = static_cast<size_t>(raw[here+i]);
			i+=1;
//cout << i << "\t\t" << crystres << endl;
			i+=crystres; //skip Ncrystallite results completelyted in Ncrystallite results

			size_t constres = static_cast<size_t>(raw[here+i]); //beginning of constitutive section
//cout << i << "\t\t" << constres << endl;
			i+=1; //skip how many results there are

			if ( constres > 0) { //read if any
				//keep skipping until found what desired
				for (unsigned int r = 0; r < owner->dlayout.constitutive.size(); ++r) { //order of dlayout has sequence as the corresponding data in the raw buffer
					if (owner->dlayout.constitutive.at(r).key.compare("edge_density") == 0) {
						for (unsigned int sys = 0; sys < constitutive.rho_e_mult; ++sys) {
							constitutive.rho_e.push_back(static_cast<real_rho>(raw[here+i+sys]));
//cout << constitutive.rho_e.back() << "__";
						}
//cout << endl << endl;
						i+= constitutive.rho_e_mult;
					}
					if (owner->dlayout.constitutive.at(r).key.compare("dipole_density") == 0) {
						for (unsigned int sys = 0; sys < constitutive.rho_d_mult; ++sys) {
							constitutive.rho_d.push_back(static_cast<real_rho>(raw[here+i+sys]));
//cout << constitutive.rho_d.back() << "__";
						}
//cout << endl;
						i+= constitutive.rho_d_mult;
					}
					//anything else we dont care
				}
			}
			else {
				continue;
			}
//break;
		} //proceed with next element aka ip
	}
}



rveAverageResults::rveAverageResults(){
	
	Vtotal = static_cast<real_xyz>(0.0);
	
	Fpavgrve = t3x3( 		static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0) );

	Favgrve = t3x3( 		static_cast<real_m33>(1.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(1.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(1.0) );

	Pavgrve = t3x3( 		static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0) );

	Strainavgrve = t3x3( 	static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0) );

	Cauchyavgrve = t3x3( 	static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
							static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0) );

	Equivavgrve = vMises( 	static_cast<real_m33>(0.0), static_cast<real_m33>(0.0) );

	loadcaseID = numeric_limits<unsigned int>::max();
	localincrID = numeric_limits<unsigned int>::max();
	globalincrID = numeric_limits<unsigned int>::max();
}

rveAverageResults::rveAverageResults(const real_xyz Vtot, t3x3 const & Fp, t3x3 const & F, t3x3 const & P,
		t3x3 const & eps, t3x3 const & cauchy, vMises const & vm,
		const unsigned int lcid, const unsigned int lincrid, const unsigned int gincrid) {
			
	Vtotal = Vtot;
	
	Fpavgrve = Fp;

	Favgrve = F;
	Pavgrve = P;
	Strainavgrve = eps;
	Cauchyavgrve = cauchy;

	/*Favgrve.a11 = F.a11;
	Favgrve.a12 = F.a12;
	Favgrve.a13 = F.a13;
	Favgrve.a21 = F.a21;
	Favgrve.a22 = F.a22;
	Favgrve.a23 = F.a23;
	Favgrve.a31 = F.a31;
	Favgrve.a32 = F.a32;
	Favgrve.a33 = F.a33;

	Pavgrve.a11 = P.a11;
	Pavgrve.a12 = P.a12;
	Pavgrve.a13 = P.a13;
	Pavgrve.a21 = P.a21;
	Pavgrve.a22 = P.a22;
	Pavgrve.a23 = P.a23;
	Pavgrve.a31 = P.a31;
	Pavgrve.a32 = P.a32;
	Pavgrve.a33 = P.a33;

	Strainavgrve.a11 = eps.a11;
	Strainavgrve.a12 = eps.a12;
	Strainavgrve.a13 = eps.a13;
	Strainavgrve.a21 = eps.a21;
	Strainavgrve.a22 = eps.a22;
	Strainavgrve.a23 = eps.a23;
	Strainavgrve.a31 = eps.a31;
	Strainavgrve.a32 = eps.a32;
	Strainavgrve.a33 = eps.a33;

	Cauchyavgrve.a11 = cauchy.a11;
	Cauchyavgrve.a12 = cauchy.a12;
	Cauchyavgrve.a13 = cauchy.a13;
	Cauchyavgrve.a21 = cauchy.a21;
	Cauchyavgrve.a22 = cauchy.a22;
	Cauchyavgrve.a23 = cauchy.a23;
	Cauchyavgrve.a31 = cauchy.a31;
	Cauchyavgrve.a32 = cauchy.a32;
	Cauchyavgrve.a33 = cauchy.a33;*/

	Equivavgrve.vMisesEquivStrain = vm.vMisesEquivStrain;
	Equivavgrve.vMisesEquivStress = vm.vMisesEquivStress;

	loadcaseID = lcid;
	localincrID = lincrid;
	globalincrID = gincrid;
}

/*
rveAverageResults::rveAverageResults(const real_xyz Vtot,
		MPI_Tensor3x3_Double const & Fp,
		MPI_Tensor3x3_Double const & F,
		MPI_Tensor3x3_Double const & P,
		MPI_Tensor3x3_Double const & eps,
		MPI_Tensor3x3_Double const & cauchy,
		vMises const & vm,
		const unsigned int lcid, const unsigned int lincrid, const unsigned int gincrid) {

	Vtotal = Vtot;

	Fpavgrve = Fp;
	Favgrve = F;
	Pavgrve = P;
	Strainavgrve = eps;
	Cauchyavgrve = cauchy;

	Equivavgrve.vMisesEquivStrain = vm.vMisesEquivStrain;
	Equivavgrve.vMisesEquivStress = vm.vMisesEquivStress;

	loadcaseID = lcid;
	localincrID = lincrid;
	globalincrID = gincrid;
}
*/

rveAverageResults::~rveAverageResults()
{
}


specOutHeader::specOutHeader()
{
	owner = NULL;

	load = "";
	wdir = "";
	geom = "";

	N = 0;
	Nip = 1;
	Ncp = N*N*N;

	matpres = 0;
	loadcases = 0;
	//freqs = 0;
	//times = 0;
	//logscales = 0;
	//nincr = 0;
	sincr = 0;
	lincr = 0;
	NX = N;
	NY = N;
	NZ = N;
	NXY = NX*NY;
	NXYZ = NXY*NZ;

	L = v3x1(0.0, 0.0, 0.0);
	DataElementsPerIncrement = 0;
	FirstByteAfterHeader = 0;
};


specOutHeader::~specOutHeader()
{
	//MK::owner is only a backreference
}



bvh_xyzm2::bvh_xyzm2()
{
	owner = NULL; //MK::owner is only a backreference

	owin_deformed = aabb3d();
	owin_final = aabb3d();
	spatialdistr = sqb();

	healthy = true;
}


bvh_xyzm2::~bvh_xyzm2()
{
	//MK::owner is only a backreference
	for (size_t mr = 0; mr < pxyz.size(); mr++) {
		if ( pxyz.at(mr) != NULL ) {
			delete pxyz.at(mr); pxyz.at(mr) = NULL;
		}
	}
	for (size_t mr = 0; mr < pid.size(); ++mr) {
		if ( pid.at(mr) != NULL ) {
			delete pid.at(mr); pid.at(mr) = NULL;
		}
	}
	for (size_t mr = 0; mr < pimage.size(); ++mr) {
		if ( pimage.at(mr) != NULL ) {
			delete pimage.at(mr); pimage.at(mr) = NULL;
		}
	}
}


void bvh_xyzm2::aabb3d_about_deformed()
{
	//compute position of all elements and their periodic images, determine axis-aligned bounding box about this volume
	aabb3d glolimits = aabb3d();
cout << "glolimits " << endl << glolimits << endl;

	#pragma omp parallel shared(glolimits)
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		aabb3d mylimits = aabb3d(); //thread-local min-max

		//###encapsulate this in a function executed by object memRegion such that memRegion returns

		//work partitioning
		memRegion* thisregion = owner->db.at(mt);

		//processing of thread-local data, other tids work on other data
		for (size_t e = 0; e < thisregion->eipid_n; ++e) {
			//##MK::do this with SIMD
			real_xyz xx = thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx; //initial position of element + deformation-induced average + fluctuant displacement
			real_xyz yy = thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy;
			real_xyz zz = thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz;

			mylimits.potentially_expand( xx, yy, zz);
		}

		//thread-individual reduction moidifying global data requires is critical
		#pragma omp critical
		{
			glolimits.potentially_expand(mylimits);
		}
	} //end of parallel region

	glolimits.scale();
	//abslimits.xsz = abslimits.xmx - abslimits.xmi;
	//abslimits.ysz = abslimits.ymx - abslimits.ymi;
	//abslimits.zsz = abslimits.zmx - abslimits.zmi;

	owin_deformed = glolimits;

	string descr = "AABB0AboutDeformed";
	owner->rveShape.add( glolimits, descr );

cout << "owin_deformed " << endl << owin_deformed << endl;
}


void bvh_xyzm2::aabb3d_add_guardzone()
{
	//strictly speaking a convex hull about the elements in the deformed configuration surplus a layer of Settings::kernelRadius
	//thickness is required only, hence the present choice works with a (slightly) to large point set
	//however, then we have to pay the convex hull construction overhead surplus the point in convex hull inclusion test
	aabb3d guarded = aabb3d(
			owin_deformed.xmi - Settings::KernelRadius,
			owin_deformed.xmx + Settings::KernelRadius,
			owin_deformed.ymi - Settings::KernelRadius,
			owin_deformed.ymx + Settings::KernelRadius,
			owin_deformed.zmi - Settings::KernelRadius,
			owin_deformed.zmx + Settings::KernelRadius );

	guarded.scale();

	owin_final = guarded;

	string descr = "AABB1AboutDeformed";
	owner->rveShape.add( guarded, descr );

cout << "owin_final " << endl << owin_final << endl;
}


void bvh_xyzm2::spatialpartitioning_init()
{
	//MK//##MK::implement functionality that grid is build only once and optimized for the its purpose!

	//##MK::heuristical spatial binning of all element positions and their periodic images inside owin_final
	//such that it requires on average approximately 27 bins only to test for finding all neighboring elements within KernelRadius to any element
	//the memory is allocated and first touched within a group of threads, these distribute the bins in xy slabs across z
	//and populate them individually

	//##MK::load-balancing is expected fairly even because the DAMASK meshes are rather regular and not adaptive
	//##MK::compromise load balancing will require number of bins along z as an integer multiple of omp_get_num_threads
	//##MK::SIMD requires number of points as an integer multiple of the SIMD width
	//##MK::in addition for strongly deformed RVEs resulting in trapezoidal shapes the aabb does not wrap the element locations very tightly
	//##MK::resulting in a considerable number of close to the boundary periodic images which will never be candidates

	real_xyz dx = owin_final.xsz;
	real_xyz dy = owin_final.ysz;
	real_xyz dz = owin_final.zsz;

	//minimum dimension of owin_final
	real_xyz dmin = min(min(dx,dy),dz);
cout << "dmin " << dmin << endl;

	//with density of elements in owin_deformed, estimate partitioning

	size_t nx = static_cast<unsigned int>( (1.0 + static_cast<real_xyz>(INCREMENTAL_KERNEL_RADIUS_STEPPING))*((dmin / Settings::KernelRadius) * dx) + 1.0 + EPSILON); //requires Settings::KernelRadius as a generalize radius
	size_t ny = static_cast<unsigned int>( (1.0 + static_cast<real_xyz>(INCREMENTAL_KERNEL_RADIUS_STEPPING))*((dmin / Settings::KernelRadius) * dy) + 1.0 + EPSILON);
	size_t nz = static_cast<unsigned int>( (1.0 + static_cast<real_xyz>(INCREMENTAL_KERNEL_RADIUS_STEPPING))*((dmin / Settings::KernelRadius) * dz) + 1.0 + EPSILON);

	spatialdistr = sqb( nx, ny, nz );

cout << "SpatialDistribution " << endl << spatialdistr << endl;

	//#pragma omp parallel
	//{
		//int tid = omp_get_thread_num();
		//int tid = MASTER;
		//##MK::generation of element position buckets in thread-local memory
		bool mehealthy = true;

		//preallocate
		pxyz.reserve(spatialdistr.nxyz);	//x,y,z coordinates of an ip point or its periodic image
		pid.reserve(spatialdistr.nxyz); //the global ID of the integration point within [0,NXYZ) //##MK::given the specific DEBUG case of 117 element single ip mesh elements so far
		//MK::also the periodic images of ips inside the container that the bvh_xyzm2 partitions spatially have the IDs of the generating ip i.e. are within [0,NXYZ)
		pimage.reserve(spatialdistr.nxyz); //a character code identifying which of the Moore periodic images it is

		for (size_t b = 0; b < spatialdistr.nxyz; ++b) {
			pxyz.push_back(NULL);	//individual arrays for coordinates to facilitate SIMD processing, ##MK may be a tradeoff as compared to bundling x,y, z in struct
			pid.push_back(NULL);
			pimage.push_back(NULL);
		}

		for (size_t b = 0; b < spatialdistr.nxyz; ++b) {
			vector<p3d>* xyzbucket = NULL;
			try { xyzbucket = new vector<p3d>; }
			catch (bad_alloc &croak) {
				stopping("Unable to allocate xyzbucket when building bounded_volume_hierarchy", owner->get_myrank(), MASTER );
				mehealthy = false; break;
			}
			pxyz.at(b) = xyzbucket;

			vector<unsigned int>* idbucket = NULL;
			try { idbucket = new vector<unsigned int>; }
			catch (bad_alloc &croak) {
				stopping("Unable to allocate idbucket when building bounded_volume_hierarchy", owner->get_myrank(), MASTER );
				mehealthy = false; break;
			}
			pid.at(b) = idbucket;

			//##MK::DEBUG
			vector<unsigned char>* imgbucket = NULL;
			try { imgbucket = new vector<unsigned char>; }
			catch (bad_alloc &croak) {
				stopping("Unable to allocate imgbucket when building bounded_volume_hierarchy", owner->get_myrank(), MASTER );
				mehealthy = false; break;
			}
			pimage.at(b) = imgbucket;
		} //all buckets initalized

		//#pragma omp critical
		//{
			if (mehealthy == false) {
				healthy = false;
			}
		//}
	//} //end of parallel region
	reporting("Spatialpartitioning initialized", owner->get_myrank(), 0, false);
}


inline unsigned int bvh_xyzm2::binning_x( const real_xyz x )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (x - owin_final.xmi) / (owin_final.xmx - owin_final.xmi) * static_cast<real_xyz>(spatialdistr.nx);
	unsigned int tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int xx = ( tmpui < spatialdistr.nx ) ? tmpui : spatialdistr.nx-1;
	return xx;
}

inline unsigned int bvh_xyzm2::binning_y( const real_xyz y )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (y - owin_final.ymi) / (owin_final.ymx - owin_final.ymi) * static_cast<real_xyz>(spatialdistr.ny);
	unsigned int tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int yy = ( tmpui < spatialdistr.ny ) ? tmpui : spatialdistr.ny-1;
	return yy;
}

inline unsigned int bvh_xyzm2::binning_z( const real_xyz z )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (z - owin_final.zmi) / (owin_final.zmx - owin_final.zmi) * static_cast<real_xyz>(spatialdistr.nz);
	unsigned int tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int zz = ( tmpui < spatialdistr.nz ) ? tmpui : spatialdistr.nz-1;
	return zz;
}



inline unsigned int bvh_xyzm2::binning_xyz( const p3d candidate )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
/*	real_xyz tmpreal = (candidate.x - owin_final.xmi) / (owin_final.xmx - owin_final.xmi) * static_cast<real_xyz>(spatialdistr.nx);
	unsigned int tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int x = ( tmpui < spatialdistr.nx ) ? tmpui : spatialdistr.nx-1;*/
	unsigned int x = binning_x( candidate.x );

/*	tmpreal = (candidate.y - owin_final.ymi) / (owin_final.ymx - owin_final.ymi) * static_cast<real_xyz>(spatialdistr.ny);
	tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int y = ( tmpui < spatialdistr.ny ) ? tmpui : spatialdistr.ny-1;*/
	unsigned int y = binning_y( candidate.y );

/*	tmpreal = (candidate.z - owin_final.zmi) / (owin_final.zmx - owin_final.zmi) * static_cast<real_xyz>(spatialdistr.nz);
	tmpui = static_cast<real_xyz>(tmpreal);
	unsigned int z = ( tmpui < spatialdistr.nz ) ? tmpui : spatialdistr.nz-1;*/
	unsigned int z = binning_z( candidate.z );

	return (x + y*spatialdistr.nx + z*spatialdistr.nxy);
}


//inline bool included( real_xyz x, real_xyz y, real_xyz z, bv3x3 const & rvedefbase, int h, int k, int l, aabb3d const & window) {}

void bvh_xyzm2::spatialpartitioning_deformedrve_perform()
{
	//forall elements and all their periodic images, check if in owin_final, if so bin them
	if ( healthy == false) {
		return;
	}

	//##MK::
	bv3x3 debug = owner->db.at(0)->grid.perOffset;
	owner->rveShape.add( debug, "TriclinicRVE");

	//thread-parallelized computing of included deformed configuration element positions and their periodic images surplus their target bin
	#pragma omp parallel
	{
		//##MK::one might optimize further this function for generating less duplicate perImages
		//as the threads which take care of ip slabs which are neither the upper nor lower slab generate multiply several points at the boundary?

		//thread-local temporary memory mybuffer to store coordinates metadata id and target bin to allow fully parallel periodic images, checking, and binning
		//with minimal concurrence when writing to global datastructure
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		bool mehealthy = true;

		vector<vector<p3dm3>*> mybuffer;
		vector<p3dm3>* s = NULL;
		try { s = new vector<p3dm3>; }
		catch (bad_alloc &exc) {
			stopping("Unable to allocate binning memory", owner->get_myrank(), mt);
			mehealthy = false;
		}
		mybuffer.push_back(s);

		if ( mehealthy == true ) {
			memRegion* thisregion = owner->db.at(mt);
			size_t eipn = thisregion->eipid_n;
			size_t lc2gl_idoffset = thisregion->eipid_start;
			bv3x3 per = thisregion->grid.perOffset; //tricline deformed rve base column vectors
			aabb3d w = owin_final;

			p3d cache[27]; //cache periodic images values of a single element
			unsigned int targetbin[27]; //cache bins into which the points get stored
			unsigned char periodicImage[27]; //cache bins which of the Moore periodic images it is, 0x00 is the center, 1 is -1-1-1, 2 0-1-1, and so forth

			unsigned int i = 0;
			real_xyz xxx = 0.0;
			real_xyz yyy = 0.0;
			real_xyz zzz = 0.0;
			real_xyz h = 0.0;
			real_xyz k = 0.0;
			real_xyz l = 0.0;

			for (size_t e = 0; e < eipn; ++e) { //MK::e is the threadlocal iterating element id, NOT the global ip ID !

				real_xyz xx = thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx;	//##MK::SIMD vectorize location of element in deformed configuration
				real_xyz yy = thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy;
				real_xyz zz = thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz;
				i = 0; //defining at which position you buffer these intermediate resutls in the caches cahce, targetpin and periodicImage

				//most points lay inside
				//<000>
				h = +0.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x00; i++; }

				//<100>
				//[+1+0+0]
				h = +1.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0e; i++; }
				//[-1+0+0]
				h = -1.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0d; i++; }
				//[+0+1+0]
				h = +0.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x10; i++; }
				//[+0-1+0]
				h = +0.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0b; i++; }
				//[+0+0+1]
				h = +0.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x16; i++; }
				//[+0+0-1]
				h = +0.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x05; i++; }

				//<110>
				//[+1+1+0]
				h = +1.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x11; i++; }
				//[-1+1+0]
				h = -1.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0f; i++; }
				//[+1-1+0]
				h = +1.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0c; i++; }
				//[-1-1+0]
				h = -1.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0a; i++; }
				//[+0+1+1]
				h = +0.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x19; i++; }
				//[+0-1+1]
				h = +0.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x13; i++; }
				//[+0+1-1]
				h = +0.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x08; i++; }
				//[+0-1-1]
				h = +0.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x02; i++; }
				//[+1+0+1]
				h = +1.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x17; i++; }
				//[-1+0+1]
				h = -1.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x15; i++; }
				//[+1+0-1]
				h = +1.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x06; i++; }
				//[-1+0-1]
				h = -1.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x04; i++; }

				//<111>
				//[+1+1+1]
				h = +1.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x1a; i++; }
				//[-1+1+1]
				h = -1.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x18; i++; }
				//[+1-1+1]
				h = +1.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x14; i++; }
				//[+1+1-1]
				h = +1.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x09; i++; }
				//[+1-1-1]
				h = +1.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x03; i++; }
				//[-1+1-1]
				h = -1.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x07; i++; }
				//[-1-1+1]
				h = -1.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x12; i++; }
				//[-1-1-1]
				h = -1.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x01; i++; }

				//now that we know how many periodic copies of e there are at all inside namely i we bin them
				for (unsigned int c = 0; c < i; ++c) {
					targetbin[c] = binning_xyz( cache[c] );
//cout << targetbin[c] << "__";
				}
//cout << endl;

				//check if still i places free in buffer construct, if so store results
				if ( mybuffer.back()->size() < (Performance::ThreadBufferAllocSize - static_cast<size_t>(i)) ) {
					vector<p3dm3>* drophere = mybuffer.back();
					unsigned int eid = static_cast<unsigned int>(e + lc2gl_idoffset);
					for (unsigned int c = 0; c < i; ++c) {
						drophere->push_back( p3dm3(cache[c].x, cache[c].y, cache[c].z,
								eid, targetbin[c], periodicImage[c] ) );
						//MK::mark m1 is the global ID of the ip,
						//m2 the ID of the bin in which the point should be stored,
						//m3 metadata telling us in which periodic BC relation the point is to its generating ip
					}
				} //existent buffer filled
				else { //get new buffer
					try {
						s = NULL;
						s = new vector<p3dm3>;
						mybuffer.push_back(s);
					}
					catch (bad_alloc &ompcroak) {
						stopping("Unable to reallocate additional threadlocal binning memory", owner->get_myrank(), mt);
						mehealthy = false;
						break; //stop working on this thread
					}

					vector<p3dm3>* dropwhere = mybuffer.back();
					unsigned int eid = static_cast<unsigned int>(e + lc2gl_idoffset); //MK::ip base and its periodic images have the same global ip ID eid
					for (unsigned int c = 0; c < i; ++c) {
						dropwhere->push_back( p3dm3(cache[c].x, cache[c].y, cache[c].z,
								eid, targetbin[c], 0x00 ) );
					}
				} //new buffer utilized for filling
			} //process next element e of thread
		} //while being healthy
		else { //upon initial allocation error
			stopping("Unable to initially allocate threadlocal binning memory", owner->get_myrank(), mt);
			mehealthy = false;
		}

		#pragma omp critical
		{
			if ( mehealthy == false ) {
				healthy = false;
				reporting("Bounding volume hierarchy construction threadlocal results faulty or incomplete", owner->get_myrank(), mt, true);
			}
			else {
				reporting("Bounding volume hierarchy construction threadlocal results ready", owner->get_myrank(), mt, true);
			}
		}

		//wait for all finishing the positioning and binning of the periodic images in w
		#pragma omp barrier


		//by now the threads have positioned, checked for inclusion, and binned all mesh element integration points surplus their periodic images points locally so fuse into global
		//##MK::for substantially large point data use incremental updating of global bvh to reduce total resistent memory
		#pragma omp critical
		{
			if ( mehealthy == false ) { //in case previous memory allocation was unsuccessful
				for (size_t b = 0; b < mybuffer.size(); ++b) {
					if (mybuffer.at(b) != NULL) {
						delete mybuffer.at(b); mybuffer.at(b) = NULL;
					}
				}
			}
			else { //valid data, so pass them into global bvh_xyzm2
				for (size_t b = 0; b < mybuffer.size(); ++b) {
					if (mybuffer.at(b) != NULL) {
						vector<p3dm3>* these = mybuffer.at(b);
						size_t nj = these->size();
						for (size_t j = 0; j < nj; ++j) {
							//MK::m1 is the global ID of the ip,
							//m2 the ID of the bin in which the point should be stored,
							//m3 metadata telling us in which periodic BC relation the point is to its generating ip
							p3dm3 piece = these->at(j);
							pxyz.at(piece.m2)->push_back( p3d(piece.x, piece.y, piece.z) );
							pid.at(piece.m2)->push_back( piece.m1 );
							//##MK::DEBUG which Moore periodic image is this?
							pimage.at(piece.m2)->push_back( piece.m3 );
						}
						//release threadlocal temporary mybuffer
						delete mybuffer.at(b); mybuffer.at(b) = NULL;
					}
				} //copy over next group of elements to process
			}
		} //end of pragma critical, thread tid has communicated local results

	} //end of parallel region

	if ( healthy == true )
		reporting("Bounding volume hierarchy construction complete", owner->get_myrank(), 0, false);
	else
		reporting("Bounding volume hierarchy construction failed", owner->get_myrank(), 0, false);
}


void bvh_xyzm2::build_bvh_xyzm2()
{
	aabb3d_about_deformed();
	aabb3d_add_guardzone();
	spatialpartitioning_init();
	spatialpartitioning_deformedrve_perform();

	if ( Settings::VisIPGridWithPerImages > 0 ) {
		bool mystatus = vtk_p3dm3( pxyz, pid, pimage, owner->thisincrement);
	}
}


void bvh_xyzm2::find_higherorder_neighbors( const p3d p, vector<dist> & results, real_xyz r )
{
	//MK:: r <= Settings::KernelRadius !
	real_xyz CurrentKernelRadius = Settings::KernelRadius;
	if ( r < CurrentKernelRadius && r > 0.0)
		CurrentKernelRadius = r;

	//clear potentially existent values
	results.clear();

	unsigned int wxmi = binning_x( p.x - CurrentKernelRadius ); //generalized coordinate in RVE =- generalized distance
	unsigned int wxmx = binning_x( p.x + CurrentKernelRadius );
	unsigned int wymi = binning_y( p.y - CurrentKernelRadius );
	unsigned int wymx = binning_y( p.y + CurrentKernelRadius );
	unsigned int wzmi = binning_z( p.z - CurrentKernelRadius );
	unsigned int wzmx = binning_z( p.z + CurrentKernelRadius );
	real_xyz dSQRKRadius = SQR(CurrentKernelRadius);

	for (unsigned int zz = wzmi; zz <= wzmx; zz++) {
		//unsigned int zoff = zz*spatialdistr.nxy;
		for (unsigned int yy = wymi; yy <= wymx; yy++) {
			//unsigned int yzoff = yy*spatialdistr.nx + zoff;
			for (unsigned int xx = wxmi; xx <= wxmx; xx++) {
				//unsigned int b = xx + yzoff;
				unsigned int b = xx + yy*spatialdistr.nx + zz*spatialdistr.nxy;

				size_t n = pxyz.at(b)->size();
				//###unroll computation SQR(dx)+SQR(dy)+SQR(dz), all have to be tested anyway
				for (size_t i = 0; i < n; ++i) {
					real_xyz dSQR = SQR(pxyz[b]->at(i).x - p.x) + SQR(pxyz[b]->at(i).y - p.y) + SQR(pxyz[b]->at(i).z - p.z);

					if ( dSQR <= dSQRKRadius ) { //inside spherical region 52%, i.e. slightly the more likely case
						results.push_back( dist(sqrt(dSQR), pid[b]->at(i)) );
					}
				} //next point bucket
			} //next xline of cells
		} //next xyslab of cells
	} //done checking intruded cells

	sort( results.begin(), results.end(), SortDistByDistanceAsc);
}


void bvh_xyzm2::destroy_bvh_xyzm2()
{
	//MK::clear datastructure for reinit and reusage
	owner = NULL;

	owin_deformed = aabb3d();
	owin_final = aabb3d();
	spatialdistr = sqb();

	for (size_t mr = 0; mr < pxyz.size(); ++mr) {
		if ( pxyz.at(mr) != NULL ) {
			delete pxyz.at(mr); pxyz.at(mr) = NULL;
		}
	}
	for (size_t mr = 0; mr < pid.size(); ++mr) {
		if ( pid.at(mr) != NULL ) {
			delete pid.at(mr); pid.at(mr) = NULL;
		}
	}
	for (size_t mr = 0; mr < pimage.size(); ++mr) {
		if ( pimage.at(mr) != NULL ) {
			delete pimage.at(mr); pimage.at(mr) = NULL;
		}
	}

	healthy = true;
}



inline aabb3d bvh_xyzm2::get_owindeformed()
{
	return owin_deformed; //#MK::was this->
}



bvh_p3dm1::bvh_p3dm1()
{
	owner = NULL; //MK::owner is only a backreference

	owin = aabb3d();
	spatialdistr = sqb();
	//p remains empty as for now
	healthy = true;
}


bvh_p3dm1::~bvh_p3dm1()
{
	//MK::owner is only a backreference
	for (size_t mr = 0; mr < pp3_points.size(); ++mr) {
		if ( pp3_points.at(mr) != NULL ) {
			delete pp3_points.at(mr); pp3_points.at(mr) = NULL;
		}
	}
}



void bvh_p3dm1::spatialpartitioning_init()
{
cout << "Building bvh_p3dm1 " << endl;
	owin.scale();
cout << "owin " << owin << endl;
	//see comments in bvh_xyzm2, owin is an AABB about all perips and uips

/*
	real_xyz dx = owin.xsz;
	real_xyz dy = owin.ysz;
	real_xyz dz = owin.zsz;

	//minimum dimension of owin
	real_xyz dmin = min(min(dx,dy),dz);
cout << "dmin " << dmin << endl;

	//estimate partitioning, we utilize this BVH space partitioning for the DBScan and the Voxelization
*/

	//MK::bin covering edge length of at least MaxKernelRadius
	real_xyz MaxKernelRadius = max(Settings::DBScanKernelRadius, Settings::PVTessKernelRadius);
cout << "maxkernelr " << MaxKernelRadius << endl;

	//size_t nx = static_cast<unsigned int>( ((dmin / MaxKernelRadius) * dx) + 1.0 + EPSILON);
	//size_t ny = static_cast<unsigned int>( ((dmin / MaxKernelRadius) * dy) + 1.0 + EPSILON);
	//size_t nz = static_cast<unsigned int>( ((dmin / MaxKernelRadius) * dz) + 1.0 + EPSILON);
	size_t nx = static_cast<unsigned int>( (owin.xsz / MaxKernelRadius) + 1.0 + EPSILON);
	size_t ny = static_cast<unsigned int>( (owin.ysz / MaxKernelRadius) + 1.0 + EPSILON);
	size_t nz = static_cast<unsigned int>( (owin.zsz / MaxKernelRadius) + 1.0 + EPSILON);

	spatialdistr = sqb( nx, ny, nz );

cout << "SpatialDistribution " << endl << spatialdistr << endl;

	//#pragma omp parallel
	//{
		//int tid = omp_get_thread_num();
		unsigned int mt = static_cast<unsigned int>(MASTER);

		//##MK::future optimization hint:
		//how it is done now will cause all ips in the fused BVH to end up in threadmemory physically closest to the core
		//executing the master threads, in ccNUMA environment may degrade performance for many threads >> 8 substantially
		//instead partition storage of the final BVH in thread-local memory sections as exemplified for the
		//storage of the DAMASK spectral geometry data in the db database (see this->db)
		bool mehealthy = true;

		//preallocate pointer to buckets holding points in 3d space with a single uint32 mark each (p3dm1)
		pp3_points.reserve(spatialdistr.nxyz);

		for (size_t b = 0; b < spatialdistr.nxyz; ++b) {
			pp3_points.push_back(NULL);
			vector<p3dm2>* pbucket = NULL;
			try {
				pbucket = new vector<p3dm2>;
				pp3_points.back() = pbucket;
			}
			catch (bad_alloc &croak) {
				stopping("Unable to allocate p3d pbucket when building bvh_p3dm1", owner->get_myrank(), mt );
				mehealthy = false; break;
			}
		} //all buckets initalized

		//#pragma omp critical
		//{
			if (mehealthy == false) {
				healthy = false;
			}
		//}
	//} //end of parallel region

	//##MK::implementing feeding back of error path in case of healthy == false
	reporting("Spatialpartitioning initialized", owner->get_myrank(), 0, false);
}


inline unsigned int bvh_p3dm1::binning_x( const real_xyz x )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (x - owin.xmi) / (owin.xmx - owin.xmi) * static_cast<real_xyz>(spatialdistr.nx);
	unsigned int tmpui = ( tmpreal >= 0.0 ) ? static_cast<real_xyz>(tmpreal) : 0;
	unsigned int xx = ( tmpui < spatialdistr.nx ) ? tmpui : spatialdistr.nx-1;
	return xx;
}

inline unsigned int bvh_p3dm1::binning_y( const real_xyz y )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (y - owin.ymi) / (owin.ymx - owin.ymi) * static_cast<real_xyz>(spatialdistr.ny);
	unsigned int tmpui = ( tmpreal >= 0.0 ) ? static_cast<real_xyz>(tmpreal) : 0;
	unsigned int yy = ( tmpui < spatialdistr.ny ) ? tmpui : spatialdistr.ny-1;
	return yy;
}

inline unsigned int bvh_p3dm1::binning_z( const real_xyz z )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	real_xyz tmpreal = (z - owin.zmi) / (owin.zmx - owin.zmi) * static_cast<real_xyz>(spatialdistr.nz);
	unsigned int tmpui = ( tmpreal >= 0.0 ) ? static_cast<real_xyz>(tmpreal) : 0;
	unsigned int zz = ( tmpui < spatialdistr.nz ) ? tmpui : spatialdistr.nz-1;
	return zz;
}

inline unsigned int bvh_p3dm1::binning_xyz( const p3d candidate )
{
	//##MK::requires optimization
	//performing binning with rectangular transfer function
	unsigned int x = binning_x( candidate.x );
	unsigned int y = binning_y( candidate.y );
	unsigned int z = binning_z( candidate.z );

	return (x + y*spatialdistr.nx + z*spatialdistr.nxy);
}


void bvh_p3dm1::build_bvh_p3dm1()
{
	spatialpartitioning_init();

	//ips (perips and uips) are copied over from grain-local sparse bvhs in fill_global_bvh_p3dm1()
}


nbp3d bvh_p3dm1::find_nearest_neighbor( const p3d p, real_xyz r )
{
	//see comments in bvh_p3dm1::find_higherorder_neighbors
	real_xyz CurrentKernelRadius = r;

	unsigned int wxmi = binning_x( p.x - CurrentKernelRadius );
	unsigned int wxmx = binning_x( p.x + CurrentKernelRadius );
	unsigned int wymi = binning_y( p.y - CurrentKernelRadius );
	unsigned int wymx = binning_y( p.y + CurrentKernelRadius );
	unsigned int wzmi = binning_z( p.z - CurrentKernelRadius );
	unsigned int wzmx = binning_z( p.z + CurrentKernelRadius );

	real_xyz dSQRKRadius = SQR(CurrentKernelRadius);
	real_xyz closestSQRDistance = numeric_limits<real_xyz>::max();
	real_xyz closestM1 = numeric_limits<unsigned int>::max();
	p3d closestPoint = p3d();
	unsigned char closestRepYesOrNo = static_cast<unsigned char>(NO);

	//pruning search space by restricting to spatially adjacent candidate buckets
	for (unsigned int zz = wzmi; zz <= wzmx; zz++) {
		for (unsigned int yy = wymi; yy <= wymx; yy++) {
			for (unsigned int xx = wxmi; xx <= wxmx; xx++) {
				//which bucketID to check all points against p?
				unsigned int b = xx + yy*spatialdistr.nx + zz*spatialdistr.nxy;

				vector<p3dm2>* thesepnts = pp3_points.at(b);

				for (size_t i = 0; i < thesepnts->size(); ++i) {

					real_xyz dSQR = SQR( (*thesepnts)[i].x - p.x) + SQR((*thesepnts)[i].y - p.y) + SQR((*thesepnts)[i].z - p.z);

					if ( dSQR <= dSQRKRadius ) { //inside spherical region 52%
						/*
						//however that it is also the closest is very unlikely ##so potential optimization potential
						//if ( dSQR > closestSQRDistance) {
						//	continue;
						//}
						////implicit else
						//closestSQRDistance = dSQR;
						//closestM1 = (*these)[i].m1;
						*/

						//##MK for now better take the safe solution
						if ( dSQR <= dSQRKRadius && dSQR <= closestSQRDistance ) {
							closestSQRDistance = dSQR;
							closestM1 = (*thesepnts)[i].m1;
							closestPoint = p3d( (*thesepnts)[i].x, (*thesepnts)[i].y, (*thesepnts)[i].z );
							closestRepYesOrNo = (*thesepnts)[i].m2;
						}
					}
				} //next point bucket
			} //next xline of cells
		} //next xyslab of cells
	} //done checking intruded cells

	//return dist(sqrt(closestSQRDistance), closestM1);
	return nbp3d( closestPoint.x, closestPoint.y, closestPoint.z, closestM1, closestRepYesOrNo );
}


void bvh_p3dm1::destroy_bvh_p3dm1()
{
	//MK::clear datastructure for reinit and reusage
	owner = NULL;
	owin = aabb3d();

	for (size_t mr = 0; mr < pp3_points.size(); ++mr) {
		if ( pp3_points.at(mr) != NULL ) {
			delete pp3_points.at(mr); pp3_points.at(mr) = NULL;
		}
	}
	pp3_points.clear();

	healthy = true;
}


inline aabb3d bvh_p3dm1::get_owin()
{
	return owin; //#MK::was this->
}


inline sqb bvh_p3dm1::get_spatialdistr()
{
	return spatialdistr; //#MK::was this->
}


grainpool::grainpool()
{
	owner = NULL;
}


grainpool::~grainpool()
{
	//do not delete owner only selfreference
	for(size_t gr = 0; gr < ipsupport.size(); ++gr) {
		if ( ipsupport.at(gr) != NULL ) {
			delete ipsupport.at(gr);
			ipsupport.at(gr) = NULL;
		}
	}
}


void grainpool::build(vector<unsigned int> const & ip2gr)
{
//cout << "Building grainpool..." << endl;
	//find unique grain ids
	set<unsigned int> unique(ip2gr.begin(), ip2gr.end());

	size_t ngrains = unique.size();

cout << "UniqueGrainIDs " << ngrains << endl;

	//##MK::thread parallel grain processing
	//##MK::pragma omp parallel shared(unique)
	//size_t tid = static_cast<size_t>(omp_get_thread_num());

	//allocate memory to store grain-resolved ip support, i.e. organize which ips support which grain
	for(size_t gr = 0; gr < ngrains; ++gr) {
		vector<unsigned int>* theseips = NULL;
		try {
			theseips = new vector<unsigned int>;
			ipsupport.push_back(theseips);
		}
		catch (bad_alloc &croak) {
			string mess = "Allocation error for grains failed!";
			stopping( mess, owner->get_myrank(), 0);
		}
	}

//cout << "IPSupport allocated" << endl;

	//identify ip to grain ownership
	unsigned int neip = static_cast<unsigned int>(ip2gr.size()); //MK::must be head.NXYZ for 117 elements
	for(unsigned int eid = 0; eid < neip; ++eid) { //ips in ip2gr are unique the unique supporting points utilize during DAMASK simulation
		ipsupport.at(ip2gr[eid])->push_back(eid);
	}


	//parallel processing of grain properties now possible
	//#pragma omp parallel
	for(size_t gr = 0; gr < ngrains; ++gr) {
		unsigned int cnts = 0;
		vector<quat> iporis;
		vector<unsigned int>* these = ipsupport.at(gr);
		for(size_t i = 0; i < these->size(); ++i) {
			cnts++;
			iporis.push_back( owner->eid2quaternion(these->at(i)) ); //##MK::memory access inefficient because random...
		}

		quatcloud res = quaternioncloud_characterize( iporis );
		quat qmean = quat( res.qm0, res.qm1, res.qm2, res.qm3 );

		//mean disorientation angle to grain mean --> aka grain reference orientation deviation (GROD)
		polyxx.push_back( grain(qmean, gr, cnts) );
	}
}


void grainpool::report()
{
	complaining( "GrainPool reporting not yet implemented", owner->get_myrank(), 0);
	//##MK::
	//###########
	//###write this->polyxx metadata to file
}


void grainpool::destroy()
{
	//MK::clear datastructure for reinit and reusage
	//do not delete owner only selfreference
	owner = NULL;
	for(size_t gr = 0; gr < ipsupport.size(); ++gr) {
		if ( ipsupport.at(gr) != NULL ) {
			delete ipsupport.at(gr);
			ipsupport.at(gr) = NULL;
		}
	}
	ipsupport.clear();
	polyxx.clear();
}


shape::shape()
{
}


shape::~shape()
{
}


void shape::add( const aabb3d & in, const string str )
{
	//compute corners of axis-aligned bounding box and add to shape library
	/*
	p3d pa = p3d( in.xmi, in.ymi, in.zmi ); //000
	p3d pb = p3d( in.xmi, in.ymi, in.zmx ); //001
	p3d pc = p3d( in.xmi, in.ymx, in.zmi ); //010
	p3d pd = p3d( in.xmi, in.ymx, in.zmx ); //011
	p3d pe = p3d( in.xmx, in.ymi, in.zmi ); //100
	p3d pf = p3d( in.xmx, in.ymi, in.zmx ); //101
	p3d pg = p3d( in.xmx, in.ymx, in.zmi ); //110
	p3d ph = p3d( in.xmx, in.ymx, in.zmx ); //111
	*/

	//##MK:: matlab patch plotting convention
	/*
	 * verts 110, 010, 011, 111, 001,101, 100, 000
	 * faces 1234,4356,6785,1287,6714,2358
	 */

	p3d pg = p3d( in.xmx, in.ymx, in.zmi ); //110
	p3d pc = p3d( in.xmi, in.ymx, in.zmi ); //010
	p3d pd = p3d( in.xmi, in.ymx, in.zmx ); //011
	p3d ph = p3d( in.xmx, in.ymx, in.zmx ); //111
	p3d pb = p3d( in.xmi, in.ymi, in.zmx ); //001
	p3d pf = p3d( in.xmx, in.ymi, in.zmx ); //101
	p3d pe = p3d( in.xmx, in.ymi, in.zmi ); //100
	p3d pa = p3d( in.xmi, in.ymi, in.zmi ); //000

	//sp.push_back( hexahedron( pa, pb, pc, pd, pe, pf, pg, ph) );
	sp.push_back( hexahedron( pg, pc, pd, ph, pb, pf, pe, pa) );
	nm.push_back( str);
}

void shape::add( const bv3x3 & in, const string str )
{
	//add corners of triclinic unit cell with column base vectors in to shape library
	//RVE1
	//to remember convention is 110, 010, 011, 111, 001,101, 100, 000
	real_xyz h, k, l;
	//h = 0.0;	k = 0.0;	l = 0.0;
	h = 1.0;	k = 1.0;	l = 0.0;
	p3d pa = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 1.0;	k = 0.0;	l = 0.0;
	h = 0.0;	k = 1.0;	l = 0.0;
	p3d pb = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 0.0;	k = 1.0;	l = 0.0;
	h = 0.0;	k = 1.0;	l = 1.0;
	p3d pc = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 0.0;	k = 0.0;	l = 1.0;
	h = 1.0;	k = 1.0;	l = 1.0;
	p3d pd = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 1.0;	k = 1.0;	l = 0.0;
	h = 0.0;	k = 0.0;	l = 1.0;
	p3d pe = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 1.0;	k = 0.0;	l = 1.0;
	h = 1.0;	k = 0.0;	l = 1.0;
	p3d pf = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 0.0;	k = 1.0;	l = 1.0;
	h = 1.0;	k = 0.0;	l = 0.0;
	p3d pg = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = 1.0;	k = 1.0;	l = 1.0;
	h = 0.0;	k = 0.0;	l = 0.0;
	p3d ph = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
					0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
					0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	string descr1 = str + "RVE1";
	sp.push_back( hexahedron( pa, pb, pc, pd, pe, pf, pg, ph) );
	nm.push_back( descr1 );

	//RVE27
	//to remember convention is 110, 010, 011, 111, 001,101, 100, 000
	//Moore periodic images of triclinic unit cell  maps 1-->2 and 0-->-1
	//h = +1.0;	k = +1.0;	l = +1.0;
	h = +2.0;	k = +2.0;	l = -1.0;
	pa = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = -1.0;	k = +1.0;	l = +1.0;
	h = -1.0;	k = +2.0;	l = -1.0;
	pb = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = +1.0;	k = -1.0;	l = +1.0;
	h = -1.0;	k = +2.0;	l = +2.0;
	pc = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = +1.0;	k = +1.0;	l = -1.0;
	h = +2.0;	k = +2.0;	l = +2.0;
	pd = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = -1.0;	k = -1.0;	l = +1.0;
	h = -1.0;	k = -1.0;	l = +2.0;
	pe = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = -1.0;	k = +1.0;	l = -1.0;
	h = +2.0;	k = -1.0;	l = +2.0;
	pf = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = +1.0;	k = -1.0;	l = -1.0;
	h = +2.0;	k = -1.0;	l = -1.0;
	pg = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	//h = -1.0;	k = +1.0;	l = +1.0;
	h = -1.0;	k = -1.0;	l = -1.0;
	ph = p3d( 	0.0 + (h*in.a11 + k* in.a12 + l*in.a13) ,
				0.0 + (h*in.a21 + k* in.a22 + l*in.a23) ,
				0.0 + (h*in.a31 + k* in.a32 + l*in.a33) );

	string descr2 = str + "RVE27";
	sp.push_back( hexahedron( pa, pb, pc, pd, pe, pf, pg, ph) );
	nm.push_back( descr2 );
}


cluster::cluster()
{
	lid = numeric_limits<unsigned int>::max();
	cnt = 0;
	aabb3d box = aabb3d();
	chosen = false;
}

cluster::~cluster()
{
	lid = numeric_limits<unsigned int>::max();
	cnt = 0;
	aabb3d box = aabb3d();
	chosen = false;
}


grGeomHdl::grGeomHdl()
{
	owner = NULL;
	cgid = numeric_limits<unsigned int>::max();
	nuip = 0;
	localbox = aabb3d();
	//let container for now empty
	grainfence = aabb3d();
	localgrid = vxlgrd();
	nuvxl = 0;

	//check through which tasks have been accomplished
	alloc_success = true;
	ips_success = true;
	sbvh_success = true;
	dbscan_success = true;
	dbpick_success = true;
	vxlinit_success = true;
	vxlfill_success = true;
	sdfinit_success = true;
	sdfsprd_success = true;

	healthy = true;
}


grGeomHdl::~grGeomHdl()
{
	//do not delete owner is only a self-reference
	cgid = numeric_limits<unsigned int>::max();
	nuip = 0;
	localbox = aabb3d();
	nuvxl = 0;

	ips.clear();

	sbvh_locator.clear();

	sbvh_points.clear();
	sbvh_uipref.clear();
	sbvh_isrepresentative.clear();

	lbl.clear();
	dbscanres.clear();
	thegrain.clear();
	grainfence = aabb3d();
	localgrid = vxlgrd();

	GrainIDField.clear();
	UIPIDField.clear();
	IsRepresentativeField.clear();
	SgnDistField1.clear();
	SgnDistField2.clear();

	alloc_success = true;
	ips_success = true;
	sbvh_success = true;
	dbscan_success = true;
	dbpick_success = true;
	vxlinit_success = true;
	vxlfill_success = true;
	sdfinit_success = true;
	sdfsprd_success = true;

	healthy = true;
}


bool grGeomHdl::debug_duplicate_check( const real_xyz x, const real_xyz y, const real_xyz z )
{
	real_xyz tol = EPSILON;
	for(size_t i = 0; i < ips.size(); ++i) { //#MK::was this->ips
		real_xyz val = SQR(ips[i].x-x) + SQR(ips[i].y-y) + SQR(ips[i].z-z);
		if ( val >= tol ) { //most likely case --> no duplicate
			continue;
		}
		else
			return true;
	}
	return false;
}


void grGeomHdl::compute_all_replica()
{
	//MK::EXECUTED FROM WITHIN OMP PARALLEL REGION!
	//scan all 27 neighboring periodic replica of the triclinic deformed point cloud considering local uip displacements
	//for all uips assigned to grain
	size_t uips_identified = ips.size();
	size_t uips_expected = static_cast<size_t>(nuip); //#MK::was this->nuip

cout << "Unique ips ident/expected/expected diff accstrategy " << uips_identified << "\t\t" << uips_expected << "\t\t" << static_cast<size_t>(owner->grains.polyxx.at(cgid).np) << endl;

	//buffer replica p3d with reference to uniqueip generator as mark m1 in ips before spatially reorganizing
	//ips into sparse bounded volume hierarchy (BVH) according to thread-synchronized globalbox
	bv3x3 per = owner->rveBaseDef; //RVE base vector in deformed configuration
	aabb3d box = aabb3d();

	//identification scheme for the Moore periodic images is, 0x00 is the center, 1 is -1-1-1, 2 0-1-1, and so forth
	real_xyz xxx = 0.0;
	real_xyz yyy = 0.0;
	real_xyz zzz = 0.0;
	real_xyz h = 0.0;
	real_xyz k = 0.0;
	real_xyz l = 0.0;

	for (size_t u = 0; u < uips_identified; ++u) { //MK::an index on this->ips[uip] valid to access only for [0,uips_identified)
		real_xyz xx = ips[u].x;
		real_xyz yy = ips[u].y;
		real_xyz zz = ips[u].z;
		unsigned int mark = ips[u].m1; //not the loop index variable u, we need the name of the original unique ip!

		//<000>
		//[000], 0x00
		h = +0.0;		k = +0.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		//unfortunately does not apply for most but runtime is scales at least O(n)
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
		//in-place appending
		//MK::do not add the image (0,0,0) as it is the uip itself and has already been added
		//<100>
		//[+1+0+0], 0x0e
		h = +1.0;		k = +0.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+0+0], 0x0d
		h = -1.0;		k = +0.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0+1+0], 0x10
		h = +0.0;		k = +1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0-1+0], 0x0b
		h = +0.0;		k = -1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0+0+1], 0x16
		h = +0.0;		k = +0.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0+0-1], 0x05
		h = +0.0;		k = +0.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//<110>
		//[+1+1+0], 0x11
		h = +1.0;		k = +1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+1+0], 0x0f
		h = -1.0;		k = +1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1-1+0], 0x0c
		h = +1.0;		k = -1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1-1+0], 0x0a
		h = -1.0;		k = -1.0;		l = +0.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0+1+1], 0x19
		h = +0.0;		k = +1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0-1+1], 0x13
		h = +0.0;		k = -1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0+1-1], 0x08
		h = +0.0;		k = +1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+0-1-1], 0x02
		h = +0.0;		k = -1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1+0+1], 0x17
		h = +1.0;		k = +0.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+0+1], 0x15
		h = -1.0;		k = +0.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1+0-1], 0x06
		h = +1.0;		k = +0.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+0-1], 0x04
		h = -1.0;		k = +0.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//<111>
		//[+1+1+1], 0x1a
		h = +1.0;		k = +1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+1+1], 0x18
		h = -1.0;		k = +1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1-1+1], 0x14
		h = +1.0;		k = -1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1+1-1], 0x09
		h = +1.0;		k = +1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[+1-1-1], 0x03
		h = +1.0;		k = -1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1+1-1], 0x07
		h = -1.0;		k = +1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1-1+1], 0x12
		h = -1.0;		k = -1.0;		l = +1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
		//[-1-1-1], 0x01
		h = -1.0;		k = -1.0;		l = -1.0;
		xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
		yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
		zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
		box.potentially_expand(xxx, yyy, zzz);
		//if (xxx <= box.xmi) box.xmi = xxx; 	if (xxx >= box.xmx)	box.xmx = xxx;	if (yyy <= box.ymi) box.ymi = yyy;	if (yyy >= box.ymx)	box.ymx = yyy;	if (zzz <= box.zmi) box.zmi = zzz;	if (zzz >= box.zmx)	box.zmx = zzz;
//if ( debug_duplicate_check(xxx,yyy,zzz) == true ) { cout << "ADDING DUPLICATE POINT FOR GRAIN " << cgid << endl; }
		ips.push_back(p3dm1(xxx, yyy, zzz, mark));
	} //all these perips support periodic replica fragments of the periodic images of grain gr

	localbox = box;
	localbox.scale();


	//##MK::DEBUG if ( Settings::VisPeriodicReplicaIPs == 0 || Settings::VisPeriodicReplicaIPs == (cgid + 1) ) {

	//##MK::thread writing file with different name than any other thread (disjoint grainIDs) therefore can be done in parallel, but of course may not be efficient when only one harddrive
	if ( Settings::VisPeriodicReplicaIPs == (cgid+1) ) { //##MK::DEBUG overwrite
		string different_than_otherthreads = "RVE27IPsGrainID" + to_string(cgid);
		bool status = vtk_p3dm1( ips, cgid, owner->thisincrement, different_than_otherthreads );
	}

	#pragma omp critical
	{
cout << "Total ips for grain " << cgid << " are " << ips.size() << endl;
	}
}


/*
void grGeomHdl::update_threadlocalbox( aabb3d & current_extremal)
{
	//MK::EXECUTED FROM WITHIN OMP PARALLEL REGION!
	if ( localbox.xmi <= current_extremal.xmi )	current_extremal.xmi = localbox.xmi;
	if ( localbox.xmx >= current_extremal.xmx )	current_extremal.xmx = localbox.xmx;
	if ( localbox.ymi <= current_extremal.ymi )	current_extremal.ymi = localbox.ymi;
	if ( localbox.ymx >= current_extremal.ymx )	current_extremal.ymx = localbox.ymx;
	if ( localbox.zmi <= current_extremal.zmi )	current_extremal.zmi = localbox.zmi;
	if ( localbox.zmx >= current_extremal.zmx )	current_extremal.zmx = localbox.zmx;
}
*/


void grGeomHdl::build_sparse_bvh()
{
	//MK::EXECUTED FROM WITHIN OMP PARALLEL REGION!

	//ips contains points which are unsorted and spatially unstructured but sampling space sparsely and inhomogeneous
	//(apart from the fact that the first [0,nuip) points are unique)
	//this is because these points sample/support the periodic images and fragments of a grain in RVE27,
	//clustering therefore obviously locally dense, showing no noise, causing most bins of the global BVH empty in ips with
	//association to grain cgid
	//except for a few, hence to improve memory locality in the subsequent grain-local DBScan algorithm
	//we reorganize the storage of the ips as follows:

	//first we identify the target bin for each ip in ips and fill the target container with dummies

	//##MK::what happens if no further memory allocatable?

	vector<unsigned int> i2bin;
	for (size_t i = 0; i < ips.size(); ++i ) {
		p3d p = p3d( ips.at(i).x, ips.at(i).y, ips.at(i).z );
		unsigned int b = owner->pp3rve27.binning_xyz(p);
		i2bin.push_back( b );
//cout << "gr/ipid/p/b\t\t" << cgid << "\t\t" << i << "\t\t" << b << "\t\t" << p << endl;
		sbvh_points.push_back(p3d());
		sbvh_uipref.push_back(numeric_limits<unsigned int>::max());
		sbvh_isrepresentative.push_back(static_cast<unsigned char>(NO));
	}

	//now we identify which binIDs have at all points and how many
	map<unsigned int, size_t> cnts;
	for(size_t i = 0; i < i2bin.size(); ++i ) {
		unsigned int key = i2bin.at(i);
		if (cnts.find(key) != cnts.end()) //key exists, this is most likely is that after a few loop iterations a key exists
			cnts[key]++;
		else
			cnts.insert(make_pair(key, 1));
	}
	//MK::the C++11 standard + the by construction disjointness of keys as binIDs guarantees
	//that the elements in cnts are strictly ordered ascendingly

	//now we build the indexing meta structure detailing the organization of the vector sbvh as follows:
	//blocks of p3d all belonging to one bin are glued together, with the bin of smallest binID first,
	//binIDs are not consecutive hence we need an indexing structure
	//that tells us where we find the p3d's in that bin

	//##MK::future optimization may eradicate some these redundant but safe maps
	//which for now aiom first to assure a correctly working algorithm
	map<unsigned int, size_t> cnts_start;
	map<unsigned int, size_t> cnts_end; //cumulate cnts
	map<unsigned int, size_t> cnts_trsfd; //how many already transferred?

	size_t startidx = 0;
	size_t endidx = 0;
	map<unsigned int, size_t>::iterator it = cnts.begin();
//cout << "sparseBinCnts on grain " << cgid << endl;
	while( it != cnts.end()) { //traversing along indexing structure by standard and see above comment sweeps ascendingly
//cout << "binID/cnts " << it->first << "\t\t" << it->second << endl;
		cnts_start.insert(make_pair(it->first, startidx));
		endidx = startidx + it->second;
		cnts_end.insert(make_pair(it->first, endidx));
		startidx = startidx + it->second;
		cnts_trsfd.insert(make_pair(it->first, 0));
		it++;
	}

	it = cnts.begin();
	while (it != cnts.end()) {
		unsigned int key = it->first;
		sbvhrange tmp = sbvhrange( cnts_start[key], cnts_end[key]);
		sbvh_locator.insert(make_pair(key, tmp));
//cout << "binID/idxs/idxe " << key << "\t\t" << sbvh_locator.find(key)->second.startidx << "\t\t" << sbvh_locator.find(key)->second.pastendidx << endl;
 		it++;
	}

//##MK::BEGIN DEBUG
	map<unsigned int, sbvhrange>::iterator jt = sbvh_locator.begin();
	while (jt != sbvh_locator.end()) {
//cout << "binID/Range[/) " << jt->first << "\t\t" << jt->second.startidx << ";" << jt->second.pastendidx << endl;
		jt++;
	}
//##MK::END DEBUG

	//fill sbvh with ips in ordered manner according to known binning
	for(size_t i = 0; i < i2bin.size(); ++i ) {
		unsigned int key = i2bin.at(i);
		//##MK::add checks if key exists...
		size_t pos = sbvh_locator[key].startidx + cnts_trsfd[key];
		sbvh_points[pos] = p3d(ips[i].x, ips[i].y, ips[i].z); //transfer point coordinate and original uniqueip mark
		sbvh_uipref[pos] = ips[i].m1;
		//we dont know yet which points of sbvh_points are in the representative cluster or not
		cnts_trsfd[key]++;
	}

#pragma omp critical
{
	cout << "SparseBinningConstructed for grain " << cgid << endl;
}

	//##MK::now at least technically ips is obsolete ips.clear();

	//MK::ips no longer needed clear them and replace the associated potentially MB memory with single page fresh vector by swapping trick, MK::because std::vector<T>::clear only deletes content but does not give allocated pages back to OS
	vector<p3dm1> dummy;
	ips.swap( dummy );
}


void grGeomHdl::spatial_range_query_sbvh_noclear_nosort(const size_t thisip, vector<unsigned int> & results, real_xyz r )
{
	//returns array indices on sbvh_points, in every case zero element is always myself i.e. ip all others not sorted by distance just inclusion in ball with radius r about point sbvh_points.at(ip)

	//does bin contain data at all?, MK::one might wonder at this point why to use at all a sbvh as the function
	//is called only in the DBScan where anyways only ips from one particular grain cluster are processed,
	//reason is that if one would store the bvh as vector<vector<...>*> construct the individual bins would
	//point not necessarily to as closely group memory pages/locations
	//of course one could quickly check if data in the bucket exists by construct.at(firstlevel) ? NULL
	//but then not only storing potentially many NULL pointer but once switching the bin one may jump farther in memory...
	//that though might be counteracted by increased control logic to process the sbvh positions, i dunno which approach is best...
	//but 500gr*100000 buckets per grain*8B per pointer is already 380MB only --- mostly --- NULL pointer of which
	//approximately only 1/500 of the total integration points will end up in the grains own sbvh...
	//instead with the sparsebvh here everything is not guaranteed (because thats the allocators choice)
	//but to be expected much better memory-aligned on sbvh thus improving data locality during processing and reducing page switches
	//a task --- spatial range query within r-ball about point location p --- is usually memorybound

	//idx, indices nbor_eip refer to positional information on sbvh[idx] on [0,size())

	//##MK::implement checks a la //if ( r < CurrentKernelRadius && r > 0.0)
	real_xyz R = r;

	//MK::performs no clear on input!!

	//take care that bins do not protrude out of the owner->pp3rve27 BVH
	sqb gspatdistr = owner->pp3rve27.get_spatialdistr();

	p3d p = sbvh_points.at(thisip);
	results.push_back(static_cast<unsigned int>(thisip));

	//identify which cuboidal region of bins to check at all
	unsigned int wxmi = owner->pp3rve27.binning_x(p.x-R); //imi
	unsigned int wxmx = owner->pp3rve27.binning_x(p.x+R); //imx
	unsigned int wymi = owner->pp3rve27.binning_y(p.y-R);
	unsigned int wymx = owner->pp3rve27.binning_y(p.y+R);
	unsigned int wzmi = owner->pp3rve27.binning_z(p.z-R);
	unsigned int wzmx = owner->pp3rve27.binning_z(p.z+R);

	//perform distance computations and checking who is in r-ball about p
	real_xyz dSQRKRadius = SQR(R);
	map<unsigned int, sbvhrange>::iterator it;
	for (unsigned int zz = wzmi; zz <= wzmx; zz++) {
		for (unsigned int yy = wymi; yy <= wymx; yy++) {
			for (unsigned int xx = wxmi; xx <= wxmx; xx++) {
				unsigned int b = xx + yy*gspatdistr.nx + zz*gspatdistr.nxy;

				//check if bin in the sparse bvh exists at all,
				//remember only bins with at least one point exist with heterogeneous spatial distribution
				it = sbvh_locator.find(b);
				if ( it != sbvh_locator.end() ) { //the bucket contains at all points
					size_t idxs = it->second.startidx;
					size_t idxp = it->second.pastendidx;

					for (size_t i = idxs; i < idxp; i++) { //MK::here we see the benefit of having sbvh_points and sbvh_uipref split, ideally the memory section of sbvh_points is shorter than when having the labels thus, quicker to load + more likely frequently reutilized
						real_xyz dSQR = SQR(sbvh_points[i].x-p.x) + SQR(sbvh_points[i].y-p.y) + SQR(sbvh_points[i].z-p.z);

						if ( dSQR <= dSQRKRadius ) { //inside spherical region 52%, i.e. slightly the more likely case
							if ( i != thisip ) { //most cases
								results.push_back( static_cast<unsigned int>(i) ); //MK::i is an index on sbvh! ''MK::should be size_t instead of unsigned int
							}
						}
					} //next (very likely) cache-aligned element in the bucket
				} //done checking an existent bucket
			} //next xline of cells
		} //next xyslab of cells
	} //done checking intruded cells

	//we do not sort saving lgN(N)
}


void grGeomHdl::dbscan_ips2grreplicates( const real_xyz eps )
{
	//CALLED FROM WITHIN PARALLEL REGION

	//performing sequential DBScan at thread scale with thread/grain-class-local data structure
	//thereby avoiding concurrency, critical sections, and reducing (if placement of grGeomHdl is in distant
	//thread-local memory) false sharing
	//other threads meanwhile do DBScanning on other grains also independently

	//ips are the periodic replications + the potentially fragmented pieces of the grain cgid which
	//in the original DAMASK simulation made contact with the boundary
	//MK::lets call these fragments and replications of a grain "replicates"
	//if the grain had no RVE domain boundary contact in total 27 replicates should exists in the Moore replication of the RVE1 uips
	//if the grain lays in a corner of the RVE domain 8 replicates exists and so forth...
	//never the less each replicated cluster ought to have the same number of supporting integrations as detected
	//for the grain in the hierarchical community detection algorithm BUT
	//MK::WARNING if grains PERIODICALLY WRAP AROUND THE ORIGINAL RVE DOMAIN (singlecrystalline, layer-stacks, etc...)
	//these ought to be multiplied integer way in the RVE27, hence their number of ips an uint multiple of nuip

	//anyhow we should now find the corresponding ips that sample the position of the replicates in 3D space
	//by performing a clustering analysis. As we do not know in advance how many replicates to expect
	//we cannot, though, apply a k-means-type clustering algorithm, instead we use DBScan...
	//... M. Ester, H.-P. Kriegel, J. Sander, X. Xu: A density-based algorithm for discovering clusters in large
	//spatial databases with noise, Proc. 2nd Int. Conf. Knowledge Discovery and Data Mining KDD-96, 1996

	//original DBScan distinguishes noise, boundary points and core points. However, we are not interested in noise,
	//there ought to be in fact no noise because all ips per grain belong to periodic fragments of that grain with ID cgid
	//also we are not interested in whether an ip bounds a cluster

	//we check so by defining a minimum ip distance eps to discern cluster
	//we do not have to worry about grains being adjacent, thanks to the partiontioning of the problem at the grain
	//scale the threadlocal ip sbvh does not contain any ips sampling/supporting adjacent grains in the polycrystal
	//as as the second DBScan parameter we assume it suffices for the cluster to find at least one neighboring
	//point in the eps, i.e. Settings::DBScanKernelRadius-ball

	//string mess = "Implementation currently uses uint32 as labels
	//thus be careful when working with >512^3 DAMASK domains";
	if ( owner->head.NXYZ > static_cast<unsigned int>(300*300*300)) {
		//complaining(mess, owner->get_myrank(), 0);
		stopping("Implementation utilizes 32-bit labels may not work for so large domains", owner->get_myrank(), omp_get_thread_num());
		dbscan_success = false; healthy = false; //#MK::was this->
		return;
	}
	//DBscan implementation:: see for instance an implementation in https://github.lcom/propanoid/DBSCAN/blob/master/dbscan.cpp
	//strictly speaking DBScan is formulated for a mean-field density instead here we utilize so to say space with Dirac spikes at the ips

	//prepare cluster labels
	vector<dbcontrol> lookup;
	for(size_t i = 0; i < sbvh_points.size(); ++i) {
		lookup.push_back( dbcontrol() );
	}

	/*
	//size_t nips = sbvh_points.size();unsigned int* labels = NULL;
	try { labels = new unsigned int[nips]; }
	catch (bad_alloc &mecroak) {
		stopping("Allocation error for label init in DBScan", owner->get_myrank(), omp_get_thread_num());
		this->dbscan_success = false; this->healthy = false;
		return;
	}
	//prepare visited bitmap, ##MK::ideally somewhen utilize also a true bitmap and not a bool array!
	bool* visited = NULL;
	try { visited = new bool[nips]; }
	catch (bad_alloc &mecroak) {
		stopping("Allocation error for visited in DBScan", owner->get_myrank(), omp_get_thread_num());
		this->dbscan_success = false; this->healthy = false;
		return;
	}
	for(size_t ip = 0; ip < nips; ++ip) {
		labels[ip] = numeric_limits<unsigned int>::max(); //flag unassigned
		visited[ip] = false;
	}
	//unsigned int clusterid = 0;
	*/

	vector<unsigned int> ipres;
	vector<unsigned int> nboripres;
	unsigned short clusterid = 0;

	for(unsigned int ip = 0; ip < sbvh_points.size(); ++ip) {
		if ( lookup.at(ip).visited == false ) {
			lookup.at(ip).visited = true;
			//p3d thisposition = sbvh_points.at(ip);

			ipres.clear();
			spatial_range_query_sbvh_noclear_nosort( ip, ipres, eps );

			//MK::given the specific problem of compact well clustered objects, there needs to be at least myself and one additional reference to an ip in ipres
			if ( ipres.size() > 1 ) { //##MK:: was 1, find at least myself, first entry is ip itself but we want at least one nearest neighbor in the DBScanKernelRadius ball, noise or grains supported by only one ip are considered noise
				//if there is at least another point in eps about me that is not myself label me
				lookup.at(ip).label = clusterid;

				//skipping myself check now successively all neighbors whether they have neighbors
				for(unsigned int j = 1; j < ipres.size(); ++j) { //skipping myself, therefore start at 1
					unsigned int nborip = ipres.at(j); //##MK::was .nbor_eipid; //MK::index of an ip according to order in sbvh

					if (lookup.at(nborip).visited == false) {
						lookup.at(nborip).visited = true;
						//p3d nborposition = sbvh_points.at(nborip);

						nboripres.clear();
						spatial_range_query_sbvh_noclear_nosort( nborip, nboripres, eps);

						//MK::nboripres will now contain nborip but ipres already contains it hence ignore
						if ( nboripres.size() > 1 ) { //only pushback those references to ipres which are not nborip
							for(size_t k = 1; k < nboripres.size(); ++k) {
								ipres.push_back( nboripres.at(k) );
							} //having added all close bys within eps about nborposition thereby climbing ape-like through space
						} //done adding if there were more than myself
					}

					if (lookup.at(nborip).label == numeric_limits<unsigned short>::max()) {
						lookup.at(nborip).label = clusterid;
					}
				} //do unless no further neigbor relations connected by distance <= eps are found

				++clusterid;
			}
		} //done processing unvisited ip
	} //proceed with next ip

	//transfer labeling
	for(size_t i = 0; i < sbvh_points.size(); ++i) {
		lbl.push_back( static_cast<unsigned int>(lookup.at(i).label) );
	}

	//##MK::thread writing file with different name than any other thread (disjoint grainIDs) therefore can be done in parallel, but of course may not be efficient when only one harddrive
	if ( Settings::VisDBScanClusterIDs == (cgid+1) ) {
		//if ( Settings::VisDBScanClusterIDs == 0 ) { //##MK::DEBUG overwritten
		string different_than_otherthreads = "RVE27DBScanGrainID" + to_string(cgid);
		bool status = vtk_p3dm1( sbvh_points, lbl, cgid, owner->thisincrement, different_than_otherthreads);
	}

#pragma omp critical
{
	cout << "DBScanPerformed for grain " << cgid << endl;
}

	//cleanup
	/*if ( visited != NULL ) { delete [] visited; visited = NULL; }
	if ( labels != NULL ) { delete [] labels; labels = NULL; }*/
}


void grGeomHdl::dbscan_pick_representative()
{
	//MK::EXECUTED FROM WITHIN OMP PARALLEL REGION!

	//how many ips per cluster?
	dbscanres.clear();
	map<unsigned int, unsigned int> relabeling;
	unsigned int newlabel = 0;
	//##MK::potentially build collection of unique clusterIDs by std::set
	//##MK::strictly speaking as the DBScan is executed locally the relabeling should not even be necessary if indexing started at newlabel
	for(size_t ip = 0; ip < lbl.size(); ++ip) {
		unsigned int oldlabel = lbl[ip];
		if (relabeling.find(oldlabel) != relabeling.end())
			continue;
		//implicit else
		relabeling.insert(make_pair(oldlabel, newlabel));
		newlabel++;
	}

	//first initialize cluster object instances
	map<unsigned int, unsigned int>::iterator it = relabeling.begin();
	while( it != relabeling.end() ) {
		dbscanres.push_back( cluster() );
		dbscanres.back().lid = it->second; //the new label of this cluster
		dbscanres.back().cnt = 0;
		it++;
	}

	//MK::MIND THAT THE labels on this->lbl are NOT REPLACED (yet)!

	//fence axis/aligned bounding boxes about the identified cluster from the DBScan
	for(size_t ip = 0; ip < sbvh_points.size(); ++ip) {
		unsigned int oldlabel = lbl[ip];
		unsigned int cid = relabeling[oldlabel];
		dbscanres[cid].cnt++;
		dbscanres[cid].box.potentially_expand( sbvh_points[ip].x, sbvh_points[ip].y, sbvh_points[ip].z );
	}

	//find center of uip point cloud and number of uip of the grain this->cgid from the Louvain labeling algorithm
	//aabb3d about uip deformed configuration with no guard zone
	p3d rvecenter = owner->pp3rve1withguard.get_owindeformed().center();

	//how many uips with the same Louvain labeling algorithm ID
	unsigned int nip_required = nuip;

cout << "Pick representative of grain " << cgid << endl;
cout << "nip_required/rvecenter " << nip_required << "\t\t" << rvecenter;
	//what is the cluster closest to rvedomain center and is it complete in ips as expected?
	unsigned int bestvalidcid = numeric_limits<unsigned int>::max();
	real_xyz bestdistance = numeric_limits<real_xyz>::max();

	//find best candidate cluster by checking all, number of cluster typically a few dozens
	for(size_t cid = 0; cid < dbscanres.size(); ++cid) {
		dbscanres[cid].box.scale();
		p3d cp = dbscanres[cid].box.center();

		real_xyz dist = sqrt(SQR(cp.x-rvecenter.x)+SQR(cp.y-rvecenter.y)+SQR(cp.z-rvecenter.z));

/*
		//##MK::BEGIN DEBUG should be in critical section
cout << "ClusterID/Center/Distance " << cid << "\t\t" << cp.x << ";" << cp.y << ";" << cp.z << "\t\t" << dist << endl;
cout << "ClusterID/AABB " << cid << dbscanres[cid].box << endl;
		//##MK::END DEBUG
*/

		if ( dist <= bestdistance && dbscanres[cid].cnt == nip_required ) {
			bestvalidcid = cid;
			bestdistance = dist;
		}
	}

	//most importantly: have we found such a cluster, meaning was the reconstruction, i.e. picking of one complete grain successful ?
	if (bestvalidcid == numeric_limits<unsigned int>::max()) {
cout << "Picking representative for grain " << cgid << " no cluster found!" << endl; //##MK::DEBUG use std:verbose functions
		dbpick_success = false; healthy = false; //#MK::was this->
		return;
	}

	//obviously something was found well then copy over the ips
	dbscanres[bestvalidcid].chosen = true;

	grainfence = dbscanres[bestvalidcid].box;
	grainfence.scale();
cout << "Grain fence for grain " << cgid << " is " << endl << grainfence << endl;
	//filter out from all ips those sampling this replicate, i.e. thegrain, positions may live in RVE27
	//for all grains during the DAMASK simulation making no domain boundary contact ought to be in RVE1 because closest to center was chosen
	//and bounding boxes frame the entire RVE27
	unsigned int targetlabel_of_thegrain = bestvalidcid;
	for(size_t ip = 0; ip < sbvh_points.size(); ++ip) {
		if (sbvh_points[ip].x < grainfence.xmi) continue; //most ips will belong to other cluster, so make use of fact that most often we will early and likely reject
		if (sbvh_points[ip].x > grainfence.xmx) continue;
		if (sbvh_points[ip].y < grainfence.ymi) continue;
		if (sbvh_points[ip].y > grainfence.ymx) continue;
		if (sbvh_points[ip].z < grainfence.zmi) continue;
		if (sbvh_points[ip].z > grainfence.zmx) continue;
		//but do not continued if ip is in fact inside the grain fence
		//but given the possibility of nonconvex shapes the AABB may still protrude in space occupied by other cluster

		//so check in addition for == with the correct cluster label
		unsigned int oldlabel = lbl[ip];
		it = relabeling.find(oldlabel);
		if ( it != relabeling.end() ) {
			if ( it->second == targetlabel_of_thegrain ) {
				thegrain.push_back( sbvh_points[ip]);
				sbvh_isrepresentative[ip] = static_cast<unsigned char>(YES);
				//MK::why is storing such boolean information necessary? because during voxelizing inside RVE27 for small RVEs periodic wraparound of
				//grains may cause that inside the local voxelgrid more voxel than building the representatively chosen cluster have the grain ID as the
				//reconstructed grain, i.e. case of lurking in periodic images of grains within localgrid
				//hence in order to pass within the voxelization process only all ips which really build the representative cluster and not any other
				//fragment of other periodic images of the same grain i.e. different cluster, but those not chosen as the representative, we store here
				//which ips in RVE local for each grain build the cluster to later filter out specifically not only for each voxel what is the closest
				//ip hence uip hence grainID but also whether this ip builds the representative cluster of cgid
			}
			else { //##MK::be safe but should really not be encountered!
				string mess = "Labeling inconsistencies during identification of grain " + to_string(cgid);
				complaining(mess, owner->get_myrank(), omp_get_thread_num() );
				dbpick_success = false; healthy = false; //#MK::was this->
				return;
			}
		} //done with checking ip assignment
	}


	//##MK::thread writing file with different name than any other thread (disjoint grainIDs) therefore can be done in parallel, but of course may not be efficient when only one harddrive
	if ( Settings::VisReconstructedGrain == (cgid+1) ) {
		//if( Settings::VisReconstructedGrain == 0 ) { //##MK::DEBUG overwrite
		string different_than_otherthread = "ReconstructedGrain" + to_string(cgid);
		bool status = vtk_p3d( thegrain, cgid, owner->thisincrement, different_than_otherthread);
	}

#pragma omp critical
{
	cout << "Picking representative for grain " << cgid << " clusterID/distance " << dbscanres[bestvalidcid].lid << "/" << bestdistance;
	cout << " clusterCnts " << dbscanres[bestvalidcid].cnt << " should be " << nuip << endl;
}

}


void grGeomHdl::build_localgrid()
{
	//EXECUTED FROM WITHIN PARALLEL REGION

	//identify at which locations the localgrid is embedded in the global one
	//MK::this->grainfence already includes blowup
	//is it numerically sound
	size_t lxmi = owner->thegrid.binning_x( grainfence.xmi );
	if ( lxmi == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; } ////#MK::was this->,   gridding was unsuccessful
	size_t lxmx = owner->thegrid.binning_x( grainfence.xmx );
	if ( lxmx == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; }

	size_t lymi = owner->thegrid.binning_y( grainfence.ymi );
	if ( lymi == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; }
	size_t lymx = owner->thegrid.binning_y( grainfence.ymx );
	if ( lymx == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; }

	size_t lzmi = owner->thegrid.binning_z( grainfence.zmi );
	if ( lzmi == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; }
	size_t lzmx = owner->thegrid.binning_z( grainfence.zmx );
	if ( lzmx == numeric_limits<size_t>::max() ) { vxlinit_success = false; healthy = false; return; }

/*
cout << "Grainfence grain " << cgid << endl << grainfence << endl;
cout << "lxyzmi = " << lxmi << ";" << lymi << ";" << lzmi << endl;
cout << "lxyzmx = " << lxmx << ";" << lymx << ";" << lzmx << endl;
*/

	localgrid.box = grainfence;
	localgrid.dcell = Settings::PVTessCubeVoxelEdgeLen;
	localgrid.origin_cntum = p3d( grainfence.xmi, grainfence.ymi, grainfence.zmi ); //##MK:: not necessarily exactly boundary but lower domain wall of voxel not center of gravity of voxel
	localgrid.origin_discr_own = vxl(0, 0, 0);
	localgrid.origin_discr_lnk = vxl( lxmi, lymi, lzmi ); //my origin with respect to the linked global thegrid

	localgrid.nx = lxmx - lxmi;
	localgrid.ny = lymx - lymi;
	localgrid.nz = lzmx - lzmi;
	localgrid.nxy = (lxmx - lxmi) * (lymx - lymi);
	localgrid.nxyz = (lxmx - lxmi) * (lymx - lymi) * (lzmx - lzmi);

cout << "Final grain-local grid for " << cgid << endl; // << " is " << localgrid << endl;

	//##MK::BEGIN DEBUG dummy GrainID field
	try {
		GrainIDField.resize( localgrid.nxyz, numeric_limits<unsigned int>::max() );
	}
	catch (exception &me_nouipid) {
cout << "Grain " << cgid << " allocation error for GrainIDField" << endl; //##MK::DEBUG should be in OMP critical
		vxlinit_success = false; healthy = false;
		return;
	}
	//##MK::END OF DEBUG

	//fill UIPIDField with dummies to assure
		//that memory is allocated when attempting storage of closest uipid during voxelization
		//that error management can prevent onthefly running out of memory
		//that the dummy value argument has already taken care of the case that there is no close uipid to the vxl from the RVE27 realization
		//	in case where the grid is so heavily distorted
	try {
		UIPIDField.resize( localgrid.nxyz, numeric_limits<unsigned int>::max() );
	}
	catch (exception &me_nouipid) {
cout << "Grain " << cgid << " allocation error for UIPIDField" << endl; //##MK::DEBUG should be in OMP critical
		vxlinit_success = false; healthy = false;
		return;
	}

	//fill IsRepresentativeField with dummies indicating that so far the voxel is not assigned to an ip included in the point cloud of the representative chosen dbscan cluster
	try {
		IsRepresentativeField.resize( localgrid.nxyz, static_cast<unsigned char>(NO) );
	}
	catch (exception &me_nouipid) {
cout << "Grain " << cgid << " allocation error for IsRepresentativeField" << endl; //##MK::DEBUG should be in OMP critical
		vxlinit_success = false; healthy = false;
		return;
	}

	real_sdf h = static_cast<real_sdf>(localgrid.dcell);
	real_sdf h_outside = static_cast<real_sdf>(-1.0) * h;
	try {
		SgnDistField1.resize( localgrid.nxyz, h_outside );
	}
	catch (exception &me_nosdf) {
cout << "Grain " << cgid << " allocation error for SgnDistField1" << endl; //##MK::DEBUG should be in OMP critical
		vxlinit_success = false; healthy = false;
		return;
	}

	real_sdf initial = static_cast<real_sdf>(FASTSWEEPING_INITIAL_VALUE);
	try {
		SgnDistField2.resize( localgrid.nxyz, initial );
	}
	catch (exception &me_nosdf) {
cout << "Grain " << cgid << " allocation error for SgnDistField2" << endl; //##MK::DEBUG should be in OMP critical
		vxlinit_success = false; healthy = false;
		return;
	}

	//##MK::thread writing file with different name than any other thread (disjoint grainIDs) therefore can be done in parallel, but of course may not be efficient when only one harddrive
	if ( Settings::VisGrainLocalVoxelGrid == (cgid+1) ) { //##MK::DEBUG overwritten
		string different_than_otherthread = "GrLocalInitialUIPVxlGridGrainID";
		string descr = "Grain local initial voxelgrid";
		bool status = vtk_vxlgrdm2( owner->thegrid, localgrid, UIPIDField, SgnDistField1,
				cgid, owner->thisincrement,
					different_than_otherthread, descr, "UIPIDMark", "SDFValue" );
	}

#pragma omp critical
{
cout << "Local voxelgrid constructed for grain " << cgid << endl;
}

}


void grGeomHdl::voxelize_via_pvtessellation()
{
	//MK::EXECUTED FROM WITHIN PARALLEL REGION

	real_xyz h = localgrid.dcell;
	real_xyz hhalf = static_cast<real_xyz>(0.5) * h;

	p3d global_conti_origin = owner->thegrid.origin_cntum;
	vxl local2global = localgrid.origin_discr_lnk;
	p3d v = p3d();

	//##MK::potential further optimization trick percompute equidistance voxel center locations

	for(size_t lz = 0; lz < localgrid.nz; ++lz) {
		size_t zoff = lz * localgrid.nxy;
		v.z = global_conti_origin.z + hhalf + (h * static_cast<real_xyz>(local2global.z + lz));

		for(size_t ly = 0; ly < localgrid.ny; ++ly) {
			size_t yzoff = ly * localgrid.nx + zoff;
			v.y = global_conti_origin.y + hhalf + (h * static_cast<real_xyz>(local2global.y + ly));

			for(size_t lx = 0; lx < localgrid.nx; ++lx) {
				//continuum coordinate exemplarily: left voxel wall boundary is a global_conti_origin.x add half cell to get to center add integer multiple of entire cells to get to local cell center
				v.x = global_conti_origin.x + hhalf + (h * static_cast<real_xyz>(local2global.x + lx));

				nbp3d nbor_if_any = owner->pp3rve27.find_nearest_neighbor( v, Settings::PVTessKernelRadius ); //returns distance and unique ip name

				if ( nbor_if_any.uipid != numeric_limits<unsigned int>::max() ) {
					//MK::cubic voxel center of gravity has nothing to do with DAMASK simulation determined location of ips (perips and uips)!

					size_t cxyz = lx + yzoff;

					UIPIDField[cxyz] = nbor_if_any.uipid;

					if ( nbor_if_any.representative == YES ) { //point is part of representative dbcluster
						IsRepresentativeField[cxyz] = YES;
					}
				}
				//nothing remains to be done because if nobody was found UIPIDField was already set to ERROR FLAG numeric_limits<unsigned int>::max() upon initialization
			} //next vxl along +x
		} //next x vxlline stacked along +y
	} //next xy vxlslab stacked along +z

#pragma omp critical
{
	cout << "VoxelizationPerformed for grain " << cgid << endl;
}

}


void grGeomHdl::compute_sgndistfun_coarse()
{
	//initial allocation of SgnDistField1 was already everything set to h_o, so only resetting to h_i for all vxl closest to uipid referencing grain == cgid required
	//MK::signed-distance function sign convention

	real_sdf h = static_cast<real_sdf>(localgrid.dcell);
	//real_sdf h_o = static_cast<real_sdf>(-1.0) * h; //outside the grain - negative
	real_sdf h_i = static_cast<real_sdf>(+1.0) * h; //inside the grain - positive

	nuvxl = 0;

	for(size_t lz = 0; lz < localgrid.nz; ++lz) {
		size_t zoff = lz * localgrid.nxy;
		for(size_t ly = 0; ly < localgrid.ny; ++ly) {
			size_t yzoff = ly * localgrid.nx + zoff;
			for(size_t lx = 0; lx < localgrid.nx; ++lx) {
				//where?
				size_t cxyz = lx + yzoff;

				//which uip is closest?
				unsigned int uipid = UIPIDField[cxyz];

				//supporting at all a grain and if so which?
				if ( uipid != numeric_limits<unsigned int>::max() ) {
					unsigned int uipid2gid = owner->eid2gid( static_cast<size_t>(uipid) );
					GrainIDField[cxyz] = uipid2gid;

					if ( uipid2gid == cgid ) { //difficult to predict what is most likely, as AABB are tight about grain targetgid, so most voxel should be of the target grain however, shape non-convexity may render lower fraction, hence then more branch mispredictions

						//to avoid having now also fragments of periodic wraparound grains inside the voxel container consider in addition whether ips is one of the representative dbscan cluster
						if ( IsRepresentativeField[cxyz] == YES ) {
							//###MK::strictly speaking an assignment should only be made

							//cout << "cgid/lx/ly/lz/uipid2gid " << cgid << "\t\t" << lx << ";" << ly << ";" << lz << "\t\t" << uipid2gid << endl;
							SgnDistField1[cxyz] = h_i;
							nuvxl++;
						}
					}
				}
			} //next vxl along +x
		} //next xline stacked along +y
	} //next xyslab stacked along +z

#pragma omp critical
{
	cout << "Coarse signed distance function initialized for grain " << cgid << endl;
}

}


void grGeomHdl::compute_sgndistfun_fsm()
{
	//EXECUTED FROM WITHIN PARALLEL REGION

	//execute in essence the fast sweeping method which was originally proposed by
	//H. Zhao, Mathematics of Computation, 74 250, 2004, 603-627
	//to compute approximate solution to Eikonal equation to get signed-distance function for all vxl
	//here we utilize implementation of C. Mie{\ss}en, N. Velinov, G. Gottstein, L. A. Barrales-Mora, MSMSE, 25, 8 , 2017,
	//https://github.com/GraGLeS/GraGLeS_3D  commit id 69753e2
	//which is a modification ##MK

	//MK::accumulate results in distance buffer 2
	real_sdf h = static_cast<real_sdf>(localgrid.dcell);
	real_sdf candidate, i_slope, distToZero;

	//init outputBuffer SgnDistField2 was already prepared

	//other than in the implementation of Mie{\ss}en in- and output buffer of same size
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax, intersec_zmin, intersec_zmax;

	//mind that size_t to int cast is unsafe as multiple wraparound can occur!
	intersec_xmin = localgrid.origin_discr_own.x;
	intersec_ymin = localgrid.origin_discr_own.y;
	intersec_zmin = localgrid.origin_discr_own.z;
	intersec_xmax = localgrid.origin_discr_own.x + localgrid.nx;
	intersec_ymax = localgrid.origin_discr_own.y + localgrid.ny;
	intersec_zmax = localgrid.origin_discr_own.z + localgrid.nz;

	//mind naming convention of indexing order differs
	//column+NX*row+depth*NXY indexing though applies the same for both

	//Mie{\ss}en indexing order i,j,k <=> row, column, depth <=> y, x, z,
	//K\"uhbach indexing order i,j,k <=> x, y, z

	//mind further that, as both buffers are the same the following identities holds::
	//m_*Distance->getMinI() => intersec_imin;
	//m_*Distance->getMaxI() => intersec_imax;
	/*
	int m_inputDistance_getMinX = intersec_xmin;
	int m_inputDistance_getMaxX = intersec_xmax;
	int m_inputDistance_getMinY = intersec_ymin;
	int m_inputDistance_getMaxY = intersec_ymax;
	int m_inputDistance_getMinZ = intersec_zmin;
	int m_inputDistance_getMaxZ = intersec_zmax;
	*/

	int m_outputDistance_getMinX = intersec_xmin;
	int m_outputDistance_getMaxX = intersec_xmax;
	int m_outputDistance_getMinY = intersec_ymin;
	int m_outputDistance_getMaxY = intersec_ymax;
	int m_outputDistance_getMinZ = intersec_zmin;
	int m_outputDistance_getMaxZ = intersec_zmax;

	// first to updates layer by layer to take advantage of the order of point in memory - there are aligned layer by layer.
	for (int k = intersec_zmin; k < intersec_zmax - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX - 1;
					j++) {
				// x-direction forward
				if (j < intersec_xmax - 1 && i < intersec_ymax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j + 1, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j + 1, k)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i,
													j + 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j + 1, k)))
						m_outputDistance_setValueAt(i, j + 1, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j + 1, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j + 1, k)))
						m_outputDistance_setValueAt(i, j + 1, k, candidate);
				}
			}
		}

		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmax - 1; j > m_outputDistance_getMinX;
					j--) {
				// x-direction outputDistanceward
				//check for sign change
				if (j > intersec_xmin && i < intersec_ymax) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i,
													j - 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
						m_outputDistance_setValueAt(i, j - 1, k, candidate);
				} else {
					//<<<<<<< HEAD
					//					candidate = m_outputDistance_getValueAt(i, j, k)
					//							+ sgn(m_outputDistance_getValueAt(i, j - 1, k))
					//									* h;
					//					if (abs(candidate)
					//							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
					//=======
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ sgn(m_outputDistance_getValueAt(i, j - 1, k))
									* h;
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
						//>>>>>>> 66c3d6672001fb0db664a7cc036f13ecf8da0d05
						m_outputDistance_setValueAt(i, j - 1, k, candidate);
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
			for (int i = intersec_ymin; i < m_outputDistance_getMaxY - 1;
					i++) {
				if (j < intersec_xmax && i < intersec_ymax - 1) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i + 1, j, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i + 1, j, k)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i + 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i + 1, j, k)))
						m_outputDistance_setValueAt(i + 1, j, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i + 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i + 1, j, k)))
						m_outputDistance_setValueAt(i + 1, j, k, candidate);
				}
			}
		}

		for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
			for (int i = intersec_ymax - 1; i > m_outputDistance_getMinY;
					i--) {
				if (j < intersec_xmax && i > intersec_ymin) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i - 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i - 1, j, k)))
						m_outputDistance_setValueAt(i - 1, j, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i - 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i - 1, j, k)))
						m_outputDistance_setValueAt(i - 1, j, k, candidate);
				}
			}
		}
	}
	// update the point into the third dimensio. the strategy has to change to avoid unneccesary cache loads
	// the idea is to compare all points in one layer first to the next and go on:
	//TODO redist into the third direction
	// z forward:
	for (int k = intersec_zmin; k < m_outputDistance_getMaxZ - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
				// x-direction forward
				if (k < intersec_zmax - 1 && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j, k + 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j, k + 1)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i, j,
													k + 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k + 1)))
						m_outputDistance_setValueAt(i, j, k + 1, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j, k + 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k + 1)))
						m_outputDistance_setValueAt(i, j, k + 1, candidate);
				}
			}
		}
	}
	// z backward:
	for (int k = intersec_zmax - 1; k > m_outputDistance_getMinZ; k--) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
				// x-direction forward
				if (k > intersec_zmin && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j, k - 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j, k - 1)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i, j,
													k - 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k - 1)))
						m_outputDistance_setValueAt(i, j, k - 1, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j, k - 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k - 1)))
						m_outputDistance_setValueAt(i, j, k - 1, candidate);
				}
			}
		}
	}

	//no clamping operation

	//##MK::thread writing file with different name than any other thread (disjoint grainIDs) therefore can be done in parallel, but of course may not be efficient when only one harddrive
	if ( Settings::VisGrainFinalReconstruction == 0 || Settings::VisGrainFinalReconstruction == (cgid+1) ) {
		bool status = vtk_grfinalrecon( owner->thegrid, localgrid, GrainIDField, UIPIDField,
				SgnDistField1, SgnDistField2, cgid, owner->thisincrement );
	}

#pragma omp critical
{
cout << "Fast sweeping method applied to grain " << cgid << endl;
}

}


void grGeomHdl::debug_sdf_sphere_seed( const real_xyz R)
{
	//allocate buffers for signed distance computation
	//allocate GrainID and place ideal sphere located at box center with continuum scale radius R
	//set all voxel with center of gravity inside to GrainID this->cgid, all others to unsigned int::max()
	GrainIDField.resize( localgrid.nxyz, numeric_limits<unsigned int>::max() ); //#MK::was this->localgrid

	real_xyz half = static_cast<real_xyz>(0.5);
	p3d center = localgrid.box.center();

cout << "DEBUGGING SDF placing sphere radius " << R << " at center " << center << endl;

	real_xyz SQRR = SQR(R);
	for (size_t z = 0; z < localgrid.nz; ++z) {
		real_xyz cz = localgrid.origin_cntum.z + (half + static_cast<real_xyz>(z)) * localgrid.dcell;
		size_t zoff = z * localgrid.nxy;

		for (size_t y = 0; y < localgrid.ny; ++y) {
			real_xyz cy = localgrid.origin_cntum.y + (half + static_cast<real_xyz>(y)) * localgrid.dcell;
			size_t yzoff = y * localgrid.nx + zoff;

			for (size_t x = 0; x < localgrid.nx; ++x) {
				real_xyz cx = localgrid.origin_cntum.x + (half + static_cast<real_xyz>(x)) * localgrid.dcell;

				real_xyz SQRdist = SQR(cx-center.x) + SQR(cy-center.y) + SQR(cz-center.z);

				if ( SQRdist <= SQRR ) {
					GrainIDField[x+yzoff] = cgid;
				}
			}
		}
	}
}


void grGeomHdl::debug_sdf_sphere_guess()
{
	//guess initial SDF by setting all vxl with GrainID == this->cgid to +h and -h otherwise
	//MK::use distance buffer 1

cout << "DEBUGGING SDF estimating initial SDF" << endl;

	real_sdf h_o = static_cast<real_sdf>(-1.0) * static_cast<real_sdf>(localgrid.dcell);
	real_sdf h_i = static_cast<real_sdf>(+1.0) * static_cast<real_sdf>(localgrid.dcell);

	//init SDF everywhere to h_o, negative on the outside
	SgnDistField1.resize( localgrid.nxyz, h_o );

	//replace SDF inside the grain to h_i, positive on the inside
	//proposed "straightforward" initialization approach for discrete Poisson-Voronoi structures

	for(size_t c = 0; c < localgrid.nxyz; ++c) {
		if ( GrainIDField[c] == cgid ) {
			SgnDistField1[c] = h_i;
		}
	}
}

//##MK::DEBUG functions emulating data access strategy of Miessen to avoid modifying too strongly the cases
inline real_sdf grGeomHdl::m_inputDistance_getValueAt( const int row, const int col, const int dep )
{
	//input buffer in our case is SgnDistField1
	int x = 0; //##MK::by construction localgrid.origin_discr_own.x;
	int y = 0; //localgrid.origin_discr_own.y;
	int z = 0; //localgrid.origin_discr_own.z;
	int nx = localgrid.nx;
	int nxy = localgrid.nxy;
	return SgnDistField1.at( (col-x) + (row-y)*nx + (dep-z)*nxy );
}

inline real_sdf grGeomHdl::m_outputDistance_getValueAt( const int row, const int col, const int dep )
{
	//output buffer in our case is SgnDistField2
	int x = 0; //##MK::by construction localgrid.origin_discr_own.x;
	int y = 0; //localgrid.origin_discr_own.y;
	int z = 0; //localgrid.origin_discr_own.z;
	int nx = localgrid.nx;
	int nxy = localgrid.nxy;
	return SgnDistField2.at( (col-x) + (row-y)*nx + (dep-z)*nxy );
}

inline void grGeomHdl::m_outputDistance_setValueAt( const int row, const int col, const int dep, const real_sdf val )
{
	//output buffer in our case is SgnDistField2
	int x = 0; //##MK::by construction localgrid.origin_discr_own.x;
	int y = 0; //localgrid.origin_discr_own.y;
	int z = 0; //localgrid.origin_discr_own.z;
	int nx = localgrid.nx;
	int nxy = localgrid.nxy;
	SgnDistField2.at( (col-x) + (row-y)*nx + (dep-z)*nxy ) = val;
}


void grGeomHdl::debug_sdf_sphere_fsm()
{
cout << "DEBUGGING SDF computing SDF using FastSweepingMethod" << endl;

	//execute in essence the fast sweeping method which was originally proposed by
	//H. Zhao, Mathematics of Computation, 74 250, 2004, 603-627
	//to compute approximate solution to Eikonal equation to get signed-distance function for all vxl
	//here we utilize implementation of C. Mie{\ss}en, N. Velinov, G. Gottstein, L. A. Barrales-Mora, MSMSE, 25, 8 , 2017,
	//https://github.com/GraGLeS/GraGLeS_3D  commit id 69753e2
	//which is a modification ##MK

	//MK::accumulate results in distance buffer 2
	real_sdf h = static_cast<real_sdf>(localgrid.dcell);
	real_sdf candidate, i_slope, distToZero;

	//init outputBuffer SgnDistField2
	SgnDistField2.resize(localgrid.nxyz, static_cast<real_sdf>(FASTSWEEPING_INITIAL_VALUE) );
	//##MK::this implicitly clamps the signed distance function outside the grain

	//other than in the implementation of Mie{\ss}en in- and output buffer of same size
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax, intersec_zmin, intersec_zmax;

	//mind that size_t to int cast is unsafe as multiple wraparound can occur!
	intersec_xmin = localgrid.origin_discr_own.x;
	intersec_ymin = localgrid.origin_discr_own.y;
	intersec_zmin = localgrid.origin_discr_own.z;
	intersec_xmax = localgrid.origin_discr_own.x + localgrid.nx;
	intersec_ymax = localgrid.origin_discr_own.y + localgrid.ny;
	intersec_zmax = localgrid.origin_discr_own.z + localgrid.nz;

	//mind naming convention of indexing order differs
	//column+NX*row+depth*NXY indexing though applies the same for both

	//Mie{\ss}en indexing order i,j,k <=> row, column, depth <=> y, x, z,
	//K\"uhbach indexing order i,j,k <=> x, y, z

	//mind further that, as both buffers are the same the following identities holds::
	//m_*Distance->getMinI() => intersec_imin;
	//m_*Distance->getMaxI() => intersec_imax;
	/*
	int m_inputDistance_getMinX = intersec_xmin;
	int m_inputDistance_getMaxX = intersec_xmax;
	int m_inputDistance_getMinY = intersec_ymin;
	int m_inputDistance_getMaxY = intersec_ymax;
	int m_inputDistance_getMinZ = intersec_zmin;
	int m_inputDistance_getMaxZ = intersec_zmax;
	*/

	int m_outputDistance_getMinX = intersec_xmin;
	int m_outputDistance_getMaxX = intersec_xmax;
	int m_outputDistance_getMinY = intersec_ymin;
	int m_outputDistance_getMaxY = intersec_ymax;
	int m_outputDistance_getMinZ = intersec_zmin;
	int m_outputDistance_getMaxZ = intersec_zmax;

	// first to updates layer by layer to take advantage of the order of point in memory - there are aligned layer by layer.
	for (int k = intersec_zmin; k < intersec_zmax - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX - 1;
					j++) {
				// x-direction forward
				if (j < intersec_xmax - 1 && i < intersec_ymax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j + 1, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j + 1, k)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i,
													j + 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j + 1, k)))
						m_outputDistance_setValueAt(i, j + 1, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j + 1, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j + 1, k)))
						m_outputDistance_setValueAt(i, j + 1, k, candidate);
				}
			}
		}

		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmax - 1; j > m_outputDistance_getMinX;
					j--) {
				// x-direction outputDistanceward
				//check for sign change
				if (j > intersec_xmin && i < intersec_ymax) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i,
													j - 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
						m_outputDistance_setValueAt(i, j - 1, k, candidate);
				} else {
					//<<<<<<< HEAD
					//					candidate = m_outputDistance_getValueAt(i, j, k)
					//							+ sgn(m_outputDistance_getValueAt(i, j - 1, k))
					//									* h;
					//					if (abs(candidate)
					//							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
					//=======
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ sgn(m_outputDistance_getValueAt(i, j - 1, k))
									* h;
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j - 1, k)))
						//>>>>>>> 66c3d6672001fb0db664a7cc036f13ecf8da0d05
						m_outputDistance_setValueAt(i, j - 1, k, candidate);
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
			for (int i = intersec_ymin; i < m_outputDistance_getMaxY - 1;
					i++) {
				if (j < intersec_xmax && i < intersec_ymax - 1) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i + 1, j, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i + 1, j, k)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i + 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i + 1, j, k)))
						m_outputDistance_setValueAt(i + 1, j, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i + 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i + 1, j, k)))
						m_outputDistance_setValueAt(i + 1, j, k, candidate);
				}
			}
		}

		for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
			for (int i = intersec_ymax - 1; i > m_outputDistance_getMinY;
					i--) {
				if (j < intersec_xmax && i > intersec_ymin) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i - 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i - 1, j, k)))
						m_outputDistance_setValueAt(i - 1, j, k, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i - 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i - 1, j, k)))
						m_outputDistance_setValueAt(i - 1, j, k, candidate);
				}
			}
		}
	}
	// update the point into the third dimensio. the strategy has to change to avoid unneccesary cache loads
	// the idea is to compare all points in one layer first to the next and go on:
	//TODO redist into the third direction
	// z forward:
	for (int k = intersec_zmin; k < m_outputDistance_getMaxZ - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
				// x-direction forward
				if (k < intersec_zmax - 1 && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j, k + 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j, k + 1)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i, j,
													k + 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k + 1)))
						m_outputDistance_setValueAt(i, j, k + 1, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j, k + 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k + 1)))
						m_outputDistance_setValueAt(i, j, k + 1, candidate);
				}
			}
		}
	}
	// z backward:
	for (int k = intersec_zmax - 1; k > m_outputDistance_getMinZ; k--) {
		for (int i = intersec_ymin; i < m_outputDistance_getMaxY; i++) {
			for (int j = intersec_xmin; j < m_outputDistance_getMaxX; j++) {
				// x-direction forward
				if (k > intersec_zmin && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance_getValueAt(i, j, k)
							* m_inputDistance_getValueAt(i, j, k - 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance_getValueAt(i, j, k - 1)
								- m_inputDistance_getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance_getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance_getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance_setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance_getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance_getValueAt(i, j,
													k - 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k - 1)))
						m_outputDistance_setValueAt(i, j, k - 1, candidate);
				} else {
					candidate = m_outputDistance_getValueAt(i, j, k)
							+ (sgn(m_outputDistance_getValueAt(i, j, k - 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance_getValueAt(i, j, k - 1)))
						m_outputDistance_setValueAt(i, j, k - 1, candidate);
				}
			}
		}
	}

	//no clamping operation
}


void grGeomHdl::debug_sdf_sphere_exact( const real_xyz R)
{
	//compare with exact analytical distance function
	//MK::in distance buffer 1

cout << "DEBUGGING SDF compute analytical SDF" << endl;
	p3d center = localgrid.box.center();
	real_xyz half = static_cast<real_xyz>(0.5);

	for (size_t z = 0; z < localgrid.nz; ++z) {
		real_xyz cz = localgrid.origin_cntum.z + (half + static_cast<real_xyz>(z)) * localgrid.dcell;
		size_t zoff = z * localgrid.nxy;

		for (size_t y = 0; y < localgrid.ny; ++y) {
			real_xyz cy = localgrid.origin_cntum.y + (half + static_cast<real_xyz>(y)) * localgrid.dcell;
			size_t yzoff = y * localgrid.nx + zoff;

			for (size_t x = 0; x < localgrid.nx; ++x) {
				real_xyz cx = localgrid.origin_cntum.x + (half + static_cast<real_xyz>(x)) * localgrid.dcell;

				real_xyz SQRdist = SQR(cx-center.x) + SQR(cy-center.y) + SQR(cz-center.z);

				if ( GrainIDField[x+yzoff] == cgid )
					SgnDistField1[x+yzoff] = static_cast<real_sdf>(+1.0) * (R - sqrt(SQRdist)); //positive inside grain
				else
					SgnDistField1[x+yzoff] = static_cast<real_sdf>(-1.0) * (sqrt(SQRdist) - R); //negative outside
			}
		}
	}
}

void grGeomHdl::debug_sdf_sphere_report()
{
	//VTK file for each voxel, field data GrainID, distance buffer 1 analytical solution, distance buffer 2 FSM solution, difference
	//writes 3d positions to VTK file
cout << "DEBUGGING SDF reporting SDF results" << endl;

	string fn = owner->get_prefix() + ".DEBUG.SDF.Sphere.vtk";

	ofstream vtk;
	vtk.open( fn.c_str(),  ofstream::out | std::ofstream::trunc  );
	if ( vtk.is_open() == true ) {
		size_t nvertices = localgrid.nxyz;
		//construct header and point coordinates
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "DAMASKPDT Debugging discretization and signeddistance approximation sphere\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "\n";
		vtk << "POINTS " << nvertices << " double\n";
		real_xyz half = static_cast<real_xyz>(0.5);
		for (size_t z = 0; z < localgrid.nz; ++z) {
			real_xyz cz = localgrid.origin_cntum.z + (half + static_cast<real_xyz>(z)) * localgrid.dcell;
			for (size_t y = 0; y < localgrid.ny; ++y) {
				real_xyz cy = localgrid.origin_cntum.y + (half + static_cast<real_xyz>(y)) * localgrid.dcell;
				for (size_t x = 0; x < localgrid.nx; ++x) {
					real_xyz cx = localgrid.origin_cntum.x + (half + static_cast<real_xyz>(x)) * localgrid.dcell;

					vtk << cx << " " << cy << " " << cz << "\n";
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

		//add field data
		vtk << "POINT_DATA " << nvertices << "\n";
		vtk << "FIELD FieldData 4\n";
		vtk << "GrainID 1 " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << GrainIDField[i] << "\n";
		}
		vtk << endl;

		vtk << "SDFana 1 " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << SgnDistField1[i] << "\n";
		}
		vtk << endl;

		vtk << "SDFfsm 1 " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			vtk << SgnDistField2[i] << "\n";
		}
		vtk << endl;

		vtk << "SDFDiff 1 " << nvertices << " double\n";
		for (size_t i = 0; i < nvertices; ++i) {
			//vtk << abs((SgnDistField1[i] - SgnDistField2[i])) << "\n"; //absolute point-wise difference analytical solution - FSM estimate
			real_sdf rel =  static_cast<real_sdf>(100.0) * abs(SgnDistField1[i]-SgnDistField2[i]) / abs(SgnDistField1[i]);
			if ( isinf(rel) == false && isnan(rel) == false )
				vtk << rel << "\n";
			else
				vtk << static_cast<real_sdf>(1.0e9) << "\n";
		}
		vtk << endl;
		vtk.flush();
		vtk.close();
	}
	else {
cout << "DEBUGGING Signed Distance Function Sphere Reporting unable to write integration point positions to VTK file" << endl;
	}
}



specOutIncr::specOutIncr()
{
	thisincrement = numeric_limits<unsigned int>::max();
	thiswrittenincrement = numeric_limits<unsigned int>::max();

	thegrid = vxlgrd();

	rveBaseIni = bv3x3();
	rveBaseDef = bv3x3();

	myrank = MASTER;
	nranks = SINGLE;
	healthy = true;
}


specOutIncr::~specOutIncr()
{
	for (size_t mr = 0; mr < db.size(); mr++) {
		if (db.at(mr) != NULL) {
			delete db.at(mr);
			db.at(mr) = NULL;
		}
	}
	thegrid = vxlgrd();
}


void specOutIncr::init_mpi_derivedtypes()
{
	/*initializes user-defined MPI datatypes*/
	MPI_Type_contiguous(9, MPI_DOUBLE, &MPI_Tensor3x3_Double_Type);
	MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_VonMises_Double_Type);
	MPI_Type_contiguous(13, MPI_DOUBLE, &MPI_DisloSpatDistr_Double_Type);

	/*
	//MPI_Datatype MPI_3DFaceIO_Type;
	int elementCounts3[2] = {1, 2};
	MPI_Aint displacements3[2] = {0, 1 * MPIIO_OFFSET_INCR_DOUBLE};
	MPI_Datatype oldTypes3[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Type_create_struct(2, elementCounts3, displacements3, oldTypes3, &MPI_3DFaceIO_Type);
	*/

	MPI_Type_commit(&MPI_Tensor3x3_Double_Type);
	MPI_Type_commit(&MPI_VonMises_Double_Type);
	MPI_Type_commit(&MPI_DisloSpatDistr_Double_Type);

	//MPI_Type_commit(&MPI_3DFaceIO_Type);
}


bool specOutIncr::specout_read_header()
{
	//read header of a spectralOut file parsing Fortran integer tags
	double tic = MPI_Wtime();
	string fn = Settings::SpectralOutFilename + ".spectralOut";
	FILE * spec;
	spec = fopen( fn.c_str() , "rb" );
	if ( spec == NULL ) {
		string s = Settings::SpectralOutFilename + " either does not exist or is inaccessible!";
		stopping( s, myrank, 0 );
		return false;
	}
	else { //start reading file
if ( get_myrank() == MASTER )
	cout << "The metadata of the spectralOut file are as follows:" << endl;

		size_t probe = 0;
		int ftag = 0; //Fortran writes integer tags whose size indicates the number of data elements following
		//closed individually by the tag again each

		//##MK::newline as single character

		//load file, leading tag
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read load tag 1!", myrank, 0); fclose(spec); return false; }
		//value itself
		head.load.resize(ftag);
		probe = fread(&head.load[0], sizeof(char), ftag, spec);
		if ( probe != static_cast<size_t>(ftag) ) { stopping( "Unable to read load value!", myrank, 0); fclose(spec); return false; }
cout << "Load filename\t\t\t" << head.load << endl;
		//trailing tag
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read load tag 2!", myrank, 0); fclose(spec); return false; }

		//working dir
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read wdir tag 1!", myrank, 0); fclose(spec); return false; }
		head.wdir.resize(ftag);
		probe = fread(&head.wdir[0], sizeof(char), ftag, spec);
		if ( probe != static_cast<size_t>(ftag) ) { stopping( "Unable to read wdir value!", myrank, 0); fclose(spec); return false; }
cout << "Working directory\t\t" << head.wdir << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read wdir tag 2!", myrank, 0); fclose(spec); return false; }

		//geometry file
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read geom tag 1!", myrank, 0); fclose(spec); return false; }
		head.geom.resize(ftag);
		probe = fread(&head.geom[0], sizeof(char), ftag, spec);
		if ( probe != static_cast<size_t>(ftag) ) { stopping( "Unable to read geom value!", myrank, 0); fclose(spec); return false; }
cout << "Geometry filename\t\t" << head.geom << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read geom tag 2!", myrank, 0); fclose(spec); return false; }

		//grid size
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read grid tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != (5+3*4) ) { stopping( "File corrupted, keyword with 5 char and three int32 expected!", myrank, 0); fclose(spec); return false; }
		//skip 5 char keyword 'grid:'
		char skip5[5];
		probe = fread(&skip5[0], sizeof(char), 5, spec);
		if ( probe != 5 ) { stopping( "Unable to extract 5 char keyword grid:", myrank, 0); fclose(spec); return false; }
		int grid[3];
		probe = fread(&grid[0], sizeof(int), 3, spec);
		if ( probe != 3 ) { stopping( "Unable to extract 3 int grid sizes", myrank, 0); fclose(spec); return false; }
		if ( grid[0] == grid[1] && grid[0] == grid[2] ) { //cubic element grid as expected for DAMASK_spectral
			head.N = static_cast<unsigned int>(grid[0]);
			head.NX = static_cast<unsigned int>(grid[0]);
			head.NY = static_cast<unsigned int>(grid[1]);
			head.NZ = static_cast<unsigned int>(grid[2]);
			head.NXY = head.NX*head.NY;
			head.NXYZ = head.NX*head.NY*head.NZ;
		}
cout << "RVE domain resolution\t\t" << head.N << endl;
cout << "RVE domain shape\t\t" << head.NX << "\t\t" << head.NY << "\t\t" << head.NZ << endl;
cout << "RVE domain area/volume\t\t" << "\t\t" << head.NXY << "\t\t" << head.NXYZ << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read grid tag 2!", myrank, 0); fclose(spec); return false; }

		//geomSize
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read geomSize tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != (5+3*8) ) { stopping( "File corrupted, keyword with 5 char and three float64 expected!", myrank, 0); fclose(spec); return false; }
		//skip 5 char keyword 'size:'
		probe = fread(&skip5[0], sizeof(char), 5, spec);
		if ( probe != 5 ) { stopping( "Unable to extract 5 char keyword size:", myrank, 0); fclose(spec); return false; }
		double geomsizes[3]; //MK::required double as Fortran KIND=8 fp value
		probe = fread(&geomsizes[0], sizeof(double), 3, spec);
		if ( probe != 3 ) { stopping( "Unable to extract 3 float64 geomSizes", myrank, 0); fclose(spec); return false; }
		head.L = v3x1(
				static_cast<real_xyz>(geomsizes[0]),
				static_cast<real_xyz>(geomsizes[1]),
				static_cast<real_xyz>(geomsizes[2]) );
cout << "RVE domain geometry\t\t" << head.L.x << "\t\t" << head.L.y << "\t\t" << head.L.z << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read geomSize tag 2!", myrank, 0); fclose(spec); return false; }

		//materialpoint_results
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read matpres tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != (26+4) ) { stopping( "File corrupted, keyword with 26 char and one int32 expected!", myrank, 0); fclose(spec); return false; }
		char skip26[26];
		//skip 26 char keyword 'materialpoint_sizeResults:'
		probe = fread(&skip26[0], sizeof(char), 26, spec);
		if ( probe != 26 ) { stopping( "Unable to extract 26 char keyword materialpoint_sizeResults:", myrank, 0); fclose(spec); return false; }
		int mpr; //MK:: Fortran int KIND=4 type
		probe = fread(&mpr, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
		head.matpres = static_cast<unsigned int>(mpr);
cout << "Materialdata per IP\t\t" << head.matpres << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read matpres tag 2!", myrank, 0); fclose(spec); return false; }

		//loadcases
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read loadcases tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != (10+4) ) { stopping( "File corrupted, keyword with 10 char and one int32 expected!", myrank, 0); fclose(spec); return false; }
		char skip10[10];
		//skip 10 char keyword 'loadcases:'
		probe = fread(&skip10[0], sizeof(char), 10, spec);
		if ( probe != 10 ) { stopping( "Unable to extract 10 char keyword loadcases:", myrank, 0); fclose(spec); return false; }
		int lc;
		probe = fread(&lc, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
		head.loadcases = static_cast<unsigned int>(lc);
cout << "Loadcases in total\t\t" << head.loadcases << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read loadcases tag 2!", myrank, 0); fclose(spec); return false; }

		//init loadcasesmetadata
		head.loadcasesmeta.clear();
		for (unsigned int l = 0; l < head.loadcases; ++l) {
			head.loadcasesmeta.push_back( lcasemeta() );
		}

		//frequencies
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read freqs tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != static_cast<int>(12+head.loadcases*4) ) { stopping( "File corrupted, keyword with 12 char and loadcases*int32 expected!", myrank, 0); fclose(spec); return false; }
		char skip12[12];
		//skip 12 char keyword 'frequencies:'
		probe = fread(&skip12[0], sizeof(char), 12, spec);
		if ( probe != 12 ) { stopping( "Unable to extract 12 char keyword frequencies:", myrank, 0); fclose(spec); return false; }
		for (unsigned int l = 0; l < head.loadcases; ++l) {
			int fq;
			probe = fread(&fq, sizeof(int), 1, spec);
			if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
			head.loadcasesmeta.at(l).freqs = static_cast<unsigned int>(fq);
cout << "Frequency loadcase " << l << "\t\t" << head.loadcasesmeta.at(l).freqs << endl;
		}
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read freqs tag 2!", myrank, 0); fclose(spec); return false; }

		//times
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read times tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != static_cast<int>(6+head.loadcases*8) ) { stopping( "File corrupted, keyword with 6 char and loadcases*double64 expected!", myrank, 0); fclose(spec); return false; }
		char skip6[6];
		//skip 6 char keyword 'times:'
		probe = fread(&skip6[0], sizeof(char), 6, spec);
		if ( probe != 6 ) { stopping( "Unable to extract 6 char keyword times:", myrank, 0); fclose(spec); return false; }
		for (unsigned int l = 0; l < head.loadcases; ++l) {
			double tms;
			probe = fread(&tms, sizeof(double), 1, spec);
			if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
			head.loadcasesmeta.at(l).times = tms;
cout << "Timing of loadcase " << l << "\t\t" << head.loadcasesmeta.at(l).times << endl;
		}
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read times tag 2!", myrank, 0); fclose(spec); return false; }

		//logscales
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read logscales tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != static_cast<int>(10+head.loadcases*4) ) { stopping( "File corrupted, keyword with 12 char and loadcases*int32 expected!", myrank, 0); fclose(spec); return false; }
		//skip 10 char keyword 'logscales:'
		probe = fread(&skip10[0], sizeof(char), 10, spec);
		if ( probe != 10 ) { stopping( "Unable to extract 10 char keyword logscales:", myrank, 0); fclose(spec); return false; }
		for (unsigned int l=0; l < head.loadcases; ++l) {
			int lgsc;
			probe = fread(&lgsc, sizeof(int), 1, spec);
			if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
			head.loadcasesmeta.at(l).logscales = static_cast<unsigned int>(lgsc);
cout << "Is loadcase " << l << " logscaled\t\t" << head.loadcasesmeta.at(l).logscales << endl;
		}
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read logscales tag 2!", myrank, 0); fclose(spec); return false; }

		//increments
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read increments tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != static_cast<int>(11+head.loadcases*4) ) { stopping( "File corrupted, keyword with 11 char and loadcases*int32 expected!", myrank, 0); fclose(spec); return false; }
		char skip11[11];
		//skip 11 char keyword 'increments:'
		probe = fread(&skip11[0], sizeof(char), 11, spec);
		if ( probe != 11 ) { stopping( "Unable to extract 11 char keyword increments:", myrank, 0); fclose(spec); return false; }
		for (unsigned int l = 0; l < head.loadcases; ++l) {
			int incrs;
			probe = fread(&incrs, sizeof(int), 1, spec);
			if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
			head.loadcasesmeta.at(l).nincr = static_cast<unsigned int>(incrs);
cout << "Increments in loadcase " << l << "\t" << head.loadcasesmeta.at(l).nincr << endl;
		}
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read increments tag 2!", myrank, 0); fclose(spec); return false; }

		//startingIncrement
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read startingIncrement tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != (18+4) ) { stopping( "File corrupted, keyword with 18 char and one int32 expected!", myrank, 0); fclose(spec); return false; }
		char skip18[18];
		//skip 18 char keyword 'startingIncrement:'
		probe = fread(&skip18[0], sizeof(char), 18, spec);
		if ( probe != 18 ) { stopping( "Unable to extract 18 char keyword startingIncrement:", myrank, 0); fclose(spec); return false; }
		int sincr;
		probe = fread(&sincr, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to extract 1 int32 count", myrank, 0); fclose(spec); return false; }
		head.sincr = static_cast<unsigned int>(sincr);
cout << "Index of starting incr\t\t" << head.sincr << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read startingIncrement tag 2!", myrank, 0); fclose(spec); return false; }

		//eoh
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read eoh tag 1!", myrank, 0); fclose(spec); return false; }
		if ( ftag != 3 ) { stopping( "File corrupted, keyword with 3 char eoh expected!", myrank, 0); fclose(spec); return false; }
		char skip3[3];
		//skip 3 char keyword 'eoh'
		probe = fread(&skip3[0], sizeof(char), 3, spec);
		if ( probe != 3 ) { stopping( "Unable to extract 3 char keyword eoh", myrank, 0); fclose(spec); return false; }
cout << "EOH" << endl;
		probe = fread(&ftag, sizeof(int), 1, spec);
		if ( probe != 1 ) { stopping( "Unable to read eoh tag 2!", myrank, 0); fclose(spec); return false; }

		//compute last incr
		unsigned int incr = head.sincr;
		for (unsigned int loadcase = 0; loadcase < head.loadcases; ++loadcase) {
			for (unsigned int i = 0; i < head.loadcasesmeta.at(loadcase).nincr; ++i) { //##MK::mind that frequencies between loadcases can differ!
				incr++;
			}
		}
		head.lincr = incr;
cout << "LastIncrement\t\t\t" << head.lincr << endl;

		//store position of current filepointer as it points to first byte past the header
		head.FirstByteAfterHeader = static_cast<size_t>(ftell(spec));
cout << "Header offset is\t\t" << head.FirstByteAfterHeader << " Bytes" << endl;
	}

	if ( fclose(spec) != 0 ) {
		string s = Settings::SpectralOutFilename + " was not successfully closed!";
		stopping( s, myrank, 0 );
		return false;
	}

	//##MK::for the moment we utilize DAMASK_spectral therefore simple cubic reduced integration point element
	head.Nip = 1;
cout << "Header Nip\t\t\t" << head.Nip << endl;
	head.Ncp = head.NXYZ;
cout << "Header Ncp\t\t\t" << head.Ncp << endl;

	head.DataElementsPerIncrement = static_cast<size_t>(head.matpres)*static_cast<size_t>(head.Nip)*static_cast<size_t>(head.Ncp);
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::SpecOutReadHeader"));

	return true;
}


bool specOutIncr::specout_read_structure_homogenization()
{
	//the spectralOut file structure is partially strongly adapted to the simulated problem, therefore
	//meta data files are kept which specify the internal layout of the spectralOut file, in particular
	//how long the containers for physical data are expected and what their internal order is
	double tic = MPI_Wtime();

	ifstream metafile;
	string metaline;
	istringstream line;
	string datapiece;
	string k;
	string v;

	metafile.open(Settings::SpectralOutFilename + ".outputHomogenization");
	if ( metafile.is_open() == true) {
		//##MK::better functionality to auto-distinguish keywords from rubbish...
		//##MK::currently kick three header lines
		if (metafile.good() == true)	getline( metafile, metaline);
		if (metafile.good() == true)	getline( metafile, metaline);
		if (metafile.good() == true)	getline( metafile, metaline);

		while( metafile.good() == true ) {
			//read a line which contains potentially a keyword k and value v
			getline( metafile, metaline);
			if ( metaline.size() > 0 ) {
				istringstream line(metaline); //get character stream in string and check first and last character if indicating a DAMASK keyword, i.e. ( and )
				line >> k >> v;

				int test1 = k.compare(0,1,"(");
				int test2 = k.compare(k.length()-1,1,")");
				if ( test1 == 0 && test2 == 0 ) { //it was a keyword so store the value in map
					dlayout.homogenization.push_back( dlayoutnode(k, v) );

cout << dlayout.homogenization.at(dlayout.homogenization.size()-1).key <<
		"__" <<  dlayout.homogenization.at(dlayout.homogenization.size()-1).valstr <<
			"__" <<  dlayout.homogenization.at(dlayout.homogenization.size()-1).valui << endl;
				}
			}
		}//next potential keyword

		metafile.close();
	}
	else {
		stopping("Unable to open spectralOut homogenization metafile", myrank, 0);
		return false;
	}
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::MetaReadHomogenization"));

	return true;
}


bool specOutIncr::specout_read_structure_crystallite()
{
	double tic = MPI_Wtime();
	ifstream metafile;
	string metaline;
	istringstream line;
	string datapiece;
	string k;
	string v;

	//string s = ;
	metafile.open(Settings::SpectralOutFilename + ".outputCrystallite");
	if ( metafile.is_open() == true) {
		//##MK::better functionality to auto-distinguish keywords from rubbish...
		//##MK::currently kick three header lines
		if (metafile.good() == true)	getline( metafile, metaline);
		if (metafile.good() == true)	getline( metafile, metaline);
		if (metafile.good() == true)	getline( metafile, metaline);

		//##MK::kick first line which currently imports (plasticity) dislotwin

		while( metafile.good() == true ) {
			//read a line which contains potentially a keyword k and value v
			getline( metafile, metaline);
			if ( metaline.size() > 0 ) {
				istringstream line(metaline); //get character stream in string and check first and last character if indicating a DAMASK keyword, i.e. ( and )
				line >> k >> v;

				//##MK::check if k and v contain at all pieces of information
				dlayout.crystallite.push_back( dlayoutnode(k, v) );

cout << dlayout.crystallite.at(dlayout.crystallite.size()-1).key <<
		"__" <<  dlayout.crystallite.at(dlayout.crystallite.size()-1).valstr <<
			"__" <<  dlayout.crystallite.at(dlayout.crystallite.size()-1).valui << endl;
			}
		}//next potential keyword

		metafile.close();
	}
	else {
		stopping("Unable to open spectralOut constitutive metafile", myrank, 0);
		return false;
	}
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::MetaReadCrystallite"));

	return true;
}


bool specOutIncr::specout_read_structure_constitutive()
{
	double tic = MPI_Wtime();

	ifstream metafile;
	string metaline;
	istringstream line;
	string datapiece;
	string k;
	string v;

	//string s = ;
	metafile.open(Settings::SpectralOutFilename + ".outputConstitutive");
	if ( metafile.is_open() == true) {
		//##MK::better functionality to auto-distinguish keywords from rubbish...
		//##MK::currently kick three header lines
		if ( metafile.good() == true )	getline( metafile, metaline);
		if ( metafile.good() == true )	getline( metafile, metaline);
		if ( metafile.good() == true )	getline( metafile, metaline);

		//##MK::currently kick also header line which specifies constitutive model
		if ( metafile.good() == true )	getline( metafile, metaline);

		while( metafile.good() == true ) {
			//read a line which contains potentially a keyword k and value v
			getline( metafile, metaline);
			if ( metaline.size() > 0 ) {
				istringstream line(metaline); //get character stream in string and check first and last character if indicating a DAMASK keyword, i.e. ( and )
				line >> k >> v;

				//##MK::check if k and v contain at all pieces of information
				dlayout.constitutive.push_back( dlayoutnode(k, v) );

cout << dlayout.constitutive.at(dlayout.constitutive.size()-1).key <<
		"__" <<  dlayout.constitutive.at(dlayout.constitutive.size()-1).valstr <<
			"__" <<  dlayout.constitutive.at(dlayout.constitutive.size()-1).valui << endl;
			}
		}//next potential keyword

		metafile.close();
	}
	else {
		stopping("Unable to open spectralOut constitutive metafile", myrank, 0);
		return false;
	}
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::MetaReadConstitutive"));

	return true;
}



int specOutIncr::map_increments2ranks()
{
	double tic = MPI_Wtime();

	//assign DAMASK_spectral out increments to processes
	//taking into account that not every converged increment has rawdata written to spectralOutfile
	//the frequency at which rawdata are dumped is in addition load-case-dependent
	//thus we distribute only existent rawdata dumps

	int incr = head.sincr;
	int writtenincr = 0; //MK::int for type consistence, not yet any known increment
	for ( unsigned int loadcase = 0; loadcase < head.loadcases; ++loadcase) {
		//with which frequency were in this loadcase data dumped to the spectralOut file?
		unsigned int loadcase_dump_freq = head.loadcasesmeta.at(loadcase).freqs;

		for ( unsigned int li = 0; li < head.loadcasesmeta.at(loadcase).nincr; ++li ) { //##MK::check if freq > 1
			if ( li % loadcase_dump_freq == 0 ) { //so if a dump was written to file
				//##MK::round robin work partitioning of existent dump data
				int thisrank = writtenincr % nranks;
				//a valid rank id?
				if (thisrank >= MASTER && thisrank < nranks) {
					incr2rank.insert( pair<int,int>(incr,thisrank) );
					incr2wincr.insert( pair<int,int>(incr,writtenincr) );

					string s = "Increment " + to_string(incr) + " in loadcase " + to_string(loadcase);
					s += " has dump data, which are mapped to rank " + to_string(thisrank) + " identified as w(ritten)incr " + to_string(writtenincr);
					reporting( s, myrank, 0, false);

					writtenincr++;
				}
				else {
					stopping("Invalid rank to map an increment on", myrank, 0);
					return 0;
				}
			}
			else { //converged DAMASK spectral increment but for which no dump was written to spectralOut
				string s = "Increment " + to_string(incr) + " in loadcase " + to_string(loadcase) + " has no written dump data in spectralOut file!";
				reporting( s, myrank, 0, false);
			}
			incr++;
		} //check next for written data for next converged DAMASK increment
	}
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "WorkScheduling::MapIncrements2Ranks"));

	//all increments have been assigned to the processes, although not guaranteeing optimal load imbalance,
	//however, for DAMASK number of elements does not change, so only moderate effect of workload, unless shape change is significant and phase transformations change
	//volume of phases which require dissimilar processing changes drastically...
	return 1;
}


void specOutIncr::map_meshelements2threads()
{
	//functionality which, based on knowing the extent of the global mesh, partitions the elements
	//into xy-slabs along the z direction by generating a number of memRegion class objects. Each of which
	//acts as an element of a database carrying all heavydata for the meshelements within a single slab
	//In effect, this is a priori static load balancing across threads, based on meshelement integration points (ips)
	//at highest granularity
	//MK::for most postprocessing functionalities the workload per meshelement ip is usually comparable
	//therefore assigning millions of ips to a thread as a first guess strategy should generate already a
	//fairly even load distribution across the threads
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		bool mehealthy = true;

		#pragma omp master
		{
			for ( unsigned int mr = 0; mr < nt; mr++ ) {
				//master thread only initializes the database, thus pointers to memRegion class instances
				//are stored in memory of the master thread
				db.push_back(NULL);
			}
		}
		//necessary, as there is no implicit barrier for the OMP master construct
		//and we need the container db prepared to populate it with pointers to allocated objects
		#pragma omp barrier

		size_t neipid = static_cast<size_t>(head.Nip) * static_cast<size_t>(head.Ncp); //##MK::size of spectralOut materialpoint_sizeResults does not matter as we not utilize all pieces of information therein
		size_t MeshElementsPerThreadTarget = neipid / static_cast<size_t>(nt);
		//##MK::such workload may not be an integer multiple of the SIMD length nor result in the same number of elements per thread
		//therefore (scalar) fallback of remainder in thread with highest mt (aka the last thread) required

		if ( MeshElementsPerThreadTarget < 1 ) {
			stopping( "Illegal meshelement partitioning", myrank, mt );
			mehealthy = false;
		}
		else {
			//allocate memRegion class object in thread-local memory by allocating locally and ! making first-touch access !
			memRegion* me = NULL;
			try { me = new memRegion; }
			catch (bad_alloc &ompcroak) {
				stopping( "Unable to allocate memRegion", myrank, mt);
				mehealthy = false;
			}

			if ( mehealthy == true ) {
				//first touch and initialize
				me->owner = this; 	//MK::here this is necessary because we want to assign a pointer to the calling object to the memRegion class object such knows who's increment read instance it
									//to allow the memRegion access to higher level data...

				//MK::WE UTILIZE THAT ORDER OF MESH ELEMENTS IS IMPLICIT x+y*NX+z*NX*NY
				size_t s = static_cast<size_t>(mt) * MeshElementsPerThreadTarget;
				size_t e = s + MeshElementsPerThreadTarget;
				if ( mt == (nt-1) ) { //the last one
					e = neipid;
				}

				me->eipid_start = s;
				me->eipid_end = e;
				me->eipid_n = e - s; //MK explicit computation necessary

				me->grid.owner = me;
				//me->grid.init();
				me->homogenization.owner = me;
				//me->homogenization.init();
				me->crystallite.owner = me;
				//me->crystallite.init();
				me->constitutive.owner = me;
				//me->constitutive.init();

				//pass pointer to memRegion class object to owner is writing to shared vector
				string mess = "MemRegion " + to_string(mt) + " eipid_s/e/n ";
				mess += to_string(static_cast<size_t>(me->eipid_start)) + ";";
				mess += to_string(static_cast<size_t>(me->eipid_end)) + ";";
				mess += to_string(static_cast<size_t>(me->eipid_n));

				#pragma omp critical
				{
					reporting(mess, myrank, mt, true);
					db.at(mt) = me;
				}
			} //linking and first touch done
		}

		//identify whether any thread failed
		#pragma omp critical
		{
			if ( mehealthy == false )
				healthy = false;
		}
	}
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshElements2Ranks" + to_string(thisincrement)));
}


void specOutIncr::grid_initial_configuration()
{
	//transform in each memRegion elem ids to x0,y0,z0 thread parallel on db
	//utilize fact that DAMASK_spectral element order is implicit x+y*NX+z*NX*NY
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		db.at(mt)->grid.eipid2xyz0();
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::InitialCoordinates" + to_string(thisincrement)));
}


real_xyz specOutIncr::rve_volume_total()
{
	real_xyz res = 0.0;
	
	unsigned int nt = db.size();
	unsigned int mr = MASTER;
	for ( 		; mr < nt; mr++ ) { //collect local averaging results
		res += db.at(mr)->crystallite.compute_total_V();
	}
	
	return res;
}


t3x3 specOutIncr::rve_volume_averaged_defpgradient()
{
	t3x3 rve_avg = t3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //initialized with identity deformation gradient tensor, hence if applied to base vectors we remain in the initial configuration

	unsigned int nt = db.size();
	unsigned int mr = MASTER;
	for (		; mr < nt; mr++ ) { //collect local averaging results
		t3x3 mr_avg = db.at(mr)->crystallite.compute_average_Fp();

		rve_avg.add( mr_avg, static_cast<real_m33>(1.0) );
	}
	if ( nt > 0 ) { //summarize into RVE global result
		real_m33 scaler = static_cast<real_m33>(nt);
		rve_avg.div( scaler );
		return rve_avg;
	}
	//implicit else
	complaining("Unable to compute RVE averaged F, instead returning identity", myrank, 0 );
	return t3x3(); //if unable return identity tensor
}


t3x3 specOutIncr::rve_volume_averaged_defgradient()
{
	t3x3 rve_avg = t3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //initialized with identity deformation gradient tensor, hence if applied to base vectors we remain in the initial configuration

	unsigned int nt = db.size();
	unsigned int mr = MASTER;
	for (		; mr < nt; mr++ ) { //collect local averaging results
		t3x3 mr_avg = db.at(mr)->crystallite.compute_average_F();

		rve_avg.add( mr_avg, static_cast<real_m33>(1.0) );
	}
	if ( nt > 0 ) { //summarize into RVE global result
		real_m33 scaler = static_cast<real_m33>(nt);
		rve_avg.div( scaler );
		return rve_avg;
	}
	//implicit else
	complaining("Unable to compute RVE averaged F, instead returning identity", myrank, 0 );
	return t3x3(); //if unable return identity tensor
}


t3x3 specOutIncr::rve_volume_averaged_piolastress()
{
	t3x3 rve_avg = t3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //initialized with zero-valued piola stress tensor, hence if applied to base vectors we remain in the initial configuration

	unsigned int nt = db.size();
	unsigned int mr = MASTER;
	for ( 		; mr < nt; mr++ ) { //collect local averaging results
		t3x3 mr_avg = db.at(mr)->crystallite.compute_average_P();

		rve_avg.add( mr_avg, static_cast<real_m33>(1.0) );
	}
	if ( nt > 0 ) { //summarize into RVE global result
		real_m33 scaler = static_cast<real_m33>(mr);
		rve_avg.div( scaler );
		return rve_avg;
	}
	//implicit else
	complaining("Unable to compute RVE averaged P, instead returning zero-valued tensor", myrank, 0 );
	return t3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

t3x3 specOutIncr::rve_volume_averaged_truestrain( t3x3 const & Fav, t3x3 const & Pav )
{
	bool mystatus = true;
	t3x3 TrueStrain = failt3x3();
	
	mystatus = computeStrainTensor( Fav, TrueStrain );

	//##MK::DEBUG
	if (mystatus == true)
cout << "RVE UStrain " << TrueStrain << endl;
	else
cout << "RVE UStrain " << "IMKL call did not converge or failed!" << endl;
	//##MK::END OF DEBUG

	return TrueStrain;
}


t3x3 specOutIncr::rve_volume_averaged_cauchystress( t3x3 const & Fav, t3x3 const & Pav )
{
	bool mystatus = true;
	t3x3 Cauchyy = failt3x3();

	mystatus = computeCauchy( Fav, Pav, Cauchyy );

	//##MK::DEBUG
	if (mystatus == true)
cout << "RVE Cauchy " << Cauchyy << endl;
	else
cout << "RVE Cauchy " << "IMKL call did not converge or failed!" << endl;
	//##MK::END OF DEBUG
	
	return Cauchyy;
}


vMises specOutIncr::rve_volume_averaged_scalars( t3x3 const & eps, t3x3 const & cau)
{
	vMises res = vMises(computeMises(eps, false), computeMises(cau, true));
cout << "RVE vMisesStrain " << res.vMisesEquivStrain << endl;
cout << "RVE vMisesStress " << res.vMisesEquivStress << endl;
	return res;
}


void specOutIncr::analyze_addStrainTensors()
{
	double tic = MPI_Wtime();

	//##MK::assure that IntelMKL is not called multithreaded here! only one should be multithreaded here SVD of 3x3 tensors is peanuts therefore run IntelMKL sequentially on every thread

	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		bool mystatus = true;
		memRegion* thisregion = db.at(mt);
		size_t eipid_n = thisregion->eipid_n;

		for ( size_t i = 0; i < eipid_n; ++i ) {
			t3x3 Fnow = thisregion->crystallite.F.at(i);

			t3x3 TrueStrain = failt3x3();
			mystatus = computeStrainTensor( Fnow, TrueStrain );

/*
 	 	 	//##MK::DEBUG
			#pragma omp critical
			{
				if (mystatus == true)
cout << "ThreadID " << tid << " Gridpoint " << i << " --> " << "UStrain " << TrueStrain << endl;
				else
cout << "ThreadID " << tid << " Gridpoint " << i << " --> " << "IMKL call did not converge or failed!" << endl;
			}
			//##MK::END OF DEBUG, in production version comment out this critical region!
*/

			thisregion->crystallite.straintensor.push_back( TrueStrain );

		} //threadlocal processing of tensors
		//MK::as the IntelMKL is threadsafe multiple calls to it dont interfere
		//MK::however for optimal performance choose either MKL_NUM_THREADS = 1 && OMP_NUM_THREADS max,
		//in this case all physical cores call in parallel the IntelMKL but make a quasi-sequential call
		//other extreme is setting MKL_NUM_THREADS max and OMP_NUM_THREADS = 1 then
		//processing of grid data proceeds sequentially but with individually potentially threadparallel IMKL calls,
		//however as the granularity in this case is fine, and for 3x3 tensor the problem size either
		//way much too small to justify in every entry the setup of a threading region within MKL, very likely
		//option 1 --- making individual sequential calls to IMKL but having OMP threads per MPI process working on data is
		//most efficient

	} //end of parallel region

/*	//perform I/O
	for (size_t mr = 0; mr < db.size(); ++mr) {
		//##MK::DEBUG for now!
		vector<vMises>* here = NULL;
		here = &db.at(mr)->crystallite.scalars;
		size_t nr = here->size();
		for ( size_t r = 0; r < nr; ++r) {
			cout << here->at(r);
		}
	}*/
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddStrainTensors" + to_string(thisincrement)));
}


void specOutIncr::analyze_addCauchy()
{
	double tic = MPI_Wtime();

	//##MK::assure that IntelMKL is not called multithreaded here! only one should be multithreaded here SVD of 3x3 tensors is peanuts therefore run IntelMKL sequentially on every thread

	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		bool mystatus = true;
		memRegion* thisregion = db.at(mt);
		size_t eipid_n = thisregion->eipid_n;

		for ( size_t i = 0; i < eipid_n; ++i ) {
			t3x3 Fnow = thisregion->crystallite.F.at(i);
			t3x3 Pnow = thisregion->crystallite.P.at(i);

			t3x3 Cauchyy = failt3x3();

			mystatus = computeCauchy( Fnow, Pnow, Cauchyy );

/*
			//##MK::DEBUG
			#pragma omp critical
			{
				if (mystatus == true)
cout << "ThreadID " << tid << " Gridpoint " << i << " --> " << "Cauchy " << Cauchyy << endl;
				else
cout << "ThreadID " << tid << " Gridpoint " << i << " --> " << "IMKL call did not converge or failed!" << endl;
			}
			//##MK::END OF DEBUG, in production version comment out this critical region!
*/

			thisregion->crystallite.stresstensor.push_back( Cauchyy );

		} //threadlocal processing of tensors
	} //end of parallel region
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddCauchy" + to_string(thisincrement)));
}


void specOutIncr::analyze_addVonMises()
{
	double tic = MPI_Wtime();

	//##MK::assure that IntelMKL is not called multithreaded here! only one should be multithreaded here SVD of 3x3 tensors is peanuts therefore run IntelMKL sequentially on every thread

	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		memRegion* thisregion = db.at(mt);

		size_t eipid_n = thisregion->eipid_n;

		for ( size_t i = 0; i < eipid_n; ++i ) {

			t3x3 Strainn = thisregion->crystallite.straintensor.at(i);
			t3x3 Cauchyy = thisregion->crystallite.stresstensor.at(i);

			vMises res = vMises(computeMises(Strainn, false), computeMises(Cauchyy, true));

/*
			//##MK::DEBUG
			#pragma omp critical
			{
cout << "ThreadID " << tid << " Gridpoint " << i << " --> " << "vMisesUStrain " << res.vMisesEquivStrain << " vMisesCauchyStress " << res.vMisesEquivStress << endl;
			}
			//##MK::END OF DEBUG, in production version comment out this critical region!
*/

			thisregion->crystallite.scalars.push_back( res );
		} //threadlocal processing of tensors
	} //end of parallel region
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddVonMises" + to_string(thisincrement)));
}


void specOutIncr::analyze_addDisplacements()
{
	//compute displacement (fluctuant and average) or fall back to zero
	double tic = MPI_Wtime();

	//get integration point grid of the RVE in initial configuration
	unsigned int NXold = head.N;
	unsigned int NYold = head.N;
	unsigned int NZold = head.N;
	unsigned int NXYZ = NXold*NYold*NZold;

	//MK::assessing addDisplacement.py with cuboidal RVEs revealed, that resolutions like 31x07x11
	//become reshaped to 11 x 7 x 31 resulting in a Fourier grid on NI <=> original NZ 11, NJ <=> orgNY 7 NK = orgNZ//2+1 16
	//this is the python command causing it: table.data[:,:9].reshape(grid[2],grid[1],grid[0],3,3)

	unsigned int NXnew = NZold;
	unsigned int NYnew = NYold;
	unsigned int NZnew = NXold;
	unsigned int grid[3] = {NXnew, NYnew, NZnew};

	//set of unique coordinates
	set<real_xyz> sx;
	set<real_xyz> sy;
	set<real_xyz> sz;

	//load point grid, F component, and on the fly detection of unique coordinates
	vector<p3d> p;
	vector<t3x3> f;

	for (size_t mr = 0; mr < db.size(); mr++) {
		size_t ni = db.at(mr)->grid.xyz0.size();
		vector<p3d>* coordinates_are_here = &db.at(mr)->grid.xyz0;
		for (size_t i = 0; i < ni; ++i ) {
			p.push_back(coordinates_are_here->at(i));
			sx.insert( coordinates_are_here->at(i).x );
			sy.insert( coordinates_are_here->at(i).y );
			sz.insert( coordinates_are_here->at(i).z );
//cout << p.back() << endl;
		}

		size_t nj = db.at(mr)->crystallite.F.size();
		vector<t3x3>* defgrads_are_here = &db.at(mr)->crystallite.F;
		for (size_t j = 0; j < nj; ++j) {
			f.push_back(defgrads_are_here->at(j));
//cout << f.back() << endl;
		}
	}
cout << "Loaded xyz and F" << endl;

	//get unique coordinates and size of the mesh
	vector<real_xyz> ux;
	ux.assign( sx.begin(), sx.end() );
	vector<real_xyz> uy;
	uy.assign( sy.begin(), sy.end() );
	vector<real_xyz> uz;
	uz.assign( sz.begin(), sz.end() );

	//aabb3d corner = get_corners( ux, uy, uz );

cout << "Unique integration point grid coordinates and AABB about these identified" << endl;

	//Intel MKL requires scaling factors consistent with precision
	double size[3] = { static_cast<double>(head.L.x), static_cast<double>(head.L.y), static_cast<double>(head.L.z) };


	vector<rfftn*> defgrad_forwardFFT;
	unsigned int nFWTransforms = 0;

	vector<irfftn*> displ_backwardFFT;
	unsigned int nBKTransforms = 0;

	vector<d3d> displ_flu; //vector carries pointer pointer to heavy data on heap
	displ_flu.reserve(NXYZ);
	vector<d3d> displ_avg;
	displ_avg.reserve(NXYZ);


	//forward 3d fourier transform of individual deformation gradient tensor components
	nFWTransforms = fourier_transform_defgradient( f, grid, size, defgrad_forwardFFT );

	if ( nFWTransforms == 9 ) {
cout << "All components of the deformation gradient tensor successfully Fourier transformed" << endl;

		//compute average displacements based on the 9 successful forwardFFTs
		displacement_average( defgrad_forwardFFT, grid, size, displ_avg );
cout << "All average displacements computed based on forward Fourier transforms" << endl;

		//compute fluctuant displacements including three inverse 3d fourier transform results for coordinates 1-->0, 2-->1, and 3-->2
		nBKTransforms = displacement_fluctuations( defgrad_forwardFFT, grid, size, displ_backwardFFT, displ_flu );

		if ( nBKTransforms == 3 ) { //if also these were successful so copy in threadlocal memory
cout << "All fluctuant displacements computed from successful inverse Fourier transform" << endl;

			#pragma omp parallel
			{
				unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

				memRegion* thisregion = db.at(mt);
				size_t eipid_s = thisregion->eipid_start;
				size_t eipid_n = thisregion->eipid_n;

				thisregion->grid.dxyz_avg.reserve(eipid_n);
				thisregion->grid.dxyz_flu.reserve(eipid_n);
				for ( size_t i = 0; i < eipid_n; ++i ) {
					thisregion->grid.dxyz_avg.push_back( displ_avg.at(eipid_s+i) );
					thisregion->grid.dxyz_flu.push_back( displ_flu.at(eipid_s+i) );
				}
			} //end of parallel region with implicit barrier

		}
		else {
			cout << "Error in inverse Fourier transform during computing fluctuant contributions" << endl;
		}
	}
	else {
		cout << "Error in forward Fourier transform during computing deformation gradient transforms" << endl;
	}

	//clear components ffts, if existent
	for (size_t i = 0; i < defgrad_forwardFFT.size(); ++i ) {
		if ( defgrad_forwardFFT.at(i) != NULL ) {
			delete defgrad_forwardFFT.at(i);
			defgrad_forwardFFT.at(i) = NULL;
		}
	}

	//clear displacement iffts, if existent
	for(size_t i = 0; i < displ_backwardFFT.size(); ++i) {
		if ( displ_backwardFFT.at(i) != NULL ) {
			delete displ_backwardFFT.at(i);
			displ_backwardFFT.at(i) = NULL;
		}
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddDisplacements" + to_string(thisincrement)));
}


void specOutIncr::analyze_ignoreDisplacements()
{
	//MK::add zero displacements to allow utilizing the same periodic grid construction irrespective of displacement handling
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		memRegion* thisregion = db.at(mt);
		size_t eipid_n = thisregion->eipid_n;

		thisregion->grid.dxyz_avg.reserve(eipid_n);
		thisregion->grid.dxyz_flu.reserve(eipid_n);
		for ( size_t i = 0; i < eipid_n; ++i ) {
			thisregion->grid.dxyz_avg.push_back( d3d(0.0, 0.0, 0.0) );
			thisregion->grid.dxyz_flu.push_back( d3d(0.0, 0.0, 0.0) );
		}
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::IgnoreDisplacements" + to_string(thisincrement)));
}


void specOutIncr::analyze_addRVEAverages( const unsigned int lc, const unsigned int li, const unsigned int glbincrid )
{
	double tic = MPI_Wtime();
	avg.push_back( rveAverageResults() );
	avg.back().loadcaseID = lc;
	avg.back().localincrID = li;
	avg.back().globalincrID = glbincrid;

	//get RVE-volume average value of F and P
	avg.back().Vtotal = rve_volume_total();
cout << "RVETotalVolume V is" << endl << avg.back().Vtotal << endl;

	avg.back().Fpavgrve = rve_volume_averaged_defpgradient();
cout << "RVEVolAveraged Fp is" << endl << avg.back().Fpavgrve << endl;

	avg.back().Favgrve = rve_volume_averaged_defgradient();
cout << "RVEVolAveraged F is" << endl << avg.back().Favgrve << endl;

	avg.back().Pavgrve = rve_volume_averaged_piolastress();
cout << "RVEVolAveraged P is" << endl << avg.back().Pavgrve << endl;

	avg.back().Strainavgrve = rve_volume_averaged_truestrain( avg.back().Favgrve, avg.back().Pavgrve );
cout << "RVEVolAveraged TrueStrain is" << endl << avg.back().Strainavgrve << endl;

	avg.back().Cauchyavgrve = rve_volume_averaged_cauchystress( avg.back().Favgrve, avg.back().Pavgrve );
cout << "RVEVolAveraged CauchyStress is" << endl << avg.back().Cauchyavgrve << endl;

	avg.back().Equivavgrve = rve_volume_averaged_scalars( avg.back().Strainavgrve, avg.back().Cauchyavgrve );
cout << "RVEVolAveraged vonMises stress " << avg.back().Equivavgrve.vMisesEquivStress << " Pa, true strain " << avg.back().Equivavgrve.vMisesEquivStrain << endl;

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddRVEAverages" + to_string(thisincrement)));
}


void specOutIncr::analyze_ipgrid_displacements()
{
	double tic = MPI_Wtime();
	//##MK::so far only for cubic domains
	rveBaseIni = bv3x3(
			//first base vector parallel to axis 1 initial undeformed configuration
			static_cast<real_xyz>(head.L.x),
			static_cast<real_xyz>(0.0),
			static_cast<real_xyz>(0.0),

			//second base column vector parallel to axis 2 initial undeformed configuration
			static_cast<real_xyz>(0.0),
			static_cast<real_xyz>(head.L.y),
			static_cast<real_xyz>(0.0),

			//third ... parallel to axis 3 initial i.e. undeformed configuration
			static_cast<real_xyz>(0.0),
			static_cast<real_xyz>(0.0),
			static_cast<real_xyz>(head.L.z) );

	//now during simulated deformation on this hexahedral RVE volume
	//there is an RVE-volume-averaged deformation gradient \bar{F} tensor, which modifies
	//the length of these base vectors and the directions into which they point ``aka the deformed configuration''
	//these averaged tensors are stored in this->avg.Favgrve for instance

	//base vectors in deformed configuration = Frveavg left multiplied bv3x3
cout << "Base column vectors in initial configuration" << endl;
cout << rveBaseIni << endl;

	rveBaseDef = leftmult( avg.back().Favgrve, rveBaseIni );

cout << "Base vectors in deformed configuration" << endl;
cout <<  rveBaseDef << endl;

	//make these distorted base column vector matrix known to all memory regions
	//because they are needed to compute the periodic replica images of the integration point grid when Settings::AnalyzeStateVarSpatialDistr == true
	for (size_t mr = 0; mr < db.size(); mr++ ) {
		db.at(mr)->grid.perOffset = rveBaseDef;
	}

	if ( Settings::AnalyzeAddDisplacements == true && Settings::ThrshldAgainstDefConfig == true ) {
		reporting("Deformation-induced integration point grid displacements/distortions are accounted for", myrank, 0, false);

		if ( head.NX == head.NY && head.NX == head.NZ ) {
			//##MK::IntelMKL-powered computing of individual integration point grid new coordinates
			//is currently supported for cubic point grid with size [1.0 1.0 1.0] only

			//##MK::implement loading of otherwise computed displacements --> inject them here!

			//MK::do not enclose in parallel region as IntelMKL is called which may parallelize internally
			//if utilizing as the threaded version of the library and MKL_NUM_THREADS > 1
			analyze_addDisplacements();
		}
		else {
			complaining("The integration point grid displacements/distortions are not accounted for", myrank, 0);

			//computation of displacements for hexahedral grid not yet possible, either feed postResults + addDisplacement data
			//##MK::if so inject code for this here

			//or stay in the initial configuration
			analyze_ignoreDisplacements();
		}
	}
	else {
		complaining("The integration point grid displacements/distortions are not accounted for", myrank, 0);

		//computation of displacements either not supported or desired, that means we stay in the initial configuration

		analyze_ignoreDisplacements();
	} //assuming distances in undeformed configuration
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::Displacements" + to_string(thisincrement)));
}

//##MK::strictly speaking these functions must be called from an memRegion object within a parallel region
//during MPI I/O such that the threads automatically assign the data to their local containers
//##potential hybrid parallel HDF5 compilation and support, albeit we make no MPI calls from within OpenMP threads





//##MK::this hybrid parallel version 1 does ot crash but MPI_File_read does not read something, for now, utilize version 2 instead
//##MK::and additionally assumes every converged increment has dump data...
/*
void specOutIncr::specout_read_heavydata1()
{
	//Fortran spectralOut writes in double precision (KIND=8), in C this is MPI_DOUBLE
	//provided the generating system and the postprocessing system have the same endianness
	//this is only portable, as the corresponding MPI I/O routines in DAMASK_spectral when using the same system
	//for x86 Intel/AMD-based workstations and clusters should be no problem, usually they are little endian
	//for supercomputers of Cray and some Sun system this must be reconsidered by implementing
	//###MK::endianness checks
	double* buf = NULL;
	size_t nalloc = 0;
	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;

	#pragma omp parallel shared(buf,nalloc,ioReadFileHdl,ioReadFileStatus)
	{
		int tid = omp_get_thread_num();
		bool mehealthy = true;

		//thread local copies of data organization constructs
		//how many to read in total, an element is in the following understood as a homog+cryst+consti data package for one integration point of one mesh element
		size_t ElementsTotal = static_cast<size_t>(head.Nip) * static_cast<size_t>(head.Ncp);
		//we do not read all data at once but in multiples of the stripe size of the underlying raid system
		size_t ElementsPerBlockTarget = 1 + (Performance::MPIIOStripeSize / static_cast<size_t>(sizeof(MPI_DOUBLE)*head.matpres) );		//how many to read with current block
		size_t ElementsPerBlockNow = (ElementsTotal > ElementsPerBlockTarget) ? ElementsPerBlockTarget : ElementsTotal;
		//for how many ips do we have heavy data read already from the file
		size_t ElementsReadAlready = 0;
		long long offset = 0;

		#pragma omp master
		{
cout << "ElementsTotal/BlockTarget/Now " << ElementsTotal << ";" << ElementsPerBlockTarget << ";" << ElementsPerBlockNow << endl;
			//MK::MPI_THREAD_FUNNELED, only master thread makes MPI calls!
			nalloc = ElementsPerBlockNow * static_cast<size_t>(head.matpres);
			try { buf = new double[nalloc]; }
			catch (std::bad_alloc &croak) {
				stopping("Unable to allocate MPI I/O read buffer at first place", myrank, 0);
				mehealthy = false;
				healthy = false;
			}
		}
		#pragma omp barrier

		if (healthy == false) {
			#pragma omp critical
			{
				stopping("Master thread was unable to allocate buffer for MPI/IO", myrank, tid); //ordered couting
				//MK::do not return or alike, because we are not allowed to simply break out of OMP parallel region
			}
		}
		else {
			//set starting point to read heavy data for particular increment
			//MK::C-style indexing, i.e. starting at 0, instead of Fortran starting at 1
			long long databytes_preceeding_incs = static_cast<long long>(thisincrement);
			databytes_preceeding_incs *= static_cast<long long>(8);
			databytes_preceeding_incs *= static_cast<long long>(head.matpres);
			databytes_preceeding_incs *= static_cast<long long>(head.Nip);
			databytes_preceeding_incs *= static_cast<long long>(head.Ncp);

			offset = static_cast<long long>(head.FirstByteAfterHeader) + databytes_preceeding_incs;
cout << "Starting offset " << offset << endl;

			#pragma omp master
			{
				MPI_File ioReadFileHdl;
				MPI_Status ioReadFileStatus;
				//##MK::MPI_COMM_SELF with very many processes on the same file may run in contention...
				string fn = Settings::SpectralOutFilename + ".spectralOut";
				MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
				//skipping header and all bytes with heavy data increments prior to the one I am interested in
				MPI_File_seek(ioReadFileHdl, offset , MPI_SEEK_SET);

				MPI_Offset now;
				MPI_File_get_position(ioReadFileHdl, &now);
cout << "MPI Offset " << static_cast<long long>(now) <<  endl;
				//##ERROR handling
			}
			#pragma omp barrier

			//remember: element in this function means a paket of 8Bytes*head.matpres, i.e. all data for one integration point of one DAMASK_spectral mesh element
			//MK::hence we do always read in complete logical packages of heavy data rather than bytes!

			//all threads will evaluate these elements!
			for ( ElementsReadAlready = 0; ElementsReadAlready < ElementsTotal; ElementsReadAlready += ElementsPerBlockNow) {
				if ( (ElementsTotal - ElementsReadAlready) >= ElementsPerBlockTarget ) { //most likely
					ElementsPerBlockNow = ElementsPerBlockTarget;
					//keep reutilizing already allocated read buffer for MPI I/O
				}
				else { //remainder portion of heavy data different number of elements requires to discard old buffer and reallocate
					ElementsPerBlockNow = ElementsTotal - ElementsReadAlready;

					#pragma omp master
					{
						if ( buf != NULL ) {
							delete [] buf; buf = NULL;
						}
						nalloc = ElementsPerBlockNow * static_cast<size_t>(head.matpres);
						try {
							buf = new double[nalloc];
						}
						catch (std::bad_alloc &croak) {
							stopping("Unable to allocate MPI I/O read buffer for remainder of increment dataset", myrank, 0);
							mehealthy = false; //break;
							healthy = false;
						}
						if ( mehealthy == true ) {
							//now that we have a buffer, clear it DEBUG to detect data inconsistencies quicker
							for ( size_t i = 0; i < nalloc; ++i ) { buf[i] = numeric_limits<real>::max(); } //0.0; }//##MK::debug - clear buffer
						}

					}
					//necessary such that non-master threads do not continue with next cycle
					#pragma omp barrier
				}

				//existing working buf which all threads know to exist and be operative
				#pragma omp master
				{
					int len = static_cast<int>(nalloc);
	cout << "On Master NAlloc sizet / int" << nalloc << "--" << len << endl;
					MPI_File_read_at( ioReadFileHdl, offset, buf, len, MPI_DOUBLE, &ioReadFileStatus); //implicit advance of fp

					int count = -1;
	cout << "count " << count << endl;
					MPI_Get_count( &ioReadFileStatus, MPI_INT, &count );
	cout << "MPIStatus get count " << count << endl;
//for (size_t i = 0; i < nalloc; ++i) { cout << buf[i] << endl; }
					MPI_Offset now;
					MPI_File_get_position(ioReadFileHdl, &now);
	cout << "MPI Offset read ElementsReadAlready " << ElementsReadAlready << ";" << static_cast<long long>(now) <<  endl;
					//##MK::data have arrived
				}
				//necessary, threads have to wait for data to arrive
				#pragma omp barrier

				//data have arrived threads can grab in parallel their portion of interest from the buf
				//all threads probe whether there are useful data in buf knowing in advance their
				//integration and mesh element id range for which they are responsible
				db.at(tid)->db_distr_homogenization( buf, ElementsReadAlready, ElementsPerBlockNow );
				db.at(tid)->db_distr_crystallite( buf, ElementsReadAlready, ElementsPerBlockNow );
				db.at(tid)->db_distr_constitutive( buf, ElementsReadAlready, ElementsPerBlockNow );

				//##MK::error handling at thread level when filling database
			} //all elements processed
		} //data reading done

		#pragma omp master
		{
			if ( buf != NULL ) {
				delete [] buf; buf = NULL;
			}

			//no MPI_Barrier required because was open MPI_COMM_SELF
			MPI_File_close(&ioReadFileHdl);
		}
		#pragma omp barrier

		string s ="Finished reading heavy data of increment " + to_string(thisincrement);
		reporting( s, myrank, 0, false);
	}//end of parallel region
}
*/


void specOutIncr::specout_read_heavydata2()
{
	//Fortran spectralOut writes in double precision (KIND=8), in C this is MPI_DOUBLE
	//provided the generating system and the postprocessing system have the same endianness
	//this is only portable, as the corresponding MPI I/O routines in DAMASK_spectral when using the same system
	//for x86 Intel/AMD-based workstations and clusters should be no problem, usually they are little endian
	//for supercomputers of Cray and some Sun system this must be reconsidered by implementing
	//###MK::endianness checks

	double tic = MPI_Wtime();
	double ltic = 0.0; //local tictoc
	double ltoc = 0.0;
	double tasks[4] = {0.0, 0.0, 0.0, 0.0}; //MPI read, homog interpret, cryst interpret, const interpret

	double* buf = NULL;
	size_t nalloc = 0;

	//MPI_File ioReadFileHdl;
	//MPI_Status ioReadFileStatus;
	long long offset = 0;

	//thread local copies of data organization constructs
	//how many to read in total, an element is in the following understood as a homog+cryst+consti data package for one integration point of one mesh element
	size_t ElementsTotal = static_cast<size_t>(head.Nip) * static_cast<size_t>(head.Ncp);
	//we do not read all data at once but in multiples of the stripe size of the underlying raid system
	size_t ElementsPerBlockTarget = 1 + (Performance::MPIIOStripeSize / static_cast<size_t>(sizeof(MPI_DOUBLE)*head.matpres) );		//how many to read with current block
	size_t ElementsPerBlockNow = (ElementsTotal > ElementsPerBlockTarget) ? ElementsPerBlockTarget : ElementsTotal;
	//for how many ips do we have heavy data read already from the file
	size_t ElementsReadAlready = 0;

cout << "ElementsTotal/BlockTarget/Now " << ElementsTotal << ";" << ElementsPerBlockTarget << ";" << ElementsPerBlockNow << endl;
	nalloc = ElementsPerBlockNow * static_cast<size_t>(head.matpres);
	try { buf = new double[nalloc]; }
	catch (std::bad_alloc &croak) {
		stopping("Unable to allocate MPI I/O read buffer at first place", myrank, 0);
		healthy = false;
	}

	if (healthy == false) {
		return;
	}
	else {
		//set starting point to read heavy data for particular increment
		//MK::C-style indexing, i.e. starting at 0, instead of Fortran starting at 1
		long long databytes_preceeding_incs = static_cast<long long>(thiswrittenincrement);
		databytes_preceeding_incs *= static_cast<long long>(8); //because DAMASKspectral writes all output in double precision floats
		databytes_preceeding_incs *= static_cast<long long>(head.matpres);
		databytes_preceeding_incs *= static_cast<long long>(head.Nip);
		databytes_preceeding_incs *= static_cast<long long>(head.Ncp);

		offset = static_cast<long long>(head.FirstByteAfterHeader);
cout << "FirstByteAfterHeader " << offset << endl;

		offset += databytes_preceeding_incs;
cout << "Starting offset " << offset << endl;

		MPI_File ioReadFileHdl;
		MPI_Status ioReadFileStatus;
		//##MK::MPI_COMM_SELF with very many processes on the same file may run in contention...
		string fn = Settings::SpectralOutFilename + ".spectralOut";
		MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
		//skipping header and all bytes with heavy data increments prior to the one I am interested in
		MPI_File_seek(ioReadFileHdl, offset , MPI_SEEK_SET);

		MPI_Offset now;
		MPI_File_get_position(ioReadFileHdl, &now);
cout << "MPI Offset " << static_cast<long long>(now) <<  endl;
		//##ERROR handling


		//remember: element in this function means a paket of 8Bytes*head.matpres, i.e. all data for one integration point of one DAMASK_spectral mesh element
		//MK::hence we do always read in complete logical packages of heavy data rather than bytes!
		int cc = 0;
		for ( ElementsReadAlready = 0; ElementsReadAlready < ElementsTotal; ElementsReadAlready += ElementsPerBlockNow) {
			cc++;
cout << cc << endl;
			if ( (ElementsTotal - ElementsReadAlready) >= ElementsPerBlockTarget ) { //most likely
				ElementsPerBlockNow = ElementsPerBlockTarget;
				//keep reutilizing already allocated read buffer for MPI I/O
			}
			else { //remainder portion of heavy data different number of elements requires to discard old buffer and reallocate
				ElementsPerBlockNow = ElementsTotal - ElementsReadAlready;

				if ( buf != NULL ) {
					delete [] buf; buf = NULL;
				}
				nalloc = ElementsPerBlockNow * static_cast<size_t>(head.matpres);
				try {
					buf = new double[nalloc];
				}
				catch (std::bad_alloc &croak) {
					stopping("Unable to allocate MPI I/O read buffer for remainder of increment dataset", myrank, 0);
					healthy = false;
				}
				if ( healthy == true ) {
					//now that we have a buffer, clear it DEBUG to detect data inconsistencies quicker
					for ( size_t i = 0; i < nalloc; ++i ) { buf[i] = numeric_limits<real>::max(); } //0.0; }//##MK::debug - clear buffer
				}
			}

			//existing working buf which all threads know to exist and be operative
			int len = static_cast<int>(nalloc);

			ltic = MPI_Wtime();
//cout << "On Master NAlloc sizet / int" << nalloc << "--" << len << endl;
			MPI_File_read( ioReadFileHdl, buf, len, MPI_DOUBLE, &ioReadFileStatus); //implicit advance of fp
			ltoc = MPI_Wtime();
			tasks[0] += (ltoc-ltic);

			//int count = -1; MPI_Get_count( &ioReadFileStatus, MPI_DOUBLE, &count ); cout << "MPIStatus get count " << count << endl;
			//for (size_t i = 0; i < nalloc; ++i) { cout << buf[i] << endl; }
			MPI_Offset now; MPI_File_get_position(ioReadFileHdl, &now);
cout << "MPI Offset read ElementsReadAlready " << ElementsReadAlready << ";" << static_cast<long long>(now) <<  endl;

			//data have arrived threads can grab in parallel their portion of interest from the buf
			//all threads probe whether there are useful data in buf knowing in advance their
			//integration and mesh element id range for which they are responsible

			//##MK::BEGIN DEBUG
			//for(size_t j = 0; j < nalloc; j = j + head.matpres) {
			//	for(unsigned int k = 0; k < head.matpres; ++k)
			//		cout << buf[j+k] << "__";
			//	cout << endl;
			//}
			//##MK::END DEBUG

			unsigned int nt = db.size();
			for (unsigned int mt = MASTER; mt < nt; mt++) {
				ltic = MPI_Wtime();
				db.at(mt)->db_distr_homogenization( buf, ElementsReadAlready, ElementsPerBlockNow );
				ltoc = MPI_Wtime();
				tasks[1] += (ltoc-ltic);

				ltic = MPI_Wtime();
				db.at(mt)->db_distr_crystallite2( buf, ElementsReadAlready, ElementsPerBlockNow );
				ltoc = MPI_Wtime();
				tasks[2] += (ltoc-ltic);
//cout << tid << " success!" << endl;

				ltic = MPI_Wtime();
				db.at(mt)->db_distr_constitutive2( buf, ElementsReadAlready, ElementsPerBlockNow );
				ltoc = MPI_Wtime();
				tasks[3] += (ltoc-ltic);
			}

			//##MK::error handling at thread level when filling database
		} //all elements processed

		//no MPI_Barrier required because was open MPI_COMM_SELF
		MPI_File_close(&ioReadFileHdl);

	} //data reading done

	if ( buf != NULL ) {
		delete [] buf; buf = NULL;
	}

	string mess ="Finished reading heavy data of increment " + to_string(thisincrement);
	reporting( mess, myrank, 0, false);

	double toc = MPI_Wtime();
	mess = "IO:SpecOutReadHeavyDataMPIFileRead" + to_string(thisincrement);
	tictoc.push_back(plog(tasks[0], mess)); //subtime spent NMPIIO File read library call
	mess = "IO:SpecOutReadHeavyDataHomogInterpret" + to_string(thisincrement);
	tictoc.push_back(plog(tasks[1], mess));
	mess = "IO:SpecOutReadHeavyDataCrystInterpret" + to_string(thisincrement);
	tictoc.push_back(plog(tasks[2], mess));
	mess = "IO:SpecOutReadHeavyDataConstInterpret" + to_string(thisincrement);
	tictoc.push_back(plog(tasks[3], mess));
	mess = "IO::SpecOutReadHeavyDataOtherTasks" + to_string(thisincrement);
	tictoc.push_back(plog(((toc-tic)-tasks[0]-tasks[1]-tasks[2]-tasks[3]), mess)); //remainder prep and mem handling
}


bool specOutIncr::specout_check_heavydata2()
{
	double tic = MPI_Wtime();

	//check that length ##MK and content of data is at expected location
	bool data_are_consistent = true;

	size_t ElementsTotal = static_cast<size_t>(head.Nip) * static_cast<size_t>(head.Ncp);

	//##MK::parallelize

	unsigned int nt = db.size();
	for (unsigned int mt = MASTER; mt < nt; mt++) {
//cout << "Databse memRegion " << t << " checking consistence..." << endl;
		memRegion* thisregion = db.at(mt);

/*
//###MK::DEBUG
cout << endl;
cout << "V size = " << thisregion->crystallite.V.size() << endl;
cout << "q size = " << thisregion->crystallite.q.size() << endl;
cout << "F size = " << thisregion->crystallite.F.size() << endl;
cout << "rhoe size = " << thisregion->constitutive.rho_e.size() << endl;
cout << "rhod size = " << thisregion->constitutive.rho_d.size() << endl;
cout << endl;
*/

		if ( thisregion->crystallite.q.size() != thisregion->eipid_n ) {
//cout << "q size different" << endl;
			data_are_consistent = false;
		}
		else {
//cout << "q size consistent" << endl;
		}
		if ( thisregion->crystallite.F.size() != thisregion->eipid_n ) {
//cout << "F size different" << endl;
			data_are_consistent = false;
		}
		else {
//cout << "F size consistent" << endl;
		}
		if ( thisregion->constitutive.rho_e.size() != (thisregion->eipid_n * thisregion->constitutive.rho_e_mult ) ) {
//cout << "rho_e size different" << endl;
			data_are_consistent = false;
		}
		else {
//cout << "rho_e size the consistent" << endl;
		}
		if ( thisregion->constitutive.rho_d.size() != (thisregion->eipid_n * thisregion->constitutive.rho_d_mult ) ) {
//cout << "rho_d size different" << endl;
			data_are_consistent = false;
		}
		else {
//cout << "rho_d size consistent" << endl;
		}

//cout << "Database memRegion " << t << " consistent" << endl;
	} //check next memRegion

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::SpecOutCheckHeavyData" + to_string(thisincrement) ));

	return data_are_consistent;
}


void specOutIncr::specout_auxiliarytasks()
{
	//##MK::implement auxiliary tasks here
	double tic = MPI_Wtime();

	//##MK::highly developmental and not optimized for speed

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::AuxiliaryTasks" + to_string(thisincrement) ));
}


void specOutIncr::write_ipgrid_textureid()
{
	double tic = MPI_Wtime();

	//##MK::highly developmental and not optimized for speed

	vector<p3d> ppp;
	vector<unsigned int> texid;
	for (size_t mr = 0; mr < db.size(); ++mr) {
		memRegion* thisregion = db.at(mr);
		for(size_t i = 0; i < thisregion->grid.xyz0.size(); ++i) {
			ppp.push_back( p3d(
					(thisregion->grid.xyz0[i].x + thisregion->grid.dxyz_avg[i].dx + thisregion->grid.dxyz_flu[i].dx),
					(thisregion->grid.xyz0[i].y + thisregion->grid.dxyz_avg[i].dy + thisregion->grid.dxyz_flu[i].dy),
					(thisregion->grid.xyz0[i].z + thisregion->grid.dxyz_avg[i].dz + thisregion->grid.dxyz_flu[i].dz) ) );
		}

		for(size_t i = 0; i < thisregion->crystallite.TextureID.size(); ++i) {
			texid.push_back( thisregion->crystallite.TextureID[i] );
		}
	}

	bool status = vtk_p3dm3( ppp, texid, thisincrement );

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::AuxiliaryTasks" + to_string(thisincrement) ));
}


void specOutIncr::bounded_volume_hierarchy()
{
	//thread parallel processing for all elements initial+displacements +- periodic images in inspection radius
	//axis-aligned bounding box about points
	//3d binning of aabb into bvh_xyzm2
	double tic = MPI_Wtime();

	pp3rve1withguard.build_bvh_xyzm2();

	//bvh_xyzm2 containers with x,y,z boostSIMD aligned
	//bvh_xyzm2 containers with element id
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::BVHXYZM2" + to_string(thisincrement)));
}



void specOutIncr::analyze_svar_closestuip_disoriangle()
{
	//we desire to compute for all meshelements the distance to their closest neighboring element
	//(considering/including periodic images) within a spherical window of Settings::KernelRadius about each element
	//whose disorientation is larger than Settings::CriticalDisoriAngle, if any such exists
	//therefore we require

	/*this->db.at(0)->constitutive.rho_d_nextfreeslot = 1;
	this->db.at(0)->homogenization.eid_nextfreeslot = 1;
	this->db.at(0)->crystallite.F_nextfreeslot = 1;*/
	double tic = MPI_Wtime();

	vector<MPI_DisloSpatDistr_Double> globalres_edge;
	vector<MPI_DisloSpatDistr_Double> globalres_dipo;

    //sanity checks
	//##MK::currently only fcc supported!
	size_t nslipsys = 12;
	if ( db.back()->constitutive.rho_e_mult != nslipsys || db.back()->constitutive.rho_d_mult != nslipsys ) {
		reporting("Attempting to analyze non-supported crystal structure!", myrank, 0, true);
		return;
	}

	//MK::in case of running sequentially utilize globalresults directly

	#pragma omp parallel shared(globalres_edge, globalres_dipo)
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		size_t myworkload = 0;

		//MK::now the threads process their elements while reading the bounded volume hierarchy pp3rve1withguard to get the neighbors
		memRegion* thisregion = db.at(mt);
		size_t eip_s = thisregion->eipid_start;
		size_t eip_n = thisregion->eipid_n;
		size_t eid = eip_s;

		vector<dist> candidates;

		//thread-local results buffer
		vector<MPI_DisloSpatDistr_Double> myres_edge;
		vector<MPI_DisloSpatDistr_Double> myres_dipo;

		for (size_t e = 0; e < eip_n; ++e, ++eid) { //only the not periodic images of each element non-periodic images are material points to probe for boundary proximity

//if ( eid == 8421504 ) {
			p3d matpoint = p3d(
					thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx,
					thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy,
					thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz );

			//get state variable values at matpoint

			dist res = dist( Settings::KernelRadius + EPSILON, numeric_limits<unsigned int>::max() );
			//size_t incrworkload = 0;
			for (int incr = INCREMENTAL_KERNEL_RADIUS_STEPPING; incr > -1; incr--) { //set stepping to begin at incr = 0 for avoid incremental increasing
				real_xyz SearchRadius = Settings::KernelRadius / (1.0 + static_cast<real_xyz>(incr));

				pp3rve1withguard.find_higherorder_neighbors( matpoint, candidates, SearchRadius);

				res = closestBoundaryInKernelIfAny( eid, candidates );
				//MK::return nbor_eipid != uint::max if a solution was found

				//incrworkload += candidates.size();

				if ( res.nbor_eipid != numeric_limits<unsigned int>::max() ) //very likely if most in early cases
					break;
				//else inspecting larger search radius, if that also does not yield, returning maximum or eventually
				//if never breaks limiting inspection KernelRadius + EPSILON is reported as distance
			}

			//pp3rve1withguard.find_higherorder_neighbors( matpoint, candidates, Settings::KernelRadius );
			//dist res = closestBoundaryInKernelIfAny( eid, candidates, myresults_suspicious ); //MK::eid gets incremented

			myworkload += candidates.size();

			//collect results for I/O
			if ( nt > SINGLE_THREADED ) {
				myres_edge.push_back( MPI_DisloSpatDistr_Double(res.d, thisregion->constitutive.rho_e, e*nslipsys) );
				myres_dipo.push_back( MPI_DisloSpatDistr_Double(res.d, thisregion->constitutive.rho_d, e*nslipsys) );
			}
			else {
				//utilize output buffer directly as there is not thread concurrency and we save memory
				globalres_edge.push_back( MPI_DisloSpatDistr_Double(res.d, thisregion->constitutive.rho_e, e*nslipsys) );
				globalres_dipo.push_back( MPI_DisloSpatDistr_Double(res.d, thisregion->constitutive.rho_d, e*nslipsys) );
			}
//}
//cout << e << "/" << eip_n << endl;
		} //next element

		//aggregation of thread-local data as the point when the threads enter the critical is not entirely predictable we generate different order in output array
		//for debugging purposes undesirable we can improve by enforcing a strict order
		#pragma omp barrier
		for (unsigned int mr = 0; mr < nt; mr++) {
			if (mt == mr) { //applies to only one thread at a time, enforced sequentialization
				if ( nt > SINGLE_THREADED ) { //only necessary when running multithreaded
					cout << "Thread " << mt << " workload " << myworkload << endl;
					reporting("Computation of higher-order-based distances locally completed", get_myrank(), mt, true);

					//copy threadlocal results critically and individually to global buffer, unsorted...
					if ( myres_edge.size() == myres_dipo.size() ) {
						size_t ni = myres_edge.size();
						for (size_t i = 0; i < ni ; ++i) {
							globalres_edge.push_back( myres_edge[i] );
							globalres_dipo.push_back( myres_dipo[i] );
						}
					}
					else {
						complaining("Computation of higher-order-based distances resulted in dissimilar long arrays", get_myrank(), mt);
					}
				}
				//else, we are working SINGLE_THREADED then data are already in globalresults_* and myres_* remain empty
				reporting("Computation of higher-order-based distances globally committed", get_myrank(), mt, true);
			}
			//mt != mr wait
			#pragma omp barrier
		}

		/*//slightly faster because using workload difference-induced gaps productively, however resulting in non-ordered file output
		#pragma omp critical
		{
			if ( nt > SINGLE_THREADED ) { //only necessary when running multithreaded
				cout << "Thread " << mt << " workload " << myworkload << endl;
				reporting("Computation of higher-order-based distances locally completed", get_myrank(), mt, true);

				//copy threadlocal results critically and individually to global buffer, unsorted...
				if ( myres_edge.size() == myres_dipo.size() ) {
					size_t ni = myres_edge.size();
					for (size_t i = 0; i < ni ; ++i) {
						globalres_edge.push_back( myres_edge[i] );
						globalres_dipo.push_back( myres_dipo[i] );
					}
				}
				else {
					complaining("Computation of higher-order-based distances resulted in dissimilar long arrays", get_myrank(), mt);
				}
			}
			//else, we are working SINGLE_THREADED then data are already in globalresults_* and myres_* remain empty
			reporting("Computation of higher-order-based distances globally committed", get_myrank(), mt, true);
		}
		*/

	} //end of parallel region, implicit barrier

	/*
	//##MK::DEBUG report suspicious cases
	vector<p3d> ppp;
	vector<unsigned int> textureid;
	for (size_t k = 0; k < globalresults_suspicious.size(); ++k) {
		ppp.push_back( eid2p3d(globalresults_suspicious[k]) );
		textureid.push_back( eid2textureid(globalresults_suspicious[k]) );
	}
	bool status = vtk_p3dm3( ppp, textureid, numeric_limits<unsigned int>::max() );
	//##MK::END OF DEBUG
	*/
	//cout << "Total number of suspicious cases " << suspicioustotal << endl;


	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddStateVarSpatialDistr" + to_string(thisincrement)));

	//write_histograms2_ascii(globalresults_d, globalresults_edge, globalresults_dipo);
	write_histograms2_binary( globalres_edge, globalres_dipo );

	//write_searchefficiency(globalresults_ncandidates);
}


void specOutIncr::compute_edge_weights( const size_t central_eid, vector<dist> const & candidates, vector<lvwtedge> & edgs )
{
	//MK::DEBUG at the moment

	//get orientation of central point
	unsigned int ceid = static_cast<unsigned int>(central_eid);
	quat central_ori = eid2quaternion( central_eid );
	real_ori qcentral[4] = { central_ori.q0, central_ori.q1, central_ori.q2, central_ori.q3 };

	for ( size_t i = 1; i < candidates.size(); ++i ) {
		//we start at 1 because we exclude self-references as edge
		unsigned int neid = candidates[i].nbor_eipid;
		quat nbor_ori = eid2quaternion( neid );
		real_ori qnbor[4] = { nbor_ori.q0, nbor_ori.q1, nbor_ori.q2, nbor_ori.q3 };

		real_ori q0 = disorientation_q0_fcc( qcentral, qnbor );

		double weight = exp( Settings::CommDetectionDisoriAngleWght * (static_cast<double>(q0) - 1.0) );

		//##MK::also weight by node distance utilizing scaling constant Settings::CommDetectionDistanceWght

		edgs.push_back( lvwtedge(ceid, neid, weight) ); //multiples will be found vice versa

//cout << edgs.back().src << "\t" << edgs.back().dest << "\t" << edgs.back().wt << endl;
	}
}


void specOutIncr::write_edge_weights( vector<lvwtedge> const & edgs )
{
	//	ofstream out( "10g10x10x10_01_tensionX.0.edges.bin", ofstream::binary);
	//	out.write( (char*)&edges[0], edges.size() * sizeof(lvwtedge) );
	//cout << "Written to 10g10x10x10_01_tensionX.0.edges.bin" << endl;
	string fn = "LOUVAIN.DEBUG.SimID." + to_string(Settings::SimID) + ".Incr." + to_string(thisincrement) + ".txt";
	ofstream ascii;
	ascii.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (ascii.is_open() == true) { //headerless
		for (size_t i = 0; i < edgs.size(); ++i)
			ascii << edgs[i].src << "\t" << edgs[i].dest << "\t" << edgs[i].wt << endl;
		//ascii.flush();
		ascii.close();
	}
cout << "Written to " << fn << endl;
}


void specOutIncr::hierarchical_community_detection( vector<lvwtedge> const & edgs, vector<unsigned int> & uip2community )
{
	cerr << setprecision(18);

	class louvainHdl* lvseqworker = NULL;
	try {
		lvseqworker = new class louvainHdl;
	}
	catch (bad_alloc &croak) {
		string mess = "Community detection with Louvain modularity cluster method was unsuccessful";
		stopping( mess, get_myrank(), 0 ); //#MK::was this->
		return;
	}

	lvseqworker->execute( edgs, uip2community );

	lvseqworker->profiling();

	if ( lvseqworker != NULL ) {
		delete lvseqworker;
		lvseqworker = NULL;
	}
}


void specOutIncr::analyze_identify_grains()
{
	double tic, toc;
	tic = MPI_Wtime();

	//##MK::own data structure for polycrystal in the future taking care of spatial placement of grains
	vector<unsigned int> uniqueip2community;

	if (Settings::GrainReconModel == E_INITIAL_TEXTUREID) {

		for(size_t eid = 0; eid < head.NXYZ; ++eid) {
			unsigned int debug_gid = eid2textureid(eid) - 1; //MK::-1 because of Fortran to C indexing...
			uniqueip2community.push_back( debug_gid );
		}

		//##MK::BEGIN DEBUG
		//ofstream bin_debug("UIP.bin", ios::out | ios::binary);
		//bin_debug.write( (char*)&uniqueip2community[0], uniqueip2community.size()*sizeof(unsigned int));
		//bin_debug.close();
		//##MK::END DEBUG

		toc = MPI_Wtime();
		tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentTextureID" + to_string(thisincrement)));
	}
	else if (Settings::GrainReconModel == E_HIERARCHICAL_CLUSTERING) {

		//testing utilization of clustering-based techniques feeding V. Blondel, J-L Guillaume, R. Lambiotte's louvain
		//method for finding communities in large networks
		///https://sites.google.com/site/findcommunities/ 2015 sourceforge code version louvain-generic

		vector<lvwtedge> edges;
		vector<dist> candidates;

		//##MK::can also be OpenMP parallelized but output order on edges is required deterministic for Louvain, hence will involve to be sequentialized aggregation part as the *_svar_* functions
		unsigned int nt = db.size();
		for (unsigned int mr = 0; mr < nt; mr++) { //analyze the grain assignment of all ips
			memRegion* thisregion = db.at(mr);
			size_t eip_s = thisregion->eipid_start;
			size_t eip_n = thisregion->eipid_n;
			size_t eid = eip_s;

			for (size_t e = 0; e < eip_n; ++e) { //only the not periodic images of each element non-periodic images are material points to probe for grain identification
				p3d matpoint = p3d(
					thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx,
					thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy,
					thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz );

				pp3rve1withguard.find_higherorder_neighbors( matpoint, candidates, Settings::GrainReconInspectionRadius);

				compute_edge_weights( eid, candidates, edges );

				eid++;
			}
	cout << "Edge weights for memRegion " << mr << " computed" << endl;
		}

		toc = MPI_Wtime();
		tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentPrepareEdgeWeights" + to_string(thisincrement)));

cout << "Total edge relations " << edges.size() << endl;

		tic = MPI_Wtime();

		//#MK::BEGIN DEBUGGING
		//write_edge_weights( edges );
		//##MK::END DEBUGGING

		hierarchical_community_detection( edges, uniqueip2community );

		toc = MPI_Wtime();
		tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentCommunityLouvain" + to_string(thisincrement)));
	}
	else {
		string mess = "No message found to identify grains";
		stopping( mess, get_myrank(), 0);
		return;
	}

	//build grain objects and data structures from the labeled ip graph
	tic = MPI_Wtime();

	grains.build( uniqueip2community );

	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentBuildPool" + to_string(thisincrement)));


	tic = MPI_Wtime();
	for (unsigned int eid = 0; eid < uniqueip2community.size(); ++eid) {
		eid_write_gid( eid, uniqueip2community.at(eid) );
	}
	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentWriteAssgnGrainIDs" + to_string(thisincrement)));


	tic = MPI_Wtime();

	write_identified_grains();

	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddGrainIdentReportMeta" + to_string(thisincrement)));

	if ( Settings::VisIPGridWithGrainID > 0 ) {
		//##MK::BEGIN DEBUG
		//render grain structure
		vector<p3d> ppp;
		for (size_t eid = 0; eid < head.NXYZ; ++eid) {
			ppp.push_back( eid2p3d(eid) );
		}
		bool status = vtk_gp3d( ppp, uniqueip2community, thisincrement );
		//##MK::END OF DEBUG
	}

	if ( Settings::VisGrainQuaternionClouds > 0 ) {
		write_grainid_and_quaternions();
	}
}


void specOutIncr::analyze_svar_grainbased_sdf()
{
	//for all successfully identified and voxelated grains check for each voxel inside its distance to the boundary
	//via utilizing the SgnDistField2 and identify corresponding unique ip and state variable value
	double tic = MPI_Wtime();

	vector<MPI_DisloSpatDistr_Double> globalres_edge;
	vector<MPI_DisloSpatDistr_Double> globalres_dipo;
	
    //sanity checks
	//##MK::currently only fcc supported!
	size_t nslipsys = 12;
	if ( db.back()->constitutive.rho_e_mult != nslipsys || db.back()->constitutive.rho_d_mult != nslipsys ) {
		reporting("Attempting to analyze non-supported crystal structure!", myrank, 0, true);
		return;
	}
	
	unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
	for(unsigned int gr = 0; gr < ngr; ++gr) {
		grGeomHdl* thisgrain = sdf.at(gr);
		memRegion* thisregion = NULL;
		size_t e = numeric_limits<size_t>::max();

		if (thisgrain != NULL) { //exists
			if (thisgrain->healthy == true) {
				size_t nv = thisgrain->localgrid.nxyz;
				unsigned int inside = thisgrain->cgid;
				for(size_t v = 0; v < nv; ++v) {
					if ( thisgrain->GrainIDField[v] == inside ) { //now we read from only 4 arrays all of which are linear and dont have to search for neighbors

						//read distance of voxel center to boundary from distance function
						real_xyz dist = static_cast<real_xyz>(thisgrain->SgnDistField2[v]);
						
						//read corresponding reference to DAMASK_spectral integration point to get material point constitutive data
						size_t uipid = static_cast<size_t>(thisgrain->UIPIDField[v]);

						//find memory region which stores the state variable values
						for(unsigned int mr = 0; mr < db.size(); mr++) {
							//##MK::further improvement potential is seen here, as the searching for the correct memory location of heavy data sweeps through a vector of pointer, which individually may result in a pointer chase through memory until finding the pointer pointing to the correct heavy data container... instead utilize temporary map of
							if( uipid < db.at(mr)->eipid_start || uipid >= db.at(mr)->eipid_end ) {
								continue;
							} //not continued, hence uipid \in [eipid_start,eipid_end)
							thisregion = db.at(mr);
							e = uipid - thisregion->eipid_start;
							break;
						}

						//found!,  now pull state variable values
						globalres_edge.push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_e, e*nslipsys) );
						globalres_dipo.push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_d, e*nslipsys) );
					}
				} //check all voxel of thisgrain.localgrid whose center is inside the grain
			}
		} //next grain multiple thread work on different grains simultaneously
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddStateVarSpatialDistr" + to_string(thisincrement)));

	write_histograms2_binary( globalres_edge, globalres_dipo );
}


/*
##MK::not working version
void specOutIncr::analyze_svar_grainbased_sdf()
{
	//for all successfully identified and voxelated grains check for each voxel inside its distance to the boundary
	//via utilizing the SgnDistField2 and identify corresponding unique ip and state variable value
	double tic = MPI_Wtime();

	vector<MPI_DisloSpatDistr_Double> globalres_edge;
	vector<MPI_DisloSpatDistr_Double> globalres_dipo;
	size_t required_cnt = 0;
	
    //sanity checks
	//##MK::currently only fcc supported!
	size_t nslipsys = 12;
	if ( db.back()->constitutive.rho_e_mult != nslipsys || db.back()->constitutive.rho_d_mult != nslipsys ) {
		reporting("Attempting to analyze non-supported crystal structure!", myrank, 0, true);
		return;
	}

	//MK::in case of running sequentially utilize globalresults directly

	#pragma omp parallel shared(globalres_edge, globalres_dipo,required_cnt)
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		//thread-local results collection in case of multi-threaded execution
		//but for the purpose of debugging we want to assure that the output order is deterministic
		vector<vector<MPI_DisloSpatDistr_Double>*> myres_edge;
		vector<vector<MPI_DisloSpatDistr_Double>*> myres_dipo;
		size_t myres_cnt = 0;

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; ++gr) {
			if ( check_if_myworkpackage(gr, mt, nt) == false ) {
				continue; //as somebody else will take care
			}

			//generate new myres buffer
			//##MK::allocate checks
			vector<MPI_DisloSpatDistr_Double>* here_e = NULL;
			vector<MPI_DisloSpatDistr_Double>* here_d = NULL;
			if ( nt > SINGLE_THREADED) {
				here_e = new vector<MPI_DisloSpatDistr_Double>;
				myres_edge.push_back( here_e );
				here_d = new vector<MPI_DisloSpatDistr_Double>;
				myres_dipo.push_back( here_d );
			}
			
			grGeomHdl* thisgrain = sdf.at(gr);
			memRegion* thisregion = NULL;
			size_t e = numeric_limits<size_t>::max();

			if (thisgrain != NULL) { //exists
				if (thisgrain->healthy == true) {
					size_t nv = thisgrain->localgrid.nxyz;
					unsigned int inside = thisgrain->cgid;
					for(size_t v = 0; v < nv; ++v) {
						if( thisgrain->GrainIDField[v] == inside ) { //now we read from only 4 arrays all of which are linear and dont have to search for neighbors

							//read distance of voxel center to boundary from distance function
							real_xyz dist = static_cast<real_xyz>(thisgrain->SgnDistField2[v]);

							//read corresponding reference to DAMASK_spectral integration point to get material point constitutive data
							size_t uipid = static_cast<size_t>(thisgrain->UIPIDField[v]);

							//find memory region which stores the state variable values
							for(unsigned int mr = 0; mr < db.size(); mr++) {
								//##MK::further improvement potential is seen here, as the searching for the correct memory location of heavy data sweeps through a vector of pointer, which individually may result in a pointer chase through memory until finding the pointer pointing to the correct heavy data container... instead utilize temporary map of
								if( uipid < db.at(mr)->eipid_start || uipid >= db.at(mr)->eipid_end ) {
									continue;
								} //not continued, hence uipid \in [eipid_start,eipid_end)
								thisregion = db.at(mr);
								e = uipid - thisregion->eipid_start;
								break;
							}

							//found!,  now pull state variable values
							if ( nt > SINGLE_THREADED ) {
								myres_edge.back()->push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_e, e*nslipsys) );
								myres_dipo.back()->push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_d, e*nslipsys) );
								myres_cnt++;
							}
							else { //SINGLE_THREADED, in this case is guaranteed because of for loop
								//utilize output buffer directly as there is not thread concurrency and we save memory
								globalres_edge.push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_e, e*nslipsys) );
								globalres_dipo.push_back( MPI_DisloSpatDistr_Double(dist, thisregion->constitutive.rho_d, e*nslipsys) );
							}
						}
					} //check all voxel of thisgrain.localgrid whose center is inside the grain
				}
			}
		} //next grain multiple thread work on different grains simultaneously

		//aggregation of thread-local data as the point when the threads enter the critical is not entirely predictable we generate different order in output array
		//for debugging purposes undesirable we can improve by enforcing a strict order
		
		//make output globalres large enough
		#pragma omp critical
		{
				required_cnt += myres_cnt;
		}
		
		//##MK::explicit flush of required_cnt
		//all have to wait until total output size is known an can be reserved by master thread
		#pragma omp barrier
		
		#pragma omp master
		{
			if ( nt > SINGLE_THREADED ) {
				globalres_edge.reserve(required_cnt);
				globalres_dipo.reserve(required_cnt);
				for(size_t i = 0; i < required_cnt; ++i) {
					globalres_edge.push_back( MPI_DisloSpatDistr_Double() );
					globalres_dipo.push_back( MPI_DisloSpatDistr_Double() );					
				}
			} //as positions are now a there we can update in parallel
		}
		
		//necessary, omp master construct has none 
		#pragma omp barrier
		
		for(unsigned int gr = 0; gr < ngr; ++gr) {
			if ( check_if_myworkpackage(gr, mt, nt) == false ) {
				continue; //as somebody else will take care
			}
			for(size_t i = 0; i < )

		}
		
		for (unsigned int mr = 0; mr < nt; mr++) {
			if (mt == mr) { //applies to only one thread at a time, enforced sequentialization
				if ( nt > SINGLE_THREADED ) { //only necessary when running multithreaded
					reporting("Computation of higher-order-based distances locally completed", get_myrank(), mt, true);

					//copy threadlocal results critically and individually to global buffer, unsorted...
					if ( myres_edge.size() == myres_dipo.size() ) {
						size_t ni = myres_edge.size();
						for (size_t i = 0; i < ni ; ++i) {
							globalres_edge.push_back( myres_edge[i] );
							globalres_dipo.push_back( myres_dipo[i] );
						}
					}
					else {
						complaining("Computation of higher-order-based distances resulted in dissimilar long arrays", get_myrank(), mt);
					}
				}
				//else, we are working SINGLE_THREADED then data are already in globalresults_* and myres_* remain empty
				reporting("Computation of higher-order-based distances globally committed", get_myrank(), mt, true);
			}
			//mt != mr wait
			#pragma omp barrier
		}


//		//aggregation of thread-local data is writing to shared variable, i.e. critical
//		#pragma omp critical
//		{
//			if ( nt != SINGLE_THREADED ) { //only necessary when running multithreaded
//				reporting("Computation of higher-order-based distances locally completed", get_myrank(), mt, true);
//
//				//copy threadlocal results critically and individually to global buffer, unsorted...
//				//size_t ni = myresults_d.size();
//				if ( myres_edge.size() == myres_dipo.size() ) {
//					size_t ni = myres_edge.size();
//					for (size_t i = 0; i < ni ; ++i) {
//						globalres_edge.push_back( myres_edge[i] );
//						globalres_dipo.push_back( myres_dipo[i] );
//					}
//				}
//				else {
//					complaining("Computation of higher-order-based distances resulted in dissimilar long arrays", get_myrank(), mt);
//				}
//			}
//			//otherwise data already in globalresults_* because single-threaded execution in quasi parallel OpenMP region
//			reporting("Computation of higher-order-based distances globally committed", get_myrank(), mt, true);
//		} //end of critical region


	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddStateVarSpatialDistr" + to_string(thisincrement)));

	//write_histograms2_ascii(globalresults_d, globalresults_edge, globalresults_dipo);
	write_histograms2_binary( globalres_edge, globalres_dipo );

	//write_searchefficiency(globalresults_ncandidates);
}
*/


bool specOutIncr::check_if_myworkpackage(unsigned int const i, unsigned int const tid, unsigned int const ntid)
{
	if ( i % ntid == tid )
		return true;
	else
		return false;
	
	/*if ( ntid > 1 ) { //only when more than one thread partition
		if ( i % ntid == tid )
			return true;
		else
			return false;
	} //else, tid is alone so of course you work on it...
	return true;*/
}


bool specOutIncr::init_threadlocalmemory_sdf()
{
	double tic = MPI_Wtime();

	//build for each member of grains a mesher and signed-distance computing instance stored in sdf
	//instantiate in thread local memory
	bool allhealthy = true;
	
	//init specific place
	sdf.clear();
	unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
	for(unsigned int gr = 0; gr < ngr; ++gr) {
		sdf.push_back(NULL);
	}
		
	
	#pragma omp parallel shared(allhealthy)
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		bool mehealthy = true;

		//##MK::assign work to threads round robin assuring that pointing to threadlocal memory
		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; ++gr) {
			//##MK::many critical regions for setting up the grain local analyzer
			//however, this will very likely pay off because threads can work in threadlocal memory
			//provided when utilizing specific allocator implementations as for instance jemalloc
			//##MK::future optimization hint, use buffer strategy, generate independently local buffer pointer first and allocate than sync
			//MK::formal overhead of this for loop should in the future better be replaced by worksharing construct
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue; //as somebody else will take care

			//round robin mine
			grGeomHdl* geomhdl = NULL;
			try { geomhdl = new grGeomHdl; }
			catch (bad_alloc &mecroak) {
				mehealthy = false;
			}
			
			geomhdl->owner = this;
			geomhdl->cgid = gr;
			geomhdl->nuip = grains.polyxx.at(gr).np;
			geomhdl->healthy = true;
				
			#pragma omp critical
			{
				//writing to globally shared array, ''MK::better in the future build threadlocal array of grgeomHdl objects and pass specOutIncrHdl only references to these
cout << "Thread " << mt << " processes grain " << gr << endl;
				sdf.at(gr) = geomhdl;
				
				/*sdf.push_back(NULL);
				sdf.back() = geomhdl;
				sdf.back()->owner = this; //MK::also here this is necessary to get address of calling higher-level class object
				sdf.back()->cgid = gr;
				sdf.back()->nuip = grains.polyxx.at(gr).np;
				sdf.back()->healthy = true;*/
			}
		} //all grains processed

		#pragma omp critical
		{
			if (mehealthy == false)
				allhealthy = false;
		}
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::InitThreadLocalGrainHdl" + to_string(thisincrement)));

	if (allhealthy == false) {
		string mess = "Unable to allocate memory for grain geometry reconstruction";
		stopping(mess, get_myrank(), 0);
		return false;
	}
	string mess = "Thread-local allocation of grain geometry reconstructor successful";
	reporting( mess, get_myrank(), 0, false );
	return true;
}


void specOutIncr::compute_perips_and_uips()
{
	double tic = MPI_Wtime();

	aabb3d globalbox = aabb3d(); //MK::global, i.e. not only one grain or only uip but all ips in the rve27 ip point cloud
	#pragma omp parallel shared(globalbox)
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		aabb3d mbox = aabb3d(); //MK::local to the thread not a grain, used to aggregate threadlocal result

cout << mbox << endl;
		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue; //as somebody else will take care

			//get first all uniqueips that build the grain
			vector<unsigned int>* passthese = grains.ipsupport.at(gr);
			for(size_t i = 0; i < passthese->size(); ++i ) {
				//get spatial data of uniqueip uip and assign name of the uip to mark m1
				unsigned int uipname = passthese->at(i);
				p3d p = eid2p3d(static_cast<size_t>(uipname)); //##MK::DEBUG before mark = loop iterator
				sdf.at(gr)->ips.push_back( p3dm1(p.x, p.y, p.z, uipname) );
			}

			//get periodic images of these uips in RVE27
			sdf.at(gr)->compute_all_replica();

			//threadlocal aggregation of bounding boxes
			mbox.potentially_expand( sdf.at(gr)->localbox );
//cout << mbox << endl;

		} //me done with all grains

		//##MK::strictly speaking there should be no barrier necessary as updating to globalbox is in critical region

		//...to identify globalbox by reducing the grain-based already thread-locally aggregated localboxes via mbox is because updating global data on my_rank
		#pragma omp critical
		{
			globalbox.potentially_expand(mbox);
cout << "Thread " << mt << " local bounding box " << mbox << endl;
		}
	} //end of parallel region

	//finally get dimensions of joint AABB
	globalbox.scale();

	pp3rve27.owin = globalbox;

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingComputePerIPs" + to_string(thisincrement)));
}


void specOutIncr::init_global_bvh_p3dm1()
{
	//knowing now the final axis-aligned extent of the ip cloud globally we define a BVH globally
	//whose dimensions and spatial partitioning we choose the same for processing each grain
	//for this task to conduct individually though only the ips of a particular grain are required therefore
	//we store a sparse representation of the bvh and once finished with each grain fuse the ips from
	//this sparse structure into the global BVH that we initialize in this function

	double tic = MPI_Wtime();

	pp3rve27.build_bvh_p3dm1(); //#MK::was this->

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::BVHP3DM1" + to_string(thisincrement)));
}


void specOutIncr::build_grainlocal_sbvh()
{
	//EXECUTED FROM WITHIN PARALLEL REGION

	//knowing now the final axis-aligned extent of the ip cloud globally and for each grain
	//we now bin into a sparsely stored bounded volume hierarchy the ips of each grain to prepare
	//for performing the thread/grain-local DBScan to identify the periodic images/replicates of each grain
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			//##MK::formal overhead of this for loop should in the future better be replaced by worksharing construct
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;

			if (sdf.at(gr) != NULL) {

				sdf.at(gr)->build_sparse_bvh();


			}
		} //me done with all grains
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingBuildLocalSparseBVHs" + to_string(thisincrement)));
}


void specOutIncr::discern_replica_via_dbscan()
{
	//now we want to identify all periodic replica in RVE27 via DBScan
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			//##MK::formal overhead of this for loop should in the future better be replaced by worksharing construct
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;
			if (sdf.at(gr) != NULL) {
				sdf.at(gr)->dbscan_ips2grreplicates( Settings::DBScanKernelRadius );
			}
		} //me done with all grains
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingDBScanIdentReplicates" + to_string(thisincrement)));
}


void specOutIncr::pick_one_replica_from_dbscan()
{
	//knowing now the cluster pick one representative of the grain closest to center of uip point cloud center
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;

			sdf.at(gr)->dbscan_pick_representative();

			if (sdf.at(gr)->healthy != true) {
				#pragma omp critical
				{
					string mess = "No periodic replica found for grain " + to_string(gr);
					complaining(mess, get_myrank(), mt);
				}
			}
		} //me done with all grains
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingDBScanPickRepresentative" + to_string(thisincrement)));

	//##MK::make optional
	write_dbscan_result();
}


void specOutIncr::write_dbscan_result()
{
	double tic = MPI_Wtime();

	string fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".DBScanResult.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "GrainID;NumberUniqueIPs;LocalClusterID;Chosen;NumberOfIPs;AABBCenterX;AABBCenterY;AABBCenterZ;XMin;XMax;YMin;YMax;ZMin;ZMax\n";
		csvlog << "1;1;1;boolean;1;1;1;1;1;1;1;1;1;1\n";
		csvlog << "GrainID;NumberUniqueIPs;LocalClusterID;Chosen;NumberOfIPs;AABBCenterX;AABBCenterY;AABBCenterZ;XMin;XMax;YMin;YMax;ZMin;ZMax\n";

		for (size_t gr = 0; gr < sdf.size(); ++gr) {
			grGeomHdl* thisgrain = sdf.at(gr);
			if (thisgrain != NULL) { //exists
				//##MK::show DBScan results for all grains
				//if (thisgrain->healthy == true) {
					for(size_t cl = 0; cl < thisgrain->dbscanres.size(); ++cl) {
						aabb3d bx = thisgrain->dbscanres.at(cl).box;
						p3d ctr = bx.center();
						csvlog << thisgrain->cgid << ";" << thisgrain->nuip << ";" << thisgrain->dbscanres.at(cl).lid << ";";
						csvlog << static_cast<unsigned int>(thisgrain->dbscanres.at(cl).chosen) << ";" << thisgrain->dbscanres.at(cl).cnt << ";";
						csvlog << ctr.x << ";" << ctr.y << ";" << ctr.z << ";";
						csvlog << bx.xmi << ";" << bx.xmx << ";" << bx.ymi << ";" << bx.ymx << ";" << bx.zmi << ";" << bx.zmx << "\n";
					}
				//}
			}
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write grain-resolved DBScanResults to file", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteDBScanResults" + to_string(thisincrement)));
}


void specOutIncr::init_global_voxelgrid_csys()
{
	double tic = MPI_Wtime();

	//identify AABB about all sdf[]->thegrains considered
	//MK::these AABB may when adding a guard zone protrude out of the Moore periodic replica
	//##MK::as of now we ignore this problem during voxeling the grains assigning in such case dummy ID uint32::max()

	//the idea of the voxelization is as follows: when voxelizing naively the entire RVE27 ip cloud a considerable number of
	//voxel may be required, given the fact that reconstructed replica require voxelizing also substantial space
	//inside the Moore periodic images it is more efficient to define a global voxel grid but voxelize/ppulate it sparsely
	//i.e. given a voxelgrid definition which encloses the entire space occupied by the Moore replica ips owin and its origin
	//we can utilize this global coordinate system of voxel and sample only those voxel laying inside a guarded AABB about each successfully identified grain
	//the approach is inspired by C. Mie{\ss}en et al. grain box patchwork approach presented in
	//C. Mie{\ss}en, Liesenjohann, Barrales-Mora, Gottstein, Shvindlerman, Acta Mat 2015

	//first find for all healthy thegrains their AABB, expand these by PVTessCubeVoxelEdgeLen*PVTessGuardZoneEdgeLen

	//grainfences initially are build from bounding boxes about ip locations now we have to distinguish center of gravity of voxel and wall of cubic voxels...
	real_xyz guard_edgelen = static_cast<real_xyz>(Settings::PVTessGuardZoneEdgeLen)*Settings::PVTessCubeVoxelEdgeLen;

cout << "Guardlen " << guard_edgelen << endl;

	pp3rve27.owin.scale(); //#MK::was this->
	aabb3d rve27 = pp3rve27.owin;
cout << "RVE27 " << rve27 << endl;

	aabb3d frame_grboxes_withguard = aabb3d();

	for(size_t gr = 0; gr < sdf.size(); ++gr) {
		if ( sdf.at(gr) != NULL ) {

			aabb3d old_grbox = sdf.at(gr)->grainfence;
	cout << "OldGrBox " << endl << old_grbox << endl;

			sdf.at(gr)->grainfence.add_guardzone(guard_edgelen);

			aabb3d new_grbox = sdf.at(gr)->grainfence;
	cout << "NewGrBox " << endl << new_grbox << endl;

			if ( new_grbox.lurking_out_of( rve27 ) == true ) {
	cout << "Grain " << gr << " protrudes out of rve27" << endl;
				sdf.at(gr)->vxlinit_success = false;
			}
			else { //inside rve27 but updateframe_grboxes_withguard
				frame_grboxes_withguard.potentially_expand( new_grbox );
				sdf.at(gr)->vxlinit_success = true;
	cout << "Updating frame_grboxes_withguard... with grain " << gr << endl << frame_grboxes_withguard << endl;
			}
		}
	}

	//blowup frame_grboxes_withguard such that it is integer multiple of Settings::PVTessCubeVoxelEdgeLen
	//blowup >= 1.0 where the numerical case of exactly 1.0 is very unlikely
	frame_grboxes_withguard.scale();

cout << "Frame_grboxes_withguard " << endl << frame_grboxes_withguard << endl;

	vxlgrd tmp = vxlgrd();
	tmp.dcell = Settings::PVTessCubeVoxelEdgeLen;
	tmp.nx = static_cast<size_t>(frame_grboxes_withguard.xsz / tmp.dcell) + 1;
	tmp.ny = static_cast<size_t>(frame_grboxes_withguard.ysz / tmp.dcell) + 1;
	tmp.nz = static_cast<size_t>(frame_grboxes_withguard.zsz / tmp.dcell) + 1;
	d3d blowup(
			(static_cast<real_xyz>(tmp.nx) * tmp.dcell) / frame_grboxes_withguard.xsz,
			(static_cast<real_xyz>(tmp.ny) * tmp.dcell) / frame_grboxes_withguard.ysz,
			(static_cast<real_xyz>(tmp.nz) * tmp.dcell) / frame_grboxes_withguard.zsz );
cout << "Blowup " << blowup << endl;

//	p3d debug_old = frame_grboxes_withguard.center(); //##MK::DEBUG

	frame_grboxes_withguard.blowup_xyz( blowup );

cout << "Blown up frame_grboxes_withguard " << endl << frame_grboxes_withguard << endl;

/*
p3d debug_new = frame_grboxes_withguard.center(); //##MK::DEBUG
cout << setprecision(32);
cout << "Frame_grboxes_withguard old center " << debug_old << endl;
cout << "Frame_grboxes_withguard new center " << debug_new << endl;
cout << setprecision(18);
*/

	//define missing values of vxlgrd
	tmp.nxy = tmp.nx * tmp.ny;
	tmp.nxyz = tmp.nx * tmp.ny * tmp.nz;
	tmp.origin_cntum = p3d( frame_grboxes_withguard.xmi,
									frame_grboxes_withguard.ymi,
										frame_grboxes_withguard.zmi ); //MK::not expected to be 0.0, 0.0, 0.0 ! can be negative, somewhere also potentially outside RVE27 !
	//MK::frame_grboxes_withguard is integer multiple of cell edge length, hence its imi locate where the domain walls of the 3d cubic voxel container is not the center of gravity locations of the 0-th voxel layer!
	tmp.origin_discr_own = vxl( 0, 0, 0);
	tmp.origin_discr_lnk = vxl( 0, 0, 0); //MK::there is no linked higher level voxelgrid for this one because it is the global vxlgrd

	tmp.box = frame_grboxes_withguard;
	tmp.box.scale();

	thegrid = tmp;

cout << "Final global voxelgrid is " << thegrid << endl;

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingInitGlobalVoxelGrid" + to_string(thisincrement)));
}


void specOutIncr::init_grainlocal_voxelgrids_csys()
{
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;
			if (sdf.at(gr) != NULL) {
				if ( sdf.at(gr)->healthy == true ) {
					sdf.at(gr)->build_localgrid();
				}
			}
		}
	} //next grain

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingInitLocalVoxelGrids" + to_string(thisincrement)));
}


void specOutIncr::write_grainlocal_vxlgrids()
{
	double tic = MPI_Wtime();

	string fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".GrainVxlGrids.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "GrainID;UIPCnt;VxlInside;EquivSphericalRadius;VxlInside/UIPCnt;ReconstructedYesNo;NX;NY;NZ;NXY;NXYZ;OrgDiscrX;OrgDiscrY;OrgDiscrZ;dcell;OrgCntX;OrgCntY;OrgCntZ;BoxXMin;BoxXMax;BoxYMin;BoxYMax;BoxZMin;BoxZMax\n";
		csvlog << "1;1;1;1;1;boolean;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1\n";
		csvlog << "GrainID;UIPCnt;VxlInside;EquivSphericalRadius;VxlInside/UIPCnt;ReconstructedYesNo;NX;NY;NZ;NXY;NXYZ;OrgDiscrX;OrgDiscrY;OrgDiscrZ;dcell;OrgCntX;OrgCntY;OrgCntZ;BoxXMin;BoxXMax;BoxYMin;BoxYMax;BoxZMin;BoxZMax\n";

		//add global grid first
		vxlgrd* v = &thegrid;
		csvlog << numeric_limits<unsigned int>::max() << ";0;0;0;0.0;0.0;";
		csvlog << v->nx << ";" << v->ny << ";" << v->nz << ";" << v->nxy << ";" << v->nxyz << ";";
		csvlog << v->origin_discr_lnk.x << ";" << v->origin_discr_lnk.y << ";" << v->origin_discr_lnk.z << ";" << v->dcell << ";";
		csvlog << v->origin_cntum.x << ";" << v->origin_cntum.y << ";" << v->origin_cntum.z << ";";
		csvlog << v->box.xmi << ";" << v->box.xmx << ";" << v->box.ymi << ";" << v->box.ymx << ";" << v->box.zmi << ";" << v->box.zmx << "\n";

		double _fgeo = static_cast<double>(3.0) / (static_cast<double>(4.0) * static_cast<double>(PI));
		double _onethird = static_cast<double>(1.0) / static_cast<double>(3.0);

		for (size_t gr = 0; gr < sdf.size(); ++gr) {
			grGeomHdl* thisgrain = sdf.at(gr);
			if (thisgrain != NULL) { //exists
				if (thisgrain->healthy == true) {
					v = &(thisgrain->localgrid);
					double vxl_vs_uip = static_cast<double>(thisgrain->nuvxl) / static_cast<double>(thisgrain->nuip);
					double ehr = pow( static_cast<double>(thisgrain->nuip) * _fgeo, _onethird );
					csvlog << thisgrain->cgid << ";" << thisgrain->nuip << ";" << thisgrain->nuvxl << ";" << ehr << ";" << vxl_vs_uip << ";" << 1 << ";";
					csvlog << v->nx << ";" << v->ny << ";" << v->nz << ";" << v->nxy << ";" << v->nxyz << ";";
					csvlog << v->origin_discr_lnk.x << ";" << v->origin_discr_lnk.y << ";" << v->origin_discr_lnk.z << ";" << v->dcell << ";";
					csvlog << v->origin_cntum.x << ";" << v->origin_cntum.y << ";" << v->origin_cntum.z << ";";
					csvlog << v->box.xmi << ";" << v->box.xmx << ";" << v->box.ymi << ";" << v->box.ymx << ";" << v->box.zmi << ";" << v->box.zmx << "\n";
				}
			}
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write grain-resolved voxel grid meta data to file", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteGrainLocalVoxelGrids" + to_string(thisincrement)));
}


void specOutIncr::fill_global_bvh_p3dm1()
{
	//MK::make use of the grain-resolved sparse bvh structures and fuse them into a global BVH, deleting the grainlocal

	//why is at all now again a global bvh required?, because as soon as
	//we voxelize we seek to find the ip closest to the center of gravity of the voxel to find the grain
	//this may be the ip supporting a particular grain for which a sparse bvh already exists but it may not,
	//we know that each uip and therefore also perip has a grainID assigned and the RVE simulation was periodic BC
	//therefore

	//##MK::future optimization hint
	//store p3dm1 objects in such a way that the pointer points to the memory of that thread which has to
	//voxelize most grains in the sptial region in which these containers lay

	double tic = MPI_Wtime();

	//##MK::do not run in parallel to avoid write concurrency, tasks seems at first glance of inherent sequential nature/requirement ?

	unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
	for(unsigned int gr = 0; gr < ngr; gr++) {
		//unsigned int grainID = sdf.at(gr)->cgid;

//cout << "Grain " << gr << " filling in ips" << endl;
		map<unsigned int, sbvhrange>::iterator it = sdf.at(gr)->sbvh_locator.begin();
		//vector<p3d>* readhere = sdf.at(gr);
		while( it != sdf.at(gr)->sbvh_locator.end() ) {

			//MK::observe that all elements on sdf.at(gr)->sbvh on the index interval [idxs,idxe)
			//go to the same bin and were assigned the same grainID ...
			unsigned int binID = it->first;
//cout << binID << " ";
			size_t idxs = it->second.startidx;
			size_t idxe = it->second.pastendidx;

			vector<p3dm2>* writehere = pp3rve27.pp3_points.at(binID);

			for(size_t i = idxs; i < idxe; ++i) {
				writehere->push_back( p3dm2(sdf.at(gr)->sbvh_points[i], sdf.at(gr)->sbvh_uipref[i], sdf.at(gr)->sbvh_isrepresentative[i]) );
			}

			it++; //enter next binID
		} //all sparse data transferred


		//##MK::swap free the sbvh_buckets as they are no longer required
		vector<p3d> dummy1;
		vector<unsigned int> dummy2;
		vector<unsigned char> dummy3;
		sdf.at(gr)->sbvh_points.swap( dummy1 );
		sdf.at(gr)->sbvh_uipref.swap( dummy2 );
		sdf.at(gr)->sbvh_isrepresentative.swap( dummy3 );

//cout << endl << endl;
	} //next grain

	//##MK::BEGIN DEBUG
	if( Settings::VisRVE27BVH == 1 ) {
		string mess = "RVE27BoundedVolHrarchy";
		string descr = "RVE27 BVH Points and BinIDs";
		bool status = vtk_bvh_p3dm1( pp3rve27.pp3_points, thisincrement, mess, descr );
	}
	//##MK::END DEBUG


	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingFillGlobalBVH" + to_string(thisincrement)));
}


void specOutIncr::voxelize_this_replica()
{
	//knowing now the cluster pick one representative of the grain closest to center of uip point cloud center
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;
			if (sdf.at(gr) != NULL) {
				if (sdf.at(gr)->healthy == true) {
					sdf.at(gr)->voxelize_via_pvtessellation();
				}
			}
		} //me done with all grains
	} //end of parallel region


	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::MeshingVoxelizingRepresentative" + to_string(thisincrement)));
}


void specOutIncr::approximate_signed_distance_function()
{
	//having now the cubic voxel grid representation of the grain we can (finally) compute the signed-distance function
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;
			if(sdf.at(gr) != NULL) {
				if (sdf.at(gr)->healthy == true) {
					sdf.at(gr)->compute_sgndistfun_coarse();
				}
			}
		} //me done with all grains
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::SgnDistFunctionInit" + to_string(thisincrement)));
}


void specOutIncr::spread_signed_distance_function()
{
	//having now the cubic voxel grid representation of the grain we can (finally) compute the signed-distance function
	double tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		unsigned int ngr = static_cast<unsigned int>(grains.polyxx.size());
		for(unsigned int gr = 0; gr < ngr; gr++) {
			if ( check_if_myworkpackage(gr, mt, nt) == false )
				continue;
			if (sdf.at(gr) != NULL ) {
				if (sdf.at(gr)->healthy == true) {
					sdf.at(gr)->compute_sgndistfun_fsm();
				}
			}
		} //me done with all grains
	} //end of parallel region

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::SgnDistFunctionFSM" + to_string(thisincrement)));

	//##MK::make optional
	write_grainlocal_vxlgrids();
}


void specOutIncr::analyze_mesh_grains()
{
	if ( init_threadlocalmemory_sdf() != true )
		return;

	compute_perips_and_uips();

	init_global_bvh_p3dm1();

	build_grainlocal_sbvh();

	//MK::from now on we do not necessarily require storage of the grGeomHdl::ips, MK::nw they are swapped free

	discern_replica_via_dbscan();

	pick_one_replica_from_dbscan();

	//MK::from now on we must check if the grain periodic image recon was successful, i.e. if sdf.at(i)->healthy == true

	init_global_voxelgrid_csys();

	init_grainlocal_voxelgrids_csys();

	//MK::even if the dbscan entirely fails for a grain, at least the periodic images of all uips to each grain were computed in
	//compute_perips_and_uips to support the PV voxelization

	fill_global_bvh_p3dm1();

	//MK::from now on we may no longer require the grain-local sbvh's

	voxelize_this_replica();

	approximate_signed_distance_function();

	spread_signed_distance_function();

	//##MK::now one would be ready to utilize the signed distance function, where it is defined to get
	//the project normal distance of a uip in the deformed configuration to the grain boundary
}


/*
void specOutIncr::analyze_reconstruct_grains()
{
	double tic = MPI_Wtime();

	//build field carrying the successively build assignment of grain ids to ips
	size_t n = head.NXYZ; //0;
cout << "Reconstructing grains with n " << n << endl;

	vector<unsigned int> pid2gid_tmp;
	for (size_t i = 0; i < n; ++i)
		pid2gid_tmp.push_back( numeric_limits<unsigned int>::max() );

	vector<grain> polyxx_tmp;
	vector<dist> candidates;

	for (unsigned int mr = 0; mr < db.size(); mr++) { //analyze the grain assignment of all ips
		memRegion* thisregion = db.at(mr);
		size_t eip_s = thisregion->eipid_start;
		size_t eip_n = thisregion->eipid_n;
		size_t eid = eip_s;

		for (size_t e = 0; e < eip_n; ++e) { //only the not periodic images of each element non-periodic images are material points to probe for grain identification
			p3d matpoint = p3d(
				thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx,
				thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy,
				thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz );

			pp3rve1withguard.find_higherorder_neighbors( matpoint, candidates, Settings::GrainReconInspectionRadius);

			//the candidates list is sorted by ascending distances to the location of matpoint

			//check if any of the candidates has a grain already assigned and if so the orientation of such grain is
			//similar to the one of the matpoint within SO3 according to LocalDisoriAngle
			//if it is return gid of that most similar grain in terms of orientation proximity
			//if not generate a new grain and return its gid
			//pid2gid_tmp[eid] = anyGrainAlreadyAssigned( eid, candidates, pid2gid_tmp, polyxx_tmp );
			pid2gid_tmp[eid] = anyGrainAlreadyExistent1( eid, candidates, pid2gid_tmp, polyxx_tmp );
			eid++;
		}
	}

	//given that the labeling of the global ips is arbitrary and therefore also the order in which the identification process
	//is performed arbitrary, grains in contact with the domain boundary will be split potentially into multiple fragments
	//the identification and labeling of these fragments is arbitrary ending up with different absolute positions for
	//fragments as aufpunkte and thus potential multiple and different reference orientations
	//simplest 1-D explanation of that issue: consider a point chain indexed from 1,2,3,... to 10
	//starting at 1 one may find a grain add ip 2 and 3 to it, because the problem is periodic lets assume ip 10 belongs also to the grain
	//however, practically, by the point we arrive at 10 coming from ip 9, 10 is not yet assigned and the immediate periodic partner ip 10 is take
	//into account but still not assigned due to the strict sequentiality of the algorithm causing that the grain built of ip 1,2, 3
	//gets not found, never tested for orientation proximity, and therefore a new grain is generated instead
	//splitting up the grain in two fragments at ip 1,2, 3 and ip 10, respectively.
	//to remedy we can rerun the labeling algorithm

cout << "Initial labeling done" << endl;
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::AddReconstructGrainsRun1" + to_string(thisincrement)));

//	tic = MPI_Wtime();
//	//a workaround is to run the algorithm twice and relabel grains, then though we need a second buffer for the final id assignment
//	//on the fly updating of pid2gid may result in data races
//	vector<unsigned int> pid2gid_final;
//	for (size_t i = 0; i < n; ++i)
//		pid2gid_final.push_back( numeric_limits<unsigned int>::max() );

//	vector<grain> polyxx_final;
//	for (unsigned int mr = 0; mr < db.size(); mr++) { //analyze the grain assignment of all ips
//		memRegion* thisregion = db.at(mr);
//		size_t eip_s = thisregion->eipid_start;
//		size_t eip_n = thisregion->eipid_n;
//		size_t eid = eip_s;

//		for (size_t e = 0; e < eip_n; ++e) { //only the not periodic images of each element non-periodic images are material points to probe for grain identification
//			p3d matpoint = p3d(
//				thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx,
//				thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy,
//				thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz );

//			pp3rve1withguard.find_higherorder_neighbors( matpoint, candidates, Settings::GrainReconInspectionRadius);

//			//the candidates list is sorted by ascending distances to the location of matpoint

//			//check if any of the candidates has a grain already assigned and if so the orientation of such grain is
//			//similar to the one of the matpoint within SO3 according to LocalDisoriAngle
//			//if it is return gid of that most similar grain in terms of orientation proximity
//			//if not generate a new grain and return its gid
//			pid2gid_final[eid] = anyGrainAlreadyExistent2( eid, candidates, pid2gid_tmp, polyxx_tmp, polyxx_final );
//			eid++;
//		}
//	}

//cout << "Final labeling done" << endl;

//	toc = MPI_Wtime();
//	tictoc.push_back(plog(tic, toc, "Processing::AddReconstructGrainsRun2" + to_string(thisincrement)));

	//##MK::BEGIN OF DEBUG
	for (size_t i = 0; i < n; ++i) {
		if ( pid2gid_tmp[i] != numeric_limits<unsigned int>::max() )
			continue;
cout << "Integration point " << i << " has no grain assigned!" << endl;
	}
	//##MK::END OF DEBUG

	//get grain metadata
	tic = MPI_Wtime();
	for (size_t i = 0; i < n; ++i) {
		polyxx_tmp.at(pid2gid_tmp[i]).np++;
	}

	//##MK::single-grain resolved average orientation, grain reference orientation deviation, or orientation spread analysis

	//##MK::DEBUG render grain structure
	vector<p3d> ppp;
	vector<unsigned int> texid;
	for (size_t i = 0; i < n; ++i) {
		ppp.push_back( eid2p3d(i) );
		texid.push_back( eid2textureid(i) );
	}
	bool status = vtk_p3dm3( ppp, texid, thisincrement );
	status = vtk_gp3d( ppp, pid2gid_tmp, thisincrement );
	//##MK::END OF DEBUG

	//##MK::transfer grain data to other analysis routines

	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::ComputeGrainMetadata" + to_string(thisincrement)));

	write_reconstructed_grains1( polyxx_tmp );
}
*/



inline void specOutIncr::eid_write_gid( const size_t eid, const unsigned int gid )
{
	unsigned int nt = db.size();
	for (unsigned int mr = 0; mr < nt; ++mr) {
		if ( eid >= db.at(mr)->eipid_start && eid < db.at(mr)->eipid_end ) {
			size_t pos = eid - db.at(mr)->eipid_start;
			db.at(mr)->crystallite.GrainID[pos] = gid;
			return;
		}
	}
}


inline quat specOutIncr::eid2quaternion( const size_t eid )
{
	//##MK::very likely this logic is more performant
	unsigned int nt = db.size();
	for (unsigned int mr = 0; mr < nt; ++mr) {
		if ( eid < db.at(mr)->eipid_start || eid >= db.at(mr)->eipid_end )
			continue;
		//implicit else
		if ( eid >= db.at(mr)->eipid_start && eid < db.at(mr)->eipid_end ) {
//##MK::DEBUGsize_t localpos = eid - db.at(t)->eipid_start; quat res = db.at(t)->crystallite.q[eid - db.at(t)->eipid_start]; cout << localpos << "\t\t" << res << endl; return res;
			return db.at(mr)->crystallite.q[eid - db.at(mr)->eipid_start];
		}
	}
	return quat();//if nothing return identity quaternion as neutral element
}


inline p3d specOutIncr::eid2p3d( const size_t eid )
{
	unsigned int nt = db.size();
	for (unsigned int mr = 0; mr < nt; ++mr) {
		if ( eid >= db.at(mr)->eipid_start && eid < db.at(mr)->eipid_end ) {
			size_t pos = eid - db.at(mr)->eipid_start;
			return p3d( (db.at(mr)->grid.xyz0[pos].x + db.at(mr)->grid.dxyz_avg[pos].dx + db.at(mr)->grid.dxyz_flu[pos].dx),
						(db.at(mr)->grid.xyz0[pos].y + db.at(mr)->grid.dxyz_avg[pos].dy + db.at(mr)->grid.dxyz_flu[pos].dy),
						(db.at(mr)->grid.xyz0[pos].z + db.at(mr)->grid.dxyz_avg[pos].dz + db.at(mr)->grid.dxyz_flu[pos].dz)  );
		}
	}
	return p3d(); //error value
}


inline unsigned int specOutIncr::eid2textureid( const size_t eid )
{
	unsigned int nt = db.size();
	for (unsigned int mr = 0; mr < nt; ++mr) {
		if ( eid >= db.at(mr)->eipid_start && eid < db.at(mr)->eipid_end ) {
			return db.at(mr)->crystallite.TextureID[eid - db.at(mr)->eipid_start];
		}
	}
	return numeric_limits<unsigned int>::max(); //error value
}


inline unsigned int specOutIncr::eid2gid( const size_t eid )
{
	unsigned int nt = db.size();
	for (unsigned int mr = 0; mr < nt; ++mr) {
		if ( eid >= db.at(mr)->eipid_start && eid < db.at(mr)->eipid_end ) {
			return db.at(mr)->crystallite.GrainID[eid - db.at(mr)->eipid_start];
		}
	}
	return numeric_limits<unsigned int>::max(); //error value
}

dist specOutIncr::closestBoundaryInKernelIfAny( const size_t central_eid, vector<dist> const & candidates ) //, vector<size_t>& suspicious )
{
	//MK::assumes that candidates is sorted in ascending order!
	//for ( vector<dist>::iterator it = candidates.begin(); it != candidates.end(); ++it ) {*it->
	//for (unsigned int c = 0; c < candidates.size(); ++c) { cout << candidates[c].d << ","; } cout << endl; return dist();

	//get orientation of central point
	quat central_ori = eid2quaternion( central_eid );
	real_ori qcentral[4] = { central_ori.q0, central_ori.q1, central_ori.q2, central_ori.q3 };


	//##MK::DEBUG
	if ( candidates.size() > 0 ) {
		if ( central_eid != candidates.at(0).nbor_eipid) { //##MK::unlikely
cout << "ERROR::Inconsistent ID\t\t" << central_eid << "\t\t" << candidates.at(0).nbor_eipid << endl;
			return dist();
		}

		if ( Settings::VisHigherOrderNeighbor == (central_eid+1) ) { //##MK::only one of all points, MK::in *.xml this index is Fortran-based while in damaskpdt in C-based
			//##MK::collect metadata of neighborhood is expensive...
			vector<p3d> ppp;
			vector<unsigned int> pppid;
			vector<real_xyz> dxyz;
			vector<real_ori> diso;
			vector<unsigned int> texid;
			//first entry in candidates is the central integration point so disori between this and itself return nan
			for (size_t j = 0; j < candidates.size(); ++j) {
				ppp.push_back( eid2p3d(candidates[j].nbor_eipid) );
				pppid.push_back( candidates[j].nbor_eipid );
				dxyz.push_back( candidates[j].d );
				quat nbori = eid2quaternion( candidates[j].nbor_eipid );
				real_ori qnb[4] = { nbori.q0, nbori.q1, nbori.q2, nbori.q3 };
				diso.push_back( RADIANT2DEGREE(disorientation_angle_fcc( qcentral, qnb )) );
				texid.push_back( eid2textureid(candidates[j].nbor_eipid) );
			}
			bool status = vtk_nbhd3d( ppp, pppid, dxyz, diso, texid, thisincrement );
		}
	}
	//##MK::END OF DEBUG


	for ( size_t i = 1; i < candidates.size(); ++i ) { //we start at one because we exclude ourselves

		quat nbor_ori = eid2quaternion( candidates[i].nbor_eipid ); //##optimization via SIMD
		real_ori qnbor[4] = { nbor_ori.q0, nbor_ori.q1, nbor_ori.q2, nbor_ori.q3 };

		//##MK::implement different thresholding algorithms in here

		real_ori disori = disorientation_angle_fcc( qcentral, qnbor ); //MK::qcentral will not change!

		if ( disori >= Settings::CriticalDisoriAngle ) {
//cout << "centralid/kickingcand/totalcand/kickingid/len/disori\t\t" << central_eid << "\t\t" << i << "\t\t" << candidates.size() << ";" << candidates[i].nbor_eipid << ";" << candidates[i].d << ";" << RADIANT2DEGREE(disori) << endl;

			//##MK::DEBUG filter out ids of suspicious
			/*if ( i == 1 ) {
				suspicious.push_back( central_eid );
			}*/

			return candidates[i];
		}
	}

	//return flagged dummy value in case we were unsuccessful, either because no neighbors (almost impossible)
	//or (the usual case) point has no higher-order neighboring point
	//with disori > critical in the predefined inspection sphere

	return dist( (Settings::KernelRadius + EPSILON), numeric_limits<unsigned int>::max() ); //MK::indicating at least//dist();
}


/*void specOutIncr::write_histograms(vector<spatdist> const & buf)
{
	double tic = MPI_Wtime();

	//##MK::sort buffer content according to distance

	//##MK::suboptimal... one file per increment/rank
	string fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".StateVarSpatDistr.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "GeneralizedDistance;StateVariableValue\n";
		//##MK::DEBUG
		csvlog<< ";m^-2\n";
		csvlog << "GeneralizedDistance;StateVariableValue\n";

		for (size_t i = 0; i < buf.size(); ++i) {
			csvlog << buf.at(i).d << ";" << buf.at(i).sval << endl;
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local processing files", this->get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteStateVarSpatialDistr" + to_string(thisincrement)));
}*/


void specOutIncr::write_searchefficiency(vector<unsigned int> const & ntested )
{
	double tic = MPI_Wtime();

	string fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".HOSearchEfficiency.csv";

	ofstream srefflog;
	srefflog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (srefflog.is_open() == true) {
		//header
		srefflog << "CandidatesProbed\n";

		for (size_t i = 0; i < ntested.size(); ++i) {
			srefflog << ntested.at(i) << "\n";
		}
		//closure
		srefflog << endl;
		srefflog.flush();
		srefflog.close();
	}
	else {
		stopping("Unable to write local search efficiency to file", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::MetaWriteSearchEfficiency" + to_string(thisincrement)));
}


void specOutIncr::write_histograms2_ascii(vector<real_xyz> const & bufd, vector<slipsysdata_fcc> const & bufedge, vector<slipsysdata_fcc> const & bufdipo)
{
	//check input length consistence
	if ( bufd.size() != bufedge.size() || bufd.size() != bufdipo.size() ) {
		reporting("Writing spatial state variable distributions input buffer are inconsistent in their lengths!", myrank, 0, false);
		return;
	}

	double tic = MPI_Wtime();

	//##MK::sort buffer content according to distance

	//##MK::suboptimal... one file per increment/rank
	string fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".EdgeDensitySpatDistr.csv";

	ofstream csvlog_e;
	csvlog_e.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog_e.is_open() == true) {
		//header
		csvlog_e << "GeneralizedDistance;EdgeDensityB2;EdgeDensityB4;EdgeDensityB5;EdgeDensityC1;EdgeDensityC3;EdgeDensityC5;";
		csvlog_e << "EdgeDensityA2;EdgeDensityA3;EdgeDensityA6;EdgeDensityD1;EdgeDensityD4;EdgeDensityD6\n";
		csvlog_e << "1;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2\n";
		csvlog_e << "GeneralizedDistance;EdgeDensityB2;EdgeDensityB4;EdgeDensityB5;EdgeDensityC1;EdgeDensityC3;EdgeDensityC5;";
		csvlog_e << "EdgeDensityA2;EdgeDensityA3;EdgeDensityA6;EdgeDensityD1;EdgeDensityD4;EdgeDensityD6\n";
		//content
		for (size_t i = 0; i < bufd.size(); ++i) {
			csvlog_e << bufd.at(i) << ";";
			csvlog_e << bufedge.at(i).b2 << ";" << bufedge.at(i).b4 << ";" << bufedge.at(i).b5 << ";";
			csvlog_e << bufedge.at(i).c1 << ";" << bufedge.at(i).c3 << ";" << bufedge.at(i).c5 << ";";
			csvlog_e << bufedge.at(i).a2 << ";" << bufedge.at(i).a3 << ";" << bufedge.at(i).a6 << ";";
			csvlog_e << bufedge.at(i).d1 << ";" << bufedge.at(i).d4 << ";" << bufedge.at(i).d6 << endl;
		}
		//closure
		csvlog_e.flush();
		csvlog_e.close();
	}
	else {
		stopping("Unable to write local processing files for edge dislocation densities", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteEdgeDensitySpatialDistr" + to_string(thisincrement)));


	tic = MPI_Wtime();
	fn = get_prefix() + ".DipoleDensitySpatDistr.csv";
    ofstream csvlog_d;
	csvlog_d.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog_d.is_open() == true) {
		//header
		csvlog_d << "GeneralizedDistance;DipoleDensityB2;DipoleDensityB4;DipoleDensityB5;DipoleDensityC1;DipoleDensityC3;DipoleDensityC5;";
		csvlog_d << "DipoleDensityA2;DipoleDensityA3;DipoleDensityA6;DipoleDensityD1;DipoleDensityD4;DipoleDensityD6\n";
		csvlog_d << "1;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2;m^-2\n";
		csvlog_d << "GeneralizedDistance;DipoleDensityB2;DipoleDensityB4;DipoleDensityB5;DipoleDensityC1;DipoleDensityC3;DipoleDensityC5;";
		csvlog_d << "DipoleDensityA2;DipoleDensityA3;DipoleDensityA6;DipoleDensityD1;DipoleDensityD4;DipoleDensityD6\n";
		//content
		for (size_t i = 0; i < bufd.size(); ++i) {
			csvlog_d << bufd.at(i) << ";";
			csvlog_d << bufdipo.at(i).b2 << ";" << bufdipo.at(i).b4 << ";" << bufdipo.at(i).b5 << ";";
			csvlog_d << bufdipo.at(i).c1 << ";" << bufdipo.at(i).c3 << ";" << bufdipo.at(i).c5 << ";";
			csvlog_d << bufdipo.at(i).a2 << ";" << bufdipo.at(i).a3 << ";" << bufdipo.at(i).a6 << ";";
			csvlog_d << bufdipo.at(i).d1 << ";" << bufdipo.at(i).d4 << ";" << bufdipo.at(i).d6 << endl;
		}
		//closure
		csvlog_d.flush();
		csvlog_d.close();
	}
	else {
		stopping("Unable to write local processing files for dipole dislocation densities", get_myrank(), 0);
	}

	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteDipoleDensitySpatialDistr" + to_string(thisincrement)));
}


void specOutIncr::write_histograms2_binary(vector<MPI_DisloSpatDistr_Double> const & bufedge, vector<MPI_DisloSpatDistr_Double> const & bufdipo)
{
	double tic, toc;
	string prefix, fn;

	tic = MPI_Wtime();
	fn = get_prefix() + ".Incr." +  to_string(thisincrement) + ".EdgeDnsSpatDistr.NR." + to_string(bufedge.size()) + ".NC." + to_string(1 + db.back()->constitutive.rho_e_mult) + ".bin";

	//utilize MPI I/O library
	//open three files volume and nfaces and see
	//two files per increment, individual specOutIncrHdl write files individually and in parallel
	MPI_File ioHdl;
	MPI_Status ioSta;

	// open the file in create and write-only mode

	//##MK::write in blocks to avoid writing more than N*sizeof(MPI_DisloSpatDistr_Double) elements &bufedge[0], static_cast
	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdl);
	long long totalOffset = 0;
	size_t eTotal = bufedge.size();
	size_t ePerBlockTarget = Performance::MPIIOStripeSize / sizeof(MPI_DisloSpatDistr_Double);
	size_t ePerBlockCurr = ePerBlockTarget;
	for( size_t eWritten = 0; eWritten < eTotal;    ) {
		ePerBlockCurr = ((eWritten + ePerBlockTarget) < eTotal) ? ePerBlockTarget : eTotal - eWritten;
cout << "eTotal/Target/Curr/Written " << eTotal << ";" << ePerBlockTarget << ";" << ePerBlockCurr << ";" << eWritten << endl;
		//get raw pointer to memory location where data are stored
		const MPI_DisloSpatDistr_Double* thissection = bufedge.data();
cout << "thissection " << thissection << endl;
		//offset pointer to section where current block to be written begin
		thissection += eWritten;
cout << "thissection " << thissection << endl;

		MPI_File_write_at( ioHdl, totalOffset, thissection, ePerBlockCurr, MPI_DisloSpatDistr_Double_Type, &ioSta);
cout << "totalOffset " << totalOffset << endl;
		totalOffset = totalOffset + (sizeof(MPI_DisloSpatDistr_Double) * ePerBlockCurr);
		eWritten += ePerBlockCurr;
cout << "totalOffset/eWritten " << totalOffset << ";" << eWritten << endl;
	}

	MPI_File_close(&ioHdl);

	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteEdgeDnsSpatDistrMPIIO" + to_string(thisincrement)));


	tic = MPI_Wtime();
	fn = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".DipoleDnsSpatDistr.NR." + to_string(bufdipo.size()) + ".NC." + to_string(1 + db.back()->constitutive.rho_d_mult) + ".bin";

	//utilize MPI library write in blocks to avoid writing more than N*sizeof(MPI_DisloSpatDistr_Double) elements &bufedge[0], static_cast

	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &ioHdl);
	totalOffset = 0;
	eTotal = bufdipo.size();
	ePerBlockTarget = Performance::MPIIOStripeSize / sizeof(MPI_DisloSpatDistr_Double);
	ePerBlockCurr = ePerBlockTarget;
	for( size_t eWritten = 0; eWritten < eTotal;    ) {
		ePerBlockCurr = ((eWritten + ePerBlockTarget) < eTotal) ? ePerBlockTarget : eTotal - eWritten;

		const MPI_DisloSpatDistr_Double* thissection = bufdipo.data();
		thissection += eWritten;

		MPI_File_write_at( ioHdl, totalOffset, thissection, ePerBlockCurr, MPI_DisloSpatDistr_Double_Type, &ioSta);

		totalOffset = totalOffset + (sizeof(MPI_DisloSpatDistr_Double) * ePerBlockCurr);

		eWritten += ePerBlockCurr;
	}

	MPI_File_close(&ioHdl);
	toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteDipoleDnsSpatDistrMPIIO" + to_string(thisincrement)));
}


/*
unsigned int specOutIncr::anyGrainAlreadyAssigned( const size_t central_eid, vector<dist> const & candidates,
	vector<unsigned int> const & p2g, vector<grain> & grpool )
{
	//MK::will and is allowed to add elements to grpool
	//MK::assumes that candidates is sorted in ascending order based on inter-point distances!

	//get orientation of central point
	quat central_ori = eid2quaternion( central_eid );
	real_ori qcentral[4] = { central_ori.q0, central_ori.q1, central_ori.q2, central_ori.q3 };

	real_xyz closestDistance = numeric_limits<real_xyz>::max();
	real_ori closestDisori = DEGREE2RADIANT(180.0);
	unsigned int closestIndex = numeric_limits<unsigned int>::max();
	unsigned int closestGID = numeric_limits<unsigned int>::max();
	bool found = false;

	//##MK::BEGIN OF DEBUG
	//unsigned int debuginspect = 4;
	//if ( central_eid == debuginspect ) {
	//	cout << "Candidate " << debuginspect << endl;
	//	cout << "Candidate " << central_eid << endl;
	//	for (size_t k = 0; k < candidates.size(); ++k)
	//		cout << candidates.at(k).d << "\t\t" << candidates.at(k).nbor_eipid << "\t\t" << p2g[candidates[k].nbor_eipid] << endl;
	//}
	//##MK::END OF DEBUG

	vector<unsigned int> matches;	//already assigned grain IDs to which central_ori has disori within acceptance range

	unsigned int ni = static_cast<size_t>(candidates.size());
	for ( size_t i = 1; i < ni; ++i ) { //we start at 1 because we exclude ourselves
		unsigned int nbor_gid = p2g[candidates[i].nbor_eipid]; //translate global ip ID to grain ID

		//##MK::BEGIN OF DEBUG
		//if ( central_eid == debuginspect )
		//	cout << i << "\t\t" << nbor_gid << endl "-->";
		//##MK::END OF DEBUG

		if ( nbor_gid != numeric_limits<unsigned int>::max() ) {
			real_ori qgr[4] = { grpool.at(nbor_gid).ori.q0, grpool.at(nbor_gid).ori.q1,
					grpool.at(nbor_gid).ori.q2, grpool.at(nbor_gid).ori.q3 };

			real_ori disori = disorientation_angle_fcc( qcentral, qgr );

			//##MK::BEGIN OF DEBUG
			//if ( central_eid == debuginspect )
			//	cout << disori << "\t\t" << Settings::GrainReconLocalDisoriAngle << endl;
			//##MK::END OF DEBUG

			if ( disori <= Settings::GrainReconLocalDisoriAngle ) {
				matches.push_back( i ); //remember matching indices of the interesting candidates
				if ( disori <= closestDisori ) { //MK::given the default this will at least be executed once if below LocalDisoriThreshold
					found = true;
					//best matching in terms of orientation in search radius results not necessarily in fusing spatial contiguous
					//grains when search radius is too larger than a few integration points
					closestDisori = disori;
					closestIndex = i;
					closestGID = grpool.at(nbor_gid).gid;
				}
			} //check for local disori angle
		} //checked whether candidate point had already grain assigned
	} //MK::we have to test all candidates for spatial contiguity

	if ( found == true ) { //most likely
cout << "centraleid/closestGID/closestDISORI\t\t" << central_eid << "\t\t" << closestGID << "\t\t" << RADIANT2DEGREE(closestDisori) << endl;
		return closestGID;
	}
	//implicit else --- either no grain yet assigned to any neighboring ip or all their assignments of oris too different in SO3
	grain agr = grain( central_ori, static_cast<unsigned int>(grpool.size()) );
	grpool.push_back( agr );

cout << "centraleid/newGrainGID/q0123\t\t" << central_eid << "\t\t" << grpool.back().gid << "\t\t";
cout << grpool.back().ori.q0 << ";" << grpool.back().ori.q1 << ";" << grpool.back().ori.q2 << ";" << grpool.back().ori.q3 << endl;

	return grpool.back().gid;
}


unsigned int specOutIncr::anyGrainAlreadyExistent1( const size_t central_eid,
		vector<dist> const & candidates, vector<unsigned int> const & p2g, vector<grain> & grpool )
{
	//MK::will and is allowed to add elements to grpool
	//MK::assumes that candidates is sorted in ascending order based on interpoint distances!
	//MK::THE SORT NEEDS TO BE STABLE !

	//get orientation of central point
	quat central_ori = eid2quaternion( central_eid );
	real_ori qcentral[4] = { central_ori.q0, central_ori.q1, central_ori.q2, central_ori.q3 };
	cout << "Central ori\t\t" << qcentral[0] << ";" << qcentral[1] << ";" << qcentral[2] << ";" << qcentral[3] << endl;

	//real_xyz closestDistance = numeric_limits<real_xyz>::max();
	real_ori closestDisori = DEGREE2RADIANT(180.0);
	unsigned int closestGID = numeric_limits<unsigned int>::max();
	bool found = false;

	//##MK::BEGIN OF DEBUG
	//unsigned int debuginspect = 4;
	//if ( central_eid == debuginspect ) {
	//	cout << "Candidate " << debuginspect << endl;
		cout << "Candidate " << central_eid << endl;
		for (size_t k = 0; k < candidates.size(); ++k)
			cout << candidates.at(k).d << "\t\t" << candidates.at(k).nbor_eipid << "\t\t" << p2g[candidates[k].nbor_eipid] << endl;
	//}
	//##MK::END OF DEBUG

	for ( size_t i = 1; i < candidates.size(); ++i ) { //we start at 1 because we exclude ourselves
		unsigned int nbor_gid = p2g[candidates[i].nbor_eipid]; //translate global ip ID to grain ID

		//##MK::BEGIN OF DEBUG
		//if ( central_eid == debuginspect )
		//	cout << i << "\t\t" << nbor_gid << endl "-->";
		//##MK::END OF DEBUG

		if ( nbor_gid != numeric_limits<unsigned int>::max() ) { //if this was possible
			real_ori qgr[4] = { grpool.at(nbor_gid).ori.q0, grpool.at(nbor_gid).ori.q1,
					grpool.at(nbor_gid).ori.q2, grpool.at(nbor_gid).ori.q3 };

			real_ori disori = disorientation_angle_fcc( qcentral, qgr );

			//##MK::BEGIN OF DEBUG
			//if ( central_eid == debuginspect )
				cout << disori << "\t\t" << Settings::GrainReconLocalDisoriAngle << "\t\t" << qgr[0] << ";" << qgr[1] << ";" << qgr[2] << ";" << qgr[3] << endl;
			//##MK::END OF DEBUG

			if ( disori <= Settings::GrainReconLocalDisoriAngle ) {
				if ( disori < closestDisori ) { //candidates[i].d <= (closestDistance+EPSILON) ) {
					//MK::given the defaults of disori, distance and the fact that the list is sorted ascendingly
					//the condition will at least once be executed ...

					//##MK:: what about multiple numerically equally distant grain options to choose from, then choose the closests in disori of the equal ones

					found = true;
					closestDisori = disori;
					closestGID = grpool.at(nbor_gid).gid;
					//##MK::for now DEBUG breaking out having found one
					break;
				}
			} //check for local disori angle
		} //checked whether candidate point had already grain assigned

	} //MK::we have to test all candidates for spatial contiguity

	if ( found == true ) { //most likely
cout << "centraleid/closestGID/closestDISORI\t\t" << central_eid << "\t\t" << closestGID << "\t\t" << RADIANT2DEGREE(closestDisori) << endl;
		return closestGID;
	}
	//implicit else --- either no grain yet assigned to any neighboring ip or all their assignments of oris too different in SO3
	grain agr = grain( central_ori, static_cast<unsigned int>(grpool.size()) );
	grpool.push_back( agr );

cout << "centraleid/newGrainGID/q0123\t\t" << central_eid << "\t\t" << grpool.back().gid << "\t\t";
cout << grpool.back().ori.q0 << ";" << grpool.back().ori.q1 << ";" << grpool.back().ori.q2 << ";" << grpool.back().ori.q3 << endl;

	return grpool.back().gid;
}
*/

/*unsigned int specOutIncr::anyGrainAlreadyExistent2( const size_t central_eid, vector<dist> const & candidates,
				vector<unsigned int> const & p2gtmp, vector<grain> const & oldgrpool, vector<grain> & newgrpool )
{
	//get grain assigned to central point
	unsigned central_gid = p2gtmp[central_eid];
	quat central_ori = oldgrpool.at(central_gid).ori;
	real_ori qcentral[4] = { central_ori.q0, central_ori.q1, central_ori.q2, central_ori.q3 };

	real_xyz closestDistance = numeric_limits<real_xyz>::max();
	real_ori closestDisori = DEGREE2RADIANT(180.0);
	unsigned int closestGID = numeric_limits<unsigned int>::max();
	bool found = false;

	for ( size_t i = 1; i < candidates.size(); ++i ) { //we start at 1 because we exclude ourselves
		unsigned int nbor_gid = p2gtmp[candidates[i].nbor_eipid]; //translate global ip ID to grain ID
		quat nbor_ori = oldgrpool.at(nbor_gid).ori;
		real_ori qgr[4] = { nbor_ori.q0, nbor_ori.q1, nbor_ori.q2, nbor_ori.q3 };

		real_ori disori = disorientation_angle_fcc( qcentral, qgr );

		if ( disori <= Settings::GrainReconLocalDisoriAngle ) {
			if ( disori < closestDisori ) { //candidates[i].d <= (closestDistance+EPSILON) ) {
				//MK::given the defaults of disori, distance and the fact that the list is sorted ascendingly
				//the condition will at least once be executed ...

				//##MK:: what about multiple numerically equally distant grain options to choose from, then choose the closests in disori of the equal ones

				found = true;
				closestDisori = disori;
				closestGID = oldgrpool.at(nbor_gid).gid;
				//##MK::for now DEBUG breaking out having found one
				break;
			}
		} //check for local disori angle
	} //MK::we have to test all candidates for spatial contiguity

	if ( found == true ) { //most likely
cout << "centraleid/closestGID/closestDISORI\t\t" << central_eid << "\t\t" << closestGID << "\t\t" << RADIANT2DEGREE(closestDisori) << endl;
		return closestGID;
	}
	//implicit else --- either no grain yet assigned to any neighboring ip or all their assignments of oris too different in SO3
	grain agr = grain( central_ori, static_cast<unsigned int>(newgrpool.size()) );
	newgrpool.push_back( agr );

cout << "centraleid/newGrainGID/q0123\t\t" << central_eid << "\t\t" << newgrpool.back().gid << "\t\t";
cout << newgrpool.back().ori.q0 << ";" << newgrpool.back().ori.q1 << ";" << newgrpool.back().ori.q2 << ";" << newgrpool.back().ori.q3 << endl;

	return newgrpool.back().gid;

	return numeric_limits<unsigned int>::max();
}*/




void specOutIncr::write_reconstructed_grains1( vector<grain> const & grpool )
{
	double tic = MPI_Wtime();

	string fn = get_prefix() + ".Incr." +  to_string(thisincrement) + ".GrainData.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "GrainID;NumberOfIPsAssigned;RefQuatQ0;RefQuat1;RefQuatQ2;RefQuatQ3\n";
		csvlog<< ";1;;;;\n";
		csvlog << "GrainID;NumberOfIPsAssigned;RefQuatQ0;RefQuat1;RefQuatQ2;RefQuatQ3\n";

		for (size_t i = 0; i < grpool.size(); ++i) {

			csvlog << grpool[i].gid << ";" << grpool[i].np << ";"
					<< grpool[i].ori.q0 << ";" << grpool[i].ori.q1 << ";" << grpool[i].ori.q2 << ";" << grpool[i].ori.q3 << "\n";
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local grain reconstruction meta data to file", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteReconstructedGrainMetaData" + to_string(thisincrement)));
}

void specOutIncr::write_identified_grains()
{
	double tic = MPI_Wtime();

	string fn = get_prefix() + ".Incr." +  to_string(thisincrement) + ".GrainData.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "GrainID;NumberOfIPsAssigned;RefQuatQ0;RefQuat1;RefQuatQ2;RefQuatQ3;GROD\n";
		csvlog<< ";1;;;;;degree\n";
		csvlog << "GrainID;NumberOfIPsAssigned;RefQuatQ0;RefQuat1;RefQuatQ2;RefQuatQ3;GROD\n";

		for (size_t i = 0; i < grains.polyxx.size(); ++i) {
			csvlog << grains.polyxx[i].gid << ";" << grains.polyxx[i].np << ";"
					<< grains.polyxx[i].ori.q0 << ";" << grains.polyxx[i].ori.q1 << ";"
						<< grains.polyxx[i].ori.q2 << ";" << grains.polyxx[i].ori.q3 << ";"
							<< "n/a" << "\n";
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local grain reconstruction meta data to file", get_myrank(), 0);
	}

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::WriteIdentifiedGrainMetaData" + to_string(thisincrement)));
}


void specOutIncr::write_grainid_and_quaternions()
{
	//##MK::DEBUG highly developmental
	unsigned int* u32_buffer = NULL;
	float* f32_buffer = NULL;
	size_t n = 0;
	for(unsigned int gr = 0; gr < grains.polyxx.size(); gr++) {
		if ( grains.ipsupport.at(gr) != NULL ) {
			n = n + grains.ipsupport.at(gr)->size();
		}
		else {
cout << "Writing of grainID and quaternion clouds was unsuccessful as uip data for at least one grain missing" << endl;
			return;
		}
	}

	if ( n != head.NXYZ ) { //#MK::was this->
cout << "Writing of grainID and quaternion clouds was unsuccessful as total number of uip data inconsistent to head.NXYZ expectation" << endl;
		return;
	}

	//##MK::check 4*n against int max
	if ( 4*n > (numeric_limits<int>::max() - 1) ) {
cout << "Writing of grainID and quaternion clouds was unsuccessful because current file buffering strategy cannot deal with so large input dataset" << endl;
		return;
	}


	try { u32_buffer = new unsigned int[n]; }
	catch (bad_alloc &croak) {
cout << "Writing of grainID and quaternion clouds was unsuccessful during u32_buffer allocation" << endl;
		return;
	}

	try { f32_buffer = new float[4*n]; }
	catch (bad_alloc &croak) {
cout << "Writing of grainID and quaternion clouds was unsuccessful during f32_buffer allocation" << endl;
		if (u32_buffer != NULL) { delete [] u32_buffer; u32_buffer = NULL; }
		return;
	}

	//set up inconsistency indicating initial values
	size_t u32_offset = 0;
	size_t f32_offset = 0;
	for(size_t i = 0; i < n;   ) { u32_buffer[i] = numeric_limits<unsigned int>::max(); ++i; }
	for(size_t i = 0; i < 4*n; ) {
		f32_buffer[i+0] = numeric_limits<float>::max();
		f32_buffer[i+1] = numeric_limits<float>::max();
		f32_buffer[i+2] = numeric_limits<float>::max();
		f32_buffer[i+3] = numeric_limits<float>::max();
		i += 4;
	}


	for(unsigned int gr = 0; gr < grains.polyxx.size(); gr++) {
		unsigned int thisgid = grains.polyxx.at(gr).gid;
		vector<unsigned int>* theseips = grains.ipsupport.at(gr);
		size_t ni = theseips->size();

		for(size_t i = 0; i < ni; ++i) {
			u32_buffer[u32_offset] = thisgid;
			u32_offset += 1;

			quat thisq = eid2quaternion( static_cast<unsigned int>(theseips->at(i)) );

			real_ori q[4] = {thisq.q0, thisq.q1, thisq.q2, thisq.q3 };

			passive2active( q );

			//because have been converted during reading in and MTex utilizes active convention as it does DAMASK
			//MK::THIS WILL CHANGE IN FUTURE VERSIONS OF DAMASK IN FAVOR FOR THE PASSIVE CONVENTION...

			f32_buffer[f32_offset+0] = static_cast<float>(q[0]);
			f32_buffer[f32_offset+1] = static_cast<float>(q[1]);
			f32_buffer[f32_offset+2] = static_cast<float>(q[2]);
			f32_buffer[f32_offset+3] = static_cast<float>(q[3]);
			f32_offset += 4;
		} //transfer all quaternions for uips from grain gr
	}

	//write to binary file in one call

	string fn1 = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".QuatCloudGID.bin";
	string fn2 = get_prefix() +  ".Incr." +  to_string(thisincrement) + ".QuatCloudORI.bin";

	//##MK::add checks
	ofstream bin_debug1(fn1.c_str(), ios::out | ios::binary);
	bin_debug1.write( (char*) u32_buffer, n * sizeof(unsigned int));
	bin_debug1.close();

	ofstream bin_debug2(fn2.c_str(), ios::out | ios::binary);
	bin_debug2.write( (char*) f32_buffer, 4 * n * sizeof(float));
	bin_debug2.close();

	//deallocate naive buffer
	if (u32_buffer != NULL) { delete [] u32_buffer; u32_buffer = NULL; }
	if (f32_buffer != NULL) { delete [] f32_buffer; f32_buffer = NULL; }

cout << "Writing of grainID and quaternion clouds successful" << endl;
}

void specOutIncr::report_rveshapes()
{
	cout << "Reporting RVE Shape corner points" << endl;
	if ( rveShape.sp.size() == rveShape.nm.size() ) {
		for(size_t i = 0; i < rveShape.sp.size(); ++i ) {
			cout << rveShape.sp.at(i) << rveShape.nm.at(i) << endl;
		}
	}
	rveShape.sp.clear();
	rveShape.nm.clear();
}

void specOutIncr::free_increment_heavydata()
{
	double tic = MPI_Wtime();

	for ( size_t t = 0; t < db.size(); ++t ) {
		if ( db.at(t) != NULL ) {
			//##MK::deallocating subunits of memRegion may not be invoked?
			delete db.at(t);
			db.at(t) = NULL;
		}
	}
	db.clear();

	rveBaseIni = bv3x3();
	rveBaseDef = bv3x3();

	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "Processing::FreeMemoryHeavyData" + to_string(thisincrement)));

	thisincrement = numeric_limits<unsigned int>::max();
	thiswrittenincrement = numeric_limits<unsigned int>::max();
}


void specOutIncr::free_increment_bvh_xyzm2()
{
	pp3rve1withguard.destroy_bvh_xyzm2();
}


void specOutIncr::free_increment_bvh_p3dm1()
{
	pp3rve27.destroy_bvh_p3dm1();
}


void specOutIncr::free_increment_grains()
{
	grains.destroy();
}


void specOutIncr::free_increment_sdf()
{
	for(size_t s = 0; s < sdf.size(); ++s) {
		if (sdf.at(s) != NULL) {
			delete sdf.at(s);
			sdf.at(s) = NULL;
		}
	}
	sdf.clear();
}




//MK::tags should be disjoint!
#define TAG_VAVG		11
#define TAG_FPAVG		12
#define TAG_FAVG		13
#define TAG_PAVG		14
#define TAG_EPSAVG		15
#define TAG_CAUAVG		16

#define TAG_VMISES		17

/*
void specOutIncr::report_flowcurve()
{
	//collect RVE-volume averaged increment data from the MPI processes at MASTER into one flow curve
	//write VTK file showing the positions of all ions in reconstructed space
	reporting( "Writing flowcurve to file", myrank, 0, true ); //##MK::true or false

	//collect from MPI_COMM_WORLD
	//##MK::allocation error handling and MPI error handling
	double tic = MPI_Wtime();
	vector<rveAverageResults> iobuffer;
	unsigned int incr = head.sincr;
	unsigned int targetincr = Settings::IncrementFirst;
	for ( unsigned int loadcase = 0; loadcase < head.loadcases; ++loadcase) {
		for ( unsigned int li = 0; li < head.loadcasesmeta.at(loadcase).nincr; ++li ) {
			if ( incr == targetincr ) {
				//not necessarily for each converged increment dump data were written
				if ( incr2rank.find(incr) != incr2rank.end() ) {
					//MK::all processes know who computes which increment
					int whom = incr2rank.find(incr)->second;
					int me = get_myrank(); //#MK::was this->

					if ( me == MASTER ) {
						if ( whom == MASTER ) { //no comm required as master has local data already, buffer directly ##MK::could also write immediately...
							//##MK::find correct data, could be optimized
							unsigned int j = 0;
							unsigned int nj = avg.size();
							for (		; j < nj; ++j) {
								if ( avg.at(j).loadcaseID != loadcase )	continue;
								if ( avg.at(j).localincrID != li )		continue;
								if ( avg.at(j).globalincrID != incr)	continue;
								//not continued, so entry j is the correct one
								break;
							}
							if ( j >= avg.size() ) { complaining("Unable to find RVE average results", myrank, 0); }
							else {
								iobuffer.push_back( rveAverageResults( avg.at(j).Vtotal, avg.at(j).Favgrve, avg.at(j).Pavgrve,
										avg.at(j).Strainavgrve, avg.at(j).Cauchyavgrve, avg.at(j).Equivavgrve, loadcase, li, incr ) );
							}
						}
						else { //comm required because different rank than MASTER processed incr
							real_xyz recvV = static_cast<real_xyz>(0.0);
							MPI_Tensor3x3_Double recvF = MPI_Tensor3x3_Double();
							t3x3 recvP = t3x3();
							t3x3 recvEpsilon = t3x3();
							t3x3 recvCauchy = t3x3();
							MPI_Recv( &recvF, 1, MPI_Tensor3x3_Double_Type, whom, TAG_FAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_VonMises_Double recvvm = MPI_VonMises_Double();
							MPI_Recv( &recvvm, 1, MPI_VonMises_Double_Type, whom, TAG_VMISES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							iobuffer.push_back( rveAverageResults(
								recvV, t3x3(static_cast<real_m33>(recvF.a11), static_cast<real_m33>(recvF.a12), static_cast<real_m33>(recvF.a13),
											static_cast<real_m33>(recvF.a21), static_cast<real_m33>(recvF.a22), static_cast<real_m33>(recvF.a23),
											static_cast<real_m33>(recvF.a31), static_cast<real_m33>(recvF.a32), static_cast<real_m33>(recvF.a33)),
												recvP, recvEpsilon, recvCauchy, vMises(static_cast<real_m33>(recvvm.eps), static_cast<real_m33>(recvvm.sigma)),
												loadcase, li, incr)  );
						}
					}  //tasks for master defined
					else {  //all other processes
						if ( whom == me ) {
							unsigned int j = 0;
							for (		; j < avg.size(); ++j) {
								if ( avg.at(j).loadcaseID != loadcase )		continue;
								if ( avg.at(j).localincrID != li )			continue;
								if ( avg.at(j).globalincrID != incr)		continue;
								//not continued, so entry j is the correct one
								break;
							}
							if ( j >= avg.size() ) {
								complaining("Unable to find RVE average results", myrank, 0);
								MPI_Tensor3x3_Double snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FAVG, MPI_COMM_WORLD);
								MPI_VonMises_Double sndvm = MPI_VonMises_Double(vMises());
								MPI_Send( &sndvm, 1, MPI_VonMises_Double_Type, MASTER, TAG_VMISES, MPI_COMM_WORLD);
							}
							else {
								MPI_Tensor3x3_Double snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Favgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FAVG, MPI_COMM_WORLD);

								MPI_VonMises_Double sndvm = MPI_VonMises_Double(avg.at(j).Equivavgrve);
								MPI_Send( &sndvm, 1, MPI_VonMises_Double_Type, MASTER, TAG_VMISES, MPI_COMM_WORLD);
							}
						}
						//me is neither MASTER nor whom, i.e. not responsible for I/O for this incr
					} //tasks for non-master processes defined

					//##MK::this barrier not necessary only for DEBUGGING !
					//MPI_Barrier(MPI_COMM_WORLD);
				} //processing of desired targetincr performed
			}

			incr++;
			targetincr += Settings::IncrementOffset;
			if ( targetincr > Settings::IncrementLast )
				targetincr = Settings::IncrementLast; //thus enforcing a stop criterion, as incr keeps increasing...
		} //find targetincrements in steps of current loadcase
	} //find targetincrements in next loadcase

	if ( get_myrank() == MASTER ) {
		//all results now on master, so pipe to file
		ofstream flowcurve_sampled;

		//each converged DAMASK spectral solver increment represent one point on the simulated flow curve
		string fn = get_prefix() + ".Incr.All.Flowcurve.csv";
		flowcurve_sampled.open( fn.c_str(), ofstream::out | std::ofstream::trunc );  //MK::thereby discarding data in already existent files of the same name!

		if ( flowcurve_sampled.is_open() == true ) {
			diary whatwasdone = diary(Settings::SimID);
			flowcurve_sampled << setprecision(32);
			flowcurve_sampled << whatwasdone << endl;

			//Origin inspired three-column header with column name tag, units, and proposal explanation
			flowcurve_sampled << "LoadCaseID;LoadCaseIncr;IncrementID;vMisesEqvTrueStrain;vMisesEqCauchyStress;Volume;F11;F12;F13;F21;F22;F23;F31;F32;F33\n";
			flowcurve_sampled << ";;;;Pa;;;;;;;;;;\n";
			flowcurve_sampled << "LoadCaseID;LoadCaseIncr;IncrementID;vMisesEqvTrueStrain;vMisesEqCauchyStress;Volume;F11;F12;F13;F21;F22;F23;F31;F32;F33\n";

			for ( unsigned int tincr = 0; tincr < iobuffer.size(); ++tincr) {
				flowcurve_sampled << iobuffer.at(tincr).loadcaseID << ";";
				flowcurve_sampled << iobuffer.at(tincr).localincrID << ";";
				flowcurve_sampled << iobuffer.at(tincr).globalincrID << ";";
				flowcurve_sampled << iobuffer.at(tincr).Equivavgrve.vMisesEquivStrain << ";";
				flowcurve_sampled << iobuffer.at(tincr).Equivavgrve.vMisesEquivStress << ";";
				flowcurve_sampled << iobuffer.at(tincr).Vtotal << ";";
				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a11 << ";" << iobuffer.at(tincr).Favgrve.a12 << ";" << iobuffer.at(tincr).Favgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a21 << ";" << iobuffer.at(tincr).Favgrve.a22 << ";" << iobuffer.at(tincr).Favgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a31 << ";" << iobuffer.at(tincr).Favgrve.a32 << ";" << iobuffer.at(tincr).Favgrve.a33 << "\n";
			}

			flowcurve_sampled.flush();
			flowcurve_sampled.close();
		}
		else {
			complaining("Flowcurve could not be written to file", myrank, 0);
		}
	} //master finished processing file

	//MK::necessary, whenever attempting MPI parallel computations requiring syncing thereafter!
	MPI_Barrier(MPI_COMM_WORLD);
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::ReportingFlowcurve"));
}
*/

void specOutIncr::report_flowcurve2()
{
	//collect RVE-volume averaged increment data from the MPI processes at MASTER into one flow curve
	//write VTK file showing the positions of all ions in reconstructed space
	reporting( "Writing flowcurve to file", myrank, 0, true ); //##MK::true or false

	//collect from MPI_COMM_WORLD
	//##MK::allocation error handling and MPI error handling
	double tic = MPI_Wtime();
	vector<rveAverageResults> iobuffer;
	unsigned int incr = head.sincr;
	unsigned int targetincr = Settings::IncrementFirst;
	for ( unsigned int loadcase = 0; loadcase < head.loadcases; ++loadcase) {
		for ( unsigned int li = 0; li < head.loadcasesmeta.at(loadcase).nincr; ++li ) {
			if ( incr == targetincr ) {
				//not necessarily for each converged increment dump data were written
				if ( incr2rank.find(incr) != incr2rank.end() ) {
					//MK::all processes know who computes which increment
					int whom = incr2rank.find(incr)->second;
					int me = get_myrank(); //#MK::was this->

					if ( me == MASTER ) {
						if ( whom == MASTER ) { //no comm required as master has local data already, buffer directly ##MK::could also write immediately...
							//##MK::find correct data, could be optimized
							unsigned int j = 0;
							unsigned int nj = avg.size();
							for (		; j < nj; ++j) {
								if ( avg.at(j).loadcaseID != loadcase )	continue;
								if ( avg.at(j).localincrID != li )		continue;
								if ( avg.at(j).globalincrID != incr)	continue;
								//not continued, so entry j is the correct one
								break;
							}
							if ( j >= avg.size() ) { complaining("Unable to find RVE average results", myrank, 0); }
							else {
								iobuffer.push_back( rveAverageResults( avg.at(j).Vtotal, avg.at(j).Fpavgrve, avg.at(j).Favgrve, avg.at(j).Pavgrve,
										avg.at(j).Strainavgrve, avg.at(j).Cauchyavgrve, avg.at(j).Equivavgrve, loadcase, li, incr ) );
							}
						}
						else { //comm required because different rank than MASTER processed incr
							real_xyz recvV = static_cast<real_xyz>(0.0);
							MPI_Recv( &recvV, 1, MPI_DOUBLE, whom, TAG_VAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							//##MK::can be optimized by packing all pieces of information into one data container
							//##MK::for now better use the safe and easy to read solution of individual messages per data element
							MPI_Tensor3x3_Double recvFp = MPI_Tensor3x3_Double();
							MPI_Recv( &recvFp, 1, MPI_Tensor3x3_Double_Type, whom, TAG_FPAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_Tensor3x3_Double recvF = MPI_Tensor3x3_Double();
							MPI_Recv( &recvF, 1, MPI_Tensor3x3_Double_Type, whom, TAG_FAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_Tensor3x3_Double recvP = MPI_Tensor3x3_Double();
							MPI_Recv( &recvP, 1, MPI_Tensor3x3_Double_Type, whom, TAG_PAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_Tensor3x3_Double recvEps = MPI_Tensor3x3_Double();
							MPI_Recv( &recvEps, 1, MPI_Tensor3x3_Double_Type, whom, TAG_EPSAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_Tensor3x3_Double recvCau = MPI_Tensor3x3_Double();
							MPI_Recv( &recvCau, 1, MPI_Tensor3x3_Double_Type, whom, TAG_CAUAVG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							MPI_VonMises_Double recvvm = MPI_VonMises_Double();
							MPI_Recv( &recvvm, 1, MPI_VonMises_Double_Type, whom, TAG_VMISES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							iobuffer.push_back( rveAverageResults(
									recvV,
									t3x3( 	static_cast<real_xyz>(recvFp.a11), static_cast<real_xyz>(recvFp.a12), static_cast<real_xyz>(recvFp.a13),
											static_cast<real_xyz>(recvFp.a21), static_cast<real_xyz>(recvFp.a22), static_cast<real_xyz>(recvFp.a23),
											static_cast<real_xyz>(recvFp.a31), static_cast<real_xyz>(recvFp.a32), static_cast<real_xyz>(recvFp.a33) ),

									t3x3( 	static_cast<real_xyz>(recvF.a11), static_cast<real_xyz>(recvF.a12), static_cast<real_xyz>(recvF.a13),
											static_cast<real_xyz>(recvF.a21), static_cast<real_xyz>(recvF.a22), static_cast<real_xyz>(recvF.a23),
											static_cast<real_xyz>(recvF.a31), static_cast<real_xyz>(recvF.a32), static_cast<real_xyz>(recvF.a33) ),

									t3x3( 	static_cast<real_xyz>(recvP.a11), static_cast<real_xyz>(recvP.a12), static_cast<real_xyz>(recvP.a13),
											static_cast<real_xyz>(recvP.a21), static_cast<real_xyz>(recvP.a22), static_cast<real_xyz>(recvP.a23),
											static_cast<real_xyz>(recvP.a31), static_cast<real_xyz>(recvP.a32), static_cast<real_xyz>(recvP.a33) ),

									t3x3( 	static_cast<real_xyz>(recvEps.a11), static_cast<real_xyz>(recvEps.a12), static_cast<real_xyz>(recvEps.a13),
											static_cast<real_xyz>(recvEps.a21), static_cast<real_xyz>(recvEps.a22), static_cast<real_xyz>(recvEps.a23),
											static_cast<real_xyz>(recvEps.a31), static_cast<real_xyz>(recvEps.a32), static_cast<real_xyz>(recvEps.a33) ),

									t3x3( 	static_cast<real_xyz>(recvCau.a11), static_cast<real_xyz>(recvCau.a12), static_cast<real_xyz>(recvCau.a13),
											static_cast<real_xyz>(recvCau.a21), static_cast<real_xyz>(recvCau.a22), static_cast<real_xyz>(recvCau.a23),
											static_cast<real_xyz>(recvCau.a31), static_cast<real_xyz>(recvCau.a32), static_cast<real_xyz>(recvCau.a33) ),
									vMises(static_cast<real_m33>(recvvm.eps), static_cast<real_m33>(recvvm.sigma) ),
											loadcase, li, incr)
									);
						}
					}  //tasks for master defined
					else {  //all other processes
						if ( whom == me ) {
							unsigned int j = 0;
							for (		; j < avg.size(); ++j) {
								if ( avg.at(j).loadcaseID != loadcase )		continue;
								if ( avg.at(j).localincrID != li )			continue;
								if ( avg.at(j).globalincrID != incr)		continue;
								//not continued, so entry j is the correct one
								break;
							}
							if ( j >= avg.size() ) {
								complaining("Unable to find RVE average results", myrank, 0);
								double snddbl = 0.0;
								MPI_Send( &snddbl, 1, MPI_DOUBLE, MASTER, TAG_VAVG, MPI_COMM_WORLD);

								MPI_Tensor3x3_Double snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FPAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_PAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_EPSAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(failt3x3());
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_CAUAVG, MPI_COMM_WORLD);

								MPI_VonMises_Double sndvm = MPI_VonMises_Double(vMises());
								MPI_Send( &sndvm, 1, MPI_VonMises_Double_Type, MASTER, TAG_VMISES, MPI_COMM_WORLD);
							}
							else {
								double snddbl = static_cast<double>(avg.at(j).Vtotal);
								MPI_Send( &snddbl, 1, MPI_DOUBLE, MASTER, TAG_VAVG, MPI_COMM_WORLD);

								MPI_Tensor3x3_Double snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Fpavgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FPAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Favgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_FAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Pavgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_PAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Strainavgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_EPSAVG, MPI_COMM_WORLD);

								snd3x3 = MPI_Tensor3x3_Double(avg.at(j).Cauchyavgrve);
								MPI_Send( &snd3x3, 1, MPI_Tensor3x3_Double_Type, MASTER, TAG_CAUAVG, MPI_COMM_WORLD);

								MPI_VonMises_Double sndvm = MPI_VonMises_Double(avg.at(j).Equivavgrve);
								MPI_Send( &sndvm, 1, MPI_VonMises_Double_Type, MASTER, TAG_VMISES, MPI_COMM_WORLD);
							}
						}
						//me is neither MASTER nor whom, i.e. not responsible for I/O for this incr
					} //tasks for non-master processes defined

					//##MK::this barrier not necessary only for DEBUGGING !
					//MPI_Barrier(MPI_COMM_WORLD);
				} //processing of desired targetincr performed
			}

			incr++;
			targetincr += Settings::IncrementOffset;
			if ( targetincr > Settings::IncrementLast )
				targetincr = Settings::IncrementLast; //thus enforcing a stop criterion, as incr keeps increasing...
		} //find targetincrements in steps of current loadcase
	} //find targetincrements in next loadcase

	if ( get_myrank() == MASTER ) {
		//all results now on master, so pipe to file
		ofstream flowcurve_sampled;

		//each converged DAMASK spectral solver increment represent one point on the simulated flow curve
		string fn = get_prefix() + ".Incr.All.Flowcurve.csv";
		flowcurve_sampled.open( fn.c_str(), ofstream::out | std::ofstream::trunc );  //MK::thereby discarding data in already existent files of the same name!

		if ( flowcurve_sampled.is_open() == true ) {
			diary whatwasdone = diary(Settings::SimID);
			flowcurve_sampled << setprecision(32);
			flowcurve_sampled << whatwasdone << endl;

			//Origin inspired three-column header with column name tag, units, and proposal explanation
			flowcurve_sampled << "LoadCaseID;LoadCaseIncr;IncrementID;Volume;";
			flowcurve_sampled << "Fp11;Fp12;Fp13;Fp21;Fp22;Fp23;Fp31;Fp32;Fp33;";
			flowcurve_sampled << "F11;F12;F13;F21;F22;F23;F31;F32;F33;";
			flowcurve_sampled << "P11;P12;P13;P21;P22;P23;P31;P32;P33;";
			flowcurve_sampled << "lnV11;lnV12;lnV13;lnV21;lnV22;lnV23;lnV31;lnV32;lnV33;";
			flowcurve_sampled << "Cauchy11;Cauchy12;Cauchy13;Cauchy21;Cauchy22;Cauchy23;Cauchy31;Cauchy32;Cauchy33;";
			flowcurve_sampled << "vMisesEqvTrueStrain;vMisesEqCauchyStress\n";

			for ( unsigned int tincr = 0; tincr < iobuffer.size(); ++tincr) {
				flowcurve_sampled << iobuffer.at(tincr).loadcaseID << ";";
				flowcurve_sampled << iobuffer.at(tincr).localincrID << ";";
				flowcurve_sampled << iobuffer.at(tincr).globalincrID << ";";
				flowcurve_sampled << iobuffer.at(tincr).Vtotal << ";";

				flowcurve_sampled << iobuffer.at(tincr).Fpavgrve.a11 << ";" << iobuffer.at(tincr).Fpavgrve.a12 << ";" << iobuffer.at(tincr).Fpavgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Fpavgrve.a21 << ";" << iobuffer.at(tincr).Fpavgrve.a22 << ";" << iobuffer.at(tincr).Fpavgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Fpavgrve.a31 << ";" << iobuffer.at(tincr).Fpavgrve.a32 << ";" << iobuffer.at(tincr).Fpavgrve.a33 << ";";

				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a11 << ";" << iobuffer.at(tincr).Favgrve.a12 << ";" << iobuffer.at(tincr).Favgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a21 << ";" << iobuffer.at(tincr).Favgrve.a22 << ";" << iobuffer.at(tincr).Favgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Favgrve.a31 << ";" << iobuffer.at(tincr).Favgrve.a32 << ";" << iobuffer.at(tincr).Favgrve.a33 << ";";

				flowcurve_sampled << iobuffer.at(tincr).Pavgrve.a11 << ";" << iobuffer.at(tincr).Pavgrve.a12 << ";" << iobuffer.at(tincr).Pavgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Pavgrve.a21 << ";" << iobuffer.at(tincr).Pavgrve.a22 << ";" << iobuffer.at(tincr).Pavgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Pavgrve.a31 << ";" << iobuffer.at(tincr).Pavgrve.a32 << ";" << iobuffer.at(tincr).Pavgrve.a33 << ";";

				flowcurve_sampled << iobuffer.at(tincr).Strainavgrve.a11 << ";" << iobuffer.at(tincr).Strainavgrve.a12 << ";" << iobuffer.at(tincr).Strainavgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Strainavgrve.a21 << ";" << iobuffer.at(tincr).Strainavgrve.a22 << ";" << iobuffer.at(tincr).Strainavgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Strainavgrve.a31 << ";" << iobuffer.at(tincr).Strainavgrve.a32 << ";" << iobuffer.at(tincr).Strainavgrve.a33 << ";";

				flowcurve_sampled << iobuffer.at(tincr).Cauchyavgrve.a11 << ";" << iobuffer.at(tincr).Cauchyavgrve.a12 << ";" << iobuffer.at(tincr).Cauchyavgrve.a13 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Cauchyavgrve.a21 << ";" << iobuffer.at(tincr).Cauchyavgrve.a22 << ";" << iobuffer.at(tincr).Cauchyavgrve.a23 << ";";
				flowcurve_sampled << iobuffer.at(tincr).Cauchyavgrve.a31 << ";" << iobuffer.at(tincr).Cauchyavgrve.a32 << ";" << iobuffer.at(tincr).Cauchyavgrve.a33 << ";";

				flowcurve_sampled << iobuffer.at(tincr).Equivavgrve.vMisesEquivStrain << ";";
				flowcurve_sampled << iobuffer.at(tincr).Equivavgrve.vMisesEquivStress << "\n";
			}

			flowcurve_sampled.flush();
			flowcurve_sampled.close();
		}
		else {
			complaining("Flowcurve could not be written to file", myrank, 0);
		}
	} //master finished processing file

	//MK::necessary, whenever attempting MPI parallel computations requiring syncing thereafter!
	MPI_Barrier(MPI_COMM_WORLD);
	double toc = MPI_Wtime();
	tictoc.push_back(plog(tic, toc, "IO::ReportingFlowcurve"));
}



void specOutIncr::spit_profiling()
{
	//##MK::further optimization aand convenience tasks: bundle all in one file, incr ID and so forth

	//##MK::suboptimal... one file per rank
	string fn = get_prefix() +  ".Incr.All.WallClock.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << setprecision(18);
		csvlog << "What;WallClock\n";
		csvlog<< ";s\n";
		csvlog << "What;MPI_Wtime\n";

		for (size_t i = 0; i < tictoc.size(); ++i) {
			csvlog << tictoc.at(i).get_what() << ";" << tictoc.at(i).get_dt() << endl;
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local processing files", get_myrank(), 0);
	}
}


void specOutIncr::debug_signed_distance_function()
{
	grGeomHdl* debugger = NULL;
	try { debugger = new grGeomHdl; }
	catch (bad_alloc &debugcroak) { return; }

	//odd-edge length cubic voxel cubus, voxel center of gravity of central voxel at 0.0, 0.0, 0.0
	debugger->owner = this; //MK::here this is necessary to get address of calling class object
	debugger->cgid = 0;
	debugger->nuip = 0;
	debugger->localbox = aabb3d();
	debugger->grainfence = aabb3d();

	vxlgrd testgrid;
	size_t oddnumber = 21; //101
	testgrid.nx = oddnumber;
	testgrid.ny = testgrid.nx;
	testgrid.nz = testgrid.nx;
	testgrid.nxy = testgrid.nx * testgrid.ny;
	testgrid.nxyz = testgrid.nx * testgrid.ny * testgrid.nz;
	testgrid.dcell = static_cast<real_xyz>(0.5)* (static_cast<real_xyz>(1.0) / static_cast<real_xyz>(oddnumber));
	real_xyz halfdistance = (static_cast<real_xyz>(0.5) + static_cast<real_xyz>((oddnumber-1)/2)) * testgrid.dcell;
	real_xyz imi = static_cast<real_xyz>(-1.0) * halfdistance;
	real_xyz imx = static_cast<real_xyz>(+1.0) * halfdistance;
	testgrid.origin_cntum = p3d( imi, imi, imi );
	testgrid.origin_discr_own = vxl(0, 0, 0);
	testgrid.origin_discr_lnk = vxl(0, 0, 0);
	testgrid.box = aabb3d( imi, imx, imi, imx, imi, imx);

	debugger->localgrid = testgrid;
	debugger->healthy = true;

	real_xyz testradius =  testgrid.dcell * static_cast<real_xyz>(5);
	debugger->debug_sdf_sphere_seed( testradius );
	debugger->debug_sdf_sphere_guess();
	debugger->debug_sdf_sphere_fsm();
	debugger->debug_sdf_sphere_exact( testradius );
	debugger->debug_sdf_sphere_report();

	if ( debugger != NULL) {
		delete debugger; debugger = NULL;
	}
}
