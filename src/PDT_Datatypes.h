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

#ifndef __PDT_DATATYPES_H__
#define __PDT_DATATYPES_H__

#include "PDT_Performance.h"

struct dbcontrol
{
	unsigned short label;
	bool visited;
	dbcontrol() : label(numeric_limits<unsigned short>::max()), visited(false) {}
	dbcontrol(const unsigned short _lbl, const bool _seen) :
		label(_lbl), visited(_seen) {}
};


struct lcasemeta
{
	unsigned int freqs;
	double times;
	unsigned int logscales;
	unsigned int nincr;
	lcasemeta() :
		freqs(1), times(numeric_limits<double>::max()), logscales(0), nincr(0) {}
	lcasemeta(const unsigned int _f, const double _t, const unsigned int _lg,
			const unsigned int _ni) : freqs(_f), times(_t), logscales(_lg), nincr(_ni) {}
	~lcasemeta() {}
};

struct vxl
{
	size_t x;
	size_t y;
	size_t z;
	vxl() : x(numeric_limits<size_t>::max()),
				y(numeric_limits<size_t>::max()),
					z(numeric_limits<size_t>::max()) {}
	vxl(const size_t _x, const size_t _y, const size_t _z) :
		x(_x), y(_y), z(_z) {}
	~vxl(){}
};

std::ostream& operator << (std::ostream& in, vxl const & val);

struct p3d
{
	real_xyz x;					//a point in 3d space
	real_xyz y;
	real_xyz z;
	p3d() : x(0.0), y(0.0), z(0.0) {}
	p3d(const real_xyz _x, const real_xyz _y, const real_xyz _z) :
		x(_x), y(_y), z(_z) {}
	~p3d(){}
};

std::ostream& operator << (std::ostream& in, p3d const & val);

struct d3d
{
	real_xyz dx;					//a point displacement in 3d space
	real_xyz dy;
	real_xyz dz;
	d3d() : dx(0.0), dy(0.0), dz(0.0) {}
	d3d(const real_xyz _dx, const real_xyz _dy, const real_xyz _dz) :
		dx(_dx), dy(_dy), dz(_dz) {}
	~d3d(){}
};

std::ostream& operator << (std::ostream& in, d3d const & val);


struct p3dm1
{
	real_xyz x;					//a point in 3d space with three property marks m
	real_xyz y;
	real_xyz z;
	unsigned int m1;

	p3dm1() : x(0.0), y(0.0), z(0.0), m1(std::numeric_limits<unsigned int>::max()) {}
	p3dm1(const real_xyz _x, const real_xyz _y, const real_xyz _z,
			const unsigned int _m1) : x(_x), y(_y), z(_z), m1(_m1) {}
	~p3dm1(){}
};

std::ostream& operator << (std::ostream& in, p3dm1 const & val);


struct p3dm2
{
	real_xyz x;					//a point in 3d space with two property marks m
	real_xyz y;
	real_xyz z;
	unsigned int m1;
	unsigned char m2;

	p3dm2() : x(0.0), y(0.0), z(0.0),
			m1(std::numeric_limits<unsigned int>::max()), m2(std::numeric_limits<unsigned char>::max()) {}
	p3dm2(const p3d _xyz, const unsigned int _m1, const unsigned char _m2) :
				x(_xyz.x), y(_xyz.y), z(_xyz.z), m1(_m1), m2(_m2) {}
	p3dm2(const real_xyz _x, const real_xyz _y, const real_xyz _z,
			const unsigned int _m1, const unsigned char _m2) :
				x(_x), y(_y), z(_z), m1(_m1), m2(_m2) {}
};

std::ostream& operator << (std::ostream& in, p3dm2 const & val);


struct p3dm3
{
	real_xyz x;					//a point in 3d space with three property marks m
	real_xyz y;
	real_xyz z;
	unsigned int m1;
	unsigned int m2;

	unsigned char m3;

	p3dm3() : x(0.0), y(0.0), z(0.0),
			m1(std::numeric_limits<unsigned int>::max()),
			m2(std::numeric_limits<unsigned int>::max()),
			m3(std::numeric_limits<unsigned char>::max())
			{}
	p3dm3(const real_xyz _x, const real_xyz _y, const real_xyz _z,
			const unsigned int _m1, const unsigned int _m2, const unsigned char _m3) :
		x(_x), y(_y), z(_z), m1(_m1), m2(_m2), m3(_m3) {}
	~p3dm3(){}
};

std::ostream& operator << (std::ostream& in, p3dm3 const & val);

struct v3x1
{
	real_xyz x;					//a column vector in 3d
	real_xyz y;
	real_xyz z;
	v3x1() : x(0.0), y(0.0), z(0.0) {}
	v3x1(const real_xyz _x, const real_xyz _y, const real_xyz _z) :
		x(_x), y(_y), z(_z) {}
	~v3x1(){}
};

std::ostream& operator << (std::ostream& in, v3x1 const & val);

struct bv3x1
{
	real_xyz a11;			//a column spatial coordinate system base vector
	real_xyz a21;
	real_xyz a31;
	bv3x1() : 	a11(static_cast<real_xyz>(0.0)),
				a21(static_cast<real_xyz>(0.0)),
				a31(static_cast<real_xyz>(0.0)) {} //initialize to neutral column vector
	bv3x1(const real_xyz _a11, const real_xyz _a21, const real_xyz _a31) :
				a11(_a11),
				a21(_a21),
				a31(_a31) {}
	~bv3x1(){}
};

std::ostream& operator << (std::ostream& in, bv3x1 const & val);

struct bv3x3
{
	real_xyz a11;				//RVE base column vector collection, 1st column || 1. axis, 2nd column || 2. axis, 3rd column || 3.
	real_xyz a12;
	real_xyz a13;
	real_xyz a21;
	real_xyz a22;
	real_xyz a23;
	real_xyz a31;
	real_xyz a32;
	real_xyz a33;
	bv3x3() :	a11(1.0), a12(0.0), a13(0.0),
				a21(0.0), a22(1.0), a23(0.0),
				a31(0.0), a32(0.0), a33(1.0) {}
	bv3x3(	const real_xyz _a11, const real_xyz _a12, const real_xyz _a13,
			const real_xyz _a21, const real_xyz _a22, const real_xyz _a23,
			const real_xyz _a31, const real_xyz _a32, const real_xyz _a33 ) :
				a11(_a11), a12(_a12), a13(_a13),
				a21(_a21), a22(_a22), a23(_a23),
				a31(_a31), a32(_a32), a33(_a33) {}
	~bv3x3(){}
};

std::ostream& operator << (std::ostream& in, bv3x3 const & val);

struct aabb3d
{
	real_xyz xmi;				//an axis-aligned bounding box
	real_xyz xmx;
	real_xyz ymi;
	real_xyz ymx;
	real_xyz zmi;
	real_xyz zmx;
	real_xyz xsz;
	real_xyz ysz;
	real_xyz zsz;
	aabb3d() : 	xmi(std::numeric_limits<real_xyz>::max()),
				xmx(std::numeric_limits<real_xyz>::lowest()),
				ymi(std::numeric_limits<real_xyz>::max()),
				ymx(std::numeric_limits<real_xyz>::lowest()),
				zmi(std::numeric_limits<real_xyz>::max()),
				zmx(std::numeric_limits<real_xyz>::lowest()),
				xsz(0.0), ysz(0.0), zsz(0.0) {}
	aabb3d( const real_xyz _xmi, const real_xyz _xmx,
			const real_xyz _ymi, const real_xyz _ymx,
			const real_xyz _zmi, const real_xyz _zmx ) :
				xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), zmi(_zmi), zmx(_zmx),
				xsz(_xmx - _xmi), ysz(_ymx - _ymi), zsz(_zmx - _zmi) {}
	~aabb3d(){}
	void scale();
	p3d center();
	void add_guardzone( const real_xyz guardlen );
	void potentially_expand(aabb3d const & cand);
	void potentially_expand(const real_xyz dx, const real_xyz dy, const real_xyz dz);
	void blowup_xyz(d3d const & factor);
	bool lurking_out_of(aabb3d const & ref);
	bool is_inside(const real_xyz x, const real_xyz y, const real_xyz z);
};

std::ostream& operator << (std::ostream& in, aabb3d const & val);

struct hexahedron
{
	p3d p1;
	p3d p2;
	p3d p3;
	p3d p4;
	p3d p5;
	p3d p6;
	p3d p7;
	p3d p8;
	hexahedron() : p1(p3d()),  p2(p3d()),  p3(p3d()),  p4(p3d()),  p5(p3d()),  p6(p3d()),  p7(p3d()),  p8(p3d()) {}
	hexahedron(const p3d _p1, const p3d _p2, const p3d _p3, const p3d _p4, const p3d _p5, const p3d _p6, const p3d _p7, const p3d _p8 ) :
		p1(_p1), p2(_p2), p3(_p3), p4(_p4), p5(_p5), p6(_p6), p7(_p7), p8(_p8) {}
};

std::ostream& operator << (std::ostream& in, hexahedron const & val);

struct sqb
{
	unsigned int nx;			//metadata to a spatial decomposition of a cuboidal 3d domain
	unsigned int ny;			//how many subdomains along ny
	unsigned int nz;

	unsigned int nxy;			//x+y*nx+z*nxy implicit indexing
	unsigned int nxyz;

	sqb() : nx(1), ny(1), nz(1), nxy(1), nxyz(1) {}
	sqb(const unsigned int _nx, const unsigned int _ny, const unsigned int _nz) :
		nx(_nx), ny(_ny), nz(_nz), nxy(_nx*_ny), nxyz(_nx*_ny*_nz) {}
	~sqb(){}
};

std::ostream& operator << (std::ostream& in, sqb const & val);


struct vxlgrd
{
	size_t nx;				//edge length of the axis-aligned voxel container
	size_t ny;
	size_t nz;
	size_t nxy;				//x + y*nx + z*nx*ny implict addressing
	size_t nxyz;

	real_xyz dcell;			//generalized distance single voxel edge represents, measure of cubic discretization
	p3d origin_cntum;		//real continuous space coordinate of voxelgrid origin
	vxl origin_discr_own;	//discrete own origin
	vxl origin_discr_lnk;	//my own origin vxl has which ID in a linked coordinate system(i.e. the global voxel grid)
	aabb3d box;

	vxlgrd() : nx(0), ny(0), nz(0), nxy(0), nxyz(0),
			dcell(numeric_limits<real_xyz>::min()),
			origin_cntum(p3d()), origin_discr_own(vxl()), origin_discr_lnk(vxl()),
			box(aabb3d()) {}
	vxlgrd(const size_t _nx, const size_t _ny, const size_t _nz,
			const real_xyz _dc, const p3d _oc, const vxl _odown, const vxl _odlnk, const aabb3d _b ) :
			nx(_nx), ny(_ny), nz(_nz), nxy(_nx*_ny), nxyz(_nx*_ny*_nz),
			dcell(_dc), origin_cntum(_oc), origin_discr_own(_odown), origin_discr_lnk(_odlnk), box(_b) {}
	~vxlgrd(){}
	size_t binning_x( const real_xyz x );
	size_t binning_y( const real_xyz y );
	size_t binning_z( const real_xyz z );
};

std::ostream& operator << (std::ostream& in, vxlgrd const & val);


struct dist
{
	real_xyz d;					//distance between two element ip coordinates
	unsigned int nbor_eipid;	//the element id of the neighbor
	dist() :
		d(numeric_limits<real_xyz>::max()),
		nbor_eipid(numeric_limits<unsigned int>::max()) {}
	dist(const real_xyz _d, const unsigned int _id) :
		d(_d), nbor_eipid(_id) {}
	~dist() {}
};


struct nbp3d
{
	real_xyz x;
	real_xyz y;
	real_xyz z;
	unsigned int uipid;			//which index in RVE1 original uip point cloud
	unsigned char representative;  //MK::one ip is constructed from only one uip, the latter has a grainID because grain Recon took care of periodicity
								//hence uipid refers to particular grain for which we will know, if recon was successful, whether the ip is part of the representative dbscan cluster or not

	nbp3d() : x(0.0), y(0.0), z(0.0),
		uipid(numeric_limits<unsigned int>::max()), representative(static_cast<unsigned char>(NO)) {}
	nbp3d(const real_xyz _x, const real_xyz _y, const real_xyz _z, const unsigned int _uip, const unsigned char _rep) :
		x(_x), y(_y), z(_z), uipid(_uip), representative(_rep) {}
	~nbp3d(){}
};


std::ostream& operator << (std::ostream& in, dist const & val);

struct spatdist
{
	real_xyz d;					//distance between location of integration point me and nearest integration point in distorted grid between along which a certain critical value is exceeded
	real_xyz sval;				//state variable value at integration point me
	spatdist() : d(numeric_limits<real_xyz>::max()),
			sval(numeric_limits<real_xyz>::max()) {}
	spatdist( const real_xyz _d, const real_xyz _sval ) :
		d(_d), sval(_sval) {}
	~spatdist() {}
};

std::ostream& operator << (std::ostream& in, spatdist const & val);


struct t3x1
{
	real_m33 a11;			//a column vector
	real_m33 a21;
	real_m33 a31;
	t3x1() : 	a11(static_cast<real_m33>(0.0)),
				a21(static_cast<real_m33>(0.0)),
				a31(static_cast<real_m33>(0.0)) {} //initialize to neutral column vector
	t3x1(const real_m33 _a11, const real_m33 _a21, const real_m33 _a31) :
				a11(_a11),
				a21(_a21),
				a31(_a31) {}
	~t3x1(){}
};

ostream& operator<<(ostream& in, t3x1 const & val);

struct t3x3
{
	real_m33 a11;				//a second order rank tensor with row-column indexing
	real_m33 a12;
	real_m33 a13;
	real_m33 a21;
	real_m33 a22;
	real_m33 a23;
	real_m33 a31;
	real_m33 a32;
	real_m33 a33;
	t3x3() :	a11(1.0), a12(0.0), a13(0.0),
				a21(0.0), a22(1.0), a23(0.0),
				a31(0.0), a32(0.0), a33(1.0) {}	//initialize to identity tensor
	t3x3( const real_m33* matrix3x3 ) :
				a11(matrix3x3[0]), a12(matrix3x3[1]), a13(matrix3x3[2]),
				a21(matrix3x3[3]), a22(matrix3x3[4]), a23(matrix3x3[5]),
				a31(matrix3x3[6]), a32(matrix3x3[7]), a33(matrix3x3[8]) {}
	t3x3(	const real_m33 _a11, const real_m33 _a12, const real_m33 _a13,
			const real_m33 _a21, const real_m33 _a22, const real_m33 _a23,
			const real_m33 _a31, const real_m33 _a32, const real_m33 _a33 ) :
				a11(_a11), a12(_a12), a13(_a13),
				a21(_a21), a22(_a22), a23(_a23),
				a31(_a31), a32(_a32), a33(_a33) {}
	~t3x3(){};
	void add( const t3x3 & increase, const real_m33 weight );
	void div( const real_m33 divisor );
};

std::ostream& operator << (std::ostream& in, t3x3 const & val);

struct vMises
{
	real_m33 vMisesEquivStrain;
	real_m33 vMisesEquivStress;	
	vMises() :
		vMisesEquivStrain(numeric_limits<real_xyz>::max()),
		vMisesEquivStress(numeric_limits<real_xyz>::max())  {}
	vMises(const real_m33 _epsilon, const real_m33 _sigma) :
		vMisesEquivStrain(_epsilon),
		vMisesEquivStress(_sigma) {}
	~vMises(){}
};

std::ostream& operator << (std::ostream& in, vMises const & val);


struct sbvhrange
{
	size_t startidx;
	size_t pastendidx;
	sbvhrange() : startidx(0), pastendidx(0) {}
	sbvhrange(const size_t _sidx, const size_t _peidx) :
		startidx(_sidx), pastendidx(_peidx) {}
	~sbvhrange(){}
};

std::ostream& operator << (std::ostream& in, sbvhrange const & val);

#endif
