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

#include "PDT_Datatypes.h"

//definition of << overloading to conventiently report values of user-defined types

ostream& operator << (ostream& in, vxl const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "\n";
	return in;
}

ostream& operator << (ostream& in, p3d const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "\n";
	return in;
}

ostream& operator << (ostream& in, d3d const & val) {
	in << val.dx << ";" << val.dy << ";" << val.dz << "\n";
	return in;
}

ostream& operator << (ostream& in, p3dm1 const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "--" << val.m1 << "\n";
	return in;
}

ostream& operator << (ostream& in, p3dm2 const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "--" << val.m1 << "--" << int(val.m2) << "\n";
	return in;
}

ostream& operator << (ostream& in, p3dm3 const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "--" << val.m1 << ";" << val.m2 << ";" << val.m3 << "\n";
	return in;
}

ostream& operator << (ostream& in, v3x1 const & val) {
	in << val.x << ";" << val.y << ";" << val.z << "\n";
	return in;
}

ostream& operator << (ostream& in, bv3x1 const & val) {
	in << val.a11 << ";" << val.a21 << ";" << val.a31 << "\n";
	return in;
}

ostream& operator << (ostream& in, bv3x3 const & val) {
	in << val.a11 << ";" << val.a12 << ";" << val.a13 << "\n";
	in << val.a21 << ";" << val.a22 << ";" << val.a23 << "\n";
	in << val.a31 << ";" << val.a32 << ";" << val.a33 << "\n";
	return in;
}


void aabb3d::scale()
{
	this->xsz = this->xmx - this->xmi;
	this->ysz = this->ymx - this->ymi;
	this->zsz = this->zmx - this->zmi;
}


p3d aabb3d::center()
{
	real_xyz half = static_cast<real_xyz>(0.5);
	return p3d(
			this->xmi + half*(this->xmx-this->xmi),
			this->ymi + half*(this->ymx-this->ymi),
			this->zmi + half*(this->zmx-this->zmi) );
}


void aabb3d::add_guardzone( const real_xyz guardlen )
{
	this->xmi -= guardlen;
	this->xmx += guardlen;
	this->ymi -= guardlen;
	this->ymx += guardlen;
	this->zmi -= guardlen;
	this->zmx += guardlen;
	this->scale();
}


void aabb3d::potentially_expand(aabb3d const & cand)
{
	//check if candidate aabb3d protrudes out of this aabb3d, if so modify this aabb3d
	if ( cand.xmi <= this->xmi ) 	this->xmi = cand.xmi;
	if ( cand.xmx >= this->xmx ) 	this->xmx = cand.xmx;

	if ( cand.ymi <= this->ymi ) 	this->ymi = cand.ymi;
	if ( cand.ymx >= this->ymx ) 	this->ymx = cand.ymx;

	if ( cand.zmi <= this->zmi ) 	this->zmi = cand.zmi;
	if ( cand.zmx >= this->zmx ) 	this->zmx = cand.zmx;

	this->scale();
}


void aabb3d::potentially_expand(const real_xyz dx, const real_xyz dy, const real_xyz dz)
{
	//check if candidate aabb3d protrudes out of this aabb3d, if so modify this aabb3d
	if ( dx <= this->xmi ) 	this->xmi = dx;
	if ( dx >= this->xmx ) 	this->xmx = dx;

	if ( dy <= this->ymi ) 	this->ymi = dy;
	if ( dy >= this->ymx ) 	this->ymx = dy;

	if ( dz <= this->zmi ) 	this->zmi = dz;
	if ( dz >= this->zmx ) 	this->zmx = dz;

	//this->scale();
}


void aabb3d::blowup_xyz(d3d const & fac)
{
	//stretch this aabb3d maintaining center by factors fac_i >= 1.0 in general dissimilarly per principal direction
	this->scale();
	p3d center_old = this->center();

	real_xyz h = static_cast<real_xyz>(0.5);
	d3d newlen = d3d( fac.dx*h*this->xsz, fac.dy*h*this->ysz, fac.dz*h*this->zsz);

	//push old domains outside, ##MK::potentially may suffer numerical issues...
	this->xmi = center_old.x - newlen.dx;
	this->xmx = center_old.x + newlen.dx;

	this->ymi = center_old.y - newlen.dy;
	this->ymx = center_old.y + newlen.dy;

	this->zmi = center_old.z - newlen.dz;
	this->zmx = center_old.z + newlen.dz;

	this->scale();
	//p3d center_new = this->center();
}


bool aabb3d::lurking_out_of(aabb3d const & ref)
{
	//true if this protrudes even numerically or partially (at least dim) out of the reference
	if ( this->xmi < ref.xmi )	return true;
	if ( this->xmx > ref.xmx )	return true;
	if ( this->ymi < ref.ymi )	return true;
	if ( this->ymx > ref.ymx )	return true;
	if ( this->zmi < ref.zmi )	return true;
	if ( this->zmx > ref.zmx )	return true;
	//not return inside
	return false;
}


bool aabb3d::is_inside(const real_xyz x, const real_xyz y, const real_xyz z)
{
	if ( x >= this->xmi && x <= this->xmx &&
			y >= this->ymi && y <= this->ymx &&
				z >= this->zmi && z <= this->zmx )
		return true;
	else
		return false;
}


std::ostream& operator << (std::ostream& in, aabb3d const & val) {
	in << val.xmi << ";" << val.xmx << "---" << val.xsz << "\n";
	in << val.ymi << ";" << val.ymx << "---" << val.ysz << "\n";
	in << val.zmi << ";" << val.zmx << "---" << val.zsz << "\n";
	return in;
}

std::ostream& operator << (std::ostream& in, hexahedron const & val)
{
	in << val.p1.x << ";" << val.p1.y << ";" << val.p1.z << ";";
	in << val.p2.x << ";" << val.p2.y << ";" << val.p2.z << ";";
	in << val.p3.x << ";" << val.p3.y << ";" << val.p3.z << ";";
	in << val.p4.x << ";" << val.p4.y << ";" << val.p4.z << ";";
	in << val.p5.x << ";" << val.p5.y << ";" << val.p5.z << ";";
	in << val.p6.x << ";" << val.p6.y << ";" << val.p6.z << ";";
	in << val.p7.x << ";" << val.p7.y << ";" << val.p7.z << ";";
	in << val.p8.x << ";" << val.p8.y << ";" << val.p8.z << ";";
	return in;
}


unsigned int sqb::where( const unsigned int _xx, const unsigned int _yy, const unsigned int _zz )
{
	//return implicit 3d coordinate
	return _xx + _yy * this->nx + _zz * this->nxy;
}


ostream& operator << (ostream& in, sqb const & val) {
	in << val.nx << ";" << val.ny << ";" << val.nz << ";" << val.nxy << ";" << val.nxyz << "\n";
	return in;
}


size_t vxlgrd::binning_x( const real_xyz x )
{
	//checks whether coordinate x \in \mathbb{R} lays in any discrete bin along the x direction, rectangular transfer function
	//##MK::requires optimization
	//real_xyz tmpreal = (x - this->box.xmi) / (this->box.xmx - this->box.xmi) * static_cast<real_xyz>(this->nx);
	//size_t tmpszi = ( tmpreal >= 0.0 ) ? static_cast<real_xyz>(tmpreal) : 0;
	//size_t xx = ( tmpszi < this->nx ) ? tmpszi : this->nx-1;

	//assuming implicit origin at 0, 0, 0 in this->origin_discr_own as expected
	real_xyz tmpreal = (x - this->box.xmi) / (this->box.xmx - this->box.xmi) * static_cast<real_xyz>(this->nx);
	if ( tmpreal >= static_cast<real_xyz>(fabs(0.0)) ) { //most likely, test only passed if tmpreal's signbit is set to positive thereby rendering the cast to size_t safe avoiding wraparound
		size_t tmpszi = static_cast<size_t>(tmpreal); //##MK::behavior of -0.0 and =0.0 ? stability?
		if ( tmpszi < this->nx )
			return tmpszi;
	}
	//
	//	else
	return numeric_limits<size_t>::max(); //indicating error because build there is no voxel with coordinate this->ni as ni points to voxel past in positive axis direction
	//similarily nothing in the negative
	//return numeric_limits<size_t>::max();
}


size_t vxlgrd::binning_y( const real_xyz y )
{
	//see vxlgrd::binning_x for comments
	real_xyz tmpreal = (y - this->box.ymi) / (this->box.ymx - this->box.ymi) * static_cast<real_xyz>(this->ny);
	if ( tmpreal >= static_cast<real_xyz>(fabs(0.0)) ) {
		size_t tmpszi = static_cast<size_t>(tmpreal);
		if ( tmpszi < this->ny )
			return tmpszi;
	}
	return numeric_limits<size_t>::max();
}


size_t vxlgrd::binning_z( const real_xyz z )
{
	//see vxlgrd::binning_x for comments
	real_xyz tmpreal = (z - this->box.zmi) / (this->box.zmx - this->box.zmi) * static_cast<real_xyz>(this->nz);
	if ( tmpreal >= fabs(static_cast<real_xyz>(fabs(0.0))) ) {
		size_t tmpszi = static_cast<size_t>(tmpreal);
		if ( tmpszi < this->nz )
			return tmpszi;
	}
	return numeric_limits<size_t>::max();
}


ostream& operator << (ostream& in, vxlgrd const & val) {
	in << "Vxlgrid/nx,y,z " << val.nx << ";" << val.ny << ";" << val.nz << "\n";
	in << "Vxlgrid/nxy/nxyz "<< val.nxy << ";" << val.nxyz << "\n";
	in << "Vxlgrid/dcell " << val.dcell << "\n";
	in << "Vxlgrid/orig continuum " << val.origin_cntum;
	in << "Vxlgrid origin discrete myown " << val.origin_discr_own;
	in << "Vxlgrdid origin discrete linked " << val.origin_discr_lnk;
	return in;
}


ostream& operator << (ostream& in, dist const & val) {
	in << val.d << ";" << val.nbor_eipid << "\n";
	return in;
}

ostream& operator << (ostream& in, spatdist const & val) {
	in << val.d << ";" << val.sval << "\n";
	return in;
}

ostream& operator<<(ostream& in, t3x1 const & val)
{
	in << val.a11 << "," << val.a21 << "," << val.a31 << "\n";
	return in;
}

void t3x3::add( const t3x3 & increase, const real_m33 weight )
{
	this->a11 += weight * increase.a11;
	this->a12 += weight * increase.a12;
	this->a13 += weight * increase.a13;

	this->a21 += weight * increase.a21;
	this->a22 += weight * increase.a22;
	this->a23 += weight * increase.a23;

	this->a31 += weight * increase.a31;
	this->a32 += weight * increase.a32;
	this->a33 += weight * increase.a33;
}

void t3x3::div( const real_m33 divisor )
{
	if (abs(divisor) > EPSILON) {
		this->a11 /= divisor;
		this->a12 /= divisor;
		this->a13 /= divisor;

		this->a21 /= divisor;
		this->a22 /= divisor;
		this->a23 /= divisor;

		this->a31 /= divisor;
		this->a32 /= divisor;
		this->a33 /= divisor;
	}
}

ostream& operator << (ostream& in, t3x3 const & val) {
	in << val.a11 << ";" << val.a12 << ";" << val.a13 << "\n";
	in << val.a21 << ";" << val.a22 << ";" << val.a23 << "\n";
	in << val.a31 << ";" << val.a32 << ";" << val.a33 << "\n";
	return in;
}

ostream& operator << (ostream& in, vMises const & val) {
	in << val.vMisesEquivStrain << " " << val.vMisesEquivStress << " Pa" << "\n";
	return in;
}

ostream& operator << (ostream& in, sbvhrange const & val) {
	in << "[" << val.startidx << ", " << val.pastendidx << "]" << "\n";
	return in;
}

ostream& operator << (ostream& in, whenTaken const & val) {
	in << "LoadcaseID " << val.loadcaseID << ", local increment " << val.localincrID << ", globalincr " << val.globalincrID << "\n";
	return in;
}



ostream& operator<<(ostream& in, voro_p3d const & val)
{
	in << "\n";
	in << val.x << ";" << val.y << ";" << val.z << "\n";
	in << "\n";
	return in;
}


void voro_p3d::normalize()
{
	voro_real len = sqrt(SQR(this->x)+SQR(this->y)+SQR(this->z));
	this->x /= len;
	this->y /= len;
	this->z /= len;
}


voro_real dot( voro_p3d const & a, voro_p3d const & b )
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}


voro_p3d cross( voro_p3d const & a, voro_p3d const & b )
{
	return voro_p3d( 	a.y*b.z - a.z*b.y,
						a.z*b.x - a.x*b.z,
						a.x*b.y - a.y*b.x );
}


ostream& operator<<(ostream& in, voro_p3dm1 const & val)
{
	in << "\n";
	in << val.x << ";" << val.y << ";" << val.z << ";" << val.m << "\n";
	in << "\n";
	return in;
}


ostream& operator<<(ostream& in, voro_p3dm2 const & val)
{
	in << "\n";
	in << val.x << ";" << val.y << ";" << val.z << "\t\tGrainID " << val.m1 << "\t\tUIPID " << val.m2 << "\n";
	in << "\n";
	return in;
}


ostream& operator<<(ostream& in, voro_p3i const & val)
{
	in << "\n";
	in << val.x << ";" << val.y << ";" << val.z << "\n";
	in << "\n";
	return in;
}


voro_p3d voro_tri3d::barycenter()
{
	return voro_p3d( 	(this->x1 + this->x2 + this->x3 ) / static_cast<voro_real>(3.0),
						(this->y1 + this->y2 + this->y3 ) / static_cast<voro_real>(3.0),
						(this->z1 + this->z2 + this->z3 ) / static_cast<voro_real>(3.0) );
}

ostream& operator<<(ostream& in, voro_tri3d const & val)
{
	in << "\n";
	in << val.x1 << ";" << val.y1 << ";" << val.z1 << "\n";
	in << val.x2 << ";" << val.y2 << ";" << val.z2 << "\n";
	in << val.x3 << ";" << val.y3 << ";" << val.z3 << "\n";
	in << "\n";
	return in;
}


void voro_aabb3d::add_epsilon_guard( const voro_real eps )
{
	this->xmi -= eps;
	this->xmx += eps;
	this->ymi -= eps;
	this->ymx += eps;
	this->ymi -= eps;
	this->ymx += eps;
}


bool voro_aabb3d::inside( voro_p3dm1 const & test )
{
	if ( test.x < this->xmi )	return false;
	if ( test.x > this->xmx )	return false;
	if ( test.y < this->ymi )	return false;
	if ( test.y > this->ymx )	return false;
	if ( test.z < this->zmi )	return false;
	if ( test.z > this->zmx )	return false;
	return true;
}


bool voro_aabb3d::inside( voro_p3d const & test )
{
	if ( test.x < this->xmi )	return false;
	if ( test.x > this->xmx )	return false;
	if ( test.y < this->ymi )	return false;
	if ( test.y > this->ymx )	return false;
	if ( test.z < this->zmi )	return false;
	if ( test.z > this->zmx )	return false;
	return true;
}


void voro_aabb3d::possibly_enlarge_me( voro_p3dm1 const & test )
{
	if ( test.x <= this->xmi )
		this->xmi = test.x;
	if ( test.x >= this->xmx )
		this->xmx = test.x;

	if ( test.y <= this->ymi )
		this->ymi = test.y;
	if ( test.y >= this->ymx )
		this->ymx = test.y;

	if ( test.z <= this->zmi )
		this->zmi = test.z;
	if ( test.z >= this->zmx )
		this->zmx = test.z;
}


void voro_aabb3d::possibly_enlarge_me( voro_p3d const & test )
{
	if ( test.x <= this->xmi )
		this->xmi = test.x;
	if ( test.x >= this->xmx )
		this->xmx = test.x;

	if ( test.y <= this->ymi )
		this->ymi = test.y;
	if ( test.y >= this->ymx )
		this->ymx = test.y;

	if ( test.z <= this->zmi )
		this->zmi = test.z;
	if ( test.z >= this->zmx )
		this->zmx = test.z;
}


voro_p3d voro_aabb3d::boxcenter( void )
{
	return voro_p3d( 	0.5*(this->xmi + this->xmx),
						0.5*(this->ymi + this->ymx),
						0.5*(this->zmi + this->zmx)  );
}


voro_p3i voro_aabb3d::blockpartitioning( const voro_uint p_total, const voro_uint p_perblock_target )
{
	voro_p3i out = voro_p3i( 1, 1, 1); //at least one block
	if ( (this->xmx - this->xmi) > EPSILON ) { //##MK::assuming cuboidal box with x direction the largest
		voro_real yrel = (this->ymx - this->ymi) / (this->xmx - this->xmi);
		voro_real zrel = (this->zmx - this->zmi) / (this->xmx - this->xmi);
		voro_real ix = pow( (static_cast<voro_real>(p_total) /
				static_cast<voro_real>(p_perblock_target)) / (yrel * zrel), (1.0/3.0) );
		out.x = ( static_cast<voro_int>(floor(ix)) > 0 ) ? static_cast<voro_int>(floor(ix)) : 1; //ceil and floor(x+0.5) ?
		voro_real iy = yrel * ix;
		out.y = ( static_cast<voro_int>(floor(iy)) > 0 ) ? static_cast<voro_int>(floor(iy)) : 1;
		voro_real iz = zrel * ix;
		out.z = ( static_cast<voro_int>(floor(iz)) > 0 ) ? static_cast<voro_int>(floor(iz)) : 1;
	}
	return out;
}


ostream& operator<<(ostream& in, voro_aabb3d const & val)
{
	in << "\n";
	in << "Axis-aligned bounding box bounds" << "\n";
	in << val.xmi << ";" << val.xmx << "\n";
	in << val.ymi << ";" << val.ymx << "\n";
	in << val.zmi << ";" << val.zmx << "\n";
	in << "\n";
	return in;
}


