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


#include "PDT_OriMath.h"

std::ostream& operator << (std::ostream& in, quat const & val)
{
	in << val.q0 << ";" << val.q1 << "i;" << val.q2 << "j;" << val.q3 << "k" << "\n";
	return in;
}


/*
void inverse_quaternion( real_ori* q )
{

//	//for general quaternion its inverse is its conjugate normalized by amount

//	real_ori norm = static_cast<real_ori>( sqrt(SQR[0]+SQR[1]+SQR[2]+SQR[3]) );
//	if ( norm > EPSILON ) { //exists
//		q[0] /= norm;
//		q[1] *= -1.0;		q[1] /= norm;
//		q[2] *= -1.0;		q[2] /= norm;
//		q[3] *= -1.0;		q[3] /= norm;
//	}
//	else {
//		//##MK::implement complaint, when inverse at least numerically not safely exists
//	}


	//for unit quaternion norm = 1.0 so
	conjugate_quaternion( q );
}


void conjugate_quaternion( real_ori* q )
{
	//q[0] remains as is
	q[1] *= -1.0;
	q[2] *= -1.0;
	q[3] *= -1.0;
}


void multiply_quaternion( real_ori* p, real_ori* q, real_ori* pq)
{
    pq[0] =  p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]; //DAMASK::ok same as in Consistent rotations tutorial, 2015 and Grimmer, 1974
    pq[1] =  p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2];
    pq[2] =  p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1];
    pq[3] =  p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0];
}


void active2passive( real_ori* q )
{
	//transforming between rotation convention
	inverse_quaternion( q );
}


void passive2active( real_ori* q )
{
	//transforming between rotation convention
	inverse_quaternion( q );
}


real_ori disorientation_angle_fcc( real_ori* p, real_ori* q )
{
	//see mentioned Grimmer, 1974 Acta Cryst A reference in declaration for details
	//passive orientation quaternions

	real_ori qconj[4];
	qconj[0] = q[0];			//in place inverse of unit quaternion q^-1
	qconj[1] = -1.0*q[1];
	qconj[2] = -1.0*q[2];
	qconj[3] = -1.0*q[3];

	real_ori qcand[4] = { 1.0, 0.0, 0.0, 0.0 };
	multiply_quaternion( p, qconj, qcand ); //misorientation, passive convention pq^-1

	//fcc fundamental quaternions
	real_ori qcandsymm[24];

	//taking absolute values as we are only interested in disorientation angle not the axis

	qcandsymm[ 0] = fabs(qcand[0]); //Grimmer 7a
	qcandsymm[ 1] = fabs(qcand[1]);
	qcandsymm[ 2] = fabs(qcand[2]);
	qcandsymm[ 3] = fabs(qcand[3]);

	real_ori _sqrt2 = static_cast<real_ori>(1.0) / sqrt(static_cast<real_ori>(2.0));
	qcandsymm[ 4] = fabs(_sqrt2*(qcand[0] + qcand[1])); //Grimmer 7b
	qcandsymm[ 5] = fabs(_sqrt2*(qcand[0] - qcand[1]));
	qcandsymm[ 6] = fabs(_sqrt2*(qcand[2] + qcand[3]));
	qcandsymm[ 7] = fabs(_sqrt2*(qcand[2] - qcand[3]));

	qcandsymm[ 8] = fabs(_sqrt2*(qcand[0] + qcand[2])); //Grimmer 7c
	qcandsymm[ 9] = fabs(_sqrt2*(qcand[0] - qcand[2]));
	qcandsymm[10] = fabs(_sqrt2*(qcand[1] + qcand[3]));
	qcandsymm[11] = fabs(_sqrt2*(qcand[1] - qcand[3]));

	qcandsymm[12] = fabs(_sqrt2*(qcand[0] + qcand[3])); //Grimmer 7d
	qcandsymm[13] = fabs(_sqrt2*(qcand[0] - qcand[3]));
	qcandsymm[14] = fabs(_sqrt2*(qcand[1] + qcand[2]));
	qcandsymm[15] = fabs(_sqrt2*(qcand[1] - qcand[2]));

	real_ori half = static_cast<real_ori>(0.5);
	qcandsymm[16] = fabs(half*(qcand[0] + qcand[1] + qcand[2] + qcand[3])); //Grimmer 7e
	qcandsymm[17] = fabs(half*(qcand[0] + qcand[1] - qcand[2] - qcand[3]));
	qcandsymm[18] = fabs(half*(qcand[0] - qcand[1] + qcand[2] - qcand[3]));
	qcandsymm[19] = fabs(half*(qcand[0] - qcand[1] - qcand[2] + qcand[3]));

	qcandsymm[20] = fabs(half*(qcand[0] + qcand[1] + qcand[2] - qcand[3])); //Grimmer 7f
	qcandsymm[21] = fabs(half*(qcand[0] + qcand[1] - qcand[2] + qcand[3]));
	qcandsymm[22] = fabs(half*(qcand[0] - qcand[1] + qcand[2] + qcand[3]));
	qcandsymm[23] = fabs(half*(qcand[0] - qcand[1] - qcand[2] - qcand[3]));

	real_ori maximum = 0.0;
	for (unsigned int i = 0; i < 24; ++i) {
		if ( qcandsymm[i] >= maximum) //##MK::potentially case < maximum hit more likely, then better if/else construct to profit from more successful branch prediction
			maximum = qcandsymm[i];
	}

	if ( maximum <= static_cast<real_ori>(1.0) ) { //p and q are the same
		return (2.0*acos(maximum)); //angle in radiant
		//cout << "Nan occurring in disori " << maximum << "-->" << p[0] << ";" << p[1] << ";" << p[2] << ";" << p[3] << "\n";
		//cout << "Nan occurring in disori " << maximum << "-->" << qconj[0] << ";" << qconj[1] << ";" << qconj[2] << ";" << qconj[3] << "\n";
		//cout << "Nan occurring in disori qcand component " << qcand[0] << ";" << qcand[1] << ";" << qcand[2] << ";" << qcand[3] << "\n";
	}
	//##MK::p and q likely the same, implicit else
	return static_cast<real_ori>(0.0);
}





real_ori disorientation_q0_fcc( real_ori* p, real_ori* q )
{
	//see mentioned Grimmer, 1974 Acta Cryst A reference in declaration for details
	//passive orientation quaternions

	real_ori qconj[4];
	qconj[0] = q[0];			//in place inverse of unit quaternion q^-1
	qconj[1] = -1.0*q[1];
	qconj[2] = -1.0*q[2];
	qconj[3] = -1.0*q[3];

	real_ori qcand[4] = { 1.0, 0.0, 0.0, 0.0 };
	multiply_quaternion( p, qconj, qcand ); //misorientation, passive convention pq^-1

	//fcc fundamental quaternions
	real_ori qcandsymm[24];

	//taking absolute values as we are only interested in disorientation angle not the axis

	qcandsymm[ 0] = fabs(qcand[0]); //Grimmer 7a
	qcandsymm[ 1] = fabs(qcand[1]);
	qcandsymm[ 2] = fabs(qcand[2]);
	qcandsymm[ 3] = fabs(qcand[3]);

	real_ori _sqrt2 = static_cast<real_ori>(1.0) / sqrt(static_cast<real_ori>(2.0));
	qcandsymm[ 4] = fabs(_sqrt2*(qcand[0] + qcand[1])); //Grimmer 7b
	qcandsymm[ 5] = fabs(_sqrt2*(qcand[0] - qcand[1]));
	qcandsymm[ 6] = fabs(_sqrt2*(qcand[2] + qcand[3]));
	qcandsymm[ 7] = fabs(_sqrt2*(qcand[2] - qcand[3]));

	qcandsymm[ 8] = fabs(_sqrt2*(qcand[0] + qcand[2])); //Grimmer 7c
	qcandsymm[ 9] = fabs(_sqrt2*(qcand[0] - qcand[2]));
	qcandsymm[10] = fabs(_sqrt2*(qcand[1] + qcand[3]));
	qcandsymm[11] = fabs(_sqrt2*(qcand[1] - qcand[3]));

	qcandsymm[12] = fabs(_sqrt2*(qcand[0] + qcand[3])); //Grimmer 7d
	qcandsymm[13] = fabs(_sqrt2*(qcand[0] - qcand[3]));
	qcandsymm[14] = fabs(_sqrt2*(qcand[1] + qcand[2]));
	qcandsymm[15] = fabs(_sqrt2*(qcand[1] - qcand[2]));

	real_ori half = static_cast<real_ori>(0.5);
	qcandsymm[16] = fabs(half*(qcand[0] + qcand[1] + qcand[2] + qcand[3])); //Grimmer 7e
	qcandsymm[17] = fabs(half*(qcand[0] + qcand[1] - qcand[2] - qcand[3]));
	qcandsymm[18] = fabs(half*(qcand[0] - qcand[1] + qcand[2] - qcand[3]));
	qcandsymm[19] = fabs(half*(qcand[0] - qcand[1] - qcand[2] + qcand[3]));

	qcandsymm[20] = fabs(half*(qcand[0] + qcand[1] + qcand[2] - qcand[3])); //Grimmer 7f
	qcandsymm[21] = fabs(half*(qcand[0] + qcand[1] - qcand[2] + qcand[3]));
	qcandsymm[22] = fabs(half*(qcand[0] - qcand[1] + qcand[2] + qcand[3]));
	qcandsymm[23] = fabs(half*(qcand[0] - qcand[1] - qcand[2] - qcand[3]));

	real_ori maximum = 0.0;
	for (unsigned int i = 0; i < 24; ++i) {
		if ( qcandsymm[i] >= maximum) //##MK::potentially case < maximum hit more likely, then better if/else construct to profit from more successful branch prediction
			maximum = qcandsymm[i];
	}

	return maximum;
}


quatcloud quaternioncloud_characterize( vector<quat> const & q )
{
	//implement mean quaternion of given quaternion set q

	//##########
	//##########

	//implement deviation from this mean ##MK::as of now fcc only and not considering Pantleon, 2005, Mat Sc A, etc.
	//i.e. this "a classical Krieger-Lassen mean"

	//##########
	//##########

	//##MK::implement further inferential statistics of the SO3 quaternion cloud as for instance
	//reported in F. Bachmann, R. Hielscher, P. E. Jupp, W. Pantleon, H. Schaeben, E. Wegert
	//Journal of Applied Crystallography, 43, 2010, 1338-1355 doi:10.1107/S002188981003027X

	return quatcloud();
}
*/



ori_real om3x3::det()
{
	ori_real row1 = +1.f*this->a11 * ((this->a22*this->a33) - (this->a32*this->a23));
	ori_real row2 = +1.f*this->a21 * ((this->a12*this->a33) - (this->a32*this->a13));
	ori_real row3 = +1.f*this->a31 * ((this->a12*this->a23) - (this->a22*this->a13));
	return row1 - row2 + row3;
}


squat om3x3::om2qu()
{
	ori_real q0 = ORI_HALF * sqrt(ORI_ONE + this->a11 + this->a22 + this->a33);
	ori_real q1 = ORI_HALF * P_IJK * sqrt(ORI_ONE + this->a11 - this->a22 - this->a33);
	ori_real q2 = ORI_HALF * P_IJK * sqrt(ORI_ONE - this->a11 + this->a22 - this->a33);
	ori_real q3 = ORI_HALF * P_IJK * sqrt(ORI_ONE - this->a11 - this->a22 + this->a33);
	if ( this->a32 < this->a23 )
		q1 = -q1;
	if ( this->a13 < this->a31 )
		q2 = -q2;
	if ( this->a21 < this->a12 )
		q3 = -q3;

	squat res = squat( q0, q1, q2, q3 );
	res.normalize();
	return res;
}


bunge om3x3::om2eu()
{
	bunge res = bunge();

	if ( fabs(this->a33) <= (ORI_ONE - ORI_EPSILON) ) { //rotation matrix component <= 1
		ori_real zeta = ORI_ONE / sqrt( ORI_ONE - SQR(this->a33) );

		res.phi1 = atan2( this->a31 * zeta, -this->a32 * zeta );
		res.Phi = acos( this->a33 );
		res.phi2 = atan2( this->a13 * zeta, this->a23 * zeta );
	}
	else {
		res.phi1 = atan2( this->a12, this->a11 );
		res.Phi = static_cast<ori_real>(MYPI) * ORI_HALF * (ORI_ONE - this->a33);
		res.phi2 = ORI_ZERO;
	}
	return res;
}


ostream& operator<<(ostream& in, om3x3 const & val)
{
	in << "\n";
	in << val.a11 << "\t\t" << val.a12 << "\t\t" << val.a13 << "\n";
	in << val.a21 << "\t\t" << val.a22 << "\t\t" << val.a23 << "\n";
	in << val.a31 << "\t\t" << val.a32 << "\t\t" << val.a33 << "\n";
	in << "\n";
	return in;
}


squat::squat( const ori_real _q0, const ori_real _q1, const ori_real _q2, const ori_real _q3 )
{
	this->q0 = _q0;
	this->q1 = _q1;
	this->q2 = _q2;
	this->q3 = _q3;
	if ( this->q0 < ORI_ZERO ) { //reverse sign for consistence, northern hemisphere
		this->q0 = -_q0;
		this->q1 = -_q1;
		this->q2 = -_q2;
		this->q3 = -_q3;
	}
}


squat::squat( const ori_real X0, const ori_real X1, const ori_real X2 )
{
	//MK: Quaternion algebra: q = -q define an equivalent rotation. for unit quaternions ||q|| = 1 so the inverse q^-1 = q*/||q||^2 simplifies to = q* with q* the conjugated q0 - (q1,q2,q3)
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134, mind quaternion order of Shoemake w  + i*v + j*x + k*z <=> 0 +  1 2 3 with
	ori_real r1 = sqrt(ORI_ONE - X0);
	ori_real r2 = sqrt(X0);
	ori_real theta1 = ORI_TWO * static_cast<ori_real>(MYPI) * X1;
	ori_real theta2 = ORI_TWO * static_cast<ori_real>(MYPI) * X2;

	this->q0 = r1*sin(theta1);
	this->q1 = r1*cos(theta1);
	this->q2 = r2*sin(theta2);
	this->q3 = r2*cos(theta2); //w,v,x,z

	if ( this->q0 < ORI_ZERO ) { //reverse sign to Northern hemisphere
		this->q0 = -this->q0;
		this->q1 = -this->q1;
		this->q2 = -this->q2;
		this->q3 = -this->q3;
	}
	//this->normalize();
}


void squat::normalize()
{
	ori_real len = sqrt(SQR(this->q0)+SQR(this->q1)+SQR(this->q2)+SQR(this->q3));
	this->q0 /= len;
	this->q1 /= len;
	this->q2 /= len;
	this->q3 /= len;
}


squat squat::invert()
{
	//for unit quaternion norm = 1.0 so
	return conjugate();
}


squat squat::conjugate()
{
	return squat( this->q0, -this->q1, -this->q2, -this->q3 );
}


om3x3 squat::qu2om()
{
	ori_real qbar = (this->q0*this->q0) -
			((this->q1*this->q1) + (this->q2*this->q2) + (this->q3*this->q3));

	om3x3 res = om3x3(
			qbar + (ORI_TWO*this->q1*this->q1),
			ORI_TWO*((this->q1*this->q2) - (P_IJK*this->q0*this->q3)),
			ORI_TWO*((this->q1*this->q3) + (P_IJK*this->q0*this->q2)),

			ORI_TWO*((this->q1*this->q2) + (P_IJK*this->q0*this->q3)),
			qbar + (ORI_TWO*this->q2*this->q2),
			ORI_TWO*((this->q2*this->q3) - (P_IJK*this->q0*this->q1)),

			ORI_TWO*((this->q1*this->q3) - (P_IJK*this->q0*this->q2)),
			ORI_TWO*((this->q2*this->q3) + (P_IJK*this->q0*this->q1)),
			qbar + (ORI_TWO*this->q3*this->q3)  );

	if ( fabs(res.det() - ORI_ONE) < ORI_EPSILON ) {
		return res;
	}
	else {
		cerr << "qu2om " << this->q0 << ";" << this->q1 << ";" << this->q2 << ";" << this->q3 << " results in improper rotation matrix " << res.det() << "\n";
		return res;
	}
}


bunge squat::qu2eu()
{
	//qu2eu
	ori_real q03 = SQR(this->q0) + SQR(this->q3); //>= 0.0
	ori_real q12 = SQR(this->q1) + SQR(this->q2); //>= 0.0
	ori_real chi = sqrt(q03*q12); //>= 0.0

	//gimbal lock
	bunge res = bunge();
	if ( chi >= ORI_EPSILON ) { //&& q12 >= ORI_EPSILON && q12 >= ORI_EPSILON ) { //most likely case
		res.phi1 = atan2( (this->q1*this->q3 - P_IJK*this->q0*this->q2) / chi,
				(-P_IJK*this->q0*this->q1 - this->q2*this->q3) / chi );
		res.Phi = atan2( ORI_TWO*chi, q03 - q12 );
		res.phi2 = atan2( (P_IJK*this->q0*this->q2 + this->q1*this->q3) / chi,
				(this->q2*this->q3 - P_IJK*this->q0*this->q1) / chi );
	}
	else {
		res.phi1 = ORI_ZERO; //error values
		res.Phi = ORI_ZERO;
		res.phi2 = ORI_ZERO;

		if ( chi < ORI_EPSILON && q12 < ORI_EPSILON ) {
			res.phi1 = atan2( -ORI_TWO*P_IJK*this->q0*this->q3, (SQR(this->q0)-SQR(this->q3)) );
			res.Phi = ORI_ZERO;
			res.phi2 = ORI_ZERO;
		}
		if ( chi < EPSILON && q03 < EPSILON ) {
			res.phi1 = atan2( +ORI_TWO*this->q1*this->q2, (SQR(this->q1)-SQR(this->q2)) );
			res.Phi = static_cast<ori_real>(MYPI);
			res.phi2 = ORI_ZERO;
		}
	}
	return res;
}


ostream& operator<<(ostream& in, squat const & val)
{
	in << "\n";
	in << val.q0 << "\t\t" << val.q1 << "\t\t" << val.q2 << "\t\t" << val.q3 << "\n";
	in << "\n";
	return in;
}


bunge::bunge( const string parseme )
{
	this->phi1 = ORI_ZERO;
	this->Phi = ORI_ZERO;
	this->phi2 = ORI_ZERO;

	stringstream parsethis;
	string datapiece;
	parsethis << parseme;

	getline( parsethis, datapiece, ';');
#ifdef ORI_SINGLE_PRECISION
	this->phi1 = DEGREE2RADIANT(stof(datapiece));
#else
	this->phi1 = DEGREE2RADIANT(stod(datapiece));
#endif

	getline( parsethis, datapiece, ';');
#ifdef ORI_SINGLE_PRECISION
	this->Phi = DEGREE2RADIANT(stof(datapiece));
#else
	this->Phi = DEGREE2RADIANT(stod(datapiece));
#endif

	getline( parsethis, datapiece, ';');
#ifdef ORI_SINGLE_PRECISION
	this->phi2 = DEGREE2RADIANT(stof(datapiece));
#else
	this->phi2 = DEGREE2RADIANT(stod(datapiece));
#endif
}


squat bunge::eu2qu()
{
	ori_real sigma = ORI_HALF*(this->phi1 + this->phi2);
	ori_real delta = ORI_HALF*(this->phi1 - this->phi2);
	ori_real c = cos(ORI_HALF*this->Phi);
	ori_real s = sin(ORI_HALF*this->Phi);

	ori_real test_northernhemi_q0 = c * cos(sigma);
	if ( test_northernhemi_q0 < ORI_ZERO ) { //reverse sign of entire quaternion to bring to Northern hemisphere
		return squat( c*cos(sigma), +P_IJK*s*cos(delta), +P_IJK*s*sin(delta), +P_IJK*c*sin(sigma) );
	}
	else {
		return squat( +c*cos(sigma), -P_IJK*s*cos(delta), -P_IJK*s*sin(delta), -P_IJK*c*sin(sigma) );
	}
}


om3x3 bunge::eu2om()
{
    ori_real c1 = cos(this->phi1);
    ori_real s1 = sin(this->phi1);
    ori_real c  = cos(this->Phi);
    ori_real s  = sin(this->Phi);
    ori_real c2 = cos(this->phi2);
    ori_real s2 = sin(this->phi2);
    
    return om3x3(    +c1*c2 - s1*c*s2,       +s1*c2 + c1*c*s2,      +s*s2,
                    -c1*s2 - s1*c*c2,       -s1*s2 + c1*c*c2,      +s*c2,
                    +s1*s,                  -c1*s,                 +c  );
}
    

ostream& operator<<(ostream& in, bunge const & val)
{
	in << "\n";
	in << val.phi1 << "\t\t" << val.Phi << "\t\t" << val.phi2 << "\n";
	in << "\n";
	return in;
}


squat multiply_quaternion( squat const & p, squat const & q)
{
	//pq = p*q
    return squat( 	p.q0*q.q0 - p.q1*q.q1 - p.q2*q.q2 - p.q3*q.q3,
    		 		p.q0*q.q1 + p.q1*q.q0 + p.q2*q.q3 - p.q3*q.q2,
    				p.q0*q.q2 - p.q1*q.q3 + p.q2*q.q0 + p.q3*q.q1,
					p.q0*q.q3 + p.q1*q.q2 - p.q2*q.q1 + p.q3*q.q0   );  //DAMASK::ok same as in Consistent rotations tutorial, 2015 and Grimmer, 1974
}


pair<ori_real,ori_real> disorientation_angle_fcc_grimmer( squat const & qcand )
{
	//fcc fundamental quaternions
	ori_real qcandsymm[24];

	//taking absolute values as we are only interested in disorientation angle not the axis
	qcandsymm[ 0] = fabs(qcand.q0); //Grimmer 7a
	qcandsymm[ 1] = fabs(qcand.q1);
	qcandsymm[ 2] = fabs(qcand.q2);
	qcandsymm[ 3] = fabs(qcand.q3);

	ori_real _sqrt2 = ORI_ONE / sqrt(static_cast<ori_real>(ORI_TWO));
	qcandsymm[ 4] = fabs(_sqrt2*(qcand.q0 + qcand.q1)); //Grimmer 7b
	qcandsymm[ 5] = fabs(_sqrt2*(qcand.q0 - qcand.q1));
	qcandsymm[ 6] = fabs(_sqrt2*(qcand.q2 + qcand.q3));
	qcandsymm[ 7] = fabs(_sqrt2*(qcand.q2 - qcand.q3));

	qcandsymm[ 8] = fabs(_sqrt2*(qcand.q0 + qcand.q2)); //Grimmer 7c
	qcandsymm[ 9] = fabs(_sqrt2*(qcand.q0 - qcand.q2));
	qcandsymm[10] = fabs(_sqrt2*(qcand.q1 + qcand.q3));
	qcandsymm[11] = fabs(_sqrt2*(qcand.q1 - qcand.q3));

	qcandsymm[12] = fabs(_sqrt2*(qcand.q0 + qcand.q3)); //Grimmer 7d
	qcandsymm[13] = fabs(_sqrt2*(qcand.q0 - qcand.q3));
	qcandsymm[14] = fabs(_sqrt2*(qcand.q1 + qcand.q2));
	qcandsymm[15] = fabs(_sqrt2*(qcand.q1 - qcand.q2));

	ori_real half = static_cast<ori_real>(ORI_HALF);
	qcandsymm[16] = fabs(half*(qcand.q0 + qcand.q1 + qcand.q2 + qcand.q3)); //Grimmer 7e
	qcandsymm[17] = fabs(half*(qcand.q0 + qcand.q1 - qcand.q2 - qcand.q3));
	qcandsymm[18] = fabs(half*(qcand.q0 - qcand.q1 + qcand.q2 - qcand.q3));
	qcandsymm[19] = fabs(half*(qcand.q0 - qcand.q1 - qcand.q2 + qcand.q3));

	qcandsymm[20] = fabs(half*(qcand.q0 + qcand.q1 + qcand.q2 - qcand.q3)); //Grimmer 7f
	qcandsymm[21] = fabs(half*(qcand.q0 + qcand.q1 - qcand.q2 + qcand.q3));
	qcandsymm[22] = fabs(half*(qcand.q0 - qcand.q1 + qcand.q2 + qcand.q3));
	qcandsymm[23] = fabs(half*(qcand.q0 - qcand.q1 - qcand.q2 - qcand.q3));

	ori_real maximum = ORI_ZERO;
	for (unsigned int i = 0; i < 24; ++i) {
		if ( qcandsymm[i] >= maximum) //##MK::potentially case < maximum hit more likely, then better if/else construct to profit from more successful branch prediction
			maximum = qcandsymm[i];
	}

	if ( maximum <= (ORI_ONE - ORI_EPSILON) ) { //p and q are not the same
		return make_pair( maximum, ORI_TWO * acos(maximum));
		//return (ORI_TWO * acos(maximum)); //angle in radiant
		//cout << "Nan occurring in disori " << maximum << "-->" << p[0] << ";" << p[1] << ";" << p[2] << ";" << p[3] << "\n";
		//cout << "Nan occurring in disori " << maximum << "-->" << qconj[0] << ";" << qconj[1] << ";" << qconj[2] << ";" << qconj[3] << "\n";
		//cout << "Nan occurring in disori qcand component " << qcand[0] << ";" << qcand[1] << ";" << qcand[2] << ";" << qcand[3] << "\n";
	}
	//maximum == ORI_ONE
	//##MK::p and q likely the same, implicit else
	return make_pair( maximum, ORI_EPSILON );

	//return ORI_EPSILON;
}

