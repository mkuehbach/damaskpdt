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
	in << val.q0 << ";" << val.q1 << "i;" << val.q2 << "j;" << val.q3 << "k" << endl;
	return in;
}


void inverse_quaternion( real_ori* q )
{
	/*
	//for general quaternion its inverse is its conjugate normalized by amount

	real_ori norm = static_cast<real_ori>( sqrt(SQR[0]+SQR[1]+SQR[2]+SQR[3]) );
	if ( norm > EPSILON ) { //exists
		q[0] /= norm;
		q[1] *= -1.0;		q[1] /= norm;
		q[2] *= -1.0;		q[2] /= norm;
		q[3] *= -1.0;		q[3] /= norm;
	}
	else {
		//##MK::implement complaint, when inverse at least numerically not safely exists
	}
	*/

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
		//cout << "Nan occurring in disori " << maximum << "-->" << p[0] << ";" << p[1] << ";" << p[2] << ";" << p[3] << endl;
		//cout << "Nan occurring in disori " << maximum << "-->" << qconj[0] << ";" << qconj[1] << ";" << qconj[2] << ";" << qconj[3] << endl;
		//cout << "Nan occurring in disori qcand component " << qcand[0] << ";" << qcand[1] << ";" << qcand[2] << ";" << qcand[3] << endl;
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
