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

#ifndef __PDT_ORIMATH_H__
#define __PDT_ORIMATH_H__


#include "PDT_TensorMath.h"
//#include "PDT_IntelMKL.h"

struct euler
{
	real_ori phi1;
	real_ori Phi;
	real_ori phi2;
	euler() : phi1(0.0), Phi(0.0), phi2(0.0) {}
	euler(const real_ori _e1, const real_ori _e2, const real_ori _e3) :
		phi1(_e1), Phi(_e2), phi2(_e3) {}
	//~euler(){}
};


struct quat
{

	real_ori q0;					//a scalar+iu+jv+kw unit quaternion SO3
	real_ori q1;
	real_ori q2;
	real_ori q3;
	quat() : q0(1.0), q1(0.0), q2(0.0), q3(0.0) {} //identity quaternion, i.e. "no rotation"
	quat(const real_ori _q0, const real_ori _q1, const real_ori _q2, const real_ori _q3) :
		q0(_q0), q1(_q1), q2(_q2), q3(_q3) {}
	//~quat(){}
};


struct quatcloud
{
	real_ori qm0;					//cloud mean in SO3
	real_ori qm1;
	real_ori qm2;
	real_ori qm3;
	real_ori grod;					//mean disorientation angle deviation to cloud mean of cloud elements
	quatcloud() : qm0(1.0), qm1(0.0), qm2(0.0), qm3(0.0), grod(0.0) {}
	quatcloud(const real_ori _qm0, const real_ori _qm1,
			const real_ori _qm2, const real_ori _qm3, const real_ori _grod) :
				qm0(_qm0), qm1(_qm1), qm2(_qm2), qm3(_qm3), grod(_grod) {}
	~quatcloud(){}
};


#define SQR(a)					((a)*(a))
#define CUBE(a)					((a)*(a)*(a))
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)				(((X) > (Y)) ? (X) : (Y))

#define MAXIMUM_MISORI_FCC		(DEGREE2RADIANT(62.8))
#define SYMMETRIES_IN_FCC		(24)
#define EPS_ENVIRONMENT			(EPSILON)

/*
//MK::Quaternion and orientation definitions
 * see: Consistent representations of and conversions between 3D rotations
 * D. Rowenhorst, A. D. Rollett, G. S. Rohrer, M. Groeber, M. Jackson, P. J. Konijnenberg, M. de Graef
 * Modelling and Simulation in Materials Science and Engineering 2015, Vol 23,
 * doi: 10.1088/0965-0393/23/8/083501
 * see in addition: Disorientations and coincidence rotations for cubic lattices
 * H. Grimmer, Acta Crystallographica Section, 1974, A30, 685-688
 * doi: 10.1107/S0567739474001719
 *
 * all quaternions are meant as unit quaternions!
 * Quaternions i^2=j^2=k^2=-1, ij=-ji=k, jk = -kj=i, ki=-ik=j,
 * q := (q_0, \textbf{q}) = (q_0, q_1, q_2, q_3)
 * Fortran-caused index shift though: 0<=>1, 1<=>2, 2<=>3, 3<=>4
 * but keeping the ijk definition w+ix+jy+kz, w<=>q0, x<=>q1, y<=>q2, z<=>q3
 * the comment "DAMASK:ok", means the function is consistent with the one in
 * DAMASK v2.0.1 under the assumption of P_EIJK=+1
 *
 * in every case we follow the rotation conventions by H.-J. Bunge,
 * Texture Analysis in Materials Science 1982
 * Bunge utilizes the 3-1-3 i.e. ZXZ convention...
 * Gimbal lock occurs when \Phi = 0
 * right-handed Cartesian reference frames, passive interpretation of rotations, rotation angle \omega \in [0,\pi]
 * counter-clockwise rotations are positive to an outbound axis
 * the range of the Bunge-Euler angles is \varphi_1 \in [0,2\pi], \Phi \in [0, \pi], \varphi_2 \in [0,2\pi]
*/
#define P_EIJK					(+1.0) 		//MK::do not change, unless reading above reference!

void inverse_quaternion( real_ori* q );
void conjugate_quaternion( real_ori* q );
void multiply_quaternion( real_ori* p, real_ori* q, real_ori* pq);
void active2passive( real_ori* q );
void passive2active( real_ori* q );
real_ori disorientation_angle_fcc( real_ori* p, real_ori* q );
real_ori disorientation_q0_fcc( real_ori* p, real_ori* q );

quatcloud quaternioncloud_characterize( vector<quat> const & q );
/*
//tensor math
bv3x3 leftmult( t3x3 const & defgrad, bv3x3 const & bvecs );

*/

#endif
