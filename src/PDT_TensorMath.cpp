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


#include "PDT_TensorMath.h"


bv3x3 leftmult( t3x3 const & defgrad, bv3x3 const & bvecs )
{
	//left multiply RVE initial configuration base vector with deformation gradient tensor
	//get column vector of RVE base vectors in deformed configuration,
	//this transforms base vectors initially aligned parallel to edges of hexahedron to triclinic cell

	//##MK:: implement mind cast from real_m33 to real_xyz
	//##MK:: implement mind cast from real_m33 to real_xyz
	//##MK:: implement mind cast from real_m33 to real_xyz
	//##MK:: implement mind cast from real_m33 to real_xyz
	bv3x3 F = bv3x3( 	static_cast<real_xyz>(defgrad.a11), static_cast<real_xyz>(defgrad.a12), static_cast<real_xyz>(defgrad.a13),
						static_cast<real_xyz>(defgrad.a21), static_cast<real_xyz>(defgrad.a22), static_cast<real_xyz>(defgrad.a23),
						static_cast<real_xyz>(defgrad.a31), static_cast<real_xyz>(defgrad.a32), static_cast<real_xyz>(defgrad.a33) );

	return bv3x3( 	F.a11*bvecs.a11 + F.a12*bvecs.a21 + F.a13*bvecs.a31,
					F.a11*bvecs.a12 + F.a12*bvecs.a22 + F.a13*bvecs.a32,
					F.a11*bvecs.a13 + F.a12*bvecs.a23 + F.a13*bvecs.a33,

					F.a21*bvecs.a11 + F.a22*bvecs.a21 + F.a23*bvecs.a31,
					F.a21*bvecs.a12 + F.a22*bvecs.a22 + F.a23*bvecs.a32,
					F.a21*bvecs.a13 + F.a22*bvecs.a23 + F.a23*bvecs.a33,

					F.a31*bvecs.a11 + F.a32*bvecs.a21 + F.a33*bvecs.a31,
					F.a31*bvecs.a12 + F.a32*bvecs.a22 + F.a33*bvecs.a32,
					F.a31*bvecs.a13 + F.a32*bvecs.a23 + F.a33*bvecs.a33 );
}


t3x3 transpose(t3x3 const &in)
{
	return t3x3(	in.a11, in.a21, in.a31,
					in.a12, in.a22, in.a32,
					in.a13, in.a23, in.a33);
}


t3x3 diag(t3x1 const &in)
{
	return t3x3(	in.a11, static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), in.a21, static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), in.a31);
}


real_m33 det(t3x3 const &in)
{
	return (in.a11*in.a22*in.a33 +
			in.a12*in.a23*in.a31 +
			in.a13*in.a21*in.a32 -
			in.a13*in.a22*in.a31 -
			in.a12*in.a21*in.a33 -
			in.a11*in.a23*in.a32);
}


real_m33 dot(t3x1 const &in1, t3x1 const &in2)
{
	//column vector dot product
	return(in1.a11*in2.a11 + in1.a21*in2.a21 + in1.a31*in2.a31);
}


void killnoise(t3x3 &m, real_m33 thrshld)
{
	//set second-order rank tensor values to zero for numerical conditioning
	m.a11 = (fabs(m.a11) > thrshld) ? m.a11 : static_cast<real_m33>(0.0);
	m.a12 = (fabs(m.a12) > thrshld) ? m.a12 : static_cast<real_m33>(0.0);
	m.a13 = (fabs(m.a13) > thrshld) ? m.a13 : static_cast<real_m33>(0.0);
	m.a21 = (fabs(m.a21) > thrshld) ? m.a21 : static_cast<real_m33>(0.0);
	m.a22 = (fabs(m.a22) > thrshld) ? m.a22 : static_cast<real_m33>(0.0);
	m.a23 = (fabs(m.a23) > thrshld) ? m.a23 : static_cast<real_m33>(0.0);
	m.a31 = (fabs(m.a31) > thrshld) ? m.a31 : static_cast<real_m33>(0.0);
	m.a32 = (fabs(m.a32) > thrshld) ? m.a32 : static_cast<real_m33>(0.0);
	m.a33 = (fabs(m.a33) > thrshld) ? m.a33 : static_cast<real_m33>(0.0);
}


real_m33 trace(t3x3 const &in)
{
	return (in.a11+in.a22+in.a33);
}


void eye(t3x3 &out)
{
	out.a11 = static_cast<real_m33>(1.0);
	out.a12 = static_cast<real_m33>(0.0);
	out.a13 = static_cast<real_m33>(0.0);

	out.a21 = static_cast<real_m33>(0.0);
	out.a22 = static_cast<real_m33>(1.0);
	out.a23 = static_cast<real_m33>(0.0);

	out.a31 = static_cast<real_m33>(0.0);
	out.a32 = static_cast<real_m33>(0.0);
	out.a33 = static_cast<real_m33>(1.0);

/*	return t3x3(	static_cast<real_m33>(1.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(1.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(1.0) );*/
}

t3x3 eye(void)
{
	return t3x3(	static_cast<real_m33>(1.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(1.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(1.0) );
}


void zerot3x1(t3x1 &out)
{
	out.a11 = static_cast<real_m33>(0.0);
	out.a21 = static_cast<real_m33>(0.0);
	out.a31 = static_cast<real_m33>(0.0);

	/*return t3x1(	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0) );*/
}

t3x1 zerot3x1(void)
{
	return t3x1(	static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0) );
}

t3x3 failt3x3(void)
{
#ifdef SINGLE_PRECISION
	return( t3x3( 	numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max(),
					numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max(),
					numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()	) );
#else
	return( t3x3(	numeric_limits<double>::max(), numeric_limits<float>::max(), numeric_limits<float>::max(),
					numeric_limits<double>::max(), numeric_limits<float>::max(), numeric_limits<float>::max(),
					numeric_limits<double>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()  ) );
#endif
}

void zerot3x3(t3x3 &out)
{
	out.a11 = static_cast<real_m33>(0.0);
	out.a12 = static_cast<real_m33>(0.0);
	out.a13 = static_cast<real_m33>(0.0);

	out.a21 = static_cast<real_m33>(0.0);
	out.a22 = static_cast<real_m33>(0.0);
	out.a23 = static_cast<real_m33>(0.0);

	out.a31 = static_cast<real_m33>(0.0);
	out.a32 = static_cast<real_m33>(0.0);
	out.a33 = static_cast<real_m33>(0.0);

	/*return t3x3( 	static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0)  );*/
}


t3x3 zerot3x3(void)
{
	return t3x3( 	static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0),
					static_cast<real_m33>(0.0), static_cast<real_m33>(0.0), static_cast<real_m33>(0.0)  );
}


t3x3 add(t3x3 const &in1, t3x3 const &in2)
{
	//add second-order rank tensor component-wise C = A+B
	return t3x3( 	in1.a11+in2.a11, in1.a12+in2.a12, in1.a13+in2.a13,
					in1.a21+in2.a21, in1.a22+in2.a22, in1.a23+in2.a23,
					in1.a31+in2.a31, in1.a32+in2.a32, in1.a33+in2.a33);
}


t3x3 substrct(t3x3 const &in1, t3x3 const &in2)
{
	return t3x3( 	in1.a11-in2.a11, in1.a12-in2.a12, in1.a13-in2.a13,
					in1.a21-in2.a21, in1.a22-in2.a22, in1.a23-in2.a23,
					in1.a31-in2.a31, in1.a32-in2.a32, in1.a33-in2.a33);
}


/*
void log( t3x3 &inout)
{
	//c++11 log both float and double
	//numpy conformant natural logarithm for all C = log(A)
	inout.a11 = log(inout.a11);	inout.a12 = log(inout.a12);	inout.a13 = log(inout.a13);
	inout.a21 = log(inout.a21);	inout.a22 = log(inout.a22);	inout.a23 = log(inout.a23);
	inout.a31 = log(inout.a31);	inout.a32 = log(inout.a32);	inout.a33 = log(inout.a33);
}
*/


t3x1 log(t3x1 const &in)
{
	return t3x1( 	log(in.a11),
					log(in.a21),
					log(in.a31) );
}


real_m33 sum(t3x3 const &in)
{
	return (	in.a11+in.a12+in.a13 +
				in.a21+in.a22+in.a23 +
				in.a31+in.a32+in.a33	);
}


t3x3 mult(real_m33 sc, t3x3 const &in)
{
	return t3x3(	sc*in.a11, sc*in.a12, sc*in.a13,
					sc*in.a21, sc*in.a22, sc*in.a23,
					sc*in.a31, sc*in.a32, sc*in.a33 );
}


t3x3 dyadic(t3x3 const &in1, t3x3 const &in2)
{
	//python a*b
	return t3x3(	in1.a11*in2.a11, in1.a12*in2.a12, in1.a13*in2.a13,
					in1.a21*in2.a21, in1.a22*in2.a22, in1.a23*in2.a23,
					in1.a31*in2.a31, in1.a32*in2.a32, in1.a33*in2.a33 );
}


t3x3 dot(t3x3 const &in1, t3x3 const &in2)
{
	//python dot(a,b) for a, b general second-order rank tensor
	return t3x3(	in1.a11*in2.a11 + in1.a12*in2.a21 + in1.a13*in2.a31,
					in1.a11*in2.a12 + in1.a12*in2.a22 + in1.a13*in2.a32,
					in1.a11*in2.a13 + in1.a12*in2.a23 + in1.a13*in2.a33,

					in1.a21*in2.a11 + in1.a22*in2.a21 + in1.a23*in2.a31,
					in1.a21*in2.a12 + in1.a22*in2.a22 + in1.a23*in2.a32,
					in1.a21*in2.a13 + in1.a22*in2.a23 + in1.a23*in2.a33,

					in1.a31*in2.a11 + in1.a32*in2.a21 + in1.a33*in2.a31,
					in1.a31*in2.a12 + in1.a32*in2.a22 + in1.a33*in2.a32,
					in1.a31*in2.a13 + in1.a32*in2.a23 + in1.a33*in2.a33 );
}


