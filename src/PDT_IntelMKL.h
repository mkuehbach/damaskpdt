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

#ifndef __PDT_INTELMKL_H__
#define __PDT_INTELMKL_H__


#include "PDT_VoroComposer.h"


#define INTELMKL_EXISTENT

#ifdef INTELMKL_EXISTENT
	//utilize Intel MKL library
	#include "mkl.h"
	//#include "mkl_dfti.h"
#endif

//singular value decomposition of general matrix
bool svd(t3x3 const &in, t3x3 &U, t3x1 &S, t3x3 &Vh);

//inverse of a general matrix
bool inv(t3x3 const &in, t3x3 &out);

//eigen decomposition
bool eig(t3x3 const &in, t3x1 &e_real, t3x1 &e_img, t3x3 &evr_real, t3x3 &evr_img);

bool computeEigDecompStretch2Strain( t3x3 const & in, const unsigned int kind, t3x3 & out);

bool computeStrainTensor( t3x3 const & in, t3x3 &eps);
bool computeStrainTensor2( t3x3 const & in,
		const unsigned int straintensor,
		const unsigned int strainkind, t3x3 & eps);



bool computeCauchy( t3x3 const & F, t3x3 const & P, t3x3 &cauchy);

real_m33 computeMises( t3x3 const & tensor, bool isstresstensor );


aabb3d get_corners( vector<real_xyz> &cx, vector<real_xyz> &cy, vector<real_xyz> &cz );


class rfftn
{
public:
	rfftn();
	~rfftn();

	void init(unsigned int const * ngridd);
	void fill(vector<t3x3> const & F, const unsigned int c);
	void forwardFFT();
	vector<char> fftTemp;	//forward transformation results carrier for real and imaginary parts
	vector<double> m_input; //input signal, i.e. tensor component

	DFTI_CONFIG_VALUE precision;
	DFTI_DESCRIPTOR_HANDLE m_f_handle;

	MKL_LONG NX;  //##MK::initialize and check if grid is zero
	MKL_LONG NY;
	MKL_LONG NZ;
	MKL_LONG NXYZ;
	MKL_LONG NI;
	MKL_LONG NJ;
	MKL_LONG NK;
	MKL_LONG NJK;
	MKL_LONG NIJK;
	MKL_LONG dimensions[3];
	MKL_LONG input_strides[4];
	MKL_LONG output_strides[4];
	MKL_Complex16* fftTempP;
};


class irfftn
{
public:
	irfftn();
	~irfftn();

	void init( unsigned int const * nfft, unsigned int const * ngridd );
	void fill(vector<double> const & in_real, vector<double> const & in_imag, const unsigned int cm);
	void backwardFFTandNormalize();

	vector<char> fftTemp;		//inverse transform input, i.e. integrated components in Fourier space for one coordinate direction
	vector<double> m_output; 	//inverse transform result

	DFTI_CONFIG_VALUE precision;
	DFTI_DESCRIPTOR_HANDLE m_b_handle;

	MKL_LONG NFFT; //##MK::initialize and check if grid is zero
	MKL_LONG NI;
	MKL_LONG NJ;
	MKL_LONG NK;
	MKL_LONG NJK;
	MKL_LONG NIJK;
	MKL_LONG NX;
	MKL_LONG NY;
	MKL_LONG NZ;
	MKL_LONG NXY;
	MKL_LONG NXYZ;
	MKL_LONG dimensions[3];
	MKL_LONG input_strides[4];
	MKL_LONG output_strides[4];
	MKL_Complex16* fftTempP;
};


unsigned int fourier_transform_defgradient( vector<t3x3> const & defgrad,
		unsigned int const * ngr, double const * nsz, vector<rfftn*> & ffts );


unsigned int displacement_fluctuations( vector<rfftn*> const & ffts,
		unsigned int const * ngr, double const * nsz, vector<irfftn*> & iffts,
		vector<d3d> & flucts );


void displacement_average( vector<rfftn*> const & ffts,
		unsigned int const * ngr, double const * nsz, vector<d3d> & flucts );


#endif
