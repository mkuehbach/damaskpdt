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


#include "PDT_IntelMKL.h"

//singular value decomposition of general matrix
bool svd(t3x3 const &in, t3x3 &U, t3x1 &S, t3x3 &Vh)
{
	//https://software.intel.com/en-us/articles/checking-correctness-of-lapack-svd-eigenvalue-and-one-sided-decomposition-routines
	//singular value decomposition of in into [U,S,Vh] consistent with [U,S,Vh] = numpy.linalg(svd(in))
	//https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesvd_ex.c.htm*/

	/* #define M 3
	#define N 3
	#define LDAA M
	#define LDU M
	#define LDVT N */

	U = zerot3x3();
	S = zerot3x1();
	Vh = zerot3x3();

#ifndef INTELMKL_EXISTENT
	#pragma omp critical
	{
		cerr << "IntelMKL commands were not compiled into the source code!" << "\n";
	}
	return false;
#else
	#ifdef SINGLE_PRECISION
		#pragma omp critical
		{
			cerr << "MKL SVD for single precision is not yet implemented!" << "\n";
		}
		return false;
	#else
		int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;
		double wkopt;
		double* work = NULL;
		double s[3], u[3*3], vt[3*3];
		double a[3*3];

		a[0]= static_cast<double>(in.a11);	a[1]= static_cast<double>(in.a12); 	a[2] = static_cast<double>(in.a13);
		a[3]= static_cast<double>(in.a21);	a[4]= static_cast<double>(in.a22); 	a[5] = static_cast<double>(in.a23);
		a[6]= static_cast<double>(in.a31);	a[7]= static_cast<double>(in.a32); 	a[8] = static_cast<double>(in.a33);

	//cout << "DGESVD Example Program Results" << "\n";

		//query and allocate the optimal workspace
		lwork = -1;
		dgesvd( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
		if ( info == 0 ) {

			lwork = static_cast<int>(wkopt);
			try {
				work = new double[lwork*sizeof(double)];
			}
			catch (bad_alloc &croak) {
				#pragma omp critical
				{
					cerr << "MKL SVD allocation of working array failed" << "\n";
				}
				return false;
			}

			dgesvd( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info );

			if ( info == 0 ) { //most likely
				/*
				//##MK::DEBUG
				cout << "MKL SVD converged" << "\n";
				*/

				//U stored column-wise by MKL
				U.a11 = static_cast<real_m33>(u[0*ldu+0]);		U.a12 = static_cast<real_m33>(u[0*ldu+1]); 	U.a13 = static_cast<real_m33>(u[0*ldu+2]);
				U.a21 = static_cast<real_m33>(u[1*ldu+0]);		U.a22 = static_cast<real_m33>(u[1*ldu+1]); 	U.a23 = static_cast<real_m33>(u[1*ldu+2]);
				U.a31 = static_cast<real_m33>(u[2*ldu+0]);		U.a32 = static_cast<real_m33>(u[2*ldu+1]); 	U.a33 = static_cast<real_m33>(u[2*ldu+2]);

				//S singular values in descending order
				S.a11 = static_cast<real_m33>(s[0*1+0]); //##MK::compiler will optimize explicitly written array position indexing away, here just for clarity against Intel example print_matrix
				S.a21 = static_cast<real_m33>(s[0*1+1]);
				S.a31 = static_cast<real_m33>(s[0*1+2]);

				//##MK::Vh stored row-wise by MKL with the C interface!
				Vh.a11 = static_cast<real_m33>(vt[0*ldvt+0]);	Vh.a12 = static_cast<real_m33>(vt[0*ldvt+1]);	Vh.a13 = static_cast<real_m33>(vt[0*ldvt+2]);
				Vh.a21 = static_cast<real_m33>(vt[1*ldvt+0]);	Vh.a22 = static_cast<real_m33>(vt[1*ldvt+1]);	Vh.a23 = static_cast<real_m33>(vt[1*ldvt+2]);
				Vh.a31 = static_cast<real_m33>(vt[2*ldvt+0]);	Vh.a32 = static_cast<real_m33>(vt[2*ldvt+1]);	Vh.a33 = static_cast<real_m33>(vt[2*ldvt+2]);

				/*
				//##MK::DEBUG
				cout << "U" << "\n";
				cout << U.a11 << "," << U.a12 << "," << U.a13 << "\n";
				cout << U.a21 << "," << U.a22 << "," << U.a23 << "\n";
				cout << U.a31 << "," << U.a32 << "," << U.a33 << "\n";

				cout << "S" << "\n";
				cout << S.a11 << "," << S.a21 << "," << S.a31 << "\n";

				cout << "Vh" << "\n";
				cout << Vh.a11 << "," << Vh.a12 << "," << Vh.a13 << "\n";
				cout << Vh.a21 << "," << Vh.a22 << "," << Vh.a23 << "\n";
				cout << Vh.a31 << "," << Vh.a32 << "," << Vh.a33 << "\n";
				*/

				delete [] work;
				return true; //successful
			}
			else {
				#pragma omp critical
				{
					cerr << "MKL SVD did not converge!" << "\n";
				}
				delete [] work;
				return false;
			}
		}
		else { //workspace allocation failed
			#pragma omp critical
			{
				cout << "MKL SVD working array size determination failed" << "\n";
			}
			return false;
		}
	#endif
#endif
}


//inverse of a general matrix
bool inv(t3x3 const &in, t3x3 &out)
{
	//inverse of general 3x3 matrix with real_m33 in
	out = zerot3x3();

#ifndef INTELMKL_EXISTENT
	#pragma omp critical
	{
		cerr << "IntelMKL commands were not compiled into the source code!" << "\n";
	}
	return false;
#else
	#ifdef SINGLE_PRECISION
		#pragma omp critical
		{//##MK::not yet implemented
			cerr << "MKL INV for single precision is not yet implemented!" << "\n";
		}
		return false;
	#else

		int SIZE = 3;
		double A[SIZE*SIZE];
		A[0] = static_cast<double>(in.a11);		A[1] = static_cast<double>(in.a12);		A[2] = static_cast<double>(in.a13);
		A[3] = static_cast<double>(in.a21);		A[4] = static_cast<double>(in.a22);		A[5] = static_cast<double>(in.a23);
		A[6] = static_cast<double>(in.a31);		A[7] = static_cast<double>(in.a32);		A[8] = static_cast<double>(in.a33);

		/*
		//##MK::DEBUG
		cout << "Original" << "\n";
		cout << A[0] << "," << A[1] << "," << A[2] << "\n";
		cout << A[3] << "," << A[4] << "," << A[5] << "\n";
		cout << A[6] << "," << A[7] << "," << A[8] << "\n";
		*/

		int info;
		int lwsp = SIZE;
		int* ipiv = NULL;
		try {
			ipiv = new int[SIZE]; //##MK::needs to be pointer?
		}
		catch (bad_alloc &croak) {
			#pragma omp critical
			{
				cerr << "MKL INV allocation of ipiv array failing" << "\n";
			}
			return false;
		}

		//LU decomposition of general matrix, double precision
		dgetrf(&SIZE, &SIZE, A, &SIZE, ipiv, &info);

		if ( info == 0 ) { //LUsuccess

			double* wsp = NULL;
			try {
				wsp = new double[lwsp*sizeof(double)];
			}
			catch (bad_alloc &croak) {
				cout << "MKL INV allocation of workspace array failing" << "\n";
				return false;
			}

			//inverse of the general matrix given its LU decomposition
			dgetri(&SIZE,A,&SIZE,ipiv,wsp,&lwsp,&info);

			if ( info == 0 ) { //inversion was successful
				/*
				//##MK::DEBUG
				cout << "Inverse" << "\n";
				cout << A[0] << "," << A[1] << "," << A[2] << "\n";
				cout << A[3] << "," << A[4] << "," << A[5] << "\n";
				cout << A[6] << "," << A[7] << "," << A[8] << "\n";
				*/

				delete [] ipiv;
				delete [] wsp;

				out.a11 = static_cast<real_m33>(A[0]);		out.a12 = static_cast<real_m33>(A[1]);		out.a13 = static_cast<real_m33>(A[2]);
				out.a21 = static_cast<real_m33>(A[3]);		out.a22 = static_cast<real_m33>(A[4]);		out.a23 = static_cast<real_m33>(A[5]);
				out.a31 = static_cast<real_m33>(A[6]);		out.a32 = static_cast<real_m33>(A[7]);		out.a33 = static_cast<real_m33>(A[8]);

				return true;
			}
			else {
				#pragma omp critical
				{
					cerr << "MKL INV inversion dgetri did not converge" << "\n";
				}
				delete [] ipiv;
				delete [] wsp;
				return false;
			}
		}
		else {
			#pragma omp critical
			{
				cerr << "MKL INV dgetrf did not converge" << "\n";
			}
			delete [] ipiv;
			return false;
		}
	#endif
#endif
}


bool eig(t3x3 const &in, t3x1 &e_real, t3x1 &e_img, t3x3 &evr_real, t3x3 &evr_img)
{
	//eigenvalues and right handed eigenvector
	e_real = zerot3x1();
	e_img = zerot3x1();
	evr_real = zerot3x3();
	evr_img = zerot3x3();

#ifndef INTELMKL_EXISTENT
	#pragma omp critical
	{
		cerr << "IntelMKL commands were not compiled into the source code!" << "\n";
	}
	return false;
#else
	#ifdef SINGLE_PRECISION
		//##MK::not yet implemented
		#pragma omp critical
		{
			cerr << "MKL EIG for single precision is not yet implemented!" << "\n";
		}
		return false;
	#else
		MKL_INT n = 3, lda = 3, ldvl = 3, ldvr = 3, info;
		/* Local arrays */
		double wr[3], wi[3], vl[3*3], vr[3*3];
		double a[3*3];
		a[0]= static_cast<double>(in.a11);	a[1]= static_cast<double>(in.a12); 	a[2] = static_cast<double>(in.a13);
		a[3]= static_cast<double>(in.a21);	a[4]= static_cast<double>(in.a22); 	a[5] = static_cast<double>(in.a23);
		a[6]= static_cast<double>(in.a31);	a[7]= static_cast<double>(in.a32); 	a[8] = static_cast<double>(in.a33);

		/*
		//##MK::DEBUG
		cout << "in original" << "\n";
		cout << a[0] << "," << a[1] << "," << a[2] << "\n";
		cout << a[3] << "," << a[4] << "," << a[5] << "\n";
		cout << a[6] << "," << a[7] << "," << a[8] << "\n";
		cout << "LAPACKE_dgeev (row-major, high-level) Example Program Results" << "\n";
		*/

		//solve eigen problem
		info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'V', n, a, lda, wr, wi, vl, ldvl, vr, ldvr ); //do not compute left eigenvectors

		if( info == 0 ) { //converged
			e_real.a11 = static_cast<real_m33>(wr[0]);	e_img.a11 = static_cast<real_m33>(wi[0]);//column vector
			e_real.a21 = static_cast<real_m33>(wr[1]);	e_img.a21 = static_cast<real_m33>(wi[1]);
			e_real.a31 = static_cast<real_m33>(wr[2]);	e_img.a31 = static_cast<real_m33>(wi[2]);

			/*
			//##MK::DEBUG
			cout << "three complex eigenvalues" << "\n";
			cout << wr[0] << " / " << wi[0] << "\n";
			cout << wr[1] << " / " << wi[1] << "\n";
			cout << wr[2] << " / " << wi[2] << "\n";
			*/

			//unroll for loop to get eigenvalues for (MKL_INT i = 0; i < n; ++i) { //only possible either all eigenvector real then three reals or one real and one complex conjugate pair
			MKL_INT	i = 0;
			MKL_INT	j = 0;
			if ( wi[j] == static_cast<double>(0.0) ) { //first eigenvalue real?
				evr_real.a11 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a11 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
					evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a12 = static_cast<real_m33>(0.0);		j++;
					if ( wi[j] == static_cast<double>(0.0) ) {
						evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a13 = static_cast<real_m33>(0.0);
					}
					else {
						#pragma omp critical
						{
							cerr << "After two real, eigenvalue inconsistent" << "\n";
						}
						return false;
					}
				}
				else { //...if not it is the complex conjugate pair
					evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a12 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
					evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); 	evr_img.a13 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); 	j+=2;
				}
			}
			else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
				evr_real.a11 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a11 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a12 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); 	j+=2;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a13 = static_cast<real_m33>(0.0);			j++;
				}
				else {
					#pragma omp critical
					{
						cerr << "After conjugate pair, third eigenvalue inconsistent" << "\n";
					}
					return false;
				}
			}
			i = 1;
			j = 0;
			if ( wi[j] == static_cast<double>(0.0) ) { //first eigenvalue real?
				evr_real.a21 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a21 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
					evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a22 = static_cast<real_m33>(0.0);		j++;
					if ( wi[j] == static_cast<double>(0.0) ) {
						evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(0.0);	//j++;
					}
					else {
						#pragma omp critical
						{
							cerr << "After two real, eigenvalue inconsistent" << "\n";
						}
						return false;
					}
				}
				else { //...if not it is the complex conjugate pair
					evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a22 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
					evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
				}
			}
			else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
				evr_real.a21 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a21 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a22 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(0.0);			j++;
				}
				else {
					#pragma omp critical
					{
						cerr << "After conjugate pair, third eigenvalue inconsistent" << "\n";
					}
					return false;
				}
			}
			i = 2;
			j = 0;
			if ( wi[j] == static_cast<double>(0.0) ) { //first eigenvalue real?
				evr_real.a31 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a31 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
					evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a32 = static_cast<real_m33>(0.0);		j++;
					if ( wi[j] == static_cast<double>(0.0) ) {
						evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(0.0);	//j++;
					}
					else {
						#pragma omp critical
						{
							cerr << "After conjugate pair, third eigenvalue inconsistent" << "\n";
						}
						return false;
					}
				}
				else { //...if not it is the complex conjugate pair
					evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a32 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
					evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
				}
			}
			else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
				evr_real.a31 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a31 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a32 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(0.0);			j++;
				}
				else {
					#pragma omp critical
					{
						cerr << "After conjugate pair, third eigenvalue inconsistent" << "\n";
					}
					return false; //cout << "After conjugate pair, third eigenvalue inconsistent" << "\n";
				}
			}

			//##MK::DEBUG
			/*
			cout << "three right-handed eigenvector" << "\n";
			//##order of these?
			cout << "first complex column vector" << "\n";
			cout << evr_real.a11 << " / " << evr_img.a11 << "\n";
			cout << evr_real.a21 << " / " << evr_img.a21 << "\n";
			cout << evr_real.a31 << " / " << evr_img.a31 << "\n";
			cout << "second complex column vector" << "\n";
			cout << evr_real.a12 << " / " << evr_img.a12 << "\n";
			cout << evr_real.a22 << " / " << evr_img.a22 << "\n";
			cout << evr_real.a32 << " / " << evr_img.a32 << "\n";
			cout << "third complex column vector" << "\n";
			cout << evr_real.a13 << " / " << evr_img.a13 << "\n";
			cout << evr_real.a23 << " / " << evr_img.a23 << "\n";
			cout << evr_real.a33 << " / " << evr_img.a33 << "\n";
			*/
			return true;

		}
		else { //not converged
			#pragma omp critical
			{
				cerr << "MKL EIG dgeev did not converge" << "\n";
			}
			return false;
		}
	#endif
#endif
}


/*
bool computeStrainTensor( t3x3 const & in, t3x3 &eps)
{
#ifndef INTELMKL_EXISTENT
	cout << "IntelMKL commands were not compiled into the source code!" << "\n";
	eps = failt3x3();
	return false;
#else
	//start with deformation gradient tensor F as input in and make polar decomposition
	t3x3 U = t3x3();
	t3x1 S = t3x1();
	t3x3 Vh = t3x3();

	if ( svd(in, U, S, Vh) == true ) {
		//get rotation of polar decomposition
		t3x3 R = dot(U,Vh);

		//inverse of this rotation
		t3x3 Rinv = t3x3();
		if ( inv(R, Rinv) == true ) {
			//for theStretch in theStretches 'U', 'V'

			//F = RU so RinvF = RinvRU = IdU
			//help = 'material strains based on right Cauchy--Green deformation, i.e., C and U')
//			t3x3 stretchU = dot(Rinv, in);
//			killnoise( stretchU, static_cast<real_m33>(1.0e-12) );

			//F = VR so FRinv = VRRinv = VId
			//help = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V')
			t3x3 stretchV = dot(in, Rinv);

			killnoise( stretchV, static_cast<real_m33>(1.0e-12) );
			
			//eigen decomposition of stretch tensor
			t3x1 Dr = t3x1(); 	//real part of eigenvalues
			t3x1 Di = t3x1();	//imaginary part
			t3x3 Vr = t3x3();	//real part of eigenvectors
			t3x3 Vi = t3x3();	//imaginary part

			if ( eig( stretchV, Dr, Di, Vr, Vi) == true ) {
			//if ( eig( stretchU, Dr, Di, Vr, Vi) == true ) {
//cout << "Eigendecomposition of stretchU was successful!" << "\n";
cout << "Eigendecomposition of stretchV was successful!" << "\n";

				//flip eigenvalue and eigenvector in case of negative (real part ##MK::is this sufficent?)eigenvalues
				if ( Dr.a11 <= 0.0 ) {
					//cout << "Flipping 1. negative eigenvalue/-vector positive" << "\n";
					Dr.a11 *= -1.0;		Di.a11 *= -1.0;
					Vr.a11 *= -1.0;		Vi.a11 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a21 *= -1.0;		Vi.a21 *= -1.0;
					Vr.a31 *= -1.0;		Vi.a31 *= -1.0;
				}
				if ( Dr.a21 <= 0.0 ) {
					//cout << "Flipping 2. negative eigenvalue/-vector positive" << "\n";
					Dr.a21 *= -1.0;		Di.a21 *= -1.0;
					Vr.a12 *= -1.0;		Vi.a12 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a22 *= -1.0;		Vi.a22 *= -1.0;
					Vr.a32 *= -1.0;		Vi.a32 *= -1.0;
				}
				if ( Dr.a31 <= 0.0 ) {
					//cout << "Flipping 3. negative eigenvalue/-vector positive" << "\n";
					Dr.a31 *= -1.0;		Di.a31 *= -1.0;
					Vr.a13 *= -1.0;		Vi.a13 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a23 *= -1.0;		Vi.a23 *= -1.0;
					Vr.a33 *= -1.0;		Vi.a33 *= -1.0;
				}

				//check for orthogonality
				//##MK::
cout << "##MK::add IntelMKL addStrainTensors add orthogonality check" << "\n";

				//##MK::implement various cases which reference (V,U based) and flavor (ln,Biot,Green)
				//U#ln
				//##MK::work only with real part of eigenvalues and eigenvector
				Dr.a11 = log(Dr.a11);
				Dr.a21 = log(Dr.a21);
				Dr.a31 = log(Dr.a31);

				t3x3 dia = diag(Dr);
				t3x3 Vtrsp = transpose(Vr);

				//eps = (np.dot(V,np.dot(np.diag(d),V.T)).real).reshape(9)
				t3x3 tmp = dot(dia, Vtrsp);

				//report strain tensor
				eps = dot(Vr,tmp);
				return true;
			}
			else {
cout << "Inversion of stretchU or V failed" << "\n";
				eps = zerot3x3();
				return false;
			}
		}
		else {
cout << "MKL addStrainTensor R inversion failed" << "\n";
			eps = zerot3x3();
			return false;
		}
	}
	else {
cout << "MKL addStrainTensor singular value decomposition of input failed" << "\n";
		eps = zerot3x3();
		return false;
	}
#endif
}
*/

//##MK::heavily annotated for debugging
//#define VERBOSE_ADDSTRAINTENSOR
bool computeEigDecompStretch2Strain( t3x3 const & in, const unsigned int kind, t3x3 & out)
{
	//eigen decomposition of stretch tensor in
	t3x1 Dr = t3x1(); 	//real part of eigenvalues
	t3x1 Di = t3x1();	//imaginary part
	t3x3 Vr = t3x3();	//real part of eigenvectors
	t3x3 Vi = t3x3();	//imaginary part

	if ( eig( in, Dr, Di, Vr, Vi) == true ) {
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Dr" << Dr << "\n";
cout << "Vr" << Vr << "\n";
cout << "Eigendecomposition of stretch was successful!" << "\n";
#endif

		//flip eigenvalue and eigenvector in case of negative (real part ##MK::is this sufficent?)eigenvalues
		if ( Dr.a11 <= 0.0 ) {
			//cout << "Flipping 1. negative eigenvalue/-vector positive" << "\n";
			Dr.a11 *= -1.0;		Di.a11 *= -1.0;
			Vr.a11 *= -1.0;		Vi.a11 *= -1.0; //matrix V contains the three column eigenvectors
			Vr.a21 *= -1.0;		Vi.a21 *= -1.0;
			Vr.a31 *= -1.0;		Vi.a31 *= -1.0;
		}
		if ( Dr.a21 <= 0.0 ) {
			//cout << "Flipping 2. negative eigenvalue/-vector positive" << "\n";
			Dr.a21 *= -1.0;		Di.a21 *= -1.0;
			Vr.a12 *= -1.0;		Vi.a12 *= -1.0; //matrix V contains the three column eigenvectors
			Vr.a22 *= -1.0;		Vi.a22 *= -1.0;
			Vr.a32 *= -1.0;		Vi.a32 *= -1.0;
		}
		if ( Dr.a31 <= 0.0 ) {
			//cout << "Flipping 3. negative eigenvalue/-vector positive" << "\n";
			Dr.a31 *= -1.0;		Di.a31 *= -1.0;
			Vr.a13 *= -1.0;		Vi.a13 *= -1.0; //matrix V contains the three column eigenvectors
			Vr.a23 *= -1.0;		Vi.a23 *= -1.0;
			Vr.a33 *= -1.0;		Vi.a33 *= -1.0;
		}
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Dr" << Dr << "\n";
cout << "Vr" << Vr << "\n";
#endif

		//check for orthogonality
		//##MK::
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "##MK::add IntelMKL addStrainTensors add orthogonality check" << "\n";
#endif

		//##MK::implement rest of the possible kinds (Biot)
		switch (kind) {
			case LNSTRAIN:
				//V#ln
				//U#ln
				Dr.a11 = log(Dr.a11);
				Dr.a21 = log(Dr.a21);
				Dr.a31 = log(Dr.a31);
				break;
			//##MK::add BIOTSTRAIN
			case GREENSTRAIN:
				//U#Green
				//##MK:: but not ! V#Green !
				Dr.a11 = (SQR(Dr.a11) - 1.0) * 0.5;
				Dr.a21 = (SQR(Dr.a21) - 1.0) * 0.5;
				Dr.a31 = (SQR(Dr.a31) - 1.0) * 0.5;
				break;
			default:
				out = zerot3x3();
				return false;
		}

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "lnDr" << Dr << "\n";
#endif
		t3x3 dia = diag(Dr);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "diag lnDr" << dia << "\n";
#endif

		t3x3 Vtrsp = transpose(Vr);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Vtrsp " << Vtrsp << "\n";
#endif

		//eps = (np.dot(V,np.dot(np.diag(d),V.T)).real).reshape(9)
		t3x3 tmp = dot(dia, Vtrsp);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "dot(dia,Vtrsp) " << tmp << "\n";
#endif

		//report strain tensor
		out = dot(Vr,tmp);
		return true;
	}
	else {
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Eigendecomposition of stretch failed!" << "\n";
#endif
		out = zerot3x3();
		return false;
	}
}


bool computeStrainTensor( t3x3 const & in, t3x3 & eps)
{
#ifndef INTELMKL_EXISTENT
	#pragma omp critical
	{
		cerr << "IntelMKL commands were not compiled into the source code!" << "\n";
	}
	eps = failt3x3();
	return false;
#else
	//start with deformation gradient tensor F as input in and make polar decomposition
	t3x3 U = t3x3();
	t3x1 S = t3x1();
	t3x3 Vh = t3x3();

	if ( svd(in, U, S, Vh) == true ) {
		//get rotation of polar decomposition
		t3x3 R = dot(U,Vh);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "U" << U << "\n";
cout << "S" << S << "\n";
cout << "Vh" << Vh << "\n";
cout << "R" << R << "\n";
#endif

		//inverse of this rotation
		t3x3 Rinv = t3x3();
		if ( inv(R, Rinv) == true ) {
cout << "Rinv" << Rinv << "\n";
			//for theStretch in theStretches 'U', 'V'

			//F = RU so RinvF = RinvRU = IdU
			//help = 'material strains based on right Cauchy--Green deformation, i.e., C and U')
			t3x3 stretchU = dot(Rinv, in);
			killnoise( stretchU, static_cast<real_m33>(EPSILON) );

			//F = VR so FRinv = VRRinv = VId
			//help = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V')
			t3x3 stretchV = dot(in, Rinv);
			killnoise( stretchV, static_cast<real_m33>(EPSILON) );

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "sU" << stretchU << "\n";
cout << "sV" << stretchV << "\n";
#endif


			//eigen decomposition of stretch tensor
			t3x1 Dr = t3x1(); 	//real part of eigenvalues
			t3x1 Di = t3x1();	//imaginary part
			t3x3 Vr = t3x3();	//real part of eigenvectors
			t3x3 Vi = t3x3();	//imaginary part

			if ( eig( stretchV, Dr, Di, Vr, Vi) == true ) {
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Dr" << Dr << "\n";
cout << "Vr" << Vr << "\n";
#endif

			//if ( eig( stretchU, Dr, Di, Vr, Vi) == true ) {
//cout << "Eigendecomposition of stretchU was successful!" << "\n";
cout << "Eigendecomposition of stretchV was successful!" << "\n";

				//flip eigenvalue and eigenvector in case of negative (real part ##MK::is this sufficent?)eigenvalues
				if ( Dr.a11 <= 0.0 ) {
					//cout << "Flipping 1. negative eigenvalue/-vector positive" << "\n";
					Dr.a11 *= -1.0;		Di.a11 *= -1.0;
					Vr.a11 *= -1.0;		Vi.a11 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a21 *= -1.0;		Vi.a21 *= -1.0;
					Vr.a31 *= -1.0;		Vi.a31 *= -1.0;
				}
				if ( Dr.a21 <= 0.0 ) {
					//cout << "Flipping 2. negative eigenvalue/-vector positive" << "\n";
					Dr.a21 *= -1.0;		Di.a21 *= -1.0;
					Vr.a12 *= -1.0;		Vi.a12 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a22 *= -1.0;		Vi.a22 *= -1.0;
					Vr.a32 *= -1.0;		Vi.a32 *= -1.0;
				}
				if ( Dr.a31 <= 0.0 ) {
					//cout << "Flipping 3. negative eigenvalue/-vector positive" << "\n";
					Dr.a31 *= -1.0;		Di.a31 *= -1.0;
					Vr.a13 *= -1.0;		Vi.a13 *= -1.0; //matrix V contains the three column eigenvectors
					Vr.a23 *= -1.0;		Vi.a23 *= -1.0;
					Vr.a33 *= -1.0;		Vi.a33 *= -1.0;
				}
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Dr" << Dr << "\n";
cout << "Vr" << Vr << "\n";
#endif

				//check for orthogonality
				//##MK::
cout << "##MK::add IntelMKL addStrainTensors add orthogonality check" << "\n";

				//##MK::implement various cases which reference (V,U based) and flavor (ln,Biot,Green)
				//U#ln
				//##MK::work only with real part of eigenvalues and eigenvector
				Dr.a11 = log(Dr.a11);
				Dr.a21 = log(Dr.a21);
				Dr.a31 = log(Dr.a31);
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "lnDr" << Dr << "\n";
#endif
				t3x3 dia = diag(Dr);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "diag lnDr" << dia << "\n";
#endif

				t3x3 Vtrsp = transpose(Vr);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Vtrsp " << Vtrsp << "\n";
#endif

				//eps = (np.dot(V,np.dot(np.diag(d),V.T)).real).reshape(9)
				t3x3 tmp = dot(dia, Vtrsp);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "dot(dia,Vtrsp) " << tmp << "\n";
#endif

				//report strain tensor
				eps = dot(Vr,tmp);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "eps " << eps << "\n";
#endif
				return true;
			}
			else {
				#pragma omp critical
				{
					cerr << "Inversion of stretchU or V failed" << "\n";
				}
				eps = zerot3x3();
				return false;
			}
		}
		else {
			#pragma omp critical
			{
				cerr << "MKL addStrainTensor R inversion failed" << "\n";
			}
			eps = zerot3x3();
			return false;
		}
	}
	else {
		#pragma omp critical
		{
			cerr << "MKL addStrainTensor singular value decomposition of input failed" << "\n";
		}
		eps = zerot3x3();
		return false;
	}
#endif
}


bool computeStrainTensor2( t3x3 const & in, const unsigned int straintensor,
				const unsigned int strainkind, t3x3 & eps)
{
#ifndef INTELMKL_EXISTENT
	cout << "IntelMKL commands were not compiled into the source code!" << "\n";
	eps = failt3x3();
	return false;
#else
	//start with deformation gradient tensor F as input in and make polar decomposition
	t3x3 U = t3x3();
	t3x1 S = t3x1();
	t3x3 Vh = t3x3();

	if ( svd(in, U, S, Vh) == true ) {

		//get rotation of polar decomposition
		t3x3 R = dot(U, Vh);

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "U" << U << "\n";
cout << "S" << S << "\n";
cout << "Vh" << Vh << "\n";
cout << "R" << R << "\n";
#endif

		//inverse of this rotation
		t3x3 Rinv = t3x3();
		if ( inv(R, Rinv) == true ) {
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Rinv" << Rinv << "\n";
#endif
			//for theStretch in theStretches 'U', 'V'

			if ( straintensor == RIGHTCAUCHYGREEN ) {
				//F = RU so RinvF = RinvRU = IdU
				//help = 'material strains based on right Cauchy--Green deformation, i.e., C and U')
				t3x3 stretchU = dot(Rinv, in);
				killnoise( stretchU, static_cast<real_m33>(EPSILON) );

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "sU" << stretchU << "\n";
#endif

				t3x3 epsU = zerot3x3();
				if ( computeEigDecompStretch2Strain( stretchU, strainkind, epsU ) == true ) {
#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "Eigendecomposition of stretchU was successful!" << "\n";
#endif
					eps = epsU;
					return true;
				}
				else {
					eps = zerot3x3();
					return false;
				}
			}
			else if ( straintensor == LEFTCAUCHYGREEN ) {
				//F = VR so FRinv = VRRinv = VId
				//help = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V')
				t3x3 stretchV = dot(in, Rinv);
				killnoise( stretchV, static_cast<real_m33>(EPSILON) );

#ifdef VERBOSE_ADDSTRAINTENSOR
cout << "sV" << stretchV << "\n";
#endif

				//eigen decomposition of stretch tensor
				t3x3 epsV = zerot3x3();
				if ( computeEigDecompStretch2Strain( stretchV, strainkind, epsV ) == true ) {
					eps = epsV;
					return true;
				}
				else {
					eps = zerot3x3();
					return false;
				}
			}
			else {
				#pragma omp critical
				{
					cerr << "Unknown straintensor type" << "\n";
				}

				eps = zerot3x3();
				return false;
			}
		}
		else {
			#pragma omp critical
			{
				cerr << "MKL addStrainTensor R inversion failed" << "\n";
			}
			eps = zerot3x3();
			return false;
		}
	}
	else {
		#pragma omp critical
		{
			cerr << "MKL addStrainTensor singular value decomposition of input failed" << "\n";
		}
		eps = zerot3x3();
		return false;
	}
#endif
}



bool computeCauchy( t3x3 const & F, t3x3 const & P, t3x3 & cauchy)
{
#ifndef INTELMKL_EXISTENT
	cout << "IntelMKL commands were not compiled into the source code!" << "\n";
	cauchy = failt3x3();
	return false;
#else
	real_m33 a = static_cast<real_m33>(1.0);
	a /= det(F);

	t3x3 Ftrsp = transpose(F);
	t3x3 PFtrsp = dot(P,Ftrsp);
	cauchy = mult(a, PFtrsp);
	return true;
#endif
}


real_m33 computeMises( t3x3 const & tensor, bool isstresstensor )
{
#ifndef INTELMKL_EXISTENT
	cout << "IntelMKL commands were not compiled into the source code!" << "\n";
	#ifdef SINGLE_PRECISION
		return numeric_limits<float>::max();
	#else
		return numeric_limits<double>::max();
	#endif
#else
	//get trace
	real_m33 trace_onethird = static_cast<real_m33>(1.0/3.0) * trace(tensor);
	t3x3 tmp = mult(trace_onethird, eye() );
	//deviator
	t3x3 dev = substrct(tensor,tmp);
	t3x3 devtrsp = transpose(dev);

	t3x3 symdev = mult(static_cast<real_m33>(0.5), add(dev,devtrsp));

//cout << symdev << "\n";
	t3x3 symdevtrsp = transpose(symdev);

//cout << symdevtrsp << "\n";
	real_m33 s = sum(dyadic(symdev,symdevtrsp));

//cout << s << "\n";

	if (isstresstensor == true)
		return sqrt(static_cast<real_m33>(3.0/2.0) * s);
	else //assume strain tensor
		return sqrt(static_cast<real_m33>(2.0/3.0) * s);
#endif
}


aabb3d get_corners( vector<real_xyz> &cx, vector<real_xyz> &cy, vector<real_xyz> &cz )
{
	aabb3d out = aabb3d();
	size_t i = 0;

	size_t ni = cx.size();
	for ( i = 0; i < ni; ++i ) {
		if ( cx[i] <= out.xmi ) { out.xmi = cx[i]; }
		if ( cx[i] >= out.xmx ) { out.xmx = cx[i]; }
	}
	ni = cy.size();
	for ( i = 0; i < ni; ++i ) {
		if ( cy[i] <= out.ymi ) { out.ymi = cy[i]; }
		if ( cy[i] >= out.ymx ) { out.ymx = cy[i]; }
	}
	ni = cz.size();
	for ( i = 0; i < ni; ++i ) {
		if ( cz[i] <= out.zmi ) { out.zmi = cz[i]; }
		if ( cz[i] >= out.zmx ) { out.zmx = cz[i]; }
	}

	cout << "Corners found" << "\n" << out << "\n";
	return out;
}



rfftn::rfftn()
{
	fftTempP = NULL;
}

rfftn::~rfftn()
{
	//do not delete fftTempP, only pointer to vector!
	DftiFreeDescriptor(&m_f_handle);
}


void rfftn::init(unsigned int const * ngridd)
{
#ifdef SINGLE_PRECISION
	cout << "MKL FFT rfftn single precision not yet implemented" << "\n";
	return;
#else

	precision = DFTI_DOUBLE;

	//prepare forward fft input buffer with real-valued deformation gradient data
	NX = static_cast<MKL_LONG>(ngridd[0]); //nfft
	NY = static_cast<MKL_LONG>(ngridd[1]); //nfft
	NZ = static_cast<MKL_LONG>(ngridd[2]); //##MK::check carefully for cuboidal geometries
	NXYZ = NX*NY*NZ;
	dimensions[0] = NX;
	dimensions[1] = NY;
	dimensions[2] = NZ;
	input_strides[0] = 0;
	input_strides[1] = NX;
	input_strides[2] = NX*NY;
	input_strides[3] = 1;

	size_t requiredSize = 	static_cast<size_t>(NX) *
							static_cast<size_t>(NY) *
							static_cast<size_t>(NZ);

	m_input.resize(requiredSize);

	//prepare forward FFT output buffer of MKL_Complex16
	NI = NX;
	NJ = NY;
	NK = static_cast<MKL_LONG>(floor(NZ/2)+1); //##MK::check carefully for cuboidal grids!
	NJK = NJ * NK;
	NIJK = NI * NJ * NK;

	output_strides[0] = 0;
	output_strides[1] = NK; //static_cast<MKL_LONG>(floor(NFFT/2))+1);
	output_strides[2] = NK*NJ; //(static_cast<MKL_LONG>(floor(NFFT/2))+1)*NFFT;
	output_strides[3] = 1;

	requiredSize = 	static_cast<size_t>(NI) *
					static_cast<size_t>(NJ) *
					static_cast<size_t>(NK) * sizeof(MKL_Complex16); //##MK::sizeof necessary as carrier buffer is vector<char>

	fftTemp.resize(requiredSize); //on empty vector resize set new elements, i.e. all to zero

	//redefine how to interpret this carrier buffer
	fftTempP = (MKL_Complex16*) &fftTemp[0];

cout << "Initialized fFFT with NXYZ = " << NX << ";" << NY << ";" << NZ << " NIJK = " << NI << ";" << NJ << ";" << NK << "\n";
#endif
}


void rfftn::fill(vector<t3x3> const & F, const unsigned int c)
{
	for (size_t eipid = 0; eipid < static_cast<size_t>(NXYZ); ++eipid) {
		if ( c == 0 ) { m_input[eipid] = static_cast<double>(F[eipid].a11); continue; }
		if ( c == 1 ) { m_input[eipid] = static_cast<double>(F[eipid].a12); continue; }
		if ( c == 2 ) { m_input[eipid] = static_cast<double>(F[eipid].a13); continue; }
		if ( c == 3 ) { m_input[eipid] = static_cast<double>(F[eipid].a21); continue; }
		if ( c == 4 ) { m_input[eipid] = static_cast<double>(F[eipid].a22); continue; }
		if ( c == 5 ) { m_input[eipid] = static_cast<double>(F[eipid].a23); continue; }
		if ( c == 6 ) { m_input[eipid] = static_cast<double>(F[eipid].a31); continue; }
		if ( c == 7 ) { m_input[eipid] = static_cast<double>(F[eipid].a32); continue; }
		if ( c == 8 ) { m_input[eipid] = static_cast<double>(F[eipid].a33); continue; }
	}
}


void rfftn::forwardFFT()
{
	DftiCreateDescriptor(&m_f_handle, precision, DFTI_REAL, 3, dimensions);
	DftiSetValue(m_f_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiSetValue(m_f_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
	DftiSetValue(m_f_handle, DFTI_INPUT_STRIDES, input_strides);
	DftiSetValue(m_f_handle, DFTI_OUTPUT_STRIDES, output_strides);

	DftiCommitDescriptor(m_f_handle);

	DftiComputeForward(m_f_handle, &m_input[0] , fftTempP);

	//##MK::DEBUG
	/* 	cout << "Done, printing results" << "\n";
	int jump = 0;
	for (int j = 0; j < elements; ++j) {
		cout << "(" << fftTempP[j].real << " , " << fftTempP[j].imag << ")" << "\n";
		jump++;
		if ( jump != 8)
			continue;
		else
			cout << "\n"; jump = 0;
	}
	*/
}


irfftn::irfftn()
{
	fftTempP = NULL;
}

irfftn::~irfftn()
{
	//do not delete fftTempP, only pointer to vector!
	DftiFreeDescriptor(&m_b_handle);
}


void irfftn::init( unsigned int const * nfft, unsigned int const * ngridd )
{
#ifdef SINGLE_PRECISION
cout << "MKL FFT rfftn single precision not yet implemented" << "\n";
	return;
#else
	precision = DFTI_DOUBLE;

	//prepare inverse FFT MKL_Complex16 input buffer
	NI = static_cast<MKL_LONG>(ngridd[0]);
	NJ = static_cast<MKL_LONG>(ngridd[1]);
	NK = static_cast<MKL_LONG>(ngridd[2]);
	NJK = NJ*NK;
	NIJK = NI*NJ*NK;
	input_strides[0] = 0;
	input_strides[1] = NK;
	input_strides[2] = NK*NJ;
	input_strides[3] = 1;
	size_t requiredSize = 	static_cast<size_t>(NI) *
							static_cast<size_t>(NJ) *
							static_cast<size_t>(NK) * sizeof(MKL_Complex16);

	fftTemp.resize(requiredSize);
	fftTempP = (MKL_Complex16*) &fftTemp[0];

	//prepare inverse FFT results buffer, one real value per 3d real space coordinate
	//NFFT = static_cast<MKL_LONG>(nfft);
	NX = static_cast<MKL_LONG>(nfft[0]);
	NY = static_cast<MKL_LONG>(nfft[1]);
	NZ = static_cast<MKL_LONG>(nfft[2]);
	NXY = NX * NY;
	NXYZ = NX * NY * NZ;
	dimensions[0] = NX; //NFFT;
	dimensions[1] = NY; //NFFT;
	dimensions[2] = NZ; //NFFT;
	output_strides[0] = 0;
	output_strides[1] = NX; //NFFT;
	output_strides[2] = NX*NY; //NFFT*NFFT;
	output_strides[3] = 1;

	requiredSize =	static_cast<size_t>(NX) *
					static_cast<size_t>(NY) *
					static_cast<size_t>(NZ); //MK::no * sizeof(precision) because is not a character but a double vector!

	m_output.resize(requiredSize);

cout << "Initialized iFFT with NIJK = " << NI << ";" << NJ << ";" << NK << " NXYZ = " << NX << ";" << NY << ";" << NZ << "\n";
#endif
}



void irfftn::fill(vector<double> const & in_real, vector<double> const & in_imag, const unsigned int cm)
{
	MKL_LONG offset = 0;
	for (MKL_LONG i = 0; i < NI; ++i) {
		for (MKL_LONG j = 0; j < NJ; ++j) {
			for (MKL_LONG k = 0; k < NK; ++k) {
				offset = i*NJK+j*NK+k;
				fftTempP[offset].real = in_real.at((offset*3)+cm);
				fftTempP[offset].imag = in_imag.at((offset*3)+cm);
			}
		}
	}
}


void irfftn::backwardFFTandNormalize()
{
	DftiCreateDescriptor(&m_b_handle, precision, DFTI_REAL, 3, dimensions);
	DftiSetValue(m_b_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiSetValue(m_b_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
	DftiSetValue(m_b_handle, DFTI_INPUT_STRIDES, input_strides);
	DftiSetValue(m_b_handle, DFTI_OUTPUT_STRIDES, output_strides);
	DftiCommitDescriptor(m_b_handle);

	DftiComputeBackward(m_b_handle, fftTempP, &m_output[0]);

	/*
	//##MK::DEBUG
cout << "Done, printing results size of m_output " << m_output.size() << "\n";
		for (size_t e = 0; e < m_output.size(); ++e) {
		cout << m_output.at(e) << "\n";
	}
	*/

	//normalization
	//double _norm = (1.0/static_cast<double>(NFFT)) * (1.0/static_cast<double>(NFFT)) * (1.0/static_cast<double>(NFFT)); //##MK::check carefully for cuboidal geometry

	double _norm = (1.0/static_cast<double>(NX)) * (1.0/static_cast<double>(NY)) * (1.0/static_cast<double>(NZ)); //##MK::check carefully for cuboidal geometry
cout << "Performing normalization with " << _norm << "\n";
	for (size_t e = 0; e < m_output.size(); ++e) {
		m_output[e] *= _norm;
	}
}


unsigned int fourier_transform_defgradient( vector<t3x3> const & defgrad, unsigned int const * ngr, double const * nsz, vector<rfftn*> & ffts )
{
	//3d discrete fft for real deformation gradient tensor components

	unsigned int component = 0;
	for (             ; component < 9; ++component) {
		rfftn* anfft = NULL;
		try { anfft = new rfftn; } //one component 3d fft per class object
		catch (bad_alloc &croak) {
			for (size_t i = 0; i < ffts.size(); ++i ) {
				if ( ffts.at(i) != NULL ) {
					delete ffts.at(i); ffts.at(i) = NULL;
				}
			}
			#pragma omp critical
			{
				cerr << "Fourier transformation of F components failed because FFT class objects were not allocateable" << "\n";
			}
			return 0;
		}

		anfft->init( ngr );
		anfft->fill(defgrad, component);
		anfft->forwardFFT();

		ffts.push_back(anfft);

cout << "FFT for component " << component << " done!" << "\n";
	}


	/*
	//##MK::DEBUG
	//organized output to check - real
cout << ffts.at(0)->NIJK <<  " component values are in total, their real parts are " << "\n";
	for ( unsigned int ffti = 0; ffti < ffts.at(0)->NIJK; ++ffti ) { //##MK::all ffts on same grid
		for (unsigned int c = 0; c < 8; ++c) {
			cout << ffts.at(c)->fftTempP[ffti].real << ";";
		}
		cout << ffts.at(8)->fftTempP[ffti].real << "\n";
	}

	cout << "\n" << "\n";
	//imag
cout << ffts.at(0)->NIJK <<  " component values are in total, their imaginary parts are " << "\n";
	for ( unsigned int ffti = 0; ffti < ffts.at(0)->NIJK; ++ffti ) {
		for (unsigned int c = 0; c < 8; ++c) {
			cout << ffts.at(c)->fftTempP[ffti].imag << ";";
		}
		cout << ffts.at(8)->fftTempP[ffti].imag << "\n";
	}
	*/

	return component;
}


unsigned int displacement_fluctuations( vector<rfftn*> const & ffts, unsigned int const * ngr,
		double const * nsz, vector<irfftn*> & iffts, vector<d3d> & flucts )
{
	//compute local position fluctuation from deformation gradient tensor ffts
	//##MK::dangerous! rather parse form fft directly or make displacement_fluctuation class with tools fft
	cout << "Now computing fluctuant displacements..." << "\n";

	//MK::do not use unsigned integer types for grd as these are passed to k_s inplace
	//they can have negative values (not the grid dimension itself but the frequency supporting points)
	int grd[3] = { 	static_cast<int>(ffts.at(0)->NI),
					static_cast<int>(ffts.at(0)->NJ),
					static_cast<int>(ffts.at(0)->NK) }; //##MK::NI old NZ, NJ old NY, NI old NK

	unsigned int ugrd[3] = { 	static_cast<unsigned int>(grd[0]),
								static_cast<unsigned int>(grd[1]),
								static_cast<unsigned int>(grd[2]) };

	int thresholds[3] = { 	static_cast<int>(floor(static_cast<double>(grd[0])/2.0)),
							static_cast<int>(floor(static_cast<double>(grd[1])/2.0)),
							static_cast<int>(floor(static_cast<double>(grd[2])/2.0)) };

	//##MK::works for cubic meshes only!
	//##MK::test for arbitrary cuboidal geometries...

	//build meshgrid of k_s inplace
	int val = 0;
	vector<int> k_s;
	k_s.reserve(grd[0]*grd[1]*grd[2]*3);
	for (int i = 0; i < grd[0]; ++i ) {
		for (int j = 0; j < grd[1]; ++j ) {
			for ( int k = 0; k < grd[2]; ++k ) {
				val = k;
				k_s.push_back(val);
				val = (j > thresholds[1]) ? j-grd[1] : j;
				k_s.push_back(val);
				val = (i > thresholds[0]) ? i-grd[0] : i;
				k_s.push_back(val);
			}
		}
	}

	/*
	cout << "k_s" << "\n";
	for ( size_t e = 0; e < k_s.size(); e=e+3 ) {
		cout << k_s[e] << "\n";
		cout << k_s[e+1] << "\n";
		cout << k_s[e+2] << "\n";
		cout << "\n" << "\n";
	}
	*/

	vector<double> k_sSquared;
	k_sSquared.reserve(ugrd[0]*ugrd[1]*ugrd[2]);

	int NI = grd[0];
	int NJ = grd[1];
	int NK = grd[2];
	int NJK = NJ * NK; //##MK::dangerous! rather parse form fft directly or make displacement_fluctuation class with tools fft
	int NIJK = NI * NJ * NK;

	long total = 0;
	int where = 0; //##MK::mind grid resolution!
	//MK::quick and dirty overflow identification only meant to indicate code revision necessary
	//however, currently DAMASK 3d simulations with 1024^3 grid points are completely out of scope!

	size_t intoverflow = static_cast<size_t>(NI) * static_cast<size_t>(NJ) * static_cast<size_t>(NK);
	if ( intoverflow >= static_cast<size_t>(numeric_limits<int>::max()) ) {
		cerr << "Potential INT overflow please revise code to change ints in this einsum to long types" << "\n";
		return 0;
	}

	//np.einsum('...l,...l',k_s,k_s)
	for (int i = 0; i < NI; ++i ) {
		for (int j = 0; j < NJ; ++j ) {
			for ( int k = 0; k < NK; ++k ) {
				//MK::innermost loop --- for ( int l = 0; l < 3; ++l) --- unrolled
				where = (i*NJK+j*NK+k)*3;
				total = static_cast<long>( SQR(k_s.at(where+0)) + SQR(k_s.at(where+1)) + SQR(k_s.at(where+2)) ); //##MK::replace at some point by [] for performance
				k_sSquared.push_back(static_cast<double>(total));
			}
		}
	}

/* 	for ( size_t e = 0; e < k_sSquared.size(); ++e)
		cout << k_sSquared.at(e) << "\n"; */

	//ignore global average frequency
	k_sSquared.at(0) = 1.0;


	//MK::integration in Fourier space involving complex numbers!
	//double integrator_real[3] = {0.0, 0.0, 0.0};
	double integrator_imag[3] = { 0.5*nsz[0]/MYPI, 0.5*nsz[1]/MYPI, 0.5*nsz[2]/MYPI };

	//-np.einsum('ijkml,ijkl,l->ijkm', F, k_s, integrator) / k_sSquared(ijk)
	vector<double> displacement_fourier_real;
	displacement_fourier_real.reserve(NIJK*3);
	vector<double> displacement_fourier_imag;
	displacement_fourier_imag.reserve(NIJK*3);

	double tmpreal = 0.0;
	double tmpimag = 0.0;
	for (int i = 0; i < NI; ++i ) {
		for (int j = 0; j < NJ; ++j ) {
			for ( int k = 0; k < NK; ++k ) {
				for ( int m = 0; m < 3; ++m ) { //##MK::potentially unroll
					tmpreal = 0.0;
					tmpimag = 0.0;

					size_t where = 	(static_cast<size_t>(i) * static_cast<size_t>(NJK)) +
									(static_cast<size_t>(j) * static_cast<size_t>(NK)) +
									static_cast<size_t>(k);

					for (int l = 0; l < 3; ++l ) {
						/*
						//fourier transformed component values of F are:
						//1 --> m=0,l=0
						//2 --> m=0,l=1 //row-wise layout
						//3 --> m=0,l=2
						//4 --> m=1,l=0
						//5 --> m=1,l=1
						//6 --> m=1,l=2
						//7 --> m=2,l=0
						//8 --> m=2,l=1
						//9 --> m=2,l=2 --> m*3+l --> e.g. 0-th entry value are kept in ffts.at(m*3+l)->fftTempP[i*NJK+j*NK+k].real
						*/
						//MK::mind that we multiply three complex numbers Fourier(i,j,k,m,l)*k_s(i,j,k,m)!
						/*
						double ar = ffts.at(m*3+l)->fftTempP[i*NJK+j*NK+k].real;
						double ai = ffts.at(m*3+l)->fftTempP[i*NJK+j*NK+k].imag;
						double br = static_cast<double>(k_s.at((i*NJK+j*NK+k)*3+l));
						double bi = 0.0;
						double cr = integrator_real[l];
						double ci = integrator_imag[l];

						double abr = ar*br; // - 0.0; //ar*br - ai*bi;
						double abi = br*ai; //0.0 + br*ai; //ar*bi + br*ai;

						//remember that double integrator_real[3] = {0.0, 0.0, 0.0} i.e. cr = 0;

						double abcr = -1.0*abi*ci; //0.0 - abi*ci; //abr*cr - abi*ci;
						double abci = abr*ci; //abr*ci + 0.0*abi; //abr*ci + cr*abi;
						*/

						tmpreal = tmpreal + ( -1.0*ffts.at(m*3+l)->fftTempP[where].imag * static_cast<double>(k_s.at((where*3)+l)) * integrator_imag[l] ); //tmpreal += abcr;
						tmpimag = tmpimag + ( ffts.at(m*3+l)->fftTempP[where].real * static_cast<double>(k_s.at((where*3)+l)) * integrator_imag[l] ); //tmpimag += abci;

						//MK::effect of k_sSquare on real+imag result ijkm from einsum is just a division of both the real and imag part

					} //l done
					/*
					double dr = tmpreal;
					double di = tmpimag;

					double er = -1.0;
					double ei = 0.0;
					*/
					//MK:: -np.einsum
					double der = -1.0*tmpreal; //tmpreal*-1.0 - 0.0; //dr*er - di*ei;
					double dei = -1.0*tmpimag; //tmpreal*0.0 + -1.0*tmpimag; //dr*ei + er*di;

					//MK:: / k_sSquared
					double fr = k_sSquared.at(where);
					//double fi = 0.0;

					double resr = der * fr / SQR(fr); //(der*fr + dei*fi) / (SQR(fr)+SQR(fi));
					double resi = fr * dei / SQR(fr); //(fr*dei - der*fi) / (SQR(fr)+SQR(fi));

					displacement_fourier_real.push_back(resr);
					displacement_fourier_imag.push_back(resi);
					//##MK::SIMD?
				} //m done
			} //k
		} //j
	} //i

	/*
	//##MK::DEBUG
	cout << "displreal.size() " << displacement_fourier_real.size() << "\n";
	for (size_t e = 0; e < displacement_fourier_real.size(); ++e) {
		cout << displacement_fourier_real.at(e) << "\n";
	}
	cout << "\n";
	cout << "displimag.size() " << displacement_fourier_imag.size() << "\n";
	for (size_t e = 0; e < displacement_fourier_imag.size(); ++e) {
		cout << displacement_fourier_imag.at(e) << "\n";
	}
	*/

	//backtransformation for each coordinate, via iFFT for each i,j,k position in frequency space there is one real result of a fluctDisplacement
	unsigned int component = 0;
	for (        ; component < 3; ++component) { //three coordinate directions
		irfftn* anifft = NULL;
		try { anifft = new irfftn; } //one component 3d fft per class object
		catch (bad_alloc &croak) {
			for (size_t i = 0; i < iffts.size(); ++i ) {
				if ( iffts.at(i) != NULL ) {
					delete iffts.at(i); iffts.at(i) = NULL;
				}
			}
			cerr << "Fluctuation displacement computation was unsuccessful because of failing allocation during iFFT" << "\n";
			return 0;
		}
		anifft->init( ngr, ugrd ); //##MK::unneccessary nsz );
		anifft->fill(displacement_fourier_real, displacement_fourier_imag, component);
		anifft->backwardFFTandNormalize();

		iffts.push_back(anifft);

cout << "iFFT for component " << component << " done!" << "\n";
	}

	//transfer displacements
	size_t neipid = iffts.at(0)->m_output.size();
	for (size_t eipid = 0; eipid < neipid; ++eipid) {
		flucts.push_back( d3d( 	static_cast<real_xyz>(iffts.at(0)->m_output.at(eipid)),
								static_cast<real_xyz>(iffts.at(1)->m_output.at(eipid)),
								static_cast<real_xyz>(iffts.at(2)->m_output.at(eipid)) )
						);
		//##MK::why is this NX*NY*NZ the same as Martin/Pratheeks?
		//##MK::potentially caused by numpy's "Normalization mode such that since numpy v.1.10 results are scaled by 1/n by default,
		//such the passing of grid[::-1] within addDisplacement.py normalizes as such
		//Intel MKL iFFT does apparently not apply such normalization hence the consistent factor of 1000.0 within machine precision
		//deviation to martins result
	}

	/*
	//##MK::DEBUG
	cout << "Computed, the fluctuant displacements per integration point are" << "\n";
	for (size_t eipid = 0; eipid < neipid; ++eipid) {
		cout << flucts.at(eipid).x << ";" << flucts.at(eipid).y << ";" << flucts.at(eipid).z << "\n";
	}
	*/

	return component;
}


void displacement_average( vector<rfftn*> const & ffts,
		unsigned int const * ngr, double const * nsz, vector<d3d> & flucts )
{
	//compute average position fluctuation from deformation gradient tensor ffts

	cout << "Now computing average displacement..." << "\n";

	//build coordinate mesh linspace(0.0,nsz[i],ngr[i],endpoint=false)
	size_t requiredSize = 	static_cast<size_t>(ngr[0]) *
							static_cast<size_t>(ngr[1]) *
							static_cast<size_t>(ngr[2]) * 3; //ijkm

	vector<real_xyz> xyz;
	xyz.resize(requiredSize);

	real_xyz scaler[3] = { 	static_cast<real_xyz>(nsz[0]) / static_cast<real_xyz>(ngr[0]),
							static_cast<real_xyz>(nsz[1]) / static_cast<real_xyz>(ngr[1]),
							static_cast<real_xyz>(nsz[2]) / static_cast<real_xyz>(ngr[2]) };

	//##MK::in place computation to save memory consider to do inplace with numpy.einsum
	/*unsigned int where = 0;
	for ( unsigned int z = 0; z < ngr[2]; ++z) {
		real_xyz zz = static_cast<real_xyz>(z) * scaler[2];
		for (unsigned int y = 0; y < ngr[1]; ++y) {
			real_xyz yy = static_cast<real_xyz>(y) * scaler[1];
			for (unsigned int x  = 0; x < ngr[0]; ++x) {
				real_xyz xx = static_cast<real_xyz>(x) * scaler[0];
				xyz.at(where+0) = xx;
				xyz.at(where+1) = yy;
				xyz.at(where+2) = zz;
				where+=3;
				cout << xx << ";" << yy << ";" << zz << "\n";
			}
		}
	}*/

	//(Favg = np.real(F_fourier[0,0,0,:,:])/ngr.prod())-eye()
	real_xyz _gridprod = 	(1.0/static_cast<real_xyz>(ngr[0])) *
							(1.0/static_cast<real_xyz>(ngr[1])) *
							(1.0/static_cast<real_xyz>(ngr[2]));

	//MK::fourier transformed component values of F are:
	//1 --> m=0,l=0
	//2 --> m=0,l=1 //row-wise layout
	//3 --> m=0,l=2
	//4 --> m=1,l=0
	//5 --> m=1,l=1
	//6 --> m=1,l=2
	//7 --> m=2,l=0
	//8 --> m=2,l=1
	//9 --> m=2,l=2 --> m*3+l --> e.g. 0-th entry value are kept in ffts.at(m*3+l)->fftTempP[i*NJK+j*NK+k].real
	real_xyz FavgEye[9] = {0.0};
	//F[0,0,0,:]
	FavgEye[0] = (static_cast<real_xyz>(ffts.at(0)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(1.0); //ffts.at(0*3+0)->fftTempP[0*NJK+0*NK+0].real;
	FavgEye[1] = (static_cast<real_xyz>(ffts.at(1)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[2] = (static_cast<real_xyz>(ffts.at(2)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[3] = (static_cast<real_xyz>(ffts.at(3)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[4] = (static_cast<real_xyz>(ffts.at(4)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(1.0);
	FavgEye[5] = (static_cast<real_xyz>(ffts.at(5)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[6] = (static_cast<real_xyz>(ffts.at(6)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[7] = (static_cast<real_xyz>(ffts.at(7)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(0.0);
	FavgEye[8] = (static_cast<real_xyz>(ffts.at(8)->fftTempP[0].real) * _gridprod) - static_cast<real_xyz>(1.0);

	//np.einsum('ml,ijkl->ijkm',Favg-np.eye(3),origCoords)
	//real_xyz total = 0.0;
	for ( unsigned int z = 0; z < ngr[2]; ++z) {
		real_xyz zz = static_cast<real_xyz>(z) * scaler[2];
		for (unsigned int y = 0; y < ngr[1]; ++y) {
			real_xyz yy = static_cast<real_xyz>(y) * scaler[1];
			for (unsigned int x  = 0; x < ngr[0]; ++x) {
				real_xyz xx = static_cast<real_xyz>(x) * scaler[0];

				d3d thisone = d3d();
				//for (unsigned int m = 0; m < 3; ++m) { //collect x,y,z
				/*total = 0.0; //m=0
				for (unsigned int l = 0; l < 3; ++l ) {
					total += FavgEye[0*3+l] * xyz.at((z*ngr[1]*ngr[0]+y*ngr[0]+x)*3+l);
				}
				thisone.x = total;*/ //MK::lessons how to turn from memory-bound to compute-bound...!
				thisone.dx = (FavgEye[0*3+0]*xx) + (FavgEye[0*3+1]*yy) + (FavgEye[0*3+2]*zz);

				/*total = 0.0; //m=1
				for (unsigned int l = 0; l < 3; ++l ) {
					total += FavgEye[1*3+l] * xyz.at((z*ngr[1]*ngr[0]+y*ngr[0]+x)*3+l);
				}
				thisone.y = total;*/
				thisone.dy = (FavgEye[1*3+0]*xx) + (FavgEye[1*3+1]*yy) + (FavgEye[1*3+2]*zz);

				/*total = 0.0; //m=2
				for (unsigned int l = 0; l < 3; ++l ) {
					total += FavgEye[2*3+l] * xyz.at((z*ngr[1]*ngr[0]+y*ngr[0]+x)*3+l);
				}
				thisone.z = total;*/
				thisone.dz = (FavgEye[2*3+0]*xx) + (FavgEye[2*3+1]*yy) + (FavgEye[2*3+2]*zz);

				flucts.push_back(thisone);
				//} //m
			} //x
		} //y
	} //z

	/*
	//##MK::DEBUG
	cout << "Computed, the average displacements per integration point are" << "\n";
	size_t neipid = flucts.size();
	for (size_t eipid = 0; eipid < neipid; ++eipid) {
		cout << flucts.at(eipid).dx << ";" << flucts.at(eipid).dy << ";" << flucts.at(eipid).dz << "\n";
	}
	*/
}
