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


#ifndef __PDT_MATH_H__
#define __PDT_MATH_H__

#include "PDT_OriMath.h"


//http://rosettacode.org/wiki/Reduced_row_echelon_form#C.2B.2B
// Matrix traits: This describes how a matrix is accessed. By
// externalizing this information into a traits class, the same code
// can be used both with native arrays and matrix classes. To use the
// default implementation of the traits class, a matrix type has to
// provide the following definitions as members:
//
// * typedef ... index_type;
//   - The type used for indexing (e.g. size_t)
// * typedef ... value_type;
//   - The element type of the matrix (e.g. double)
// * index_type min_row() const;
//   - returns the minimal allowed row index
// * index_type max_row() const;
//   - returns the maximal allowed row index
// * index_type min_column() const;
//   - returns the minimal allowed column index
// * index_type max_column() const;
//   - returns the maximal allowed column index
// * value_type& operator()(index_type i, index_type k)
//   - returns a reference to the element i,k, where
//     min_row() <= i <= max_row()
//     min_column() <= k <= max_column()
// * value_type operator()(index_type i, index_type k) const
//   - returns the value of element i,k
//
// Note that the functions are all inline and simple, so the compiler
// should completely optimize them away.
template<typename MatrixType> struct matrix_traits
{
  typedef typename MatrixType::index_type index_type;
  typedef typename MatrixType::value_type value_type;
  static index_type min_row(MatrixType const& A)
  { return A.min_row(); }
  static index_type max_row(MatrixType const& A)
  { return A.max_row(); }
  static index_type min_column(MatrixType const& A)
  { return A.min_column(); }
  static index_type max_column(MatrixType const& A)
  { return A.max_column(); }
  static value_type& element(MatrixType& A, index_type i, index_type k)
  { return A(i,k); }
  static value_type element(MatrixType const& A, index_type i, index_type k)
  { return A(i,k); }
};

// specialization of the matrix traits for built-in two-dimensional
// arrays
template<typename T, std::size_t rows, std::size_t columns>
 struct matrix_traits<T[rows][columns]>
{
  typedef std::size_t index_type;
  typedef T value_type;
  static index_type min_row(T const (&)[rows][columns])
  { return 0; }
  static index_type max_row(T const (&)[rows][columns])
  { return rows-1; }
  static index_type min_column(T const (&)[rows][columns])
  { return 0; }
  static index_type max_column(T const (&)[rows][columns])
  { return columns-1; }
  static value_type& element(T (&A)[rows][columns],
                             index_type i, index_type k)
  { return A[i][k]; }
  static value_type element(T const (&A)[rows][columns],
                            index_type i, index_type k)
  { return A[i][k]; }
};

// Swap rows i and k of a matrix A
// Note that due to the reference, both dimensions are preserved for
// built-in arrays
template<typename MatrixType>
 void swap_rows(MatrixType& A,
                 typename matrix_traits<MatrixType>::index_type i,
                 typename matrix_traits<MatrixType>::index_type k)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  // check indices
  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    std::swap(mt.element(A, i, col), mt.element(A, k, col));
}

// divide row i of matrix A by v
template<typename MatrixType>
 void divide_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::value_type v)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(v != 0);

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) /= v;
}

// in matrix A, add v times row k to row i
template<typename MatrixType>
 void add_multiple_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::index_type k,
                  typename matrix_traits<MatrixType>::value_type v)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) += v * mt.element(A, k, col);
}

// convert A to reduced row echelon form
template<typename MatrixType>
 void to_reduced_row_echelon_form(MatrixType& A)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  index_type lead = mt.min_row(A);

  for (index_type row = mt.min_row(A); row <= mt.max_row(A); ++row)
  {
    if (lead > mt.max_column(A))
      return;
    index_type i = row;
    while (mt.element(A, i, lead) == 0)
    {
      ++i;
      if (i > mt.max_row(A))
      {
        i = row;
        ++lead;
        if (lead > mt.max_column(A))
          return;
      }
    }
    swap_rows(A, i, row);
    divide_row(A, row, mt.element(A, row, lead));
    for (i = mt.min_row(A); i <= mt.max_row(A); ++i)
    {
      if (i != row)
        add_multiple_row(A, i, row, -mt.element(A, i, lead));
    }
  }
}



#endif
