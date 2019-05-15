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

#ifndef __PDT_ABINITIO_H__
#define __PDT_ABINITIO_H__
//natural constants
#define MYPI							(3.141592653589793238462643383279502884197169399375105820974)
#define kboltzmann						(1.3806488e-23)		//J/Kmol
#define echarge							(1.602176565e-19)	//Coulomb
#define Navogadro						(6.022140857e+23)	//1/mol
#define RGAS							(8.31446154) 		//(Navogadro)*(kboltzmann)) //J/K/mol

//natural beauty
#define SQR(a)							((a)*(a))
#define CUBE(a)							((a)*(a)*(a))
#define MIN(X,Y)						(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)						(((X) < (Y)) ? (Y) : (X))
#define CLAMP(x,lo,hi)					(MIN(MAX((x), (lo)), (hi)))

/*
//logical
#define NO								0x00
#define YES								0x01
*/

//naming definitions for more intuitive code reading
#define TWO_DIMENSIONS					2
#define THREE_DIMENSIONS				3

#endif
