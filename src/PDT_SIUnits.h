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

#ifndef __PDT_SIUNITS_H__
#define __PDT_SIUNITS_H__

//SI unit conversions
#define CELSIUS2KELVIN(T)				((273.15) + (T))
#define KELVIN2CELSIUS(T)				((T) - (273.15))

#define DEGREE2RADIANT(theta)			((theta)/(180.0)*(PI))
#define RADIANT2DEGREE(rad)				((rad)/(PI)*(180.0))

//scaling conversions
#define METER2CENTIMETER(x)				((x)*(1.0e2))
#define CENTIMETER2METER(x)				((x)*(1.0e-2))

#define METER2MILLIMETER(x)				((x)*(1.0e3))
#define MILLIMETER2METER(x)				((x)*(1.0e-3))

#define METER2MICRON(x)					((x)*(1.0e6))
#define MICRON2METER(x)					((x)*(1.0e-6))

#define METER2NANOMETER(x)				((x)*(1.0e9))
#define NANO2METER(x)					((x)*(1.0e-9))

#define METER2PICOMETER(x)				((x)*(1.0e12))
#define PICO2METER(x)					((x)*(1.0e-12))


#endif
