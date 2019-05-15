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


#include "PDT_HDF5.h"


void debug_hdf5( void )
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* Copyright by The HDF Group.                                               *
	* Copyright by the Board of Trustees of the University of Illinois.         *
	* All rights reserved.                                                      *
	*                                                                           *
	* This file is part of HDF5.  The full HDF5 copyright notice, including     *
	* terms governing use, modification, and redistribution, is contained in    *
	* the COPYING file, which can be found at the root of the source code       *
	* distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
	* If you do not have access to either file, you may request a copy from     *
	* help@hdfgroup.org.                                                        *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	#define FILE "dset.h5"
	//identifiers
	hid_t       file_id, dataset_id, dataspace_id;
	hsize_t     dims[2];
	herr_t      status;

	//create a new file using default properties
	file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create the data space for the dataset
	dims[0] = 4;
	dims[1] = 6;
	dataspace_id = H5Screate_simple(2, dims, NULL);

	//create the dataset
	dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id,
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//end access to the dataset and release resources used by it
	status = H5Dclose(dataset_id);

	//terminate access to the data space
	status = H5Sclose(dataspace_id);

	//close the file
	status = H5Fclose(file_id);
}
