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


#include "PDT_GrainPerImgRecon.h"

perImgReconHdl::perImgReconHdl()
{
//	cgid = numeric_limits<unsigned int>::max();
//	nuip = numeric_limits<unsigned int>::max();
//	//ips and localbox construct themselves
}

perImgReconHdl::~perImgReconHdl()
{
//	ips.clear();
}

void perImgReconHdl::compute_all_images()
{
/*	//forall elements and all their periodic images, check if in owin_final, if so bin them
	if ( healthy == false) {
		return;
	}

	//thread-parallelized computing of included deformed configuration element positions and their periodic images surplus their target bin
	#pragma omp parallel
	{
		//##MK::one might optimize further this function for generating less duplicate perImages
		//as the threads which take care of ip slabs which are neither the upper nor lower slab generate multiply several points at the boundary?

		//thread-local temporary memory mybuffer to store coordinates metadata id and target bin to allow fully parallel periodic images, checking, and binning
		//with minimal concurrence when writing to global datastructure
		int tid = omp_get_thread_num();
		bool mehealthy = true;

		vector<vector<p3dm3>*> mybuffer;
		vector<p3dm3>* s = NULL;
		try { s = new vector<p3dm3>; }
		catch (bad_alloc &exc) {
			stopping("Unable to allocate binning memory", owner->get_myrank(), tid);
			mehealthy = false;
		}
		mybuffer.push_back(s);

		if ( mehealthy == true ) {
			memRegion* thisregion = owner->db.at(tid);
			size_t eipn = thisregion->eipid_n;
			size_t lc2gl_idoffset = thisregion->eipid_start;
			bv3x3 per = thisregion->grid.perOffset; //tricline deformed rve base column vectors
			aabb3d w = owin_final;

			p3d cache[27]; //cache periodic images values of a single element
			unsigned int targetbin[27]; //cache bins into which the points get stored
			unsigned char periodicImage[27]; //cache bins which of the Moore periodic images it is, 0x00 is the center, 1 is -1-1-1, 2 0-1-1, and so forth

			unsigned int i = 0;
			real_xyz xxx = 0.0;
			real_xyz yyy = 0.0;
			real_xyz zzz = 0.0;
			real_xyz h = 0.0;
			real_xyz k = 0.0;
			real_xyz l = 0.0;

			for (size_t e = 0; e < eipn; ++e) { //MK::e is the threadlocal iterating element id, NOT the global ip ID !

				real_xyz xx = thisregion->grid.xyz0[e].x + thisregion->grid.dxyz_avg[e].dx + thisregion->grid.dxyz_flu[e].dx;	//##MK::SIMD vectorize location of element in deformed configuration
				real_xyz yy = thisregion->grid.xyz0[e].y + thisregion->grid.dxyz_avg[e].dy + thisregion->grid.dxyz_flu[e].dy;
				real_xyz zz = thisregion->grid.xyz0[e].z + thisregion->grid.dxyz_avg[e].dz + thisregion->grid.dxyz_flu[e].dz;
				i = 0; //defining at which position you buffer these intermediate resutls in the caches cahce, targetpin and periodicImage

				//most points lay inside
				//<000>
				h = +0.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x00; i++; }

				//<100>
				//[+1+0+0]
				h = +1.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0e; i++; }
				//[-1+0+0]
				h = -1.0;		k = +0.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0d; i++; }
				//[+0+1+0]
				h = +0.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x10; i++; }
				//[+0-1+0]
				h = +0.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0b; i++; }
				//[+0+0+1]
				h = +0.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x16; i++; }
				//[+0+0-1]
				h = +0.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x05; i++; }

				//<110>
				//[+1+1+0]
				h = +1.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x11; i++; }
				//[-1+1+0]
				h = -1.0;		k = +1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0f; i++; }
				//[+1-1+0]
				h = +1.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0c; i++; }
				//[-1-1+0]
				h = -1.0;		k = -1.0;		l = +0.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x0a; i++; }
				//[+0+1+1]
				h = +0.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x19; i++; }
				//[+0-1+1]
				h = +0.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x13; i++; }
				//[+0+1-1]
				h = +0.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x08; i++; }
				//[+0-1-1]
				h = +0.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x02; i++; }
				//[+1+0+1]
				h = +1.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x17; i++; }
				//[-1+0+1]
				h = -1.0;		k = +0.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x15; i++; }
				//[+1+0-1]
				h = +1.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x06; i++; }
				//[-1+0-1]
				h = -1.0;		k = +0.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x04; i++; }

				//<111>
				//[+1+1+1]
				h = +1.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x1a; i++; }
				//[-1+1+1]
				h = -1.0;		k = +1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x18; i++; }
				//[+1-1+1]
				h = +1.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x14; i++; }
				//[+1+1-1]
				h = +1.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x09; i++; }
				//[+1-1-1]
				h = +1.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x03; i++; }
				//[-1+1-1]
				h = -1.0;		k = +1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x07; i++; }
				//[-1-1+1]
				h = -1.0;		k = -1.0;		l = +1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x12; i++; }
				//[-1-1-1]
				h = -1.0;		k = -1.0;		l = -1.0;
				xxx = xx + (h*per.a11 + k*per.a12 +	l*per.a13);
				yyy = yy + (h*per.a21 + k*per.a22 + l*per.a23);
				zzz = zz + (h*per.a31 + k*per.a32 +	l*per.a33);
				if (xxx >= w.xmi && xxx <= w.xmx && yyy >= w.ymi && yyy <= w.ymx && zzz >= w.zmi && zzz <= w.zmx ) { cache[i] = p3d(xxx,yyy,zzz); periodicImage[i] = 0x01; i++; }

				//now that we know how many periodic copies of e there are at all inside namely i we bin them
				for (unsigned int c = 0; c < i; ++c) {
					targetbin[c] = binning_xyz( cache[c] );
//cout << targetbin[c] << "__";
				}
//cout << endl;

				//check if still i places free in buffer construct, if so store results
				if ( mybuffer.back()->size() < (Performance::ThreadBufferAllocSize - static_cast<size_t>(i)) ) {
					vector<p3dm3>* drophere = mybuffer.back();
					unsigned int eid = static_cast<unsigned int>(e + lc2gl_idoffset);
					for (unsigned int c = 0; c < i; ++c) {
						drophere->push_back( p3dm3(cache[c].x, cache[c].y, cache[c].z,
								eid, targetbin[c], periodicImage[c] ) );
						//MK::mark m1 is the global ID of the ip,
						//m2 the ID of the bin in which the point should be stored,
						//m3 metadata telling us in which periodic BC relation the point is to its generating ip
					}
				} //existent buffer filled
				else { //get new buffer
					try {
						s = NULL;
						s = new vector<p3dm3>;
						mybuffer.push_back(s);
					}
					catch (bad_alloc &ompcroak) {
						stopping("Unable to reallocate additional threadlocal binning memory", owner->get_myrank(), tid);
						mehealthy = false;
						break; //stop working on this thread
					}

					vector<p3dm3>* dropwhere = mybuffer.back();
					unsigned int eid = static_cast<unsigned int>(e + lc2gl_idoffset); //MK::ip base and its periodic images have the same global ip ID eid
					for (unsigned int c = 0; c < i; ++c) {
						dropwhere->push_back( p3dm3(cache[c].x, cache[c].y, cache[c].z,
								eid, targetbin[c], 0x00 ) );
					}
				} //new buffer utilized for filling
			} //process next element e of thread
		} //while being healthy
		else { //upon initial allocation error
			stopping("Unable to initially allocate threadlocal binning memory", owner->get_myrank(), tid);
			mehealthy = false;
		}

		#pragma omp critical
		{
			if ( mehealthy == false ) {
				healthy = false;
			}
			reporting("Bounding volume hierarchy construction threadlocal results ready", owner->get_myrank(), tid, true);
		}

		//wait for all finishing the positioning and binning of the periodic images in w
		#pragma omp barrier


		//by now the threads have positioned, checked for inclusion, and binned all mesh element integration points surplus their periodic images points
		#pragma omp critical
		{
			if ( mehealthy == false ) { //in case previous memory allocation was unsuccessful
				for (size_t b = 0; b < mybuffer.size(); ++b) {
					if (mybuffer.at(b) != NULL) {
						delete mybuffer.at(b); mybuffer.at(b) = NULL;
					}
				}
			}
			else { //valid data, so pass them into global bvh
				for (size_t b = 0; b < mybuffer.size(); ++b) {
					if (mybuffer.at(b) != NULL) {
						vector<p3dm3>* these = mybuffer.at(b);
						size_t nj = these->size();
						for (size_t j = 0; j < nj; ++j) {
							//MK::m1 is the global ID of the ip,
							//m2 the ID of the bin in which the point should be stored,
							//m3 metadata telling us in which periodic BC relation the point is to its generating ip
							p3dm3 piece = these->at(j);
							px.at(piece.m2)->push_back( piece.x );
							py.at(piece.m2)->push_back( piece.y );
							pz.at(piece.m2)->push_back( piece.z );
							pid.at(piece.m2)->push_back( piece.m1 );
							//##MK::DEBUG which Moore periodic image is this?
							pimage.at(piece.m2)->push_back( piece.m3 );
						}
						//release threadlocal temporary mybuffer
						delete mybuffer.at(b); mybuffer.at(b) = NULL;
					}
				} //copy over next group of elements to process
			}
		} //end of pragma critical, thread tid has communicated local results
	} //end of parallel region

	if ( healthy == true )
		reporting("Bounding volume hierarchy construction complete", owner->get_myrank(), 0, false);
	else
		reporting("Bounding volume hierarchy construction failed", owner->get_myrank(), 0, false);
*/
}
