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

#ifndef __PDT_PROFILER_H__
#define __PDT_PROFILER_H__

#include "PDT_Verbose.h"


//program profiling should use double precision in general as
//MPI_Wtime() and omp_get_wtime() fires in double precision

//type of computational operations
#define APT_XX				0		//default, unspecified
#define APT_UTL				1		//utility
#define APT_SPE_IO			2		//SpectralOut I/O
#define APT_MPI_IO			3		//MPI I/O
#define APT_VIS_IO			4
#define APT_GRA				5		//grain recon
#define APT_DIS				6		//disori
#define APT_SDF				7		//signed distance
#define APT_TES				8		//tessellation
#define APT_COR				9		//actual distance correlations


#define APT_IS_UNKNOWN		-1
#define APT_IS_SEQ			0
#define APT_IS_PAR			1		//information telling that parallelism is used

#define MEMORY_NOSNAPSHOT_TAKEN		-1


struct memsnapshot
{
	size_t virtualmem;
	size_t residentmem;
	memsnapshot() : virtualmem(MEMORY_NOSNAPSHOT_TAKEN),
			residentmem(MEMORY_NOSNAPSHOT_TAKEN) {}
	memsnapshot(const size_t _vm, const size_t _rm) :
		virtualmem(_vm), residentmem(_rm) {}
};


struct plog
{
	double dt;
	double tstart;
	double tend;
	size_t virtualmem;		//virtual memory consumption in bytes
	size_t residentmem;		//resident set size in bytes, i.e. number of pages process as in real memory times system specific page size
	string what;
	unsigned short typ;		//task identifier
	unsigned short pll;		//parallelism identifier
	unsigned int i;			//running number to identify the which-th snapshot
							//used because for easier utilizability of the result
							//we sort in ascending processing time thereby however having the mem data not as a time trajectory
	unsigned int tskid;		//to distinguish for the individual increments

	plog() : dt(0.0), tstart(0.0), tend(0.0), virtualmem(-1), residentmem(-1),
			what(""), typ(APT_XX), pll(APT_IS_SEQ), i(0), tskid(0) {}
	plog(const double _dt, const size_t _vm, const size_t _rm, const string _s,
			const unsigned short _t, const unsigned short _p, const unsigned int _i, const unsigned int _tskid) :
				dt(_dt), tstart(0.0), tend(0.0), virtualmem(_vm), residentmem(_rm),
				what(_s), typ(_t), pll(_p), i(_i), tskid(_tskid) {}
	plog(const double _ts, const double _te, const string _s,
			const unsigned short _t, const unsigned short _p, const unsigned int _i, const unsigned int _tskid) :
				dt(_te - _ts), tstart(_ts), tend(_te), virtualmem(-1), residentmem(-1),
				what(_s), typ(_t), pll(_p), i(_i), tskid(_tskid) {}	//version not tracking memory consumption
	plog(const double _ts, const double _te, const size_t _vm, const size_t _rm,
			const string _s, const unsigned short _t, const unsigned short _p,
			const unsigned int _i, const unsigned int _tskid) :
				dt(_te - _ts), tstart(_ts), tend(_te), virtualmem(_vm), residentmem(_rm),
				what(_s), typ(_t), pll(_p), i(_i), tskid(_tskid) {}	//version tracking memory consumption
};

class profiler
{
public:
	profiler() {};
	~profiler() {};

	//void prof(const string whichenv, const unsigned short category, const double st, const double en);
	void prof_elpsdtime_and_mem(const string whichenv, const unsigned short category,
			const unsigned short parallelism, memsnapshot const & mem, const double st, const double en,
			const unsigned int taskid );
	void prof_elpsdtime_and_mem(const string whichenv, const unsigned short category,
				const unsigned short parallelism, memsnapshot const & mem, const double dt,	const unsigned int taskid );
		memsnapshot get_memoryconsumption( void );
	size_t get_nentries( void );
	void report_memory( pair<size_t,size_t> const & in );
	void spit_profiling( const unsigned int simid, const int rank );

private:
	vector<plog> evn;
};


#endif
