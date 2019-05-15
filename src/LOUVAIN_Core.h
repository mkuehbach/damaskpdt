//File: louvain.h
//-- community detection header file
//-----------------------------------------------------------------------------
//Community detection
//Based on the article "Fast unfolding of community hierarchies in large networks"
//Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
//And based on the article "A Generalized and Adaptive Method for Community Detection"
//Copyright (C) 2014 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
//This file is part of Louvain algorithm.
//
//Louvain algorithm is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Louvain algorithm is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
//Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
//Email    : jean-loup.guillaume@lip6.fr
//Location : Paris, France
//Time	    : July 2014
//-----------------------------------------------------------------------------
//see readme.txt for more details

#ifndef __LOUVAIN_CORE_H__
#define __LOUVAIN_CORE_H__


#include "LOUVAIN_GraphBinary.h"

//using namespace std;


struct lvlogiter
{
	//##MK::alignment is poor
	int levelid;
	size_t mvcnt;
	long double cur_qual;
	long double new_qual;
	double dt;
	lvlogiter() : levelid(static_cast<int>(0)), mvcnt(static_cast<size_t>(0)),
			cur_qual(static_cast<long double>(0.0)), new_qual(static_cast<long double>(0.0)),
			dt(static_cast<double>(0.0)) {}
	lvlogiter(const int _lid, const size_t _mvs, const long double _cq, const long double _nq, const double _dt) :
		levelid(_lid), mvcnt(_mvs), cur_qual(_cq), new_qual(_nq), dt(_dt) {}
};

struct lvloglevel
{
	int level;
	int nbnodes;
	unsigned long long nblinks;
	long double tweight;
	lvloglevel() : level(static_cast<int>(0)), nbnodes(static_cast<int>(0)),
			nblinks(static_cast<unsigned long long>(0)), tweight(static_cast<long double>(0.0)) {}
	lvloglevel(const int _lvl, const int _nbn, const unsigned long long _nbl, const long double _twght) :
		level(_lvl), nbnodes(_nbn), nblinks(_nbl), tweight(_twght) {}
};


class Louvain
{
public:
	vector<real_louvain> neigh_weight;
	vector<int> neigh_pos;
	int neigh_last;

	//number of pass for one level computation
	//if -1, compute as many pass as needed to increase quality
	int nb_pass;

	//a new pass is computed if the last one has generated an increase 
	//better than eps_impr
	//if 0.0L even a minor increase is enough to go for one more pass
	real_louvain eps_impr;

	//Quality functions used to compute communities
	Quality* qual;


	//constructors:
	//reads graph from file using graph constructor
	//type defined the weighted/unweighted status of the graph file
	Louvain (int nb_pass, real_louvain eps_impr, Quality* q);

	//initiliazes the partition with something else than all nodes alone
	//void init_partition(char *filename_part);

	//compute the set of neighboring communities of node
	//for each community, gives the number of links from node to comm
	void neigh_comm(int node);

	//displays the graph of communities as computed by one_level
	void partition2graph();

	//displays the current partition (with communities renumbered from 0 to k-1)
	void display_partition();

	vector<int>* pass_partition();

	//generates the binary graph of communities as computed by one_level
	Graph partition2graph_binary();

	//compute communities of the graph for one level
	//return true if some nodes have been moved
	bool one_level( const int thislevel, vector<lvlogiter> & tctc );

	//vector<lvlogiter> tctc1;
	//vector<lvloglevel> tctc2;
};


class louvainHdl
{
public:
	louvainHdl();
	~louvainHdl();

	void init_quality(Graph *g, unsigned short nbc);
	void execute(vector<lvwtedge> const & these_edges, vector<unsigned int> & uip2comm);
	void get_tctc( Louvain const & fromhere);
	void profiling( const unsigned int thisrank, const unsigned int thisincr );

	real_louvain precision;
	real_louvain alpha;
	real_louvain sum_se;
	real_louvain sum_sq;
	real_louvain max_w;

	Quality *q;

	int wtype;
	int nb_pass;
	int display_level;
	int kmin;

	unsigned short id_qual;
	bool verbose;

	//mt19937 myrandom;

	vector<lvlogiter> tictoc1;
	vector<lvloglevel> tictoc2;
};

#endif
