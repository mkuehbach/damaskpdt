//File: graph_binary.h
//-- graph handling header file
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
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
//see README.txt for more details


#ifndef __LOUVAIN_GRAPH_BINARY_H__
#define __LOUVAIN_GRAPH_BINARY_H__

#include "PDT_GrainObject.h"
/*
#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <vector>

#include <unistd.h>

#include<random>
#include<functional>
using namespace std;
*/


#define WEIGHTED   0
//#define UNWEIGHTED 1




//test utilizing double
typedef long double	real_louvain;

//##MK::at the moment does not work...
//typedef double 		real_louvain;

/*
struct lvwtedge
{
	unsigned int src;		//ip i
	unsigned int dest;		//ip j
	double wt;				//weight of edge ij

	lvwtedge() : src(numeric_limits<unsigned int>::max()), dest(numeric_limits<unsigned int>::max()), wt(0.0) {}
	lvwtedge(const unsigned int _src, const unsigned int _dst, const double _wt) :
		src(_src), dest(_dst), wt(_wt) {}
	~lvwtedge(){}
};
*/


class Graph
{
public:
	vector<vector<pair<int, real_louvain> > > lnks;

	int nb_nodes;
	unsigned long long nb_lnks;
	unsigned long long nb_links;

	real_louvain total_weight;
	int sum_nodes_w;

	vector<unsigned long long> degrees;
	vector<int> links;
	vector<real_louvain> weights;

	vector<int> nodes_w;

	Graph();

	//binary file format is
	//4 bytes for the number of nodes in the graph
	//8*(nb_nodes) bytes for the cumulative degree for each node:
	// deg(0)=degrees[0]
	// deg(k)=degrees[k]-degrees[k-1]
	//4*(sum_degrees) bytes for the links
	//IF WEIGHTED, 10*(sum_degrees) bytes for the weights in a separate file
	Graph(vector<lvwtedge> edg);
	//Graph(char *filename, char *filename_w, int type);

	//return the biggest weight of links in the graph
	//real_louvain max_weight();

	//assign a weight to a node (needed after the first level)
	void assign_weight(int node, int weight);

	//add selfloop to each vertex in the graph
	void add_selfloops();

	void display(void);
	void display_reverse(void);
	void display_binary(char *outfile);
	//bool check_symmetry();

	//return the number of neighbors (degree) of the node
	inline int nb_neighbors(int node);

	//return the number of self loops of the node
	inline real_louvain nb_selfloops(int node);

	//return the weighted degree of the node
	inline real_louvain weighted_degree(int node);

	//return pointers to the first neighbor and first weight of the node
	inline pair<vector<int>::iterator, vector<real_louvain>::iterator > neighbors(int node);
};


inline int Graph::nb_neighbors(int node)
{
	//assert(node>=0 && node<nb_nodes);

	if (node==0)
		return degrees[0];
	else
		return (int)(degrees[node]-degrees[node-1]);
}

inline real_louvain Graph::nb_selfloops(int node)
{
	//assert(node>=0 && node<nb_nodes);

	pair<vector<int>::iterator, vector<real_louvain>::iterator > p = neighbors(node);
	for (int i=0 ; i<nb_neighbors(node) ; i++) {
		if (*(p.first+i)==node) {
			if (weights.size()!=0)
				return (real_louvain)*(p.second+i);
			else 
				return 1.0L;
		}
	}
	return 0.0L;
}

inline real_louvain Graph::weighted_degree(int node)
{
	//assert(node>=0 && node<nb_nodes);

	if (weights.size()==0)
		return (real_louvain)nb_neighbors(node);
	else {
		pair<vector<int>::iterator, vector<real_louvain>::iterator > p = neighbors(node);
		real_louvain res = 0.0L;
		for (int i=0 ; i<nb_neighbors(node) ; i++) {
			res += (real_louvain)*(p.second+i);
		}
		return res;
	}
}

inline pair<vector<int>::iterator, vector<real_louvain>::iterator > Graph::neighbors(int node)
{
	//assert(node>=0 && node<nb_nodes);

	if (node==0)
		return make_pair(links.begin(), weights.begin());
	else if (weights.size()!=0)
		return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
	else
		return make_pair(links.begin()+degrees[node-1], weights.begin());
}


class Quality
{
public:
	Graph & g; //network to compute communities for
	int size; //nummber of nodes in the network and size of all vectors
	string name;

	vector<int> n2c; //community to which each node belongs
	Quality(Graph &gr, const std::string& n) :g(gr), size(g.nb_nodes), name(n) {}

	virtual ~Quality();

	//remove the node from its current community with which it has dnodecomm links
	virtual void remove(int node, int comm, real_louvain dnodecomm) = 0;

	//insert the node in comm with which it shares dnodecomm links
	virtual void insert(int node, int comm, real_louvain dnodecomm) = 0;

	//compute the gain of quality by adding node to comm
	virtual real_louvain gain(int node, int comm, real_louvain dnodecomm, real_louvain w_degree) = 0;

	//compute the quality of the current partition
	virtual real_louvain quality() = 0;
};


class Modularity : public Quality
{
public:

	vector<real_louvain> in, tot; //used to compute the quality participation of each community

	Modularity(Graph & gr);
	~Modularity();

	inline void remove(int node, int comm, real_louvain dnodecomm);

	inline void insert(int node, int comm, real_louvain dnodecomm);

	inline real_louvain gain(int node, int comm, real_louvain dnodecomm, real_louvain w_degree);

	real_louvain quality();
};


inline void Modularity::remove(int node, int comm, real_louvain dnodecomm)
{
	//assert(node >= 0 && node<size);

	in[comm] -= 2.0L*dnodecomm + g.nb_selfloops(node);
	tot[comm] -= g.weighted_degree(node);

	n2c[node] = -1;
}

inline void Modularity::insert(int node, int comm, real_louvain dnodecomm)
{
	//assert(node >= 0 && node<size);

	in[comm] += 2.0L*dnodecomm + g.nb_selfloops(node);
	tot[comm] += g.weighted_degree(node);

	n2c[node] = comm;
}

inline real_louvain Modularity::gain(int node, int comm, real_louvain dnc, real_louvain degc)
{
	//assert(node >= 0 && node<size);

	real_louvain totc = tot[comm];
	real_louvain m2 = g.total_weight;

	return (dnc - totc*degc / m2);
}

#endif
