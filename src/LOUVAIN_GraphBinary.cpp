//File: graph_binary.cpp
//-- graph handling source
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



#include "LOUVAIN_GraphBinary.h"


Graph::Graph()
{
	nb_nodes = 0;
	nb_lnks = 0ULL;
	nb_links = 0ULL;

	total_weight = 0.0L;
	sum_nodes_w = 0;
}


Graph::Graph(vector<lvwtedge> edg)
{
	//##MK::by definition working on weighted graph
	int wtype = WEIGHTED;

	//transfer step from PDT_DAMASK data structure into louvain
	//takes care of the step reading in a text file into converter help app of louvain
	unsigned long long nb_lnks = 0ULL;

	for (size_t i = 0; i < edg.size(); ++i) {
		unsigned int src = edg[i].src;
		unsigned int dest = edg[i].dest;
		real_louvain weight = static_cast<real_louvain>(edg[i].wt);

		if (lnks.size() <= max(src,dest)+1) {
			lnks.resize(max(src,dest)+1);
		}

		lnks[src].push_back(make_pair(dest,weight));
		if (src!=dest)
			lnks[dest].push_back(make_pair(src,weight));

		nb_lnks += 1ULL;
	}

	//##MK::no renumbering

	//##MK::cleaning procedure
	for (unsigned int i=0 ; i < lnks.size() ; i++) {
		map<int, real_louvain> m;
		map<int, real_louvain>::iterator it;

		for (unsigned int j=0 ; j<lnks[i].size() ; j++) {
			it = m.find(lnks[i][j].first);
			if (it==m.end())
				m.insert(make_pair(lnks[i][j].first, lnks[i][j].second));
			else if (wtype == WEIGHTED)
				it->second+=lnks[i][j].second;
		}

		vector<pair<int, real_louvain> > v;
		for (it = m.begin(); it != m.end(); it++)
			v.push_back(*it);

		lnks[i].clear();
		lnks[i] = v;
	} //for each key src

	//##MK::no displaying, i.e. writing of file from converter but in-place passing to louvain app graph class object
	//transfer to graph_binary structure class object
	int s = lnks.size();
	nb_nodes = s;

//cout << "nb_nodes/s" << endl;
//cout << nb_nodes << "\t\t" << s << endl;

	//transfer cumulative degree sequence
	//Read cumulative degree sequence: 8 bytes for each node
	//cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.

//cout << "degrees" << endl;
//cout << "initial degrees.size()/lnks.size() " << degrees.size() << "\t\t" << lnks.size() << endl;

	//degrees.reserve(nb_nodes);
	unsigned long long tot = 0ULL;
	for (int i = 0; i < s; i++) {
		tot += (unsigned long long)lnks[i].size();
		degrees.push_back(tot);
	}
//for(unsigned int k = 0; k < degrees.size(); k++) { cout << k << "--->" << degrees[k] << endl; }
//cout << "degrees.size() " << degrees.size() << endl;

	//transfer links
	//Read links: 4 bytes for each link (each link is counted twice)
//cout << "nb_links" << endl;

	nb_links = degrees[nb_nodes-1];
	unsigned int k = 0;
	for(int i = 0; i < s; i++) {
		for (unsigned int j = 0; j < lnks[i].size(); j++) {
			links.push_back( lnks[i][j].first );
//cout << k << "--->" << links.back() << endl;
			k++;
		}
	}
//cout << "links.size() " << links.size() << endl;

	//##MK::always IF WEIGHTED, read weights: 10 bytes for each link (each link is counted twice)
//cout << "weights" << endl;
	//weights.resize(0);
	//if (wtype == WEIGHTED) {
	//weights.reserve(nb_links);
	k = 0;
	for(int i = 0; i < s; i++) {
		for (unsigned int j = 0; j < lnks[i].size(); j++) {
			weights.push_back( lnks[i][j].second );
//cout << k << "--->" << weights.back() << endl;
			k++;
		}
	}
//cout << "weights.size() " << weights.size() << endl;

	//Compute total weight
//cout << "total_weight" << endl;
	total_weight = 0.0L;
	for (int i = 0; i < nb_nodes; i++) {
		total_weight += (real_louvain)weighted_degree(i); //##MK::cast is unnecessary...
	}
//cout << total_weight << endl;
//cout << endl << endl << endl;

	nodes_w.assign(nb_nodes, 1);
	sum_nodes_w = nb_nodes;
}

/*
Graph::Graph(char *filename, char *filename_w, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, sizeof(int));
  if (finput.rdstate() != ios::goodbit) {
    cerr << "The file " << filename << " is not a valid graph" << endl;
    exit(EXIT_FAILURE);
  }
  
  //Read cumulative degree sequence: 8 bytes for each node
  //cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);
  finput.read((char *)&degrees[0], nb_nodes*sizeof(unsigned long long));

  //Read links: 4 bytes for each link (each link is counted twice)
  nb_links = degrees[nb_nodes-1];
  links.resize(nb_links);
  finput.read((char *)(&links[0]), nb_links*sizeof(int));

  //IF WEIGHTED, read weights: 10 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight = 0.0L;
  if (type==WEIGHTED) {
    ifstream finput_w;
    finput_w.open(filename_w,fstream::in | fstream::binary);
    if (finput_w.is_open() != true) {
      cerr << "The file " << filename_w << " does not exist" << filename << endl;
      exit(EXIT_FAILURE);
    }

    weights.resize(nb_links);
    finput_w.read((char *)(&weights[0]), nb_links*sizeof(real_louvain));
    if (finput_w.rdstate() != ios::goodbit) {
      cerr << "The file " << filename_w << " does not correspond to valid weights for the graph" << filename << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Compute total weight
  for (int i=0 ; i<nb_nodes ; i++)
    total_weight += (real_louvain)weighted_degree(i);

  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}

real_louvain
Graph::max_weight() {
  real_louvain max = 1.0L;

  if (weights.size()!=0)
    max = *max_element(weights.begin(),weights.end());
  
  return max;
}
*/


void Graph::assign_weight(int node, int weight)
{
	sum_nodes_w -= nodes_w[node];

	nodes_w[node] = weight;

	sum_nodes_w += weight;
}

void Graph::add_selfloops()
{
	vector<unsigned long long> aux_deg;
	vector<int> aux_links;

	unsigned long long sum_d = 0ULL;

	for (int u=0 ; u < nb_nodes ; u++) {
		pair<vector<int>::iterator, vector<real_louvain>::iterator> p = neighbors(u);
		int deg = nb_neighbors(u);

		for (int i=0 ; i < deg ; i++) {
			int neigh = *(p.first+i);
			aux_links.push_back(neigh);
		}

		sum_d += (unsigned long long)deg;

		if (nb_selfloops(u) == 0.0L) {
			aux_links.push_back(u); //add a selfloop
			sum_d += 1ULL;
		}

		aux_deg.push_back(sum_d); //add the (new) degree of vertex u
	}

	links = aux_links;
	degrees = aux_deg;

	nb_links += (unsigned long long)nb_nodes;
}

void Graph::display()
{
	for (int node=0 ; node<nb_nodes ; node++) {
		pair<vector<int>::iterator, vector<real_louvain>::iterator > p = neighbors(node);
		cout << node << ":" ;
		for (int i=0 ; i<nb_neighbors(node) ; i++) {
			if (true) {
				if (weights.size()!=0)
					cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
				else
					cout << " " << *(p.first+i);
			}
		}
		cout << endl;
	}
}

void Graph::display_reverse()
{
	for (int node=0 ; node<nb_nodes ; node++) {
		pair<vector<int>::iterator, vector<real_louvain>::iterator > p = neighbors(node);
		for (int i=0 ; i<nb_neighbors(node) ; i++) {
			if (node>*(p.first+i)) {
				if (weights.size()!=0)
					cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
				else
					cout << *(p.first+i) << " " << node << endl;
			}
		}
	}
}

/*
bool
Graph::check_symmetry() {
  int error = 0;
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<real_louvain>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      int neigh = *(p.first+i);
      real_louvain weight = *(p.second+i);

      pair<vector<int>::iterator, vector<real_louvain>::iterator > p_neigh = neighbors(neigh);
      for (int j=0 ; j<nb_neighbors(neigh) ; j++) {
	int neigh_neigh = *(p_neigh.first+j);
	real_louvain neigh_weight = *(p_neigh.second+j);

	if (node==neigh_neigh && weight!=neigh_weight) {
	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
	  if (error++==10)
	    exit(0);
	}
      }
    }
  }
  return (error==0);
}
*/


void Graph::display_binary(char *outfile)
{
	ofstream foutput;
	foutput.open(outfile ,fstream::out | fstream::binary);

	foutput.write((char *)(&nb_nodes),sizeof(int));
	foutput.write((char *)(&degrees[0]),sizeof(unsigned long long)*nb_nodes);
	foutput.write((char *)(&links[0]),sizeof(int)*nb_links);
}



Quality::~Quality()
{
	n2c.clear();
}



Modularity::Modularity(Graph & gr) :Quality(gr, "Newman-Girvan Modularity")
{
	n2c.resize(size);

	in.resize(size);
	tot.resize(size);

	//initialization
	for (int i = 0; i<size; i++) {
		n2c[i] = i;
		in[i] = g.nb_selfloops(i);
		tot[i] = g.weighted_degree(i);
	}
}

Modularity::~Modularity()
{
	in.clear();
	tot.clear();
}

real_louvain Modularity::quality()
{
	real_louvain q = 0.0L;
	real_louvain m2 = g.total_weight;

	for (int i = 0; i<size; i++) {
		if (tot[i] > 0.0L)
			q += in[i] - (tot[i] * tot[i]) / m2;
	}

	q /= m2;

	return q;
}

