//File: louvain.cpp
//-- community detection source file
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

#include "LOUVAIN_Core.h"

//using namespace std;


Louvain::Louvain(int nbp, real_louvain epsq, Quality* q) {
	qual = q;

	neigh_weight.resize(qual->size,-1);
	neigh_pos.resize(qual->size);
	neigh_last = 0;

	nb_pass = nbp;
	eps_impr = epsq;
}

/*
void Louvain::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  //read partition
  while (!finput.eof()) {
    int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = qual->n2c[node];
      neigh_comm(node);

      qual->remove(node, old_comm, neigh_weight[old_comm]);
      
      int i=0;
      for (i=0 ; i<neigh_last ; i++) {
	int best_comm = neigh_pos[i];
	real_louvain best_nblinks = neigh_weight[neigh_pos[i]];
	if (best_comm==comm) {
	  qual->insert(node, best_comm, best_nblinks);
	  break;
	}
      }
      if (i==neigh_last)
	qual->insert(node, comm, 0);
    }
  }
  finput.close();
}
*/

void Louvain::neigh_comm(int node)
{
	for (int i=0 ; i<neigh_last ; i++)
		neigh_weight[neigh_pos[i]]=-1;

	neigh_last = 0;

	pair<vector<int>::iterator, vector<real_louvain>::iterator> p = (qual->g).neighbors(node);
	int deg = (qual->g).nb_neighbors(node);

	neigh_pos[0] = qual->n2c[node];
	neigh_weight[neigh_pos[0]] = 0;
	neigh_last = 1;

	for (int i=0 ; i<deg ; i++) {
		int neigh  = *(p.first+i);
		int neigh_comm = qual->n2c[neigh];
		real_louvain neigh_w = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);

		if (neigh!=node) {
			if (neigh_weight[neigh_comm]==-1) {
				neigh_weight[neigh_comm] = 0.0L;
				neigh_pos[neigh_last++] = neigh_comm;
			}
			neigh_weight[neigh_comm] += neigh_w;
		}
	}
}

void Louvain::partition2graph()
{
	vector<int> renumber(qual->size, -1);
	for (int node=0 ; node<qual->size ; node++) {
		renumber[qual->n2c[node]]++;
	}

	int end=0;
	for (int i=0 ; i< qual->size ; i++)
		if (renumber[i]!=-1)
			renumber[i]=end++;

	for (int i=0 ; i< qual->size ; i++) {
		pair<vector<int>::iterator, vector<real_louvain>::iterator> p = (qual->g).neighbors(i);
		int deg = (qual->g).nb_neighbors(i);
		for (int j=0 ; j<deg ; j++) {
			int neigh = *(p.first+j);
			cout << renumber[qual->n2c[i]] << " " << renumber[qual->n2c[neigh]] << endl;
		}
	}
}

void Louvain::display_partition()
{
	vector<int> renumber(qual->size, -1);
	for (int node=0 ; node < qual->size ; node++) {
		renumber[qual->n2c[node]]++;
	}

	int end=0;
	for (int i=0 ; i < qual->size ; i++)
	if (renumber[i]!=-1)
		renumber[i] = end++;

	for (int i=0 ; i < qual->size ; i++)
		cout << i << " " << renumber[qual->n2c[i]] << endl;
}


vector<int>* Louvain::pass_partition()
{
	vector<int>* bucket = NULL;
	bucket = new vector<int>;
	//##MK::testing

	vector<int> renumber(qual->size, -1);
	for (int node = 0; node < qual->size; node++) {
		renumber[qual->n2c[node]]++;
	}

	int end = 0;
	for (int i = 0; i < qual->size; i++)
		if (renumber[i] != -1)
			renumber[i] = end++;

	bucket->reserve(qual->size);
	for (int i = 0; i < qual->size; i++) {
		//cout << i << " " << renumber[qual->n2c[i]] << endl;
		bucket->push_back(renumber[qual->n2c[i]]);
	}
	return bucket;
}


Graph Louvain::partition2graph_binary()
{
	//Renumber communities
	vector<int> renumber(qual->size, -1);
	for (int node=0 ; node < qual->size ; node++)
		renumber[qual->n2c[node]]++;

	int last=0;
	for (int i=0 ; i < qual->size ; i++) {
		if (renumber[i]!=-1)
			renumber[i] = last++;
	}

	//Compute communities
	vector<vector<int> > comm_nodes(last);
	vector<int> comm_weight(last, 0);

	for (int node = 0 ; node < (qual->size) ; node++) {
		comm_nodes[renumber[qual->n2c[node]]].push_back(node);
		comm_weight[renumber[qual->n2c[node]]] += (qual->g).nodes_w[node];
	}

	//Compute weighted graph
	Graph g2;
	int nbc = comm_nodes.size();

	g2.nb_nodes = comm_nodes.size();
	g2.degrees.resize(nbc);
	g2.nodes_w.resize(nbc);

	for (int comm=0 ; comm<nbc ; comm++) {
		map<int,real_louvain> m;
		map<int,real_louvain>::iterator it;

		int size_c = comm_nodes[comm].size();

		g2.assign_weight(comm, comm_weight[comm]);

		for (int node=0 ; node<size_c ; node++) {
			pair<vector<int>::iterator, vector<real_louvain>::iterator> p = (qual->g).neighbors(comm_nodes[comm][node]);
			int deg = (qual->g).nb_neighbors(comm_nodes[comm][node]);
			for (int i=0 ; i<deg ; i++) {
				int neigh = *(p.first+i);
				int neigh_comm = renumber[qual->n2c[neigh]];
				real_louvain neigh_weight = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);

				it = m.find(neigh_comm);
				if (it==m.end())
					m.insert(make_pair(neigh_comm, neigh_weight));
				else
					it->second += neigh_weight;
			}
		}

		g2.degrees[comm] = (comm==0)?m.size():g2.degrees[comm-1]+m.size();
		g2.nb_links += m.size();

		for (it = m.begin() ; it!=m.end() ; it++) {
			g2.total_weight += it->second;
			g2.links.push_back(it->first);
			g2.weights.push_back(it->second);
		}
	}

	return g2;
}

bool Louvain::one_level( const int thislevel, vector<lvlogiter> & tctc )
{
	bool improvement=false;
	int nb_moves;
	int nb_pass_done = 0;
	real_louvain new_qual = qual->quality();
	real_louvain cur_qual = new_qual;

	vector<int> random_order(qual->size);
	for (int i=0 ; i < qual->size ; i++)
		random_order[i]=i;
	for (int i=0 ; i < qual->size-1 ; i++) {
		//##MK::class structure design of Louvain is very poor
		//int rand_pos = bind(uniform_int_distribution<int>(0, qual->size - i), myrandom()) + i;
		int rand_pos = rand()%(qual->size-i)+i;
		int tmp = random_order[i];
		random_order[i] = random_order[rand_pos];
		random_order[rand_pos] = tmp;
	}

	//repeat while
	//there is an improvement of quality
	//or there is an improvement of quality greater than a given epsilon
	//or a predefined number of pass have been done

	//MK::empirical observation of the labeling algorithm showed very slow convergence because of a much too low default eps_impr for the actual convergence rate
	//after an initial power law convergence behavior in the initial approx 20 iterations
	//hence we store instead the
	real_louvain new_eps_impr = 0.0L;
	do {
		double tic = MPI_Wtime();

		cur_qual = new_qual;
		nb_moves = 0;
		nb_pass_done++;

		//for each node: remove the node from its community and insert it in the best community
		for (int node_tmp = 0 ; node_tmp < qual->size ; node_tmp++) {
			int node = random_order[node_tmp];
			int node_comm = qual->n2c[node];
			real_louvain w_degree = (qual->g).weighted_degree(node);

			//computation of all neighboring communities of current node
			neigh_comm(node);
			//remove node from its current community
			qual->remove(node, node_comm, neigh_weight[node_comm]);

			//compute the nearest community for node
			//default choice for future insertion is the former community
			int best_comm = node_comm;
			real_louvain best_nblinks  = 0.0L;
			real_louvain best_increase = 0.0L;
			for (int i=0 ; i<neigh_last ; i++) {
				real_louvain increase = qual->gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
				if (increase>best_increase) {
					best_comm = neigh_pos[i];
					best_nblinks = neigh_weight[neigh_pos[i]];
					best_increase = increase;
				}
			}

			//insert node in the nearest community
			qual->insert(node, best_comm, best_nblinks);

			if (best_comm!=node_comm)
				nb_moves++;
//cout << "Louvain processing node_tmp " << node_tmp << " of " << qual->size << endl;
		}

		new_qual = qual->quality();

		if (nb_moves>0)
			improvement=true;

		double toc = MPI_Wtime();
		tctc.push_back( lvlogiter( thislevel, nb_moves, cur_qual, new_qual, toc-tic ) );

cout << "Thislevel, NumberOfMoves, new qual-curr qual, eps_impr, new_eps_impr " << thislevel << ";" << nb_moves << "\t\t" << (new_qual - cur_qual) << ";" << eps_impr << "\t\t" << new_eps_impr << endl;

		if ( nb_pass_done == 1) { //MK::store quality improvement of first iteration
			new_eps_impr = (new_qual - cur_qual) * static_cast<real_louvain>(Settings::GrainReconLouvainPrecision);
cout << "Setting dynamic precision to precision fraction, value " << static_cast<real_louvain>(Settings::GrainReconLouvainPrecision) << "\t\t" << new_eps_impr << endl;
		}

	//##MK::original stop criterion was with eps_impr by default 1.0e-6L
	//} while (nb_moves>0 && new_qual-cur_qual > eps_impr);
	} while( nb_moves > 0 && (new_qual - cur_qual) > new_eps_impr );

	return improvement;
}


//File: main_louvain.cpp
//-- community detection, sample main file
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


//#include "LOUVAIN_graph_binary.h"
//#include "LOUVAIN_louvain.h"



//#include "LOUVAIN_modularity.h"
/*
#include "zahn.h"
#include "owzad.h"
#include "goldberg.h"
#include "condora.h"
#include "devind.h"
#include "devuni.h"
#include "dp.h"
#include "shimalik.h"
#include "balmod.h"
*/

//using namespace std;


louvainHdl::louvainHdl()
{
	/*
	char *filename = NULL;
	char *filename_w = NULL;
	char *filename_part = NULL;
	*/
	precision = static_cast<long double>(Settings::GrainReconLouvainPrecision); //0.000001L;
	alpha = 0.5L;
	sum_se = 0.0L;
	sum_sq = 0.0L;
	max_w = 1.0L;

	q = NULL;

	wtype = WEIGHTED;
	nb_pass = 0;
	display_level = -1; //-2;
	kmin = 1;

	id_qual = 0;
	verbose = true;

	//myrandom(1234);
	srand(-1234);
}


louvainHdl::~louvainHdl()
{
	//if (q != NULL) {
	//	delete q; q = NULL;
	//}

	//do not delete owner is only backreference
}



/*
void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [-q id_qual] [-c alpha] [-k min] [-w weight_file] [-p part_file] [-e epsilon] [-l display_level] [-v] [-h]" << endl << endl;
  cerr << "input_file: file containing the graph to decompose in communities" << endl;

  cerr << "-q id\tthe quality function used to compute partition of the graph (modularity is chosen by default):" << endl << endl;

  cerr << "\tid = 0\t -> the classical Newman-Girvan criterion (also called \"Modularity\")" << endl;
  cerr << "\tid = 1\t -> the Zahn-Condorcet criterion" << endl;
  cerr << "\tid = 2\t -> the Owsinski-Zadrozny criterion (you should specify the value of the parameter with option -c)" << endl;
  cerr << "\tid = 3\t -> the Goldberg Density criterion" << endl;
  cerr << "\tid = 4\t -> the A-weighted Condorcet criterion" << endl;
  cerr << "\tid = 5\t -> the Deviation to Indetermination criterion" << endl;
  cerr << "\tid = 6\t -> the Deviation to Uniformity criterion" << endl;
  cerr << "\tid = 7\t -> the Profile Difference criterion" << endl;
  cerr << "\tid = 8\t -> the Shi-Malik criterion (you should specify the value of kappa_min with option -k)" << endl;
  cerr << "\tid = 9\t -> the Balanced Modularity criterion" << endl;

  cerr << endl;

  cerr << "-c al\tthe parameter for the Owsinski-Zadrozny quality function (between 0.0 and 1.0: 0.5 is chosen by default)" << endl;
  cerr << "-k min\tthe kappa_min value (for Shi-Malik quality function) (it must be > 0: 1 is chosen by default)" << endl;

  cerr << endl;

  cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)" << endl;
  cerr << "-p file\tstart the computation with a given partition instead of the trivial partition" << endl;
  cerr << "\tfile must contain lines \"node community\"" << endl;
  cerr << "-e eps\ta given pass stops when the quality is increased by less than epsilon" << endl;
  cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure" << endl;
  cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level" << endl;
  cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and quality" << endl;
  cerr << "-h\tshow this usage message" << endl;

  exit(0);
}

void set_args()
{
	type = WEIGHTED;
	//filename_w = NULL;
	id_qual = 0;
	alpha = 0.5L;
	kmin = 1;
	//filename_part = NULL;
	precision = 0.000001L;
	display_level = -2;
	verbose = true;
	//filename = NULL;
}

void parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
	type = WEIGHTED;
        filename_w = argv[i+1];
	i++;
	break;
      case 'q':
	id_qual = (unsigned short)atoi(argv[i+1]);
	i++;
	break;
      case 'c':
	alpha = atof(argv[i+1]);
	i++;
	break;
      case 'k':
	kmin = atoi(argv[i+1]);
	i++;
	break;
      case 'p':
        filename_part = argv[i+1];
	i++;
	break;
      case 'e':
	precision = atof(argv[i+1]);
	i++;
	break;
      case 'l':
	display_level = atoi(argv[i+1]);
	i++;
	break;
      case 'v':
	verbose = true;
	break;
      case 'h':
	usage(argv[0], "");
	break;
      default:
	usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
  if (filename == NULL)
    usage(argv[0], "No input file has been provided\n");
}

void
display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << ": " << ctime (&rawtime);
}
*/


void louvainHdl::init_quality(Graph *g, unsigned short nbc)
{
	if (nbc > 0)
		delete q;

	q = new Modularity(*g);

/*
	switch (id_qual) {
		case 0:
			q = new Modularity(*g);
			break;
		case 1:
			if (nbc == 0)
				max_w = g->max_weight();
			q = new Zahn(*g, max_w);
			break;
		case 2:
			if (nbc == 0)
				max_w = g->max_weight();
			if (alpha <= 0. || alpha >= 1.0)
				alpha = 0.5;
			q = new OwZad(*g, alpha, max_w);
			break;
		case 3:
			if (nbc == 0)
				max_w = g->max_weight();
			q = new Goldberg(*g, max_w);
			break;
		case 4:
			if (nbc == 0) {
				g->add_selfloops();
				sum_se = CondorA::graph_weighting(g);
			}
			q = new CondorA(*g, sum_se);
			break;
		case 5:
			q = new DevInd(*g);
			break;
		case 6:
			q = new DevUni(*g);
			break;
		case 7:
			if (nbc == 0) {
				max_w = g->max_weight();
				sum_sq = DP::graph_weighting(g);
			}
			q = new DP(*g, sum_sq, max_w);
			break;
		case 8:
			if (kmin < 1)
				kmin = 1;
			q = new ShiMalik(*g, kmin);
			break;
		case 9:
			if (nbc == 0)
				max_w = g->max_weight();
			q = new BalMod(*g, max_w);
			break;
		default:
			q = new Modularity(*g);
			break;
	}
 */
}


void louvainHdl::execute(vector<lvwtedge> const & these_edges, vector<unsigned int> & uip2comm )
{
	//use fixed seed fixed seed and MersenneTwister of louvainHdl class object instance instead of rand()
	//srand(-1234); //srand(time(NULL)+getpid());

	/*
	set_args();
	parse_args(argc, argv);

	time_t time_begin, time_end;
	time(&time_begin);
	*/

	//if (verbose)
	//	display_time("Begin");
	unsigned short nb_calls = 0;

	Graph g(these_edges);

	//Graph g(filename, filename_w, type);
	init_quality(&g, nb_calls);
	nb_calls++;

	if (verbose)
		cout << endl << "Computation of communities with the " << q->name << " quality function" << endl << endl;

	Louvain c(-1, precision, q);
	//if (filename_part!=NULL)
	//	c.init_partition(filename_part);

	bool improvement = true;

	real_louvain quality = (c.qual)->quality();
	real_louvain new_qual;

	int level = 0;

	vector<vector<int>*> ipassignment;

	do {
		if (verbose) {
			cout << "level " << level << ":\n";
			//display_time("  start computation");
			cout << "  network size: "
				<< (c.qual)->g.nb_nodes << " nodes, "
				<< (c.qual)->g.nb_links << " links, "
				<< (c.qual)->g.total_weight << " weight" << endl;
		}
		tictoc2.push_back( lvloglevel(level, (c.qual)->g.nb_nodes, (c.qual)->g.nb_links, (c.qual)->g.total_weight) );

		improvement = c.one_level( level, this->tictoc1 );

		new_qual = (c.qual)->quality();

		if (++level == display_level)
			(c.qual)->g.display();
		if (display_level == -1) {
			//c.display_partition();	//MK::here the assignment is written to cout data to be processed with hierarchy

			ipassignment.push_back(NULL);
			ipassignment.back() = c.pass_partition();
		}

		g = c.partition2graph_binary();
		init_quality(&g, nb_calls);
		nb_calls++;

		c = Louvain(-1, precision, q);

		if (verbose)
			cout << "  quality increased from " << quality << " to " << new_qual << endl;

		//quality = (c.qual)->quality();
		quality = new_qual;

		/*
		if (verbose)
			display_time("  end computation");
		if (filename_part!=NULL && level==1) //do at least one more computation if partition is provided
			improvement=true;
		*/
	} while (improvement);

	/*
	time(&time_end);
	if (verbose) {
		display_time("End");
		cerr << "Total duration: " << (time_end-time_begin) << " sec" << endl;
	}
	cerr << new_qual << endl;
	*/

	//traverse along hierarchical partitioning to find which community the initial ip ended up in
	//##MK::check with main_hierarchy
	unsigned int nip = static_cast<unsigned int>(ipassignment.at(0)->size());
	unsigned int nlevels = static_cast<unsigned int>(ipassignment.size());
//cout << "nip/nlevels " << nip << "\t\t" << nlevels << endl;
	for (unsigned int ip = 0; ip < nip; ip++) {
		unsigned int where = ip;
		unsigned int here = ipassignment.at(0)->at(where);
//cout << ip << "->" << here;
		for (unsigned int level = 1; level < nlevels; level++) {
			where = here;
			here = ipassignment.at(level)->at(where);
			//cout << "-->" << here;
		}
		//cout << endl;

		uip2comm.push_back( here );
//cout << ip << "--->" << here << endl;
	}
cout << "Grain segmentation completed" << endl;

	for (size_t i = 0; i < ipassignment.size(); i++) {
		if (ipassignment.at(i) != NULL) {
			delete ipassignment.at(i);
			ipassignment.at(i) = NULL;
		}
	}
cout << "Hierarchy structure deleted" << endl;

cout << "Profiling data transferred" << endl;

	if ( q != NULL ) {
		delete q;
	}
}


void louvainHdl::get_tctc( Louvain const & fromhere )
{

}


void louvainHdl::profiling()
{
	//##MK::further optimization and convenience tasks: bundle all in one file, incr ID and so forth

	//##MK::suboptimal... one file per rank
	string fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".LouvainProfilingIterations.csv";

	ofstream csvlog1;
	csvlog1.open(fn.c_str(), ofstream::out | ofstream::trunc);
	csvlog1 << setprecision(32);
	if (csvlog1.is_open() == true) {
		//header
		csvlog1 << "CurrentLevel;NumberOfMoves;CurrentQuality;NewQuality;NewMinusCurrentQual;SeqWallClock\n";
		csvlog1 << "1;1;1;1;1;s\n";
		csvlog1 << "CurrentLevel;NumberOfMoves;CurrentQuality;NewQuality;NewMinusCurrentQual;SeqWallClock\n";
		for (size_t i = 0; i < tictoc1.size(); ++i) {
			csvlog1 << tictoc1.at(i).levelid << ";" << tictoc1.at(i).mvcnt << ";";
			csvlog1 << tictoc1.at(i).cur_qual << ";" << tictoc1.at(i).new_qual << ";";
			csvlog1 << tictoc1.at(i).new_qual - tictoc1.at(i).cur_qual << ";" << tictoc1.at(i).dt << "\n";
		}
		csvlog1.flush();
		csvlog1.close();
	}

	fn = "DAMASKPDT.SimID." + to_string(Settings::SimID) + ".LouvainProfilingLevels.csv";
	ofstream csvlog2;
	csvlog2.open(fn.c_str(), ofstream::out | ofstream::trunc);
	csvlog2 << setprecision(32);
	if (csvlog2.is_open() == true) {
		//header
		csvlog2 << "Level;NumberOfLinks;NumberOfNodes;TotalWeight;SeqWallClock\n";
		csvlog2 << "1;1;1;1;s\n";
		csvlog2 << "Level;NumberOfLinks;NumberOfNodes;TotalWeight;SeqWallClock\n";
		for (size_t i = 0; i < tictoc2.size(); ++i) {
			csvlog2 << tictoc2.at(i).level << ";" << tictoc2.at(i).nblinks << ";" << tictoc2.at(i).nbnodes << ";" << tictoc2.at(i).tweight << ";n/a\n";
		}
		csvlog2.flush();
		csvlog2.close();
	}
}

/*
int main(int argc, char **argv)
{
	cerr << setprecision(18);

	louvainHdl* lvseqworker = NULL;
	lvseqworker = new louvainHdl; //MK::invokes the default constructor

	vector<lvwtedge> edges;

	//MK::BEGIN OF DEBUG build edges to test for successful code stripping and refactoring yet functioning of louvain original code
	ifstream finput;
	finput.open( "LOUVAIN_MOD_TEST.txt", fstream::in);
	while (!finput.eof()) {
		unsigned int ss, dd;
		real_louvain wt = 1.0L;
		finput >> ss >> dd >> wt;
		edges.push_back( lvwtedge( ss, dd, wt) );
//cout << edges.back().src << "\t\t" << edges.back().dest << "\t\t" << edges.back().wt << endl;
	}
	finput.close();
	//##MK::END OF DEBUG

	lvseqworker->execute( edges );

	if ( lvseqworker != NULL ) {
		delete lvseqworker;
	}
}
*/
