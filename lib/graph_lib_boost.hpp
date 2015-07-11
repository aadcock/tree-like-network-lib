/*
This header file defines some graph related classes and functions and
is used in conjunction with the Boost Graph Library

By Aaron Adcock, PhD Candidate, Stanford University
2011-2013
*/

#ifndef GRAPH_H_GUARD
#define GRAPH_H_GUARD

#define CHUNKSIZE 512

#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
// #include "aaron/boost/graph/property_maps.hpp"

using namespace boost;
using namespace std;

/****************************************
 * Structures & Classes 
 ****************************************/

/* 
 Used to sort an index vector based on comparisons of the personalized
 Page Rank vector's indices.  This allows the algorithm to keep track
 of which nodes are kept during the sweep cut.x  
*/

template <class T> struct compareIndgtl
{
  compareIndgtl(const T arr):
    arr(arr){}
  
  bool operator() (const size_t a, const size_t b) const
  {
    return arr[a] > arr[b];
  }

  const T arr;
};

template <class T> struct compareIndltg
{
  compareIndltg(const T arr):
    arr(arr){}
  
  bool operator() (const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }

  const T arr;
}
 ;

/* 
 Used to keep track of the k_core of the vertex (color) for
 visualization algorithms
*/

class TreeDecomp;

struct bfs_node
{
  size_t parent;
  size_t level;

  bfs_node() {};
  bfs_node(size_t a, size_t b) {parent = a; level = b;};
};

struct vert_prop
{
  int core_num;
  double color;
  vector<size_t> prevIndices;
  vector<unsigned long> distances;
  vector<vector<unsigned long> > predecessors;
  vector<float> distribution;
};

struct edge_prop
{
  // double weight;
  // edge_prop() { weight = 1.0; };
};

struct graph_prop
{
  bool core_calc, tree_calc;
  TreeDecomp * Tp;

  graph_prop() { core_calc = false; tree_calc = false; };
};

/* Type Definitions */
typedef adjacency_list< vecS, vecS, directedS, vert_prop, edge_prop, graph_prop > Graph;
typedef adjacency_list< vecS, vecS, undirectedS, vert_prop, edge_prop, graph_prop > UGraph;

typedef graph_traits<Graph>::vertex_descriptor Vert;
typedef graph_traits<UGraph>::vertex_descriptor UVert;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<UGraph>::edge_descriptor UEdge;
typedef pair<int,int> intPair;
typedef graph_traits<Graph>::vertices_size_type v_size_t;
typedef graph_traits<UGraph>::vertices_size_type v_Usize_t;
typedef graph_traits<Graph>::edges_size_type e_size_t;
typedef graph_traits<UGraph>::edges_size_type e_Usize_t;

// Used with boost graphViz function to generate visualization

template<class vertProp>
class vertex_writer 
{
public:
  vertex_writer(vertProp _vprop) : vprop(_vprop) {}
  void operator()(ostream& out, const Vert& v) const
  {
    vert_prop vp       = vprop[v];
    unsigned long pind = vp.prevIndices[vp.prevIndices.size()-1];
      
    out<<"[color="<<vp.color<<"]";// <<" label=\""<<pind<<"\""<<"]";
  }
private:
  vertProp vprop;
};

// struct vertex_writer
// {
  //  void operator()(ostream& out, const Vert& v) const
  // {
  //   out<<"[fillcolor="<<v.color<<" label=\"c "<<v.core_num<<"\"]";
//   }
// }

struct edge_writer
{
  void operator()(ostream& out, const Edge& e) const
  {
    out<<"";
  }
};

struct graph_writer
{
  void operator()(ostream& out) const
  {
    out<<"node [colorscheme=\"rdylgn11\" shape=\"circle\"]\n";
  }
};

// Loads a color file
template <class T>
void load_color_file(string & color_file, vector<T> & color_vec)
{
  ifstream color;
  color.open(color_file.c_str());

  if(color.is_open())
    {
      while(color.good())
	{
	  string line;
	  bool flag = true;

	  getline(color,line);
	  if(line.substr(0,1)=="#")
	    flag = false;

	  if(flag)
	    {
	      istringstream line_stream(line);

	      while(line_stream.good())
		{
		  T val = 0;
		  double dummy;
		  if(!(line_stream>>dummy))
		    cerr<<"Color file contains non-numeric data.\n";
		  else
		    {
		      val = (T) dummy;
		      color_vec.push_back(val);		     
		    }
		}
	    }
	}
    }
  else
    {
      cerr<<"Failed to open color file\n";
      exit(EXIT_FAILURE);
    }
}

/* Function Headers */

// loadGraph functions take a graph defined from a vertex pairs
// representing edges These lists can come from a file or from an STL
// list container of intPair
Graph loadGraph(string filename, string delimiter);
Graph loadGraph(string filename);
Graph loadGraph(vector<intPair> edgeList);
Graph loadGraph(string filename, string delimiter,const bool directed);

// loads a graph in dimacs format
Graph loadDimacsGraph(string filename);
Graph loadDimacsGraph(string filename, const bool directed);

// Contained within the same source file as loadGraph, these are
// functions for searching through a vector of sorted integers
int binarySearch(vector<int>& sortedVector, int key);
int binarySearch(vector<v_size_t>& sortedVector, v_size_t key);

// in file approx_Page_Rank, this function approximates the
// personalized page rank vector using vector s as the seed vector.  s
// generally sums to 1
void approxPR(Graph& G, vector<double>& r, const double& alpha, const double& eps,vector<double>& p, const int ind);

// This function computes and returns the step by step probabilities
// associated with a random walk of num_steps and initial walker
// distribution initial distribution.
vector< vector<double> > random_walk(Graph& G, vector<double>& initial_distribution, const int num_steps);

// Calculates the ncp plots using the best community found by
// calculating page rank vectors (using approxPR) and calculating
// sweep cuts The function then keeps the cut with the best
// conductance at each size scale
vector<double> ncp_calc(Graph& G, const v_size_t maxC, const int step_size, const double alpha); // , const float eps);

vector<double> ncp_calc(Graph& G, const v_size_t maxC, const int step_size, const double alpha, bool display_output);

vector<double> ncp_calc(Graph& G, const v_size_t max_community, const int step_size, const double alpha, vector< vector<Vert> > & best_communities, bool display_output);

// vector<double> ncp_calc_snap(Graph& G, const int maxC, const int
// iter, const double alpha, const double eps);

// Calculates the ncp plot for each k_core
void k_core_ncp(Graph& G, const string prefix);
void k_core_ncp(Graph& G, const string prefix,vector<int>& core,const double alpha,const int step);

// Takes G and outputs the largest connected component
Graph connected(Graph& G,int s);
Graph connected(Graph& G);

// Takes G and labels nodes by component.  Note, the component labels
// are unordered
vector<v_size_t> components(Graph& G);

// Just provides a list of nodes in giant component
vector<Vert> list_giant_component(Graph& G);

// Takes G and outputs the subgraph induced on the subset S of G
Graph subset(Graph& G, vector<Vert>& S);
//Graph subset(Graph& G, vector<v_size_t>& Sindex);
UGraph undirected_subset(Graph& G, vector<Vert>& S);

// Converts an index vector to a vertex vector
vector<Vert> convert_indices_to_vertices(Graph & G, vector<v_size_t> & index);
// Produces k_core decomposition
// void k_core(Graph& G);
vector<int> k_core(Graph& G);

// Produces statistics on graph, using k_core decomposition
void k_core_stat(vector<vector<float> > & core_stat,Graph& G);
void k_core_stat(vector<vector<float> >& core_stat, Graph & G, const string output_file_prefix);


// Produces snapshots of graph using integers labels on nodes
void labeled_snapshot(Graph& G, vector<int>& node_labels, int label, string filename, const double alpha, const double eps);

// Produces distribution of labels in nodes touched by personalized
// PageRank diffusion
vector<double> community_snapshot(Graph& G, vector<int>& node_labels, int label, const double alpha, const double eps);

map<int, double> community_snapshot(Graph& G, map<v_size_t, int>& node_labels, int label, const double alpha, const double eps);

// Produces snapshots of the graph using a seed node and the approxPR
// function to generate neighborhoods
void spectral_snapshot(Graph& G, vector<Vert>& snapshot, vector<double>& page_rank, Vert seed_node, const double alpha, const double eps);

// Produces plots of k_core statistics
void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string> label, string title, string xlabel, string ylabel, string xtics); 

void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string> label, string title, string xlabel, string ylabel, string xtics,bool suppressOutput); 

void produce_plot(string infile, string outfilename,int x_col, vector<int> columns_to_plot, vector<string> label, string title, string xlabel, string ylabel, string xtics,bool suppressOutput); 

void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xstart,string xtic,bool suppressOutput);

void produce_3d_plot(string infile, string outfilename, int z_col, vector<string>label, string title, string xlabel, string ylabel, string zlabel, string xtics, string ytics);

void produce_3d_plot(string infile, string outfilename, int x_col, int y_col, int z_col, vector<string>label, string title, string xlabel, string ylabel, string zlabel, string xtics, string ytics, bool suppressOutput);

void produce_loglog_plot(string infile, string outfilename, string outplotname, vector<int> columns_to_plot, vector<string> label, string title, string xlabel, string ylabel, string xtics);

void produce_loglog_plot(string infile, string outfilename,string outplotname, vector<int> columns_to_plot, vector<string> label, string title, string xlabel, string ylabel, string xtics, bool suppressOutput);

// Calculates and stores distances/predecessors
void BFS_distance(Graph& G, string outputFileName);
void BFS_distance(Graph& G);

// Calculates and stores gromov hyperbolicity
vector< vector<double> > calc_gromov(Graph& G, v_size_t diam);

// Calculates gromov hyperbolicity AND quads_to_keep quadrupeds from
// each delta value
vector< vector<double> > calc_gromov_quads(Graph& G, v_size_t diam, int quads_to_keep, int dist_keep, string quad_output,bool single);

// Calculates delta slim triangles
vector<vector<double> > calc_delta_slim(Graph& G, v_size_t diam);

// Calculates delta fat triangles (another definition of hyperbolicity)
vector<vector<double> > calc_delta_fat(Graph& G, v_size_t diam);

// Breadth-first search of graph from source node s to all pairs,
// keeps distances and predecessors on (all) shortest paths or finds
// all vertices in connected component
void BFS_vertices_found(const Vert& s, const Graph& G, vector<bool>& searched, const v_size_t & current_comp, vector<v_size_t> & component, v_size_t & comp_size);

void BFS_source_all(const Vert& s,const Graph& G,vector<v_size_t>& distances, vector<vector<Vert> >& predecessors);
void BFS_source_all(const Vert& s,const Graph& G,vector<v_size_t>& distances);
void BFS_source_all(const Vert& s,const Graph& G,vector<v_size_t>& distances, vector<v_size_t>& num_paths);
void BFS_source_all(const Vert& s, const Graph& G, vector<bfs_node> & bfs_tree, vector<size_t> & tree_label, size_t & height);

// Breadth-first search of graph from source node s to target node t.
// Returns all vertices on all shortest paths.
void BFS_source_target(const Vert& s, const Vert& t, Graph& G, vector<Vert>& path, vector< vector<Vert> >& predecessors);

// BFS of graph from source node s to target node to, but returns only
// a single shortest path.
void BFS_source_target_single(const Vert& s, const Vert& t, Graph& G, vector<Vert>& path, vector< vector<Vert> >& predecessors);

// Breadth-first search of graph that terminates after two targets are
// found, returns all vertices on all shortest paths to both.
void BFS_two_target(const Vert& s, const Vert& t1, const Vert& t2, Graph& G, vector<Vert>& path1, vector<Vert>& path2, vector< vector<Vert> >& predecessors);

// Breadth-first search of graph that terminates after three targets
// are found, returns all vertices on all shortest paths to all three.
void BFS_three_target(const Vert& s, const Vert& t1, const Vert& t2, const Vert& t3, Graph& G, vector<Vert>& path1, vector<Vert>& path2, vector<Vert>& path3, vector< vector<Vert> >& predecessors);

// Breadth first search to determine the eccentricity of source node
int BFS_eccentricity(const Vert& s, const Graph& G);

// Fast algorithm for calculating tree eccentricity of all nodes in
// linear time.
vector<v_size_t> tree_eccentricity(Graph& T);

// forms a vector containing all vertices on ANY shortest path between
// the source s and the target t.
void form_path(const Vert& s, const Vert& t, const Graph& G, const vector< vector<Vert> >& predecessors, vector<Vert>& path, queue<Vert>& Q, vector<bool>& onPath,property_map<Graph,vertex_index_t>::type& index);

// writes a dimacs format version of the graph to a text (.dimacs)
// file
void write_dimacs(Graph& G, string output_file_name);

// writes out a color file for use with Blair Sullivan's tree
// decomp. code
void write_color_file(string output_file, string input_file, vector<double> color);

void write_color_file(string output_file, string input_file, vector<int> color);

// Writes a graphViz .dot file
void write_graphviz(Graph & G, string output_file, vector<double> color, vector<int> square_nodes, bool directed);

// Writes a size scaled graphViz .dot file
void write_scaled_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, bool directed);

// With square nodes
void write_scaled_square_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, vector<int> square_nodes, bool directed);

// Writes a size scaled, labeled, graphViz
void write_scaled_labeled_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, bool directed);

Graph get_k_shells(vector<v_size_t> & shells, Graph & G);
Graph get_k_shells(vector<v_size_t> & shells, vector<v_size_t> & core, Graph & G);

Graph collapse_core(v_size_t k_core, Graph & G);

// Machine learning routines
bool label_propagation(Graph & G, vector<Vert> & labeled_vertices, vector<int> & known_labels, vector<int> & predicted_labels, size_t num_iter);
double label_accuracy(vector<int> & predicted_labels, vector<int> & labels);

/* Distance oracle types */

// Basic distance oracle object
class distanceOracle
{
public:
  distanceOracle() { }
  
  virtual size_t query(size_t i, size_t j) = 0;
};

// Distance oracle based on a BFS spanning tree, uses square root LCA
// algorithm. 
class rootBFSDistanceOracle: public distanceOracle 
{
private:
  vector<size_t> tree_label;
  vector<bfs_node> bfs_tree;
  vector<size_t> section_ancestor;
  size_t tree_height;

  // Note, this requires the nodes to be sorted by level in bfs_tree,
  // with bfs_tree[0].level = 0;  This function finds the ancestor of
  // each node in the tree, which is in the previous section.
  void preprocess_tree()
  {
    size_t levels_per_section = (size_t) sqrt((double) tree_height);

    section_ancestor.assign(bfs_tree.size(), 0);
    for (size_t i = 0; i < bfs_tree.size(); ++i)
      {
	bfs_node current = bfs_tree[i];

	// Everywhere else arrays start with zero, so same here.
	if (current.level < levels_per_section)
	  section_ancestor[i] = 0;  
	else if (current.level % levels_per_section == 0)
	  section_ancestor[i] = current.parent;
	else
	  section_ancestor[i] = section_ancestor[current.parent];
      }
  }

public:
  rootBFSDistanceOracle(Vert& s, Graph& G):distanceOracle()
  {
    BFS_source_all(s, G, bfs_tree, tree_label, tree_height);
    preprocess_tree();
  }
  rootBFSDistanceOracle(vector<bfs_node> & input_bfs_tree, vector<size_t> input_tree_label, size_t height):distanceOracle()
  {
    bfs_tree = input_bfs_tree;
    tree_label = input_tree_label;
    tree_height = height;

    preprocess_tree();
  }

  void set_tree(Vert s, Graph & G)
  {
    BFS_source_all(s, G, bfs_tree, tree_label, tree_height);
    preprocess_tree();
  }

  void set_tree(vector<bfs_node> & input_bfs_tree, vector<size_t> input_tree_label, size_t height)
  {
    bfs_tree = input_bfs_tree;
    tree_label = input_tree_label;
    tree_height = height;

    preprocess_tree();
  }

  size_t query(size_t i, size_t j)
  {
    size_t tree_label_i = tree_label[i];
    size_t tree_label_j = tree_label[j];

    size_t x_label = tree_label_i;
    size_t y_label = tree_label_j;
    
    // First find common section
    while (section_ancestor[x_label] != section_ancestor[y_label])
      {
	bfs_node x = bfs_tree[x_label];
	bfs_node y = bfs_tree[y_label];

	// Move node further from root towards root
	if (x.level > y.level)
	  x_label = section_ancestor[x_label];
	else
	  y_label = section_ancestor[y_label];
      }

    // Now, since having the same section ancestor means that both tree
    // nodes are in the next section from the ancestor (i.e. x and y
    // are in the same section and they, at WORST, have the section
    // ancestor as the least common ancestor, we can just step until
    // they match!
    while (x_label != y_label)
      {
	bfs_node x = bfs_tree[x_label];
	bfs_node y = bfs_tree[y_label];

	if (x.level > y.level)
	  x_label = x.parent;
	else
	  y_label = y.parent;
      }
    
    // Finally, calculate tree distance
    bfs_node tree_node_i = bfs_tree[tree_label_i];
    bfs_node tree_node_j = bfs_tree[tree_label_j];
    bfs_node tree_node_lca = bfs_tree[y_label];

    size_t bfs_distance = tree_node_i.level + tree_node_j.level - 2 * tree_node_lca.level;

    return bfs_distance;
  }
};

#endif
