/*
  This header file defines tree decomposition related classes and
  functions and is used in conjunction with the Boost Graph Library
  and the graph_lib_boost.hpp header file produced by Aaron Adcock.

  In more detail, this library provides classes/structures involved
  with reading in (maybe in the future computing TD's, but for now it
  assumes that the tree decompostion is computed using the INDDGO
  software) tree decompositions, manipulating/using them, and will
  contain function headers which use the tree decompositions to
  produce features/models for machine learning algorithms.

  By Aaron Adcock, PhD Candidate, Stanford University
  2013
 */

#ifndef TREE_H_GUARD
#define TREE_H_GUARD

#include "graph_lib_boost.hpp"

/*
  
 */

class TreeDecomp 
{
private:
  Graph T;
  vector< vector<v_size_t> > bags;
  vector< vector<v_size_t> > vertex_locations;
  Graph * Gp;
  
  // Returning T is slow, though as a tree, the size is linear O(n)
  // where n is the number of vertices...work on this later...could
  // make it a pointer to a tree...
public:
  Graph get_tree() { return T; };
  int get_num_bags() { return num_vertices(T); };
  inline vector<size_t> get_bag_sizes();
  vector< vector<v_size_t> > get_bags() { return bags; };
  vector<v_size_t> get_bag(size_t bag) 
  { 
    if (bag < bags.size()) 
      return bags[bag]; 
    else
      {
	cout<<"Requested bag is out of range, returning empty set\n";
	vector<v_size_t> dummy;
	return dummy;
      }
  };
  vector< vector<v_size_t> > get_vertex_locations() { return vertex_locations; };
  Graph * get_graph_pointer() { return Gp; }
  inline TreeDecomp() { };
  inline TreeDecomp(Graph &, Graph &, vector< vector<v_size_t> > &, vector< vector< size_t> > &);
  vector<v_size_t> get_tree_eccentricity() { return tree_eccentricity(T); };
};

// Is this->Gp = &G correct?
TreeDecomp::TreeDecomp(Graph & T, Graph & G, vector< vector<v_size_t> > & bags, vector< vector<v_size_t> > & vertex_locations)
{
  this->T = T;
  this->Gp = &G;
  this->bags = bags;
  this->vertex_locations = vertex_locations;

  if (num_vertices(this->T) != this->bags.size())
    {
      cerr<<"The number of bags in the tree decomposition graph structure does not equal the number of bags in the bags vector.  Exiting. \n";
      exit(EXIT_FAILURE);
    }
}

vector<size_t> TreeDecomp::get_bag_sizes()
{
  size_t num_bags = bags.size();
  vector<size_t> bag_sizes;
  bag_sizes.reserve(num_bags);

  for (size_t i = 0; i < num_bags; ++i)
    bag_sizes.push_back((bags[i]).size());
    
  return bag_sizes;
}

/* Function Declarations */

// Loads a tree decomposition from file
TreeDecomp loadTreeDecomp(Graph & G, string filename);

// Creates a binary vector for each node in underlying graph that
// represents which bags each node is in
void createBagVectors(TreeDecomp TD, vector< vector<v_size_t> > & feature_vector, vector<size_t> & bags_to_keep);

// Writes a boolean feature vector to file
void writeFeatureVectors(string filename, vector< vector<v_size_t> > & feature_vector);

// Calculates minimum conductance communities by building them from tree decompositions along eccentricity lines
vector<double> td_eccentricity_communities(Graph& G, TreeDecomp& td, size_t max_community, bool large_comp_only);
vector<double> td_eccentricity_communities(Graph& G, TreeDecomp & td, size_t max_community, vector<double> conductance,  bool large_comp_only);
vector<double> td_eccentricity_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, vector< vector<Vert> > & best_communities,  bool large_comp_only);

// Calculates minimum conductance communities from bags of tree decompositions
vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, bool large_comp_only);
vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, bool large_comp_only);
vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, vector< vector<Vert> > & best_communities, bool large_comp_only);

// Calculate stats on a community set in a tree
void compute_community_tree_stats(vector< vector<Vert> > & communities, Graph & G, TreeDecomp & TD, vector<size_t> & num_bags, vector<size_t> & min_ecc, vector<size_t> & max_card, vector<size_t> & med_card, vector<size_t> & surface_area);
#endif
