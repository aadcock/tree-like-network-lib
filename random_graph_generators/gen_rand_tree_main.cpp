/*
This compiles into an executable that will generate a random tree
network.  By 'random tree', I mean that it will generate a b-ary tree,
and then add edges randomly (though randomly will be based on distance
in the tree).  If no randomness is desired, it will produce a basic tree.

By Aaron Adcock, PhD Candidate, Stanford University Jan. 2014
 */

#include "graph_lib_boost.hpp"

#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <queue>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace boost;

int get_tree_distance(v_size_t x, v_size_t y, int b);

int main(int argc, char *argv[])
{
  // Initialize variables, defaults are binary tree with no randomness
  string output_file = "tree";
  int b = 2;
  double avg_degree = 3;
  unsigned int num_nodes = 63;
  double locality_exponent = 0;
  bool noisy_tree = false;

  if (argc < 2)
    {
      cout<<"-o <output_file> \n-b <create_b-ary_tree> \n-n <number_of_nodes_in_tree> \n-noisy <flag present indicates adding edges probabilistically based on underlying tree.  No flag produces b-ary tree> \n-d <avg_degree of output tree> \n-locality <locality_exponent>\n";
      return(-1);
    }

  // Parse command line
  for(int i=1; i < argc; i++)
    {
      stringstream ss;
      if (strcmp(argv[i], "-o") == 0)
	output_file = argv[i+1];
      
      if (strcmp(argv[i], "-n") == 0)
	{
	  ss<<argv[i+1];
	  ss>>num_nodes;
	}

      if (strcmp(argv[i], "-b") == 0)
	{
	  ss<<argv[i+1];
	  ss>>b;
	}

      if (strcmp(argv[i], "-noisy") == 0)
	noisy_tree = true;

      if (strcmp(argv[i], "-d") == 0)
	{
	  ss<<argv[i+1];
	  ss>>avg_degree;
	}

      if (strcmp(argv[i], "-locality") == 0)
	{
	  ss<<argv[i+1];
	  ss>>locality_exponent;
	}
    }

  if(b <  1)
    {
      cout<<"b must be greater than or equal to 1.  Default of 2 is used. \n";
      b = 5;
    }

  if (num_nodes <= 0)
    {
      cout<<"The number of nodes must be positive. Default value of b^6 being used\n";
      num_nodes = pow(b, 6);
    }

  if (locality_exponent < 0)
    {
      cout<<"The locality exponent must be positive. Default value of 0 being used. \n";
      locality_exponent = 0;
    }

  cout<<"Output File Name: "<<output_file<<"\n";

  queue<int> next_edge;
  Graph G(num_nodes);

  if (noisy_tree)
    {
      // First construct distance array
      vector<double> cum_sum;
      vector< pair<v_size_t, v_size_t> > edges;
      cum_sum.reserve(num_nodes * num_nodes / 2 + 1);
      edges.reserve(num_nodes * num_nodes / 2 + 1);

      srand(time(NULL));
      // Build cdf vector of edge probabilities, might take too long.
      // We will see.
      for (size_t i = 0; i < num_nodes; ++i)
	for (size_t j = i+1; j < num_nodes; ++j)
	  {
	    pair<v_size_t, v_size_t> e(i, j);
	    edges.push_back(e);
	    cum_sum.push_back(cum_sum.back() + 1.0 / pow(get_tree_distance(i, j, b), locality_exponent));
	  }

      // Normalize
      for (size_t i = 0; i < cum_sum.size(); ++i)
	cum_sum[i] /= cum_sum.back();

      size_t num_edges = num_nodes * avg_degree;

      // Draw the correct number of edges from distribution decided above
      for (size_t i = 0; i < num_edges; ++i)
	{
	  // draw edge from distribution.  Note, rand() is not the
	  // best way to do this, may need to switch to a better
	  // random number generator
	  double r = double(rand()) / RAND_MAX;
	  bool not_found = true;
	  int edge_index = 0;
	  while (not_found)
	    {
	      if (cum_sum[edge_index] > r)
		{
		  Vert v1 = vertex(edges[edge_index].first, G);
		  Vert v2 = vertex(edges[edge_index].second, G);
		  bool edge_exists = (edge(v1, v2, G)).second;

		  if (edge_exists)
		    {
		      edge_index = 0;
		      r = double(rand()) / RAND_MAX;
		    }
		  else
		    {
		      add_edge(v1, v2, G);
		      add_edge(v2, v1, G);
		      not_found = false;
		    }
		}
	      else
		edge_index++;
	    }
	}
    }
  else
    {
      for (int i = 0; i < b; ++i)
	next_edge.push(0);
    
      for (int curr_vert = 1; curr_vert < num_nodes; ++curr_vert)
	{
	  // Get vertices for new edge
	  int prev_vert = next_edge.front();
	  Vert v1 = vertex(curr_vert, G);
	  Vert v2 = vertex(prev_vert, G);

	  // Must add edge in both directions as G is a directed graph
	  // representation
	  add_edge(v1, v2, G);
	  add_edge(v2, v1, G);
      
	  // Update edge queue
	  next_edge.pop();

	  for (int i = 0; i < b; ++i)
	    next_edge.push(curr_vert);
	}
    }

  // Need to write out graph here.
  write_dimacs(G, output_file);
}

/* 
   This function takes an integers x, y which correspond to a b-ary
   tree with node 0 at the root.  It then computes the distance
   between x and y.

   This function thinks of trees like this
                                          Order
                      0          Level 0: 0:0
                    /   \
                  1       2      Level 1: 1:0, 2:1
                /   \   /   \
               3     4 5     6   Level 2  3:0, 4:1, 5:2, 6:3
 */

int get_tree_distance(v_size_t x, v_size_t y, int b)
{
  int distance = 0;

  // if equal, distance is zero.
  if (x == y)
    return distance;

  // First need tree level of both x and y
  int tree_level_x = floor(log(x + 1) / log(b));
  int tree_level_y = floor(log(y + 1) / log(b));

  int level_diff = abs(tree_level_x - tree_level_y);

  // Now need tree order
  int order_x = x + 1 - pow(b, tree_level_x);
  int order_y = y + 1 - pow(b, tree_level_y);

  // Move the leafier node up to the rootier node's level
  if (tree_level_x > tree_level_y)
    order_x = floor(double(order_x) * pow(b, tree_level_y) / pow(b, tree_level_x));
  else if (tree_level_x < tree_level_y)
    order_y = floor(double(order_y) * pow(b, tree_level_x) / pow(b, tree_level_y));

  // Calculate distance
  distance = level_diff + 2.0 * floor(log(abs(order_x - order_y) + 1) / log(b));

  return distance;
}
