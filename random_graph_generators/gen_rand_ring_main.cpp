/*
This compiles into an executable that will generate a random ring
network.  By 'random ring', I mean that it will generate a ring,
and then add edges randomly (though randomly will be based on distance
in the ring).  If no randomness is desired, it will produce a basic ring.

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

int get_ring_distance(v_size_t x, v_size_t y, int num_nodes);

int main(int argc, char *argv[])
{
  // Initialize variables, defaults are binary ring with no randomness
  string output_file = "ring";
  double avg_degree = 3;
  unsigned int num_nodes = 63;
  double locality_exponent = 0;
  bool noisy_ring = false;

  if (argc < 2)
    {
      cout<<"-o <output_file> \n-n <number_of_nodes_in_ring> \n-noisy <flag present indicates adding edges probabilistically based on underlying ring.  No flag produces basic ring> \n-d <avg_degree of output ring> \n-locality <locality_exponent>\n";
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

      if (strcmp(argv[i], "-noisy") == 0)
	noisy_ring = true;

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

  if (num_nodes <= 0)
    {
      cout<<"The number of nodes must be positive. Default value of b^6 being used\n";
      num_nodes = 100;
    }

  if (locality_exponent < 0)
    {
      cout<<"The locality exponent must be positive. Default value of 0 being used. \n";
      locality_exponent = 0;
    }

  cout<<"Output File Name: "<<output_file<<"\n";

  Graph G(num_nodes);

  if (noisy_ring)
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
	    cum_sum.push_back(cum_sum.back() + 1.0 / pow(get_ring_distance(i, j, num_nodes), locality_exponent));
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
      for (int curr_vert = 0; curr_vert < num_nodes; ++curr_vert)
	{
	  Vert v1 = vertex(curr_vert, G);
	  Vert v2;
	  if (curr_vert == num_nodes-1)
	    v2 = vertex(0, G);
	  else
	    v2 = vertex(curr_vert + 1, G);
	  
	  // Must add edge in both directions as G is a directed graph
	  // representation
	  bool edge_exists = (edge(v1, v2, G)).second;
	  
	  if (!edge_exists)
	    add_edge(v1, v2, G);

	  edge_exists = (edge(v2, v1, G)).second;
	  if (!edge_exists)
	    add_edge(v2, v1, G);
	}
    }

  // Need to write out graph here.
  write_dimacs(G, output_file);
}

/* 
   This function takes an integers x, y which correspond to a ring
   with the node indices ordered sequentially.  
 */

int get_ring_distance(v_size_t x, v_size_t y, int num_nodes)
{
  int distance = 0;

  // if equal, distance is zero.
  if (x == y)
    return distance;

  // Distance on ring 
  if (abs(x - y) <= num_nodes / 2)
    distance = abs(x-y);
  else
    distance = num_nodes / 2 - ((int) abs(x - y) % (num_nodes / 2));

  return distance;
}
