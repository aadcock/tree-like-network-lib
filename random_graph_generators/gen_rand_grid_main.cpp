/*
This compiles into an executable that will generate a random Euclidean
grid network.  By 'random grid', I mean that it will generate a b
dimensional grid, with NxNx...xN sides, and then add edges randomly
(though randomly will be based on distance in the grid).  If no
randomness is desired, it will produce a basic grid.

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

int get_grid_distance(v_size_t x, v_size_t y, int b, v_size_t num_nodes_side);

int main(int argc, char *argv[])
{
  // Initialize variables, defaults are binary tree with no randomness
  string output_file = "grid";
  int b = 2;
  double avg_degree = 3;
  unsigned int num_nodes_side = 10;  // Number of nodes on a side
  double locality_exponent = 0;
  bool noisy_grid = false;

  if (argc < 2)
    {
      cout<<"-o <output_file> \n-b <create_b_dimensional_grid> \n-n <number_of_nodes_on_side> \n-noisy <flag present indicates adding edges probabilistically based on underlying grid.  No flag produces b dimensional grid> \n-d <avg_degree of output noisy grid, does nothing if noisy flag is not used> \n-locality <locality_exponent>\n";
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
	  ss>>num_nodes_side;
	}

      if (strcmp(argv[i], "-b") == 0)
	{
	  ss<<argv[i+1];
	  ss>>b;
	}

      if (strcmp(argv[i], "-noisy") == 0)
	noisy_grid = true;

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
      b = 2;
    }

  if (num_nodes_side <= 0)
    {
      cout<<"The number of nodes must be positive. Default value of 10 being used\n";
      num_nodes_side = 10;
    }

  if (locality_exponent < 0)
    {
      cout<<"The locality exponent must be positive. Default value of 0 being used. \n";
      locality_exponent = 0;
    }

  cout<<"Output File Name: "<<output_file<<"\n";


  // Calculate total number of nodes in graph
  v_size_t total_nodes = 1;
  for (int i = 0; i < b; ++i)
    total_nodes *= num_nodes_side;

  Graph G(total_nodes);

  if (noisy_grid)
    {
      // First construct distance array
      vector<double> cum_sum;
      vector< pair<v_size_t, v_size_t> > edges;
      cum_sum.reserve(total_nodes * total_nodes / 2 + 1);
      edges.reserve(total_nodes * total_nodes / 2 + 1);

      srand(time(NULL));
      // Build cdf vector of edge probabilities, might take too long.
      // We will see.
      for (size_t i = 0; i < total_nodes; ++i)
	for (size_t j = i+1; j < total_nodes; ++j)
	  {
	    pair<v_size_t, v_size_t> e(i, j);
	    edges.push_back(e);
	    cum_sum.push_back(cum_sum.back() + 1.0 / pow(get_grid_distance(i, j, b, num_nodes_side), locality_exponent));
	  }

      // Normalize
      for (size_t i = 0; i < cum_sum.size(); ++i)
	cum_sum[i] /= cum_sum.back();

      size_t num_edges = total_nodes * avg_degree;

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
      for (v_size_t curr_vert = 0; curr_vert < total_nodes; ++curr_vert)
	{
	  v_size_t next_node = 1;
	  Vert v1 = vertex(curr_vert, G);
	  while (next_node < total_nodes)
	    {
	      // If statement makes sure not to add edges to the last
	      // nodes in the grid
	      if ((curr_vert + next_node) % (next_node * num_nodes_side) >= next_node)
		{
		  Vert v2 = vertex(curr_vert + next_node, G);
		  add_edge(v1, v2, G);
		  add_edge(v2, v1, G);
		}
	      next_node *= num_nodes_side;
	    }
	}
    }

  // Need to write out graph here.
  write_dimacs(G, output_file);
}

/* 
   This function takes an integers x, y which correspond to a
   b-dimensional grid.  It then computes the distance between x and y.

   This function thinks of grids numbered like this (it gets harder to
   imagine in higher dimensions, but the numbering just repeats for
   each layer.

              0 - 1 - 2 - 3
	      |   |   |   |
	      4 - 5 - 6 - 7
	      |   |   |   |
	      8 - 9 - 10- 11
	      |   |   |   |
	      12- 13- 14- 15
 */

int get_grid_distance(v_size_t x, v_size_t y, int b, v_size_t num_nodes_side)
{
  int distance = 0;

  // if equal, distance is zero.
  if (x == y)
    return distance;

  // First, translate each node into coordinates.  Node zero is at (0,0, ... ,0)
  vector< int > coordinates_x;
  vector< int > coordinates_y;

  v_size_t next_dimension = pow(num_nodes_side, b - 1);
  for (int i = 0; i < b; ++i)
    {
      coordinates_x.push_back(x / next_dimension);
      x = x % next_dimension;
      coordinates_y.push_back(y / next_dimension);
      y = y % next_dimension;
      next_dimension /= num_nodes_side;
    }

  for (int i = 0; i < b; ++i)
    distance += abs(coordinates_x[i] - coordinates_y[i]);

  return distance;
}
