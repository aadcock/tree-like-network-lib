/*
This compiles into an executable that will generate a random line
graph in one of the following ways: 

1. Planted Community line.  This is similar to the planted community
model (intra-cluster edge probability > extra-cluster edge probability
for k communities) except that the clusters are placed on a line and
are connected with adjacent clusters with some probability.

2. Layered line.  The idea of this is to make a simple model of the
core-periphery networks.  There will be a dense set of nodes in the
center of the line.  Then one each side, there will be a set of nodes
which are slightly less dense and so on.  Thus we get a kind of
gradation out to the periphery of the line where the sparsest node
sets are.  Each layer will have same number of nodes and an expected
number of edges to be generated when that node is added to the
network. Edges will only occur between nodes and adjacent layers.

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


int main(int argc, char *argv[])
{
  // Initialize variables, defaults are binary tree with no randomness
  string output_file = "tree";
  int k = 2;
  unsigned int num_nodes = 63;
  double cluster_density = 0.8;
  double line_density = 0.2;
  bool layered_line = false;
  srand(time(NULL));
  if (argc < 2)
    {
      cout<<"-o <output_file> \n-layered <Uses the layered line model instead> \n-n <number_of_nodes> \n-k <number_of_clusters> \n-cluster_density <probabilty of adding edge inside community> \n-line_density <probability of adding an edge to an adjacent cluster>\n";
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

      if (strcmp(argv[i], "-k") == 0)
	{
	  ss<<argv[i+1];
	  ss>>k;
	}

      if (strcmp(argv[i], "-layered") == 0)
	layered_line = true;

      if (strcmp(argv[i], "-cluster_density") == 0)
	{
	  ss<<argv[i+1];
	  ss>>cluster_density;
	}

      if (strcmp(argv[i], "-line_density") == 0)
	{
	  ss<<argv[i+1];
	  ss>>line_density;
	}
    }

  if(k <  1)
    {
      cout<<"k must be greater than or equal to 1.  Default of 2 is used. \n";
      k = 2;
    }

  if (num_nodes <= 0)
    {
      cout<<"The number of nodes must be positive. Default value of b^6 being used\n";
      num_nodes = 100;
    }

  if (cluster_density < 0.0 || cluster_density > 1.0)
    {
      cout<<"The community density must be between 0 and 1. Default value of 0.8 being used. \n";
      cluster_density = 0.8;
    }

  if (line_density < 0.0 || line_density > 1.0)
    {
      cout<<"The line density must be between 0 and 1. Default value of 0.2 being used. \n";
      cluster_density = 0.2;
    }

  if (layered_line)
    {
      cout<<"Not implemented yet, sorry. \n";
      exit(EXIT_FAILURE);
    }

  cout<<"Output File Name: "<<output_file<<"\n";

  Graph G(num_nodes);

  // Place each node in a community
  int community_size = num_nodes / k;
  vector<int> community_membership;

  int community = 0;
  for (int i = 0; i < num_nodes; ++i)
    {
      community_membership.push_back(community);

      if ((i + 1) % community_size == 0)
	community = min(k-1, community+1);
    }

  // This is slow, but for now just look at every edge and add it with
  // the appropriate probability.
  if (!layered_line)
    {
      for (int i = 0; i < num_nodes; ++i)
	for (int j = i+1; j < num_nodes; ++j)
	  {
	    if (community_membership[i] == community_membership[j])
	      {
		double r = double(rand()) / RAND_MAX;
		if (r < cluster_density)
		  {
		    Vert v1 = vertex(i, G);
		    Vert v2 = vertex(j, G);
		    add_edge(v1, v2, G);
		    add_edge(v2, v1, G);
		  }
	      }
	    else if (abs(community_membership[i] - community_membership[j]) == 1)
	      {
	        double r = double(rand()) / RAND_MAX;
		if (r < line_density)
		  {
		    Vert v1 = vertex(i, G);
		    Vert v2 = vertex(j, G);
		    add_edge(v1, v2, G);
		    add_edge(v2, v1, G);
		  }
	      }
	  }
    }

  write_dimacs(G, output_file);
}
