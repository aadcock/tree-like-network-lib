/*
  This file defines the oracle class which serves as an approximate
  distance oracle for large networks.  The basic ideas are from the
  paper: "A Geometric Distance Oracle for Large Real-World Graphs", by
  Ajwani et al (on the ArXiv at the writing of this code).  We want to
  improve on their method of picking the root node by using the
  $k$-core decomposition.

  By Aaron Adcock, Ph.D. Candidate, 
  Stanford University 2014
 */

#include "../lib/graph_lib_boost.hpp"
#include <cmath>

int main(int argc, char *argv[])
{
  string input_file = "";
  bool use_degree = false;
  //User input

  if(argc < 2)
    {
      cout<<"-i <input_file> \n-deg <use degree instead of k-core for selecting root \n";
      return(-1);
    }

  for(int i = 1; i < argc; i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i") == 0)
	input_file = argv[i + 1];
      else if (strcmp(argv[i], "-deg") == 0)
	use_degree = true;
    }

  if (input_file.length() == 0)
    {
      cerr<<"Must provide an input file! \n";
      exit(EXIT_FAILURE);
    }
  
  // Get input
  string delimiter = "\t";
  srand(time(NULL));
  bool directed = false;
  Graph G = loadGraph(input_file, delimiter, directed);
  G = connected(G);
  size_t size = num_vertices(G);

  vector<int> ranking(size, 0);

  if (use_degree)
    {
      for (size_t i = 0; i < ranking.size(); ++i)
	ranking[i] = out_degree(vertex(i, G), G);
    }
  else
    ranking = k_core(G);
  
  vector<int> vert_sorted(ranking.size(), 0);
  for (size_t i = 0; i < vert_sorted.size(); ++i)
    vert_sorted[i] = i;

  sort(vert_sorted.begin(), vert_sorted.end(), compareIndgtl<vector<int> >(ranking));

  cout<<"Graph loaded.\n";
  
  Vert root = vertex(vert_sorted[0], G);

  // Create oracle
  rootBFSDistanceOracle oracle(root, G);
  
  BFS_distance(G);

  // Compare the oracle to actual distances
  graph_traits<Graph>::vertex_iterator vi, vie;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  vector<double> absolute_error;
  vector<size_t> num_pairs;
  
  for (tie(vi, vie) = vertices(G); vi != vie; ++vi)
    {
      Vert curr = *vi;
      size_t j  = index[curr];
      vector<unsigned long> distances = G[curr].distances;
      for (size_t i = 0; i < distances.size(); ++i)
	{
	  int oracle_dist = (int) oracle.query(j, i);
	  size_t real_dist = (size_t) distances[i];
	  
	  double abs_error;
	  abs_error = abs(real_dist - oracle_dist);

	  // Note, resizing vector only works because the diameter
	  // should be small relative to the size of the graph for
	  // most graphs we are interested in.
	  if (real_dist > absolute_error.size())
	    {
	      size_t over = real_dist - absolute_error.size();
	      for (size_t k = 0; k < over; ++k)
		{
		  absolute_error.push_back(0);
		  num_pairs.push_back(0);
		}
	    }

	  if (i != j)
	    {
	      absolute_error[real_dist - 1] += abs_error;
	      ++num_pairs[real_dist - 1];
	    }
	}
    }

  for (size_t i = 0; i < absolute_error.size(); ++i)
    {
      if (num_pairs[i] != 0)
	absolute_error[i] = absolute_error[i] / num_pairs[i];

      cout<<i + 1<<": "<<absolute_error[i]<<" over "<<num_pairs[i]<<" pairs.\n";
    }

  return 0;
}
