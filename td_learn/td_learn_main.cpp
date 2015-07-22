/*
  This main file is used to implement and test learning algorithms
  using a network and an associated tree decomposition.  Currently (as
  of Nov. 2013), the tree decompositions are generated using the
  INDDGO open source software and are read in from a modified dimacs
  file. The graph associated with the tree decomposition must also be
  specified (currently, there aren't any identifiers in the tree file
  which name the graph file that it was generated from).

  TODO: Need to get a more appropriate random function, for testing
  purposes the rand() is usable, but not for actual results.

  By Aaron Adcock
  PhD Candidate Stanford University
  2013
 */
#include "../lib/graph_lib_boost.hpp"
#include "../lib/tree_lib_boost.hpp"

int main(int argc, char * argv[])
{
  // Input variables!
  string input_treename = "";
  string input_graphname = "";
  string output_filename = "a.out";
  size_t lower_threshold = 0;
  size_t upper_threshold = (size_t) -1;
  bool take_rand_sample = false;
  double sample_size = 0.1;
  
  // Command line flags
  for (int i = 1; i < argc; ++i)
    {
      if (strcmp(argv[i], "-t") == 0)
	input_treename = argv[i+1];
      else if (strcmp(argv[i], "-g") == 0)
	input_graphname = argv[i+1];
      else if (strcmp(argv[i], "-o") == 0)
	output_filename = argv[i+1];
      else if (strcmp(argv[i], "-lt") == 0)
	lower_threshold = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-ut") == 0)
	upper_threshold = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-rand") == 0)
	{
	  take_rand_sample = true;
	  sample_size = atof(argv[i+1]);
	}
    }

  if (sample_size > 1 || sample_size <= 0)
    {
      cerr<<"You must provide a sample size between 0 and 1.  Exiting.\n";
      exit(EXIT_FAILURE);
    }

  if (input_treename.length() == 0 || input_graphname.length() == 0)
    {
      cerr<<"You must provide a tree file and a graph file. Exiting.\n";
      exit(EXIT_FAILURE);
    }
  
  // Load Graph, find connected component (INDDGO only produces
  // tree decompositions from connected components)
  Graph G = loadGraph(input_graphname);
  G = connected(G);

  v_size_t num_nodes = num_vertices(G);

  cout<<"Graph loaded. The connected component has "<<num_nodes<<" nodes.\n";
  
  // Load TD from INDDGO tree file
  TreeDecomp TD(loadTreeDecomp(G, input_treename));
  
  cout<<"Tree decomposition has "<<TD.get_num_bags()<<" bags.\n";
  
  vector<size_t> bag_sizes = TD.get_bag_sizes();
  vector<size_t> histogram(10, 0);

  // Create a bag size histogram
  cout<<"Bag Histogram: \n";
  int max = 0;
  int min = 30000;
  vector<size_t> bags_to_keep;
  for (int i = 0; i < bag_sizes.size(); ++i)
    {
      size_t bag_size = bag_sizes[i];
      if (bag_size > max)
	max = bag_size;

      if (bag_size < min)
	min = bag_size;

      if (!take_rand_sample)
	{
	  if (bag_size >= lower_threshold && bag_size <= upper_threshold)
	    bags_to_keep.push_back(i);
	}
      else
	{
	  double r = (double) rand() / (double) RAND_MAX;
	  if (r <= sample_size)
	    bags_to_keep.push_back(i);
	}

    }
  
  cout<<"Using given threshold values, we are producing features from "<<bags_to_keep.size()<<" bags.\n";

  for (int i = 0; i < bag_sizes.size(); ++i)
    {
      int index = 9 * (double) (bag_sizes[i] - min) / (double) (max - min);
      ++histogram[index];
    }

  int prev_val = min;
  for (int i = 0; i < histogram.size(); ++i)
    {
      int curr_val = (i + 1) * (double) (max - min) / histogram.size() + min;

      cout<<prev_val<<"-"<<curr_val<<": "<<histogram[i]<<"\n";
      prev_val = curr_val;
    }

  // Create feature vectors.  Two options: 
  // 1. Basic vectors just keep
  // track of whether or not bag is in one of the landmark set (either
  // chosen randomly or through a set of thresholds).
  // 2. Keep track of distance in tree (in graph?) from landmark points <-
  vector< vector<v_size_t> > feature_vector;

  createBagVectors(TD, feature_vector, bags_to_keep);

  cout<<"Feature vector contains features for "<<feature_vector.size()<<" nodes.\n";
  
  int count = 0;
  max = 0;
  for (int i = 0; i < feature_vector.size(); ++i)
    {
      if (feature_vector[i].size() == 0)
	++count;
      
      if (feature_vector[i].size() > max)
	max = feature_vector[i].size();
    }

  cout<<"Vertices missing a bag: "<<count<<".\n";
  cout<<"Maximum number of simultaneous bags by a vertex: "<<max<<".\n";

  writeFeatureVectors(output_filename, feature_vector);
}
