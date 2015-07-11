/*
  This file provides functions associated with learning on tree
  decompositions. The basic functions simply convert tree
  decompositions into input for other machine learning algorithms (ie
  feature vectors, etc).  Eventually I plan to add algorithms which
  learn "natively" on the tree decompositions.

  Aaron Adcock
  PhD Candidate Stanford University
  2013
 */

#include "tree_lib_boost.hpp"
#include "../graph_lib_boost.hpp"

/* Given bags_to_keep, it produces a modified vertex_location list which removes 
   all bags which are not in the list.
   TODO: If we have memory issues, might want to consider removing
   elements from feature_vector.  If this results in a performance
   loss, then maybe change to feature_list.
*/

void createBagVectors(TreeDecomp TD, vector< vector<v_size_t> > & feature_vector, vector<size_t> & bags_to_keep)
{
  size_t num_bags = bags_to_keep.size();
  size_t num_nodes = num_vertices(*TD.get_graph_pointer());
  vector< vector<size_t> > vertex_locations = TD.get_vertex_locations();
  vector<v_size_t> current_vertex_bags;

  // Clear feature vector and assign empty vectors to each slot
  feature_vector.clear();
  feature_vector.assign(num_nodes, current_vertex_bags);

  sort(bags_to_keep.begin(), bags_to_keep.end());
  
  for (int i = 0; i < feature_vector.size(); ++i)
    {
      current_vertex_bags.swap(vertex_locations[i]);

      if (current_vertex_bags.size() == 0)
	cout<<"WARNING: Node "<<i+1<<" does not have a bag location listed. Check TD file.\n";

      for (int j = 0; j < current_vertex_bags.size(); ++j)
	{
	  v_size_t bag = current_vertex_bags[j];
	  int index = binarySearch(bags_to_keep, bag);

	  if (index >= 0)
	    feature_vector[i].push_back(bag);
	}
    }
}

/* The files are large and mostly sparse, so we only write the bag indices */
void writeFeatureVectors(string filename, vector< vector<v_size_t> > & feature_vector)
{
  ofstream output_file;
  output_file.open(filename.c_str());

  size_t num_nodes = feature_vector.size();
  vector<v_size_t> current_node_bags;
  
  for (int i = 0; i < num_nodes; ++i)
    {
      current_node_bags.swap(feature_vector[i]);

      // Check if each node is listed in at least one bag.  If not, send warning
      if (current_node_bags.size() == 0)
	cout<<"WARNING: Node "<<i+1<<" is not listed in any bag.\n";

      for (int j = 0; j < current_node_bags.size(); ++j)
	{
	  output_file<<current_node_bags[j]<<" ";
	}
      output_file<<"\n";
    }
  output_file.close();
}
