/*
  This file generates an executable which takes a community file (each
  node is assigned an integer community) and uses this to find
  statistics on the community's localization in the tree.
  
  By Aaron Adcock, PhD candidate at Stanford University
  Oct. 2013 
*/

#include "graph_lib_boost.hpp"
#include "tree_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
#include <sstream>
#include <string.h>

template <class T>
double vector_max(vector<T> values);

vector<size_t> size_of_component(vector<size_t> & components);

pair< pair<double, double>, pair<double, double> > max_coverage_by_branch(TreeDecomp & TD, vector<v_size_t> & component, int community_label, size_t community_size, vector<int> & community);

void write_data_file (string outputFileName, vector<size_t> data) 
{
  ofstream output;
  output.open(outputFileName.c_str());

  for(v_size_t i = 0; i < data.size(); i++)
    output<<(i+1)<<"\t"<<data[i]<<"\n";
  
  output.close();
}

int main(int argc, char *argv[])
{
  string graph_file = "";
  string tree_file = "";
  string community_file = "";
  string output_prefix = "community";

  if(argc<2)
    {
      cout<<"The options are (-g -t -community are required): \n-g <input_graph_file> \n-t <input_tree_decomposition> \n-o <output_file_name_prefix> \n-community <community_file>\n";
      return 0;
    }
  
  for (int i = 1; i < argc; i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-g") == 0)
	{
	  ss<<argv[i+1];
	  ss>>graph_file;
	}
      else if (strcmp(argv[i], "-t") == 0)
	{
	  ss<<argv[i+1];
	  ss>>tree_file;
	}	
      else if (strcmp(argv[i],"-o") == 0)
	{
	  ss<<argv[i+1];
	  ss>>output_prefix;
	}
      else if (strcmp(argv[i], "-community") == 0)
	{
	  ss<<argv[i+1];
	  ss>>community_file;
	}
    }
  ////

  if(graph_file.length() == 0 || tree_file.length() == 0 || community_file.length() == 0)
    {
      cout<<"You must provide an input graph/tree/community file\n";
      return(0);
    }
  else
    cout<<"Graph: "<<graph_file<<"\nTree: "<<tree_file<<"\nCommunity: "<<community_file<<"\n";

  cout.flush();
 
  // Load Graph
  Graph G = loadGraph(graph_file, "\t", false);
  G       = connected(G);
  v_size_t size = num_vertices(G);
  cout<<"Connected Network has "<<size<<" vertices.\n";

  // Load Tree decomposition
  TreeDecomp TD = loadTreeDecomp(G, tree_file);
  vector< vector<v_size_t> > bags = TD.get_bags();

  // Load Communities
  vector<int> community;
  load_color_file(community_file, community);
  cout<<"Color file loaded.\n";

  if (community.size() != size)
    {
      cerr<<"Community vector length: "<<community.size()<<"\n";
      cerr<<"Graph size: "<<size<<"\n";
      cerr<<"Community vector not the right size. Exiting. \n";
      exit(EXIT_FAILURE);
    }

  // Get sorted list of community labels.  I assume one is the
  // "unknown" label, in my data this is typically 0.
  vector<int> community_list = community;

  sort(community_list.begin(), community_list.end());
  vector<int>::iterator it = unique(community_list.begin(), community_list.end());
  community_list.resize(distance(community_list.begin(), it));
 
  // Map from communities to relative frequency in network
  map<int, double> community_frequency;
  map<int, size_t> community_sizes;

  cout<<"Community frequencies.\n";
  for (size_t i = 0; i < community.size(); ++i)
    {
      if (community_frequency.find(community[i]) != community_frequency.end())
	{
	  community_frequency[community[i]] += 1.0 / community.size();
	  community_sizes[community[i]] += 1;
	}
      else
	{    
	  community_frequency[community[i]] = 1.0 / community.size();
	  community_sizes[community[i]] = 1;
	}
    }
  // Three ways to find bags which contain the community in a
  // meaningful way.  First, the community forms a plurality of nodes.
  // Second, the community forms a majority of nodes.  Third, the
  // community is more frequent in the bag than in the network.  The
  // first two definitions will typically place a bag in a unique
  // location excepting ties (there may be many ties in the small
  // bags).  The frequency definition could potentially place any bag
  // in multiple communities.
  vector<size_t> dummy;
  map<int, vector<size_t> > bag_plural_community;
  map<int, vector<size_t> > bag_major_community;
  map<int, vector<size_t> > bag_unlikely_community;
  
  // Initialize maps
  for (size_t i = 0; i < community_list.size(); ++i)
    {
      bag_plural_community[community_list[i]] = dummy;
      bag_major_community[community_list[i]] = dummy;
      bag_unlikely_community[community_list[i]] = dummy;
    }

  // Place bags.
  cout<<"Begin processing bags\n";
  for (size_t i = 0; i < bags.size(); ++i)
    {
      vector<v_size_t> current_bag = bags[i];
      map<int, double> community_bag_freq;
      // For each bag, calculate frequency of each community
      for (size_t j = 0; j < current_bag.size(); ++j)
	{
	  int comm = community[current_bag[j]];
	  if (community_bag_freq.find(comm) != community_bag_freq.end())
	    community_bag_freq[comm] += 1.0 / current_bag.size();
	  else
	    community_bag_freq[comm] = 1.0 / current_bag.size();
	}

      map<int, double>::iterator mit;
      vector<int> max_communities;
      double max_frequency = 0;
      // Find plural frequency
      for (mit = community_bag_freq.begin(); mit != community_bag_freq.end(); ++mit)
	{
	  if (mit->second > max_frequency)
	    {
	      max_communities.clear();
	      max_communities.push_back(mit->first);
	      max_frequency = mit->second;
	    }
	  else if (mit->second == max_frequency)
	    max_communities.push_back(mit->first);	    
	} 
      // If plural community, add bag to plural community vector
      for (size_t j = 0; j < max_communities.size(); ++j)
	bag_plural_community[max_communities[j]].push_back(i);

      // If majority, add bag to major community vector
      if (max_frequency >= 0.5)
	for (size_t j = 0; j < max_communities.size(); ++j)
	  bag_major_community[max_communities[j]].push_back(i); 

      // If bag frequency is larger than expected, add to expected vector
      for (mit = community_frequency.begin(); mit != community_frequency.end(); ++mit)
	{
	  if (mit->second < community_bag_freq[mit->first])
	    bag_unlikely_community[mit->first].push_back(i);
	}
    }

  // What do I want to output here?  Number of components.  Largest
  // components and their contents (how much of the community do they
  // cover, list of members.  Structure of the subtree?  Persistence
  // on the community?  Store file of all components.

  // Loop through community labels. Form subgraph for each community.
  // Write out files.  Then display summary to terminal
  for (size_t i = 0; i < community_list.size(); ++i)
    {
      vector<Vert> plural_bags, major_bags, likely_bags;
      int curr_comm = community_list[i];

      Graph T = TD.get_tree();
      plural_bags = convert_indices_to_vertices(T, bag_plural_community[curr_comm]);
      major_bags = convert_indices_to_vertices(T, bag_major_community[curr_comm]);
      likely_bags = convert_indices_to_vertices(T, bag_unlikely_community[curr_comm]);
      
      Graph G_plural = subset(T, plural_bags); 
      Graph G_major  = subset(T, major_bags);
      Graph G_likely = subset(T, likely_bags);
      
      vector <v_size_t> plural_comp = components(G_plural);
      vector <v_size_t> major_comp = components(G_major);
      vector <v_size_t> likely_comp = components(G_likely);
      
      size_t num_plural_comp = vector_max(plural_comp) + 1;
      size_t num_major_comp = vector_max(major_comp) + 1;
      size_t num_likely_comp = vector_max(likely_comp) + 1;

      vector<size_t> size_plural_comps = size_of_component(plural_comp);
      vector<size_t> size_major_comps = size_of_component(major_comp);
      vector<size_t> size_unlikely_comps = size_of_component(likely_comp);

      size_t num_plural_greater_5 = 0;
      size_t num_major_greater_5 = 0;
      size_t num_likely_greater_5 = 0;

      for (size_t j = 0; j < size_plural_comps.size(); ++j)
	if (size_plural_comps[j] >= 5)
	  num_plural_greater_5 += 1;

      for (size_t j = 0; j < size_major_comps.size(); ++j)
	if (size_major_comps[j] >= 5)
	  num_major_greater_5 += 1;

      for (size_t j = 0; j < size_unlikely_comps.size(); ++j)
	if (size_unlikely_comps[j] >= 5)
	  num_likely_greater_5 += 1;

      pair< pair<double, double>, pair<double, double> > plural_coverage = max_coverage_by_branch(TD, plural_comp, community_list[i], community_sizes[community_list[i]], community);
      pair< pair<double, double>, pair<double, double> > major_coverage = max_coverage_by_branch(TD, major_comp, community_list[i], community_sizes[community_list[i]], community);
      pair< pair<double, double>, pair<double, double> > likely_coverage = max_coverage_by_branch(TD, likely_comp, community_list[i], community_sizes[community_list[i]], community);

      cout<<"Community "<<community_list[i]<<"\n";
      cout<<"Community Frequency: "<<community_frequency[community[i]]<<"\n";
      cout<<"Number of tree components:\n";
      cout<<"Plural: "<<num_plural_comp<<"\n";
      cout<<"Major: "<<num_major_comp<<"\n";
      cout<<"Likely: "<<num_likely_comp<<"\n";      
      cout<<"Number of tree components with more than 5 bags:\n";
      cout<<"Plural: "<<num_plural_greater_5<<"\n";
      cout<<"Major: "<<num_major_greater_5<<"\n";
      cout<<"Likely: "<<num_likely_greater_5<<"\n";      
      cout<<"Max coverage by single component:\n";
      cout<<"Plural (largest): "<<plural_coverage.first.first<<" coverage of community, "<<plural_coverage.second.first<<" of nodes in component.\n";
      cout<<"Plural (2nd largest): "<<plural_coverage.first.second<<" coverage of community, "<<plural_coverage.second.second<<" of nodes in component.\n";
      cout<<"Major (largest): "<<major_coverage.first.first<<" coverage of community, "<<major_coverage.second.first<<" of nodes in component.\n";
      cout<<"Major (2nd largest): "<<major_coverage.first.second<<" coverage of community, "<<major_coverage.second.second<<" of nodes in component.\n";
      cout<<"Likely (largest): "<<likely_coverage.first.first<<" coverage of community, "<<likely_coverage.second.first<<" of nodes in component.\n\n";
      cout<<"Likely (2nd largest): "<<likely_coverage.first.second<<" coverage of community, "<<likely_coverage.second.second<<" of nodes in component.\n\n";
    }
}

// Finds max of a vector of value types
template <class T> 
double vector_max(vector<T> values)
{
  double max = values[0];
  for (size_t i = 1; i < values.size(); ++i)
    if (values[i] > max)
      max = (double) values[i];

  return max;
}

// Takes a vector of labels and returns the size of each label.
// Assumes that components has labels which are nonnegative integers,
// sequential, and begin at zero.
vector<size_t> size_of_component(vector<size_t> & components)
{
  vector<size_t> unique_components = components;
  sort(unique_components.begin(), unique_components.end());
  vector<size_t>::iterator it = unique(unique_components.begin(), unique_components.end());
  unique_components.resize(distance(unique_components.begin(), it));

  unique_components.assign(unique_components.size(), 0);

  for (size_t i = 0; i < components.size(); ++i)
    unique_components[components[i]] += 1;

  return unique_components;
}

// Looks at all the tree components (branches) of the tree for the
// given community label.  Finds the branch with the most coverage,
// returns that value (a number between 0 and 1).  Note that this
// algorithm assumes that the component vector is tailored
// appropriately to both the community label and the community vector.
// It does not check to see if this is true.  GIGO.
pair<pair<double, double>, pair<double, double> > max_coverage_by_branch(TreeDecomp & TD, vector<v_size_t> & component, int community_label, size_t community_size, vector<int> & community)
{
  // First get unique components
  vector<size_t> unique_components;
  unique_components.insert(unique_components.end(), component.begin(), component.end());
  sort(unique_components.begin(), unique_components.end());
  vector<size_t>::iterator it = unique(unique_components.begin(), unique_components.end());
  unique_components.resize(distance(unique_components.begin(), it)); 

  vector<double> coverage (unique_components.size(), 0);
  vector<double> fraction_of_branch (unique_components.size(), 0);
  double max_coverage = 0;
  double second_coverage = 0;
  double max_coverage_fraction = 0;
  double second_coverage_fraction = 0;
  for (size_t i = 0; i < unique_components.size(); ++i)
    {
      vector<v_size_t> nodes_in_branch;
      for (size_t j = 0; j < component.size(); ++j)
	{
	  if (component[j] == unique_components[i])
	    {
	      vector<v_size_t> bag = TD.get_bag(j);
	      nodes_in_branch.insert(nodes_in_branch.end(), bag.begin(), bag.end());
	    }
	}
      sort(nodes_in_branch.begin(), nodes_in_branch.end());
      vector<v_size_t>::iterator vec_it = unique(nodes_in_branch.begin(), nodes_in_branch.end());
      nodes_in_branch.resize(distance(nodes_in_branch.begin(), vec_it));
      for (size_t j = 0; j < nodes_in_branch.size(); ++j)
	{
	  if (community[nodes_in_branch[j]] == community_label)
	    coverage[i] += 1.0;
	}
      fraction_of_branch[i] = coverage[i] / nodes_in_branch.size();
      coverage[i] /= community_size;

      if (coverage[i] > max_coverage)
	{	 
	  second_coverage = max_coverage;
	  second_coverage_fraction = max_coverage_fraction;
	  max_coverage = coverage[i];
	  max_coverage_fraction = fraction_of_branch[i];
	}
      else if (coverage[i] > second_coverage)
	{
	  second_coverage = coverage[i];
	  second_coverage_fraction = fraction_of_branch[i];
	}
    }
  
  pair<double, double> coverage_stats(max_coverage, second_coverage);
  pair<double, double> coverage_fractions(max_coverage_fraction, second_coverage_fraction);
  pair<pair<double, double>, pair<double, double> > ret_val(coverage_stats, coverage_fractions);
  return  ret_val;
}
