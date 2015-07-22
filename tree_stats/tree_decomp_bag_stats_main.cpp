/*
  The main purpose of this function is to produce relevant statistics
  on the structure of tree decomposition bags.  It will read in a tree
  decomposition, the original graph, and produces several files.  One
  file will contain statistics that are printed in a format that is
  ready for LaTex.  The second and third files will contain stats on
  individual bags which can be used to produce scatterplots, perhaps
  correlating the width of a bag or the eccentricity of a bag, and the
  internal bag structure statistics as well as the overlap statistics
  of adjacent bags.  The final files will be visualization files of
  internal bag structures and overlaps.  The will be colored by the
  coloring provided by a color file or by k-core or degree
  distributions which can be calculated easily.

  By Aaron Adcock
  PhD Candidate Stanford University
  2014
 */
#include "../lib/graph_lib_boost.hpp"
#include "../lib/tree_lib_boost.hpp"
#include <set>

template <class T>
double vector_min(vector<T> values);
template <class T>
double vector_max(vector<T> values);
template <class T>
double vector_var(vector<T> values);
template <class T>
double vector_median(vector<T> values);
template <class T> 
double vector_avg(vector<T> values);

int main(int argc, char * argv[])
{
  // Input variables!
  string input_treename = "";
  string input_graphname = "";
  string output_file_prefix = "graph_stats";
  string color_file = "";
  string bag_list_file = "";
  bool calc_stats = false;
  bool visualize_bags = false;
  bool visualize_branch = false;
  bool directed = false;
  size_t ecc_thresh = (size_t) - 1;
  bool visualize_graph = false;

  if(argc < 2)
    {
      cout<<"-g <input_graph> \n-t <input_tree> \n-o <output_file_prefix> \n-directed <directed graph input> \n-color <color_file, must match graph size> \n-calc_stats <Calculates stats on tree> \n-visualize_bags <outputs a graph_viz file for each bag provided by bag_list> \n-visualize_branch <outputs a graph_viz file for the subgraph induced by all nodes in all bags provided by bag_list> \n-bag_list <file containing single line of space seperated integers which correspond to the indices of bags in the tree file that you want visualized> \n-visualize_graph \n";
      return(-1);
    }
  
  // Command line flags
  for (int i = 1; i < argc; ++i)
    {
      if (strcmp(argv[i], "-t") == 0)
	input_treename = argv[i+1];
      else if (strcmp(argv[i], "-g") == 0)
	input_graphname = argv[i+1];
      else if (strcmp(argv[i], "-o") == 0)
	output_file_prefix = argv[i+1];
      else if (strcmp(argv[i], "-directed") == 0)
	directed = true;
      else if (strcmp(argv[i], "-color") == 0)
	color_file = argv[i+1];
      else if (strcmp(argv[i], "-calc_stats") == 0)
	calc_stats = true;
      else if (strcmp(argv[i], "-visualize_bags") == 0)
	visualize_bags = true;
      else if (strcmp(argv[i], "-visualize_branch") == 0)
	visualize_branch = true;
      else if (strcmp(argv[i], "-bag_list") == 0)
	bag_list_file = argv[i+1];
      else if (strcmp(argv[i], "-ecc_thresh") == 0)
	{
	  stringstream ss;
	  ss<<argv[i+1];
	  ss>>ecc_thresh; 
	}
      else if (strcmp(argv[i], "-visualize_graph") == 0)
	visualize_graph = true;
    }

  if (!calc_stats && !visualize_bags && !visualize_branch)
    {
      cerr<<"You must pick at least one of -calc_stats, -visualize_bags, -visualize_branch. \n";
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

  vector<double> color;
  vector<double> cores;
  if (color_file.length() != 0)
    {
      load_color_file(color_file, color);
      vector<int> temp  = k_core(G);
      cores.assign(temp.begin(), temp.end());
    }
  else
    {
      vector<int> temp = k_core(G); 
      color.assign(temp.begin(), temp.end());
      cores.assign(temp.begin(), temp.end());
    }
  v_size_t num_nodes = num_vertices(G);

  cout<<"Graph loaded. The connected component has "<<num_nodes<<" nodes.\n";
  
  // Load TD from INDDGO tree file
  TreeDecomp TD(loadTreeDecomp(G, input_treename));
  
  cout<<"Tree decomposition has "<<TD.get_num_bags()<<" bags.\n";
  
  size_t num_tree_nodes = TD.get_num_bags();
  vector<size_t> bags_to_viz;

  // Calculate average k-core, median k-core
  

  // load list of bags to visualize, if command line asks for it
  if ((visualize_bags || visualize_branch) && bag_list_file.length() != 0)
    {
      // Reads first noncommented line of bags file and puts contents
      // in input_distribution
      ifstream bags_file;

      cout<<"Reading bag list\n";
      bags_file.open(bag_list_file.c_str());
      bool line_not_found = true;
      while (bags_file.good() && line_not_found)
	{
	  string line;
	  getline(bags_file, line);

	  // Skip lines beginning with #
	  if (line.substr(0,1) != "#")
	    {
	      line_not_found = false;

	      stringstream bags_string(line);
	      size_t count = 0;
	      // Only keep the first num_tree_nodes values
	      while (count < num_tree_nodes && bags_string.good())
		{
		  size_t temp;
		  bags_string>>temp;
		  bags_to_viz.push_back(temp);
		  count++;
		}
	    }
	}
      bags_file.close();
    }

  vector <v_size_t> dummy;
  vector< vector<v_size_t> > vertices_in_bag(bags_to_viz.size(), dummy);
  size_t total_verts_in_branch = 0;
  // If visualization is set, find vertices associated with each bag in list
  if (visualize_bags || visualize_branch)
    {
      cout<<"Putting vertices into bags\n";
      for (int i = 0; i < bags_to_viz.size(); ++i)
	{
	  vertices_in_bag[i] = TD.get_bag(bags_to_viz[i]);
	  total_verts_in_branch += vertices_in_bag[i].size();
	}
    }

  // Run through list of bags and produce a subgraph
  if (visualize_bags)
    {
      cout<<"Finding subgraph associated with each bag\n";
      Graph G_bag;
      for (int i = 0; i < vertices_in_bag.size(); ++i)
	{
	  G_bag = subset(G, vertices_in_bag[i]);
	  string output_bag_file = output_file_prefix;
	  string bag_num;
	  output_bag_file.append("_bag_");
	  stringstream ss;
	  ss<<bags_to_viz[i];
	  ss>>bag_num;
	  output_bag_file.append(bag_num);
	  output_bag_file.append(".dot");
	  vector<int> square_nodes;
	  write_graphviz(G_bag, output_bag_file, color, square_nodes, false);
	}
    }

  if (visualize_branch)
    {
      cout<<"Finding subgraph associated with branch\n";
      Graph G_branch;
      vector<v_size_t> branch;
      branch.reserve(total_verts_in_branch);
      for (int i = 0; i < vertices_in_bag.size(); ++i)
	branch.insert(branch.end(), vertices_in_bag[i].begin(), vertices_in_bag[i].end());
      
      G_branch = subset(G, branch);
      string output_branch_file = output_file_prefix;
      string bag_num;
      output_branch_file.append("_branch.dot");
      vector<int> square_nodes;
      write_graphviz(G_branch, output_branch_file, color, square_nodes, false);
    }

  if (visualize_graph)
    {
      cout<<"Producing graph visualization of bags/branches\n";
      string output_graph_file = output_file_prefix;
      output_graph_file.append("_graph.dot");
      vector<double> nodes_in_branch (num_nodes, 0); 
      vector<int> square_nodes;

      for (size_t i = 0; i < vertices_in_bag.size(); ++i)
	for (size_t j = 0; j < vertices_in_bag[i].size(); ++j)
	  nodes_in_branch[vertices_in_bag[i][j]] = 1.0;
      
      write_graphviz(G, output_graph_file, nodes_in_branch, square_nodes, false);
    }
  
  if (calc_stats)
    {
      cout<<"Calculating stats\n";
      // Calculate clique-iness/adjacent overlap of all bags
      Graph T = TD.get_tree();
      vector< vector<v_size_t> > bags = TD.get_bags();
      vector< double > avg_k_core(bags.size(), 0);
      vector< double > med_k_core(bags.size(), 0);
      vector< double > bag_density(bags.size(), 0);
      vector< double > bag_overlap(bags.size(), 0);
      vector< v_size_t > bag_ecc(bags.size(),0);
      property_map<Graph, vertex_index_t>::type index = get(vertex_index, T);

      bag_ecc = TD.get_tree_eccentricity();

      for (int i = 0; i < bags.size(); ++i)
	{
	  Graph G_bag = subset(G, bags[i]);
	  Vert bag = vertex(i, T);
	  size_t num_links = num_edges(G_bag);
	  vector<v_size_t> bag_nodes = bags[i];
	  size_t bag_size = bag_nodes.size();

	  graph_traits<Graph>::adjacency_iterator ait, aite;

	  vector<double> bag_cores;
	  for (int j = 0; j < bags[i].size(); ++j)
	    bag_cores.push_back(cores[bags[i][j]]);

	  avg_k_core[i] = vector_avg(bag_cores);
	  med_k_core[i] = vector_median(bag_cores);	    
	  
	  double avg_intersection = 0;
	  for (tie(ait, aite) = adjacent_vertices(bag, T); ait != aite; ++ait)
	    {
	      int size_intersection = 0;
	      Vert adj_bag = *ait;
	      size_t bag_index = index[adj_bag];
	      vector<v_size_t> adj_bag_nodes = bags[bag_index];
	      
	      // Not an easy hash table implementation for pre c++11,
	      // so using set (a binary search tree). Put verts in
	      // adjacent bag in set, then check if verts in bag are
	      // in set.
	      set<v_size_t> test_set;
	      test_set.insert(adj_bag_nodes.begin(), adj_bag_nodes.end());

	      for (int j = 0; j < bag_size; ++j)
		if (test_set.find(bag_nodes[j]) != test_set.end())
		  size_intersection++;

	      if (bag_size <= bags[bag_index].size())
		avg_intersection += (double) size_intersection / bag_size;
	      else
		avg_intersection += (double) size_intersection / bags[bag_index].size();
	    }
	  
	  // Averaged over all neighboring bags in tree.
	  bag_overlap[i] = (double) avg_intersection / (out_degree(bag, T));

	  // may need a divide by two here, Clique-iness!
	  bag_density[i] = (double) num_links / (bags[i].size() * (bags[i].size() - 1));
	}

      vector<size_t> bag_sizes = TD.get_bag_sizes();

      // Produce output First, produce bag color files (currently
      //unclear how to use them) and visualization files
      string overlap_file = output_file_prefix;
      overlap_file.append("_overlap.color");
      string density_file = output_file_prefix;
      density_file.append("_density.color");

      write_color_file(overlap_file, input_graphname, bag_overlap);
      write_color_file(density_file, input_graphname, bag_density);

      overlap_file = output_file_prefix;
      overlap_file.append("_overlap.dot");
      density_file = output_file_prefix;
      density_file.append("_density.dot");
      // Colors must be between 0 and 1 (ie bag_overlap/bag_density)
      write_scaled_graphviz(T, overlap_file, bag_overlap, bag_sizes, directed);
      write_scaled_graphviz(T, density_file, bag_density, bag_sizes, directed);

      // Produce gnuplot scatter_plot files
      
      string output_scatter_plot = output_file_prefix;
      output_scatter_plot.append("_scatter_plot.data");

      ofstream output_file;
      output_file.open(output_scatter_plot.c_str());

      output_file<<"# This contains data from "<<input_graphname<<" and its tree decomposition.  The first column gives the index of the bag.  The second column gives the size of the bag.  The third column gives the eccentricity (as a node in the tree decomp.) of the bag.  The fourth column gives the average overlap among the bag and its neighbors.  The fifth column gives the relative density of the bag, 6th column degree of bag in tree, 7th column avg k-core of bag, 8th column med k-core of bag. \n";
      output_file<<"# Basic stats in format: max min avg med var \n";
      output_file<<"# Basic bag size stats:  \n";
      output_file<<"# "<<vector_max(bag_sizes)<<" "<<vector_min(bag_sizes)<<" "<<vector_avg(bag_sizes)<<" "<<vector_median(bag_sizes)<<" "<<vector_var(bag_sizes)<<"\n";
      output_file<<"# Basic bag eccentricity stats:  \n";
      output_file<<"# "<<vector_max(bag_ecc)<<" "<<vector_min(bag_ecc)<<" "<<vector_avg(bag_ecc)<<" "<<vector_median(bag_ecc)<<" "<<vector_var(bag_ecc)<<"\n";
      output_file<<"# Basic overlap stats:  \n";
      output_file<<"# "<<vector_max(bag_overlap)<<" "<<vector_min(bag_overlap)<<" "<<vector_avg(bag_overlap)<<" "<<vector_median(bag_overlap)<<" "<<vector_var(bag_overlap)<<"\n";
      output_file<<"# Basic density stats:  \n";
      output_file<<"# "<<vector_max(bag_density)<<" "<<vector_min(bag_density)<<" "<<vector_avg(bag_density)<<" "<<vector_median(bag_density)<<" "<<vector_var(bag_density)<<"\n";
    
      for (int i = 0; i < bag_sizes.size(); ++i)
	output_file<<index[i]<<"\t"<<bag_sizes[i]<<"\t"<<bag_ecc[i]<<"\t"<<bag_overlap[i]<<"\t"<<bag_density[i]<<" "<<out_degree(vertex(i, T), T)<<"\t"<<avg_k_core[i]<<"\t"<<med_k_core[i]<<"\n";

      output_file.close();
    }	     

  

  return -1;
}

template <class T>
double vector_avg(vector<T> values)
{
  double avg = 0;
  for (int i = 0; i < values.size(); ++i)
    avg += (double) values[i];

  return avg / (double) values.size();
}

// double vector_avg(vector<int> values)
// {
//   double avg = 0;
//   for (int i = 0; i < values.size(); ++i)
//     avg += (double) values[i];

//   return avg / (double) values.size();
// }

template <class T>
double vector_median(vector<T> values)
{
  sort(values.begin(), values.end());
  size_t size = values.size();
  if (size % 2 == 0)
    return (double) (values[size / 2] + values[size / 2 - 1]) / 2.0;
  else
    return (double) values[size / 2];
}

// double vector_median(vector<int> values)
// {
//   sort(values.begin(), values.end());
//   size_t size = values.size();
//   if (size % 2 == 0)
//     return (double) (values[size / 2] + values[size / 2 - 1]) / 2.0;
//   else
//     return (double) values[size / 2];
// }

template <class T>
double vector_var(vector<T> values)
{
  double square_sum = 0;
  double sum = 0;

  for (int i = 0; i < values.size(); ++i)
    {
      square_sum += (double) values[i] * values[i];
      sum += (double) values[i];
    }

  return (square_sum / values.size()) - (sum * sum) / (values.size() * values.size());
}

// double vector_var(vector<int> values)
// {
//   double square_sum = 0;
//   double sum = 0;

//   for (int i = 0; i < values.size(); ++i)
//     {
//       square_sum += (double) values[i] * values[i];
//       sum += (double) values[i];
//     }

//   return (square_sum / values.size()) - (sum * sum) / (values.size() * values.size());
// }

template <class T> 
double vector_max(vector<T> values)
{
  double max = values[0];
  for (int i = 1; i < values.size(); ++i)
    if (values[i] > max)
      max = (double) values[i];

  return max;
}

// double vector_max(vector<int> values)
// {
//   double max = (double) values[0];
//   for (int i = 1; i < values.size(); ++i)
//     if (values[i] > max)
//       max = (double) values[i];

//   return max;
// }

template <class T>
double vector_min(vector<T> values)
{
  double min = values[0];
  for (int i = 1; i < values.size(); ++i)
    if (values[i] < min)
      min = (double) values[i];

  return min;
}

// double vector_min(vector<int> values)
// {
//   double min = (double) values[0];
//   for (int i = 1; i < values.size(); ++i)
//     if (values[i] < min)
//       min = (double) values[i];

//   return min;
// }
