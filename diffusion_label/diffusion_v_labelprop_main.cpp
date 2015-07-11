/*
  This function runs a diffusion at every node and checks to see how
  much of the neighborhood is the same label as the seed node.

  Aug. 20, 2014
  By Aaron Adcock
  Ph.D. EE, Stanford University
*/

#include "graph_lib_boost.hpp"
#include "tree_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <string.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>

class compare_second_greater
{
public:
  bool operator()(pair<Vert, int> first_pair, pair<Vert, int> second_pair)
  {
    return first_pair.second > second_pair.second;
  }
};

double f1_score(double p, double r);

void read_label_file(string input_label_file, map<v_size_t, int> & labels, map<v_size_t, bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> & integer_to_label);

vector<int> tree_core(Graph & G, TreeDecomp & T);

vector<int> tree_d_core(Graph G);

Graph get_labeled_cores(Graph & G, int max_core);

Graph get_labeled_periphery(Graph & G, int min_core);

Graph get_labeled_tree_cores(Graph & G, TreeDecomp & T, int max_core);

Graph get_labeled_tree_d_cores(Graph & G, int max_d);

map<v_size_t, int> convert_labels(Graph & G, map<v_size_t, int> labels);

vector<int> tree_d_core(Graph G, int d);

template <class T>
T vector_max(vector<T> values);

// Main function for driving snapshots
int main(int argc, char *argv[])
{
  // Variable declarations for user input
  string input_file;
  string output_file_prefix = "graph";
  string label_file = "";
  string tree_file = "";
  bool use_tree_core = false;
  bool directed = true;
  bool suppressOutput = false;
  double diffusion_size = .001;
  double alpha  = .001; 
  double sample = 0.1;
  int cluster_label = -1;
  double num_iterations = 10.0;
  bool core_d = false;

  // User input
  if(argc < 2)
    {
      cout<<"-i <input_file> \n-o <output_file_prefix> \n-labels <label file, if blank uses k_core>\n-a <alpha_value for PageRank> \n-e <eps, error value for PageRank, 1 / ~community size>\n-sample <percentage labels to use>\n-d <y=directed graph input> \n-s <y=suppress terminal output> \n-cluster_label <cluster stats to output>\n-tree_file <tree_decomp_file>\n-core_d <Use cores from core decompositions>\n";
      return(-1);
    }

  for (int i = 1; i < argc; ++i)
    {
      stringstream ss;
      if (strcmp(argv[i], "-i") == 0)
	{
	  ss<<argv[i + 1];
	  ss>>input_file;
	}
      else if (strcmp(argv[i], "-o") == 0)
	{
	  ss<<argv[i + 1];
	  ss>>output_file_prefix;
	}	
      else if (strcmp(argv[i], "-labels") == 0)
	{
	  ss<<argv[i + 1];
	  ss>>label_file;
	}
      else if (strcmp(argv[i], "-labels") == 0)
	sample = atof(argv[i + 1]);
      else if (strcmp(argv[i], "-a") == 0)
	alpha = atof(argv[i + 1]);
      else if (strcmp(argv[i], "-e") == 0)
	diffusion_size = atof(argv[i + 1]);
      else if (strcmp(argv[i], "-d") == 0)
	{
	  if (strcmp(argv[i + 1], "n") == 0)
	    directed = false;
	}
      else if (strcmp(argv[i], "-s") == 0)
	{
	  if(strcmp(argv[i + 1], "y") == 0)
	    suppressOutput = true;
	}
      else if (strcmp(argv[i], "-cluster_label") == 0)
	{
	  cluster_label = atoi(argv[i + 1]);
	}
      else if (strcmp(argv[i], "-tree_file") == 0)
	{
	  ss<<argv[i + 1];
	  ss>>tree_file;
	  use_tree_core = true;
	}
      else if (strcmp(argv[i], "-core_d") == 0)
	{
	  core_d = true;
	}
    }

  // Default values
  if (input_file.length() == 0)
    {
      cout<<"You must provide an input file\n";
      return(0);
    }
  else
    cout<<input_file<<"\n";

  if (alpha <= 0)
    {
      cout<<"Alpha must be greater than 0, using default .001\n";
      alpha = .001;
    }

  if (diffusion_size <= 0)
    {
      cout<<"Diffusion size must be greater than 0, using default .005\n";
      diffusion_size = .005;
    }

  // Load graph
  Graph orig_G = loadGraph(input_file, "\t");
  orig_G = connected(orig_G, -1);
  v_size_t size = num_vertices(orig_G);
  
  // If tree file is set, load tree
  TreeDecomp T;
  if (use_tree_core)
    {
      T = loadTreeDecomp(orig_G, tree_file);
    }

  property_map<Graph, vertex_index_t>::type index = get(vertex_index, orig_G);

  if(!suppressOutput)
    cout<<"Size of connected component of original G "<<size<<"\n";

  // Calculate k-core, note that the core # is also stored as a property on each node

  map<v_size_t, int> orig_labels;
  map<v_size_t, bool>  is_labeled;
  map<int, int> label_to_integer;
  map<int, int> integer_to_label;

  vector<int> core;
  if (core_d)
    core = tree_d_core(orig_G);
  else if (use_tree_core) 
    core = tree_core(orig_G, T);
  else
    core = k_core(orig_G);

  if (label_file.length() == 0)
    {	
      graph_traits<Graph>::vertex_iterator vi, vie;
      for (tie(vi, vie) = vertices(orig_G); vi != vie;  ++vi)
	{
	  orig_labels[index[*vi]] = core[index[*vi]];
	  integer_to_label[core[index[*vi]]] = core[index[*vi]];
	}
    }
  else
    read_label_file(label_file, orig_labels, is_labeled, label_to_integer, integer_to_label);

  int max_core = vector_max(core);

  vector<double> avg_precision(max_core + 1, 0.0);
  vector<double> avg_recall(max_core + 1, 0.0);

  int min_core = 30;

  if (use_tree_core)
    min_core = 30;

  for (int nums = 0; nums < num_iterations; ++nums) {
    for (int c = max_core + 1; c > min_core; --c)
      {
	Graph G_cored; 
      
	if (core_d)
	  G_cored = get_labeled_tree_d_cores(orig_G, c);
	else if (use_tree_core) 
	  G_cored = get_labeled_tree_cores(orig_G, T, c);
	else
	  G_cored = get_labeled_periphery(orig_G, c);
      
	property_map<Graph, vertex_index_t>::type cored_index = get(vertex_index, G_cored);
	map<v_size_t, int> cored_labels = convert_labels(G_cored, orig_labels);
	map<v_size_t, int>::iterator map_it;
	// Need to make a map of labels to the associated vertices, then use
	// this to pick the "labeled" data (in the ML sense of labeled) and
	// "unlabeled" data or test set.

	map<int, vector<v_size_t> > communities;
	map<v_size_t, vector<double> > node_to_pagerank;
	// Fill communities map with community => set of nodes in community
	for (map_it = cored_labels.begin(); map_it != cored_labels.end(); ++map_it)
	  {
	    if (!suppressOutput)
	      cout<<"Vert: "<<map_it->first<<" community: "<<map_it->second<<"\n";
	    if (communities.find(map_it->second) == communities.end())
	      {
		vector<v_size_t> empty_vector;
		communities[map_it->second] = empty_vector;
	      }

	    communities[map_it->second].push_back(map_it->first);
	  }

	int num_labels = communities.size();
	// Initialize the node_to_pagerank with empty vector, eventually
	// this will contain the results of the personalized page rank 
	// started at that node
	for (map_it = cored_labels.begin(); map_it != cored_labels.end(); ++map_it)
	  {
	    vector<double> dummy(num_labels, 0);
	    node_to_pagerank[map_it->first] = dummy;
	  }
  
	// Pick random sample.
	// cout<<"Get sample.\n";
	srand(time(NULL));
	map<int, vector<v_size_t> >::iterator it;
	vector<v_size_t> sample_nodes;
	vector<v_size_t> sample_labels;
	int count = 0;
	for (it = communities.begin(); it != communities.end(); ++it)
	  {
	    if (it->first >= 0)
	      {
		if (!suppressOutput)
		  cout<<++count<<"\n";
		vector<v_size_t> vert_ids = it->second;
		size_t comm_size = vert_ids.size();
		size_t num_samples = sample * comm_size;
		if (num_samples == 0)
		  num_samples = 1;

		vector<v_size_t> comm_sample;
		while (comm_sample.size() < num_samples)
		  {
		    int new_sample = rand() % vert_ids.size();
	      
		    bool unique_flag = true;
		    for (size_t i = 0; i < comm_sample.size(); ++i)
		      if (comm_sample[i] == vert_ids[new_sample])
			{
			  unique_flag = false;
			  break;
			}

		    if (unique_flag)
		      comm_sample.push_back(vert_ids[new_sample]);
		  }

		// TODO need to adjust for the "labeled" vertices;

		for (size_t i = 0; i < comm_sample.size(); ++i)
		  {
		    sample_nodes.push_back(comm_sample[i]);
		    sample_labels.push_back(cored_labels[comm_sample[i]]);
		    cored_labels[comm_sample[i]] = -2;
		  }
	      }
	  }
      
	// Now run a diffusion at each node.  Need structures for holding
	// results.  Precision, recall, etc.
	// cout<<"Starting diffusions\n";
	for (size_t i = 0; i < sample_nodes.size(); ++i)
	  {
	    if (!suppressOutput) 
	      cout<<"Node "<<sample_nodes[i]<<" diffusion.\n";
	    Vert seed = vertex(sample_nodes[i], G_cored);
	    vector<Vert> touched_nodes;
	    vector<double> page_rank;
	    spectral_snapshot(G_cored, touched_nodes, page_rank, seed, alpha, diffusion_size);
	    if (!suppressOutput)
	      cout<<"Finished snapshot, "<<touched_nodes.size()<<" nodes touched.\n";
	  
	    for (size_t j = 0; j < touched_nodes.size(); ++j)
	      {
		Vert touched = touched_nodes[j];
		v_size_t v_ind = cored_index[touched];
		// cout<<v_ind<<", ";
		// cout.flush();
		int v_label = sample_labels[i];
		node_to_pagerank[v_ind][v_label] += page_rank[v_ind];
	      }
	    if (!suppressOutput)
	      cout<<"\nFinished finding touched nodes\n";
	  }

	map<v_size_t, vector<double> >::iterator pr_it;
	vector<size_t> label_num_correct(num_labels, 0);
	vector<size_t> label_num_predict(num_labels, 0);

	for (pr_it = node_to_pagerank.begin(); pr_it != node_to_pagerank.end(); ++pr_it)
	  {
	    if (!suppressOutput)
	      cout<<"Node "<<pr_it->first<<" label: "<<cored_labels[pr_it->first]<<"\n";
	    vector<double> ranks = pr_it->second;
	    double max_label_rank = vector_max(ranks);
	    if (max_label_rank > 0)
	      {
		size_t ind = -1;
		for (size_t i = 0; i < ranks.size(); ++i)
		  {
		    if (ranks[i] > 0  && !suppressOutput)
		      cout<<"("<<i<<","<<ranks[i]<<"), ";
		    if (ranks[i] == max_label_rank)
		      ind = i;
		  }
		if (!suppressOutput)
		  cout<<"\nPredicted label "<<ind<<"\n";
		++label_num_predict[ind];
		if (ind == cored_labels[pr_it->first])
		  ++label_num_correct[ind];
	      }
	  }

	vector<double> precision(num_labels, 0);
	vector<double> recall(num_labels, 0);

	if (avg_precision.size() == 0.0) {
	  avg_precision.assign(num_labels, 0.0);
	  avg_recall.assign(num_labels, 0.0);
	}

	for (size_t i = 0; i < precision.size(); ++i)
	  {
	    if (label_num_predict[i] > 0)
	      precision[i] = (double) label_num_correct[i] / (double) label_num_predict[i];

	    if (communities[i].size() > 0)
	      recall[i] = (double) label_num_correct[i] / (double) communities[i].size();
	  
	    if (cluster_label == -1 || cluster_label == integer_to_label[i]) 
	      {
		avg_precision[max_core + 1 - c] += precision[i];
		avg_recall[max_core + 1 - c] += recall[i];
		// cout<<"#Max Core "<<c<<"\n#Label: "<<integer_to_label[i]<<"\n#Precision: "<<precision[i]<<"\n#Recall: "<<recall[i]<<"\nF1 Score: "<<f1_score(precision[i], recall[i])<<"\nPredicted: "<<label_num_predict[i]<<" Total nodes with this label: "<<communities[i].size()<<"\n";
		// cout<<integer_to_label[i]<<"\t"<<c<<"\t"<<precision[i]<<"\t"<<recall[i]<<"\t"<<f1_score(precision[i], recall[i])<<"\t"<<label_num_predict[i]<<"\t"<<communities[i].size()<<"\n";
	      }
	  }
      }
  }
  for (size_t i = 0; i < avg_precision.size(); ++i) {
    if (num_iterations != 0 && avg_recall[i] != 0 && avg_precision[i] != 0) {
      cout<<cluster_label<<"\t"<<i<<"\t"<<avg_precision[i] / num_iterations<<"\t"<<avg_recall[i] / num_iterations<<"\t"<<f1_score(avg_precision[i] / num_iterations, avg_recall[i] / num_iterations)<<"\n";
    }
  }

  return 0;
}

/* Compute f1 score given precision / recall */
double f1_score(double p, double r)
{
  return 2.0 * (p * r) / (p + r);
}

/* 
 Reads in labels from file.  Currently does not check vertex indices
 to see if they are valid.  Will throw error if they are larger than
 the size of is_labeled.  Assumes file of following format:

 map meanings:
 labels: a map from vertex_id to an integer, sequential label
 label_to_integer: a map from the original integer label to sequential integer label
 integer_to_label: a map from the sequential integer label to the original integer label

 vertex_id vertex_label
 vertex_id vertex_label
 ...
*/
void read_label_file(string input_label_file, map<v_size_t, int> & labels, map<v_size_t, bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> & integer_to_label)
{
    ifstream label_file(input_label_file.c_str());

    if (label_file.is_open())
    {
      string line;
      int line_number = 0;
      int labels_so_far = 0;
      while (getline(label_file, line))
	{
	  line_number++;
	  //	  cout<<line_number<<"\n";
	  // Comment character is #
	  if (line.substr(0,1) != "#")
	    {
	      istringstream line_stream(line);
	      int vertex_id, vertex_label;

	      if (line_stream>>vertex_id>>vertex_label)
		{
		  vertex_id = vertex_id - 1;
		  is_labeled[vertex_id] = true;
		  if (label_to_integer.find(vertex_label) == label_to_integer.end())
		    {
		      label_to_integer[vertex_label] = labels_so_far;
		      integer_to_label[labels_so_far] = vertex_label;
		      labels[vertex_id] = label_to_integer[vertex_label];
		      labels_so_far++;
		    }
		  else
		    labels[vertex_id] = label_to_integer[vertex_label]; 
		}
	    }
	}
    }
}

// Returns periphery of graph.  The prevIndices vertex property will
// contain the original indices of each vertex, i.e. prevIndices[0]
// should have the original index, which can be used for finding
// labels.
Graph get_labeled_cores(Graph & G, int max_core) 
{
  vector<int> core = k_core(G);

  // First find periphery (index off-by-one??)
  vector<Vert> periphery;
  for (int j = 0; j < core.size(); ++j)
    if (core[j] <= max_core)
      periphery.push_back(vertex(j, G));
      
  Graph G_periph = subset(G, periphery);
  return connected(G_periph);
}

// Returns core of graph.  The prevIndices vertex property will
// contain the original indices of each vertex, i.e. prevIndices[0]
// should have the original index, which can be used for finding
// labels.
Graph get_labeled_periphery(Graph & G, int min_core) 
{
  vector<int> core = k_core(G);

  // First find core (index off-by-one??)
  vector<Vert> graph_core;
  for (int j = 0; j < core.size(); ++j)
    if (core[j] > min_core)
      graph_core.push_back(vertex(j, G));
      
  Graph G_core = subset(G, graph_core);
  return connected(G_core);
}

Graph get_labeled_tree_cores(Graph & G, TreeDecomp & T, int max_core) 
{
  vector<int> core = tree_core(G, T);

  vector<Vert> periphery;
  for (int j = 0; j < core.size(); ++j)
    if (core[j] <= max_core)
      periphery.push_back(vertex(j, G));
  
  Graph G_periph = subset(G, periphery);
  return connected(G_periph);
}

Graph get_labeled_tree_d_cores(Graph & G, int max_d)
{
  vector<int> core = tree_d_core(G);

  vector<Vert> periphery;
  for (int j = 0; j < core.size(); ++j)
    if (core[j] <= max_d)
      periphery.push_back(vertex(j, G));
  
  Graph G_periph = subset(G, periphery);
  return connected(G_periph);
}

vector<int> tree_core(Graph & G, TreeDecomp & T)
{
  vector< vector<v_size_t> > vertex_locations = T.get_vertex_locations();

  vector<int> core(vertex_locations.size(), 0);
  for (int i = 0; i < core.size(); ++i) 
    {
      core[i] = (int) vertex_locations[i].size();
    }

  return core;
}

vector<int> tree_d_core(Graph G) 
{
  cout<<"Priority queue...";
  priority_queue<pair<Vert, int>, vector<pair<Vert,int> >, compare_second_greater> q;
  cout<<"done.\n";
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  cout<<"Initialize...";
  graph_traits<Graph>::vertex_iterator vi, vie;
  for (tie(vi, vie) = vertices(G); vi != vie; ++vi) 
    q.push(make_pair(*vi, out_degree(*vi, G)));
  cout<<"done.\n";
  
  vector<int> core(num_vertices(G), 51);
  while (!q.empty())
    {
      pair<Vert, int> v_deg_pair = q.top();
      q.pop();
      cout<<"Vert: "<<v_deg_pair.first<<" deg: "<<v_deg_pair.second<<"\n";

      if (v_deg_pair.second == out_degree(v_deg_pair.first, G))
	{
	  size_t ind = index[v_deg_pair.first];
	  core[ind] = v_deg_pair.second;
	  graph_traits<Graph>::adjacency_iterator ai, aie;
	  vector<Vert> neighbors;
	  for (tie(ai, aie) = adjacent_vertices(v_deg_pair.first, G); ai != aie; ++ai)
	    neighbors.push_back(*ai);
	      	  
	  for (size_t i = 0; i < neighbors.size(); ++i) 
	    for (size_t j = i + 1; j < neighbors.size(); ++j)
	      {
		Vert u = neighbors[i];
		Vert v = neighbors[j];
		if (u != v)
		  {
		    pair<Edge, bool> edge_exists = edge(u, v, G);
		    if (!edge_exists.second)
		      {
			add_edge(u, v, G);
			add_edge(v, u, G);
		      }
		  }
	      }
	  clear_vertex(v_deg_pair.first, G);
	  remove_vertex(v_deg_pair.first, G);
	}
      else
	q.push(make_pair(v_deg_pair.first, out_degree(v_deg_pair.first, G)));
    }

  cout<<"Done.\n";
  return core;
}

// Takes the basic label vector and converts it to the subgraph
map<v_size_t, int> convert_labels(Graph & G, map<v_size_t, int> labels)
{
  map<v_size_t, int> new_labels;

  v_size_t n = num_vertices(G);

  // INDEX?!? OFF_BY_ONE?>?!?
  for (v_size_t i = 0; i < n; ++i)
    {
      new_labels[i] = labels[G[vertex(i, G)].prevIndices[0]];
    }

  return new_labels;
}

// Finds max of a vector of value types
template <class T> 
T vector_max(vector<T> values)
{
  T max = values[0];
  for (size_t i = 1; i < values.size(); ++i)
    if (values[i] > max)
      max = values[i];

  return max;
}
