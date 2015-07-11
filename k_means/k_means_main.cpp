/*
  This file generates an executable which takes a community file (each
  node is assigned an integer community), uses k-means to find 
  
  By Aaron Adcock, PhD (Finishing research performed at Stanford University)
  Oct. 2014
*/

#include "graph_lib_boost.hpp"
#include "tree_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
#include <sstream>
#include <string.h>
#include <cmath>
#include <utility>

template <class Type> 
bool vector_find(Type key, vector<Type> values);

double k_means_step(Graph& G, vector<size_t> & clusters, vector< vector<double> > & centroids);
size_t total_in_community(int comm, vector<size_t> & nodes, vector<int> & communities);

vector<size_t> network_k_means(Graph & G, size_t k);
vector<size_t> size_of_component(vector<size_t> & components);
vector<size_t> convert_clusters_to_vertices(vector<size_t> & vert_to_clusters, vector<size_t> & clusters);

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
  size_t k = 0;

  if(argc<2)
    {
      cout<<"The options are (-g -community are required): \n-g <input_graph_file> \n-o <output_file_name_prefix> \n-community <community_file>\n-num_clusters <default is number of communities in community file>";
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
      else if (strcmp(argv[i], "-num_clusters")== 0)
	{
	  ss<<argv[i+1];
	  ss>>k;
	}
    }
  ////

  if(graph_file.length() == 0 || community_file.length() == 0)
    {
      cout<<"You must provide an input graph/tree/community file\n";
      return(0);
    }
  else
    cout<<"Graph: "<<graph_file<<"\nTree: "<<tree_file<<"\nCommunity: "<<community_file<<"\n";

  cout.flush();
 
  // Load Graph
  Graph G = loadGraph(graph_file, "\t", false);
  v_size_t size = num_vertices(G);
  cout<<"Connected Network has "<<size<<" vertices.\n";

  // Load Communities
  vector<int> communities;
  load_color_file(community_file, communities);
  cout<<"Color file loaded.\n";

  if (communities.size() != size)
    {
      cerr<<"Community vector length: "<<communities.size()<<"\n";
      cerr<<"Graph size: "<<size<<"\n";
      cerr<<"Community vector not the right size. Exiting. \n";
      exit(EXIT_FAILURE);
    }

  // Get sorted list of community labels.  I assume one is the
  // "unknown" label, in my data this is typically 0.
  vector<int> community_list = communities;

  sort(community_list.begin(), community_list.end());
  vector<int>::iterator it = unique(community_list.begin(), community_list.end());
  community_list.resize(distance(community_list.begin(), it));
 
  // Map from communities to relative frequency in network
  map<int, double> community_frequency;
  map<int, size_t> community_sizes;

  cout<<"Community frequencies.\n";
  for (size_t i = 0; i < communities.size(); ++i)
    {
      if (community_frequency.find(communities[i]) != community_frequency.end())
	{
	  community_frequency[communities[i]] += 1.0 / communities.size();
	  community_sizes[communities[i]] += 1;
	}
      else
	{    
	  community_frequency[communities[i]] = 1.0 / communities.size();
	  community_sizes[communities[i]] = 1;
	}
    }

  // Three ways to assign clusters to communities in a
  // meaningful way.  First, the community forms a plurality of nodes.
  // Second, the community forms a majority of nodes.  Third, the
  // community is more frequent in the cluster than in the network.  The
  // first two definitions will typically place a bag in a unique
  // location excepting ties (there may be many ties in the small
  // bags).  The frequency definition could potentially place any bag
  // in multiple communities.
  vector<size_t> dummy;
  map<int, vector<size_t> > cluster_plural_community;
  map<int, vector<size_t> > cluster_major_community;
  map<int, vector<size_t> > cluster_unlikely_community;
  
  // Initialize maps
  for (size_t i = 0; i < community_list.size(); ++i)
    {
      cluster_plural_community[community_list[i]] = dummy;
      cluster_major_community[community_list[i]] = dummy;
      cluster_unlikely_community[community_list[i]] = dummy;
    }

  /* Run through each label and find best cluster */
  // Compute clusters 
  if (k <= 0)
    k = community_list.size();

  vector<size_t> vert_to_cluster = network_k_means(G, k);

  cout<<"Creating cluster data structures\n";
  vector< vector<size_t> > clusters;
  for (size_t i = 0; i < k; ++i)
    {
      vector<size_t> dummy;
      clusters.push_back(dummy);
    }
 
  for (size_t i = 0; i < vert_to_cluster.size(); ++i)
    {
      cout<<"Vert "<<i<<" is in cluster "<<vert_to_cluster[i]<<"\n";
      size_t clust = vert_to_cluster[i];
      clusters[clust].push_back(i);
    }

  for (size_t i = 0; i < clusters.size(); ++i)
    cout<<"Cluster "<<i<<" of size "<<clusters[i].size()<<"\n";

  cout<<"Calculating community-cluster frequencies\n";
  for (size_t i = 0; i < clusters.size(); ++i)
    {
      map<int, double> community_cluster_frequency;
      vector<size_t> current_cluster = clusters[i];
      for (size_t j = 0; j < current_cluster.size(); ++j)
	{
	  // Find relative frequency of each community in the cluster
	  size_t vert = current_cluster[j];
	  int comm = communities[vert];
	  if (community_cluster_frequency.find(comm) == community_cluster_frequency.end())
	    community_cluster_frequency[comm] = 0;
	  community_cluster_frequency[comm] += 1.0 / current_cluster.size();
	}

      cout<<"Associating with cluster to community\n";
       // Need to use each test type
      map<int, double>::iterator mit;
      vector<int> max_communities;
      double max_frequency = 0;
      for (mit = community_cluster_frequency.begin(); mit != community_cluster_frequency.end(); ++mit)
	{
	  cout<<"Community: "<<mit->first<<" frequency :"<<mit->second<<"\n";
	  if (mit->second > max_frequency)
	    {
	      max_communities.clear();
	      max_communities.push_back(mit->first);
	      max_frequency = mit->second;
	    }
	  else if (mit->second == max_frequency)
	    max_communities.push_back(mit->first);
	}

      // Plural
      for (size_t j = 0; j < max_communities.size(); ++j)
	cluster_plural_community[max_communities[j]].push_back(i);

      // Major
      if (max_frequency >= 0.5)
	for (size_t j = 0; j < max_communities.size(); ++j)
	  cluster_major_community[max_communities[j]].push_back(i);

      // Likely
      for (mit = community_frequency.begin(); mit != community_frequency.end(); ++mit)
	{
	  if (mit->second <= community_cluster_frequency[mit->first])
	    cluster_unlikely_community[mit->first].push_back(i);
	}
    }

  cout<<"Calculating stats\n";
  // Should have a map with the "best' set of clusters for each
  // community given the different definitions of best.  Unlike with
  // TD, here we just need to calculate the stats (recall precision)
  map<int, pair<double, double> > plural_stats;
  map<int, pair<double, double> > major_stats;
  map<int, pair<double, double> > unlikely_stats;
  for (size_t i = 0; i < community_list.size(); ++i)
    {
      // First convert into vertices (should be disjoint, we won't
      // check for repeats
      int curr_comm = community_list[i];
      vector<size_t> plural_nodes = convert_clusters_to_vertices(vert_to_cluster, cluster_plural_community[curr_comm]);
      vector<size_t> major_nodes = convert_clusters_to_vertices(vert_to_cluster, cluster_major_community[curr_comm]);
      vector<size_t> likely_nodes = convert_clusters_to_vertices(vert_to_cluster, cluster_unlikely_community[curr_comm]);

      size_t plural_total = total_in_community(curr_comm, plural_nodes, communities);
      size_t major_total = total_in_community(curr_comm, major_nodes, communities);
      size_t likely_total = total_in_community(curr_comm, likely_nodes, communities);

      if (plural_total != 0)
	plural_stats[curr_comm] = make_pair((double) plural_total / (double) community_sizes[curr_comm], (double) plural_total / (double) plural_nodes.size());
      else
	plural_stats[curr_comm] = make_pair(0.0, 0.0);

      if (major_total != 0)
	major_stats[curr_comm] = make_pair((double) major_total / (double) community_sizes[curr_comm], (double) major_total / (double) major_nodes.size());
      else
	major_stats[curr_comm] = make_pair(0.0, 0.0);

      if (likely_total != 0)
	unlikely_stats[curr_comm] = make_pair((double) likely_total / (double) community_sizes[curr_comm], (double) likely_total / (double) likely_nodes.size());
      else
	unlikely_stats[curr_comm] = make_pair(0.0, 0.0);
    }

  cout<<"Stats for the communities:\n";
  for (size_t i = 0; i < community_list.size(); ++i)
    {
      int curr_comm = community_list[i];
      cout<<"Community "<<curr_comm<<", size = "<<community_sizes[curr_comm]<<"\n";
      cout<<"Plural recall: "<<plural_stats[curr_comm].first<<" precision: "<<plural_stats[curr_comm].second<<"\n";
      cout<<"Major recall: "<<major_stats[curr_comm].first<<" precision: "<<major_stats[curr_comm].second<<"\n";
      cout<<"Unlikely recall: "<<unlikely_stats[curr_comm].first<<" precision: "<<unlikely_stats[curr_comm].second<<"\n";
    }
}

vector<size_t> convert_clusters_to_vertices(vector<size_t> & vert_to_clusters, vector<size_t> & clusters)
{
  vector<size_t> nodes;
  for (size_t i = 0; i < vert_to_clusters.size(); ++i)
    if (vector_find(vert_to_clusters[i], clusters))
      nodes.push_back(i);

  return nodes;
}

// Returns (recall,precision) of community in given cluster
size_t total_in_community(int comm, vector<size_t> & nodes, vector<int> & communities)
{
  size_t total = 0;
  for (size_t i = 0; i < nodes.size(); ++i)
    if (communities[nodes[i]] == comm)
      total++;

  return total;
}

// Finds max of a vector of value types
template <class Type> 
bool vector_find(Type key, vector<Type> values)
{
  for (size_t i = 0; i < values.size(); ++i)
    if (values[i] == key)
      return true;

  return false;
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

/*
  Takes a graph and finds k-clusters using k-means (Lloyd's
  algorithm).  Vectors are produced using the adjacency matrix.  The
  initial conditions are found by randomly selecting k nodes from the
  networks.  Proceeds until convergence or max_iter is reached.

 */
vector<size_t> network_k_means(Graph & G, size_t k)
{
  // Constants
  const double CONVERGED = 0.001;
  const size_t MAX_ITER = 100;

  // First pick starting nodes
  size_t n = num_vertices(G);
  
  srand(time(NULL));

  // Pick vertex indices to start with (random)
  vector<Vert> start_nodes;
  start_nodes.reserve(k);
  for (size_t i = 0; i < k; ++i)
    {
      v_size_t ind = rand() % n;
      start_nodes.push_back(vertex(ind, G));
    }

  // Convert nodes into adjacency vectors
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  
  vector< vector<double> > centroids;
  for (size_t i = 0; i < k; ++i)
    {
      Vert node = start_nodes[i];
      vector<double> adj_vector(n, 0);
      graph_traits<Graph>::adjacency_iterator ai, aie;
      for (tie(ai, aie) = adjacent_vertices(node, G); ai != aie; ++ai)
	{
	  // Convert node object into index
	  v_size_t ind = index[*ai];
	  adj_vector[ind] = 1.0;
	}
      centroids.push_back(adj_vector);
    }
  // Iterations: 1. Associate each node to a centroid
  //             2. Update centroids
  double eps = 2.0;
  size_t num_iter = 0;
  vector<size_t> clusters(n, 0);
  while (eps > CONVERGED && num_iter < MAX_ITER)
    {
      num_iter++;
      cout<<"On iteration "<<num_iter<<"\n";
      eps = k_means_step(G, clusters, centroids);
      cout<<eps<<"\n";
    }

  return clusters;
}

/*
  Runs a single step of the k-means algorithm.
 **/
double k_means_step(Graph& G, vector<size_t> & clusters, vector< vector<double> > & centroids)
{
  // Basic values
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  double eps = 0;
  v_size_t n = num_vertices(G);
  vector<int> empty_flag(centroids.size(), 1);

  // First update clusters
  for (size_t i = 0; i < clusters.size(); ++i)
    {
      double min_distance = -1.0;
      int min_cluster = -1;
      Vert node = vertex(i, G);

      for (size_t j = 0; j < centroids.size(); ++j)
	{
	  vector<double> centroid = centroids[j];
	  
	  double sum = 0.0;
	  graph_traits<Graph>::adjacency_iterator ai, aie;
	  for (tie(ai, aie) = adjacent_vertices(node, G); ai != aie; ++ai)
	    {
	      size_t ind = index[*ai];
	      
	      sum += (centroid[ind] - 1) * (centroid[ind] - 1);
	    }
	  
	  if (j == 0 || sum < min_distance)
	    {
	      min_distance = sum;
	      min_cluster = j;
	    }
	}
      if (min_cluster >= 0)
	{
	  clusters[i] = min_cluster;
	  empty_flag[min_cluster] = 0;
	}
      else
	cerr<<"Error.  No cluster chosen\n";
    }

  for (size_t i = 0; i < empty_flag.size(); ++i)
    {
      if (empty_flag[i] == 1)
	{
	  size_t rand_vert = rand() % clusters.size();
	  clusters[rand_vert] = i;
	}
    }

  // Clusters updated, Moving on to centroids
  vector<double> dummy(n, 0);
  vector< vector<double> > new_centroids(centroids.size(), dummy);
  vector<double> cluster_sizes(centroids.size(), 0);
  for (size_t i = 0; i < clusters.size(); ++i)
    {
      size_t cluster = clusters[i];
      cluster_sizes[cluster] += 1;
      graph_traits<Graph>::adjacency_iterator ai, aie;
      Vert node = vertex(i, G);
      for (tie(ai, aie) = adjacent_vertices(node, G); ai != aie; ++ai)
	{
	  size_t ind = index[*ai];
	  new_centroids[cluster][ind] += 1.0;
	}
    }

  for (size_t i = 0; i < new_centroids.size(); ++i)
    {
      double tmp_eps = 0;
      for (size_t j = 0; j < new_centroids[i].size(); ++j)
	{
	  new_centroids[i][j] /= cluster_sizes[i]; 
	  tmp_eps += pow((centroids[i][j] - new_centroids[i][j]), 2);
	  centroids[i][j] = new_centroids[i][j];
	}
      eps += sqrt(tmp_eps);
    }

  return eps / (sqrt(n) * centroids.size());
}
