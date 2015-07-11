/*
This function aims to take snapshots of a graph based on the k-core
decomposition.  One option is to use a local diffusion in different
cores of the network. The local diffusion used is the Andersen, Chung,
Lu method.

By Aaron Adcock, PhD Candidate, Stanford University (~2011)
modified May, 2014
*/

#include "../graph_lib_boost.hpp"
#include <cstdlib>
#include <ctime>
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <boost/graph/copy.hpp>

/*
  This function runs an approximate personalized PageRank using the
  provided seed node and the provided parameters (alpha, eps).
  Returned is the page_rank vector (loosely, the approximate
  probability of visiting each node in the network) and the set of
  nodes which were touched (snapshot).  These vectors are intended for
  use in the write_graphviz(~) function.
 */
void spectral_snapshot(Graph& G, vector<Vert>& snapshot, vector<double>& page_rank, Vert seed_node, const double alpha, const double eps)
{
  // Vertex index map
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  
  // Basic values
  long size = num_vertices(G);
  long r = index[seed_node];
  graph_traits<Graph>::vertex_iterator vi, vi_end;
  if (size <= r || r < 0) {
    cout<<"Seed Node: "<<seed_node<<" r: "<<r<<" size: "<<size<<"\n";
    for (tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi) 
      cout<<index[*vi]<<" ";
    cout.flush();
  }
   
  // Probability (PageRank!) vector
  page_rank.assign(size, 0);
  // Seed vector
  vector<double> seed(size, 0);
  snapshot.clear();
  seed[r] = 1;
  // Personalized Page Rank using Andersen, Chung, Lu method
  approxPR(G, seed, alpha, eps, page_rank, r);
  // Find non-zero nodes (i.e. nodes which the diffusion reaches)
  vector<double> new_page_rank;
  for (tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi)
    {
      Vert v = *vi;
      int ind = index[v];

      if (page_rank[ind] != 0)
	snapshot.push_back(v);
    }
}

/*
  This function uses the spectral snapshot function to run a diffusion
  on every node with a particular integer value assigned to it, in a
  network G.  It writes out the resulting subgraph to graph_viz files.
  The intent is to label the nodes using the core number (i.e. looking
  at k-shells in the network), but this method is actually more
  generally useful.  The scaling/coloring are chosen using the k-core
  value and page rank value

  TODO: Get some sort of stats on the labelling

  INPUT: Graph G
         
 */

void labeled_snapshot(Graph& G, vector<int>& node_labels, int label, string filename, const double alpha, const double eps)
{
  // Make sure that node_labels and size of graph is the same.
  long size = num_vertices(G);
  if (size > node_labels.size())
    {
      cerr<<"Not enough node labels were provided. Exiting\n";
      exit(EXIT_FAILURE);
    }

  // Vectors and diffusion parameters
  vector<double> page_rank(size, 0);
  vector<Vert> touched_nodes;

  for (int i = 0; i < size; ++i)
    {
      if (node_labels[i] == label)
	{
	  Vert seed_node = vertex(i, G);
	  spectral_snapshot(G, touched_nodes, page_rank, seed_node, alpha, eps); 
	  
	  // Writeout a graphviz file.  Need to update filename.
	  stringstream convert_node;
	  convert_node << i;
	  stringstream convert_label;
	  convert_label << label;
	  string full_filename = filename;
	  full_filename.append("_label_");
	  full_filename.append(convert_label.str());
	  full_filename.append("_node_");
	  full_filename.append(convert_node.str());
	  full_filename.append(".dot");
	  vector<size_t> node_sizing;
	  
	  int min_label = node_labels[0];
	  for (int j = 0; j < node_labels.size(); ++j)
	    {
	      if (min_label > node_labels[j])
		min_label = node_labels[j];
	    }

	  // Need to scale colors for this function
	  double max_page_rank = page_rank[0];
	  double min_page_rank = page_rank[0];
	  vector<double> coloring;
	  int num_nodes_in_subgraph = 0;
	  vector<int> seed_node_ind;
	  for (int j = 0; j < page_rank.size(); ++j)
	    {
	      if (max_page_rank < page_rank[j])
		max_page_rank = page_rank[j];
	      
	      if (min_page_rank > page_rank[j])
		min_page_rank = page_rank[j];

	      if (page_rank[j] != 0)
		{
		  if (j == i)
		    seed_node_ind.push_back(num_nodes_in_subgraph);
		  
		  node_sizing.push_back(node_labels[j] - min_label);
		  coloring.push_back(page_rank[j]);
		  num_nodes_in_subgraph++;
		}
	    }
	  
	  for (int j = 0; j < coloring.size(); ++j)
	    coloring[j] = (coloring[j] - min_page_rank) / (max_page_rank - min_page_rank);

	  Graph diffusion_neighborhood = subset(G, touched_nodes);

	  write_scaled_square_graphviz(diffusion_neighborhood, full_filename, coloring, node_sizing,seed_node_ind, false);
	}
    }
}


/*
  This function uses the spectral snapshot function to run a diffusion
  on every node with a particular integer value assigned to it in the
  network G.  It then returns a vector of statistics showing the
  labels of each node touched by the diffusion.  It is called
  community snapshot because the original intent of the function is to
  use labels which represent communities in the network (i.e. sets of
  nodes with more connections within than without).  However, it may
  be used for any set of labelings (such as the k-core or degree of
  the node).

  INPUT: Graph G 
  OUTPUT: A vector of statistics across the labels.
         
 */

vector<double> community_snapshot(Graph& G, vector<int>& node_labels, int label, const double alpha, const double eps)
{
  // Make sure that node_labels and size of graph is the same.
  long size = num_vertices(G);
  if (size > node_labels.size())
    {
      cerr<<"Not enough node labels were provided. Exiting\n";
      exit(EXIT_FAILURE);
    }

  // Vectors and diffusion parameters
  vector<double> page_rank(size, 0);
  vector<Vert> touched_nodes;

  // Vertex index map
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  map<int, int> label_to_index;
  int count = 0;
  for (size_t i = 0; i < size; ++i)
    {
      if (label_to_index.find(node_labels[i]) == label_to_index.end())
	{
	  label_to_index[node_labels[i]] = count;
	  //ind = count;
	  count++;
	}
    }

  vector<double> distribution(count, 0);
  for (int i = 0; i < size; ++i)
    { 
      if (node_labels[i] == label)
	{
	  Vert seed_node = vertex(i, G);
	  spectral_snapshot(G, touched_nodes, page_rank, seed_node, alpha, eps); 
	  
	  // add +1 to correct distribution slot.  Consider writing single touch version
	  for (size_t j = 0; j < touched_nodes.size(); ++j)
	    {
	      Vert v = touched_nodes[j];
	      size_t v_ind = index[v];
	      int v_label = node_labels[v_ind];

	      ++distribution[label_to_index[v_label]];
	    }
	}
    }

  return distribution;
}

/*
  This function uses the spectral snapshot function to run a diffusion
  on every node with a particular integer value assigned to it in the
  network G.  It then returns a vector of statistics showing the
  labels of each node touched by the diffusion.  It is called
  community snapshot because the original intent of the function is to
  use labels which represent communities in the network (i.e. sets of
  nodes with more connections within than without).  However, it may
  be used for any set of labelings (such as the k-core or degree of
  the node).

  INPUT: Graph G 
  OUTPUT: A vector of statistics across the labels.
         
 */

map<int, double> community_snapshot(Graph& G, map<v_size_t, int> & node_labels, int label, const double alpha, const double eps)
{
  // Make sure that node_labels and size of graph is the same.
  long size = num_vertices(G);
  
  // Vectors and diffusion parameters
  vector<double> page_rank(size, 0);
  vector<Vert> touched_nodes;

  // Vertex index map
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  map<int, double> distribution;
  map<v_size_t, int>::iterator map_it;
  for (map_it = node_labels.begin(); map_it != node_labels.end(); ++map_it)
    { 
      if (map_it->second == label)
	{
	  Vert seed_node = vertex(map_it->first - 1, G);
	  spectral_snapshot(G, touched_nodes, page_rank, seed_node, alpha, eps); 
	  
	  // add +1 to correct distribution slot.  Consider writing single touch version
	  for (size_t j = 0; j < touched_nodes.size(); ++j)
	    {
	      Vert v = touched_nodes[j];
	      v_size_t v_ind = index[v];
	      int v_label = node_labels[v_ind+1];
	      
	      if (distribution.find(v_label) == distribution.end())
		{
		  distribution[v_label] = 1.0;
		}
	      else
		{
		  ++distribution[v_label];
		}
	    }
	}
    }

  return distribution;
}
