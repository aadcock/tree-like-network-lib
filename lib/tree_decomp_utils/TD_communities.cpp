/* 
   The purpose of the functions in this file are to generate candidate
   communities using a tree decomposition.  The theory behind this is
   that since tree decompositions of low width provide small graph
   seperators, then they must be able to find good communities or at
   the very least, important bottlenecks in the network.

   By Aaron Adcock, PhD Candidate at Stanford University
   Feb 2014 
*/

#include "../graph_lib_boost.hpp"
#include "../tree_lib_boost.hpp"
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>

double compute_conductance(vector<v_size_t> & list_vertices, Graph & G);
double compute_conductance(set<Vert> & verts, Graph & G);

// Function for finding minimum conductance communities across various
// size scales.  This function just keeps track of the minimum
// conductance, not the actual communities and produces an NCP like
// output.  The communities are built up using a provided tree
// decomposition.  We do this by calculating the eccentricity of each
// bag in the tree.  We then progressively build communities by
// growing communities on the tree (ie, connected components of trees
// by growing the community from the leaves in the tree up to the
// center)

/*
  INPUT:
  G:                   Graph under investigation
  td:                  Associated tree decomposition
  max_community:       Maximum community size to find 

  OUTPUT:
  A double vector 
*/
vector<double> td_eccentricity_communities(Graph& G, TreeDecomp& td, size_t max_community, bool large_comp_only)
{
  vector<double> conductance(max_community, 1.0);
  return td_eccentricity_communities(G, td, max_community, conductance, large_comp_only);
}

vector<double> td_eccentricity_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, bool large_comp_only)
{
  // First calculate tree eccentricty
  Graph T = td.get_tree();
  v_size_t n = num_vertices(T);

  if (n == 0 || max_community == 0)
    {
      cerr<<"Tree decomposition is empty.\n";
      vector<double> dummy;
      return dummy;
    }

  // Calculate eccentricity
  vector< v_size_t > td_eccentricity = tree_eccentricity(T);

  // Calculate range
  v_size_t max_ecc = td_eccentricity[0];
  v_size_t min_ecc = td_eccentricity[0];

  for (int i = 1; i < td_eccentricity.size(); ++i)
    {
      if (max_ecc < td_eccentricity[i])
	max_ecc = td_eccentricity[i];
      
      if (min_ecc > td_eccentricity[i])
	min_ecc = td_eccentricity[i];
    }

  cout<<"Max: "<<max_ecc<<" Min: "<<min_ecc<<"\n";

  // Begin making branch communities
  for (size_t ecc = max_ecc; ecc > min_ecc; --ecc)
    {
      vector<Vert> valid_bags;

      for (size_t i = 0; i < n; ++i)
	if (td_eccentricity[i] >= ecc)
	  valid_bags.push_back(vertex(i, T));

      // Get tree components
      Graph subtrees = subset(T, valid_bags);
      vector<v_size_t> component = components(subtrees);
      
      // Get number of components
      size_t num_comps = 1;
      for (size_t i = 0; i < component.size(); ++i)
	if ((component[i] + 1) > num_comps)
	  num_comps = component[i] + 1;
      
      for (size_t i = 0; i < num_comps; ++i)
	{
	  vector<v_size_t> vertices_in_subgraph;
	  for (size_t j = 0; j < component.size(); ++j)
	    {
	      size_t index_in_td = subtrees[j].prevIndices.back();
	      if (component[j] == i && td_eccentricity[index_in_td] != ecc)
		{
		  vector<v_size_t> nodes = td.get_bag(index_in_td);
		  vertices_in_subgraph.insert(vertices_in_subgraph.end(), nodes.begin(), nodes.end());
		}
	    }

	  sort(vertices_in_subgraph.begin(), vertices_in_subgraph.end());
	  vector<v_size_t>::iterator it = unique(vertices_in_subgraph.begin(), vertices_in_subgraph.end());
	  vertices_in_subgraph.resize(distance(vertices_in_subgraph.begin(), it));

	  if (large_comp_only)
	    {
	      Graph subgraph = subset(G, vertices_in_subgraph);
	      vector<Vert> big_comp = list_giant_component(subgraph);

	      vertices_in_subgraph.clear();
	      for(size_t big_it = 0; big_it < big_comp.size(); ++big_it)
		vertices_in_subgraph.push_back(subgraph[big_comp[big_it]].prevIndices.back());
	    }
	  
	  if (vertices_in_subgraph.size() - 1 < max_community && vertices_in_subgraph.size() > 0)
	    {
	      double c = compute_conductance(vertices_in_subgraph, G);
	      if (c < conductance[vertices_in_subgraph.size() - 1])
		conductance[vertices_in_subgraph.size() - 1] = c;
	    }
	}
    }
  return conductance;
}

vector<double> td_eccentricity_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, vector< vector<Vert> > & best_communities,  bool large_comp_only)
{
  // First calculate tree eccentricty
  Graph T = td.get_tree();
  v_size_t n = num_vertices(T);

  if (n == 0 || max_community == 0)
    {
      cerr<<"Tree decomposition is empty.\n";
      vector<double> dummy;
      return dummy;
    }

  if (best_communities.size() < max_community)
    best_communities.resize(max_community);

  // Calculate eccentricity
  vector< v_size_t > td_eccentricity = tree_eccentricity(T);

  // Calculate range
  v_size_t max_ecc = td_eccentricity[0];
  v_size_t min_ecc = td_eccentricity[0];

  for (int i = 1; i < td_eccentricity.size(); ++i)
    {
      if (max_ecc < td_eccentricity[i])
	max_ecc = td_eccentricity[i];
      
      if (min_ecc > td_eccentricity[i])
	min_ecc = td_eccentricity[i];
    }

  cout<<"Max: "<<max_ecc<<" Min: "<<min_ecc<<"\n";

  // Begin making branch communities
  for (size_t ecc = max_ecc; ecc > min_ecc; --ecc)
    {
      vector<Vert> valid_bags;

      for (size_t i = 0; i < n; ++i)
	if (td_eccentricity[i] >= ecc)
	  valid_bags.push_back(vertex(i, T));

      // Get tree components
      Graph subtrees = subset(T, valid_bags);
      vector<v_size_t> component = components(subtrees);
      
      // Get number of components
      size_t num_comps = 1;
      for (size_t i = 0; i < component.size(); ++i)
	if ((component[i] + 1) > num_comps)
	  num_comps = component[i] + 1;
      
      for (size_t i = 0; i < num_comps; ++i)
	{
	  vector<v_size_t> vertices_in_subgraph;
	  for (size_t j = 0; j < component.size(); ++j)
	    {
	      size_t index_in_td = subtrees[j].prevIndices.back();
	      if (component[j] == i && td_eccentricity[index_in_td] != ecc)
		{
		  vector<v_size_t> nodes = td.get_bag(index_in_td);
		  vertices_in_subgraph.insert(vertices_in_subgraph.end(), nodes.begin(), nodes.end());
		}
	    }

	  sort(vertices_in_subgraph.begin(), vertices_in_subgraph.end());
	  vector<v_size_t>::iterator it = unique(vertices_in_subgraph.begin(), vertices_in_subgraph.end());
	  vertices_in_subgraph.resize(distance(vertices_in_subgraph.begin(), it));

	  if (large_comp_only)
	    {
	      Graph subgraph = subset(G, vertices_in_subgraph);
	      vector<Vert> big_comp = list_giant_component(subgraph);

	      vertices_in_subgraph.clear();
	      for(size_t big_it = 0; big_it < big_comp.size(); ++big_it)
		vertices_in_subgraph.push_back(subgraph[big_comp[big_it]].prevIndices.back());
	    }
	  
	  if (vertices_in_subgraph.size() - 1 < max_community && vertices_in_subgraph.size() > 0)
	    {
	      double c = compute_conductance(vertices_in_subgraph, G);
	      if (c < conductance[vertices_in_subgraph.size() - 1])
		{
		  conductance[vertices_in_subgraph.size() - 1] = c;
		  vector <Vert> vert_subgraph;
		  for (size_t ind = 0; ind < vertices_in_subgraph.size(); ++ind)
		    vert_subgraph.push_back(vertex(vertices_in_subgraph[ind], G));

		  best_communities[vertices_in_subgraph.size() - 1] = vert_subgraph;
		}
	    }
	}
    }
  return conductance;
}

vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, bool large_comp_only)
{
  vector<double> conductance(max_community, 1.0);
  return td_bag_communities(G, td, max_community, conductance, large_comp_only);
}



vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, bool large_comp_only)
{
  size_t num_bags = td.get_num_bags();
  
  for (size_t i = 0; i < num_bags; ++i)
    {
      vector< v_size_t > bag = td.get_bag(i);
      
      if (large_comp_only)
	{
	  Graph subgraph = subset(G, bag);
	  vector<Vert> big_comp = list_giant_component(subgraph);

	  bag.clear();
	  for(size_t big_it = 0; big_it < big_comp.size(); ++big_it)
	    bag.push_back(subgraph[big_comp[big_it]].prevIndices.back());
	}

      if (bag.size() - 1 < max_community)
	{
	  double c = compute_conductance(bag, G);
	  if (c < conductance[bag.size() - 1])
	    conductance[bag.size() - 1] = c;
	}
    }
  return conductance;
}

vector<double> td_bag_communities(Graph& G, TreeDecomp& td, size_t max_community, vector<double> conductance, vector< vector<Vert> > & best_communities, bool large_comp_only)
{
  size_t num_bags = td.get_num_bags();
  
  for (size_t i = 0; i < num_bags; ++i)
    {
      vector< v_size_t > bag = td.get_bag(i);
      
      if (large_comp_only)
	{
	  Graph subgraph = subset(G, bag);
	  vector<Vert> big_comp = list_giant_component(subgraph);

	  bag.clear();
	  for(size_t big_it = 0; big_it < big_comp.size(); ++big_it)
	    bag.push_back(subgraph[big_comp[big_it]].prevIndices.back());
	}

      if (bag.size() - 1 < max_community)
	{
	  double c = compute_conductance(bag, G);
	  if (c < conductance[bag.size() - 1])
	    {
	      conductance[bag.size() - 1] = c;
	      vector<Vert> bag_vert;
	      for (size_t ind = 0; ind < bag.size(); ++ind)
		bag_vert.push_back(vertex(bag[ind], G));
	      
	      best_communities[bag.size() - 1] = bag_vert;
	    }
	}
    }
  return conductance;
}

void compute_community_tree_stats(vector< vector<Vert> > & communities, Graph & G, TreeDecomp & TD, vector<size_t> & num_bags, vector<size_t> & min_ecc, vector<size_t> & max_card, vector<size_t> & med_card, vector<size_t> & surface_area)
{
  // First need indexing function
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  // Intialize
  size_t num_communities = communities.size();
  num_bags.clear();
  min_ecc.clear();
  max_card.clear();
  med_card.clear();
  surface_area.clear();

  num_bags.assign(num_communities, 0);
  min_ecc.assign(num_communities, 0);
  max_card.assign(num_communities, 0);
  med_card.assign(num_communities, 0);
  surface_area.assign(num_communities, 0);

  vector<size_t> bag_card = TD.get_bag_sizes();
  vector<size_t> td_eccentricity = TD.get_tree_eccentricity();
  vector< vector<v_size_t> > vertex_in_bag = TD.get_vertex_locations();

  // Iterate over communities, calculate value for each stats vector

  for (int i = 0; i < communities.size(); ++i)
    {
      set<v_size_t> community_bag_set;
      vector<Vert> current_community = communities[i];
      vector<size_t> med_card_vector;
      size_t max_size_bag = 0;
      size_t min_ecc_bag = (size_t) - 1;
      for (int j = 0; j < current_community.size(); ++j)
	{
	  Vert curr_vert = current_community[j];
	  vector<v_size_t> vertex_bags = vertex_in_bag[index[curr_vert]];
	  community_bag_set.insert(vertex_bags.begin(), vertex_bags.end());
	  for (int k = 0; k < vertex_bags.size(); ++k)
	    {
	      med_card_vector.push_back(bag_card[vertex_bags[k]]);
	      v_size_t bag_ind = vertex_bags[k];
	      if (bag_card[bag_ind] > max_size_bag)
		max_size_bag = bag_card[bag_ind];

	      if (td_eccentricity[bag_ind] < min_ecc_bag)
		min_ecc_bag = td_eccentricity[bag_ind];	      
	    }
	}

      if (current_community.size() != 0)
	{
	  set<Vert> community_verts;
	  community_verts.insert(current_community.begin(), current_community.end());

	  graph_traits<Graph>::adjacency_iterator ait, aite;

	  for (int j = 0; j < current_community.size(); ++j)
	    for (tie(ait, aite) = adjacent_vertices(current_community[j], G); ait != aite; ++ait)
	      if (community_verts.find(*ait) == community_verts.end())
		surface_area[i]++;

	  num_bags[i] = community_bag_set.size();
	  min_ecc[i] = min_ecc_bag;
	  max_card[i] = max_size_bag;
	  sort(med_card_vector.begin(), med_card_vector.end());
	  med_card[i] = med_card_vector[med_card_vector.size() / 2];
	}
    }
}

double compute_conductance(vector<v_size_t> & list_vertices, Graph & G)
{
  set<Vert> verts;
  for (int i = 0; i < list_vertices.size(); ++i)
    verts.insert(vertex(list_vertices[i], G));
  
  return compute_conductance(verts, G);
}

// This function computes the conductance of the set of vertices verts
// in the graph G.
double compute_conductance(set<Vert> & verts, Graph & G)
{
  set<Vert>::iterator sit;
  graph_traits<Graph>::adjacency_iterator ait, aite;
  
  size_t set_volume = 0;
  size_t cut_volume = 0;
  for (sit = verts.begin(); sit != verts.end(); ++sit)
    {
      Vert v = *sit;
      
      graph_traits<Graph>::degree_size_type deg = out_degree(v, G);
	  
      set_volume += deg;
      cut_volume += deg;
      for (tie(ait, aite) = adjacent_vertices(v, G); ait != aite; ++ait)
	{
	  Vert u = *ait;
	  
	  if (verts.find(u) != verts.end())
	    cut_volume--;
	}
    }

  if (set_volume > num_edges(G) - set_volume)
    set_volume = (num_edges(G) - set_volume);

  return (double) cut_volume / (double) set_volume;
}
