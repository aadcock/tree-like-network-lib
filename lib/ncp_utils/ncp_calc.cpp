/* This function tries to repeat the ncp_profile plots from the
   "Empirical Comparison of Algorithms for Network Community Detection"
   by Leskovec, Lang, Mahoney.  Uses a spectral method based on
   personalized Page Rank vectors and sweep cuts (see
   approx_Page_Rank.cpp for more information.

   By Aaron Adcock, PhD Candidate at Stanford University
   July 2011 
*/


#include "../graph_lib_boost.hpp"
#include <cstdlib>
#include <ctime>
#include <vector>


// Main function for finding minimum conductance communities across
// various size scales.  This function just keeps track of the minimum
// conductance, not the actual communities.  This could be changed
// using the index vectors below.

/*
  INPUT:
  G:                   Graph under investigation
  max_community:       Maximum community size to look at (ie how far to take sweep cut) 
  alpha:               Teleport probability for personalized Page Rank
  step_size:           Determines how much to reduce epsilon size by between runs.  Min epsilon = 1 / (max_community * 40)

  OUTPUT:
  A double vector 
*/

vector<double> ncp_calc(Graph& G, const v_size_t max_community, const int step_size, const double alpha)
{
  return(ncp_calc(G,max_community,step_size,alpha,false));
}


vector<double> ncp_calc(Graph& G, const v_size_t max_community, const int step_size, const double alpha, bool display_output)
{
  double eps = 1;
  // Pick a seed for random number generator (just used standard c
  // generator, could use better random generator)
  srand(time(NULL));
  int size = num_vertices(G);
  long step = step_size;
  vector<double> finalCond(max_community,1.0);
  vector<int> ind(size,0);
  e_size_t edgeCount = num_edges(G);

  // Create index map for Graph
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  // Page Rank iterations, each time a sweep cut is performed to look at different sizes of communities
  for (unsigned long i = 5; i <= 40 * max_community; i += step)
    {
       
      vector<unsigned v_size_t> original_index (size, 0);       // Keep track of nodes after
      vector<bool> sampled(size, false);
      
      eps = 1.0 / double(i);

      if(display_output)
	cout<<"i: "<<i<<" eps: "<<eps<<"\n";
	
      int numIter = 8 * size / (i) + 3;

      // Need step to get larger, or else running time takes FOREVER.
      long tempStep = step;
      step *= 1.25;
      
      if (step == tempStep)
	step++;

      // cout<<i<<"\n";
      for (int j = 0; j < numIter; j++)
	{
	  int r = rand() % size;              // index of seed node
	  vector<double> s (size,0);
      
	  sampled[r] = true;

	  s[r] = 1;

	  vector<double> p (size,0);
	  approxPR(G,s,alpha,eps,p,r);

	  // Normalize PR vector by outdegree of each vertex
	  for (int k = 0; k < size; k++)
	    { 
	      Vert v = vertex(k, G);
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v, G);
	      p[k] = p[k] / double(outDegree);

	      ind[k] = k;
	    }

	  // Sort ind vector based on normalized Page Rank values
	  sort(ind.begin(), ind.end(), compareIndgtl<vector<double>&>(p));

	  // original_index is now indexed by the original vertex index and points to the rank of the vertex
	  for (int k = 0; k < size; k++)
	    original_index[ind[k]] = k;

	  // Calculate conductance of each cut along PR vector
	  long vol = 0;
	  long out = 0;
	  // 	cout<<"Iter: "<<i<<"\n";

	  // Do Sweep cut up to maximum community size, calculate conductance
	  // as we add vertices to cut
	  finalCond[0] = 1;
	  unsigned int kIter = 0;

	  while (p[ind[kIter]] != 0 && kIter < max_community)
	    {

	      unsigned int k = kIter;
	      vector<int>::iterator it;
	      Vert v    = vertex(ind[k],G);

	      // Out edge iterators
	      graph_traits<Graph>::adjacency_iterator out_i, out_e;

	      // Count outgoing edges for v
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);

	      // update edges leaving cut set (out)
	      out += outDegree;
	      // update internal volume of cut set (vol)
	      vol += outDegree;

	      // Iterate over out going edges of v
	      for (tie(out_i, out_e) = adjacent_vertices(v, G); out_i != out_e; ++out_i)
		{
		  Vert vt = *out_i;

		  // Check if edge crosses set boundary, if so, update
		  // out and vol accordingly
		  if (original_index[index[vt]] < k && original_index[index[vt]] >=0)
		    {
		      // note that 2 is due to outgoing/incoming edge
		      // as G is a directed graph representation of an
		      // undirected graph.
		      out -= 2;  
		      // vol--;
		    }
		}

	      
	      // Pick minimum of cut set and complement volume
	      if(vol > (edgeCount-vol))
		vol = edgeCount - vol;

	    
	    
	      // conductance c is now given by: c = out/vol
	      double c = double(out)/vol;

	      // Check if this beats previous minimum
	      if(c < finalCond[k])
		finalCond[k] = c;
	      
	      kIter++;
	      
	    }
	}
      
    }

  return finalCond;
}

/*
  This version keeps the minimum conductance community found at each
  community size in the vector of vectors best_communities.
 */
vector<double> ncp_calc(Graph& G, const v_size_t max_community, const int step_size, const double alpha, vector< vector<Vert> > & best_communities, bool display_output)
{
  double eps = 1;
  // Pick a seed for random number generator (just used standard c
  // generator, could use better random generator)
  srand(time(NULL));
  int size = num_vertices(G);
  long step = step_size;
  vector<double> finalCond(max_community,1.0);
  vector<int> ind(size,0);
  e_size_t edgeCount = num_edges(G);

  // Initialize best_communities vector
  vector<Vert> dummy;
  best_communities.assign(max_community, dummy);

  // Create index map for Graph
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  // Page Rank iterations, each time a sweep cut is performed to look at different sizes of communities
  for (unsigned long i = 5; i <= 40 * max_community; i += step)
    {
       
      vector<unsigned v_size_t> original_index (size, 0);       // Keep track of nodes after
      vector<bool> sampled(size, false);
      
      eps = 1.0 / double(i);

      if(display_output)
	cout<<"i: "<<i<<" eps: "<<eps<<"\n";
	
      int numIter = 8 * size / (i) + 3;

      // Need step to get larger, or else running time takes FOREVER.
      long tempStep = step;
      step *= 1.25;
      
      if (step == tempStep)
	step++;

      // cout<<i<<"\n";
      for (int j = 0; j < numIter; j++)
	{
	  // index of seed node
	  int seed_node = rand() % size;              
	  vector<double> starting_distribution (size,0);
      
	  sampled[seed_node] = true;

	  starting_distribution[seed_node] = 1;

	  vector<double> page_rank (size,0);
	  approxPR(G, starting_distribution, alpha, eps, page_rank, seed_node);

	  // Normalize PR vector by outdegree of each vertex
	  for (int k = 0; k < size; k++)
	    { 
	      Vert v = vertex(k, G);
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v, G);
	      page_rank[k] = page_rank[k] / double(outDegree);

	      ind[k] = k;
	    }

	  // Sort ind vector based on normalized Page Rank values
	  sort(ind.begin(), ind.end(), compareIndgtl<vector<double>&>(page_rank));

	  // original_index is now indexed by the original vertex
	  // index and points to the personalized page ranking of the
	  // vertex
	  for (int k = 0; k < size; k++)
	    original_index[ind[k]] = k;

	  // Calculate conductance of each cut along PR vector
	  long vol = 0;
	  long out = 0;
	  // 	cout<<"Iter: "<<i<<"\n";

	  // Do Sweep cut up to maximum community size, calculate conductance
	  // as we add vertices to cut
	  finalCond[0] = 1;
	  vector<Vert> temp_community;
	  unsigned int kIter = 0;

	  while (page_rank[ind[kIter]] != 0 && kIter < max_community)
	    {
	      unsigned int k = kIter;
	      vector<int>::iterator it;
	      Vert v    = vertex(ind[k], G);
	      temp_community.push_back(v);

	      // Out edge iterators
	      graph_traits<Graph>::adjacency_iterator out_i, out_e;

	      // Count outgoing edges for v
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);

	      // update edges leaving cut set (out)
	      out += outDegree;
	      // update internal volume of cut set (vol)
	      vol += outDegree;

	      // Iterate over out going edges of v
	      for (tie(out_i, out_e) = adjacent_vertices(v, G); out_i != out_e; ++out_i)
		{
		  Vert vt = *out_i;

		  // Check if edge crosses set boundary, if so, update
		  // out and vol accordingly
		  if (original_index[index[vt]] < k && original_index[index[vt]] >=0)
		    {
		      // note that 2 is due to outgoing/incoming edge
		      // as G is a directed graph representation of an
		      // undirected graph.
		      out -= 2;  
		      // vol--;
		    }
		}
	      
	      // Pick minimum of cut set and complement volume
	      if (vol > (edgeCount - vol))
		vol = edgeCount - vol;

	      // conductance c is now given by: c = out/vol
	      double c = double(out)/vol;

	      // Check if this beats previous minimum
	      if (c < finalCond[k])
		{
		  finalCond[k] = c;
		  best_communities[k] = temp_community;
		}
	      
	      kIter++;  
	    }
	}     
    }

  return finalCond;
}
