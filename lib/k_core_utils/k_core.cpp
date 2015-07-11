/*
This function is used to compute the k-core of a given graph.  The k-core is a maximal subgraph where each vertex has a degree greater than k WITHIN the subgraph (or k-core). Loosely follows "An O(m) Algorithm for Cores Decomposition of Networks", by Bategelj, Zaversnik 2003 (this paper is the source of much of the notation, so I recommend reading this paper)

By Aaron Adcock, PhD Candidate, Stanford University
August 2011
*/

#include "../graph_lib_boost.hpp"

using namespace std;
using namespace boost;

// void k_core(Graph& G)
// {
//   vector<int> dummy_var = k_core(G);
// }

vector<int> k_core(Graph& G)
{
  int num_vert = num_vertices(G);
  int max_degree = 0;
  vector<int> vert(num_vert, 0), pos(num_vert, 0);
  vector<int> deg(num_vert, 0);

  graph_traits<Graph>::adjacency_iterator ait, aitend;
  graph_traits<Graph>::vertex_iterator vit, vitend;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  for (tie(vit, vitend) = vertices(G); vit != vitend; ++vit)
    {
      int vind = index[(*vit)];
      //cout<<vind<<"\n";
      deg[vind] = out_degree(*vit, G);
      if (deg[vind] > max_degree)
	max_degree = deg[vind];
    }

  vector<int> bin(max_degree + 1, 0);
  vector<int>::iterator it;

  // Put number of nodes of degree in each slot of bin
  for (it = deg.begin(); it<deg.end(); ++it)
    bin[(*it)]++;
  
  // cumsum on bin
  int start = 0;
  for (int i = 0; i <= max_degree; i++)
    {
      int num = bin[i];
      bin[i]  = start;
      start  += num;
    }
  
  // Figure out which bin each vertex is in, store in pos
  // Reorder using vert
  for (int i = 0; i < num_vert; i++)
    {
      pos[i]       = bin[(deg[i])];
      vert[(pos[i])] = i;
      bin[(deg[i])]++;
    }
  // for(int i=0;i<num_vert;i++)
  //   cout<<"Node: "<<i+1<<" Deg: "<<deg[i]<<" Pos: "<<pos[i]<<"\n";

  // shift bin up
  for (int i = max_degree; i >= 1; i--)
    bin[i] = bin[i - 1];
  
  bin[0] = 0;

  for (int i = 0; i < num_vert; i++)
    {
      int v = vert[i];
      for (tie(ait, aitend) = adjacent_vertices(vertex(v, G), G); ait != aitend; ++ait)
	{
	  int u = index[(*ait)];
	  if (deg[u] > deg[v])
	    {
	      // u is the neighbor, w is the first node with same
	      // degree as u. Since we are decrementing u, we want to
	      // swap their spots (w will still be in the deg_u group,
	      // but now u is adjacent to the deg_u - 1 group.
	      int deg_u = deg[u];
	      int pos_u = pos[u];
	      int pos_w = bin[deg_u];
	      int  w = vert[pos_w];
	      
	      if (u != w)
		{
		  pos[u]   = pos_w;
		  vert[pos_u] = w;
		  pos[w]   = pos_u;
		  vert[pos_w] = u;
		}
	      bin[deg_u]++;
	      deg[u]--;
	    }
	}
    }
  for (tie(vit, vitend) = vertices(G); vit != vitend; ++vit)
    {
      Vert v  = *vit;
      int ind = index[v];
      
      G[v].core_num  = deg[ind];
      G[graph_bundle].core_calc = true;
    }

  return deg;
}

