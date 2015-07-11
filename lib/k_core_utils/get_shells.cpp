/*
This function takes a list of k-core numbers (shells) and returns the subgraph induced by keeping all nodes in those k-shells.

By Aaron Adcock, PhD Candidate, Stanford University
August 2012
 */

#include "../graph_lib_boost.hpp"


Graph get_k_shells(vector<v_size_t> & shells, Graph & G)
{
  v_size_t n = num_vertices(G);
  vector<int> core;
  vector<Vert> vert_to_keep;

  graph_traits<Graph>::vertex_iterator vit, vite;

  //Sort shell vector
  sort(shells.begin(),shells.end());

  //Check if k_core exists for all vertices
  core = k_core(G);


  for(tie(vit,vite)=vertices(G);vit!=vite;++vit)
    {
      Vert v      = *vit;
      v_size_t shell_num = G[v].core_num;

      int found = binarySearch(shells,shell_num);

      if(found>=0)
	vert_to_keep.push_back(v);
    } 

  return(subset(G,vert_to_keep));
}

Graph get_k_shells(vector<v_size_t> & shells, vector<v_size_t> & core, Graph & G)
{

  v_size_t n = num_vertices(G);
  vector<Vert> vert_to_keep;

  graph_traits<Graph>::vertex_iterator vit, vite;
  property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

  for(tie(vit,vite)=vertices(G);vit!=vite;++vit)
    {
      Vert v      = *vit;
      v_size_t vi = index[v];
      v_size_t shell_num = core[vi];

      int found = binarySearch(shells,shell_num);

      if(found>=0)
	vert_to_keep.push_back(v);
    } 

  return(subset(G,vert_to_keep));

}
