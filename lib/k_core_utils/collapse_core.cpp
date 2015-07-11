/*
This function takes a k-core of a node and collapses it into a single node.  All edges connected to the k-core will now be connected to this node.

By Aaron Adcock, PhD Candidate, Stanford University
October 2012
 */

#include "../graph_lib_boost.hpp"

using namespace boost;


Graph collapse_core(v_size_t core_to_collapse, Graph & G)
{
  Graph G_collapse(G);
  
  v_size_t n = num_vertices(G_collapse);
  vector<int> core;
  vector<Vert> vert_to_keep;

  graph_traits<Graph>::vertex_iterator vit, vite;

  //Begin with G_collapse = G
  

  //Check if k_core exists for all vertices, need to improve this by making it a graph property
  core = k_core(G_collapse);
	  
  int max_core = 0;
  for(int i=0;i<core.size();++i)
    {
      if(core[i] > max_core)
	max_core = core[i];
    }
  
  //collapse k-core into core_node
  Vert core_node = add_vertex(G_collapse);

  graph_traits<Graph>::vertex_iterator vi, vie;

  //Iterate over nodes, to look for nodes in k-core
  for(tie(vi,vie)=vertices(G_collapse); vi!=vie; ++vi)
    {
      Vert v = *vi;

      v_size_t v_core = G_collapse[v].core_num;
      
      //if node is in k-core, add edges to core_node, and then remove edges
      //Note!  This assumes that I am using a directed graph with only reciprocal edges to 
      //represent an undirected graph.  This is currently what I do...but I think I should change it
      if(v_core >= core_to_collapse && v!=core_node)
	{

	  graph_traits<Graph>::out_edge_iterator oei,oeie;
	  for(tie(oei,oeie)=out_edges(v,G_collapse);oei!=oeie;++oei)
	    {
	      Edge e = *oei;

	      Vert s = source(e,G_collapse);
	      Vert t = target(e,G_collapse);

	      // cout<<"Vertex v: "<<v<<"\n";
	      // cout<<"Vertex s: "<<s<<"\n";
	      // cout<<"Vertex t: "<<t<<"\n";

	      v_size_t t_core = G_collapse[t].core_num;

	      if(t_core < core_to_collapse)
		{

		  pair<Edge,bool> check;

		  check = edge(core_node,t,G_collapse);
		  if(!check.second)
		    add_edge(core_node,t,G_collapse);

		  check = edge(t,core_node,G_collapse);
		  if(!check.second)
		    add_edge(t,core_node,G_collapse);
		}
	    }
     
	  
	  
	  clear_vertex(v,G_collapse);
	}
    }

  //Get a connected component of G_collapse, removes all cleared vertices
  G_collapse = connected(G_collapse);



  return(G_collapse);
}

