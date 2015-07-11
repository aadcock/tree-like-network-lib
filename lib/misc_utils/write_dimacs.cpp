/*
This function writes a graph from my library into a dimacs file:
p edge #vertices #edges
e id1 id2
...

It was originally written for use in a script with Blair Sullivan's tree decomposition executable to produce a tree decomposition with coloring based on the k-core decomposition. Her code uses dimacs format rather than edge lists

The main file below is more general.

By Aaron Adcock, Feb. 2012, PhD Candidate, Stanford University
*/

#include "../graph_lib_boost.hpp"

void write_dimacs(Graph& G, string output_file_name)
{
  graph_traits<Graph>::edge_iterator eit, eite;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  const v_size_t num_vert = num_vertices(G);
  const e_size_t num_edge = num_edges(G)/2;

  ofstream output_file;
  output_file.open(output_file_name.c_str());

  output_file<<"p edge "<<num_vert<<" "<<num_edge<<"\n";

  //Loop writes file
  for(tie(eit,eite)=edges(G);eit!=eite;eit++)
    {
      const Edge e = *eit;
      const Vert s = source(e,G);
      const Vert t = target(e,G);

      const v_size_t si = index[s]+1;
      const v_size_t ti = index[t]+1;
      
      if(si < ti)
	output_file<<"e "<<si<<" "<<ti<<"\n";
    }
  
  output_file.close();
}
