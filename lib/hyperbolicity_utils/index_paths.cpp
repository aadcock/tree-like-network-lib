/*
This file contains a function which takes a source node and a target
node and finds indexes ALL shortest paths between them.  This is done
using the breadth first search predecessors (precomputed) and a stack.
When completed, every node on a shortest path will have the index of
each path it is on associated with it.  This function is intended for
use with calc_gromov.cpp and BFS.cpp.  Note that the size of the
graph, the index map, and the vertices on any shortest path from s to
t can be provided to function (rather than just G) if precalculated

By Aaron Adcock, PhD candidate at Stanford University, Jan. 2012
 */

#include "../graph_lib_boost.hpp"


unsigned int index_paths(const Vert& s, const Vert& t, const property_map<Graph,vertex_index_t>::type& index, const v_size_t& n, const vector< vector<Vert> >& p, vector< vector<unsigned int> >& ind)
{
  //property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);
  //v_size_t n = num_vertices(G);

  if(ind.empty())
    {
      const vector<unsigned int> dummy;
      const vector< vector<unsigned int> > ind_temp(n,dummy);
      ind = ind_temp;
    }
  else if(ind.size()<n)
    {
      const vector<unsigned int> dummy;
      for(size_t i=0;i<ind.size();i++)
  	ind[i].clear();

      for(size_t i=0;i<n-ind.size();i++)
  	ind.push_back(dummy);
    }
  else 
    {
      for(v_size_t i=0;i<n;i++)
	ind[i].clear();
    }
  
  // for(size_t i=0;i<n;i++)
  //   {
  //     ind_meta[i][0]=0;
  //   }

  
  unsigned int path = 0;

  //vector hay_stack is being used as a stack, with the exception that I will iterate through the nodes whenever a path is found
  vector<path_node> hay_stack;
  hay_stack.reserve(CHUNKSIZE);

  //  v_size_t si = index[s];
  const v_size_t ti = index[t];

  path_node init;
  path_node top;
  path_node next;

  Vert pv;
  v_size_t pvi;
  v_size_t vi;

  init.v = t;
  init.predecessors = ti;
  init.count = 0;

  hay_stack.push_back(init);

  while(!hay_stack.empty())
    {

      top = hay_stack.back();
      hay_stack.pop_back();

      if(top.count < p[top.predecessors].size())
	{
	  if(top.v == s)
	    {
	      hay_stack.push_back(top);
	      for(size_t i=0;i<hay_stack.size(); i++)
		{
		  pv  = hay_stack[i].v;
		  pvi = index[pv];
	      
		  ind[pvi].push_back(path);

		
		}
	      path++;
	      hay_stack.pop_back();
	      
	      //cout<<"Source pop!\n";
	    }
	  else
	    {
	      next.v = p[top.predecessors][top.count];
	      vi = index[next.v];
	      top.count++;
	      next.predecessors = vi;
	      next.count = 0;
	      hay_stack.push_back(top);
	      hay_stack.push_back(next);
	      //cout<<"Push!\n";
	    }
	}

	
    }
  return(path);
}
