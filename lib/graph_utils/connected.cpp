
/*
  Originally used Boost graph library to find connected components,
  does not any longer.  Uses the BFS.cpp function BFS_vertices_found,
  a breadth-first-search, to find the vertices in connected
  components.  This modification was made in ~Feb. 2012

  By Aaron Adcock, PhD Candidate at Stanford University
  July 2011
*/

#include "../graph_lib_boost.hpp"
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <iostream>
//#include <boost/graph/strong_components.hpp>

using namespace boost;

Graph connected(Graph& G, int s)
{
  return(connected(G));
}

Graph connected(Graph& G)
{
  // Component will be used to store which connected component each
  // vertex belongs to
  v_size_t n = num_vertices(G);
  vector<v_size_t> component(n);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  // Calculates connected components, stores in component
  // v_size_t num = strong_components(G,&component[0]);
  // vector<v_size_t> comp_size(num);
  // Calculates connected components without using Boost
  
  graph_traits<Graph>::vertex_iterator vi, vie;
  vector<bool> found(n,false);
  v_size_t current_comp = 0;
  vector<v_size_t> comp_size;
  for(tie(vi,vie) = vertices(G); vi != vie; ++vi)
    {
      Vert v = *vi;
      v_size_t vind = index[v];
      v_size_t current_comp_size = 0;

      if(!found[vind])
	{
	  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);
	  comp_size.push_back(current_comp_size);
	  ++current_comp;
	}
    }

  //Find largest component
  v_size_t maxComp = 0;
  v_size_t maxInd = 0;
  for(unsigned int i = 0; i < component.size(); ++i)
    {
      int ind = component[i];
      //comp_size[ind]++;
      if(comp_size[ind]>maxComp)
	{
	  maxComp = comp_size[ind];
	  maxInd  = ind;
	}
    }

  //Use subset function to get subgraph of G which consists of largest connected component in G
  vector<Vert> S;
  for(unsigned int i = 0; i<component.size(); i++)
    {
      if(component[i] == maxInd)
	S.push_back(vertex(i,G));
    }  
  
  
  //subset function, S is the subset of vertex indices in subgraph
  return subset(G, S);
}

vector<Vert> list_giant_component(Graph& G)
{
  // Component will be used to store which connected component each
  // vertex belongs to
  v_size_t n = num_vertices(G);
  vector<v_size_t> component(n);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  // Calculates connected components, stores in component
  // v_size_t num = strong_components(G,&component[0]);
  // vector<v_size_t> comp_size(num);
  // Calculates connected components without using Boost
  
  graph_traits<Graph>::vertex_iterator vi, vie;
  vector<bool> found(n,false);
  v_size_t current_comp = 0;
  vector<v_size_t> comp_size;
  for(tie(vi,vie) = vertices(G); vi != vie; ++vi)
    {
      Vert v = *vi;
      v_size_t vind = index[v];
      v_size_t current_comp_size = 0;

      if(!found[vind])
	{
	  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);
	  comp_size.push_back(current_comp_size);
	  ++current_comp;
	}
    }


  //Find largest component
  v_size_t maxComp = 0;
  v_size_t maxInd = 0;
  for(unsigned int i = 0; i < component.size(); ++i)
    {
      int ind = component[i];
      //comp_size[ind]++;
      if(comp_size[ind]>maxComp)
	{
	  maxComp = comp_size[ind];
	  maxInd  = ind;
	}
    }

  //Use subset function to get subgraph of G which consists of largest connected component in G
  vector<Vert> S;
  for(unsigned int i = 0; i<component.size(); i++)
    {
      if(component[i] == maxInd)
	S.push_back(vertex(i,G));
    }  
  
  return S;
}

vector<v_size_t> components(Graph& G)
{
  // Component will be used to store which connected component each
  // vertex belongs to
  v_size_t n = num_vertices(G);
  vector<v_size_t> component(n);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  // Calculates connected components, stores in component
  // v_size_t num = strong_components(G,&component[0]);
  // vector<v_size_t> comp_size(num);
  // Calculates connected components without using Boost
  
  graph_traits<Graph>::vertex_iterator vi, vie;
  vector<bool> found(n,false);
  v_size_t current_comp = 0;
  vector<v_size_t> comp_size;
  for(tie(vi,vie) = vertices(G); vi != vie; ++vi)
    {
      Vert v = *vi;
      v_size_t vind = index[v];
      v_size_t current_comp_size = 0;

      if(!found[vind])
	{
	  BFS_vertices_found(v, G, found, current_comp, component, current_comp_size);
	  comp_size.push_back(current_comp_size);
	  ++current_comp;
	}
    }

  return component;
}
