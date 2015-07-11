/*This function takes a graph G and a set of vertices S and returns
  the induced subgraph of the nodes whose indices belong to S. Assumes use of Boost Graph Library

  by Aaron Adcock, PhD candidate at Stanford Universisty
  July 2011*/


#include "../graph_lib_boost.hpp"
#include <boost/graph/copy.hpp>
#include <list>

Graph subset(Graph& G, vector<Vert>& S)
{
  //Edge iterator class

  if(S.size() == 0)
    {
      Graph Gnew;
      return(Gnew);
    }

  graph_traits<Graph>::edge_iterator it, itend;
  graph_traits<Graph>::vertex_iterator vit, vitend;
  int size = S.size();
  vector<v_size_t> Sindex;

  //This set of code defines an index set for a graph of type 'Graph',
  //this can then be used to reference different vertices in the graph
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  property_map<Graph, vertex_bundle_t>::type vertProp = get(vertex_bundle,G);

  vector<intPair> edgeList;
  vector<vert_prop> vertPropList;
  vector<bool> inGraph(size,false);
  vector<Vert> Snew;
  Graph Gnew;
  //Vert is a wrapper for an integer type, so it can be sorted, I think

  for(int i = 0; i < size; i++)
    Sindex.push_back(index[(S[i])]);

  // sort(S.begin(),S.end(),compareIndltg<vector<v_size_t>&>(Sindex));
  sort(Sindex.begin(),Sindex.end());
  vector<v_size_t>::iterator vec_it;

  vec_it = unique(Sindex.begin(),Sindex.end());
  Sindex.resize(distance(Sindex.begin(),vec_it));
  size = Sindex.size();

  for(int i = 0; i < size; i++)
    {
      Vert nv = add_vertex(Gnew);
      Vert ov = vertex(Sindex[i],G);
      v_size_t ind = index[ov];

      Snew.push_back(nv);

      Gnew[nv].core_num = G[ov].core_num;
      Gnew[nv].color = G[ov].color;
      Gnew[nv].distribution = G[ov].distribution;
      Gnew[nv].prevIndices = G[ov].prevIndices;
      Gnew[nv].prevIndices.push_back(ind);
    }

  // Iterate over the edges of G

  for(int i = 0; i < size; i++)
    {
      Vert v = vertex(Sindex[i], G);
      graph_traits<Graph>::adjacency_iterator ait, aitend;

      for(tie(ait, aitend) = adjacent_vertices(v, G); ait != aitend; ++ait)
	{
	  Vert u      = *ait;	  
	  v_size_t ui = index[u];
	  int ind     = binarySearch(Sindex,ui);
	  
	  if(ind>=0)
	    {
	      bool edgeExist;
	      Edge e;
	      
	      tie(e,edgeExist) = edge(Snew[i],Snew[ind],Gnew);
	      
	      if(!edgeExist)
		add_edge(Snew[i],Snew[ind],Gnew);

	      tie(e,edgeExist) = edge(Snew[ind],Snew[i],Gnew);
	      
	      if(!edgeExist)
		add_edge(Snew[ind],Snew[i],Gnew);
	    }
	}
      
    }

  // for(tie(it,itend)=edges(G);it!=itend;++it)
  //   {
  //     Edge currEdge = *it;
  //     Vert v1       = source(currEdge,G);
  //     v_size_t v1i  = index[v1];
  //     Vert v2       = target(currEdge,G);
  //     v_size_t v2i  = index[v2];
  //     int ind1; 
  //     int ind2;


  //     //Are both edges in the new graph?
  //     ind1 = binarySearch(Sindex,v1i);
  //     ind2 = binarySearch(Sindex,v2i);

  //     if(ind1>=0 && ind2>=0)
  // 	{
  // 	  // intPair newEdge;
  // 	  // newEdge.first  = ind1;
  // 	  // newEdge.second = ind2;

  // 	  // edgeList.push_back(newEdge);

  // 	  inGraph[ind1] = true;
  // 	  inGraph[ind2] = true;

  // 	}
  //   } 
 
  return Gnew;
}

// Index of nodes is provided directly to function
// Graph subset(Graph& G, vector<v_size_t>& Sindex)
// {
//   //Edge iterator class

//   if(Sindex.size() == 0)
//     {
//       Graph Gnew;
//       return(Gnew);
//     }

//   graph_traits<Graph>::edge_iterator it, itend;
//   graph_traits<Graph>::vertex_iterator vit, vitend;
//   int size = Sindex.size();

//   property_map<Graph, vertex_bundle_t>::type vertProp = get(vertex_bundle,G);

//   vector<intPair> edgeList;
//   vector<vert_prop> vertPropList;
//   vector<bool> inGraph(size,false);
//   vector<Vert> Snew;
//   Graph Gnew;

//   // sort(S.begin(),S.end(),compareIndltg<vector<v_size_t>&>(Sindex));
//   sort(Sindex.begin(),Sindex.end());
//   vector<v_size_t>::iterator vec_it;

//   vec_it = unique(Sindex.begin(),Sindex.end());
//   Sindex.resize(distance(Sindex.begin(),vec_it));
//   size = Sindex.size();

//   for(int i = 0; i < size; i++)
//     {
//       Vert nv = add_vertex(Gnew);
//       Vert ov = vertex(Sindex[i], G);
//       v_size_t ind = index[ov];

//       Snew.push_back(nv);

//       Gnew[nv].core_num = G[ov].core_num;
//       Gnew[nv].color = G[ov].color;
//       Gnew[nv].distribution = G[ov].distribution;
//       Gnew[nv].prevIndices = G[ov].prevIndices;
//       Gnew[nv].prevIndices.push_back(ind);
//     }

//   // Iterate over the edges of G, keep those in subgraph
//   for(int i = 0; i < size; i++)
//     {
//       Vert v = vertex(Sindex[i], G);
//       graph_traits<Graph>::adjacency_iterator ait, aitend;

//       for(tie(ait, aitend) = adjacent_vertices(v, G); ait != aitend; ++ait)
// 	{
// 	  Vert u      = *ait;	  
// 	  v_size_t ui = index[u];
// 	  int ind     = binarySearch(Sindex,ui);
	  
// 	  // If target node is in subset, add edge if not already
// 	  // added (note this currently adds edges in both directions
// 	  if(ind >= 0)
// 	    {
// 	      bool edgeExist;
// 	      Edge e;
	      
// 	      tie(e,edgeExist) = edge(Snew[i],Snew[ind],Gnew);
	      
// 	      if(!edgeExist)
// 		add_edge(Snew[i],Snew[ind],Gnew);

// 	      tie(e,edgeExist) = edge(Snew[ind],Snew[i],Gnew);
	      
// 	      if(!edgeExist)
// 		add_edge(Snew[ind],Snew[i],Gnew);
// 	    }
// 	}
      
//     } 
//   return Gnew;
// }


//This function is used to create an undirected, unweighted graph from a directed graph.  Each directed edge between two nodes is converted into a single undirected edge
UGraph undirected_subset(Graph& G, vector<Vert>& S)
{
  //Edge iterator class
  graph_traits<Graph>::edge_iterator it, itend;
  graph_traits<Graph>::vertex_iterator vit, vitend;
  int size = S.size();
  vector<v_size_t> Sindex;

  //This set of code defines an index set for a graph of type 'Graph',
  //this can then be used to reference different vertices in the graph
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  property_map<Graph, vertex_bundle_t>::type vertProp = get(vertex_bundle,G);

  vector<intPair> edgeList;
  vector<vert_prop> vertPropList;
  vector<bool> inGraph(size,false);
  vector<Vert> Snew;
  Graph Gnew;
  //Vert is a wrapper for an integer type, so it can be sorted, I think

  for(int i=0;i<size;i++)
    Sindex.push_back(index[(S[i])]);

  // sort(S.begin(),S.end(),compareIndltg<vector<v_size_t>&>(Sindex));
  sort(Sindex.begin(),Sindex.end());
  
  //Iterate over the edges of G

  for(int i=0;i<size;i++)
    {
      Vert nv = add_vertex(Gnew);
      Vert ov = vertex(Sindex[i],G);
      v_Usize_t ind = index[ov];

      Snew.push_back(nv);
       
      Gnew[nv].core_num = G[ov].core_num;
      Gnew[nv].color = G[ov].color;
      Gnew[nv].distribution = G[ov].distribution;
      Gnew[nv].prevIndices.push_back(ind);
    }

  for(int i=0; i<size; i++)
    {
      Vert v = vertex(Sindex[i],G);
      graph_traits<Graph>::adjacency_iterator ait, aitend;

      for(tie(ait,aitend)=adjacent_vertices(v,G);ait!=aitend;ait++)
	{
	  Vert u      = *ait;	  
	  v_size_t ui = index[u];
	  int ind     = binarySearch(Sindex,ui);
	  
	  if(ind>=0)
	    {
	      bool edgeExist;
	      Edge e;
	      
	      tie(e,edgeExist) = edge(Snew[i],Snew[ind],Gnew);
	      
	      if(!edgeExist)
		{
		  tie(e,edgeExist) = edge(Snew[ind],Snew[i],Gnew);
		  if(!edgeExist)
		    add_edge(Snew[i],Snew[ind],Gnew);
		}
	      
	    }
	}
      
    }

  UGraph UG;

  copy_graph(Gnew, UG);

  return UG;
}

vector<Vert> convert_indices_to_vertices(Graph & G, vector<v_size_t> & index)
{
  vector <Vert> vertex_list;
  for (size_t i = 0; i < index.size(); ++i)
    vertex_list.push_back( vertex(index[i], G) );

  return vertex_list;
}
