/*
This function aims to take snapshots of a graph based on the k-core decomposition.  One option is to use a local diffusion.

By Aaron Adcock, PhD Candidate, Stanford University
 */

#include "../graph_lib_boost.hpp"
#include <cstdlib>
#include <ctime>
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <boost/graph/copy.hpp>



void k_shell_snapshot(Graph& G, vector<int>& core, int core_num, string filename)
{
  srand(time(NULL));
  vector<int> keep_core;
  int size = num_vertices(G);
  //  long edgeCount = num_edges(G);
  vector<float> s(size,0);
  vector<float> p;
  vector<int> count(size,0);
  vector<int> indOrig(size,0);
  int ind;
  //float alpha = .2;
  //float eps   = .001;

  for(int i=0;i<core.size();i++)
    {
      if(core[i]==core_num)
	keep_core.push_back(i);
    }

  int r = 0;

  if(keep_core.size()>0)
    r = rand() % keep_core.size();
  else
    return;

  ind = keep_core[r];
  s[ind] = 1;
  /*
  p = approxPR(G,s,alpha,eps);

  for(int j=0; j<size; j++)
    {
      Vert v = vertex(j,G);
      count[j] = j;
      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);
      p[j] = p[j]/outDegree;
    }
  

  sort(count.begin(),count.end(),compareInd<vector<float>&>(p));

  for(int j=0;j<size;j++)
    indOrig[count[j]] = j;

  long vol = 0;
  long out = 0;
  int maxC = 200;


  float c = 1.0;
  int cutInd = 0;
  for(int j=0; j<maxC; j++)
    {
      //   cout<<"Comm: "<<j<<"\n";
      vector<int>::iterator it;
	   
      Vert v    = vertex(count[j],G);
      //  cout<<"index: "<<index<<"\n";

      //Out edge iterators
      graph_traits<Graph>::out_edge_iterator out_i, out_e;

      //Create index map for Graph
      property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

      //Count outgoing edges for v
      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);

      //update edges leaving cut set (out)
      out += outDegree;
      //update internal volume of cut set (vol)
      vol += outDegree;

      //Iterate over out going edges of v
      for(tie(out_i,out_e)=out_edges(v,G);out_i!=out_e;++out_i)
	{
	  Edge e = *out_i;
	  //	cout<<"node: "<<node1<<"\n";

	  Vert vt = target(e,G);

	  //Check if edge crosses set boundary, if so, update out and vol accordingly
	  if(indOrig[index[vt]]<j && indOrig[index[vt]]>=0)
	    {
	      out -= 2;  //note that 2 is due to outgoing/incoming edge as G is a directed graph representation of an undirected graph.
	      vol--;
	    }
	}
      //cout<<edgeCount<<"\n";
	
      //Pick minimum of cut set and complement volume
      if(vol>(edgeCount-vol))
	vol = edgeCount - vol;
      // cout<<"Comm size "<<j<<": "<<ctime2-ctime1<<"\n";
	    
	    
      //conductance c is now given by: c = out/vol
      float tempC = float(out)/vol;
	    
      //Check if this beats previous minimum
      if(tempC<c)
	{
	  c = tempC;
	  cutInd = j;
	}
    }

for(int j=0;j<=cutInd;j++)
    S.push_back(vertex(count[j],G));
  */
  vector<Vert> S;

  for(int j=0;j<=core.size();j++)
    {
      if(core[j]==core_num)
	S.push_back(vertex(j,G));
    }
  Graph localG = subset(G,S);
  //localG = connected(localG,-1);

  adjacency_list<vecS, vecS, undirectedS> localUGraph;

  copy_graph(localG,localUGraph);

  ofstream outViz;
  
  outViz.open(filename.c_str());


  write_graphviz(outViz,localUGraph);
  outViz.close();
}
