/*This file runs through a directory of .dot files produced by my boost
  k_core algorithm and produces the layouts using GraphViz. 

By Aaron Adcock, PhD Candidate, Stanford University
*/

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
#include <string>
#include <fstream>
#include <iostream>;


using namespace boost;
using namespace std;
typedef adjacency_list<vecS, vecS, directedS> Graph;
typedef adjacency_list<vecS, vecS, undirectedS> UGraph;
typedef graph_traits<Graph>::vertex_descriptor Vert;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<UGraph>::edge_descriptor Edge;

int main(int argc, char *argv[])
{
  
  int max_core = 81;
  dynamic_properties dp;
  ifstream graphFile;
  
  for(int i=0;i<=max_core;i++)
    {
      string s = "graphViz";
      stringstream ss;
      Graph G;
      UGraph UG;

      ss<<i;
      s.append(ss.str());
      s.append(".dot");
      graphFile.open(s);
      read_graphviz(graphFile,G,dp,"node_id");
      graphFile.close();

      copy_graph(G,UG);

      ofstream outViz;
      outViz.open(s);
      write_graphviz(outViz,UG);
      outViz.close();
    }

}
