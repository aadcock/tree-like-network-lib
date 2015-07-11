/*
This compiles into an executable that will generate Erdos-Renyi G_np
random graphs.  The Erdos-Renyi generator used is from the Boost Graph
Library.  See the Boost documentation.

By Aaron Adcock, PhD Candidate, Stanford University Feb. 2012
 */

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/graph_traits.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>


using namespace std;
using namespace boost;

typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
typedef boost::sorted_erdos_renyi_iterator<boost::mt19937, Graph> ERGen;


int main(int argc, char *argv[])
{
  string outputFilePrefix = "er";
  int N = -1;
  int n = -1;
  double p = -1;

  cout<<RAND_MAX<<"\n";

  for(int i=1;i<argc;i++)
    {
      stringstream ss;

      if(strcmp(argv[i], "-o")==0)
	outputFilePrefix = argv[i+1];
      
      if(strcmp(argv[i], "-n")==0)
	{
	  ss<<argv[i+1];
	  ss>>n;
	}

      if(strcmp(argv[i], "-N")==0)
	{
	  ss<<argv[i+1];
	  ss>>N;
	}

      if(strcmp(argv[i], "-p")==0)
	{
	  ss<<argv[i+1];
	  ss>>p;
	}
    }

  if(N<=0)
    {
      cout<<"The number of graphs must be positive. Default value of 5 being used\n";
      N = 5;
    }

  if(n<=0)
    {
      cout<<"The number of vertices must be positive. Default value of 10 being used\n";
      n = 10;
    }

  if(p>1 || p<0)
    {
      cout<<"The edge probability must be between 0 and 1. Default value of .5 being used\n";
      p = .5;
    }

  cout<<"Output File Prefix: "<<outputFilePrefix<<"\n";

  cout<<"# of graphs: "<<N<<"\n# of nodes: "<<n<<"\nEdge probability"<<p<<"\n";


  //The code in the outer loop is based on the erdos_renyi_generator page in the Boost 1.48.0 graph library manual
  boost::mt19937 gen;

  
  for(int i=0;i<N;i++)
    {


      Graph G(ERGen(gen, n, p),ERGen(),n);
      
      graph_traits<Graph>::edge_iterator eit, eite;

      stringstream ss;
      string outputFile = outputFilePrefix;
      
      ss<<"_"<<n<<"_"<<p<<"_"<<i<<".txt";

      outputFile.append(ss.str());

      ofstream output;
      output.open(outputFile.c_str());

      output<<"#This file contains an erdos renyi Gnp graph with n="<<n<<" and p="<<p<<".\n";

      for(tie(eit,eite)=edges(G);eit!=eite;eit++)
	{
	  graph_traits<Graph>::edge_descriptor e = *eit;
	  graph_traits<Graph>::vertex_descriptor s;
	  graph_traits<Graph>::vertex_descriptor t;
	  property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

	  s = source(e,G);
	  t = target(e,G);

	  output<<index[s]<<"\t"<<index[t]<<"\n";
	}
      
      output.close();
    }

  
}
