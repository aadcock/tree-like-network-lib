/*
This file generates N preferential attachment graphs, with n nodes and m new edges added with each node.

gen_pa -o "output file prefix" -N N -n n -m m

By Aaron Adcock, PhD Candidate, Stanford University April. 2012
 */

#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/graph/graph_traits.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>


using namespace std;
using namespace boost;

typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vert;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertices_size_type v_size_t;
//typedef boost::sorted_erdos_renyi_iterator<boost::random::mt19937, Graph> ERGen;

boost::random::mt19937 gen(time(0));
  

v_size_t roll_weighted_die(double weights[],v_size_t n)
{
  boost::random::discrete_distribution<> dist(weights, weights+n);
  return dist(gen);
}



int main(int argc, char *argv[])
{
  string outputFilePrefix = "pa";
  int N = -1;
  int n = -1;
  int m = 2;
  //  double p = 0;
  

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

      if(strcmp(argv[i], "-m")==0)
	{
	  ss<<argv[i+1];
	  ss>>m;
	}

      // if(strcmp(argv[i], "-p")==0)
      // 	{
      // 	  ss<<argv[i+1];
      // 	  ss>>p;
      // 	}
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

  if(m<=0)
    {
      cout<<"The number of edges to add each iteration must be positive. Default value of 2 being used\n";
      m=2;
    }

  // if(p>1 || p<0)
  //   {
  //     cout<<"The edge probability must be between 0 and 1. Default value of .5 being used\n";
  //     p = .5;
  //   }

  cout<<"Output File Prefix: "<<outputFilePrefix<<"\n";

  cout<<"# of graphs: "<<N<<"\n# of nodes: "<<n<<"\n";


  //The code in the outer loop is based on the preferential attachment algorithm described by Barabasi & Albert
 

  
  for(int i=0;i<N;i++)
    {

      //Set up data structures

      Graph G;
      //vector<int> degrees;
      //property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);
     
      //Add start vertices, begin with the vertices connected to each other.
      Vert v1 = add_vertex(G);
      Vert v2 = add_vertex(G);
      add_edge(v1,v2,G);
      add_edge(v2,v1,G);

      //Degree list contains degree of each node, this is fed to roll_weighted_die as the die weights
      double * degree_list;
      degree_list = new double [n];

      for(int j=0;j<n;++j)
	degree_list[j] = 0;

      degree_list[0] = 1;
      degree_list[1] = 1;

      for(int j=2;j<n;++j)
	{
	  Vert v = add_vertex(G);

	  //vector<v_size_t> new_edges;
	  vector<bool> edge_exist(n,false);

	  for(int k=0;k<m;++k)
	    {
	      //Use weighted die to get new edge, check to see if edge already selected (no double edges)
	      v_size_t n_e = roll_weighted_die(degree_list,n);
	      while(edge_exist[n_e])
		n_e = roll_weighted_die(degree_list,n);
	      
	      //add new vertex/edge, update degree_list
	      Vert n_v = vertex(n_e,G);
	      add_edge(v,n_v,G);
	      add_edge(n_v,v,G);
	      ++degree_list[n_e];
	      edge_exist[n_e] = true;
	    }
	  
	  degree_list[j] = m;	  
	}
      
      delete[] degree_list;
      boost::graph_traits<Graph>::edge_iterator eit, eite;

      //Write graph
      stringstream ss;
      string outputFile = outputFilePrefix;
      
      ss<<"_"<<n<<"_"<<m<<"_"<<i<<".txt";

      outputFile.append(ss.str());

      ofstream output;
      output.open(outputFile.c_str());

      output<<"#This file contains an preferential attachment graph with n="<<n<<" and m="<<m<<".\n";

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
