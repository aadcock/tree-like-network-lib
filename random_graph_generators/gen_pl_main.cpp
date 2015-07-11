/*
This file generates a power law graph according to a paper by Fan Chung.  It takes in a few different parameters and produces N power law graphs with n vertices, a power law exponent gamma and ~average degree degree.

By Aaron Adcock, PhD Candidate, Stanford University Feb. 2012
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
#include <cmath>

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
  string outputFilePrefix = "er";
  int N = -1;
  int n = -1;

  double gamma = -2.1;
  int mu_ne    = 1;
  int mu_ee    = 1;
  double mu_d  = 2;

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

      if(strcmp(argv[i], "-gamma")==0)
	{
	  ss<<argv[i+1];
	  ss>>gamma;
	}

      if(strcmp(argv[i], "-degree")==0)
	{
	  ss<<argv[i+1];
	  ss>>mu_d;
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

  if(mu_d<=0)
    {
      cout<<"Average degree must be positive. Default value of 2 being used\n";
      mu_d = 2;
    }

  if(gamma<=0)
    {
      cout<<"Gamma must be positive. Default value of 2.5 being used\n";
      gamma = 2.5;
      mu_ne = 1;
      mu_ee = 1;

      //implicitly mu_d is 2 from above.
    }
  else
    {
      cout<<"Input gamma = "<<gamma<<"\n";
      mu_ne = floor(mu_d*(gamma-2)+.5);
      mu_ee = floor(((3-gamma)/(gamma-2))*mu_ne+.5);      

      cout<<"Gamma attempt: "<<(2+double(mu_ne)/(mu_ne+mu_ee))<<"\n"; 
    }

  // if(p>1 || p<0)
  //   {
  //     cout<<"The edge probability must be between 0 and 1. Default value of .5 being used\n";
  //     p = .5;
  //   }

  cout<<"Output File Prefix: "<<outputFilePrefix<<"\n";

  cout<<"# of graphs: "<<N<<"\n# of nodes: "<<n<<"\n";
 

  
  for(int i=0;i<N;i++)
    {


      Graph G;
      //vector<int> degrees;
      //property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);
     

      Vert v1 = add_vertex(G);
      Vert v2 = add_vertex(G);
      add_edge(v1,v2,G);
      add_edge(v2,v1,G);

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

	  for(int k=0;k<mu_ne;++k)
	    {
	      v_size_t n_e = roll_weighted_die(degree_list,n);
	      while(edge_exist[n_e])
		n_e = roll_weighted_die(degree_list,n);
	      
	      Vert n_v = vertex(n_e,G);
	      add_edge(v,n_v,G);
	      add_edge(n_v,v,G);
	      ++degree_list[n_e];
	      edge_exist[n_e] = true;
	    }

	  for(int k=0;k<mu_ee;++k)
	    {
	      v_size_t v_1 = vertex(roll_weighted_die(degree_list,n),G);

	      bool flag = false;
	      int count = 0;

	      while(flag != true && count < 10)
		{
		  v_size_t v_2 = vertex(roll_weighted_die(degree_list,n),G);

		  flag = true;
	      
		  graph_traits<Graph>::adjacency_iterator ait, aite;

		  for(tie(ait,aite)=adjacent_vertices(v_1,G);ait!=aite;++ait)
		    {
		      v_size_t v_a = *ait;
		      
		      if(v_a == v_2)
			flag = false;
		    }
			   
		  if(v_1==v_2)
		    flag = false;
	      
		  if(flag)
		    {
		      add_edge(v_1,v_2,G);
		      add_edge(v_2,v_1,G);
		  
		      ++degree_list[v_1];
		      ++degree_list[v_2];
		    }

		  ++count;
	      
		}
	      
	      degree_list[j] = mu_ne;	  
	    }
	}
      delete[] degree_list;
      boost::graph_traits<Graph>::edge_iterator eit, eite;

      stringstream ss;
      string outputFile = outputFilePrefix;
      
      ss<<"_"<<n<<"_"<<mu_d<<"_"<<2+mu_ne/(mu_ne+mu_ee)<<"_"<<i<<".txt";

      outputFile.append(ss.str());

      ofstream output;
      output.open(outputFile.c_str());

      output<<"#This file contains an power law graph with n="<<n<<" avg degree="<<mu_d<<" and attempted gamma = "<<(2 + double(mu_ne)/(mu_ne+mu_ee))<<".\n";

      property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

      for(tie(eit,eite)=edges(G);eit!=eite;eit++)
	{
	  graph_traits<Graph>::edge_descriptor e = *eit;
	  graph_traits<Graph>::vertex_descriptor s;
	  graph_traits<Graph>::vertex_descriptor t;

	  s = source(e,G);
	  t = target(e,G);

	  output<<index[s]<<"\t"<<index[t]<<"\n";
	}
      
      output.close();
    }

  
}
