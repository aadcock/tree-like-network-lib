/*
This file generates a hierarchical graph according to the paper:
Hierarchical Organization in Complex Networks,
E. Ravasz and A.L. Barabasi,
Physical Review E 67

Depending on the values given, the output can be the 'stochastic' or deterministic version
of the graph described in the paper.

gen_hier -o output -n n -m m -p p

By Aaron Adcock, PhD Candidate, Stanford University June 2012
 */

#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>


#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <vector>


using namespace std;
using namespace boost;

typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
typedef boost::adjacency_list<vecS,vecS,undirectedS> UGraph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vert;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertices_size_type v_size_t;
//typedef boost::sorted_erdos_renyi_iterator<boost::random::mt19937, Graph> ERGen;

boost::random::mt19937 gen(time(0));
  

v_size_t roll_weighted_die(double weights[],v_size_t n)
{
  //Rolls a weighted die of size/probabilities given by weights
  boost::random::discrete_distribution<> dist(weights, weights+n);
  return dist(gen);
}



int main(int argc, char *argv[])
{
  string outputFilePrefix = "hier";
  int n = -1;
  int m = 3;
  double p = 1;
  
  bool stoch_flag = false;

  //Parsing the inputs
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

      if(strcmp(argv[i], "-m")==0)
	{
	  ss<<argv[i+1];
	  ss>>m;
	}

      if(strcmp(argv[i], "-p")==0)
	{
	  ss<<argv[i+1];
	  ss>>p;

	  stoch_flag = true;
	  cout<<p<<"\n";
	}

    }

  if(n<=0)
    {
      cout<<"The number of vertices must be positive. Default value of 10 being used\n";
      n = 10;
    }

  if(m<=2)
    {
      cout<<"The hierarchy base must be greater than 2, Default value of 3 being used\n ";
      m=3;
    }

  if(stoch_flag && (p>1 || p<0))
    {
      cout<<"The node fraction must be between 0 and 1. Default of .5 being used\n";
      p = .5;
    }


  cout<<"Output File Prefix: "<<outputFilePrefix<<"\n";

  cout<<"# of nodes: "<<n<<"\n";
  cout<<"Base clique: "<<m<<"\n";
  
  if(stoch_flag)
    cout<<"Probability of adding edge: "<<p<<"\n";

  //Finding the largest, without exceeding n, deterministic hierarchical graph
  double log_level;
  int num_levels;
  log_level = log(n)/log(m);
  num_levels = floor(log_level);

  n = pow(m,num_levels);

  //The code in the outer loop is based on the preferential attachment algorithm described by Barabasi & Albert
 
  Graph G(n);
  graph_traits<Graph>::vertex_iterator vi, vie;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  
  tie(vi,vie) = vertices(G);

  //Determines whether a node is a periphery node or not, ie will we connect it deeper in
  vector<bool> periphery(n,true);

  if(!stoch_flag)
    {
      //Deterministic version of graph

      //Make Cliques out of every m nodes
      for(int i=0;i<n;i+=m)
	{
	  for(int j=0;j<m;++j)
	    {
	      Vert s = *(vi+i+j);
	      if(j==0)
		periphery[i+j] = false;

	      for(int k=j+1;k<m;++k)
		{		    
		  Vert t = *(vi+i+k);
		  add_edge(s,t,G);
		}
	    }
	}
      
      //Make larger 'modules', m^2, m^3, etc.
      for(int i=1;i<num_levels;++i)
	{
	  double step = pow(m,i+1);
	  double prev_step = pow(m,i);

	  for(int j=0;j<n;j+=step)
	    {
	      int center = j;
	      Vert c = *(vi+center);
	      for(int k=0;k<step;k+=prev_step)
		{
		  for(int l=0;l<prev_step;++l)
		    {
		      if(k==0)
			{
			  periphery[j+k+l] = false;
			}
		      else if(periphery[j+k+l])
			{
			  Vert s = *(vi+j+k+l);
			  add_edge(s,c,G);
			}
		    }
		}
	    }
	}
    }
  else
    {
      //The Stochastic version of the graph
      double * degree_list;
      degree_list = new double [n];

      for(int i=0;i<n;++i)
	degree_list[i] = 0;
	        
      
      for(int i=0;i<n;i+=m)
	{
	  for(int j=0;j<m;++j)
	    {
	      Vert s = *(vi+i+j);
	      if(j==0)
		periphery[i+j] = false;

	      for(int k=j+1;k<m;++k)
		{		    
		  Vert t = *(vi+i+k);
		  add_edge(s,t,G);
		  v_size_t s_i = index[s];
		  v_size_t t_i = index[t];

		  ++degree_list[s_i];
		  ++degree_list[t_i];
		}
	    }
	}
      //make larger modules

	  double * p_die;
	  p_die = new double [2];
	  p_die[0] = p;
	  p_die[1] = 1-p;
    
      for(int i=1;i<num_levels;++i)
	{

	  p_die[0] = pow(p,i);
	  p_die[1] = 1-pow(p,i);
	  double step = pow(m,i+1);
	  double prev_step = pow(m,i);

	  int p_step = int(prev_step);
	  //cout<<p_step<<"\n";
	  //Preferential attachment to central nodes
	  double * central_degree_list;
	  central_degree_list = new double [p_step];
	  for(int j=0;j<p_step;++j)
	    {
	      central_degree_list[j] = degree_list[j];
	    }

	  vector<int> repeat_edges(int(step),-1);
	  
	  for(int j=prev_step;j<step;++j)
	    {
	      int new_edge = roll_weighted_die(p_die,2);

	      if(new_edge==0)
		{
		  v_size_t n_e = roll_weighted_die(central_degree_list,p_step);
		  //cout<<"Node: "<<j<<" Rolled: "<<n_e<<"\n";

		  Vert v   = *(vi+j);
		  Vert n_v = vertex(n_e,G);
		  add_edge(v,n_v,G);
		  ++degree_list[j];
		  ++degree_list[n_e];
		  ++central_degree_list[n_e];
		  
		  repeat_edges[j] = n_e;
		}
	    }
	  
	  //Replicating edges in each copy of the central m^i nodes
	  for(int j=step;j<n;j+=step)
	    {
	      for(int k=prev_step;k<step;++k)
		{
		  if(repeat_edges[k]>=0)
		    {
		      v_size_t n_e = repeat_edges[k]+j;

		      Vert v = *(vi+j+k);
		      Vert n_v = vertex(n_e,G);

		      add_edge(v,n_v,G);
		      ++degree_list[j+k];
		      ++degree_list[n_e];
		    }
		}
	    }


	  delete[] central_degree_list;
	}
      
      delete[] p_die;
      delete[] degree_list;
    }  
  

  //Create file
  boost::graph_traits<Graph>::edge_iterator ei, eie;


  stringstream ss;
  string outputFile = outputFilePrefix;
      
  if(stoch_flag)
    ss<<"_"<<n<<"_"<<m<<"_"<<p<<".txt";
  else
    ss<<"_"<<n<<"_"<<m<<"_d.txt";

  outputFile.append(ss.str());

  ofstream output;
  output.open(outputFile.c_str());

  output<<"#This file contains a Barabasi Ravasz hierarchical graph with  n="<<n<<" base m ="<<m<<" and p = "<<p<<".\n";


  for(tie(ei,eie)=edges(G);ei!=eie;ei++)
    {
      Edge e = *ei;
      Vert s;
      Vert t;

      s = source(e,G);
      t = target(e,G);

      output<<index[s]<<"\t"<<index[t]<<"\n";
    }
      
  output.close();

  // UGraph UG;
  // copy_graph(G,UG);

  // ofstream outViz;
  // outViz.open("hier_gviz.dot");

  // write_graphviz(outViz,UG);
  // outViz.close();
}
