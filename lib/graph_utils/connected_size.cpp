/*
A short utility for getting basic graph stats quickly
Main file for calculating connected component size, edges, etc.  Output options need updated.

-i input file...REQUIRED
-o output directory...default = "../Graphs/Results/"
-n output files prefix...default = "graph"
-d if y, indicates a directed graph...default = undirected graph
-g if f, indicates only gromov four point, if s indicates only slim delta, if blank both will be calculated
-s if y, suppresses terminal output...default = false
*/

#include "../graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cmath>

int main(int argc, char *argv[])
{
  string inputFile = "../Graphs/test_graph.txt";  

  
  //User input

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	inputFile = argv[i+1];
    }

  ////Default values

  if(inputFile.length()==0)
    {
      cerr<<"You must provide an input file\n";
      exit(EXIT_FAILURE);
    }


  // if(outputFilePrefix.length()==0)
  //   {
  //     cerr<<"Invalid output prefix.  Default output 'graph' will be used"<<"\n";
  //     outputFilePrefix = "graph";
  //   }
  // //Make directories

  // int l = outputDirectory.length();
  // int c = outputDirectory.compare(l-1,1,"/");
  // if(c != 0)
  //   outputDirectory.append("/");


  // int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  // if(direct_made!=0 && errno!=EEXIST)
  //   {
  //     cerr<<"Failed to make output directory\n";
  //     exit(EXIT_FAILURE);
  //   }
  // outputDirectory.append(outputFilePrefix);
  // outputDirectory.append("/");
  // direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);  ////
  
  // if(direct_made!=0 && errno!=EEXIST)
  //   {
  //     cerr<<"Failed to make output directory folder\n";
  //     exit(EXIT_FAILURE);
  //   }

  ////////////////////////////////////
  //----------------------------------
  ////////////////////////////////////

  string delimiter = "\t";
  
  int N, E, N_c, E_c;
  double d, d_t, C, D, D_b, D_p;
  

  Graph G = loadGraph(inputFile,delimiter,false);  
  



  N = num_vertices(G);
  E = num_edges(G)/2;

  G = connected(G,-1);

  N_c = num_vertices(G);
  E_c = num_edges(G)/2;

  d = (1.0*num_edges(G))/(1.0*num_vertices(G));
  
  D   = 0;
  D_b = 0;
  D_p = 0;
  C   = 0;
  d_t = 0;

  // BFS_distance(G);  Too much memory to store whole d-matrix on most machines
  //   cout<<"made it here!\n";
  graph_traits<Graph>::vertex_iterator vi,vie;
  graph_traits<Graph>::adjacency_iterator ai, aie;
  graph_traits<Graph>::adjacency_iterator ai2, aie2;
  vector<Vert> visited;
  property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

  std::pair<Edge,bool> check;
  int count_low_deg = 0;

  vector<v_size_t> distances;
  vector<v_size_t> num_paths;
  double N_paths = 0;


  for(tie(vi,vie) = vertices(G); vi != vie; ++vi)
    {
      Vert v = *vi;

      BFS_source_all(v,G,distances,num_paths);

      size_t size = distances.size();
      double loc_clust = 0;
      double ord_two_deg = 0;
      visited.clear();
      

      for(size_t i=0;i<size;++i)
	{

	  if(distances[i] > D)
	    D = distances[i];

	  if(distances[i] > 0 && distances[i] <= 2)
	    ord_two_deg++;

	  D_p += ((double) num_paths[i])*distances[i];
	  D_b += distances[i];
	  N_paths += num_paths[i];
	}

      d_t += pow(out_degree(v,G),2);
      distances.clear();
      //clustering, NOTE, assumes that this is an undirected graph

      for(tie(ai,aie) = adjacent_vertices(v,G); ai != aie; ++ai)
      	{
      	  Vert a  = *ai;

      	  for(int i=0;i<visited.size();++i)
      	    {
	      check = edge(visited[i],a,G);
	      
      	      if(check.second)
      		loc_clust++;
      	    }

	  visited.push_back(a);
      	}

      int k_i = out_degree(v,G);

      if(k_i>1)
	C += 2.0*loc_clust/(1.0*k_i*(k_i-1.0));
      else
	count_low_deg++;
      
}
  
  D_p /= N_paths;
  D_b /= (pow(N_c,2)-N_c);
  d_t /= (2.0*E_c);
  C   /= (N_c-count_low_deg);

  int pos = inputFile.find(".txt");

  inputFile = inputFile.substr(0,pos);

  cout<<"Input file & N & E & N_c & E_c & deg & d_t & C & Diam & D_b & D_p & # paths\n";
  cout<<inputFile<<" & "<<N<<" & "<<E<<" & "<<N_c<<" & "<<E_c<<" & "<<d<<" & "<<d_t<<" & "<<C<<" & "<<D<<" & "<<D_b<<" & "<<D_p<<" & "<<N_paths<<" & "<<"Description"<<" \\\\\n";
  // for(tie(vi,vie)=vertices(G);vi!=vie;++vi)
  //   {
  //     Vert v = *vi;
  //     v_size_t d = out_degree(v,G);
  //     cout<<v<<"\t"<<d<<"\n";
  //   }

  return(0);
}
