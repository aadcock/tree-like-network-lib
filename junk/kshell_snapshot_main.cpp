/*
This file generates an executable that calculates the Network Community Profile plot of the input graph.  The intended use is for shell scripts in calculating various statistics of graphs.
  
By Aaron Adcock, PhD candidate at Stanford University
Oct. 2011 
*/

#include "graph_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <string.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>

//Simple color assignment function, partitions core numbers naively (size/# of colors per partition)

void color_core(Graph& G)
{
  graph_traits<Graph>::vertex_iterator vit, vitend;
  int max_core = 0;

  for(tie(vit,vitend)=vertices(G);vit!=vitend;vit++)
    {
      Vert v = *vit;

      if(G[v].core_num > max_core)
	max_core = G[v].core_num;
    }
  
  ////Color graph (rdylgn11)
  int num_colors  = 11;
  vector<int> color_range;
  int interval = max_core/num_colors; 

  if(interval == 0)
    interval++;

  for(int i=0;i<num_colors-1;i++)
    {
      color_range.push_back(interval*(i+1));
    }

  for(tie(vit,vitend)=vertices(G);vit!=vitend;vit++)
    {
      Vert v = *vit;

      if(G[v].core_num < color_range[0])
	G[v].color = 11;
      else if(G[v].core_num < color_range[1])
	G[v].color = 10;
      else if(G[v].core_num < color_range[2])
	G[v].color = 9;
      else if(G[v].core_num < color_range[3])
	G[v].color = 8;
      else if(G[v].core_num < color_range[4])
	G[v].color = 7;
      else if(G[v].core_num < color_range[5])
	G[v].color = 6;
      else if(G[v].core_num < color_range[6])
	G[v].color = 5;
      else if(G[v].core_num < color_range[7])
	G[v].color = 4;
      else if(G[v].core_num < color_range[8])
	G[v].color = 3;
      else if(G[v].core_num < color_range[9])
	G[v].color = 2;
      else
	G[v].color = 1;
    }
	  
	  
  ////

}


int main(int argc, char *argv[])
{
  //Variable declarations for user input
  string inputFile;
  string outputDirectory  = "./Graphs/Results/";
  string outputFilePrefix = "graph";
  bool directed = true;
  bool suppressOutput = false;
  double diffusion_size = .001;
  double alpha  = .001;
  
  //User input

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	{
	  ss<<argv[i+1];
	  ss>>inputFile;
	}
      else if(strcmp(argv[i], "-n")==0)
	{
	  ss<<argv[i+1];
	  ss>>outputFilePrefix;
	}	
      else if(strcmp(argv[i],"-o")==0)
	{
	  ss<<argv[i+1];
	  ss>>outputDirectory;
	}
      else if(strcmp(argv[i],"-a")==0)
	alpha = atof(argv[i+1]);
      else if(strcmp(argv[i],"-e")==0)
	diffusion_size = atof(argv[i+1]);
      else if(strcmp(argv[i],"-d")==0)
	{
	  if(strcmp(argv[i+1], "n")==0)
	    directed = false;
	}
      else if(strcmp(argv[i],"-s")==0)
	{
	  if(strcmp(argv[i+1], "y")==0)
	    suppressOutput = true;
	}
    }

  ////Default values

  if(inputFile.length()==0)
    {
      cout<<"You must provide an input file\n";
      return(0);
    }
  else
    cout<<inputFile<<"\n";

  if(alpha <= 0)
    {
      cout<<"Alpha must be greater than 0, using default .001\n";
      alpha = .001;
    }

  if(diffusion_size <=0)
    {
      cout<<"Diffusion size must be greater than 0, using default .005\n";
      diffusion_size = .005;
    }

  //Make directories

  int l = outputDirectory.length();
  int c = outputDirectory.compare(l-1,1,"/");
  if(c != 0)
    outputDirectory.append("/");


  int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  outputDirectory.append(outputFilePrefix);
  outputDirectory.append("/");
  direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);  ////

  //Create Graph
  Graph G       = loadGraph(inputFile,"\t",directed);
  G             = connected(G,-1);
  v_size_t size = num_vertices(G);
  
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);


  if(!suppressOutput)
    cout<<"Size of connected component of G "<<size<<"\n";

  //Calculate k-core, note that the core # is also stored as a property on each node

  vector<int> core = k_core(G);
    
  int maxCore = 0;  
  for(int i=0;i<core.size();i++)
    {
      if(maxCore<core[i])
	maxCore = core[i];
    }

  //Construct a vector that arranges nodes into k-shells (core # groups)

  vector< vector<int> > core_groups(maxCore+1);
  for(int i=0; i<core.size();i++)
    core_groups[core[i]].push_back(i);
  
  srand(time(NULL));
  
  vector< vector<double> > snapshot_dist;
  //currently vector of distributions for each core, could make 'for each snapshot'

  

  for(int i=0;i<=maxCore;i++)
    {
      vector<int> shell = core_groups[i];

      int size_shell = shell.size();
      int num_iter   = size_shell/10;

      if(size_shell > 1 && size_shell <20)
	num_iter = 2;

      if(size_shell == 1)
	num_iter = 1;

      vector<bool> check(size_shell,false);

      stringstream ss;
      ss<<i;

      string snapshotDir = outputDirectory;
      snapshotDir.append("snapshot/");
      direct_made = mkdir(snapshotDir.c_str(),S_IRWXU);
      snapshotDir.append("shell_");
      snapshotDir.append(ss.str());
      snapshotDir.append("/");
      direct_made = mkdir(snapshotDir.c_str(),S_IRWXU);

      vector< double > spectral_core_dist(maxCore+1,0);
      int sum = 0;

      for(int j=0; j<num_iter;j++)
	{
	  vector<Vert> snapshot;
	  int r = rand() % size_shell;
	  stringstream ss2;

	  while(check[r])
	    r = (r+1) % size_shell;

	  check[r] = true;

	  Vert rv = vertex(shell[r],G);
	  string snapshotFile = snapshotDir;

	  ss2<<shell[r];
	  snapshotFile.append("/");
	  string snapshotAppend = outputFilePrefix;
	  snapshotAppend.append("_spectral_snapshot");
	  snapshotAppend.append(ss2.str());
	  snapshotFile.append(snapshotAppend);
	  
	  string snapshotGviz = snapshotDir;
	  direct_made = mkdir(snapshotGviz.c_str(),S_IRWXU);
	  snapshotGviz.append(snapshotAppend);
	  string snapshotViz = snapshotFile;
	  snapshotFile.append(".txt");

	  spectral_snapshot(G,snapshot,rv,alpha,diffusion_size);

	  for(int k=0;k<snapshot.size();k++)
	    {
	      Vert sv = snapshot[k];
	      int ind = index[sv];
	  
	      spectral_core_dist[(core[ind])]++;
	      sum++;
	    }

	  //Color the vertices of G using core #
	  color_core(G);
	 
	  
	  //Undirected graph makes less confusing visualization in GraphViz
	  UGraph UG_s = undirected_subset(G,snapshot);

	  ofstream outViz;

	  property_map<UGraph , vertex_bundle_t >::type vBundle = get(vertex_bundle,UG_s);
	  
	  vertex_writer<property_map< UGraph, vertex_bundle_t >::type> vw(vBundle);

	  default_writer dw;
	  graph_writer gw;

	  outViz.open(snapshotGviz.c_str());
	  write_graphviz(outViz, UG_s, vw, dw,gw);
	  
	  outViz.close();
	  
	  string gVizCommand = "neato -Tps ";
	  gVizCommand.append(snapshotViz);
	  gVizCommand.append(" -o ");
	  gVizCommand.append(snapshotGviz);
	  gVizCommand.append(".ps");

	  FILE *gp = popen(gVizCommand.c_str(),"w");
	  pclose(gp);
	}

      if(sum!=0)
	{
	  for(int j=0;j<spectral_core_dist.size();j++)
	    spectral_core_dist[j] /= sum;
	}
      snapshot_dist.push_back(spectral_core_dist);
    }
  
  string distOut = outputDirectory;
  distOut.append(outputFilePrefix);
  distOut.append("_spectral_core_dist.txt");

  ofstream outDistFile;
  outDistFile.open(distOut.c_str());
  vector<int> columns_to_plot;
  vector<string> label;

  for(int i=0; i<=maxCore;i++)
    {
      outDistFile<<i<<"\t";
      for(int j=0;j<snapshot_dist.size();j++)
	{
	  vector<double> C;
	  C.assign(snapshot_dist[j].begin(),snapshot_dist[j].end());
	  if(C.size()>i)
	    outDistFile<<C[i]<<"\t";
	  else
	    outDistFile<<"\t";
	}
      outDistFile<<"\n";
      if(i%5==0)
	{
	  columns_to_plot.push_back((i+2));
	  stringstream ss;
	  ss<<i;
	  string column_label = "Shell: ";
	  column_label.append(ss.str());
	  label.push_back(column_label);
	}
    }
  outDistFile.close();

  string distOutPlot = outputDirectory;
  distOutPlot.append(outputFilePrefix);
  distOutPlot.append("_spec_core_dist_plot");
  produce_plot(distOut, distOutPlot, columns_to_plot, label,"Estimate of Diffusion Distribution using k-shell decomposition", "in k-shell #","Percent of Nodes","",suppressOutput);
}
