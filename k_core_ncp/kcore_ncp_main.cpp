
/*This file generates an executable that calculates the Network Community
  Profile plot of the input graph.  The intended use is for shell scripts
  in calculating various statistics of graphs.
  
  By Aaron Adcock, PhD candidate at Stanford University
  Oct. 2011 */

#include "graph_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <string.h>

int main(int argc, char *argv[])
{
  string inputFile;
  string outputDirectory  = "./Graphs/Results/";
  string outputFilePrefix = "graph";
  bool directed = true;
  bool suppressOutput = false;
  v_size_t maxC = 0;
  int step_size = 10;
  double alpha  = .001;
  

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
      else if(strcmp(argv[i],"-mc")==0)
	maxC = atoi(argv[i+1]);
      else if(strcmp(argv[i],"-ss")==0)
	step_size = atoi(argv[i+1]);
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

  ////

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

  if(step_size <=0)
    {
      cout<<"Step size must be greater than 0, using default 10\n";
      step_size = 10;
    }

  int l = outputDirectory.length();
  int c = outputDirectory.compare(l-1,1,"/");
  if(c != 0)
    outputDirectory.append("/");

  int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  ////

  Graph G = loadGraph(inputFile,"\t",directed);
  G       = connected(G,-1);
  v_size_t size = num_vertices(G);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

    if(maxC <=0)
    {
      cout<<"Maximum Community Size must be greater than 0, using default graph size/2\n";
      maxC = size/2;
    }

  if(!suppressOutput)
    cout<<"Size of connected component of G "<<size<<"\n";

  vector<int> core = k_core(G);

 
    
  int maxCore = 0;  
  for(int i=0;i<core.size();i++)
    {
      if(maxCore<core[i])
	maxCore = core[i];

      // if(!suppressOutput)
      // cout<<core[i]<<"\t";
    }

  cout<<"\nMax Core is: "<<maxCore<<"\n";
  vector<string> label;
  vector<int> ncpColumn;
  string ncpOutputFile = outputDirectory;
  ncpOutputFile.append(outputFilePrefix);
  ncpOutputFile.append("/");

  string outputPlotPrefix = ncpOutputFile;

  direct_made = mkdir(ncpOutputFile.c_str(),S_IRWXU);

  ncpOutputFile.append(outputFilePrefix);
  ncpOutputFile.append("_kcore_ncp.txt");
  
  k_core_ncp(G,ncpOutputFile,core,alpha,step_size);

  
  outputPlotPrefix.append("kcore_ncp_plots/");
  direct_made = mkdir(outputPlotPrefix.c_str(),S_IRWXU);
  outputPlotPrefix = "./";
  outputPlotPrefix.append(outputFilePrefix);
  outputPlotPrefix.append("_kcore_ncp_maxcore");

  string title = outputFilePrefix;
  title.append(" NCP Plot of k-Cores");

  if(!suppressOutput)
    cout<<"Creating Plot Scripts\n";

  for(int i=0; i<maxCore;i++)
    {
      string tempLabel = outputFilePrefix;
      stringstream ns;
      
      tempLabel.append(" core: ");
      ns<<i;
      string core_num = ns.str();
      tempLabel.append(core_num);
      
      ncpColumn.push_back(i+2);
      label.push_back(tempLabel);

      
      if((i+1)%7 == 0 || (i+1) == maxCore)
	{
	  string outputPlotFile = outputPlotPrefix;
	  outputPlotFile.append(core_num);

	  string outputfilename = outputDirectory;
	  outputfilename.append("/");
	  outputfilename.append(outputFilePrefix);
	  outputfilename.append("/kcore_ncp_plots/");
	  outputfilename.append(outputFilePrefix);
	  outputfilename.append("_kcore_ncp_maxcore");
	  outputfilename.append(core_num);

	  ncpOutputFile = "../";
	  ncpOutputFile.append(outputFilePrefix);
	  ncpOutputFile.append("_kcore_ncp.txt");

	  //cout<<ncpOutputFile<<"\n"<<outputfilename<<"\n"<<outputPlotFile<<"\n";

	  produce_loglog_plot(ncpOutputFile,outputfilename,outputPlotFile,ncpColumn,label,title,"Community Size", "Conductance","",suppressOutput);

	  label.clear();
	  ncpColumn.clear();
	}
    }
  
}
