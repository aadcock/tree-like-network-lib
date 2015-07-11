/*This file is used for testing my graph functions as I build them using the boost library, it is also used for calculating ncp_plots using approxPR 

By Aaron Adcock, PhD candidate at Stanford University
July 2011*/

#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "graph_lib_boost.hpp"
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>


int main(int argc, char *argv[])
{

  string inputFile;
  string outputDirectory  = "./Graphs/Results/";
  string outputFilePrefix = "graph";
  bool directed = false;
  bool suppressOutput = false;
  

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

  int l = outputDirectory.length();
  int c = outputDirectory.compare(l-1,1,"/");
  if(c != 0)
    outputDirectory.append("/");

  int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  ////

  Graph  G = loadGraph(inputFile,"\t",directed);

  int size = num_vertices(G);

  cout<<"Size of G "<<size<<"\n";

  //Size of connected component
  
  Graph Gfull = G;
  G = connected(G,-1);
  cout<<"Size of Gconn "<<num_vertices(G)<<"\n";

  size = num_vertices(G);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  vector<int> core = k_core(G);
  
  // cout<<"Core numbers: \n";
  // for(int i=0; i<core.size();i++)
  //   cout<<core[i]<<"\t";

  // cout<<"\n";


  int max_core = 0;
 
  for(int i=0; i<core.size();i++)
    {      
      if(core[i]>max_core)
	max_core = core[i];
    }
  
	////
  
  string outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);

  vector<vector<float> > core_stat;

  k_core_stat(core_stat,G,outputFile);
  
  if(!suppressOutput)
    {
      vector<float> core_dist(max_core+1,0);
  
      for(int i=0;i<core.size();i++)
	core_dist[(core[i])]++;

      
      // for(int iter=0;iter<max_core;iter++)
      //   {
      //     //cout<<iter+1<<" ";
      //     string filename = "./Temp/graphViz";
      //     stringstream ss;
      //     ss<<(iter+1);
      //     filename.append(ss.str());
      //     filename.append(".dot");
      //     snapshot(Gfull,core,iter+1,filename);
      //   }
      float sum = 0;
      for(int i=1;i<=max_core+1;i++)
	{
	  core_dist[i-1] = core_dist[i-1]/size;
	  sum += core_dist[i-1];
	  cout<<(i-1)<<" :"<<core_dist[i-1]<<"\n";
	}

      cout<<"Sum: "<<sum<<"\n";

      cout<<"Shell #\tEdges\tEdgesin\tAvg\tVar\t0\t1-5in\t1-5out\t5-10in\t5-10out\t10-20in\t10-20ou\t>20in\t>20ou\tdistr\n";
      cout<<fixed;
      for(int i=0;i<core_stat.size();i++)
  	{
  	  vector<float> stat = core_stat[i];
  	  cout<<i<<": \t";
  	  for(int j=0;j<stat.size();j++)
  	    {
  	      cout<<setprecision(1)<<stat[j]<<"\t";
  	    }
  	  cout<<"\n";
  	}
    }

  string core_stat_file = outputFile;
  core_stat_file.append("_core_stat_output.txt");

  string plot_avg_output = outputFilePrefix;
  plot_avg_output.append("_kcore_avg");

  vector<int> columns_to_plot;
  columns_to_plot.push_back(4);

  vector<string> label;
  label.push_back(outputFilePrefix);

  string xtics = "";

  for(int i=0;i<=max_core;++i)
    {
      if(i==max_core)
	{
	  stringstream ss1;
	  ss1<<"\""<<i<<"\" "<<i;
	  string next_tic = ss1.str();
	  xtics.append(next_tic);
	}
      else
	{
	  stringstream ss1;
	  ss1<<"\""<<i<<"\" "<<i<<", ";
	  string next_tic = ss1.str();
	  xtics.append(next_tic);
	}
    }
  string xstart = "0";
  string xtic = "1";
  if(max_core>180)
    xtic = "20";
  else if(max_core>90)
    xtic = "10";
  else if(max_core>18)
    xtic = "5";
    
  produce_plot(core_stat_file, plot_avg_output, columns_to_plot, label, "Average k-shell Edge Jump", "Start Shell", "Avg. shell Change",xstart,xtic,false);

  // vector<int> cores_to_plot;
  // cores_to_plot.push_back(1);
  // cores_to_plot.push_back(5);
  // cores_to_plot.push_back(10);
  // cores_to_plot.push_back(20);
  // cores_to_plot.push_back(40);
  // cores_to_plot.push_back(80);

  // vector<int> columns_to_plot;
  // vector<string> label;
  // for(int i=0;i<cores_to_plot.size();i++)
  //   {
  //     stringstream ss1;
  //     ss1<<"Core "<<cores_to_plot[i];
  //     columns_to_plot.push_back(cores_to_plot[i]+2);
  //     label.push_back(ss1.str());
  //   }

  // string xtics = "\"0\" 0, \"1-5\" 1, \"5-10\" 2, \"10-20\" 3, \">20\" 4";

  // string outfilename1 = prefix;
  // outfilename1.append("_out_plot");

  // string outfilename2 = prefix;
  // outfilename2.append("_in_plot");

  // produce_plot(filein, outfilename2,columns_to_plot,label,"Edges to k-core Interior", "Size of Jump", "Percent of Edges", xtics);

  // produce_plot(fileout,outfilename1,columns_to_plot,label,"Edges to k-core Exterior", "Size of Jump", "Percent of Edges", xtics);


    
}
