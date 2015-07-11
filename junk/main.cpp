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
  
  string directory;
  cout<<"Location of Data file: ";
  getline(cin,directory);

  string prefix;
  cout<<"Data File Prefix (.txt will be appended)?: ";
  getline(cin,prefix);

  string strDirected;
  bool directed;
  cout<<"Is graph file in directed format? (y or n): ";
  getline(cin,strDirected);
  if(strDirected.compare(0,1,"y")==0 || strDirected.compare(0,1,"Y")==0)
    directed = true;
  else
    directed = false;

  string strSuppress;
  bool suppressOutput;
  cout<<"Suppress terminal output? (y or n): ";
  getline(cin,strDirected);
  if(strSuppress.compare(0,1,"y")==0 || strSuppress.compare(0,1,"Y")==0)
    suppressOutput = true;
  else
    suppressOutput = false;

  string outputDirectory = "../Graphs/Results/";
  //  cout<<"Output Directory: ";
  // getline(cin,outputDirectory

  //Load the graph from file
 
  string direct = outputDirectory;
  direct.append(prefix);
  direct.append("_results/");

  int directory_made =  mkdir(direct.c_str(),S_IRWXU);

  string fileUse=directory;
  fileUse.append(prefix);
  fileUse.append(".txt");

  Graph  G = loadGraph(fileUse,"\t",directed);

  int size = num_vertices(G);

  cout<<"Size of G "<<size<<"\n";

  //Size of connected component
  
  Graph Gfull = G;
  G = connected(G,-1);
  cout<<"Size of Gconn "<<num_vertices(G)<<"\n";

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  //Calculating the ncp plot (see 'Empirical Comparison of Algorithms for Network Community detection
  //Current implementation (as of July 2011) uses Spectral sweep cuts with random seeds to find best community structure at different size scales
  //Page Rank/ncp Parameters  
  //alpha gives teleportation probability
  //by setting one component of s to 1, you choose the seed node
  //step gives the starting step size that the ncp algorithm works with

  vector<double> s (size,0);
  vector<double> p;
  const double alpha = .001;
  const int step = 10;

  vector<int> core = k_core(G);
  for(int i=0; i<core.size();i++)
    cout<<core[i]<<"\t";

  cout<<"\n";

  string ncpName = direct;
  ncpName.append("ncp_plot/");
  directory_made =  mkdir(ncpName.c_str(),S_IRWXU);

  ncpName.append(prefix);
  
  //Calculate iterated ncp
  //k_core_ncp(G,ncpName,core,alpha,step);
  
  ncpName.append("_kcore_ncp.txt");

  vector<int> ncpColumn;
  vector<string> ncpData;

  for(int i=1;i<11;i++)
    {
      string data = prefix;
      stringstream ns;

      data.append(" core: ");
      ns<<i;
      data.append(ns.str());
      
      ncpColumn.push_back(i+1);
      
      ncpData.push_back(data);
    }  
  ////  


  int max_core = 0;
 
  for(int i=0; i<core.size();i++)
    {      
      if(core[i]>max_core)
	max_core = core[i];
    }
  
  vector<vector<int> > core_groups(max_core+1);
  for(int i=0;i<core.size();i++)
      core_groups[core[i]].push_back(i);
  
  srand(time(NULL));

  vector< vector<double> > snapshot_dist;

  for(int i=0;i<=max_core; i++)
    {
      vector<int> shell = core_groups[i];

      int size_shell = shell.size();
      int num_iter   = size_shell/10;

      if(size_shell > 1 && size_shell < 20)
	num_iter = 2;
      
      if(size_shell == 1)
	num_iter = 1;
      
      vector<bool> check(size_shell,false);

      stringstream ss;
      ss<<i;

      string snapshotDir = direct;
      snapshotDir.append("snapshot/");
      directory_made =  mkdir(snapshotDir.c_str(),S_IRWXU);
      snapshotDir.append("core_");
      snapshotDir.append(ss.str());
      directory_made =  mkdir(snapshotDir.c_str(),S_IRWXU);

      vector< double > spectral_core_dist(max_core+1,0);
      int sum = 0;

      for(int j=0;j<num_iter;j++)
	{
	  
	  double diffusion_size = .01;
	  vector<Vert> snapshot;
	  int r = rand() % size_shell;
	  stringstream ss2;

	  while(check[r])
	    r = (r+1) % size_shell;
	  
	  check[r] = true;

	  Vert rv = vertex(shell[r],G);
	  string snapshotFile = snapshotDir;

	  ss2<<j;
	  snapshotFile.append("/");
	  snapshotFile.append(prefix);
	  snapshotFile.append("_spectral_snapshot");
	  snapshotFile.append(ss2.str());
	  
	  string snapshotGviz = snapshotFile;
	  snapshotGviz.append("_gviz.txt");
	  snapshotFile.append(".txt");


	  spectral_snapshot(G,snapshot,rv,alpha,diffusion_size);

	  for(int k=0;k<snapshot.size();k++)
	    {
	      Vert sv = snapshot[k];
	      int ind = index[sv];

	      spectral_core_dist[(core[ind])]++;
	      sum++;
	    }

	  Graph G_s = subset(G,snapshot);

	  adjacency_list<vecS,vecS,undirectedS> localUGraph;

	  copy_graph(G_s,localUGraph);
	  
	  ofstream outViz;

	  outViz.open(snapshotGviz.c_str());
	  write_graphviz(outViz,localUGraph);

	  outViz.close();
	}

      if(sum!=0)
	{
	  for(int j=0;j<spectral_core_dist.size();j++)
	    spectral_core_dist[j] /= sum;
	}
      snapshot_dist.push_back(spectral_core_dist);
 
    }
  string distOut = direct;
  distOut.append(prefix);
  distOut.append("_spectral_core_dist.txt");
      
  ofstream outDistFile;
  outDistFile.open(distOut.c_str());
  vector<int> columns_to_plot;
  vector<string> label;

  for(int i=0;i<=max_core;i++)
    {
      outDistFile<<i<<"\t";
      for(int j=0;j<snapshot_dist.size();j++)
	{
	  vector<double> C;
	  C.assign(snapshot_dist[j].begin(),snapshot_dist[j].end());
	  if(C.size()>i)
	    outDistFile<<C[i]<<"\t";
	  else
	    outDistFile<<"\t";         //C[C.size()-1]<<"\t";
	}
      outDistFile<<"\n";
      columns_to_plot.push_back((i+2));
      stringstream ss;
      ss<<i;
      string column_label = "Shell: ";
      column_label.append(ss.str());
      label.push_back(column_label);
    }
  outDistFile.close();

  string distOutPlot = direct;
  distOutPlot.append(prefix);
  distOutPlot.append("_spec_core_dist_plot");
  produce_plot(distOut, distOutPlot, columns_to_plot, label, "Random Seed Diffusion Distribution using k-shell decomposition", "in k-Shell #", "\% of Nodes","");




	////
  vector<vector<float> > core_stat;
  string core_stat_file = direct;
  core_stat_file.append("/");
  core_stat_file.append(prefix);

  k_core_stat(core_stat,Gfull,core_stat_file);
  
  if(!suppressOutput)
    {
      vector<int> core_dist(max_core+1,0);
  
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

      for(int i=1;i<=max_core;i++)
  	cout<<i<<" :"<<core_dist[i-1]<<"\n";

      cout<<"Core #\tEdges\tEdgesin\tAvg\tVar\t0\t1-5in\t1-5out\t5-10in\t5-10out\t10-20in\t10-20ou\t>20in\t>20ou\n";
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
