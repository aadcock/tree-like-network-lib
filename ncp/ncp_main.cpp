/*
  This file generates an executable that calculates the Network Community
  Profile plot of the input graph.  The intended use is for shell scripts
  in calculating various statistics of graphs.
  
  By Aaron Adcock, PhD candidate at Stanford University
  Oct. 2011 
*/

#include "graph_lib_boost.hpp"
#include "tree_lib_boost.hpp"
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
#include <sstream>
#include <string.h>

void write_data_file (string outputFileName, vector<size_t> data) 
{
  ofstream output;
  output.open(outputFileName.c_str());

  for(v_size_t i = 0; i < data.size(); i++)
    output<<(i+1)<<"\t"<<data[i]<<"\n";
  
  output.close();
}

int main(int argc, char *argv[])
{
  string inputFile;
  string outputDirectory  = "./";
  string outputFilePrefix = "graph";
  bool directed = true;
  bool suppressOutput = false;
  v_size_t maxC = 0;
  int step_size = 10;
  double alpha  = .001;
  bool produce_graph_viz = false;
  bool tree_provided = false;
  bool tree_calc = false;
  string tree_file = "";

  if(argc<2)
    {
      cout<<"The options are (-i and -n and -o are required):\n-i <input-file>\n-n <name_prefix> \n-o <output_directory> \n-alpha <alpha value, default .001> \n-maxcomm <max community size> \n-step <target community step size, default 10>\n-produce_graph_viz\n-tree_file <tree_decomp_file>\n-calc_tree_community <calculates communities using TD> \n";
    }
  

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
      else if(strcmp(argv[i],"-alpha")==0)
	alpha = atoi(argv[i+1]);
      else if(strcmp(argv[i],"-maxcomm")==0)
	maxC = atoi(argv[i+1]);
      else if(strcmp(argv[i],"-step")==0)
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
      else if (strcmp(argv[i], "-produce_graph_viz") == 0)
	produce_graph_viz = true;
      else if (strcmp(argv[i], "-tree_file") == 0)
	{
	  tree_provided = true;
	  ss << argv[i+1];
	  ss >> tree_file;
	}
      else if (strcmp(argv[i], "-calc_tree_community") == 0)
	tree_calc = true;
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
 cout.flush();
  if(step_size <=0)
    {
      cout<<"Step size must be greater than 0, using default 10\n";
      step_size = 10;
    }
  cout.flush(); 
  cout<<"directory\n";
  cout.flush();
  int l = outputDirectory.length();
  int c = outputDirectory.compare(l-1,1,"/");
  if(c != 0)
    outputDirectory.append("/");

  ////
  Graph G = loadGraph(inputFile,"\t",false);
  G       = connected(G);
  v_size_t size = num_vertices(G);

  if(maxC <= 0)
    {
      cout<<"Maximum Community Size must be greater than 0, using default graph size/2\n";
      maxC = size/2;
    }

  if(!suppressOutput)
    cout<<"Size of connected component of G "<<size<<"\n";

  cout.flush();
  vector<double> community_Conductance (maxC, 1.0);
  vector<Vert> dummy;
  vector< vector<Vert> > best_communities (maxC, dummy);
  if (produce_graph_viz)
    {
      if (tree_provided && tree_calc)
	{
	  TreeDecomp td = loadTreeDecomp(G, tree_file);
	  cout<<community_Conductance.size()<<" "<<best_communities.size()<<"\n";
	  community_Conductance = td_bag_communities(G, td, maxC, community_Conductance, best_communities, true);
	  community_Conductance = td_eccentricity_communities(G, td, maxC, community_Conductance, best_communities, true);
	  
	}
      else
	community_Conductance = ncp_calc(G, maxC, step_size, alpha, best_communities, true);
      
      vector <int> temp = k_core(G);
      vector <double> color;
      color.assign(temp.begin(), temp.end());
      for (int i = 0; i < best_communities.size(); ++i)
	{
	  Graph H = subset(G, best_communities[i]);
	  // Need to produce color vector specific to this subset.
	  vector<int> square_nodes;
	  string output_graph_viz = outputDirectory;
	  output_graph_viz.append(outputFilePrefix);
	  stringstream ss;
	  string num;
	  ss<<i+1;
	  ss>>num;
	  output_graph_viz.append("_community_size_");
	  output_graph_viz.append(num);
	  string output_color_file = output_graph_viz;
	  output_graph_viz.append(".dot");
	  write_graphviz(H, output_graph_viz, color, square_nodes, false);

	  vector<double> in_community (size, 0);
	  for (int j = 0; j < best_communities[i].size(); ++j)
	    in_community[best_communities[i][j]] = 1.0;

	  output_color_file.append(".color");
	  write_color_file(output_color_file, inputFile, in_community);
	}

      if (tree_provided)
	{
	  TreeDecomp td = loadTreeDecomp(G, tree_file);
	  vector<size_t> num_bags, min_ecc, max_card, med_card, surface_area;
	  
	  compute_community_tree_stats(best_communities, G, td, num_bags, min_ecc, max_card, med_card, surface_area);
	  
	  string output_file_long_prefix = outputDirectory;
	  output_file_long_prefix.append(outputFilePrefix);

	  string out_num_bags = output_file_long_prefix;
	  out_num_bags.append("_num_bags.data");
	  string out_min_ecc = output_file_long_prefix;
	  out_min_ecc.append("_min_ecc.data");
	  string out_max_card = output_file_long_prefix;
	  out_max_card.append("_max_card.data");
	  string out_med_card = output_file_long_prefix;
	  out_med_card.append("_med_card.data");
	  string out_surface_area = output_file_long_prefix;
	  out_surface_area.append("_surface_area.data");

	  write_data_file(out_num_bags, num_bags);
	  write_data_file(out_min_ecc, min_ecc);
	  write_data_file(out_max_card, max_card);
	  write_data_file(out_med_card, med_card);
	  write_data_file(out_surface_area, surface_area);
	}
    }
  else
    {
      if (tree_provided && tree_calc)
	{
	  TreeDecomp td = loadTreeDecomp(G, tree_file);
	  community_Conductance = td_bag_communities(G, td, maxC, community_Conductance, true);
	  community_Conductance = td_eccentricity_communities(G, td, maxC, community_Conductance, true);	  
	}
      else
	community_Conductance = ncp_calc(G, maxC, step_size, alpha, true);
    }
  cout<<"Conductance Computed\n";

  ofstream outputNCP;
  string outputFileName = outputDirectory;
  outputFileName.append(outputFilePrefix);
  outputFileName.append("_ncp.txt");
  outputNCP.open(outputFileName.c_str());

  for(v_size_t i=0; i<community_Conductance.size();i++)
      outputNCP<<(i+1)<<"\t"<<community_Conductance[i]<<"\n";
  
  outputNCP.close();

  string outputPlotFile = outputDirectory;
  outputPlotFile.append(outputFilePrefix);
  outputPlotFile.append("_ncp_plot");

  string label = outputFilePrefix;
  label.append(" ncp plot");

  string title = outputFilePrefix;
  title.append(" NCP Plot");

  vector<int> int_to_vec;
  int_to_vec.push_back(2);

  vector<string> string_to_vec;
  string_to_vec.push_back(label);

  string outputPNG = "./";
  outputPNG.append(outputFilePrefix);
  outputPNG.append("_ncp_plot");

  //  produce_loglog_plot(outputFileName,outputPlotFile,outputPNG,int_to_vec,string_to_vec,title,"Community size", "Conductance","",suppressOutput);

  
}


