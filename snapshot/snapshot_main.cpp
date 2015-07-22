/*
This file generates an executable that calculates the Network
Community Profile plot of the input graph.  The intended use is for
shell scripts in calculating various statistics of graphs.
  
By Aaron Adcock, PhD candidate at Stanford University
Oct. 2011 (modified May 2014)
*/

#include "../lib/graph_lib_boost.hpp"
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

void read_label_file(string input_label_file, map<v_size_t, int> & labels, map<v_size_t, bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> & integer_to_label);

// Main function for driving snapshots
int main(int argc, char *argv[])
{
  //Variable declarations for user input
  string input_file;
  string output_file_prefix = "graph";
  string label_file = "";
  bool directed = true;
  bool suppressOutput = false;
  double diffusion_size = .001;
  double alpha  = .001; 
 
  //User input
  if(argc < 2)
    {
      cout<<"-i <input_file> \n-o <output_file_prefix> \n-labels <label file, if blank uses k_core>\n-a <alpha_value for PageRank, currently not in use> \n-e <eps, error value for PageRank, currently not in use> \n-d <y=directed graph input> \n-s <y=suppress terminal output> \n";
      return(-1);
    }

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	{
	  ss<<argv[i+1];
	  ss>>input_file;
	}
      else if(strcmp(argv[i], "-o")==0)
	{
	  ss<<argv[i+1];
	  ss>>output_file_prefix;
	}	
      else if(strcmp(argv[i], "-labels") == 0)
	{
	  ss<<argv[i + 1];
	  ss>>label_file;
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

  if (input_file.length() == 0)
    {
      cout<<"You must provide an input file\n";
      return(0);
    }
  else
    cout<<input_file<<"\n";

  if (alpha <= 0)
    {
      cout<<"Alpha must be greater than 0, using default .001\n";
      alpha = .001;
    }

  if (diffusion_size <= 0)
    {
      cout<<"Diffusion size must be greater than 0, using default .005\n";
      diffusion_size = .005;
    }

  // Load graph
  Graph G       = loadGraph(input_file, "\t");
  G             = connected(G,-1);
  v_size_t size = num_vertices(G);
  
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);


  if(!suppressOutput)
    cout<<"Size of connected component of G "<<size<<"\n";

  //Calculate k-core, note that the core # is also stored as a property on each node

  map<v_size_t, int> labels;
  map<v_size_t, bool>  is_labeled;
  map<int, int> label_to_integer;
  map<int, int> integer_to_label;

  if (label_file.length() == 0)
    {
      vector<int> core = k_core(G);
	
      graph_traits<Graph>::vertex_iterator vi, vie;
      for (tie(vi, vie) = vertices(G); vi != vie; ++vi)
	{
	  labels[index[*vi]] = core[index[*vi]];
	  integer_to_label[core[index[*vi]]] = core[index[*vi]];
	}
    }
  else
    read_label_file(label_file, labels, is_labeled, label_to_integer, integer_to_label);

  // // Find maximum label value
  // int max_label = 0;  
  // for (map_it = labels.begin(); map_it != labels.end(); ++map_it)
  //   {
  //     pair<v_size_t, int> item = *map_it;
  //     if (max_label < item.second)
  // 	max_label = item.second;
  //   }
  cout<<"Finished getting labels.\n";
  map<v_size_t, int>::iterator map_it;

  // Produce snapshots of each label
  map<int, bool> label_done;
  for (map_it = labels.begin(); map_it != labels.end(); ++map_it)  
    {
      if (label_done.find(map_it->second) == label_done.end())
	{
	  map<int, double> distribution = community_snapshot(G, labels, map_it->second, alpha, diffusion_size);
	  label_done[map_it->second] = true;
	
	  cout<<"Processing "<<integer_to_label[map_it->second]<<"-label:\n";
	  cout.flush();
	  double total = 0;
	  map<int, double>::iterator map_jt;
	  for (map_jt = distribution.begin(); map_jt != distribution.end(); ++map_jt)
	    total += map_jt->second;

	  if (total != 0)
	    {
	      for (map_jt = distribution.begin(); map_jt != distribution.end(); ++map_jt)
		cout<<integer_to_label[map_jt->first]<<": "<<map_jt->second / total <<"\n";
	    }
	  else
	    cout<<"No nodes in this label.\n";
	}
    }
}


/* 
 Reads in labels from file.  Currently does not check vertex indices
 to see if they are valid.  Will throw error is they are larger than
 the size of is_labeled.  Assumes file of following format:

 map meanings:
 labels: a map from vertex_id to an integer, sequential label
 label_to_integer: a map from the original integer label to sequential integer label
 integer_to_label: a map from the sequential integer label to the original integer label

 vertex_id vertex_label
 vertex_id vertex_label
 ...
*/
void read_label_file(string input_label_file, map<v_size_t, int> & labels, map<v_size_t, bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> & integer_to_label)
{
    ifstream label_file(input_label_file.c_str());

    if (label_file.is_open())
    {
      string line;
      int line_number = 0;
      int labels_so_far = 0;
      while (getline(label_file, line))
	{
	  line_number++;
	  //	  cout<<line_number<<"\n";
	  // Comment character is #
	  if (line.substr(0,1) != "#")
	    {
	      istringstream line_stream(line);
	      int vertex_id, vertex_label;

	      if (line_stream>>vertex_id>>vertex_label)
		{
		  is_labeled[vertex_id] = true;
		  if (label_to_integer.find(vertex_label) == label_to_integer.end())
		    {
		      label_to_integer[vertex_label] = labels_so_far;
		      integer_to_label[labels_so_far] = vertex_label;
		      labels[vertex_id] = label_to_integer[vertex_label];
		      labels_so_far++;
		    }
		  else
		    labels[vertex_id] = label_to_integer[vertex_label]; 
		}
	    }
	}
    }
}
