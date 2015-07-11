/*
Main file for performing "learning" on a graph. Basic idea, take a
graph and a some labels on the nodes.  Split nodes into a training and
test set, run the ML algorithm, produce results.  This replaces my
previous, "run graph algorithms in C, ML in Matlab framework".
Currently only uses label propagation. Later, may use as a general ML
framework.

By Aaron Adcock
Stanford University, PhD Candidate
Dec. 2013

Executable options:

-i input graph...REQUIRED
-labels vertices for which labels are known...REQUIRED
-o output file...default = "./ml.out"
-core k_min...run label prop on nodes in cores greater than k_min
-periph k_max...run label prop on nodes only in cores less than k_max
-d if y, indicates a directed graph...default = undirected graph
-s if y, suppresses terminal output...default = false
*/

#include "graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <map>
#include <limits>

// Function used to read label file
void read_label_file(string input_label_file, map<v_size_t, int> & labels, vector<bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> integer_to_label);

int main(int argc, char *argv[])
{
  string input_graph = "";
  string input_label_file = "";
  string output_file = "./label_prop.out";
  int k_min = 0; 
  int k_max = numeric_limits<int>::max();
  bool directed = false;
  bool suppress_output = false;
  int max_iterations = 10;
  //User input

  if(argc < 2)
    {
      cout<<"-i <input_graph_file> \n-labels <input_labels_file> \n-o <output_file> \n-max_iter <max_iterations> \n-core <min_core_num> \n -periph <max_core_num> \n-d <y=directed graph input> \n-s <y=suppress terminal output>\n";
      return(-1);
    }

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if (strcmp(argv[i], "-i") == 0)
	input_graph = argv[i + 1];
      else if (strcmp(argv[i], "-labels") == 0)
	input_label_file = argv[i + 1];
      else if (strcmp(argv[i],"-o") == 0)
	output_file = argv[i + 1];
      else if (strcmp(argv[i], "-max_iter") == 0)
	max_iterations = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-core") == 0)
	k_min = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-periph") == 0)
	k_max = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-d") == 0)
	{	
	  if (strcmp(argv[i+1], "y") == 0)
	    directed = true;
	}
      else if (strcmp(argv[i], "-s") == 0)
	{	  
	  if (strcmp(argv[i+1], "y") == 0)
	    suppress_output = true;
	}

    }

  ////Default values

  if (input_graph.length()==0)
    {
      cerr<<"You must provide an input graph file\n";
      exit(EXIT_FAILURE);
    }
  else if (input_label_file.length()==0)
    {
      cerr<<"You must provide an input labels file\n";
      exit(EXIT_FAILURE);
    }
  else if (!suppress_output)
    cout<<"Input File: "<<input_graph<<"\n";

  // Get graph, edge direction reversed will get connected component
  // later as we want to make sure labels are applied beforehand
  Graph G = loadGraph(input_graph);
  cout<<"Graph loaded\n";

  map<v_size_t, int> labels;
  vector<bool> is_labeled(num_vertices(G), false);
  // I will assume that the labels come as ints, but not necessarily
  // ints that begin at zero.  Thus two converters.
  map<int,int> label_to_integer;
  map<int,int> integer_to_label;

  read_label_file(input_label_file, labels, is_labeled, label_to_integer, integer_to_label);

  size_t num_labels = label_to_integer.size();
  cout<<"Labels loaded\n";

  // Need to get train/test set Should I do this and make sure all
  // labels are represented?  This shouldn't be an issue if each label
  // sufficiently represented in the data...I'll ignore this for now
  
}
  
/* 
 Reads in labels from file.  Currently does not check vertex indices
 to see if they are valid.  Will throw error is they are larger than
 the size of is_labeled.  
*/
void read_label_file(string input_label_file, map<v_size_t, int> & labels, vector<bool> & is_labeled, map<int, int> & label_to_integer, map<int, int> integer_to_label)
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
