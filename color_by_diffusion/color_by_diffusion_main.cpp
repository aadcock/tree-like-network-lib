/*
Main file for calculating diffusion colorings on a graph.  Executable options:

-i input file...REQUIRED
-o output files prefix...default = "graph"
-directed if flagged, the graph is directed. Undirected otherwise
-num_steps if flagged, gives number of steps.  default is 10
-start_node if flagged, places a value of 1 at node given
-distribution file containing a row of numbers representing the starting distribution values

TODO: bad input distribution (ie, what if there are characters rather
than doubles in first line) needs to throw an error

*/

#include "../lib/graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

int main(int argc, char *argv[])
{
  string input_file = "";
  string output_file = "";
  bool directed = false;
  int num_steps = 10;
  int input_node = -1;
  string input_distribution_file = "";
  vector<double> input_distribution;
  
  //User input

  if(argc < 2)
    {
      cout<<"-i <input_file> \n-o <output_file_prefix> \n-directed <directed graph input> \n-num_steps <number_of_diffusion_steps> \n-start_node <index_of_diffusion_start_node> \n-distribution <input_distribution_file>\n";
      return(-1);
    }

  for(int i = 1; i < argc; i++)
    {
      stringstream ss;
      if (strcmp(argv[i], "-i") == 0)
	input_file = argv[i+1];
      else if (strcmp(argv[i],"-o") == 0)
	output_file = argv[i+1];
      else if (strcmp(argv[i],"-directed") == 0)
	directed = true;
      else if (strcmp(argv[i],"-start_node") == 0)
	{
	  ss<<argv[i+1];
	  ss>>input_node;
	}
      else if (strcmp(argv[i], "-num_steps") == 0)
	{
	  ss<<argv[i+1];
	  ss>>num_steps;
	}
      else if (strcmp(argv[i], "-distribution") == 0)
	input_distribution_file = argv[i+1];
    }

  ////Default values

  if (input_file.length()==0)
    {
      cerr<<"You must provide an input file\n";
      exit(EXIT_FAILURE);
    }
  cout<<"Input File: "<<input_file<<"\n";

  if (output_file.length()==0)
    {
      cerr<<"Default output 'graph_diffusion' will be used"<<"\n";
      output_file = "graph_diffusion";
    }

  if (input_distribution_file.length() == 0 && input_node < 0)
    {
      cerr<<"You must provide either a starting distribution or an input node\n";
      exit(EXIT_FAILURE);
    }


  ////////////////////////////////////
  //**********************************
  ////////////////////////////////////

  string delimiter = "\t";


  Graph G = loadGraph(input_file, delimiter, directed);
  G = connected(G);
  v_size_t size = num_vertices(G);
  input_distribution.assign(size, 0);

  if (input_distribution_file.length() == 0)
    input_distribution[input_node] = 1.0;
  else
    {
      // Reads first noncommented line of graph file and puts contents
      // in input_distribution
      ifstream graph_file;
      graph_file.open(input_distribution_file.c_str());

      if (graph_file.is_open())
	{
	  bool flag = true;
	  while (graph_file.good() && flag)
	    {
	      string line;
	      getline(graph_file, line);

	      // Skip lines beginning with #
	      if (line.substr(0,1) != "#")
		{
		  flag = false;

		  stringstream distribution_string(line);
		  size_t count = 0;
		  // Only keep the first size values
		  while (count < size && distribution_string.good())
		    {
		      distribution_string>>input_distribution[count];
		      count++;
		    }
		}
	    }
	  graph_file.close();
	}
      else
	{
	  cerr<<"Failed to open distribution file.\n";
	  exit(EXIT_FAILURE);
	}
    }
    
  vector< vector<double> > diffusions = random_walk(G, input_distribution, num_steps);
  string input_output_file = output_file;
  input_output_file.append("_");
  input_output_file.append("0.color");

  write_color_file(input_output_file, input_file, input_distribution);
  for (size_t i = 0; i < diffusions.size(); ++i)
    {
      string output_file_name = output_file;
      output_file_name.append("_");
      stringstream int_convert;
      int_convert<<i+1;
      string step;
      int_convert>>step;
      output_file_name.append(step);
      output_file_name.append(".color");

      write_color_file(output_file_name, input_file, diffusions[i]);
    }
}
