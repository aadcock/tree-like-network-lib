/*
Writes input file to dimacs format.  Assumes undirected.

By Aaron Adcock, Feb 2012, PhD Candidate, Stanford University
 */

#include "../lib/graph_lib_boost.hpp"

int main(int argc, char* argv[])
{
  string input_file_name;
  string output_file_name = "graph";
  string delimiter = "\t";

  bool input_flag  = false;
  bool output_flag = false;

  for(int i=1;i<argc;i++)
    {
      if(strcmp(argv[i], "-i")==0)
	{
	  input_file_name = argv[i+1];
	  input_flag      = true;
	}
      else if(strcmp(argv[i],"-o")==0)
	{
	  output_file_name = argv[i+1];
	  output_flag      = true;
	}

      if(strcmp(argv[i], "-d")==0)
	delimiter = argv[i+1];
    }

  if(!input_flag)
    {
      cerr<<"You must provide an input file\n";
      cout<<"Use -i <file> for input file name\nUse -o <file> for output file name\n";
      exit(EXIT_FAILURE);
    }

  if(!output_flag)
    {
      cerr<<"Warning: No output file provided, default of \"graph\" was used\n";
    }

  Graph G = loadGraph(input_file_name,delimiter,0);
  cout<<"Vertices: "<<num_vertices(G)<<" Edges: "<<num_edges(G)<<"\n";
  G       = connected(G,-1);
  
  write_dimacs(G, output_file_name);

  return(0);
}
