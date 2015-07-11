/*
This makes an executable to call the write_graphviz function.

Produces a graph_viz .dot file for laying out the network.  The color
vector is a set of numbers which produce a heat map of colors from
blue (0, or min value) to red(1, max value)

 By Aaron Adcock, May, 2013 PhD Candidate at Stanford University
 */

#include "../lib/graph_lib_boost.hpp"
#include <boost/graph/copy.hpp>

int main(int argc, char *argv[])
{
   
  if(argc < 2)
    {
      cout<<"-i <input_file> \n-o <output_file_prefix> \n-d <y=directed graph output> \n-color <color_file> \n-index <index_file> <-- This function just displays the corresponding indices and colors from a file for reference\n-square <vertex to be highlighted>\n-labeled <node_index_displayed>\n";
      return(-1);
    }

  string input_file;
  string output_file = "graph_gviz.dot";
  string color_file;
  string index_file;

  vector<int> square_nodes;
  bool directed = false;
  bool labeled = false;

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	input_file = argv[i+1];
      else if(strcmp(argv[i],"-o")==0)
	output_file = argv[i+1];
      else if(strcmp(argv[i],"-d")==0)
	{	
	  if(strcmp(argv[i+1], "y")==0)
	    directed = true;
	}
      else if(strcmp(argv[i],"-color")==0)
	color_file = argv[i+1];
      else if(strcmp(argv[i],"-index")==0)
	index_file = argv[i+1];
      else if(strcmp(argv[i],"-square")==0)
	{
	  istringstream ss(argv[i+1]);
	  int val;
	  ss>>val;
	  square_nodes.push_back(val);
	}
      else if (strcmp(argv[i], "-labeled") == 0)
	labeled = true;
    }

  sort(square_nodes.begin(),square_nodes.end());

  if(input_file.length()==0)
    {
      cerr<<"You must provide an input file\n";
      exit(EXIT_FAILURE);
    }

  if(index_file.length()!=0)
    {
      cout<<"Indices in provided index file: \n";
      vector<double> temp; 
      load_color_file(index_file, temp);
      vector<int> index;

      for(size_t i = 0; i < temp.size(); ++i)
	{
	  index.push_back((int) temp[i]);
	  cout<<i<<" "<<index[i]<<"\n";
	}
      cout<<"Program done.\n";

      return 0;
    }

  vector<double> color;
  load_color_file(color_file, color);
      
  Graph orig_G = loadGraph(input_file);

  v_size_t n = num_vertices(orig_G);

  cout<<"Graph size: "<<n<<"\n";
  cout<<"Graph edges: "<<num_edges(orig_G)<<"\n";

  if(n!=color.size())
    {
      cout<<"Vector size: "<<color.size()<<"\n";
      cerr<<"Color vector not sized properly. replacing with min value\n";
      color.assign(n,1);
    }

  vector<size_t> sizes;
  sizes.reserve(n);

  sizes.assign(color.begin(), color.end());

  if (labeled)
    {
      double max = color[0];
      for (size_t i = 0; i < color.size(); ++i)
	{
	  if (max < color[i])
	    max = color[i];
	}

      for (size_t i = 0; i < color.size(); ++i)
	color[i] = color[i] / max;

      write_scaled_labeled_graphviz(orig_G, output_file, color, sizes, directed);
    }
  else
    write_graphviz(orig_G, output_file, color, square_nodes, directed);
  
  return 0;
}





