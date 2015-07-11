/*
This main function produces a file containing the k_core information
of the graph in the form of a list of core numbers for each vertex.
The order is determined by the index of the node and is the same as
the indices produced by the write dimacs program.  Thus, this k-core
file is a viable color file for the INDDGO tree decomposition
software.

It can also be used to take an existing color file of the giant
component of a network (note that if the color file is not for the
giant component, this code will terminate) and produce a color file
for the k-periphery of the network.

COMMAND LINE INPUT OPTIONS:

-i input_file
-color color_file
-write_core 
-write_degree
-periph k-periphery to produce color file for
-o output_file

By Aaron Adcock, PhD Candidate, Stanford University, 2012
 */

#include "../lib/graph_lib_boost.hpp"

Graph get_periphery_colors(Graph & G, vector<double> & color_vector, vector<int> & core_vector, vector<double> & periphery_vector, int threshold);

int main(int argc, char *argv[])
{

  // Input variables. Empty strings indicate no input. If color_file
  // is provided, then core_to_write indicates the core to write out a
  // color file for or a nonpositive core_to_write value indicates
  // writing all core color files. 
  string input_filename  = "";
  string output_filename = "";
  string color_filename  = "";
  int periphery_to_write  = -1; 
  bool write_core    = false;
  bool write_degree  = false;

  // Parse command line inputs
  for(int i = 1; i < argc; ++i)
    {
      if (strcmp(argv[i], "-i") == 0)
	input_filename = argv[i+1];
      else if (strcmp(argv[i], "-o") == 0)
	output_filename = argv[i+1];
      else if (strcmp(argv[i], "-color") == 0)
	color_filename = argv[i+1];
      else if (strcmp(argv[i], "-write_core") == 0)
	write_core = true;
      else if (strcmp(argv[i], "-write_degree") == 0)
	write_degree = true;
      else if (strcmp(argv[i], "-periph") == 0)
	periphery_to_write = atoi(argv[i+1]);
    }

  if(input_filename.length()==0)
    {
      cerr<<"You must provide an input file! All other commands are optional\n";
      cout<<"Use -i <file> for input file name\n"
	    "Use -color <file> for the input color file name\n"
	    "Use -write_core to write out core numbers as the colors\n"
	    "Use -write_degree to write out degree numbers as the colors\n"
	    "Use -periph to indicate the k-periphery, if no value is provided "
	    "then the color files for each periphery will be printed if no color "
	    "file is provided, then no k-periphery color files will be printed\n"
	    "Use -o <file> for output file name prefix\n";
      exit(EXIT_FAILURE);
    }

  // Load network from input file and display basic info
  Graph G = loadGraph(input_filename,"\t",0);
  cout<<"Size of G: "<<num_vertices(G)<<"\n";

  // Calculate giant component and display basic info
  G = connected(G,-1);

  cout<<"Size of connected G: "<<num_vertices(G)<<"\n";
  
  //Calculate k core of G
  vector<int> core_vector = k_core(G);
  vector<int> degree_vector;

  cout<<"k-core calculated\n";

  graph_traits<Graph>::vertex_iterator vit, vite;

  for(tie(vit,vite)=vertices(G);vit!=vite;++vit)
    {
      Vert v  = *vit;
      int deg = out_degree(v,G);
      degree_vector.push_back(deg);
    }
  
  cout<<"Degrees calculated\n";

  string output_filename_k,output_filename_d;

  if(output_filename.length() == 0)
    {
      output_filename_k = input_filename;
      size_t pos  = output_filename_k.find_last_of(".");
      output_filename_k = output_filename_k.substr(0,pos);
      output_filename_k.append("_k_core.color");

      output_filename_d = input_filename;
      pos  = output_filename_d.find_last_of(".");
      output_filename_d = output_filename_d.substr(0,pos);
      output_filename_d.append("_deg.color");
    }
  else
    {
      output_filename_k = output_filename;
      output_filename_k.append("_k_core.color");
      
      output_filename_d = output_filename;
      output_filename_d.append("_deg.color");
    }

  // Makes color file
  if (write_core)
    write_color_file(output_filename_k, input_filename, core_vector);
  if (write_degree)
    write_color_file(output_filename_d, input_filename, degree_vector);

  // Uses independent color file for whole network to produce a color
  // file for just a k-periphery of the network.  The k-periphery is
  // defined as all nodes with a core number less than or equal to
  // k. If k is greater than or equal to k_max for the network, then
  // this is the same as writing the k-core for a giant component of
  // the network.
  if (color_filename.length() != 0)
    {
      vector<double> color_vector;
      ifstream color_file;
      string line;

      color_file.open(color_filename.c_str());
      getline(color_file, line);

      while (line[0] == '#')
	getline(color_file, line);

      cout<<"Color file opened\n";
      stringstream linestream(line);
      double num;
      while (linestream >> num)
	color_vector.push_back(num);
      
      cout<<"Color vector created\n";
      if (periphery_to_write > 0)
	{
	  vector<double> periphery_color_vector;
	  get_periphery_colors(G, color_vector, core_vector, periphery_color_vector, periphery_to_write);
	  
	  string output_filename_color = color_filename;
	  size_t pos = output_filename_color.find_last_of(".");
	  output_filename_color = output_filename_color.substr(0,pos);
	  stringstream temp_stream;

	  temp_stream << output_filename_color << "_periphery_" << periphery_to_write << ".color";

	  output_filename_color = temp_stream.str();

	  write_color_file(output_filename_color, input_filename, periphery_color_vector);
	}
      else
	{
	  int max_core = -1;
	  for (size_t i = 0; i < core_vector.size(); ++i)
	    if (max_core < core_vector[i])
	      max_core = core_vector[i];
         
	  cout<<"Max Core calculated\n";
	  vector<double> periphery_color_vector;
	  for (int i = max_core; i > 0; --i)
	    {
	      cout<< "Periphery " << i << "\n";
	      Graph G_periph = get_periphery_colors(G, color_vector, core_vector, periphery_color_vector, i);
	      
	      cout<<"Calculated periphery colors\n";
	      string output_filename_color = color_filename;
	      size_t pos = output_filename_color.find_last_of(".");
	      output_filename_color = output_filename_color.substr(0,pos);
	      stringstream color_stream;

	      color_stream << output_filename_color << "_periphery_" << i << ".color";
	      output_filename_color = color_stream.str();

	      cout<< "Color file produced is: " << output_filename_color << "\n";
	      cout<< "Size of periphery color vector: "<<periphery_color_vector.size()<<"\n";
	      write_color_file(output_filename_color, input_filename, periphery_color_vector);
	      cout<< "Color file written\n";
	      
	      string output_dimacs = input_filename;
	      pos = output_dimacs.find_last_of(".");
	      output_dimacs = output_dimacs.substr(0,pos);
	      stringstream output_stream;
	      output_stream << output_dimacs << "_periphery_" << i << ".dimacs";
	      output_dimacs = output_stream.str();

	      write_dimacs(G_periph, output_dimacs);
	    }
	}
    }

  return(0);
}

/*
  This function is a helper function which takes an input color
  vector, a vector for storing the output color vector, a core vector,
  and a threshold for keeping a node in the periphery color vector.
  For example, if the threshold is set to 2, then the periphery color
  vector will contain all colors connected to a node with a core
  number less than or equal to 2.

  INPUT:
  color_vector  The input color vector (should be same length as core_vector)
  core_vector   The input core vector
  periphery_vector The output vector containing only colors above the threshold
  threshold A double representing the color vector cutoff
*/

Graph get_periphery_colors(Graph & G, vector<double> & color_vector, vector<int> & core_vector, vector<double> & periphery_vector, int threshold)
{
  periphery_vector.clear();
  vector<Vert> periphery_vertices;

  // Boost graph library traits
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  graph_traits<Graph>::vertex_iterator vit, vitend;

  for (tie(vit, vitend) = vertices(G); vit != vitend; ++vit)
    {
      Vert v = *vit;
      G[v].color = color_vector[index[v]];
      if (G[v].core_num <= threshold)
	periphery_vertices.push_back(v);
    }
  
  // This assumes that index order is preserved at all steps.
  Graph G_periph = subset(G, periphery_vertices);
  G_periph = connected(G_periph);

  for (tie(vit, vitend) = vertices(G_periph); vit != vitend; ++vit)
    {
      Vert v = *vit;
      periphery_vector.push_back(G_periph[v].color);
    }
  
  return G_periph;
}
