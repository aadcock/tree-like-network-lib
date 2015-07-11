/*
This function takes in a tree decomposition, written in the
tree-decomposition dimacs format used in the INDDGO software, for the
purpose of writing out the number of bags containing nodes with
different values computed on the nodes.

By Aaron Adcock, Stanford University, 
PhD Candidate, 2013
 */

#include "graph_lib_boost.hpp"

//Loads all values in a file into a vector.  In this main, we assume
//that these numbers are ordered by the same vertex index as used in
//constructing the tree decomposition
//vector<int> load_color_file(string & color_file);

int main(int argc, char* argv[])
{

  //Command line input
  string input_file_name = "";
  string output_file_name     = "stats.txt";
  
  vector<string> color_file_names;

  for(int i=1;i<argc;++i)
    {
      if(strcmp(argv[i], "-i")==0)
	input_file_name = argv[i+1];
      else if (strcmp(argv[i], "-o")==0)
	output_file_name = argv[i+1];
      else if (strcmp(argv[i], "-color")==0)
	{
	  string temp; 
	  temp = argv[i+1];
	  color_file_names.push_back(temp);
	}
    }

  if(input_file_name=="")
    {
      cerr<<"You must provide an input file\n";
      cout<<"Use -i <file> for input file name, -o <file> for output file name, and -color <file> for color files\n";
      exit(EXIT_FAILURE);
    }
  
  //load tree decomposition from file
  Graph G = loadDimacsTreeDecomp(input_file_name);
  cout<<"tree loaded\n";
  vector < vector<int> > color_values;

  bool first = true;
  int n=0;
  
  //load color files
  for(int i=0;i<color_file_names.size();++i)
    {
      vector<int> temp = load_color_file(color_file_names[i]);

      if(first)	
	{
	  n = temp.size();	
	  first = false;
	}
      else
	{
	  if(n!=temp.size())
	    {
	      cerr<<"Color vectors of different length, exiting "<<n<<" "<<temp.size()<<"\n";
	      exit(EXIT_FAILURE);
	    }
	}

      color_values.push_back(temp);
    }

  //  v_size_t tree_size = num_vertices(G);

  graph_traits<Graph>::vertex_iterator vit, vitend;

  vector<int> num_bags;

  //Fine size of underlying network
  if(n!=0)
    {
      cout<<"Network size: "<<n<<"\n";
      for(int i=0;i<n;++i)
	num_bags.push_back(0);
//num_bags.assign(n,0);
    }
  else
    {
      cout<<"Determining size of network\n";
      unsigned long max = 0;
      for(tie(vit,vitend)=vertices(G);vit!=vitend;++vit)
	{
	  Vert v = *vit;

	  vector<unsigned long> bag = G[v].bag;

	  for(int i=0;i<bag.size();++i)
	    {
	      if(bag[i]>=max)
		max = bag[i];
	    }
	}
      n=max;
      num_bags.assign(max,0);
      cout<<"Network size: "<<max<<"\n";
    }


  //  cout<<"here?\n";

  for(tie(vit,vitend)=vertices(G);vit!=vitend;++vit)
    {
      Vert v = *vit;
      vector<unsigned long> bag = G[v].bag;

      for(int i=0;i<bag.size();++i)
	{
	  //	  cout<<i<<" "<<bag[i]<<" "<<num_bags.size()<<"\n";
	  int id = bag[i]-1;
	  ++num_bags[id];
	}
    }

  ofstream output_file;
  output_file.open(output_file_name.c_str());

  output_file<<"#Counts number of bags in "<<input_file_name<<" containing each vertex.  These values are in the last and first column.  Additional color files (originally containing graph indices and/or k-core number) may be in the other columns.\n";
  
  for(int i=0;i<num_bags.size();++i)
    {
      
      output_file<<i<<"\t";

      for(int j=0;j<color_values.size();++j)
	output_file<<color_values[j][i]<<"\t";

      output_file<<num_bags[i]<<"\n";

    }
  
  output_file.close();

  cout<<"Finished file\n";


  return 0;
  
}


// vector<int> load_color_file(string & color_file)
// {

//   vector<int> color_vec;

//   ifstream color;
//   color.open(color_file.c_str());

//   if(color.is_open())
//     {
//       while(color.good())
// 	{
// 	  string line;
// 	  bool flag = true;

// 	  getline(color,line);
// 	  if(line.substr(0,1)=="#")
// 	    flag = false;

// 	  if(flag)
// 	    {
 
// 	      istringstream line_stream(line);

// 	      while(line_stream.good())
// 		{
// 		  int val;
// 		  if(line_stream>>val)
// 		    color_vec.push_back(val);
// 		}
// 	    }
// 	}
   
//       color.close();
   
//       return color_vec;
//     }
//   else
//     {
//       cerr<<"Failed to open color file\n";
//       exit(EXIT_FAILURE);
//     }

  
//   return color_vec;
// }
