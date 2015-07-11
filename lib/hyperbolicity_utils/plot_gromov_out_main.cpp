/*
Main file for plotting gromov hyperbolicity output.  Options need updated.

-i input file...REQUIRED
-o output directory...default = "../Graphs/Results/"
-n output files prefix...default = "graph"
-d if y, indicates a directed graph...default = undirected graph
-g if f, indicates only gromov four point, if s indicates only slim delta, if blank both will be calculated
-s if y, suppresses terminal output...default = false
*/

#include "../graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

int main(int argc, char *argv[])
{
  string inputFile = "../Graphs/test_graph.txt";
  string outputDirectory  = "../Graphs/Results/";
  string outputFilePrefix = "graph";
  bool four_point = true;
  bool calc_slim  = true;
  bool directed = false;
  bool suppressOutput = false;
  
  //User input

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	inputFile = argv[i+1];
      else if(strcmp(argv[i], "-n")==0)
	outputFilePrefix = argv[i+1];
      else if(strcmp(argv[i],"-o")==0)
	outputDirectory = argv[i+1];
      else if(strcmp(argv[i],"-d")==0)
	{	
	  if(strcmp(argv[i+1], "y")==0)
	    directed = true;
	}
      else if(strcmp(argv[i],"-s")==0)
	{	  
	  if(strcmp(argv[i+1], "y")==0)
	    suppressOutput = true;
	}
      else if(strcmp(argv[i], "-g")==0)
	{
	  if(strcmp(argv[i+1], "f")==0)
	    calc_slim = false;
	  
	  if(strcmp(argv[i+1], "s")==0)
	    four_point = false;
	}
    }

  ////Default values

  if(inputFile.length()==0)
    {
      cerr<<"You must provide an input file\n";
      exit(EXIT_FAILURE);
    }
  else if(!suppressOutput)
    cout<<"Input File: "<<inputFile<<"\n";

  if(outputFilePrefix.length()==0)
    {
      cerr<<"Invalid output prefix.  Default output 'graph' will be used"<<"\n";
      outputFilePrefix = "graph";
    }
  //Make directories

  int l = outputDirectory.length();
  int c = outputDirectory.compare(l-1,1,"/");
  if(c != 0)
    outputDirectory.append("/");


  // int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  // if(direct_made!=0 && errno!=EEXIST)
  //   {
  //     cerr<<"Failed to make output directory\n";
  //     exit(EXIT_FAILURE);
  //   }
  // outputDirectory.append(outputFilePrefix);
  // outputDirectory.append("/");
  // direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);  ////
  
  // if(direct_made!=0 && errno!=EEXIST)
  //   {
  //     cerr<<"Failed to make output directory folder\n";
  //     exit(EXIT_FAILURE);
  //   }

  ////////////////////////////////////
  //**********************************
  ////////////////////////////////////

  ifstream input_data;

  string line;
  double max_delta = -1;

  vector<string> dummy;
  vector< vector<string> > new_file(1,dummy);

  input_data.open(inputFile.c_str());

  if(input_data.is_open())
    {

      while(input_data.good())
	{
	  getline(input_data,line);
	  if(line.substr(0,1) != "#" && line.substr(0,1) != "\n")
	    {
	      size_t pos1 = line.find("\t");
	      if(pos1 != string::npos)
		{
		  size_t pos2 = line.find("\n");
		  string delta = line.substr(0,pos1);
		  stringstream ss;
		  ss<<delta;
		  double new_delta;
		  ss>>new_delta;
		  
		  //cout<<delta<<"\t"<<new_delta<<"\n";
		  string percent = line.substr(pos1+1,pos2-pos1);
		  if(new_delta>max_delta)
		    {
		      max_delta = new_delta;
		      new_file.push_back(dummy);
		      new_file[2*new_delta].push_back(delta);
		    }

		  if(four_point)
		    new_delta = 2*new_delta;
		  
		  new_file[new_delta].push_back(percent);
		}
	    }
	}
    } 

  input_data.close();

  ofstream temp_plot;

  string temp_file = inputFile;
  //temp_file.append(outputFilePrefix);
  //temp_file.append("_gromov.txt");
  temp_plot.open(temp_file.c_str());
  
  for(int i=0;i<new_file.size();i++)
    {
      vector<string> new_line = new_file[i];

      for(int j=0;j<new_line.size();j++)
	{
	  temp_plot<<new_line[j]<<"\t";
	}

      temp_plot<<"\n";
    }

  temp_plot.close();

  return 0;
}
