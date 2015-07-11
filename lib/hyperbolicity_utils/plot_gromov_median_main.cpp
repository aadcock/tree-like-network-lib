/*
Main file for plotting median or mean of gromov hyperbolicity. output options need updated.

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
  
  
  //********************************


  ifstream input_data;

  string line;
  double max_delta = -1;

  //first index gives the delta value, the second index gives the quadruplet size
  vector< vector<double> > gromov_output;

  input_data.open(inputFile.c_str());

  if(input_data.is_open())
    {

      while(input_data.good())
	{
	  getline(input_data,line);
	  if(line.substr(0,1) != "#")
	    {
	      vector<double> temp;
	      stringstream ss(line);
	      double num;
	      while(ss >> num)
		temp.push_back(num);
	      gromov_output.push_back(temp);
	    }
	}
    } 

  input_data.close();
  size_t n = gromov_output[1].size();
  vector<double> medians(n-1,0);
  vector<double> sums(n-1,0);
  vector<double> maxs(n-1,0);
  vector<double> avg(n-1,0);
  
  for(size_t i=0;i<gromov_output.size();i++)
    {
      for(size_t j=1;j<gromov_output[i].size();j++)
	{
	  double gromovOut = gromov_output[i][j];
	  //cout<<i<<"\n"<<j<<"\n"<<gromov_output[i].size()<<"\n"<<n<<"\n Done with that set.\n";
	  bool sum_flag = false;
	  if(sums[j-1] < .5)
	    sum_flag = true;

	  sums[j-1] = sums[j-1] + gromovOut;
	  double delta = i;

	  if(four_point)
	    delta = double(i)/2.0;

	  if(sums[j-1]>.5 && sum_flag)
	    medians[j-1] = delta;

	  if(gromovOut!=0)
	    {
	      if(maxs[j-1] < delta)
		maxs[j-1] = delta;
	    }

	  avg[j-1] += gromovOut*delta;
	}
    }

  ofstream output_file;

  string outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);
  outputFile.append("_median.txt");
  
  output_file.open(outputFile.c_str());

  for(size_t i=0;i<medians.size();i++)
    {
      output_file<<i+1<<"\t"<<medians[i]<<"\n";
    }

  output_file.close();
  outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);
  outputFile.append("_max.txt");

  output_file.open(outputFile.c_str());
  
  for(size_t i=0;i<maxs.size();i++)
    {
      output_file<<i+1<<"\t"<<maxs[i]<<"\n";
    }
  output_file.close();

  outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);
  outputFile.append("_max_scale.txt");

  output_file.open(outputFile.c_str());
  
  for(size_t i=0;i<maxs.size();i++)
    {
      double scaled_gro = double(maxs[i])/double(i+1);
      output_file<<i+1<<"\t"<<scaled_gro<<"\n";
    }
  output_file.close();


  outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);
  outputFile.append("_avg_scale.txt");

  output_file.open(outputFile.c_str());
  
  for(size_t i=0;i<avg.size();i++)
    {
      double scaled_gro = double(avg[i])/double(i+1);
      output_file<<i+1<<"\t"<<scaled_gro<<"\n";
    }
  output_file.close();
  
  return 0;
}
