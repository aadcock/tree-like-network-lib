/*
Produces plot of gromov hyperbolicity stats.  output options need updated.

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

  string outfilename = outputDirectory;
  outfilename.append(outputFilePrefix);
  

  ifstream temp;
  temp.open(inputFile.c_str());
  string line;
  int count = 0;
  bool flag = true;

  vector<long> quad_size;
  if(temp.is_open())
    {
      stringstream ss;
      while(flag)
	{
	  getline(temp,line);
	  string first_two_letters="";

	  if(line.size()>=2)
	    first_two_letters=line.substr(0,2);
	  
	  if(line.substr(0,1) != "#")
	    {
	      line.erase(line.find_last_not_of(" \n\r\t")+1);
	      for(int i=0;i<line.length();i++)
		{
		  char ch = line[i];
		  string s;
		  s = ch;
		  if(s.compare("\t")==0)
		    {
		      //cout<<line<<"\n";
		      //cout<<s<<"\n";
		      //cout<<count<<"\n";
		      count++;
		    }
		}
	      flag = false;
	    }
	  else if(strcmp("#~",first_two_letters.c_str())==0)
	    {
	      int pos = line.find_last_of("\t");
	      stringstream snum;
	      //	      cout<<pos<<" "<<line<<"  "<<first_two_letters<<" here?\n";
	      snum<<line.substr(pos,line.size()-pos);
	      //	      cout<<"nope\n";
	      long tmp;
	      snum>>tmp;
	      cout<<tmp<<"\n";
	      quad_size.push_back(tmp);
	    }
	}

    }

  temp.close();
  //cout<<"count: "<<count<<"\n";

  ////////////////////////////////////////////
  //Create delta distribution file
  ////////////////////////////////////////////

  ofstream out_file_stream;
  string out_file_quad = outputDirectory;

  out_file_quad.append(outputFilePrefix);
  out_file_quad.append("_quad_distr.txt");
  
  out_file_stream.open(out_file_quad.c_str());
  temp.open(inputFile.c_str());
  vector<double> sum_frac_delta;
  vector<double> frac_size(quad_size.size(),0);
  int index_sfd = -1;
  long sum = 0;

  for(int i=0;i<quad_size.size();++i)
    sum+=quad_size[i];
  
  for(int i=0;i<frac_size.size();++i)
    frac_size[i] = double(quad_size[i])/double(sum);
  
  if(temp.is_open())
    {
      stringstream ss;
      while(temp.good())
	{
	  
	  getline(temp,line);
	  if(line.substr(0,1) != "#")
	    {
	      sum_frac_delta.push_back(0);
	      bool flag_delt = false;
	      istringstream buffer(line);
	      for(int i=0;i<=count;i++)
		{
		  double z;

		  if(buffer>>z)
		    {		  
		      if(i!=0)
			sum_frac_delta[index_sfd] += z*frac_size[i-1];
		      else
			{
			  flag_delt = true;
			  index_sfd=2*z;
			}
		    }
		  
		}
	      
	      if(flag_delt)
		out_file_stream<<double(index_sfd)/2.0<<"\t"<<sum_frac_delta[index_sfd]<<"\n";
	    }

	}
    }

  temp.close();
  out_file_stream.close();
  
  ////////////////////////////////////////////
   vector<int> columns_to_plot;
   vector<string> label;
  // for(int i=2;i<=count+1;i++)
  //   {
  //     columns_to_plot.push_back(i);
  //     stringstream ss;
  //     ss<<(i-1);
  //     string s = "Quads with max Side: ";
  //     s.append(ss.str());

  //     label.push_back(s);
  //   }

  columns_to_plot.push_back(2);
  label.push_back("Delta Distribution");

  string title = inputFile;
  title.append(" Delta Four Point Distributions");

  string xlabel = "Delta";
  string ylabel = "Percent of Quadruplets";
  string xtics = "";

  produce_plot(out_file_quad, out_file_quad,columns_to_plot, label, title, xlabel, ylabel, xtics, false);
  
  return 0;
}
