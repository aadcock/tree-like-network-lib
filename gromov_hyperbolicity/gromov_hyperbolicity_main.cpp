/*
Main file for calculating gromov hyperbolicity on a graph.  Executable options:

-i input file...REQUIRED
-o output directory...default = "../Graphs/Results/"
-n output files prefix...default = "graph"
-d if y, indicates a directed graph...default = undirected graph
-g if f, indicates only gromov four point, if s indicates only slim delta, if blank both will be calculated
-s if y, suppresses terminal output...default = false
*/

#include "graph_lib_boost.hpp"
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
  bool calc_fat   = true;
  bool quadruped  = false;
  bool directed = false;
  bool suppressOutput = false;
  bool single = false;

  int quads_to_keep = 1;
  int dist_keep     = 0;

  
  //User input

  if(argc < 2)
    {
      cout<<"-i <input_file> \n-n <output_file_prefix> \n-o <output_file_directory> \n-d <y=directed graph input> \n-s <y=suppress terminal output> \n-g gromov options: \n<four> performs four-point calculation only \n<slim> performs slim calculation only \n<fat> performs fat calculation only \n<four_quad> performs four-point calculations and saves some number of geodesic quadrupeds in dimacs format with a index file to original graph \nif this option is not specified fat/slim/four-point calculations are performed without quadrupeds > \n-quads <num_to_keep>\n-dist <lower_bound_on_quad_dist>\n-single returns only a single quadruped, if four_quad is enabled";
      return(-1);
    }

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
	  if(strcmp(argv[i+1], "four")==0)
	    {
	      calc_slim = false;
	      calc_fat  = false;
	    }

	  if(strcmp(argv[i+1], "four_quad")==0)
	    {
	      calc_slim  = false;
	      calc_fat   = false;
	      four_point = false;
	      quadruped  = true;
	    }
	  
	  if(strcmp(argv[i+1], "slim")==0)
	    {
	      four_point = false;
	      calc_fat   = false;
	    }

	  if(strcmp(argv[i+1], "fat")==0)
	    {	    
	      four_point = false;
	      calc_slim  = false;
	    }

	}
      else if(strcmp(argv[i],"-quads")==0)
	{
	  istringstream ss(argv[i+1]);
	  ss>>quads_to_keep;
	}
      else if(strcmp(argv[i],"-dist")==0)
	{
	  istringstream ss(argv[i+1]);
	  ss>>dist_keep;
	}
      else if(strcmp(argv[i],"-single")==0)
	{
	  single = true;
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


  int direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);
  if(direct_made!=0 && errno!=EEXIST)
    {
      cerr<<"Failed to make output directory\n";
      exit(EXIT_FAILURE);
    }
  outputDirectory.append(outputFilePrefix);
  outputDirectory.append("/");
  direct_made = mkdir(outputDirectory.c_str(),S_IRWXU);  ////
  
  if(direct_made!=0 && errno!=EEXIST)
    {
      cerr<<"Failed to make output directory folder\n";
      exit(EXIT_FAILURE);
    }

  ////////////////////////////////////
  //**********************************
  ////////////////////////////////////

  string delimiter = "\t";

  graph_traits<Graph>::vertex_iterator vi,vie;
  
  double t1, t2, t3, t4, t5, t6, t7;
  
  t1 = time(NULL);
  Graph G = loadGraph(inputFile,delimiter,directed);  
  t2 = time(NULL);
  cout<<"Load time: "<<t2-t1<<"\n";

  G = connected(G,-1);
  t3 = time(NULL);
  cout<<"Connected time: "<<t3-t2<<"\n";

  string outputFile = outputDirectory;
  outputFile.append(outputFilePrefix);
  BFS_distance(G);
  t4 = time(NULL);
  cout<<"BFS time: "<<t4-t3<<"\n";

  unsigned int dmax = 0;
  for(tie(vi,vie) = vertices(G); vi != vie; vi++)
    {
      Vert v = *vi;
      size_t size = G[v].distances.size();
      
      for(size_t i=0;i<size;i++)
	{
	  if(G[v].distances[i] > dmax)
	    dmax = G[v].distances[i];
	}
    }
  if(four_point)
    {
      vector< vector<double> > delta = calc_gromov(G,dmax);

      unsigned long numQuad = 0;
      for(size_t i=0;i<delta.size();i++)
	numQuad += delta[i][0];

      ofstream outGrom;
      string outputFileGrom = outputFile;
      outputFileGrom.append("_gromov.txt");
      outGrom.open(outputFileGrom.c_str());

      ofstream outAvg;
      string outputFileAvg = outputFile;
      outputFileAvg.append("_avg_distr.txt");
      outAvg.open(outputFileAvg.c_str());
      outGrom<<"#Column format: \n";
      outGrom<<"#Max Length in Quadruplet\tNumber of Quadruplets\n";

      for(size_t i=0;i<delta.size();i++)
	{
	  outGrom<<"#~\t"<<i+1<<"\t"<<delta[i][0]<<"\n";
	}

      outGrom<<"#Column format: \n";
      outGrom<<"#Delta\tPercent of quadruplets for each length\n";

      for(size_t j=1;j<delta[0].size();j++)
	{
	  //size_t size = delta[i].size();	

	  
	  outGrom<<(j-1)/2.0<<"\t";
	  for(size_t i=0;i<delta.size();i++)
	    {
	      if(delta[i][0]!=0)
		(delta[i])[j] = (delta[i])[j]/(delta[i])[0];
	  
	      
	      outGrom<<delta[i][j]<<"\t";
	    }
	  outGrom<<"\n";
	  
	}

      for(size_t i=0;i<delta.size();i++)
	{
	  double avg = 0;
	  for(size_t j=1;j<delta[i].size();j++)
	    {
	      avg += delta[i][j]*(j-1)/2.0;
	    }
	  outAvg<<i+1<<"\t"<<avg<<"\n";
	}

      outGrom.close();
      outAvg.close();

      t5 = time(NULL);
      if(!suppressOutput)
	{
	  cout<<"Gromov time: "<<t5-t4<<"\n";

	  cout<<"Results for "<<inputFile<<":\n";
	  cout<<"Diameter of Graph: "<<dmax<<"\n";
	  cout<<"# of quadruplets: "<<numQuad<<"\n";
	  cout<<"Delta distribution of Graph: "<<"\n";
	  cout.setf(ios::fixed,ios::floatfield);


	  for(size_t i=0;i<delta.size();i++)
	    {
	      cout<<"Max. Length in Quadruple: "<<i+1<<"\nPercent of quadruplets: "<<delta[i][0]/numQuad<<"\n";
	      double avg = 0;
	      for(size_t j=1;j<delta[i].size();j++)
		{	
		  cout.precision(1);
		  cout<<(j-1)/2.0;
		  cout.precision(10);
		  cout<<": "<<delta[i][j]<<"\n";
	  
		  avg += delta[i][j]*(j-1)/2.0;
		}
	      cout<<"Average: "<<avg<<"\n";
      

	    }

	  cout<<"Total time: "<<t5-t1<<"\n";
	}
    }
  // ////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////

  if(quadruped)
    {
      string quad_output = outputDirectory;
      quad_output.append(outputFilePrefix);

      vector< vector<double> > delta = calc_gromov_quads(G,dmax,quads_to_keep,dist_keep,quad_output,single);

      unsigned long numQuad = 0;
      for(size_t i=0;i<delta.size();i++)
	numQuad += delta[i][0];

      ofstream outGrom;
      string outputFileGrom = outputFile;
      outputFileGrom.append("_gromov.txt");
      outGrom.open(outputFileGrom.c_str());

      ofstream outAvg;
      string outputFileAvg = outputFile;
      outputFileAvg.append("_avg_distr.txt");
      outAvg.open(outputFileAvg.c_str());
      outGrom<<"#Column format: \n";
      outGrom<<"#Max Length in Quadruplet\tNumber of Quadruplets\n";

      for(size_t i=0;i<delta.size();i++)
	{
	  outGrom<<"#~\t"<<i+1<<"\t"<<delta[i][0]<<"\n";
	}

      outGrom<<"#Column format: \n";
      outGrom<<"#Delta\tPercent of quadruplets for each length\n";

      for(size_t j=1;j<delta[0].size();j++)
	{
	  //size_t size = delta[i].size();	

	  
	  outGrom<<(j-1)/2.0<<"\t";
	  for(size_t i=0;i<delta.size();i++)
	    {
	      if(delta[i][0]!=0)
		(delta[i])[j] = (delta[i])[j]/(delta[i])[0];
	  
	      
	      outGrom<<delta[i][j]<<"\t";
	    }
	  outGrom<<"\n";
	  
	}

      for(size_t i=0;i<delta.size();i++)
	{
	  double avg = 0;
	  for(size_t j=1;j<delta[i].size();j++)
	    {
	      avg += delta[i][j]*(j-1)/2.0;
	    }
	  outAvg<<i+1<<"\t"<<avg<<"\n";
	}

      outGrom.close();
      outAvg.close();

      t5 = time(NULL);
      if(!suppressOutput)
	{
	  cout<<"Gromov time: "<<t5-t4<<"\n";

	  cout<<"Results for "<<inputFile<<":\n";
	  cout<<"Diameter of Graph: "<<dmax<<"\n";
	  cout<<"# of quadruplets: "<<numQuad<<"\n";
	  cout<<"Delta distribution of Graph: "<<"\n";
	  cout.setf(ios::fixed,ios::floatfield);


	  for(size_t i=0;i<delta.size();i++)
	    {
	      cout<<"Max. Length in Quadruple: "<<i+1<<"\nPercent of quadruplets: "<<delta[i][0]/numQuad<<"\n";
	      double avg = 0;
	      for(size_t j=1;j<delta[i].size();j++)
		{	
		  cout.precision(1);
		  cout<<(j-1)/2.0;
		  cout.precision(10);
		  cout<<": "<<delta[i][j]<<"\n";
	  
		  avg += delta[i][j]*(j-1)/2.0;
		}
	      cout<<"Average: "<<avg<<"\n";
      

	    }

	  cout<<"Total time: "<<t5-t1<<"\n";
	}
    }


  // ///////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  if(calc_slim)
    {
      t5 = time(NULL);
      vector< vector<double> > delta_slim = calc_delta_slim(G,dmax);
      t6 = time(NULL);
      unsigned long numTri = 0;
      for(size_t i=0;i<delta_slim.size();i++)
	numTri += delta_slim[i][0];

      ofstream outSlim;
      string outputFileSlim = outputFile;
      outputFileSlim.append("_slim.txt");
      outSlim.open(outputFileSlim.c_str());

      ofstream outSlimAvg;
      string outputFileSlimAvg = outputFile;
      outputFileSlimAvg.append("_slim_avg_distr.txt");
      outSlimAvg.open(outputFileSlimAvg.c_str());
      outSlim<<"#Column format: \n";
      outSlim<<"#Max side\tNumber of Triplets\n\n";


      for(size_t i=0;i<delta_slim.size();i++)
	{
	  outSlim<<"#~\t"<<i+1<<"\t"<<delta_slim[i][0]<<"\n";
	}

      outSlim<<"#Column format: \n";
      outSlim<<"#Delta\tPercent of Triplets\n\n";      

      for(size_t j=1;j<delta_slim[0].size();j++)
	{
	  //size_t size = delta_slim[i].size();	

	  
	  outSlim<<(j-1)<<"\t";
	  for(size_t i=0;i<delta_slim.size();i++)
	    {
	      if(delta_slim[i][0]!=0)
		(delta_slim[i])[j] = (delta_slim[i])[j]/(delta_slim[i])[0];
	  
	      
	      outSlim<<delta_slim[i][j]<<"\t";
	    }
	  outSlim<<"\n";
	  
	}

      for(size_t i=0;i<delta_slim.size();i++)
	{
	  double avg = 0;
	  for(size_t j=1;j<delta_slim[i].size();j++)
	    {
	      avg += delta_slim[i][j]*(j-1);
	    }
	  outSlimAvg<<i+1<<"\t"<<avg<<"\n";
	}

      outSlim.close();
      outSlimAvg.close();

      if(!suppressOutput)
	{
	  cout<<"Slim time: "<<t6-t5<<"\n";

	  cout<<"Results for "<<inputFile<<":\n";
	  cout<<"Diameter of Graph: "<<dmax<<"\n";
	  cout<<"# of Triplets: "<<numTri<<"\n";
	  cout<<"Delta slim distribution of Graph: "<<"\n";
	  cout.setf(ios::fixed,ios::floatfield);


	  for(size_t i=0;i<delta_slim.size();i++)
	    {
	      cout<<"Max. Length in Triplet: "<<i+1<<"\nPercent of triplets: "<<delta_slim[i][0]/numTri<<"\n";
	      double avg = 0;
	      for(size_t j=1;j<delta_slim[i].size();j++)
		{	
		  cout.precision(1);
		  cout<<(j-1);
		  cout.precision(10);
		  cout<<": "<<delta_slim[i][j]<<"\n";
	  
		  avg += delta_slim[i][j]*(j-1);
		}
	      cout<<"Average: "<<avg<<"\n";
      

	    }
	  t7 = time(NULL);
	  cout<<"Total time: "<<t7-t5<<"\n";
	}
    }

// ////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////

  if(calc_fat)
    {
      t5 = time(NULL);
      vector< vector<double> > delta_fat = calc_delta_fat(G,dmax);
      t6 = time(NULL);
      unsigned long numTri = 0;
      for(size_t i=0;i<delta_fat.size();i++)
	numTri += delta_fat[i][0];

      ofstream outFat;
      string outputFileFat = outputFile;
      outputFileFat.append("_fat.txt");
      outFat.open(outputFileFat.c_str());

      ofstream outFatAvg;
      string outputFileFatAvg = outputFile;
      outputFileFatAvg.append("_fat_avg_distr.txt");
      outFatAvg.open(outputFileFatAvg.c_str());
      outFat<<"#Column format: \n";
      outFat<<"#Max side\tNumber of Triplets\n\n";


      for(size_t i=0;i<delta_fat.size();i++)
	{
	  outFat<<"#~\t"<<i+1<<"\t"<<delta_fat[i][0]<<"\n";
	}

      outFat<<"#Column format: \n";
      outFat<<"#Delta\tPercent of Triplets\n\n";      

      for(size_t j=1;j<delta_fat[0].size();j++)
	{
	  //size_t size = delta_fat[i].size();	

	  
	  outFat<<(j-1)<<"\t";
	  for(size_t i=0;i<delta_fat.size();i++)
	    {
	      if(delta_fat[i][0]!=0)
		(delta_fat[i])[j] = (delta_fat[i])[j]/(delta_fat[i])[0];
	  
	      
	      outFat<<delta_fat[i][j]<<"\t";
	    }
	  outFat<<"\n";
	  
	}

      for(size_t i=0;i<delta_fat.size();i++)
	{
	  double avg = 0;
	  for(size_t j=1;j<delta_fat[i].size();j++)
	    {
	      avg += delta_fat[i][j]*(j-1);
	    }
	  outFatAvg<<i+1<<"\t"<<avg<<"\n";
	}

      outFat.close();
      outFatAvg.close();

      if(!suppressOutput)
	{
	  cout<<"Fat time: "<<t6-t5<<"\n";

	  cout<<"Results for "<<inputFile<<":\n";
	  cout<<"Diameter of Graph: "<<dmax<<"\n";
	  cout<<"# of Triplets: "<<numTri<<"\n";
	  cout<<"Delta fat distribution of Graph: "<<"\n";
	  cout.setf(ios::fixed,ios::floatfield);


	  for(size_t i=0;i<delta_fat.size();i++)
	    {
	      cout<<"Max. Length in Triplet: "<<i+1<<"\nPercent of triplets: "<<delta_fat[i][0]/numTri<<"\n";
	      double avg = 0;
	      for(size_t j=1;j<delta_fat[i].size();j++)
		{	
		  cout.precision(1);
		  cout<<(j-1);
		  cout.precision(10);
		  cout<<": "<<delta_fat[i][j]<<"\n";
	  
		  avg += delta_fat[i][j]*(j-1);
		}
	      cout<<"Average: "<<avg<<"\n";
      

	    }
	  t7 = time(NULL);
	  cout<<"Total time: "<<t7-t5<<"\n";
	}
    }

  return 0;
}
