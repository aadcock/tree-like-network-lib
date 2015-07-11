/*
  A short utility for getting basic graph stats quickly
  Main file for calculating connected component size, edges, etc.  Output options need updated.

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
#include <cmath>

int main(int argc, char *argv[])
{
  string inputFile = "./dissect/as/as20000102.dimacs";  
  string outputFile = "./graph";
  bool display_output = true;
  bool write_graph_files = true;
  bool gro_calc = false;
  bool ncp_calc = false;
  int core_to_start = -1;
  //User input

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if(strcmp(argv[i], "-i")==0)
	inputFile = argv[i+1];

      if(strcmp(argv[i], "-o")==0)
	outputFile = argv[i+1];

      if(strcmp(argv[i], "-silent")==0)
	display_output = false;

      if(strcmp(argv[i], "-nofile")==0)
	write_graph_files = false;

      if(strcmp(argv[i], "-gromov")==0)
	gro_calc = true;
   
      if(strcmp(argv[i], "-ncp")==0)
	ncp_calc = true;

      if(strcmp(argv[i], "-core")==0)
	core_to_start = atoi(argv[i+1]);
    }

  ////Default values

  if(inputFile.length()==0)
    {
      cerr<<"You must provide an input file\n";
      exit(EXIT_FAILURE);
    }

  //////////////////////////////////////////////////////////////
  //-----------------------------------------------------------
  //////////////////////////////////////////////////////////////

  Graph G_o = loadGraph(inputFile,"\t",false);
  Graph G = connected(G_o);
 
  vector<int> core;
  vector<v_size_t> distances;
  vector<v_size_t> num_paths;

  v_size_t N_paths,D,n;
  int max_core = 0;
  double A_d, A_p;

  A_d = 0;
  A_p = 0;

  D       = 0;
  N_paths = 0;
  n       = num_vertices(G);

  graph_traits<Graph>::vertex_iterator vi, vie;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  core = k_core(G);
  

  for(int i=0;i<core.size();++i)
    {
      if(core[i] > max_core)
	max_core = core[i];
    }
  
  if(core_to_start<=0)
    core_to_start = max_core;

  for(tie(vi,vie) = vertices(G); vi != vie; ++vi)
    {
      Vert v = *vi;
      //cout<<G[v].core_num<<" ";
      BFS_source_all(v,G,distances,num_paths);

      size_t size = distances.size();

      for(size_t i=0;i<size;++i)
	{

	  if(distances[i] > D)
	    D = distances[i];

	  A_p += num_paths[i] * distances[i];
	  A_d += distances[i];
	  N_paths += num_paths[i];
	}
    }

  A_p /= N_paths;
  A_d /= (pow(n, 2) - n);

  cout<<"Max Core: "<<max_core<<" The diameter: "<<D<<" The avg. path: "<<A_p<<" The avg. dist: "<<A_d<<" The size: "<<n<<"\n";
  cout.flush();

  string periph_diam   = outputFile;
  string collapse_diam = outputFile;
  string core_diam     = outputFile;

  ofstream file_periph;
  ofstream file_collapse;
  ofstream file_core;

  periph_diam.append("_periph_diam.out");
  collapse_diam.append("_collapse_diam.out");
  core_diam.append("_core_diam.out");

  file_periph.open(periph_diam.c_str());
  file_collapse.open(collapse_diam.c_str());
  file_core.open(core_diam.c_str());
  
  file_periph<<(max_core+1)<<"\t"<<D<<"\t"<<A_d<<"\t"<<A_p<<"\n";
  file_collapse<<(max_core+1)<<"\t"<<D<<"\t"<<A_d<<"\t"<<A_p<<"\n";
  file_core<<(max_core+1)<<"\t"<<0<<"\t"<<0<<"\t"<<"\n";

  Graph G_core, G_periph, G_collapse;
  vector<v_size_t> core_list;
  vector<v_size_t> periph_list;

  for(int i=1;i<=max_core;i++)
    periph_list.push_back(i);


  for(int i=max_core; i>0; --i)
    {

      core_list.push_back(i);
      periph_list.pop_back();

      if(i<=core_to_start)
	{
	  G_core     = get_k_shells(core_list,G);
	  G_periph   = get_k_shells(periph_list,G);
	  G_collapse = collapse_core(i,G); 

	  stringstream ssperiph;
	  ssperiph<<i;
	  string outputPeriph   = outputFile;
	  string outputCore     = outputFile;
	  string outputCollapse = outputFile;
	  outputPeriph.append("_periph_");
	  outputCore.append("_core_");
	  outputCollapse.append("_collapse_");
	  outputPeriph.append(ssperiph.str());
	  outputCore.append(ssperiph.str());
	  outputCollapse.append(ssperiph.str());
	  outputPeriph.append(".dimacs");
	  outputCore.append(".dimacs");
	  outputCollapse.append(".dimacs");

	  if(write_graph_files)
	    {
	      write_dimacs(G_periph,outputPeriph);
	      write_dimacs(G_core,outputCore);
	      write_dimacs(G_collapse,outputCollapse);
	    }
	  //Need to make something that checks number of connected components, and then runs this on each one.  Just because the whole graph is connected doesn't mean the cores have to be.

	  D = 0;
	  v_size_t n_core = num_vertices(G_core);
	  A_p = 0;
	  A_d = 0;
	  N_paths = 0;
	  for(tie(vi,vie) = vertices(G_core); vi != vie; ++vi)
	    {
	      Vert v = *vi;

	      BFS_source_all(v,G_core,distances,num_paths);

	      size_t size = distances.size();

	      for(size_t j=0;j<size;++j)
		{

		  if(distances[j] > D)
		    D = distances[j];

		  A_p += num_paths[j]*distances[j];
		  A_d += distances[j];
		  N_paths += num_paths[j];
		}
	    }

	  A_p /= N_paths;

	  vector<bool> found(n_core,false);
	  v_size_t current_comp = 0;
	  vector<v_size_t> comp_size;
	  vector<v_size_t> component(n_core);

	  index = get(vertex_index,G_core);

	  for(tie(vi,vie)=vertices(G_core);vi!=vie;++vi)
	    {
	      Vert v = *vi;
	      v_size_t vind = index[v];
	      v_size_t current_comp_size = 0;

	      if(!found[vind])
		{
		  BFS_vertices_found(v,G_core,found,current_comp,component,current_comp_size);
		  comp_size.push_back(current_comp_size);
		  ++current_comp;
		}
	    }
	  if(display_output)
	    cout<<"\nCore: "<<i<<"\n";

	  double sum_of_components = 0;
      
	  if(display_output)
	    cout<<"# of components in core: "<<comp_size.size()<<" Size of core: "<<num_vertices(G_core)<<"\n";

      
	  for(size_t j=0;j<comp_size.size();++j)
	    {
	      if(comp_size[j] > (num_vertices(G_core)/10.0) && display_output)
		cout<<"Component "<<j<<" Size: "<<comp_size[j]<<"\n";

	      sum_of_components += pow(comp_size[j],2)-comp_size[j];
	    }
      

	  A_d /= sum_of_components;

	  if(display_output)
	    cout<<"The diameter: "<<D<<" The avg. path: "<<A_p<<" The avg. dist: "<<A_d<<" The size: "<<n_core<<"\n";

	  file_core<<(i)<<"\t"<<D<<"\t"<<A_d<<"\t"<<A_p<<"\n";
	  
	  cout.flush();

	  //////////////////////////////////////////////	
	  if(gro_calc && n_core>4)
	    {
	      G_core = connected(G_core);
	      
	      //cout<<"Vertices in new core: "<<num_vertices(G_core)<<"\n";

	      BFS_distance(G_core);
	    
	      vector< vector<double> > delta = calc_gromov(G_core,D);
	      stringstream ss,tt;

	      ss<<i;
	      unsigned long numQuad = 0;
	      for(size_t j=0;j<delta.size();++j)
		numQuad += delta[j][0];

	      ofstream outGrom;
	      string outputFileGrom = outputFile;
	      outputFileGrom.append("_gromov_core_");
	      outputFileGrom.append(ss.str());
	      outputFileGrom.append(".txt");
	      outGrom.open(outputFileGrom.c_str());

	      tt<<i;

	      ofstream outAvg;
	      string outputFileAvg = outputFile;
	      outputFileAvg.append("_avg_distr_core_");
	      outputFileAvg.append(tt.str());
	      outputFileAvg.append(".txt");
	      outAvg.open(outputFileAvg.c_str());
	      outGrom<<"#Note: Hyperbolicity only calculated on largest component of core\n";
	      outGrom<<"#Column format: \n";
	      outGrom<<"#Max Length in Quadruplet\tNumber of Quadruplets\n";

	      for(size_t j=0;j<delta.size();++j)
		{
		  outGrom<<"#~\t"<<j+1<<"\t"<<delta[j][0]<<"\n";
		}

	      outGrom<<"#Column format: \n";
	      outGrom<<"#Delta\tPercent of quadruplets for each length\n";

	      for(size_t j=1;j<delta[0].size();j++)
		{
		  //size_t size = delta[i].size();	

	  
		  outGrom<<(j-1)/2.0<<"\t";
		  for(size_t k=0;k<delta.size();k++)
		    {
		      if(delta[k][0]!=0)
			(delta[k])[j] = (delta[k])[j]/(delta[k])[0];
	  
	      
		      outGrom<<delta[k][j]<<"\t";
		    }
		  outGrom<<"\n";
	  
		}

	      for(size_t k=0;k<delta.size();k++)
		{
		  double avg = 0;
		  for(size_t j=1;j<delta[k].size();j++)
		    {
		      avg += delta[k][j]*(j-1)/2.0;
		    }
		  outAvg<<k+1<<"\t"<<avg<<"\n";
		}

	      outGrom.close();
	      outAvg.close();
	    }

	  ////////////////////////////////////////////////////////



	  v_size_t n_periph = num_vertices(G_periph);

	  found.assign(n_periph,false);
	  current_comp = 0;
	  comp_size.clear();
	  component.assign(n_periph,0);

	  index = get(vertex_index,G_periph);

	  for(tie(vi,vie)=vertices(G_periph);vi!=vie;++vi)
	    {
	      Vert v = *vi;
	      v_size_t vind = index[v];
	      v_size_t current_comp_size = 0;

	      if(!found[vind])
		{
		  BFS_vertices_found(v,G_periph,found,current_comp,component,current_comp_size);
		  comp_size.push_back(current_comp_size);
		  ++current_comp;
		}
	    }


	  if(display_output)
	    cout<<"# of components in periphery: "<<comp_size.size()<<" Size of periphery: "<<num_vertices(G_periph)<<"\n";

	  for(size_t j=0;j<comp_size.size();++j)
	    if(comp_size[j] > (num_vertices(G_periph)/10.0) && display_output)
	      cout<<"Component "<<j<<" Size: "<<comp_size[j]<<"\n";



	  G_periph = connected(G_periph);
	  n_periph = num_vertices(G_periph);

	  /////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////
	  //ncp
	  //////////////////////////////////////////////////////////////////////

	  // v_size_t max_community = n_periph/2;
	  // vector<double> ncp_periph = ncp_calc(G_periph,max_community,1,.001,false);

	  // stringstream ssncp;
	  // ssncp<<i;
	  // string outputNcp = outputFile;
	  // ofstream output_file_ncp;
	  // outputNcp.append("_ncp_periph_");
	  // outputNcp.append(ssncp.str());
	  // outputNcp.append(".txt");

	  // output_file_ncp.open(outputNcp.c_str());

	  // for(v_size_t k=0; k<max_community;++k)
	  // 	output_file_ncp<<(k+1)<<"\t"<<ncp_periph[k]<<"\n";

	  // output_file_ncp.close();


	  // string label = outputFile;
	  // label.append(" ncp plot");

	  // string title = outputFile;
	  // title.append(" NCP Plot");

	  // vector<int> int_to_vec;
	  // int_to_vec.push_back(2);

	  // vector<string> string_to_vec;
	  // string_to_vec.push_back(label);

	  // string outputPNG = "";
	  // outputPNG.append(outputFile);
	  // outputPNG.append("_ncp_periph_");
	  // outputPNG.append(ssncp.str());
	  // outputPNG.append("_ncp_plot");

	  // produce_loglog_plot(outputNcp,outputPNG,outputPNG,int_to_vec,string_to_vec,title,"Community size", "Conductance","",false);


	  // /////////////////////////////////////////////////////////////////////////
	  // /////////////////////////////////////////////////////////////////////////
	  D = 0;
	  A_p = 0;
	  A_d = 0;
	  N_paths = 0;
	  for(tie(vi,vie) = vertices(G_periph); vi != vie; ++vi)
	    {
	      Vert v = *vi;

	      BFS_source_all(v,G_periph,distances,num_paths);

	      size_t size = distances.size();

	      for(size_t j=0;j<size;++j)
		{

		  if(distances[j] > D)
		    D = distances[j];

		  A_p += num_paths[j]*distances[j];
		  A_d += distances[j];
		  N_paths += num_paths[j];
		}
	    }

	  A_p /= N_paths;
	  A_d /= (pow(n_periph,2)-n_periph);
	


	  if(gro_calc && n_periph>4)
	    {
	      BFS_distance(G_periph);
	    
	      vector< vector<double> > delta = calc_gromov(G_periph,D);
	      stringstream ss,tt;

	      ss<<i;
	      unsigned long numQuad = 0;
	      for(size_t j=0;j<delta.size();++j)
		{
		  numQuad += delta[j][0];
		}

	      

	      ofstream outGrom;
	      string outputFileGrom = outputFile;
	      outputFileGrom.append("_gromov_periph_");
	      outputFileGrom.append(ss.str());
	      outputFileGrom.append(".txt");
	      outGrom.open(outputFileGrom.c_str());
	      
	      tt<<i;

	      ofstream outAvg;
	      string outputFileAvg = outputFile;
	      outputFileAvg.append("_avg_distr_periph_");
	      outputFileAvg.append(tt.str());
	      outputFileAvg.append(".txt");
	      outAvg.open(outputFileAvg.c_str());
	      outGrom<<"#Column format: \n";
	      outGrom<<"#Max Length in Quadruplet\tNumber of Quadruplets\n";

	      
	      for(size_t j=0;j<delta.size();++j)
		{
		  outGrom<<"#~\t"<<j+1<<"\t"<<delta[j][0]<<"\n";
		}

	      outGrom<<"#Column format: \n";
	      outGrom<<"#Delta\tPercent of quadruplets for each length\n";
	     
	      for(size_t j=1;j<delta[0].size();j++)
		{
		  //size_t size = delta[i].size();	

	  
		  outGrom<<(j-1)/2.0<<"\t";
		  for(size_t k=0;k<delta.size();k++)
		    {
		      if(delta[k][0]!=0)
			(delta[k])[j] = (delta[k])[j]/(delta[k])[0];
	  
	      
		      outGrom<<delta[k][j]<<"\t";
		    }
		  outGrom<<"\n";	  
		}

	      for(size_t k=0;k<delta.size();k++)
		{
		  double avg = 0;
		  for(size_t j=1;j<delta[k].size();j++)
		    {
		      avg += delta[k][j]*(j-1)/2.0;
		    }
		  outAvg<<k+1<<"\t"<<avg<<"\n";
		}

	      outGrom.close();
	      outAvg.close();
	    }

	  //cout<<"\nCore: "<<i<<"\n";
	  if(display_output)      
	    cout<<"The diameter: "<<D<<" The avg. path: "<<A_p<<" The avg. dist: "<<A_d<<" The size: "<<n_periph<<"\n";

	  file_periph<<(i)<<"\t"<<D<<"\t"<<A_d<<"\t"<<A_p<<"\n";

	  cout.flush();
	  //////////////////////////////////////////////////
	  //////////////////////////////////////////////////

	  D = 0;
	  v_size_t n_collapse = num_vertices(G_collapse);
	  A_p = 0;
	  A_d = 0;
	  N_paths = 0;
	  for(tie(vi,vie) = vertices(G_collapse); vi != vie; ++vi)
	    {
	      Vert v = *vi;

	      BFS_source_all(v,G_collapse,distances,num_paths);

	      size_t size = distances.size();

	      for(size_t j=0;j<size;++j)
		{

		  if(distances[j] > D)
		    D = distances[j];

		  A_p += num_paths[j]*distances[j];
		  A_d += distances[j];
		  N_paths += num_paths[j];
		}
	    }

	  A_p /= N_paths;

	  found.assign(n_collapse,false);
	  current_comp = 0;
	  comp_size.clear();
	  component.assign(n_collapse,0);

	  index = get(vertex_index,G_collapse);

	  for(tie(vi,vie)=vertices(G_collapse);vi!=vie;++vi)
	    {
	      Vert v = *vi;
	      v_size_t vind = index[v];
	      v_size_t current_comp_size = 0;

	      if(!found[vind])
		{
		  BFS_vertices_found(v,G_collapse,found,current_comp,component,current_comp_size);
		  comp_size.push_back(current_comp_size);
		  ++current_comp;
		}
	    }


	  sum_of_components = 0;
      
	  if(display_output)
	    cout<<"# of components in collapse: "<<comp_size.size()<<" Size of collapse: "<<num_vertices(G_collapse)<<"\n";

	  cout.flush();
	  for(size_t j=0;j<comp_size.size();++j)
	    {
	      if(comp_size[j] > (num_vertices(G_collapse)/10.0) && display_output)
		cout<<"Component "<<j<<" Size: "<<comp_size[j]<<"\n";

	      sum_of_components += pow(comp_size[j],2)-comp_size[j];
	    }
      

	  A_d /= sum_of_components;

	  if(display_output)
	    cout<<"The diameter: "<<D<<" The avg. path: "<<A_p<<" The avg. dist: "<<A_d<<" The size: "<<n_collapse<<"\n";

	  file_collapse<<(i)<<"\t"<<D<<"\t"<<A_d<<"\t"<<A_p<<"\n";

	  cout.flush();
	  //////////////////////////////////////////////	
	  if(gro_calc && n_collapse>4)
	    {
	      BFS_distance(G_collapse);
	    
	      vector< vector<double> > delta = calc_gromov(G_collapse,D);
	      stringstream ss,tt;

	      ss<<i;
	      unsigned long numQuad = 0;
	      for(size_t j=0;j<delta.size();j++)
		numQuad += delta[j][0];

	      ofstream outGrom;
	      string outputFileGrom = outputFile;
	      outputFileGrom.append("_gromov_collapse_");
	      outputFileGrom.append(ss.str());
	      outputFileGrom.append(".txt");
	      outGrom.open(outputFileGrom.c_str());

	      tt<<i;

	      ofstream outAvg;
	      string outputFileAvg = outputFile;
	      outputFileAvg.append("_avg_distr_collapse_");
	      outputFileAvg.append(tt.str());
	      outputFileAvg.append(".txt");
	      outAvg.open(outputFileAvg.c_str());
	      outGrom<<"#Column format: \n";
	      outGrom<<"#Max Length in Quadruplet\tNumber of Quadruplets\n";

	      for(size_t j=0;j<delta.size();j++)
		{
		  outGrom<<"#~\t"<<j+1<<"\t"<<delta[j][0]<<"\n";
		}

	      outGrom<<"#Column format: \n";
	      outGrom<<"#Delta\tPercent of quadruplets for each length\n";

	      for(size_t j=1;j<delta[0].size();j++)
		{
		  //size_t size = delta[i].size();	

	  
		  outGrom<<(j-1)/2.0<<"\t";
		  for(size_t k=0;k<delta.size();k++)
		    {
		      if(delta[k][0]!=0)
			(delta[k])[j] = (delta[k])[j]/(delta[k])[0];
	  
	      
		      outGrom<<delta[k][j]<<"\t";
		    }
		  outGrom<<"\n";
	  
		}

	      for(size_t k=0;k<delta.size();k++)
		{
		  double avg = 0;
		  for(size_t j=1;j<delta[k].size();j++)
		    {
		      avg += delta[k][j]*(j-1)/2.0;
		    }
		  outAvg<<k+1<<"\t"<<avg<<"\n";
		}

	      outGrom.close();
	      outAvg.close();
	    }
	}
      ////////////////////////////////////////////////////////

    }
    

  file_periph.close();
  file_collapse.close();
  file_core.close();
    
    
  return(0);
}
