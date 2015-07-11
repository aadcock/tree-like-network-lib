/*This function calculates the statistics of nodes using the
  k_core_decomposition of a graph

By Aaron Adcock, PhD candidate at Stanford University
August 2011*/

#include <vector>
#include <utility>
#include "../graph_lib_boost.hpp"
#include <iomanip>
#include <stdio.h>
#include <fstream>


void k_core_stat(vector<vector<float> >& core_stat, Graph & G)
{
  k_core_stat(core_stat,G,"");
}

void k_core_stat(vector<vector<float> >& core_stat, Graph & G, const string output_file_prefix)
{
  //basic parameters/graph stats for algorithm
  int size = num_vertices(G);
  int e_size = num_edges(G);
  vector<int> core = k_core(G);
  int max_core = 0;
  int num_bins = 11;

  //statistic vectors
  vector< vector<int> > jumpStat;  //edges,edgesin,edges0,edges5,edges10,edges20,edges20+
  vector< pair<int,int> > stat; //running totals, (avg,var)

  //maximum k_core value
  for(size_t i=0; i<core.size();i++)
    {      
      if(core[i]>max_core)
	max_core = core[i];
    }

  int min_core = max_core;

  for(size_t i=0; i<core.size();i++)
    {
      if(core[i]<min_core)
	min_core = core[i];
    }

  //core distribution (ie, how many nodes at each core)
  vector<int> core_dist(max_core+1,0);
  for(int i=0;i<core.size();i++)
    core_dist[(core[i])]++;


  //initialize vectors
  for(int i=0; i<=max_core; i++)
    {
      vector<int> dummy(num_bins,0);
      pair<int,int> p;

      p.first  = 0;
      p.second = 0;
      jumpStat.push_back(dummy);
      stat.push_back(p);
    }


  //Get index vector and create vertex iterators
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  graph_traits<Graph>::vertex_iterator vit, vitend;

  double avg_degree = 0;

  //Iterate through vertices
  for(tie(vit,vitend)=vertices(G);vit!=vitend;++vit)
    {
      Vert v = *vit;
      int ind   = index[v];
      int c_val = core[ind];
      
      avg_degree += out_degree(v,G);

      graph_traits<Graph>::out_edge_iterator out_i, out_e;

      //iterate over outgoing edges of current node
      for(tie(out_i,out_e)=out_edges(v,G);out_i!=out_e;++out_i)
	{
	  Edge e  = *out_i;                           //current edge
	  Vert vt = target(e,G);                      //target node (neighbor of v)
	  
	  ind = index[vt];
	  int vc_val = core[ind];                    //core value of vt
	  int jump   = vc_val-c_val;                 //jump from current core

	  jumpStat[c_val][0]++;

	  stat[c_val].first  += jump;
	  stat[c_val].second += jump*jump;
	  
	  //jump>0 => a jump 'in',  larger core number
	  //jump<0 => a jump 'out', smaller core number

	  if(jump>0)
	    {
	      jumpStat[c_val][1]++;
	     
	      if(jump <= 5)
		jumpStat[c_val][3]++;
	      else if(jump <= 10)
		jumpStat[c_val][5]++;
	      else if(jump <= 20)
		jumpStat[c_val][7]++;
	      else
		jumpStat[c_val][9]++;	      
	    }

	  if(jump<0)
	    {
	      if(jump >= -5)
		jumpStat[c_val][4]++;
	      else if(jump >= -10)
		jumpStat[c_val][6]++;
	      else if(jump >= -20)
		jumpStat[c_val][8]++;
	      else
		jumpStat[c_val][10]++;	      
	    }
	  
	  if(jump == 0)
	    jumpStat[c_val][2]++;
	  
	}
    }
  
  avg_degree = avg_degree/double(size);
  vector<vector<float> > final_stat;

  for(int i=0; i<=max_core; i++)
    {
      vector<float> temp(14,0);
      
      float total = float((jumpStat[i][0]));

      if(total!=0)
	{
	  temp[0]  = total;
	  temp[1]  = 100*jumpStat[i][1]  /total;
	  temp[2]  = stat[i].first       /total;
	  temp[3]  = stat[i].second      /total - temp[3]*temp[3];
	  temp[4]  = 100*jumpStat[i][2]  /total;
	  temp[5]  = 100*jumpStat[i][3]  /total;
	  temp[6]  = 100*jumpStat[i][4]  /total;
	  temp[7]  = 100*jumpStat[i][5]  /total;
	  temp[8]  = 100*jumpStat[i][6]  /total;
	  temp[9]  = 100*jumpStat[i][7]  /total;
	  temp[10] = 100*jumpStat[i][8]  /total;
	  temp[11] = 100*jumpStat[i][9]  /total;
	  temp[12] = 100*jumpStat[i][10] /total;
	  temp[13] = core_dist[i];
	}
      final_stat.push_back(temp);
    }
  core_stat.clear();
  core_stat.insert(core_stat.begin(),final_stat.begin(),final_stat.end());

  if(output_file_prefix != "")
    {
      //Print Results
      ofstream graph2File;
      string file2 = output_file_prefix;
      file2.append("_core_dist.txt");
      graph2File.open(file2.c_str());

      

      if(graph2File.is_open())
	{
	  for(int i=0; i<core_dist.size(); i++)
	    {
	      graph2File <<i<<"\t"<<core_dist[i]<<"\n";
	    }
	}
      graph2File.close();

      file2 = output_file_prefix;
      file2.append("_core_latex_dist.txt");
      graph2File.open(file2.c_str());

      if(graph2File.is_open())
	{
	  int num_cores   = core_dist.size() - 1;
	  int core_25_per_ind = floor(num_cores/4.0 + 0.5);
	  int core_50_per_ind = floor(num_cores/2.0 + 0.5);
	  int core_75_per_ind = floor(3.0*num_cores/4.0 + 0.5);

	  int core_25_per = 0;
	  int core_50_per = 0;
	  int core_75_per = 0;

	  if(core_25_per_ind == 0)
	    ++core_25_per_ind;
	  
	  if(core_50_per_ind == 0)
	    ++core_50_per_ind;

	  if(core_75_per_ind == 0)
	    ++core_75_per_ind;

	  for(int i=num_cores;i>core_25_per_ind;--i)
	    core_25_per += core_dist[i];

	  for(int i=num_cores;i>core_50_per_ind;--i)
	    core_50_per += core_dist[i];

	  for(int i=num_cores;i>core_75_per_ind;--i)
	    core_75_per += core_dist[i];

	  
	  graph2File<<output_file_prefix<<" & "<<size<<" & "<<setprecision(4)<<avg_degree<<" & "<<min_core<<" & "<<max_core<<" & "<<100*double(core_dist[min_core])/double(size)<<" & "<<100*double(core_dist[num_cores])/double(size)<<" \\\\ \n";
	}
      else
	{
	  cerr<<"Error opening "<<file2<<", file not produced\n";
	}

      ofstream graph3File;
      ofstream graphJumpInFile;
      ofstream graphJumpOutFile;

      string filein  = output_file_prefix;
      filein.append("_jump_in_output.txt");

      string fileout = output_file_prefix;
      fileout.append("_jump_out_output.txt");

      string file3   = output_file_prefix;
      file3.append("_core_stat_output.txt");

      graph3File.open(file3.c_str());
      graphJumpInFile.open(filein.c_str());
      graphJumpOutFile.open(fileout.c_str());

      if(graph3File.is_open() && graphJumpInFile.is_open() && graphJumpOutFile.is_open())
	{
	  graph3File<<"#Nodes: "<<size<<"\n";
	  graph3File<<"#Edges: "<<e_size<<"\n";
	  graph3File<<"#Core #\tEdges\tEdgesin\tAvg\tVar\t0\t1-5in\t1-5out\t5-10in\t5-10out\t10-20in\t10-20ou\t>20in\t>20ou\tdistr\n";
	  graph3File<<fixed;
	  for(int i=0;i<core_stat.size();i++)
	    {
	      vector<float> stat = core_stat[i];
	      graph3File<<i<<": \t";
	      for(int j=0;j<stat.size();j++)
		{
		  graph3File<<setprecision(1)<<stat[j]<<"\t";
		}
	      graph3File<<"\n";
	    }

	  for(int i=0;i<5;i++)
	    {
	      graphJumpInFile<<i;
	      graphJumpOutFile<<i;
	      for(int j=0; j<core_stat.size();j++)
		{
		  graphJumpOutFile<<"\t"<<core_stat[j][4+2*i];
		  if(i==0)
		    graphJumpInFile<<"\t"<<core_stat[j][4];
		  else
		    graphJumpInFile<<"\t"<<core_stat[j][5+(i-1)*2];

		}
	      graphJumpInFile<<"\n";
	      graphJumpOutFile<<"\n";
	    }
	}
      graph3File.close();
      graphJumpInFile.close();
      graphJumpOutFile.close();
    }
  else
    {
      cout<<"Output file prefix was empty. \n";
    }

}
