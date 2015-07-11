/*
This function calculates the ncp plot over successive k_core
decompositions.  If the k_core has not been provided, the function will
find the k_core decomposition of the graph.
By Aaron Adcock, PhD Candidate, Stanford University
Sept, 2011
 */

#include "../graph_lib_boost.hpp"

void k_core_ncp(Graph& G,const string ncpName)
{
  vector<int> core = k_core(G);
  const double alpha = .001;
  const int step = 10;
  k_core_ncp(G,ncpName,core,alpha,step);
}

void k_core_ncp(Graph& G,const string ncpName,vector<int>& core,const double alpha,const int step)
{
  vector<Vert> ind(core.size(),0);
  int max_core = 0;

  for(size_t i=0;i<core.size();i++)
    ind[i]=i;

  sort(ind.begin(),ind.end(),compareIndltg<vector<int>&>(core));

  for(size_t i=0; i<core.size();i++)
    {      
      if(core[i]>max_core)
	max_core = core[i];
    }

  long maxClustSize = 0;
  vector< vector<double> > Cdata;
  for(int i=1; i<=max_core; i++)
    {
      vector<double> C;
      vector<Vert> S;
      stringstream ns;    

      while(core[ind[0]]<i)
	{
	  //	  cout<<core[ind[0]]<<"\n";
	  ind.erase(ind.begin());
	}

      S = ind;

      Graph G_k = subset(G,S);
      G_k = connected(G_k,-1);
      
      C = ncp_calc(G_k,num_vertices(G_k)/2,step,alpha);
      // for(int j=1;j<C.size();j++)
      // 	 cout<<C[j]<<"\t";
      // cout<<"ncp calculated\n";
      // cout<<"\n"<<"Vector of Vector\n";
      Cdata.push_back(C);
  
      // for(int j=0; j<Cdata[i-1].size();j++)
      // 	cout<<Cdata[i-1][j]<<"\t";
      // cout<<"NCP plotted\n";
      if(C.size()>maxClustSize)
	maxClustSize = C.size();
      
    }
  

  ofstream graphFile;
  graphFile.open(ncpName.c_str());
  graphFile.precision(10);
  graphFile<<scientific;

  for(int i=0;i<maxClustSize;i++)
    {
      graphFile<<i+1<<"\t";
      for(size_t j=0;j<Cdata.size();j++)
	{
	  vector<double> C;
	  C.assign(Cdata[j].begin(),Cdata[j].end());
	  if(C.size()>i)
	    graphFile<<C[i]<<"\t";
	  else
	    graphFile<<"\t";//C[C.size()-1]<<"\t";
	}
      graphFile<<"\n";
    }
  graphFile.close();
}
