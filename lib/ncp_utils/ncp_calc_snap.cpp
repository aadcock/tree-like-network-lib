/*This function tries to repeat the ncp_profile plots from the "Empirical Comparison of Algorithms for Network Community Detection" by Leskovec, Lang, Mahoney.  Uses a spectral method based on personalized Page Rank vectors and sweep cuts (see approx_Page_Rank_boost.cpp for more information.

By Aaron Adcock, PhD Candidate at Stanford University
July 2011 */


#include "../graph_lib_boost.hpp"
#include <cstdlib>
#include <ctime>
#include <vector>



//Main function for finding minimum conductance communities across various size scales.  This function just keeps track of the minimum conductance, not the actual communities.  This could be changed using the index vectors below.

/*
G:          Graph under investigation
maxC:       Maximum community size to look at (ie how far to take sweep cut)
iter:       Number of Page Rank iterations with random seed nodes to do
alpha:      Teleport probability for personalized Page Rank
eps:        Approximation factor for personalized Page Rank (smaller = more accuracy)
 */
vector<float> ncp_calc_snap(Graph& G, const int maxC, const int iter, const float alpha, const float eps)
{

    //Pick a seed for random number generator (just used standard c generator, could use better random generator)
    srand(time(NULL));
    int size = num_vertices(G);
    vector<float> finalCond(maxC,1.0);
    vector<int> ind(size,0);
    graph_traits<Graph>::edges_size_type edgeCount = num_edges(G)/2;  //The 2 is because this is represented as a directed graph, each undirected edge represented by 2 directed edges

    //Page Rank iterations, each time a sweep cut is performed to look at different sizes of communities

    for(int csize=1; csize<maxC+1;csize++){
      float c = 1;
      for(int i=0; i<size-csize; i++)
	{
       
	  vector<int> indOrig (size,0);       //Keep track of nodes after sort
	  int r = rand() % size;              //index of seed node
	  vector<float> s (size,0);           //Seed vector (a 1 at node r)
	  vector<float> p;                    //Store final page rank vector
	
	  //Index vector that will be sorted by Page Rank vector
	  for(int j=0;j<size;j++)
	    ind[j] = j;

	  s[r] = 1;
	
	  //Calculate Page Rank (time variable for roughly assessing time of algorithm)
	  long time1 = time(NULL);
	  p = approxPR_boost(G,s,alpha,1/double(csize+1));
	  long time2 = time(NULL);
	  //	cout<<"\nIter "<<i<<": "<<time2-time1<<"\n";

	  //Normalize PR vector by outdegree of each vertex
	  for(int j=0; j<size; j++)
	    { 
	      Vert v = vertex(j,G);
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);
	      p[j] = p[j]/outDegree;
	    }

	  //Sort ind vector based on normalized Page Rank values
	  time1 = time(NULL);
	  sort(ind.begin(),ind.end(),compareInd<vector<float>&>(p));
	  time2 = time(NULL);
	
	  //indOrig is now indexed by the original vertex index and points to the rank of the vertex
	  for(int j=0;j<size;j++)
	    indOrig[ind[j]] = j;
	  //	cout<<"Time of Sort: "<<time2-time1<<"\n";



	  //Calculate conductance of each cut along PR vector
	  long vol = 0;
	  long out = 0;
	  //	cout<<"Iter: "<<i<<"\n";

	  //Do Sweep cut up to maximum community size, calculate conductance as we add vertices to cut
	 
	  for(int j=0; j<maxC; j++)
	    {
	      //   cout<<"Comm: "<<j<<"\n";
	      vector<int>::iterator it;
	   
	      Vert v    = vertex(ind[j],G);
	      //  cout<<"index: "<<index<<"\n";

	      //Out edge iterators
	      graph_traits<Graph>::out_edge_iterator out_i, out_e;

	      //Create index map for Graph
	      property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

	      //Count outgoing edges for v
	      graph_traits<Graph>::degree_size_type outDegree = out_degree(v,G);

	      //update edges leaving cut set (out)
	      out += outDegree;
	      //update internal volume of cut set (vol)
	      vol += outDegree;

	      //Iterate over out going edges of v
	      for(tie(out_i,out_e)=out_edges(v,G);out_i!=out_e;++out_i)
		{
		  Edge e = *out_i;
		  //	cout<<"node: "<<node1<<"\n";

		  Vert vt = target(e,G);

		  //Check if edge crosses set boundary, if so, update out and vol accordingly
		  if(indOrig[index[vt]]<j && indOrig[index[vt]]>=0)
		    {
		      out -= 2;  //note that 2 is due to outgoing/incoming edge as G is a directed graph representation of an undirected graph.
		      //		    vol--;
		    }
		}
	      //cout<<edgeCount<<"\n";
	
	      //Pick minimum of cut set and complement volume
	      if(vol>(edgeCount-vol))
		vol = edgeCount - vol;
	      // cout<<"Comm size "<<j<<": "<<ctime2-ctime1<<"\n";
	    
	    
	      //conductance c is now given by: c = out/vol
	      float ctemp = float(out)/vol;
	    
	      if(ctemp<finalCond[j])
		finalCond[j] = ctemp;
	      //Check if this beats previous minimum
	   
	    }
	  long ctime2 = time(NULL);
	  //	cout<<"Conductance Calc for 200: "<<ctime2-ctime1<<"\n";
	}

    }

    return finalCond;
  }
