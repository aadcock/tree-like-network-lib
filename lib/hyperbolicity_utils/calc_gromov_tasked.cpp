/*
This file contains a function for computing the Gromov 'delta' hyperbolicity of the given, connected, graph G.  It is meant to be easily parallelizable.

By Aaron Adcock, PhD candidate at Stanford University, Dec. 2011
 */

#include "../graph_lib_boost.hpp"
#include <omp.h>


//Calculates max delta of triangles with sides uv, uw, vw
//double triangle_delta(Graph& G,vector<Vert>& uv, vector<Vert>& uw, vector<Vert>& vw);


struct path_node
{
  Vert v;
  v_size_t predecessors;
  unsigned int count;
};


//double pair_delta(Graph& G, const vector<Vert>& side,const vector<Vert>& other_side1,const vector<Vert>& other_side2,const vector<vector<unsigned int> >& path_ind1,const vector<vector<unsigned int> >& path_ind2,const unsigned int& num_paths1, const unsigned int& num_paths2,const property_map<Graph,vertex_index_t>::type& index);

double pair_delta(vector<vector<v_size_t> >& d, const vector<Vert>& side, const vector<Vert>& other_side1, const vector<Vert>& other_side2, const vector<vector<unsigned int> >& path_ind1, const vector<vector<unsigned int> >& path_ind2, const unsigned int& num_paths1, const unsigned int& num_paths2, const property_map<Graph,vertex_index_t>::type& index, vector<unsigned int>& delta_v, vector<unsigned int>& k1, vector<unsigned int> k2, const size_t& n);

//index_paths indexes ALL shortest paths found from a BFS.  p is a list of predecessor vertices from the BFS
unsigned int index_paths(const Vert& s, const Vert& t,  const property_map<Graph,vertex_index_t>::type& index, const v_size_t& n, const vector< vector<Vert> >& p, vector< vector<unsigned int> >& ind, vector<path_node>& node_stack);

//Forms paths using predecessors, uses a vector for the queue to save on free() calls.
void form_path(const Vert& s, const Vert& t, const Graph& G, const vector< vector<Vert> >& predecessors, vector<Vert>& path, vector<Vert>& Q, vector<bool>& onPath,property_map<Graph,vertex_index_t>::type& index);

//Calculates the gromov four_point condition.  Excuse the large sections of commented code, this function's final form is still in flux.
vector< vector<double> > calc_gromov(Graph& G, v_size_t diam)
{
  ///////////////////////////
  //Delta data structures
  ///////////////////////////

  vector< vector<double> > delta;
  double maxDelta     = 0;

  if(num_vertices(G) < 4)
    {
      cerr<<"Graph is too small for hyperbolicity, returning zero\n";
      vector<double> delta_vec(1,0);
      delta.push_back(delta_vec);
      return delta;
    }

  Vert v = vertex(1,G);

  vector<bool> found(num_vertices(G),false);
  v_size_t current_comp = 0;
  v_size_t current_comp_size = 0;
  vector<v_size_t> component(num_vertices(G));
  
  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);

  vector<bool>::iterator it;
  for(it = found.begin();it<found.end();++it)
    {
      bool v_found = *it;
      if(!v_found)
	{
	  cerr<<"Graph passed to calc_gromov has multiple components\n";
	  G = connected(G);
	}
    }
  //This code was being used to make a global d matrix
  //  vector<Vert> vert_vec;
  //size_t size;
  //vector< vector<v_size_t> > d;

  // for(tie(vit,vite) = vertices(G); vit!= vite; vit++)
  //   {
  //     Vert v = *vit;
  //     vert_vec.push_back(v);
  //     //d.push_back(G[v].distances);
  //   }

  

  for(v_size_t i=0;i<diam;i++)
    {
      vector<double> zero_dist(2*diam+1,0);
      
      delta.push_back(zero_dist);
      //delta_loc.push_back(zero_dist);
    }
  //  cout<<"delta_loc: "<<delta_loc.size()<<"\n";
  //  size = vert_vec.size();


  //  vector< vector<double> > delta_loc_tmp;
  //double temp_d;

  //Determine number of threads in use
  int num_threads = omp_get_max_threads();

  //  vector< vector< vector<double> > > shared_delta(num_threads,delta_loc_tmp);
  //vector<double*> shared_max_delta(num_threads,&temp_d);

  //*****************
  //These are the shared variables between threads.  Each thread has access to a certain pointer.
  //The pointer should then point to a copy of the variable for the thread's private (local) use.
  //Shared max_delta and shared_delta are updated by each thread, must be careful to avoid contention.
  //*****************

  int size_of_array = 64/sizeof(double*);  //64 is cacheline size, needs to be updated on different machines for performance
  int type_size     = sizeof(double*);
  if(type_size<sizeof(vector<vector<double> >*))
    type_size = sizeof(vector<vector<double> >*);

  int counter = 1;

  while(num_threads*type_size>size_of_array)
    {
      ++counter;
      size_of_array = counter*64/sizeof(double*);
    }
  ////////////////////////////////////////////////////////////
  //Pointers to shared data arrays that have copy of the delta information for each thread,
  //ie, if there are 8 threads, there will be an 8 (pointers to) arrays each accessed by only one thread,
  //At the end of the parallel section of code, the 8 arrays will be combined into a single array
  ///////////////////////////////////////////////////////////

  double ** shared_max_delta;
  shared_max_delta = new double * [size_of_array];

  vector< vector<double> > ** shared_delta;
  shared_delta = new vector< vector<double> > * [size_of_array];

  //vector< vector<v_size_t> > ** shared_d;
  //  shared_d = new vector< vector<v_size_t> > * [size_of_array];

  //  property_map<Graph,vertex_index_t>::type ** shared_index;
  //shared_index = new property_map<Graph,vertex_index_t>::type * [size_of_array];

  //vector<Vert> ** shared_vert;
  //shared_vert = new vector<Vert> * [size_of_array];


  //OMP parallel region, this current implementation generates tasks at two of the four levels of the nested for loop, these tasks are then passed to threads which execute them and update their private copy of the results array
  //After the tasks have all been generated and executed the results are collated
#pragma omp parallel shared(delta,maxDelta,shared_max_delta,diam,shared_delta) //private(i,j,k,l,u,v,x,y,ui,vi,xi,yi,s1,s2,s3,currDelta,maxPair,delta_loc,maxDelta_loc)
  {

    //    size_t i,j,k,l;
    //Vert u,v,x,y;
    //v_size_t ui,vi,xi,yi;
    //v_size_t uv,ux,vx,uy,vy,xy;
    //v_size_t s1,s2,s3;
    //double currDelta = 0;
    //unsigned long maxPair = 1;
    double maxDelta_loc = 0;

    vector<double> zero_dist(2*diam+1,0);
    zero_dist.reserve(256/sizeof(double));

    vector< vector<double> > delta_loc(diam,zero_dist);

    int thread_id = omp_get_thread_num();

    //shared_delta[thread_id].swap(delta_loc);
    shared_max_delta[thread_id] = &maxDelta_loc;

    graph_traits<Graph>::vertex_iterator vit, vite;
    property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

    vector<Vert> vert_vec;
    vector< vector<v_size_t> > d;

    for(tie(vit,vite) = vertices(G); vit!= vite; vit++)
      {
	Vert v = *vit;
	vert_vec.push_back(v);
	d.push_back(G[v].distances);
      }

    shared_delta[thread_id] = &delta_loc;
    //    shared_d[thread_id] = &d;
    //shared_index[thread_id] = &index;
    //shared_vert[thread_id] = &vert_vec;

    size_t size = vert_vec.size();
    //#pragma omp for schedule(dynamic,1) //collapse(4) 

#pragma omp single
    {
      for(size_t i=0;i<size;++i)
	{
#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec,index)
	  {
	    for(size_t j=i+1;j<size;++j)
	      {
#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec,index)
		{
		  int thread_num = omp_get_thread_num();
		  
		  //		  vector<vector<double> > * delta_point = shared_delta[thread_num];
		  vector<vector<double> > & loc_delta = *shared_delta[thread_num];
		   
		  double & max_delta_loc = *shared_max_delta[thread_num];
		  
		  //		  vector<vector<v_size_t> > * d_point = shared_d[thread_num];
		  //vector<vector<v_size_t> > & loc_d = *d_point;
		  
		  //vector<Vert> * vert_point = shared_vert[thread_num];
		  //vector<Vert> & loc_vert_vec = *vert_point;
		  
		  //property_map<Graph,vertex_index_t>::type * index_point = shared_index[thread_num];
		  //property_map<Graph,vertex_index_t>::type & loc_index   = *index_point;
		  
		  for(size_t k=j+1;k<size;++k)
		    {
		      //#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec)
		      {
			for(size_t l=k+1;l<size;++l)
			  {
			    
			    
			    const Vert u = (vert_vec)[i];
			    const v_size_t ui = (index)[u];
			    //cerr<<"Vertex: "<<u<<"\n";
			    
			    const Vert v = (vert_vec)[j];
			    const v_size_t vi = (index)[v];
			    const v_size_t uv = (d)[ui][vi];//.distances[vi];
			    
			    const Vert x = (vert_vec)[k];
			    const v_size_t xi = (index)[x];
			    const v_size_t ux = (d)[ui][xi];//G[u].distances[xi];
			    const v_size_t vx = (d)[vi][xi];//G[v].distances[xi];
			    
			    const Vert y = (vert_vec)[l];		  
			    const v_size_t yi = (index)[y];
			    const v_size_t uy = (d)[ui][yi];//G[u].distances[yi];
			    const v_size_t vy = (d)[vi][yi];//G[v].distances[yi];
			    const v_size_t xy = (d)[xi][yi];//G[x].distances[yi];
			    const v_size_t s1 = uv + xy;
			    const v_size_t s2 = ux + vy;
			    const v_size_t s3 = uy + vx;
			    
			    double currDelta=0;
			    unsigned long maxPair=1;
			    
			    //find max pair distance
			    maxPair = uv;
			    
			    if(xy>maxPair)
			      maxPair = xy;
			    
			    if(ux>maxPair)
			      maxPair = ux;
			    
			    if(vy>maxPair)
			      maxPair = vy;
			    
			    if(uy>maxPair)
			      maxPair = uy;
			    
			    if(vx>maxPair)
			      maxPair = vx;
			    
			    //Delta calculations
		  
			    if(s1 >= s2 && s3 >= s2)
			      {
				if(s1 >= s3)
				  currDelta = (s1 - s3)/2.0;
				else
				  currDelta = (s3 - s1)/2.0;
			      }
			    else if(s2 >= s1 && s3 >= s1)
			      {
				if(s2 >= s3)
				  currDelta = (s2 - s3)/2.0;
				else
				  currDelta = (s3 - s2)/2.0;
			      }
			    else
			      {
				if(s1 >= s2)
				  currDelta = (s1 - s2)/2.0;
				else
				  currDelta = (s2 - s1)/2.0;
			      }

			
			    //++shared_delta[thread_num][maxPair-1][2*currDelta+1];
			    ++(loc_delta)[maxPair-1][2*currDelta+1];
			
			
			    if(currDelta>max_delta_loc)
			      max_delta_loc = currDelta;

			    //++shared_delta[thread_num][maxPair-1][0];
			    ++(loc_delta)[maxPair-1][0];
			
			  }      		  
		      }
		    }
		}
	      }
	  }
	}
    }
#pragma omp critical (collate_delta) 
    {
      int thread_num = omp_get_thread_num();
      vector<vector<double> > * delta_point = shared_delta[thread_num];
      vector<vector<double> > & loc_delta = *delta_point;

      size_t i_size  = loc_delta.size();
      //cout<<delta_loc.size()<<"\n";
      for(size_t i=0;i<i_size;++i)
	{
	  size_t j_size = loc_delta[i].size();
	  for(size_t j=0;j<j_size;++j)
	    {
	      //    cout<<i<<" "<<j<<"\n";
	      delta[i][j] += loc_delta[i][j];
	    }
	}

      if(maxDelta_loc > maxDelta)
	maxDelta = maxDelta_loc;
    }
  } 

  //take care of pointer!
  delete [] shared_delta;
  //delete [] shared_d;
  //  delete [] shared_index;
  // delete [] shared_vert;
  delete [] shared_max_delta;

  //Needs removed before parallel code is implemented
  for(v_size_t i=0;i<diam;++i)
    {
      delta[i].erase(delta[i].begin()+2*maxDelta+2,delta[i].end());
    }

  return(delta);
}


//Calculates the gromov four_point condition.  Excuse the large sections of commented code, this function's final form is still in flux.
vector< vector<double> > calc_gromov_quads(Graph& G, v_size_t diam, int quads_to_keep, int dist_keep, string quad_output)
{
  srand(time(NULL));
  vector<int> core = k_core(G);
  ///////////////////////////
  //Delta data structures
  ///////////////////////////

  vector< vector<double> > delta;
  double maxDelta     = 0;

  if(num_vertices(G) < 4)
    {
      cerr<<"Graph is too small for hyperbolicity, returning zero\n";
      vector<double> delta_vec(1,0);
      delta.push_back(delta_vec);
      return delta;
    }

  Vert v = vertex(1,G);

  vector<bool> found(num_vertices(G),false);
  v_size_t current_comp = 0;
  v_size_t current_comp_size = 0;
  vector<v_size_t> component(num_vertices(G));
  
  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);

  vector<bool>::iterator it;
  for(it = found.begin();it<found.end();++it)
    {
      bool v_found = *it;
      if(!v_found)
	{
	  cerr<<"Graph passed to calc_gromov has multiple components\n";
	  G = connected(G);
	}
    }
  //This code was being used to make a global d matrix
  //  vector<Vert> vert_vec;
  //size_t size;
  //vector< vector<v_size_t> > d;

  // for(tie(vit,vite) = vertices(G); vit!= vite; vit++)
  //   {
  //     Vert v = *vit;
  //     vert_vec.push_back(v);
  //     //d.push_back(G[v].distances);
  //   }

  

  for(v_size_t i=0;i<diam;i++)
    {
      vector<double> zero_dist(2*diam+1,0);
      
      delta.push_back(zero_dist);
      //delta_loc.push_back(zero_dist);
    }
  //  cout<<"delta_loc: "<<delta_loc.size()<<"\n";
  //  size = vert_vec.size();


  //  vector< vector<double> > delta_loc_tmp;
  //double temp_d;

  //Determine number of threads in use
  int num_threads = omp_get_max_threads();

  vector<int> zero_threads(num_threads,0);
  vector< vector<int> > quads_kept(diam,zero_threads);
    
  //  vector< vector< vector<double> > > shared_delta(num_threads,delta_loc_tmp);
  //vector<double*> shared_max_delta(num_threads,&temp_d);

  //*****************
  //These are the shared variables between threads.  Each thread has access to a certain pointer.
  //The pointer should then point to a copy of the variable for the thread's private (local) use.
  //Shared max_delta and shared_delta are updated by each thread, must be careful to avoid contention.
  //*****************

  int size_of_array = 64/sizeof(double*);  //64 is cacheline size, needs to be updated on different machines for performance
  int type_size     = sizeof(double*);
  if(type_size<sizeof(vector<vector<double> >*))
    type_size = sizeof(vector<vector<double> >*);

  int counter = 1;

  while(num_threads*type_size>size_of_array)
    {
      ++counter;
      size_of_array = counter*64/sizeof(double*);
    }
  ////////////////////////////////////////////////////////////
  //Pointers to shared data arrays that have copy of the delta information for each thread,
  //ie, if there are 8 threads, there will be an 8 (pointers to) arrays each accessed by only one thread,
  //At the end of the parallel section of code, the 8 arrays will be combined into a single array
  ///////////////////////////////////////////////////////////

  double ** shared_max_delta;
  shared_max_delta = new double * [size_of_array];

  vector< vector<double> > ** shared_delta;
  shared_delta = new vector< vector<double> > * [size_of_array];

  //vector< vector<v_size_t> > ** shared_d;
  //  shared_d = new vector< vector<v_size_t> > * [size_of_array];

  //  property_map<Graph,vertex_index_t>::type ** shared_index;
  //shared_index = new property_map<Graph,vertex_index_t>::type * [size_of_array];

  //vector<Vert> ** shared_vert;
  //shared_vert = new vector<Vert> * [size_of_array];


  //OMP parallel region, this current implementation generates tasks at two of the four levels of the nested for loop, these tasks are then passed to threads which execute them and update their private copy of the results array
  //After the tasks have all been generated and executed the results are collated
#pragma omp parallel shared(delta,maxDelta,shared_max_delta,diam,shared_delta) //private(i,j,k,l,u,v,x,y,ui,vi,xi,yi,s1,s2,s3,currDelta,maxPair,delta_loc,maxDelta_loc)
  {

    //    size_t i,j,k,l;
    //Vert u,v,x,y;
    //v_size_t ui,vi,xi,yi;
    //v_size_t uv,ux,vx,uy,vy,xy;
    //v_size_t s1,s2,s3;
    //double currDelta = 0;
    //unsigned long maxPair = 1;
    double maxDelta_loc = 0;

    vector<double> zero_dist(2*diam+1,0);
    zero_dist.reserve(256/sizeof(double));

    vector< vector<double> > delta_loc(diam,zero_dist);

    int thread_id = omp_get_thread_num();

    //shared_delta[thread_id].swap(delta_loc);
    shared_max_delta[thread_id] = &maxDelta_loc;

    graph_traits<Graph>::vertex_iterator vit, vite;
    property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

    vector<Vert> vert_vec;
    vector< vector<v_size_t> > d;

    for(tie(vit,vite) = vertices(G); vit!= vite; vit++)
      {
	Vert v = *vit;
	vert_vec.push_back(v);
	d.push_back(G[v].distances);
      }

    shared_delta[thread_id] = &delta_loc;
    //    shared_d[thread_id] = &d;
    //shared_index[thread_id] = &index;
    //shared_vert[thread_id] = &vert_vec;

    size_t size = vert_vec.size();
    //#pragma omp for schedule(dynamic,1) //collapse(4) 
    

#pragma omp single
    {
      for(size_t i=0;i<size;++i)
	{
#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec,index)
	  {
	    for(size_t j=i+1;j<size;++j)
	      {
#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec,index)
		{
		  int thread_num = omp_get_thread_num();
		  int delta_count = 0;
		  //		  vector<vector<double> > * delta_point = shared_delta[thread_num];
		  vector<vector<double> > & loc_delta = *shared_delta[thread_num];
		   
		  double & max_delta_loc = *shared_max_delta[thread_num];
		  
		  //		  vector<vector<v_size_t> > * d_point = shared_d[thread_num];
		  //vector<vector<v_size_t> > & loc_d = *d_point;
		  
		  //vector<Vert> * vert_point = shared_vert[thread_num];
		  //vector<Vert> & loc_vert_vec = *vert_point;
		  
		  //property_map<Graph,vertex_index_t>::type * index_point = shared_index[thread_num];
		  //property_map<Graph,vertex_index_t>::type & loc_index   = *index_point;
		  
		  for(size_t k=j+1;k<size;++k)
		    {
		      //#pragma omp task shared(shared_delta,shared_max_delta,d,vert_vec)
		      {
			for(size_t l=k+1;l<size;++l)
			  {
			    
			    
			    const Vert u = (vert_vec)[i];
			    const v_size_t ui = (index)[u];
			    //cerr<<"Vertex: "<<u<<"\n";
			    
			    const Vert v = (vert_vec)[j];
			    const v_size_t vi = (index)[v];
			    const v_size_t uv = (d)[ui][vi];//.distances[vi];
			    
			    const Vert x = (vert_vec)[k];
			    const v_size_t xi = (index)[x];
			    const v_size_t ux = (d)[ui][xi];//G[u].distances[xi];
			    const v_size_t vx = (d)[vi][xi];//G[v].distances[xi];
			    
			    const Vert y = (vert_vec)[l];		  
			    const v_size_t yi = (index)[y];
			    const v_size_t uy = (d)[ui][yi];//G[u].distances[yi];
			    const v_size_t vy = (d)[vi][yi];//G[v].distances[yi];
			    const v_size_t xy = (d)[xi][yi];//G[x].distances[yi];
			    const v_size_t s1 = uv + xy;
			    const v_size_t s2 = ux + vy;
			    const v_size_t s3 = uy + vx;
			    
			    double currDelta=0;
			    unsigned long maxPair=1;
			    
			    //find max pair distance
			    maxPair = uv;
			    
			    if(xy>maxPair)
			      maxPair = xy;
			    
			    if(ux>maxPair)
			      maxPair = ux;
			    
			    if(vy>maxPair)
			      maxPair = vy;
			    
			    if(uy>maxPair)
			      maxPair = uy;
			    
			    if(vx>maxPair)
			      maxPair = vx;
			    
			    //Delta calculations
		  
			    if(s1 >= s2 && s3 >= s2)
			      {
				if(s1 >= s3)
				  currDelta = (s1 - s3)/2.0;
				else
				  currDelta = (s3 - s1)/2.0;
			      }
			    else if(s2 >= s1 && s3 >= s1)
			      {
				if(s2 >= s3)
				  currDelta = (s2 - s3)/2.0;
				else
				  currDelta = (s3 - s2)/2.0;
			      }
			    else
			      {
				if(s1 >= s2)
				  currDelta = (s1 - s2)/2.0;
				else
				  currDelta = (s2 - s1)/2.0;
			      }


			    int delta_ind   = int(currDelta*2);
			    double r; 

			    if(delta_count<quads_to_keep && maxPair>=dist_keep)
			      {
				r = 0.0;//(double)rand()/RAND_MAX;
				delta_count = 0;
				for(int thr_id=0;thr_id < quads_kept[delta_ind].size();++thr_id)
				  {
				    delta_count += quads_kept[delta_ind][thr_id];
				  }
			      }
			    
			      
			    vector<Vert> quadruped_nodes;
			    if(delta_count<quads_to_keep && maxPair>=dist_keep && r<0.5)
			      {
				vector<Vert> u_v;
				vector<Vert> u_x;
				vector<Vert> u_y;
				vector<Vert> v_x;
				vector<Vert> v_y;
				vector<Vert> x_y;
				vector< vector<Vert> > pred;

				BFS_three_target(u,v,x,y,G,u_v,u_x,u_y,pred);
				BFS_two_target(v,x,y,G,v_x,v_y,pred);
				BFS_source_target(x,y,G,x_y,pred);

				vector<Vert> master = u_v;
				
				master.insert(master.begin(),u_x.begin(),u_x.end());
				master.insert(master.begin(),u_y.begin(),u_y.end());
				master.insert(master.begin(),v_x.begin(),v_x.end());
				master.insert(master.begin(),v_y.begin(),v_y.end());
				master.insert(master.begin(),x_y.begin(),x_y.end());

				Graph Gnew = subset(G,master);

				string quad_output_curr = quad_output;


				char delt_char[40];
				quad_output_curr.append("_delt_");
				snprintf(delt_char,sizeof(delt_char),"%d",delta_ind);
				quad_output_curr.append(delt_char);

				char dist_char[40];
				quad_output_curr.append("_dist_");
				snprintf(dist_char,sizeof(dist_char),"%d",maxPair);
				quad_output_curr.append(dist_char);

				char u_char [40];
				quad_output_curr.append("_");
				snprintf(u_char, sizeof(u_char),"%d",ui);
				quad_output_curr.append(u_char);
				
				char v_char[40];
				quad_output_curr.append("_");
				snprintf(v_char, sizeof(v_char),"%d",vi);
				quad_output_curr.append(v_char);
				
				char x_char[40];
				quad_output_curr.append("_");
				snprintf(x_char, sizeof(x_char),"%d",xi);
				quad_output_curr.append(x_char);

				char y_char[40];
				quad_output_curr.append("_");
				snprintf(y_char, sizeof(y_char),"%d",yi);
				quad_output_curr.append(y_char);

				string quad_orig_indices = quad_output_curr;
				quad_orig_indices.append("_orig_indices.txt");
				
				string quad_orig_core    = quad_output_curr;
				quad_orig_core.append("_orig_k_core.color");
				quad_output_curr.append(".dimacs");
				vector<int> prev_indices;
				
				graph_traits<Graph>::vertex_iterator newit, newitend;

				for(tie(newit,newitend)=vertices(Gnew);newit!=newitend;++newit)
				  {
				    Vert new_v = *newit;
				    int ind_ind = Gnew[new_v].prevIndices.size()-1;
				    int old_ind = Gnew[new_v].prevIndices[ind_ind];
				    prev_indices.push_back(old_ind);
				    
				  }
				vector <int> core_vec = prev_indices;

				for(int ii=0;ii<core_vec.size();++ii)
				  {
				    core_vec[ii] = core[prev_indices[ii]];
				  }
				++quads_kept[delta_ind][thread_num];
				write_dimacs(Gnew,quad_output_curr);
				write_color_file(quad_orig_indices,quad_output,prev_indices);
				write_color_file(quad_orig_core,quad_output,core_vec);
			      }
			      

			    //++shared_delta[thread_num][maxPair-1][2*currDelta+1];
			    ++(loc_delta)[maxPair-1][2*currDelta+1];
			
			
			    if(currDelta>max_delta_loc)
			      max_delta_loc = currDelta;

			    //++shared_delta[thread_num][maxPair-1][0];
			    ++(loc_delta)[maxPair-1][0];
			
			  }      		  
		      }
		    }
		}
	      }
	  }
	}
    }
#pragma omp critical (collate_delta) 
    {
      int thread_num = omp_get_thread_num();
      vector<vector<double> > * delta_point = shared_delta[thread_num];
      vector<vector<double> > & loc_delta = *delta_point;

      size_t i_size  = loc_delta.size();
      //cout<<delta_loc.size()<<"\n";
      for(size_t i=0;i<i_size;++i)
	{
	  size_t j_size = loc_delta[i].size();
	  for(size_t j=0;j<j_size;++j)
	    {
	      //    cout<<i<<" "<<j<<"\n";
	      delta[i][j] += loc_delta[i][j];
	    }
	}

      if(maxDelta_loc > maxDelta)
	maxDelta = maxDelta_loc;
    }
  } 

  //take care of pointer!
  delete [] shared_delta;
  //delete [] shared_d;
  //  delete [] shared_index;
  // delete [] shared_vert;
  delete [] shared_max_delta;

  //Needs removed before parallel code is implemented
  for(v_size_t i=0;i<diam;++i)
    {
      delta[i].erase(delta[i].begin()+2*maxDelta+2,delta[i].end());
    }

  return(delta);
}




vector< vector<double> > calc_delta_slim(Graph& G, v_size_t diam)
{



  //Variable initializations
  const vector<double> zero_dist(diam+2,0);
  double maxDelta = 0;


  // for(tie(vit,vite)=vertices(G); vit!=vite; vit++)
  //   {
  //     Vert v =*vit;
  //     vert_vec.push_back(v);
  //   }


  vector< vector<double> > delta(diam,zero_dist);


  if(num_vertices(G) < 3)
    {
      cerr<<"Graph is too small for hyperbolicity, returning zero\n";
      return delta;
    }

  Vert v = vertex(1,G);

  vector<bool> found(num_vertices(G),false);
  v_size_t current_comp = 0;
  v_size_t current_comp_size = 0;
  vector<v_size_t> component(num_vertices(G));
  
  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);

  vector<bool>::iterator it;
  for(it = found.begin();it<found.end();++it)
    {
      bool v_found = *it;
      if(!v_found)
	{
	  cerr<<"Graph passed to calc_delta_slim has multiple components\n";
	  G = connected(G);
	}
    }

#pragma omp parallel shared(delta,maxDelta,diam)
  {

    //Variable assignements,
    //Notation: puv = "Path indices on side uv"
    //           uv = "Vertices on side uv"
    //vertex iterators
    graph_traits<Graph>::vertex_iterator vit, vite;
  
    //Property maps
    property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

    vector<Vert> vert_vec;
    vector< vector<v_size_t> > d;

    for(tie(vit,vite) = vertices(G); vit!= vite; ++vit)
      {
	Vert v = *vit;
	vert_vec.push_back(v);
	d.push_back(G[v].distances);
      }

    const size_t size = vert_vec.size();

    vector<unsigned int> dummy;
    vector< vector<unsigned int> > puv(size,dummy);
    vector< vector<unsigned int> > puw(size,dummy);
    vector< vector<unsigned int> > pvw(size,dummy);

    vector<path_node> node_stack;
    node_stack.reserve(CHUNKSIZE);

    for(size_t i=0;i<size;i++)
      {
    	puv[i].reserve(CHUNKSIZE);
    	puw[i].reserve(CHUNKSIZE);
    	pvw[i].reserve(CHUNKSIZE);
	
      }
    
    vector<unsigned int> delta_v;
    vector<unsigned int> k1;
    vector<unsigned int> k2;

    delta_v.reserve(diam+1);
    k1.reserve(CHUNKSIZE);
    k2.reserve(CHUNKSIZE);

    vector<Vert> uv,uw,vw;
    vector<bool> on_path;
    queue<Vert> Q;

    on_path.reserve(size);
    uv.reserve(diam+1);
    uw.reserve(diam+1);
    vw.reserve(diam+1);

    vector<vector<Vert> > vPredecessors;
    vector<vector<Vert> > uPredecessors;

    vector< vector<double> > delta_loc(diam,zero_dist);
    double maxDelta_loc = 0;

    double uvDelta, uwDelta, vwDelta, currDelta;
    unsigned int nuv, nuw, nvw;
    
    //distribute outer for loops
#pragma omp for schedule(dynamic,1)
    for(size_t i=0;i<size;i++)
      {
	const Vert u = vert_vec[i];
	const v_size_t ui = index[u];
	
// #pragma omp single nowait
// 	{
// 	  cerr<<"Vertex: "<<u<<"\n";
// 	}

	uPredecessors = G[u].predecessors;
	//cerr<<"Vertex: "<<u<<"\n";

	for(size_t j=i+1;j<size;++j)
	  {
	    const Vert v = vert_vec[j];
	    const v_size_t vi = index[v];
	    const v_size_t uvd = d[ui][vi];
	  
	    vPredecessors = G[v].predecessors;
	    form_path(u,v,G,uPredecessors, uv, Q, on_path, index);

	    // if(j==(size-1))
	    //   {
	    // 	for(size_t k=0;k<size;k++)
	    // 	  {
	    // 	    puv[k].clear();
	    // 	    puw[k].clear();
	    // 	    pvw[k].clear();
		    
	    // 	    uv_meta[k][1] = 0;
	    // 	    uw_meta[k][1] = 0;
	    // 	    vw_meta[k][1] = 0;
	    // 	  }
	    //   }
	    
	    for(size_t k=j+1;k<size;++k)
	      {
		const Vert w = vert_vec[k];
		const v_size_t wi = index[w];

		const v_size_t uwd = d[ui][wi];
		const v_size_t vwd = d[vi][wi];

		v_size_t maxSide = uvd;

		if(uwd>maxSide)
		  maxSide = uwd;

		if(vwd>maxSide)
		  maxSide = vwd;
		
		  form_path(u,w,G,uPredecessors,uw,Q,on_path,index);
		  form_path(v,w,G,vPredecessors,vw,Q,on_path,index);
		
		//**Code commented out for BFS (less memory, MUCH slower version)**
		//BFS_two_target(u,v,w,G,uv,uw,predecessors);
	      
		// vector< vector<unsigned int> > puv;
		// vector< vector<unsigned int> > puw;
		// vector< vector<unsigned int> > pvw;


		  //#pragma omp parallel sections num_threads(3)
		    //{
		  //#pragma omp section
		  nuv = index_paths(u,v,index,size,uPredecessors,puv,node_stack);
		  //#pragma omp section		
		  nuw = index_paths(u,w,index,size,uPredecessors,puw,node_stack);
		  //#pragma omp section		
		  nvw = index_paths(v,w,index,size,vPredecessors,pvw,node_stack);	      
		  //}
		//**Debugging commented out**
		// cout<<"u: "<<u<<" v: "<<v<<"\n";
		// for(int l=0;l<uv.size();l++)
		// 	{
		// 	  Vert c = uv[l];
		// 	  v_size_t ci = index[c];
		// 	  cout<<c<<"\nOn Paths: ";
		// 	  for(int z=0;z<puv[ci].size();z++)
		// 	    {
		// 	      cout<<puv[ci][z]<<" ";
		// 	    }
		// 	  cout<<"\n";
		// 	}
		// char a;
		// cin>>a;	  
  
		//**Again, BFS code**
		//BFS_source_target(v,w,G,vw,predecessors);


		//Calculates delta for each side, ie returns the maximum minimum distance for 
		//any vertex on a shortest path between the pair 

		  //#pragma omp parallel sections shared(uvDelta,uwDelta,vwDelta) num_threads(3)
		    //{
		  //#pragma omp section
		  uvDelta = pair_delta(d,uv,uw,vw,puw,pvw,nuw,nvw,index,delta_v,k1,k2,size);
		  //#pragma omp section
		  uwDelta = pair_delta(d,uw,uv,vw,puv,pvw,nuv,nvw,index,delta_v,k1,k2,size);
		  //#pragma omp section
		  vwDelta = pair_delta(d,vw,uv,uw,puv,puw,nuv,nuw,index,delta_v,k1,k2,size);
		  //}
		//The worst of all pairs is the delta of that triplet
		currDelta = uvDelta;
		if(uwDelta>currDelta)
		  currDelta = uwDelta;
		if(vwDelta>currDelta)
		  currDelta = vwDelta;

		++delta_loc[maxSide-1][currDelta+1];
		if(currDelta>maxDelta_loc)
		  maxDelta_loc = currDelta;

		++delta_loc[maxSide-1][0]; 
	      }
	  }
      }

    #pragma omp critical
    {
      //cout<<delta_loc.size()<<"\n";
      for(size_t i=0;i<delta_loc.size();++i)
	for(size_t j=0;j<delta_loc[i].size();++j)
	  {
	    //cout<<i<<" "<<j<<"\n";
	    delta[i][j] += delta_loc[i][j];
	  }

      if(maxDelta_loc > maxDelta)
	maxDelta = maxDelta_loc;      
    }

  }
  for(int i=0;i<diam;++i)
    {
      delta[i].erase(delta[i].begin()+maxDelta+2,delta[i].end());
    }

  return(delta);
}




vector< vector<double> > calc_delta_fat(Graph& G, v_size_t diam)
{
  
  //Variable initializations
  v_size_t size = num_vertices(G);
  v_size_t biggest_delta = 3*diam;

  const vector<double> zero_dist(biggest_delta+2,0);
  double maxDelta = 0;


  // for(tie(vit,vite)=vertices(G); vit!=vite; vit++)
  //   {
  //     Vert v =*vit;
  //     vert_vec.push_back(v);
  //   }



  vector< vector<double> > delta(diam,zero_dist);

  if(num_vertices(G) < 3)
    {
      cerr<<"Graph is too small for hyperbolicity, returning zero\n";
   
      return delta;
    }

  Vert v = vertex(1,G);

  vector<bool> found(num_vertices(G),false);
  v_size_t current_comp = 0;
  v_size_t current_comp_size = 0;
  vector<v_size_t> component(num_vertices(G));
  
  BFS_vertices_found(v,G,found,current_comp,component,current_comp_size);

  vector<bool>::iterator it;
  for(it = found.begin();it<found.end();++it)
    {
      bool v_found = *it;
      if(!v_found)
	{
	  cerr<<"Graph passed to calc_gromov has multiple components\n";
	  G = connected(G);
	}
    }


#pragma omp parallel shared(delta,maxDelta,diam,size)
  {

    //Variable assignements,
    //Notation: puv = "Path indices on side uv"
    //           uv = "Vertices on side uv"
    //vertex iterators
    graph_traits<Graph>::vertex_iterator vit, vite;
  
    //Property maps
    property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);

    vector<Vert> vert_vec;
    vector< vector<v_size_t> > d;

    for(tie(vit,vite) = vertices(G); vit!= vite; ++vit)
      {
	Vert v = *vit;
	vert_vec.push_back(v);
	d.push_back(G[v].distances);
      }

    //    const size_t size = vert_vec.size();

    vector<unsigned int> dummy;
    
    vector<Vert> uv,uw,vw;
    vector<bool> on_path;
    vector<Vert> Q;

    Q.reserve(2*diam);
    on_path.reserve(size);
    uv.reserve(diam+1);
    uw.reserve(diam+1);
    vw.reserve(diam+1);

    vector<vector<Vert> > vPredecessors;
    vector<vector<Vert> > uPredecessors;

    vector< vector<double> > delta_loc(diam,zero_dist);
    double maxDelta_loc = 0;
    
    double currDelta = size;
    size_t uvi,uwi,vwi;
    
    Vert x,y,z;
    v_size_t xi,yi,zi;
    v_size_t xyd,yzd,xzd;
    v_size_t sum;

    //distribute outer for loops
#pragma omp for schedule(dynamic,1)
    for(size_t i=0;i<size;++i)
      {
	const Vert u = vert_vec[i];
	const v_size_t ui = index[u];
	
// #pragma omp single nowait
// 	{
// 	  cerr<<"Vertex: "<<u<<"\n";
// 	}

	uPredecessors = G[u].predecessors;
	//cerr<<"Vertex: "<<u<<"\n";

	for(size_t j=i+1;j<size;++j)
	  {
	    const Vert v = vert_vec[j];
	    const v_size_t vi = index[v];
	    const v_size_t uvd = d[ui][vi];
	  
	    vPredecessors = G[v].predecessors;
	    form_path(u,v,G,uPredecessors, uv, Q, on_path, index);

	    // if(j==(size-1))
	    //   {
	    // 	for(size_t k=0;k<size;k++)
	    // 	  {
	    // 	    puv[k].clear();
	    // 	    puw[k].clear();
	    // 	    pvw[k].clear();
		    
	    // 	    uv_meta[k][1] = 0;
	    // 	    uw_meta[k][1] = 0;
	    // 	    vw_meta[k][1] = 0;
	    // 	  }
	    //   }
	    
	    for(size_t k=j+1;k<size;++k)
	      {
		const Vert w = vert_vec[k];
		const v_size_t wi = index[w];

		const v_size_t uwd = d[ui][wi];
		const v_size_t vwd = d[vi][wi];

		v_size_t maxSide = uvd;

		if(uwd>maxSide)
		  maxSide = uwd;

		if(vwd>maxSide)
		  maxSide = vwd;
		
		  form_path(u,w,G,uPredecessors,uw,Q,on_path,index);
		  form_path(v,w,G,vPredecessors,vw,Q,on_path,index);

		  //Finding the worst delta_fat triangle in each triplet
		  currDelta = size;
		  for(uvi=0;uvi<uv.size();++uvi)
		    {
		      x  = uv[uvi];
		      xi = index[x];
		      for(uwi=0;uwi<uw.size();++uwi)
			{
			  y   = uw[uwi];
			  yi  = index[y];
			  xyd = d[xi][yi];
			  for(vwi=0;vwi<vw.size();++vwi)
			    {
			      z   = vw[vwi];
			      zi  = index[z];
			      yzd = d[yi][zi];
			      xzd = d[xi][zi];
			      sum = xyd+yzd+xzd;

			      //Is there a way to confirm the load balancing issue on this?  For example, 
			      //I could set the predecessors to a fixed set of nodes and see what happens each time?
			      //I could give everyone a set of the predecessor matrix and see how it scales in that case.

			      if(currDelta > sum)
				currDelta = sum;
			      
			    }
			}
		    }

		++delta_loc[maxSide-1][currDelta+1];
		if(currDelta>maxDelta_loc)
		  maxDelta_loc = currDelta;

		++delta_loc[maxSide-1][0]; 
	      }
	  }
      }

    #pragma omp critical
    {
      //cout<<delta_loc.size()<<"\n";
      for(size_t i=0;i<delta_loc.size();++i)
	for(size_t j=0;j<delta_loc[i].size();++j)
	  {
	    //cout<<i<<" "<<j<<"\n";
	    delta[i][j] += delta_loc[i][j];
	  }

      if(maxDelta_loc > maxDelta)
	maxDelta = maxDelta_loc;      
    }

  }
  for(size_t i=0;i<diam;++i)
    {
      delta[i].erase(delta[i].begin()+maxDelta+2,delta[i].end());
    }

  return(delta);
}




double pair_delta(vector<vector<v_size_t> >& d, const vector<Vert>& side, const vector<Vert>& other_side1, const vector<Vert>& other_side2, const vector<vector<unsigned int> >& path_ind1, const vector<vector<unsigned int> >& path_ind2, const unsigned int& num_paths1, const unsigned int& num_paths2, const property_map<Graph,vertex_index_t>::type& index, vector<unsigned int>& delta_v, vector<unsigned int>& k1, vector<unsigned int> k2, const size_t& n)
{
  //const v_size_t n = num_vertices(G);
  Vert x,y;
  v_size_t xi,yi;
  unsigned int dxy;
  unsigned long max1,max2;
  size_t i,j,k;
  unsigned path_max=num_paths2;
  const size_t other_side1_size = other_side1.size();
  const size_t other_side2_size = other_side2.size();
  const size_t side_size = side.size();
  unsigned int delta = 0;

  //delta_v.clear();
  k1.clear();
  k2.clear();

  //delta_v.assign(side.size(),0);
  k1.assign(num_paths1,n);
  k2.assign(num_paths2,n);


  if(num_paths1 > num_paths2)
    path_max = num_paths1;

  //Outer loop running through all points on all shortest paths (given by side)
  for(i=0;i<side_size;++i)
    {
      x = side[i];
      xi = index[x];
      
      //const v_size_t xi = index[x];
      //const vector<unsigned long> distance = G[x].distances;
      //k1.assign(num_paths1,n);o
      //k2.assign(num_paths2,n);
      // for(j=0;j<path_max;j++)
      // 	{
      // 	  if(j<num_paths1)
      // 	    k1[j] = n;
	  
      // 	  if(j<num_paths2)
      // 	    k2[j] = n;
	    
      // 	}
	
      //vector<unsigned int> k1(num_paths1, n);
      //vector<unsigned int> k2(num_paths2, n);

      //Look at each point on other sides (indexed by path using path_ind* for each vertex) 
      for(j=0;j<other_side1_size;++j)
	{
	  y = other_side1[j];
	  yi = index[y];
	  //vector<unsigned int> y_paths = path_ind1[yi];

	  dxy  = d[xi][yi];
	  
	  for(k=0;k<path_ind1[yi].size();++k)
	    {
	      if(k1[(path_ind1[yi][k])] > dxy)
		k1[(path_ind1[yi][k])] = dxy;
	    }
	  
	}
      
      for(j=0;j<other_side2_size;++j)
	{
	  y      = other_side2[j];
	  yi = index[y];
	  //vector<unsigned int> y_paths = path_ind2[yi];

	  dxy  = d[xi][yi];
	  
	  for(k=0;k<path_ind2[yi].size();++k)
	    {
	      if(k2[(path_ind2[yi][k])] > dxy)
		k2[(path_ind2[yi][k])] = dxy;
	    }
	  
	}

      max1 = 0;
      max2 = 0;
      for(j=0;j<path_max;++j)
	{
	  if(j<num_paths1)
	    {
	      if(k1[j]>max1)
		max1 = k1[j];

	      k1[j] = n;
	    }
	  
	  if(j<num_paths2)
	    {
	      if(k2[j]>max2)
		max2 = k2[j];

	      k2[j] = n;
	    }
   
	}

       if(max1>max2)
	 {
	   //delta_v[i] = max2;
	   if(max2>delta)
	     delta = max2;
	 }
       else
	 {
	   //delta_v[i] = max1;
	   if(max1>delta)
	     delta = max1;
	 }

    }

  // for(i=0;i<delta_v.size();i++)
  //   {
  //     if(delta_v[i]>delta)
  // 	delta = delta_v[i];
  //   }
  
  return(delta);
}




unsigned int index_paths(const Vert& s, const Vert& t, const property_map<Graph,vertex_index_t>::type& index, const v_size_t& n, const vector< vector<Vert> >& p, vector< vector<unsigned int> >& ind, vector<path_node>& node_stack)
{
  //property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);
  //v_size_t n = num_vertices(G);

  size_t i;

  if(ind.empty())
    {
      const vector<unsigned int> dummy;
      ind.assign(n,dummy);
    }
  else if(ind.size()<n)
    {
      const vector<unsigned int> dummy;
      for(size_t i=0;i<ind.size();++i)
  	ind[i].clear();

      for(size_t i=0;i<n-ind.size();++i)
  	ind.push_back(dummy);
    }
  else 
    {
      for(i=0;i<n;++i)
	ind[i].clear();
    }
  
  // for(size_t i=0;i<n;i++)
  //   {
  //     ind_meta[i][0]=0;
  //   }

  
  unsigned int path = 0;

  node_stack.clear();
  //vector node_stack is being used as a stack, with the exception that I will iterate through the nodes whenever a path is found
  //vector<path_node> node_stack;
  //  node_stack.reserve(CHUNKSIZE);

  //  v_size_t si = index[s];
  const v_size_t ti = index[t];

  path_node init;
  path_node top;
  path_node next;

  Vert pv;
  v_size_t pvi;
  v_size_t vi;

  init.v = t;
  init.predecessors = ti;
  init.count = 0;

  node_stack.push_back(init);

  while(!node_stack.empty())
    {

      top = node_stack.back();
      node_stack.pop_back();

      if(top.count < p[top.predecessors].size())
	{
	  if(top.v == s)
	    {
	      node_stack.push_back(top);
	      for(i=0;i<node_stack.size(); ++i)
		{
		  pv  = node_stack[i].v;
		  pvi = index[pv];
	      
		  ind[pvi].push_back(path);

		
		}
	      path++;
	      node_stack.pop_back();
	      
	      //cout<<"Source pop!\n";
	    }
	  else
	    {
	      next.v = p[top.predecessors][top.count];
	      vi = index[next.v];
	      ++top.count;
	      next.predecessors = vi;
	      next.count = 0;
	      node_stack.push_back(top);
	      node_stack.push_back(next);
	      //cout<<"Push!\n";
	    }
	}

	
    }
  return(path);
}

void form_path(const Vert& s, const Vert& t, const Graph& G, const vector< vector<Vert> >& predecessors, vector<Vert>& path, vector<Vert>& Q, vector<bool>& onPath,property_map<Graph,vertex_index_t>::type& index)
{

  // property_map<Graph,vertex_index_t>::type index = get(vertex_index,G);
  const v_size_t n  = num_vertices(G);
  const v_size_t si = index[s];
  size_t q_front = 0;
  size_t q_end   = 0;
  
  Q.clear();

  onPath.clear();
  path.clear();

  onPath.assign(n,false);
  
  //Q push
  Q.push_back(t);
  ++q_end;
  
  if(n!=predecessors.size())
    cerr<<"Warning, in form_path: Size of predecessor vector does not match number of vertices in G\n";

  if(predecessors[si].size()!=1 || predecessors[si][0]!=s)
    cerr<<"Warning, in form_path: The predecessor of source s is not s, possible invalid predecessors vector\n";

  while(q_end!=q_front)
    {
      const Vert current = Q[q_front];
      const v_size_t ci  = index[current];

      const vector<Vert> preds = predecessors[ci];
      
      //Q pop()
      ++q_front;

      path.push_back(current);

      if(current!=s)
	{
	  for(unsigned int i=0;i<preds.size();++i)
	    {
	      const Vert p      = preds[i];
	      const v_size_t pi = index[p];
	  
	      if(!onPath[pi])
		{
		  Q.push_back(p);
		  ++q_end;

		  onPath[pi] = true;
		}
	    }
	}
    }
}
