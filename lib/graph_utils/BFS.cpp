/*
This function is intended to perform a breadth first search of a graph, beginning at source node s, storing distances and predecessors (all predecessors), as it reaches each node in the graph.  It is intended to be used with ../graph_lib_boost.hpp.  It is used instead of the Boost breadth-first-search because it took less time to write a breadth-first-search with less functionality than to understand and modify the Boost algorithm.

By Aaron Adcock, PhD Candidate, Stanford University
Jan 2012
*/

#include "../graph_lib_boost.hpp"

//Takes the predecessors for a given source s, with target t,
//and returns all vertices on shortest paths between them.
void form_path(const Vert& s, const Vert& t, const Graph& G, const vector< vector<Vert> >& predecessors, vector<Vert>& path, queue<Vert>& Q, vector<bool>& onPath,property_map<Graph,vertex_index_t>::type& index);

// Finds all vertices in a given connected component of a graph G,
// returns vertices in component and component size in comp_size,
// vertices searched in searched
void BFS_vertices_found(const Vert& s, const Graph& G, vector<bool>& searched, const v_size_t & current_comp, vector<v_size_t> & component, v_size_t & comp_size)
{
  //Data structure assignment and basic graph statistics (size, vertex index map)
  comp_size = 0;

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  //searched.assign(n,false);
  queue<Vert> unsearched;
  v_size_t si = index[s];

  searched[si] = true;
  ++comp_size;
  component[si] = current_comp;

  //seed node for queue (ie where search begins)
  unsearched.push(s);

  while (!unsearched.empty())
    {
      //next node in queue
      const Vert current = unsearched.front();
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate over adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie;++ai)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];

	  //check if already visited, if not mark and add to queue
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      ++comp_size;
	      component[ni] = current_comp;
	      
	      unsearched.push(next);
	    }

	}
    }
}


// Finds distances and shortest paths (where they exist) to all
// vertices in G beginning at vertex s Returns distances in vector
// distances (distance from s to vertex with index i is distances[i])
// Returns shortest paths in predecessors (predecessors[i] = vector of
// vertices adjacent to i that were found immediately before i, ie are
// predecessors to i on any shortest path from s to i.
void BFS_source_all(const Vert& s, const Graph& G, vector<v_size_t>& distances, vector< vector<Vert> >& predecessors)
{

  //Data structures set up
  distances.clear();
  predecessors.clear();

  v_size_t n = num_vertices(G);
  vector<v_size_t> dist(n,0);
  const vector<Vert> dummy;

  for (unsigned int i=0;i<n;i++)
    predecessors.push_back(dummy);
  
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  vector<bool> searched(n,false);
  queue<Vert> unsearched;

  //source (seed) node
  searched[(index[s])] = true;
  
  unsearched.push(s);

  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate through adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie;++ai)
	{
	  Vert next  = *ai;
	  v_size_t ni = index[next];
	  
	  //if adjacent vertex unsearched, mark, update distance, add
	  //to queue, if searched, check to see if it has the current
	  //distance of interest and then add to predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      unsearched.push(next);
	    }
	  else if(dist[ni]==(dist[ci]-1))
	    {
	      predecessors[ci].push_back(next);
	    }
	}
    }
  predecessors[(index[s])].push_back(s);
  distances = dist;
}

// Finds distances and shortest paths (where they exist) to all
// vertices in G beginning at vertex s.  This version of the function
// returns the vector bfs_tree, a set of bfs_node's sorted by level, a
// vector of labels (maps graph node index to tree node index), and
// returns the height of the tree.  A more standard return.
void BFS_source_all(const Vert& s, const Graph& G, vector<bfs_node> & bfs_tree, vector<size_t> & tree_label, size_t & height)
{
  v_size_t n = num_vertices(G);
  bfs_node root_node(0, 0);

  bfs_tree.clear();
  bfs_tree.assign(n, root_node); 
  tree_label.clear();
  tree_label.assign(n, 0);
  size_t node_to_label = 0;
  height = 0;
  
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  vector<int> searched(n, 0);
  queue<Vert> unsearched;

  // As we add nodes to the queue, they need to be added to bfs_tree
  // in bfs_node form and to the tree_label.
  searched[(index[s])] = 1;
  unsearched.push(s);
  bfs_tree[node_to_label] = root_node;
  tree_label[index[s]] = node_to_label++;

  while (!unsearched.empty())
    {
      // Vist next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];

      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      // Iterate through adjacent vertices
      for (tie(ai, aie) = adjacent_vertices(current, G); ai != aie; ++ai)
	{
	  const Vert next = *ai;
	  const v_size_t ni = index[next];

	  if (searched[ni] == 0)
	    {
	      searched[ni] = 1;
	      unsearched.push(next);
	      
	      // Add node to bfs_tree and tree_label, keep track of
	      // max level.
	      bfs_node next_node;
	      next_node.parent = tree_label[ci];
	      next_node.level = bfs_tree[next_node.parent].level + 1;
	      bfs_tree[node_to_label] = next_node;
	      tree_label[ni] = node_to_label++;
	      if (next_node.level > height)
		height = next_node.level;
	    }
	}
    }
}


void BFS_source_all(const Vert& s, const Graph& G, vector<v_size_t>& distances)
{

  //Data structures set up
  distances.clear();

  v_size_t n = num_vertices(G);
  vector<v_size_t> dist(n,0);
  const vector<Vert> dummy;

  
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  vector<bool> searched(n,false);
  queue<Vert> unsearched;

  //source (seed) node
  searched[(index[s])] = true;
  
  unsearched.push(s);

  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate through adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie;++ai)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];
	  
	  //if adjacent vertex unsearched, mark, update distance, add
	  //to queue, if searched, check to see if it has the current
	  //distance of interest and then add to predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      unsearched.push(next);
	    }
	}
    }
  distances = dist;
}


//Finds distances between s and all other vertices, also keeps track
//of number of shortest paths between the two.
void BFS_source_all(const Vert& s, const Graph& G, vector<v_size_t>& distances, vector<v_size_t>& num_paths)
{

  /* Data structures set up */
  
  // Clear structures to store distance to each vertex and the number
  // of paths to each vertex
  distances.clear();
  num_paths.clear();

  v_size_t n = num_vertices(G);
  // Temporary structure
  vector<v_size_t> dist(n,0);

  // Probably superfluous
  num_paths.assign(n,0);
  
  // Get index map of vertices
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  vector<bool> searched(n,false);
  queue<Vert> unsearched;

  // Source (seed) node is searched
  searched[(index[s])] = true;

  // Set to one, but then set to zero at end of function
  num_paths[(index[s])] = 1;

  unsearched.push(s);

  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate through adjacent vertices
      for (tie(ai,aie) = adjacent_vertices(current,G); ai != aie; ++ai)
	{
	  const Vert next   = *ai;
	  const v_size_t ni = index[next];
	  
	  //if adjacent vertex unsearched, mark, update distance, add
	  //to queue, if searched, check to see if it has the current
	  //distance of interest and then add to predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      unsearched.push(next);
	    }
	  else if(dist[ni] == (dist[ci]-1))
	    {
	      num_paths[ci] += num_paths[ni];
	    }
	}
    }
  num_paths[(index[s])] = 0;
  distances = dist;
}


//Finds shortest path/distance between source s and target t in
//G. Returns a shortest path and all predecessor nodes (which form all
//shortest paths).
void BFS_source_target(const Vert& s, const Vert& t, Graph& G, vector<Vert>& path, vector< vector<Vert> >& predecessors)
{
  //Set up data structures
  path.clear();
  predecessors.clear();

  v_size_t n = num_vertices(G);

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  
  vector<bool> searched(n,false);
  queue<Vert> unsearched;
  vector<v_size_t> dist(n,0);
  vector<Vert> dummy;

 
  //Source node
  for (unsigned int i=0;i<n;i++)
    predecessors.push_back(dummy);

  searched[(index[s])] = true;
  
  unsearched.push(s);

  bool notFound = true;
  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;
      
      //iterate over adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie; ai++)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];

	  //if not searched, mark, add to queue, check if it is target
	  //node, end while loop if t found Else, update predecessor
	  //vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      
	      if(notFound)
		unsearched.push(next);
	      
	      if(next == t)
		notFound = false;
		

	    }
	  else if(dist[ni] == (dist[ci]-1))
	    {
	      predecessors[ci].push_back(next);
	    }
	}
    }
  predecessors[(index[s])].push_back(s);

  //Form path from predecessor node
  queue<Vert> Q;
  vector<bool> on_path;

  form_path(s,t,G,predecessors,path,Q,on_path,index);
  // const v_size_t tind = index[t];
  // graph_traits<Graph>::adjacency_iterator ai, aie;
  // for (tie(ai,aie)=adjacent_vertices(t,G); ai!=aie; ai++)
  //   {
  //     const Vert a = *ai;
  //     const v_size_t aind = index[a];
      
  //     if(searched[aind] && dist[aind]==(dist[tind]-1))
  // 	predecessors[tind].push_back(a);
  //   }
  

}


//Same as source target, but returns only one geodesic.
void BFS_source_target_single(const Vert& s, const Vert& t, Graph& G, vector<Vert>& path, vector< vector<Vert> >& predecessors)
{
  //Set up data structures
  path.clear();
  predecessors.clear();

  v_size_t n = num_vertices(G);

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  
  vector<bool> searched(n,false);
  queue<Vert> unsearched;
  vector<v_size_t> dist(n,0);
  vector<Vert> dummy(1,0);

 
  //Source node
  for (unsigned int i=0;i<n;i++)
    predecessors.push_back(dummy);

  searched[(index[s])] = true;
  
  unsearched.push(s);

  bool notFound = true;
  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;
      
      //iterate over adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie; ai++)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];

	  // if not searched, mark, add to queue, check if it is
	  // target node, end while loop if t found Else, update
	  // predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      
	      if(notFound)
		unsearched.push(next);
	      
	      if(next == t)
		notFound = false;
		

	    }
	  else if(dist[ni] == (dist[ci]-1))
	    {
	      predecessors[ci][0] = next;
	    }
	}
    }
  predecessors[(index[s])][0] = s;

  //Form path from predecessor node
  queue<Vert> Q;
  vector<bool> on_path;

  form_path(s, t, G, predecessors, path, Q, on_path, index);
  // const v_size_t tind = index[t];
  // graph_traits<Graph>::adjacency_iterator ai, aie;
  // for (tie(ai,aie)=adjacent_vertices(t,G); ai!=aie; ai++)
  //   {
  //     const Vert a = *ai;
  //     const v_size_t aind = index[a];
      
  //     if(searched[aind] && dist[aind]==(dist[tind]-1))
  // 	predecessors[tind].push_back(a);
  //   }
  

}
    
//Find shortest path to two target nodes.  Essentially the same as
//above algorithm, except it finds two paths.  Created to save time in
//finding geodesics in graph triangles
void BFS_two_target(const Vert& s, const Vert& t1, const Vert& t2, Graph& G, vector<Vert>& path1, vector<Vert>& path2, vector< vector<Vert> >& predecessors)
{
  //set up data structures
  path1.clear();
  path2.clear();
  predecessors.clear();
  
  //path1.push_back(s);
  //  path2.push_back(s);

  v_size_t n = num_vertices(G);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  
  vector<bool> searched(n,false);
  queue<Vert> unsearched;
  vector<v_size_t> dist(n,0);
  const vector<Vert> dummy;

  for (unsigned int i=0;i<n;i++)
    predecessors.push_back(dummy);

  //Source vertex
  searched[(index[s])] = true;
  
  unsearched.push(s);

  bool notFound1 = true;
  bool notFound2 = true;
  
  while (!unsearched.empty())
    {
      //next unsearched vertex
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate over adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie; ai++)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];

	  //Check if next vertex already searched, update
	  //distance/mark/add to queue if so If not, update
	  //predecessor vectors
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;

	      if(notFound1 || notFound2)
		unsearched.push(next);

	      if(next == t1)
		notFound1 = false;
		
	      if(next == t2)		
		notFound2 = false;
		
	    }
	  else if(dist[ni] == (dist[ci]-1))
	    {
	      predecessors[ci].push_back(next);
	    }
	}
    }
  
  predecessors[(index[s])].push_back(s);

  //form 2 paths from predecessors/G
  queue<Vert> Q;
  vector<bool> on_path;

  form_path(s,t1,G,predecessors,path1,Q,on_path,index);
  form_path(s,t2,G,predecessors,path2,Q,on_path,index);

}

// Find shortest path to three target nodes.  Essentially the same as
// above algorithm, except it finds two paths.  Created to save time
// in finding geodesics in graph triangles
void BFS_three_target(const Vert& s, const Vert& t1, const Vert& t2, const Vert& t3, Graph& G, vector<Vert>& path1, vector<Vert>& path2, vector<Vert>& path3, vector< vector<Vert> >& predecessors)
{
  //set up data structures
  path1.clear();
  path2.clear();
  path3.clear();
  predecessors.clear();
  
  //path1.push_back(s);
  //  path2.push_back(s);

  v_size_t n = num_vertices(G);
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  
  vector<bool> searched(n,false);
  queue<Vert> unsearched;
  vector<v_size_t> dist(n,0);
  const vector<Vert> dummy;

  for (unsigned int i=0;i<n;i++)
    predecessors.push_back(dummy);

  //Source vertex
  searched[(index[s])] = true;
  
  unsearched.push(s);

  bool notFound1 = true;
  bool notFound2 = true;
  bool notFound3 = true;
  
  while (!unsearched.empty())
    {
      //next unsearched vertex
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate over adjacent vertices
      for (tie(ai,aie)=adjacent_vertices(current,G); ai!=aie; ai++)
	{
	  const Vert next  = *ai;
	  const v_size_t ni = index[next];

	  //Check if next vertex already searched, update
	  //distance/mark/add to queue if so If not, update
	  //predecessor vectors
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;

	      if(notFound1 || notFound2 || notFound3)
		unsearched.push(next);

	      if(next == t1)
		notFound1 = false;
		
	      if(next == t2)		
		notFound2 = false;

	      if(next == t3)
		notFound3 = false;
		
	    }
	  else if(dist[ni] == (dist[ci]-1))
	    {
	      predecessors[ci].push_back(next);
	    }
	}
    }
  
  predecessors[(index[s])].push_back(s);

  //form 2 paths from predecessors/G
  queue<Vert> Q;
  vector<bool> on_path;

  form_path(s,t1,G,predecessors,path1,Q,on_path,index);
  form_path(s,t2,G,predecessors,path2,Q,on_path,index);
  form_path(s,t3,G,predecessors,path3,Q,on_path,index);

}

// In the target BFS functions, not all vertices in predecessors are
// predecessors on the path of interest (we essentially have the
// vertices for a ball around s, form_path finds the vertices on ALL
// shortest paths (s, t) and returns them
void form_path(const Vert& s, const Vert& t, const Graph& G, const vector< vector<Vert> >& predecessors, vector<Vert>& path, queue<Vert>& Q, vector<bool>& onPath,property_map<Graph,vertex_index_t>::type& index)
{

  // property_map<Graph,vertex_index_t>::type index =
  // get(vertex_index,G); Sets up data structures, most passed in to
  // increase efficiency Note that predecessors contains all
  // predecessors to each vertex when graph is searched starting at s.
  const v_size_t n  = num_vertices(G);
  const v_size_t si = index[s];
  
  if(!Q.empty())
    {
      cerr<<"Warning, in form_path: Q passed to form_path that is not empty\n";
      while (!Q.empty())
  	Q.pop();
    }

  onPath.clear();
  path.clear();

  onPath.assign(n,false);
  Q.push(t);
  
  if(n!=predecessors.size())
    cerr<<"Warning, in form_path: Size of predecessor vector does not match number of vertices in G\n";

  if(predecessors[si].size()!=1 || predecessors[si][0]!=s)
    cerr<<"Warning, in form_path: The predecessor of source s is not s, possible invalid predecessors vector\n";

  //move back through the graph from target vertex t
  while (!Q.empty())
    {
      const Vert current = Q.front();
      const v_size_t ci  = index[current];

      const vector<Vert> preds = predecessors[ci];
      
      Q.pop();
      path.push_back(current);

      //Add to path/queue all vertices that are predecessors to
      //current (which at some point were predecessors to t)
      //Terminates when queue is empty, ie all paths have reached back
      //to s.
      if(current != s)
	{
	  for (unsigned int i = 0; i < preds.size(); i++)
	    {
	      const Vert p      = preds[i];
	      const v_size_t pi = index[p];
	  
	      if(!onPath[pi])
		{
		  Q.push(p);
		  onPath[pi] = true;
		}
	    }
	}
    }
}

int BFS_eccentricity(const Vert& s, const Graph& G)
{

  /* Data structures set up */
  v_size_t n = num_vertices(G);
  // Temporary structure
  vector<v_size_t> dist(n,0);
  
  // Get index map of vertices
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);
  vector<bool> searched(n,false);
  queue<Vert> unsearched;

  size_t max_dist = 0;

  // Source (seed) node is searched
  searched[(index[s])] = true;

  unsearched.push(s);

  while (!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate through adjacent vertices
      for (tie(ai,aie) = adjacent_vertices(current,G); ai != aie; ++ai)
	{
	  const Vert next   = *ai;
	  const v_size_t ni = index[next];
	  
	  //if adjacent vertex unsearched, mark, update distance, add
	  //to queue, if searched, check to see if it has the current
	  //distance of interest and then add to predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      unsearched.push(next);
	      if (dist[ni] > max_dist)
		max_dist = dist[ni];
	    }
	}
    }
  return max_dist;
}
