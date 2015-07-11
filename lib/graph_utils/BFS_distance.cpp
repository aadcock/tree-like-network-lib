/*
This file contains a function which takes a graph G, calculates all of
the pairwise distances (and eventually, keeps all of the shortest
paths) of G.  It is implemented using a breadth-first search (BFS).
For weighted graphs, Djikstra's algorithm should be used.

By Aaron Adcock, PhD candidate at Stanford University, Dec. 2011
 */

#include "../graph_lib_boost.hpp"
#include <fstream>
#include <queue>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>

void find_furthest_vertex(Graph & T, Vert source, Vert& furthest_vertex, v_size_t & max_distance);

void BFS_distance(Graph& G, string outputFileName)
{

  graph_traits<Graph>::vertex_iterator vi, vie;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);


  //Matrices for storing short path predecessors/shortest distance
  //between vertices.
  std::vector< std::vector<vector<Vert> > > short_path;
  std::vector< std::vector<v_size_t> > short_dist;
  std::vector< v_size_t > labels;

  for(tie(vi,vie)=vertices(G);vi!=vie;vi++)
    {
      Vert s = *vi;  //s is the source for this BFS
      

      labels.push_back(index[s]);
      //p is for keeping the predecessors
      vector<vector<Vert> > p(num_vertices(G));
      //v_size_t d[num_vertices(G)];
      //fill_n(d,num_vertices(G),0);

      //p[s] = s;

      //The visitor concept lets actions be performed as each node is 
      //visited throughout the search tree.  In this case, I am using two 
      //BGL defined visitors record_distances and record_predecessors
      //which keep the distances and the path predecessors.  I need to modify
      //this to keep ALL shortest paths.
      // breadth_first_search(G,s,visitor(make_bfs_visitor(std::make_pair(record_distances(d,on_tree_edge()),record_predecessors(&p[0],on_tree_edge() ) ) ) ) );



      //d_v keeps the distances to other vertices
      //      vector<v_size_t> d_v(d,d+num_vertices(G));
      vector<v_size_t> d;
      BFS_source_all(s,G,d,p);
      G[s].distances.clear();
      G[s].distances.insert(G[s].distances.begin(),d.begin(),d.end());

      G[s].predecessors.clear();
      G[s].predecessors = p;

      //Aaron, you need to modify this to create the short matrices
      //based on vertex index, that way you can be sure what matrix
      //row corresponds to what vertex.

      short_path.push_back(p);
      short_dist.push_back(d);
    }

  //distance iterator and predecessor iterator
  std::vector< vector<v_size_t> >::iterator dit;
  std::vector< vector<vector<Vert> > >::iterator pit;

  int count = 0;

  ofstream outputDistance;
  string distanceFile = outputFileName;
  distanceFile.append("_distances.txt");
  outputDistance.open(distanceFile.c_str());
  
  ////////////////////////////////////////////////////
  //Write the distance file
  ////////////////////////////////////////////////////
  for(dit=short_dist.begin();dit<short_dist.end();dit++)
    {
      std::vector<v_size_t>::iterator sit;
      std::vector<v_size_t> dt = *dit;

      outputDistance<<labels[count]<<"\t";

      for(sit = dt.begin();sit<dt.end();sit++)
	{
	  v_size_t x = *sit;
	  outputDistance<<x<<"\t";
	}
      count++;
      outputDistance<<"\n";
    } 
  outputDistance.close();

  // distanceFile = outputFileName;
  // distanceFile.append("_distances.txt");
  // outputDistance.open(distanceFile);

  // counter = 0;
  // count = 0;
  // graph_traits<Graph>::vertex_iterator vit,vitend;
  // for(tie(vit,vitend)=vertices(G); vit!=vitend; ++vit)
  //   {
  //     Vert v = *vit;
  //     std::vector<unsigned long>::iterator sit;
  //     std::vector<unsigned long> dt = G[v].distances;
  //     for(sit = dt.begin();sit<dt.end();sit++)
  // 	{
  // 	  unsigned long x = *sit;
  // 	  outputDistance<<x<<"\t";
  // 	}
  //     count++;
  //     outputDistance<<"\n";
  //   }
  // outputDistance.close();

  string pathFile = outputFileName;
  pathFile.append("_short_paths.txt");

  outputDistance.open(pathFile.c_str());

  /////////////////////////////////////////////////////
  //Write the predecessor file
  /////////////////////////////////////////////////////
   count = 0;
  for(pit=short_path.begin();pit<short_path.end();pit++)
    {
      std::vector< vector<Vert> >::iterator sit;
      std::vector< vector<Vert> > dt = *pit;

      outputDistance<<labels[count]<<"\t";
      count++;

      for(sit = dt.begin();sit<dt.end();sit++)
	{
	  vector<Vert> x = *sit;
	  std::vector<Vert>::iterator xit;

	  outputDistance<<"(";
	  for(xit = x.begin();xit<x.end();xit++)
	    {
	      Vert w = *xit;
	      if(xit == x.end()-1)
		outputDistance<<index[w];
	      else
		outputDistance<<index[w]<<",";
	    }
	  outputDistance<<")\t";
	}
      outputDistance<<"\n";
    }
   
  outputDistance.close();
}



//This version does NOT print output files and does NOT keep shortest
//path info.  It just attaches the associated distance information to
//each node.  This is more memory efficient, especially for large
//graphs.  This function uses n^2 in memory.
void BFS_distance(Graph& G)
{
  graph_traits<Graph>::vertex_iterator vi, vie;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  //Matrices for storing short path predecessors/shortest distance
  //between vertices.  std::vector< std::vector<v_size_t> >
  //short_dist;
  std::vector< v_size_t > labels;

  for(tie(vi, vie) = vertices(G); vi != vie; vi++)
    {
      Vert s = *vi;  //s is the source for this BFS
      

      labels.push_back(index[s]);
      //p is for keeping the predecessors

      //v_size_t d[num_vertices(G)];
      //fill_n(d,num_vertices(G),0);

      //p[s] = s;

      //The visitor concept lets actions be performed as each node is 
      //visited throughout the search tree.  In this case, I am using two 
      //BGL defined visitors record_distances and record_predecessors
      //which keep the distances and the path predecessors.  I need to modify
      //this to keep ALL shortest paths.
      // breadth_first_search(G,s,visitor(make_bfs_visitor(std::make_pair(record_distances(d,on_tree_edge()),record_predecessors(&p[0],on_tree_edge() ) ) ) ) );



      //d_v keeps the distances to other vertices
      //      vector<v_size_t> d_v(d,d+num_vertices(G));
      vector<v_size_t> d;
      BFS_source_all(s,G,d);
      G[s].distances.clear();
      G[s].distances.insert(G[s].distances.begin(),d.begin(),d.end());

      //Aaron, you need to modify this to create the short matrices
      //based on vertex index, that way you can be sure what matrix
      //row corresponds to what vertex.

      //short_dist.push_back(d);
    }

}

// This function calculates the eccentricity of a tree.  It runs in
// O(n) time.  First, find the center/radius of the tree, (2n
// operations).  Then we calculate the distance between that node and
// all others.  The input is assumed to be a connected tree.  If it is not a
// connected tree, returns empty vector.
vector<v_size_t> tree_eccentricity(Graph& T)
{
  // Basic properties
  v_size_t n_v = num_vertices(T);
  e_size_t n_e = num_edges(T);

  // Test of tree
  if (n_v == (n_e + 1))
    {
      cout<<"Tree not passed to eccentricity algorithm\n";
      vector<size_t> empty_vector;
      return empty_vector;
    }

  // Pick a vertex
  graph_traits<Graph>::vertex_iterator vi, vie;
  tie(vi, vie) = vertices(T);
  Vert source = *vi;
  Vert diameter_vertex_1, diameter_vertex_2;
  v_size_t diameter = 0;

  // Have to find furthest vertex twice to get diameter/diameter vertices
  find_furthest_vertex(T, source, diameter_vertex_1,  diameter);
  find_furthest_vertex(T, diameter_vertex_1, diameter_vertex_2, diameter);

  // Get vertex path between diameter points
  vector<Vert> path;
  vector< vector<Vert> > temp;
  BFS_source_target_single(diameter_vertex_1, diameter_vertex_2, T, path, temp);

  // The path should contain vertices in order from diameter 2 to
  // diameter 1, starting with diameter 2 and excluding diameter 1.
  // Just need to grab the appropriate vertices, must figure out if
  // tree is bicentered

  bool bicentered = (diameter % 2) != 0;
  vector<Vert> central_nodes;

  if (bicentered)
    {
      central_nodes.push_back(path[diameter / 2]);
      central_nodes.push_back(path[diameter / 2 + 1]);
    }
  else
    central_nodes.push_back(path[diameter / 2]);

  // To calculate eccentricity, find distance to closest central node
  // and add that to floor(diameter / 2).  If current code is slow,
  // could write a custom BFS_source_all that updates distances as the
  // BFS occurs.
  vector <v_size_t> eccentricity(n_v, (v_size_t) - 1);
  vector < vector<v_size_t> > distances (central_nodes.size(), eccentricity);
  for (int i = 0; i < central_nodes.size(); ++i)
    BFS_source_all(central_nodes[i], T, distances[i]);

  for (int i = 0; i < n_v; ++i)
    for (int j = 0; j < central_nodes.size(); ++j)
      if (eccentricity[i] > distances[j][i] + diameter / 2)
	eccentricity[i] = distances[j][i] + diameter / 2;

  return eccentricity;
}

// Returns the vertex at furthest distance as well as the distance.
// Touches 2n vertices in a tree.  In a nontree, has to check the
// searched vector more often.
void find_furthest_vertex(Graph & T, Vert source, Vert& furthest_vertex, v_size_t & max_distance)
{
  /* Data structures set up */

  v_size_t n = num_vertices(T);
  // Temporary structure
  vector<v_size_t> dist(n, 0);
  
  // Get index map of vertices
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, T);
  vector<bool> searched(n, false);
  std::queue<Vert> unsearched;

  max_distance = 0;
  furthest_vertex = source; 

  // Source (seed) node is searched
  searched[(index[source])] = true;

  unsearched.push(source);

  while(!unsearched.empty())
    {
      //next unsearched node
      const Vert current = unsearched.front();
      const v_size_t ci  = index[current];
      
      unsearched.pop();
      graph_traits<Graph>::adjacency_iterator ai, aie;

      //iterate through adjacent vertices
      for(tie(ai,aie) = adjacent_vertices(current,T); ai != aie; ++ai)
	{
	  const Vert next   = *ai;
	  const v_size_t ni = index[next];
	  
	  // if adjacent vertex unsearched, mark, update distance, add
	  // to queue, if searched, check to see if it has the current
	  // distance of interest and then add to predecessor vector
	  if(!searched[ni])
	    {
	      searched[ni] = true;
	      dist[ni] = dist[ci]+1;
	      unsearched.push(next);
	      if (dist[ni] > max_distance)
		{
		  max_distance = dist[ni];
		  furthest_vertex = next;
		}
	    }
	}
    }
}
