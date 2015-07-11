
/* Approximates the local PageRank vector according to eps and alpha
using the procedure found in Anderson, Chung, Lang's "Local Graph
Partitioning using PageRank Vectors".  Jure Leskovec's SNAP code was
referenced as well.

eps is the approximation factor (which is used in the ncp calculations
as 1/approximate community size), alpha is the teleport probability, r
is the starting distribution (ie residual), p is the output
approximate page rank, ind is an index that eliminates the need to run
through the entire r vector at the beginning.  Probably should be
removed

The local PageRank vector is calculated by repeatedly pushing the
probability from each node to another.

Jan. 2014, I also added a function which runs a vanilla diffusion from
a provided initial mass vector.  ie If you use the indicator vector,
it shows how the mass spreads around the network.  It currently runs
for a set number of iterations.

This function uses the Boost Graph Library

by Aaron Adcock, PhD Candidate at Stanford University 
July 2011  */



#include "../graph_lib_boost.hpp"
#include <string>
#include <vector>
#include <queue>
#include <iostream>
#include <float.h>

using namespace boost;

// The function push is used to push probability from one node to the
// surrounding nodes.  alpha of the probability stays at the current
// node, the rest is split evenly among the neighbors.  Main workhorse
// of algorithm

void push(int ind, vector<double> & p, vector<double> & r, const double& alpha, Graph& G, const double& eps, queue<int> & Q);

// Main Function: Sets up queue of nodes based on how much probability
// is 'on' each node.  This queue is then updated by the push function
// and nodes are added or removed based on the degree normalized
// 'probability' at each node of r.

void approxPR(Graph& G, vector<double>& r, const double& alpha, const double& eps, vector<double>& p, const int ind)
{
  int size = r.size();
  // vector<double> p (size,0);
  queue<int> epsQueue;
   
  // s gives the amount of probability at each node
  // Add i to queue if it contains any probability

  if(ind < 0)
    {
      for (size_t i = 0; i < size; ++i)
	{
	  if (r[i] > 0.0)
	    epsQueue.push(i);
	}
    }
  else if (r[ind] >= 0)
    {
      epsQueue.push(ind);
    }

  // Run through nodes in queue, note that push adds additional nodes
  // to queue epsQueue.pop();
  while (!epsQueue.empty())
    {
      // Get node at front of queue
      int rindex = epsQueue.front();
      epsQueue.pop();
      // Call push to update p, r using node with index rindex, and epsQueue
      push(rindex, p, r, alpha, G, eps, epsQueue);

    }

 //return(p);

}

void push(int ind, vector<double> & p, vector<double> & r, const double& alpha, Graph& G, const double& eps, queue<int> & Q)
{

  // Creates index set for the vertices of G
  
  // Used to store push value
  double val;

  // Graph iterator
  graph_traits<Graph>::adjacency_iterator it, itend;  

  // Assign current vertex to v
  Vert v = vertex(ind, G);

  // Find number of outgoing edges of G (Note: this is assumed to be a
  // directed graph)
  graph_traits<Graph>::degree_size_type outDegree = out_degree(v, G);

  // Probability pushed from r vector to the approximation vector p
  p[ind] += alpha * (r[ind] - 0.5 * eps * outDegree);

  // Store value to be distributed among neighbors
  val = (1 - alpha) * (r[ind] - 0.5 * eps * outDegree) / outDegree;
 
  // Probability that remains on vertex ind in r vector
  r[ind] = 0.5 * eps * outDegree;
  
  // iterate adjacent vertices
  for (tie(it, itend) = adjacent_vertices(v, G); it != itend; ++it)
    { 
      Vert neighbor = *it;

      graph_traits<Graph>::degree_size_type neighbor_degree = out_degree(neighbor, G);

      const double rPrev = r[neighbor];
      
      // Store new value and add to Q if amount of probability in r
      // exceeds threshold, and neighbor wasn't already in Q
      r[neighbor] += val;
       if(rPrev < eps * neighbor_degree && r[neighbor] > eps * neighbor_degree)
	 Q.push(neighbor);
    }

}

/*
  Runs a random walk where the walker is equally likely to take any
  one of the outgoing links or stay on the node they are at.

  The vector distribution sets the starting distribution of places the
  walker can start from. If values in distribution do not sum to one,
  then it will be normalized.  If values are negative, exit with
  error.  num_steps is the number of steps the walker takes

  Returns vector of probability that walker is at each node after num_steps.

  Added 2014
 */
vector< vector<double> > random_walk(Graph& G, vector<double>& initial_distribution, const int num_steps)
{
  // First make sure initial_distribution is a distribution 
  
  // Need to initialize a few variables, like number of vertices in graph
  double sum_val = 0;
  v_size_t size = num_vertices(G);

  // Note, if initial_distribution is larger than number of vertices, values
  // at indices past the number of vertices are ignored
  initial_distribution.resize(size, 0.0);

  for (size_t i = 0; i < size; ++i)
    {
      if (initial_distribution[i] < 0)
	{
	  cout<<"Distribution must be nonnegative!\n";
	  exit(EXIT_FAILURE);
	}
      sum_val += initial_distribution[i];
    }

  if (sum_val == 0)
    cout<<"Warning. Distribution is all zeros, using uniform distribution\n";
    
  for (size_t i = 0; i < size; ++i)
    {
      if (sum_val == 0)
	initial_distribution[i] = 1.0 / double(size);
      else
	initial_distribution[i] /= sum_val;
    }

  // Now perform walk

  vector< vector<double> > diffusion;
  vector<double> prev_walk = initial_distribution;
  graph_traits<Graph>::degree_size_type walk_options;
  graph_traits<Graph>::adjacency_iterator ait, aite;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);

  for (int i = 0; i < num_steps; ++i)
    {
      vector<double> walk(size, 0);
      for (size_t j = 0; j < size; ++j)
	{
	  // Only need to distribute from nodes with nonzero mass
	  if (prev_walk[j] != 0)
	    {
	      // Find out degree + 1
	      Vert v = vertex(j, G);
	      walk_options = out_degree(v, G) + 1;
	      
	      // First add probability of staying put.
	      double diffusion_value = prev_walk[j] / double(walk_options);
	      walk[j] += diffusion_value;
	      
	      // Now add probability of moving to adjacent nodes
	      for (tie(ait, aite) = adjacent_vertices(v, G); ait != aite; ++ait)
		{
		  Vert u = *ait;
		  walk[index[u]] += diffusion_value;
		}
	    }
	}
      prev_walk = walk;
      diffusion.push_back(walk);
    }
  return diffusion;
}
