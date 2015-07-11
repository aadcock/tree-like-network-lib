/*
Loads a graph from an edgelist in the provided text file The loaded
graph's nodes are reindexed from 1 to the size of the graph.  Assumes
the graph is undirected unless 'directed' flag is provided and set to
true.

Also contains load_color_file loads a vector of ints or doubles from
file.  These can be used to color a graph in a visualization.

Latest version includes functions which parse arbitrary whitespace,
dimacs files, and the INDDGO dimacs-tree-decomposition format --2013

Code to load INDDGO tree decompositions has been moved to
loadTreeDecomposition.cpp -- Nov. 2013

Aaron Adcock
PhD Candidate Stanford University
2011
*/

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include "../graph_lib_boost.hpp"

using namespace boost;

/////////////////////////////////////
//Binary Search
//Searches through a sorted vector for a specific key
//Returns -1 if key not found
/////////////////////////////////////

int binarySearch(vector<int>& sortedVector, int key)
{
  //Check if sortedVector is empty
  if(sortedVector.size()==0)
    return(-1);

  //Declare variables, l is length of vector
  //k is midpoint of current range [klower,kupper)
  const int l = sortedVector.size();
  int k = l/2;
  bool flag = true;
  int kupper = l;
  int klower = 0;

   
  while(sortedVector[k] != key && (klower < kupper))
    {
      flag = false;

      if(sortedVector[k] > key)
	kupper = k;
      else if(sortedVector[k] < key)
	klower = k;
       
      if((kupper-klower)<=1 && k==klower)
	{
	  klower = kupper;
	}
      else
	{
	  k = klower + (kupper-klower)/2;

	  if(sortedVector[k]==key)
	    flag = true;
	}
    }

  if(flag==true)
    {
      return(k);
    }
  else
    {
      return(-1);
    }

}

////////////////////////////
//Searches through a set of vertices
//v_size_t declared in ../graph_lib_boost.hpp
///////////////////////////

int binarySearch(vector<v_size_t>& sortedVector, v_size_t key)
{
  if(sortedVector.size()==0)
    return(-1);
 
  long l = sortedVector.size();
  long k = l/2;
  bool flag = true;
  long kupper = l;
  long klower = 0;

   
  while(sortedVector[k] != key && (klower < kupper))
    {
      flag = false;

      if(sortedVector[k] > key)
	kupper = k;
      else if(sortedVector[k] < key)
	klower = k;
       
      if((kupper-klower)<=1 && k==klower)
	{
	  klower = kupper;
	}
      else
	{
	  k = klower + (kupper-klower)/2;

	  if(sortedVector[k]==key)
	    flag = true;
	}
    }

  if(flag==true)
    {
      return(k);
    }
  else
    {
      return(-1);
    }

}

////////////////////////////////////////////////////////////////////
//Loads a graph from a vector of intPairs which represent
//a edgelist
///////////////////////////////////////////////////////////////////

Graph loadGraph(vector<intPair> edgeList)
{
  vector<intPair>::iterator it;
  list<int> uniqueList;
  vector<int> sortedVector;
  
  // Make a list of vertices (starts off with repeats) by running through edges
  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge=*it;
      uniqueList.push_back(currEdge.first);
      uniqueList.push_back(currEdge.second);
    }

  // Sort, remove repeated vertices, put in vector
  uniqueList.sort();
  uniqueList.unique();
  sortedVector.insert(sortedVector.end(),uniqueList.begin(),uniqueList.end());

  // Create graph of appropriate size
  Graph G(uniqueList.size());

  // Add edges to graph
  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge =*it;
      intPair newEdge;
      int newid1;
      int newid2;
      bool flag = false;

      // Check to see if id is in the unique vector of nodes,
      // Just a sanity check
      newid1 = binarySearch(sortedVector,currEdge.first);
      if(newid1<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      newid2 = binarySearch(sortedVector,currEdge.second);
      if(newid2<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      if(flag)
	{
	  Vert v1 = vertex(newid1,G);
	  Vert v2 = vertex(newid2,G);
	  
	  add_edge(v1,v2,G);
	}
    }

  return(G);
}

//////////////////////////////////////////////////// 
// The next two functions are shorthand for reading
// an edgelist from file (should work on any whitespace
// delimited edgelist
///////////////////////////////////////////////////
Graph loadGraph(string filename)
{
  return loadGraph(filename,"\t",0);
}

Graph loadGraph(string filename,string delimiter)
{
  return loadGraph(filename,delimiter,0);
}

/////////////////////////////////////////////////////////////////////
// This is the main function for loading a graph from an edgelist.
// The first argument gives the full filepath, the second gives the delimiter
// the third argument tells whether to load this as a directed/undirecte graph
/////////////////////////////////////////////////////////////////////

Graph loadGraph(string filename, string delimiter, const bool directed)
{

  ifstream graphFile;  
  string line;
  vector<int> sortedVector;
  list<intPair> edgeList;
  int l = delimiter.length();
  list<int> vertices;
  list<int> uniqueList;
  int num_vert = 0;

  // This assumes that the graph file has # for unused lines and that
  // it is in the edge list format (lists of node *delimiter* node)
  // Winter '12, This should now also work for dimacs functions,
  // though you can use an explicit call to the loadDimacs function
  // also included in this file.


  // Open file
  graphFile.open(filename.c_str());

  // Check if open succeeded
  if (graphFile.is_open())
    {
      // flag keeps track if file is dimacs file
      bool flag = true;
      while (graphFile.good())
	{
	  getline(graphFile,line);
	  
	  // checking for dimacs (starts file with 'p num_nodes num_edges'
	  if(line.substr(0,1)=="p")
	      flag = false;

	  // checking for dimacs (has an 'e' in front of each edge)
	  if(line.substr(0,1) == "e")
	    {
	      flag = false;
	      
	      // look for spaces (delimiter in dimacs file)
	      size_t pos1 = line.find(" ");
	      size_t pos2 = line.find(" ",pos1+1);

	      if(pos1 != string::npos)
		{
		  // Find end of line, plus first two numbers
		  size_t pos3 = line.find("\n");
		  string node1id = line.substr(pos1+1,pos2);
		  string node2id = line.substr(pos2+1,pos3-pos2);
		  intPair currEdge;
		  istringstream id1(node1id,istringstream::in);
		  istringstream id2(node2id,istringstream::in);

		  // get nodeid's ([1,n] in dimacs format)
		  id1>>currEdge.first;
		  id2>>currEdge.second;
		  vertices.push_back(currEdge.first);
		  vertices.push_back(currEdge.second);

		  if(currEdge.first != currEdge.second)
		    edgeList.push_back(currEdge);

		  if(!directed && (currEdge.first != currEdge.second))
		    {
		      int temp = currEdge.second;

		      currEdge.second = currEdge.first;
		      currEdge.first  = temp;
		      
		      edgeList.push_back(currEdge);
		    }
		}
	    } 
	  // not dimacs
	  else if(line.substr(0,1) != "#" && flag)
	    {

	      // This looks for an arbitrary number of spaces because of the output of a Fortran code
	      string::iterator sit;
	      int nodeid1;
	      int nodeid2;
	      
	      bool one_int_found  = false;
	      bool read_an_int = false;
	      bool two_int_found  = false;

	      string node_string = "";
	      
	      // Read until two tokens are found 
	      for(int i=0; i<line.length(); ++i)
		{
		  char c = line[i];

		  if(isdigit(c))
		    {
		      node_string.append(line.substr(i,1));
		      read_an_int = true;
		    }
		  else if(read_an_int)
		    {

		      if(!one_int_found)
			{
			  nodeid1 = atoi(node_string.c_str());
			  node_string = "";
			  one_int_found=true;
			}
		      else if(!two_int_found)
			{
			  nodeid2 = atoi(node_string.c_str());
			  node_string = "";
			  two_int_found=true;
			}
		      read_an_int = false;
		    }
		}

	      if(read_an_int && !two_int_found)
		{
		  nodeid2 = atoi(node_string.c_str());
		  node_string = "";
		  two_int_found=true;
		}

	      if(!one_int_found || !two_int_found)
		{
		  cerr<<"Warning, Line in graph file does not contain at least two integers, moving to next line\n";
		}
	      else
		{
		  intPair currEdge;
		  currEdge.first  = nodeid1;
		  currEdge.second = nodeid2;
		  vertices.push_back(currEdge.first);
		  vertices.push_back(currEdge.second);

		  if(currEdge.first != currEdge.second)
		    edgeList.push_back(currEdge);
		  
		  if(!directed && (currEdge.first != currEdge.second))
		    {
		      int temp = currEdge.second;
		      currEdge.second = currEdge.first;
		      currEdge.first  = temp;
		      
		      edgeList.push_back(currEdge);
		    }
		}

	      //  size_t pos1 = line.find(delimiter);
	      //  if(pos1 != string::npos)
	      //  	{
	      //  	  size_t pos2 = line.find("\n");
	      //  	  string node1id = line.substr(0,pos1);
	      //  	  string node2id = line.substr(pos1+l,pos2-pos1);
	      //  	  intPair currEdge;
	      //  	  istringstream id1(node1id,istringstream::in);
	      //  	  istringstream id2(node2id,istringstream::in);
		  
	      //  	  id1>>currEdge.first;
	      //  	  id2>>currEdge.second;
	      //  	  vertices.push_back(currEdge.first);
	      //  	  vertices.push_back(currEdge.second);
		  
	      //  	  edgeList.push_back(currEdge);
		  
	      //  	  if(!directed)
	      //  	    {
	      //  	      int temp = currEdge.second;

	      //  	      currEdge.second = currEdge.first;
	      //  	      currEdge.first  = temp;
		      
	      //  	      edgeList.push_back(currEdge);
	      //  	    }
	      //  	}
	    }

	}
    }
  else
    {
      std::cerr<<"Error opening graph file \n"<<filename.c_str();
      exit(EXIT_FAILURE);
    }
  graphFile.close();

  vertices.sort();
  vertices.unique();
  num_vert = vertices.size();
  // Sorted vector contains all nodeids for sorting
  list<intPair>::iterator it;
  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge=*it;
      uniqueList.push_back(currEdge.first);
      uniqueList.push_back(currEdge.second);
    }
  
  uniqueList.sort();
  uniqueList.unique();

  sortedVector.insert(sortedVector.end(),uniqueList.begin(),uniqueList.end());

  // Create graph of appropriate size
  Graph G(uniqueList.size());

  // Add edges
  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge =*it;
      intPair newEdge;
      int newid1;
      int newid2;
      bool flag = false;

      // Search for node in nodelist as sanity check
      newid1 = binarySearch(sortedVector,currEdge.first);
      if(newid1<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      newid2 = binarySearch(sortedVector,currEdge.second);
      if(newid2<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      // If the searches are succesful (real nodes) then 
      // check if edges already exists.  If it doesn't, add the edge
      if(flag)
	{
	  Vert v1 = vertex(newid1,G);
	  Vert v2 = vertex(newid2,G);
	  bool not_edge = true;

	  graph_traits<Graph>::adjacency_iterator ai, aie;
	  for(tie(ai,aie) = adjacent_vertices(v1,G);ai!=aie;ai++)
	    {
	      Vert u = *ai;
	      if(v2==u)
		not_edge = false;
	    }
	  if(not_edge)
	    add_edge(v1,v2,G);
	}
    }
  return(G);
}

//////////////////////////////////////////////// 
// The following functions explicitly look for a dimacs format
///////////////////////////////////////////////
Graph loadDimacsGraph(string filename)
{
  return loadDimacsGraph(filename,0);
}


Graph loadDimacsGraph(string filename, const bool directed)
{

  ifstream graphFile;  
  string line;
  vector<int> sortedVector;
  list<intPair> edgeList;
  list<int> vertices;
  list<int> uniqueList;
  int num_vert = 0;
  const string delimiter = " ";
  const int l = delimiter.length();

  // This assumes that the graph file has a letter at the beginning of each line and
  // in accordance with the dimacs graph format

  graphFile.open(filename.c_str());
  if (graphFile.is_open())
    {
      while (graphFile.good())
	{
	  getline(graphFile,line);
	  if(line.substr(0,1) == "e")
	    {
	      size_t pos1 = line.find(delimiter);
	      size_t pos2 = line.find(delimiter,pos1+1);
	      if(pos1 != string::npos)
		{
		  size_t pos3 = line.find("\n");
		  string node1id = line.substr(pos1+1,pos2);
		  string node2id = line.substr(pos2+l,pos3-pos2);
		  intPair currEdge;
		  istringstream id1(node1id,istringstream::in);
		  istringstream id2(node2id,istringstream::in);

		  id1>>currEdge.first;
		  id2>>currEdge.second;
		  vertices.push_back(currEdge.first);
		  vertices.push_back(currEdge.second);

		  edgeList.push_back(currEdge);

		  if(!directed)
		    {
		      int temp = currEdge.second;

		      currEdge.second = currEdge.first;
		      currEdge.first  = temp;
		      
		      edgeList.push_back(currEdge);
		    }
		}
	    }

	}
    }
  else
    {
      std::cerr<<"Error opening graph file \n", filename.c_str();
      exit(EXIT_FAILURE);
    }
  graphFile.close();

  vertices.sort();
  vertices.unique();
  num_vert = vertices.size();
  // Sorted vector contains all nodeids for sorting
  list<intPair>::iterator it;
  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge=*it;
      uniqueList.push_back(currEdge.first);
      uniqueList.push_back(currEdge.second);
    }
  
  uniqueList.sort();
  uniqueList.unique();

  sortedVector.insert(sortedVector.end(),uniqueList.begin(),uniqueList.end());

  Graph G(uniqueList.size());

  for(it=edgeList.begin();it!=edgeList.end();it++)
    {
      intPair currEdge =*it;
      intPair newEdge;
      int newid1;
      int newid2;
      bool flag = false;

      newid1 = binarySearch(sortedVector,currEdge.first);
      if(newid1<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      newid2 = binarySearch(sortedVector,currEdge.second);
      if(newid2<0)
	std::cerr<<"Error in GraphLoad finding id\n";
      else
	  flag = true;

      if(flag)
	{
	  Vert v1 = vertex(newid1,G);
	  Vert v2 = vertex(newid2,G);
	  bool not_edge = true;

	  graph_traits<Graph>::adjacency_iterator ai, aie;
	  for(tie(ai,aie) = adjacent_vertices(v1,G);ai!=aie;ai++)
	    {
	      Vert u = *ai;
	      if(v2==u)
		not_edge = false;
	    }
	  if(not_edge)
	    add_edge(v1,v2,G);
	}
    }
  return(G);
}
