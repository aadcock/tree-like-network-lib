/*
  Loads a tree decompostion from a modified dimacs file (such as the one
  produced by the INDDGO open source software). Loads the tree
  decomposition into the TreeDecomp format defined in
  tree_lib_boost.hpp.

  Aaron Adcock
  PhD Candidate Stanford University
  2013
*/

#include "../graph_lib_boost.hpp"
#include "tree_lib_boost.hpp"

using namespace boost;

/*
  This function loads a tree decomposition located at the location of
  filename.  The tree decomposition must be stored in a modified
  dimacs file.
*/

TreeDecomp loadTreeDecomp(Graph & G, string filename)
{
  ifstream treeFile;  
  string line;
  vector<int> sortedVector;
  list<intPair> edgeList;
  list<v_size_t> vertices;
  list<v_size_t> uniqueList;
  v_size_t num_bags = 0;
  v_size_t listed_num_bags = 0;
  const string delimiter = " ";
  const int l = delimiter.length();
  
  vector< vector<v_size_t> > bags;
  vector<size_t> size_t_dummy;
  vector< vector<size_t> > vertex_locations(num_vertices(G), size_t_dummy);

  // This assumes that the graph file has a letter at the beginning of
  // each line and in accordance with the dimacs graph format.  In the
  // modified dimacs, tree decomposition format, B indicates a list of
  // nodes in a bag.  The bag format is
  //
  // B n u v w ...
  //
  // where B is the indicator character, n is the bag number, u v w
  // ... are the nodes in the bag.  The first line must contain the
  // number of nodes, or else an error is shown.

  bool num_bags_found = false;
  int count_bags = 0;
  treeFile.open(filename.c_str());
  if (treeFile.is_open())
    {
      while (treeFile.good())
	{
	  getline(treeFile,line);
	  if(line.substr(0,1) == "e")
	    {
	      if (!num_bags_found)
		{
		  cerr<<"The first line did not begin with p treed num_bags num_edges\n";
		  exit(EXIT_FAILURE);
		}
	      size_t delim1_pos = line.find(delimiter);
	      size_t delim2_pos = line.find(delimiter,delim1_pos+1);
	      if(delim1_pos != string::npos)
		{
		  size_t end_of_line = line.find("\n");
		  string node1id = line.substr(delim1_pos+1,delim2_pos);
		  string node2id = line.substr(delim2_pos+l,end_of_line-delim2_pos);
		  intPair currEdge;
		  istringstream id1(node1id,istringstream::in);
		  istringstream id2(node2id,istringstream::in);
		  int temp1;
		  int temp2;
		  id1>>temp1;
		  id2>>temp2;
		  
		  currEdge.first  = temp1;
		  currEdge.second = temp2;
		  vertices.push_back(currEdge.first);
		  vertices.push_back(currEdge.second);

		  edgeList.push_back(currEdge);

		  int temp = currEdge.second;

		  currEdge.second = currEdge.first;
		  currEdge.first  = temp;
		      
		  edgeList.push_back(currEdge);
		    
		}
	    }
	  else if(line.substr(0,1) == "p")
	    {
	      istringstream line_stream(line);
	      string token_1,token_2;
	      int token_3, token_4;

	      try
		{
		  line_stream>>token_1>>token_2>>token_3>>token_4;
		}catch (...)
		{
		  cerr<<"Does not follow INDDGO tree format.  p treed bags edges not first line\n";
		  exit(EXIT_FAILURE);
		}
	      
	      

	      // 	      size_t end_of_line = line.find("\n");
	      // string node1id = line.substr(delim1_pos+1,delim2_pos);
	      // string node2id = line.substr(delim2_pos+l,end_of_line-delim2_pos);
	      // intPair currEdge;
	      // istringstream id1(node1id,istringstream::in);
	      // istringstream id2(node2id,istringstream::in);
	      num_bags_found = true;
	      listed_num_bags = token_3;  
	    }
	  else if(line.substr(0,1) == "B")
	    {
	      //  cout<<"Begin reading bags, on bag "<<++count_bags<<"\n";
	      istringstream bag_line(line,istringstream::in);
		  
	      char dummy[50];
	      v_size_t bag_id;
	      v_size_t vertex_id;
	      int count = 0;
	      vector<v_size_t> curr_bag;
	      while(bag_line.good())
		{
		  count++;
		  if(count==2)
		    {
		      bag_line>>bag_id;
		      bag_id = bag_id - 1;
		    }
		  else if(count>3)
		    {
		      bag_line>>vertex_id;
		      vertex_id = vertex_id - 1;
		      curr_bag.push_back(vertex_id);
		      vertex_locations[vertex_id].push_back(bag_id);
		    }
		  else if (count == 3)
		    {
		      size_t bag_size;
		      bag_line>>bag_size;
		      curr_bag.reserve(bag_size);
		    }
		  else
		    bag_line>>dummy;
		}
		  
	      bags.push_back(curr_bag);
	    }
	}

    }
  else
    {
      std::cerr<<"Error opening graph file \n", filename.c_str();
      exit(EXIT_FAILURE);
    }
  treeFile.close();
  
  cout<<"Finished loading edge list/bag list from file\n";
  vertices.sort();
  vertices.unique();
  num_bags = vertices.size();

  if (listed_num_bags == 1 && num_bags == 0)
    num_bags = 1;

  if (num_bags != listed_num_bags)
    {
      cerr<<"Number of vertices found does not match number of vertices listed\n";
      exit(EXIT_FAILURE);
    }

  cout<<"Number of bags in tree: "<<num_bags<<"\n";
  Graph T(num_bags);
  // Sorted vector contains all nodeids for sorting
  list<v_size_t>::iterator v_it;
  
  //   uniqueList.sort();
  // uniqueList.unique();

  sortedVector.insert(sortedVector.end(),vertices.begin(),vertices.end());

  // cout<<sortedVector[0]<<"\n";
  list<intPair>::iterator it;
  // cout<<"Finished sorting/initializing tree\n";
  for (it = edgeList.begin(); it != edgeList.end(); ++it)
    {
      intPair currEdge = *it;
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

      // cout<<"Add edge "<<newid1<<","<<newid2<<"\n";
      if (flag)
	{
	  Vert v1 = vertex(newid1,T);
	  Vert v2 = vertex(newid2,T);

	  bool not_edge = true;

	  graph_traits<Graph>::adjacency_iterator ai, aie;
	  for (tie(ai,aie) = adjacent_vertices(v1,T); ai != aie; ++ai)
	    {
	      Vert u = *ai;
	      if (v2 == u)
		not_edge = false;
	    }

	  if(not_edge)
	    add_edge(v1,v2,T);
	}
    }
 
  TreeDecomp treeDecomp(T, G, bags, vertex_locations);
  // cout<<(treeDecomp.get_bags()).size()<<"\n";
  return(treeDecomp);
}
