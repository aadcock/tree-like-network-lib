/*
Main file for performing label propagation on a graph. Currently only
allows input labels. Later, input distributions may be allowed.

By Aaron Adcock
Stanford University, PhD Candidate
Dec. 2013
Updated June, 2014

Executable options:

-i input graph...REQUIRED
-labels vertices for which labels are known...REQUIRED
-o output file...default = "./label_prop.out"
-core k_min...run label prop on nodes in cores greater than k_min
-periph k_max...run label prop on nodes only in cores less than k_max
-d if y, indicates a directed graph...default = undirected graph
-s if y, suppresses terminal output...default = false
*/

#include "../lib/graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <map>
#include <limits>

int main(int argc, char *argv[])
{
  string input_graph = "";
  string input_label_file = "";
  string output_file = "./label_prop.out";
  string seed_type = "random";
  string sort_type = "descend";
  float train_percentage = .10;
  bool directed = false;
  bool suppress_output = false;
  int max_iterations = 10;
  //User input

  if(argc < 2)
    {
      cout<<"-i <input_graph_file> \n-labels <input_labels_file> \n-o <output_file> \n-max_iter <max_iterations> \n-seed <deg for degrees, core for cores, rand for random> \n -sort <ascend: sort by ascending seed value, descend: sort by descending seed value> \n-train <fraction of labels to keep for training> \n-d <y=directed graph input> \n-s <y=suppress terminal output>\n";
      return(-1);
    }

  for(int i=1;i<argc;i++)
    {
      stringstream ss;
      if (strcmp(argv[i], "-i") == 0)
	input_graph = argv[i + 1];
      else if (strcmp(argv[i], "-labels") == 0)
	input_label_file = argv[i + 1];
      else if (strcmp(argv[i],"-o") == 0)
	output_file = argv[i + 1];
      else if (strcmp(argv[i], "-max_iter") == 0)
	max_iterations = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-seed") == 0)
	seed_type = argv[i + 1];
      else if (strcmp(argv[i], "-sort") == 0)
	sort_type = argv[i + 1];
      else if (strcmp(argv[i], "-train") == 0)
	train_percentage = atof(argv[i + 1]);
      else if (strcmp(argv[i], "-d") == 0)
	{	
	  if (strcmp(argv[i+1], "y") == 0)
	    directed = true;
	}
      else if (strcmp(argv[i], "-s") == 0)
	{	  
	  if (strcmp(argv[i+1], "y") == 0)
	    suppress_output = true;
	}

    }

  ////Default values

  if (input_graph.length()==0)
    {
      cerr<<"You must provide an input graph file\n";
      exit(EXIT_FAILURE);
    }
  else if (input_label_file.length()==0)
    {
      cerr<<"You must provide an input labels file\n";
      exit(EXIT_FAILURE);
    }
  else if (!suppress_output)
    cout<<"Input File: "<<input_graph<<"\n";

  /* Label propagation algorithm (I may break this section of code
     into a separate file so that this just contains the input/output)
     
     How this code should work:

     For the labeled vertices, assign a fixed value distribution (ie
     p(label) = 1, p(!label) = 0).

     For the unlabeled vertices, assign a uniform distribution (ie
     p(x) = 1/n where n is the number of labels)

     while (not converged, or max_iterations not reached)
     {
          for each vertex
	  {
               Get neighbor's distributions;
	       Average them; 
	       Set vertex distribution to average;
          }
     }

     Output distribution?
     Output top label?
     
   */

  // Get graph, edge direction reversed will get connected component
  // later as we want to make sure labels are applied beforehand
  Graph G = loadGraph(input_graph);
  cout<<"Graph loaded\n";

  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  
  ifstream label_file(input_label_file.c_str());

  vector<int> labels(num_vertices(G), -1);
  vector<Vert> labeled_vertices;
  map<int, vector<Vert> > vertices_by_label;
  labeled_vertices.reserve(num_vertices(G));
  if (label_file.is_open())
    {
      string line;
      int line_number = 0;
      int labels_so_far = 0;
      while (getline(label_file, line))
	{
	  line_number++;
	  // Comment character is #
	  if (line.substr(0,1) != "#")
	    {
	      istringstream line_stream(line);
	      int vertex_id, vertex_label;

	      if (line_stream>>vertex_id>>vertex_label)
		{
		  if (vertex_label >= 0)
		    labeled_vertices.push_back(vertex(vertex_id - 1, G));
		  
		  labels[vertex_id - 1] = vertex_label;
		  
		  if (vertices_by_label.find(vertex_label) == vertices_by_label.end())
		    {
		      vector<Vert> initial;
		      initial.push_back(vertex(vertex_id - 1, G));
		      vertices_by_label[vertex_label] = initial;
		    }
		  else
		    vertices_by_label[vertex_label].push_back(vertex(vertex_id - 1, G));
		}
	    }
	}
    }

  srand(time(NULL));
  size_t num_train_labels = (size_t) (labeled_vertices.size() * train_percentage);
  vector<int> train_labels(num_train_labels, -1);
  vector<Vert> train_vertices;

  // Pick labels to keep based on seed string/randomness
  if (strcmp(seed_type.c_str(), "deg") == 0)
    {
      // Need to order nodes based on degree, first need to calculate
      // degrees.
      //      vector<size_t> degrees_labeled_vertices;
      //
      //degrees_labeled_vertices.reserve(labeled_vertices.size());
      //for (size_t i = 0; i < labeled_vertices.size(); ++i)
      //degrees_labeled_vertices.push_back(out_degree(labeled_vertices[i], G));

      // if (sort_type == "ascend")
      // 	sort(labeled_vertices.begin(), labeled_vertices.end(), compareIndltg<vector<size_t> >(degrees_labeled_vertices));
      // else
      // 	sort(labeled_vertices.begin(), labeled_vertices.end(), compareIndgtl<vector<size_t> >(degrees_labeled_vertices));

      // // load top of sorted, labeled vertices
      // for (size_t i = 0; i < num_train_labels; ++i)
      // 	{
      // 	  v_size_t vind = labeled_vertices[i];
      // 	  train_vertices.push_back(vertex(labeled_vertices[i], G));
      // 	  train_labels[i] = labels[vind];
      // 	}
      //      map<int, vector<size_t> > degrees_by_label;
      vector<int> labels;
      for (map<int, vector<size_t> >::iterator it = vertices_by_label.begin(); it != vertices_by_label.end(); ++it)
	{
	  int label = it->first;
	  labels.push_back(label);
	  vector<Vert> vertices = it->second;
	  vector<size_t> temp_degrees(vertices.size(), 0);
	  //degrees_by_label[label] = degrees;

	  for (size_t i = 0; i < vertices.size(); ++i)
	    temp_degrees[i] = out_degree(vertices[i], G);

	  if (sort_type == "ascend")
	    sort(vertices_by_label[label].begin(), vertices_by_label[label].end(), compareIndltg<vector<size_t> >(temp_degrees));
	  else
	    sort(vertices_by_label[label].begin(), vertices_by_label[label].end(), compareIndgtl<vector<size_t> >(temp_degrees));
	}

      int num_seeds_per_label = num_train_labels / labels.size();
      size_t count = -1;
      size_t adjust = 0;
      cout<<labels.size()<<"\n";
      for (size_t i = 0; i < num_train_labels; ++i)
	{
	  if (i % num_seeds_per_label == 0)
	    count++;

	  if (count == labels.size())
	    {
	      count = 0;
	      adjust += num_seeds_per_label;
	    }

	  size_t ind = i % num_seeds_per_label + adjust;
	  v_size_t vind = vertices_by_label[labels[count]][ind];
	  train_vertices.push_back(vertex(vind, G));
	  train_labels[i] = labels[count];
	}
    }
  else if (strcmp(seed_type.c_str(), "core") == 0)
    {
      // vector<int> core = k_core(G);
      // vector<double> labeled_core;

      // for (size_t i = 0; i < labeled_vertices.size(); ++i)
      // 	labeled_core.push_back(core[labeled_vertices[i]]);

      // vector<size_t> degrees_labeled_vertices;
      // size_t max_degree = 0;
      // degrees_labeled_vertices.reserve(labeled_vertices.size());
      // for (size_t i = 0; i < labeled_vertices.size(); ++i)
      // 	{
      // 	  size_t deg = out_degree(labeled_vertices[i], G);
      // 	  degrees_labeled_vertices.push_back(deg);
      // 	  if (deg > max_degree)
      // 	    max_degree = deg;
      // 	}

      // // adjust core by degree, thus ordering by core then degree
      // for (size_t i = 0; i < labeled_vertices.size(); ++i)
      // 	labeled_core[i] += (double) degrees_labeled_vertices[i] / ((double) 2*max_degree);

      // if (sort_type == "ascend")
      // 	sort(labeled_vertices.begin(), labeled_vertices.end(), compareIndltg<vector<double> >(labeled_core));
      // else
      // 	sort(labeled_vertices.begin(), labeled_vertices.end(), compareIndgtl<vector<double> >(labeled_core));

      // // load top of sorted, labeled vertices
      // for (size_t i = 0; i < num_train_labels; ++i)
      // 	{
      // 	  v_size_t vind = labeled_vertices[i];
      // 	  train_vertices.push_back(vertex(labeled_vertices[i], G));
      // 	  train_labels[i] = labels[vind];
      // 	}
      vector<int> labels;
      vector<int> core = k_core(G);
      for (map<int, vector<size_t> >::iterator it = vertices_by_label.begin(); it != vertices_by_label.end(); ++it)
	{
	  int label = it->first;
	  labels.push_back(label);
	  vector<Vert> vertices = it->second;
	  vector<size_t> temp_cores(vertices.size(), 0);
	  //degrees_by_label[label] = degrees;

	  for (size_t i = 0; i < vertices.size(); ++i)
	    temp_cores[i] = (size_t) core[index[vertices[i]]];

	  if (temp_cores.size() > 850)
	    cout<<temp_cores[850]<<"\n";

	  if (sort_type == "ascend")
	    sort(vertices_by_label[label].begin(), vertices_by_label[label].end(), compareIndltg<vector<size_t> >(temp_cores));
	  else
	    sort(vertices_by_label[label].begin(), vertices_by_label[label].end(), compareIndgtl<vector<size_t> >(temp_cores));
	}

      int num_seeds_per_label = num_train_labels / labels.size();
      size_t count = -1;
      size_t adjust = 0;
      cout<<labels.size()<<"\n";
      for (size_t i = 0; i < num_train_labels; ++i)
	{
	  if (i % num_seeds_per_label == 0)
	    count++;

	  if (count == labels.size())
	    {
	      count = 0;
	      adjust += num_seeds_per_label;
	    }

	  size_t ind = i % num_seeds_per_label + adjust;
	  v_size_t vind = vertices_by_label[labels[count]][ind];
	  train_vertices.push_back(vertex(vind, G));
	  train_labels[i] = labels[count];
	}
    }
  else
    {
      for (int i = 0; i < num_train_labels; ++i)
	{
	  size_t rand_ind = rand() % labeled_vertices.size();
	  train_vertices.push_back(vertex(labeled_vertices[rand_ind], G));
	  train_labels[i] = labels[labeled_vertices[rand_ind]];
	}
    }

  for (int i = 0; i < train_vertices.size(); ++i)
    cout<<"Vertex: "<<train_vertices[i]<<" label: "<<train_labels[i]<<"\n";

  vector<int> predicted_labels;
  label_propagation(G, train_vertices, train_labels, predicted_labels, max_iterations);

  cout<<"Classification accuracy: "<<label_accuracy(predicted_labels, labels)<<"\n";
}
