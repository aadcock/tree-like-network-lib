/*
  File containing machine learning routines for use in this graph
  library.  First routines to be added are label propagation and a
  classification accuracy routine.  Later may use cross-validation
  routine.

  By Aaron Adcock
  Stanford University, PhD Candidate (defended)
  June 2014
*/

#include "../graph_lib_boost.hpp"
#include <algorithm>
/*
  Function outputs a classification vector using label propagation
  seeded with the labels in known_labels.  The known_labels vector
  must be of the same length as the labeled_vertices vector.  The
  labeled vertices vector gives the associated vertex for each label
  in known_labels.  Output is a vector of labels.  Currently, we make
  use of a vertex property called distribution, which is a vector of
  floats representing the probability of each label for that vector.
  Currently we have no modeling, we just assume each vector is equally
  likely.  In the future, we could set the distribution vector for
  each node.

  Currently output is placed in predicted labels and in the
  distribution vector of G.  The distribution vector of G worries me,
  I think that should be changed (ie, what if I want to keep
  something in G[v].distribution)
*/

bool label_propagation(Graph & G, vector<Vert> & labeled_vertices, vector<int> & known_labels, vector<int> & predicted_labels, size_t num_iter)
{
  if (labeled_vertices.size() != known_labels.size())
    {
      cerr<<"The labeled vertex vector and the known labels vector must be the same size.\n";
      return false;
    }

  if (num_iter == 0)
    {
      cerr<<"Number of iterations for label propagation must be greater than 0.  Setting to 50.\n";
      num_iter = 50;
    }
  // Find how many labels are in vector.
  // I will assume that the labels come as ints, but not necessarily
  // ints that begin at zero.  Thus two converters.
  map<int, size_t> label_to_integer;
  map<size_t, int> integer_to_label;

  vector<int> unique_labels = known_labels;
  sort(unique_labels.begin(), unique_labels.end());
  vector<int>::iterator it = unique(unique_labels.begin(), unique_labels.end());
  unique_labels.resize(distance(unique_labels.begin(), it));

  size_t num_labels = 0;
  for (size_t i; i < unique_labels.size(); ++i)
    {
      label_to_integer[unique_labels[i]] = num_labels;
      integer_to_label[num_labels] = unique_labels[i];
      ++num_labels;
    }
  
  // Run through vertices and update each distribution to be uniform
  // (unlabeled) or deterministic (labeled)
  graph_traits<Graph>::vertex_iterator vit, vite;
  property_map<Graph, vertex_index_t>::type index = get(vertex_index, G);
  vector<float> uniform_distribution(num_labels, 1.0 / (float) num_labels);
  vector<float> zero_distribution(num_labels, 0);
  vector<int> is_labeled(num_vertices(G), 0);
  
  // set distributions
  for (tie(vit, vite) = vertices(G); vit != vite; ++vit)
    {
      Vert v = *vit;
      v_size_t ind = index[v];
      G[v].distribution = uniform_distribution;      
    }

  for (size_t i = 0; i < labeled_vertices.size(); ++i)
    {
      Vert v = labeled_vertices[i];
      is_labeled[index[v]] = 1;
      size_t ind = label_to_integer[known_labels[i]];
      G[v].distribution = zero_distribution;
      G[v].distribution[ind] = 1.0;
    }
  
  cout<<"Initial Distributions set\n";
    
  // Begin iterations
  // For now, use 50 iterations, later use a convergence criteria +
  // user specified or default max iterations

  graph_traits<Graph>::adjacency_iterator ait, aite;
  vector<float> new_distribution(num_labels, 0);
  for (int i = 0; i < num_iter; ++i)
    {
      cout<<"On iteration: "<<i<<".\n";
      for (tie(vit, vite) = vertices(G); vit != vite; ++vit)
	{
	  Vert v = *vit;
	  v_size_t ind = index[v];
	  if (is_labeled[ind] != 1)
	    {
	      new_distribution.assign(num_labels, 0);
	      int degree_v = out_degree(v, G);
	      for (tie(ait, aite) = adjacent_vertices(v, G); ait != aite; ++ait)
		{
		  Vert u = *ait;
		  for (int j = 0; j < new_distribution.size(); ++j)
		    new_distribution[j] += G[u].distribution[j] / degree_v;
		}
	      G[v].distribution = new_distribution;
	    }
	}
    }

  predicted_labels.reserve(num_vertices(G));
  predicted_labels.clear();
  predicted_labels.assign(num_vertices(G), -1);
  
  for (tie(vit, vite) = vertices(G); vit != vite; ++vit)
    {
      Vert v = *vit;
      v_size_t ind = index[v];
      vector<float> & label_distribution = G[v].distribution;
      
      float max_prob = 0;
      int max_label = -1;
      for (size_t i = 0; i < label_distribution.size(); ++i)
	{
	  if (label_distribution[i] > max_prob)
	    {
	      max_prob = label_distribution[i];
	      max_label = i;
	    }
	}

      predicted_labels[ind] = integer_to_label[max_label];
    }

  return true;
}

/*
  Gives accuracy of a labeling in predicted_labels as compared to the
  real labels.  Note, labels must be positive.  Negative values are
  discounted.  Thus, for every vertex a predicted and real entry must
  be provided, though if it is negative it will be ignored.
 */
double label_accuracy(vector<int> & predicted_labels, vector<int> & labels)
{
  if (predicted_labels.size() != labels.size())
    {
      cerr<<"There must be the same number of predicted labels as actual labels.\n";
      return -1;
    }

  double total_predictions = 0;
  double total_correct = 0;
  for (size_t i = 0; i < predicted_labels.size(); ++i)
    {
      if (predicted_labels[i] >= 0 && labels[i] >= 0)
	{
	  total_predictions++;
	  if (predicted_labels[i] == labels[i])
	    total_correct++;
	}
    }

  return total_correct / total_predictions;
}
