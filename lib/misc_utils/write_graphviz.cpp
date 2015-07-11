/*
This function is used to take a network (either from an edge list or a
dimacs file) and produce a graphViz .dot file that lays the network
out and colors it as desired.  As desired means we add a vector of
numbers that give the relative color values.  The colors are then
determined using a blue to green to red heat map as the node fill.  We
give the option of highlighting certain nodes (doubling size and
adding a border).

By Aaron Adcock, May, 2013
PhD Candidate at Stanford University
*/

/*
  Converts a value into an rgb char array (ie, a set of 3 hex values).
  The result is a heat map style color function with blue being 0 and
  red being 1.
 */
#include "../graph_lib_boost.hpp"

void get_rgb_value(char rgb[], double val);


/*
  This function produces a graphviz file of the graph passed to it,
  but it colors the nodes based on the colors provided in color.  All
  nodes are same size, though the nodes in the vector square_nodes are
  doubled in size and square as a way of highlighting important nodes.
 */
void write_graphviz(Graph & G, string output_file, vector<double> color, vector<int> square_nodes, bool directed)
{
  ofstream output;
  output.open(output_file.c_str());

  double min = color[0];
  double max = color[0];

  v_size_t n = num_vertices(G);

  for(size_t i = 0; i < color.size(); ++i)
    {
      if(color[i]>max)
	max = color[i];

      if(color[i]<min)
	min = color[i];
    }


  output<<"Graph color {\n";
  output<<"splines=false;\n";

  graph_traits<Graph>::edge_iterator eit, eitend;
  graph_traits<Graph>::vertex_iterator vit, vitend;

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  for (tie(vit,vitend) = vertices(G); vit != vitend; ++vit)
    {

      Vert      v = *vit;
      v_size_t vi = index[v];
      char rgb[8];
      double width=0.2;      
      double val;
      
      if (color.size() != n && G[v].prevIndices.size() != 0)
	val = ((double) color[G[v].prevIndices[0]] - min) / (double) (max - min);
      else
	val = ((double) color[vi] - min) / (double) (max - min);

      get_rgb_value(rgb,val);

      int found = binarySearch(square_nodes,vi);
      
      if(found==-1)
	output<<vi<<"[label=\"\", width="<<width<<", shape=circle, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";
      else
	output<<vi<<"[label=\"\", width="<<2.0*width<<", shape=square, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";

    }
  
  if(directed)
    {
      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert u,v;
	  u = source(e,G);
	  v = target(e,G);

	  v_size_t ui = index[u];
	  v_size_t vi = index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }
  else
    {
      cout<<"Graph viz undirected\n";
      UGraph UG((v_Usize_t) n );

      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert s = source(e,G);
	  Vert t = target(e,G);

	  //	  cout<<"s "<<s<<" t "<<t<<"\n";

	  long ui = index[s];
	  long vi = index[t];

	  //cout<<"ui "<<ui<<" vi "<<vi<<"\n";

	  UVert u = vertex(ui,UG);
	  UVert v = vertex(vi,UG);

	  //cout<<"u "<<u<<" v "<<v<<"\n";

	  UEdge exist;
	  bool yes1,yes2;

	  tie(exist,yes1) = edge(u,v,UG);
	  tie(exist,yes2) = edge(v,u,UG);
	  if(!yes1 && !yes2)
	    add_edge(u,v,UG);
	}

      graph_traits<UGraph>::edge_iterator uit, uitend;
      property_map<UGraph, vertex_index_t>::type u_index = get(vertex_index,UG);

      for(tie(uit,uitend)=edges(UG);uit!=uitend;++uit)
	{
	  UEdge e = *uit;
	  UVert u,v;
	  u = source(e,UG);
	  v = target(e,UG);

	  long ui = u_index[u];
	  long vi = u_index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }

  output<<"}";
  output.close();
}


/*
  This function produces a graphviz file of the graph passed to it,
  but it colors the nodes based on the colors provided in color and
  sizes the nodes based on the node size parameter.  To prevent the
  files from getting out of hand, we make the largest bag 1.1 and the
  smallest bag 0.1.
 */
void write_scaled_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, bool directed)
{
  ofstream output;
  output.open(output_file.c_str());

  double color_min = color[0];
  double color_max = color[0];
  size_t size_min = sizes[0];
  size_t size_max = sizes[0];

  v_size_t n = num_vertices(G);

  if (color.size() != sizes.size() || color.size() != n)
    {
      cout<<"Color vector, size vector must all be of length num_nodes in graph. Exiting. \n";
      exit(EXIT_FAILURE);
    }

  for (size_t i = 0; i < color.size(); ++i)
    {
      if(color[i] > color_max)
	color_max = color[i];

      if(color[i] < color_min)
	color_min = color[i];

      if(sizes[i] > size_max)
	size_max = sizes[i];

      if(sizes[i] < size_min)
	size_min = sizes[i];
    }

  if (color_max > 1 or color_min < 0)
    {
      cout<<"Color vector must be in range [0,1]. Exiting.\n";
      exit(EXIT_FAILURE);
    }

  output<<"Graph color {\n";
  output<<"splines=false;\n";

  graph_traits<Graph>::edge_iterator eit, eitend;
  graph_traits<Graph>::vertex_iterator vit, vitend;

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  for(tie(vit,vitend)=vertices(G);vit!=vitend;++vit)
    {
      Vert      v = *vit;
      v_size_t vi = index[v];
      char rgb[8];
      double width = .75;
      if (size_max != size_min)
	width = 0.1 + ((double) sizes[vi] - size_min) / (double) (size_max - size_min);
      
      double val = 0.5;

      val = color[vi];
      
      get_rgb_value(rgb, val);
     
      output<<vi<<"[label=\"\", width="<<width<<", shape=circle, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";

    }
  
  if(directed)
    {
      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert u,v;
	  u = source(e,G);
	  v = target(e,G);

	  v_size_t ui = index[u];
	  v_size_t vi = index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }
  else
    {
      cout<<"Graph viz undirected\n";
      UGraph UG((v_Usize_t) n );

      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert s = source(e,G);
	  Vert t = target(e,G);

	  //	  cout<<"s "<<s<<" t "<<t<<"\n";

	  long ui = index[s];
	  long vi = index[t];

	  //cout<<"ui "<<ui<<" vi "<<vi<<"\n";

	  UVert u = vertex(ui,UG);
	  UVert v = vertex(vi,UG);

	  //cout<<"u "<<u<<" v "<<v<<"\n";

	  UEdge exist;
	  bool yes1,yes2;

	  tie(exist,yes1) = edge(u,v,UG);
	  tie(exist,yes2) = edge(v,u,UG);
	  if(!yes1 && !yes2)
	    add_edge(u,v,UG);
	}

      graph_traits<UGraph>::edge_iterator uit, uitend;
      property_map<UGraph, vertex_index_t>::type u_index = get(vertex_index,UG);

      for(tie(uit,uitend)=edges(UG);uit!=uitend;++uit)
	{
	  UEdge e = *uit;
	  UVert u,v;
	  u = source(e,UG);
	  v = target(e,UG);

	  long ui = u_index[u];
	  long vi = u_index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }

  output<<"}";
  output.close();
}

/*
  This function produces a graphviz file of the graph passed to it,
  but it colors the nodes based on the colors provided in color and
  sizes the nodes based on the node size parameter.  To prevent the
  files from getting out of hand, we make the largest bag 1.1 and the
  smallest bag 0.1.
 */
void write_scaled_square_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, vector<int> square_nodes, bool directed)
{
  ofstream output;
  output.open(output_file.c_str());

  double color_min = color[0];
  double color_max = color[0];
  size_t size_min = sizes[0];
  size_t size_max = sizes[0];
 
  v_size_t n = num_vertices(G);
  vector<double> node_weight(n, 0);

  if (color.size() != sizes.size() || color.size() != n)
    {
      cout<<"Color vector, size vector must all be of length num_nodes in graph. Exiting. \n";
      exit(EXIT_FAILURE);
    }

  for (size_t i = 0; i < color.size(); ++i)
    {
      if(color[i] > color_max)
	color_max = color[i];

      if(color[i] < color_min)
	color_min = color[i];

      if(sizes[i] > size_max)
	size_max = sizes[i];

      if(sizes[i] < size_min)
	size_min = sizes[i];
    }

  if (color_max > 1 or color_min < 0)
    {
      cout<<"Color vector must be in range [0,1]. Exiting.\n";
      exit(EXIT_FAILURE);
    }

  output<<"Graph color {\n";
  output<<"splines=false;\n";

  graph_traits<Graph>::edge_iterator eit, eitend;
  graph_traits<Graph>::vertex_iterator vit, vitend;

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  for (tie(vit,vitend) = vertices(G); vit != vitend; ++vit)
    {
      Vert      v = *vit;
      v_size_t vi = index[v];
      char rgb[8];
      double width = .75;
      if (size_max != size_min)
	width = 0.1 + ((double) sizes[vi] - size_min) / (double) (size_max - size_min);
      
      double val = 0.5;
      node_weight[vi] = 1.0 / width;

      val = color[vi];
      
      get_rgb_value(rgb, val);

      int found = binarySearch(square_nodes,vi);
      
      if (found == -1)
	output<<vi<<"[label=\"\", width="<<width<<", shape=circle, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";
      else
	output<<vi<<"[label=\"\", width="<<width<<", shape=square, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";
     


    }
  
  if(directed)
    {
      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert u,v;
	  u = source(e,G);
	  v = target(e,G);

	  v_size_t ui = index[u];
	  v_size_t vi = index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }
  else
    {
      cout<<"Graph viz undirected\n";
      UGraph UG((v_Usize_t) n );

      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert s = source(e,G);
	  Vert t = target(e,G);

	  //	  cout<<"s "<<s<<" t "<<t<<"\n";

	  long ui = index[s];
	  long vi = index[t];

	  //cout<<"ui "<<ui<<" vi "<<vi<<"\n";

	  UVert u = vertex(ui,UG);
	  UVert v = vertex(vi,UG);

	  //cout<<"u "<<u<<" v "<<v<<"\n";

	  UEdge exist;
	  bool yes1,yes2;

	  tie(exist,yes1) = edge(u,v,UG);
	  tie(exist,yes2) = edge(v,u,UG);
	  if(!yes1 && !yes2)
	    add_edge(u,v,UG);
	}

      graph_traits<UGraph>::edge_iterator uit, uitend;
      property_map<UGraph, vertex_index_t>::type u_index = get(vertex_index,UG);

      for(tie(uit,uitend)=edges(UG);uit!=uitend;++uit)
	{
	  UEdge e = *uit;
	  UVert u,v;
	  u = source(e,UG);
	  v = target(e,UG);

	  long ui = u_index[u];
	  long vi = u_index[v];

	  long max_size = ui;
	  if (sizes[vi] > max_size)
	    max_size = vi;

	  output<<ui<<" -- "<<vi<<" [weight="<<node_weight[max_size]<<"];\n";
	}
    }

  output<<"}";
  output.close();
}

void write_scaled_labeled_graphviz(Graph & G, string output_file, vector<double> color, vector<size_t> sizes, bool directed)
{
  ofstream output;
  output.open(output_file.c_str());

  double color_min = color[0];
  double color_max = color[0];
  size_t size_min = sizes[0];
  size_t size_max = sizes[0];

  v_size_t n = num_vertices(G);

  if (color.size() != sizes.size() || color.size() != n)
    {
      cout<<"Color vector, size vector must all be of length num_nodes in graph. Exiting. \n";
      exit(EXIT_FAILURE);
    }

  for (size_t i = 0; i < color.size(); ++i)
    {
      if(color[i] > color_max)
	color_max = color[i];

      if(color[i] < color_min)
	color_min = color[i];

      if(sizes[i] > size_max)
	size_max = sizes[i];

      if(sizes[i] < size_min)
	size_min = sizes[i];
    }

  if (color_max > 1 or color_min < 0)
    {
      cout<<"Color vector must be in range [0,1]. Exiting.\n";
      exit(EXIT_FAILURE);
    }

  output<<"Graph color {\n";
  output<<"splines=false;\n";

  graph_traits<Graph>::edge_iterator eit, eitend;
  graph_traits<Graph>::vertex_iterator vit, vitend;

  property_map<Graph, vertex_index_t>::type index = get(vertex_index,G);

  for(tie(vit,vitend)=vertices(G);vit!=vitend;++vit)
    {
      Vert      v = *vit;
      v_size_t vi = index[v];
      char rgb[8];
      double width = .75;
      if (size_max != size_min)
	width = 0.1 + ((double) sizes[vi] - size_min) / (double) (size_max - size_min);
      
      double val = 0.5;

      val = color[vi];
      
      get_rgb_value(rgb, val);
      stringstream ss;

      string label;
      ss<<vi;
      ss>>label;

      output<<vi<<"[label=\""<<label<<"\", width="<<width<<", shape=circle, style=\"filled\", fillcolor=\""<<rgb<<"\"];\n";

    }
  
  if(directed)
    {
      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert u,v;
	  u = source(e,G);
	  v = target(e,G);

	  v_size_t ui = index[u];
	  v_size_t vi = index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }
  else
    {
      cout<<"Graph viz undirected\n";
      UGraph UG((v_Usize_t) n );

      for(tie(eit,eitend)=edges(G);eit!=eitend;++eit)
	{
	  Edge e = *eit;
	  Vert s = source(e,G);
	  Vert t = target(e,G);

	  //	  cout<<"s "<<s<<" t "<<t<<"\n";

	  long ui = index[s];
	  long vi = index[t];

	  //cout<<"ui "<<ui<<" vi "<<vi<<"\n";

	  UVert u = vertex(ui,UG);
	  UVert v = vertex(vi,UG);

	  //cout<<"u "<<u<<" v "<<v<<"\n";

	  UEdge exist;
	  bool yes1,yes2;

	  tie(exist,yes1) = edge(u,v,UG);
	  tie(exist,yes2) = edge(v,u,UG);
	  if(!yes1 && !yes2)
	    add_edge(u,v,UG);
	}

      graph_traits<UGraph>::edge_iterator uit, uitend;
      property_map<UGraph, vertex_index_t>::type u_index = get(vertex_index,UG);

      for(tie(uit,uitend)=edges(UG);uit!=uitend;++uit)
	{
	  UEdge e = *uit;
	  UVert u,v;
	  u = source(e,UG);
	  v = target(e,UG);

	  long ui = u_index[u];
	  long vi = u_index[v];

	  output<<ui<<" -- "<<vi<<";\n";
	}
    }

  output<<"}";
  output.close();
}

/*Calculates rgb heat map stored result in rgb, uses val between 0-1,
  0 is blue and 1 is red*/
void get_rgb_value(char rgb[], double val)
{
  int red, green, blue;
  double temp_color_value;

  //Red
  temp_color_value = 0;
  if(val > 0.875)
    temp_color_value = -4.0*val + 4.5;
  else if(val > .625)
    temp_color_value = 1.0;
  else if(val > .375)
    temp_color_value = 4.0*val - 1.5;

  red = 255*temp_color_value;

  //green
  temp_color_value = 0;
  if(val > 0.875)
    temp_color_value = 0;
  else if(val > .625)
    temp_color_value = -4.0*val + 3.5;
  else if(val > .375)
    temp_color_value = 1.0;
  else if(val>.125)
    temp_color_value = 4.0*val - 0.5;

  green = 255*temp_color_value;

  //blue
  temp_color_value = 0;
  if(val > 0.625)
    temp_color_value = 0;
  else if(val > .375)
    temp_color_value = -4.0*val + 2.5;
  else if(val > .125)
    temp_color_value = 1.0;
  else if(val>=0.0)
    temp_color_value = 4.0*val + .5;
  
  blue = 255*temp_color_value;

  sprintf(rgb,"#%02x%02x%02x",red,green,blue);
}
