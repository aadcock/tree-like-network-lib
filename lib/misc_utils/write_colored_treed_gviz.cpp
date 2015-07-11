/* 
   The original functions in this file were provided by Blair Sullivan
   from her tree-decomposition code.  I have modified them to color
   the labels based on the k_core (or any coloring) of the
   node. Everything else is her original code.

   Modified By Aaron Adcock, PhD Candidate, Stanford University
*/





/*
 * Writes the graph to the provided dimacs file.
 */
void GraphVizGraphWriter::write_graph(Graph *g)
{
  int i, j;
  FILE *out;
  string graphviz_file = out_file_name;
  if ((out = fopen(graphviz_file.c_str(), "w")) == NULL)
    fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, graphviz_file.c_str());

  fprintf(out, "Graph G{\n");
  fprintf(out, "overlap=false;\n");
  list<int>::iterator it;
  vector<Node> nodes = g->get_nodes();
  list<int> nbrs;

  for (i = 0; i < g->get_num_nodes(); i++)
    {
      nbrs = nodes[i].get_nbrs();
      it = nbrs.begin();
      while (it != nbrs.end())
	{
	  j = *it;
	  if (i < j)
	    // Make this 1-based in true DIMACS spirit
	    fprintf(out, "%d -- %d;\n", i + 1, j + 1);
	  ++it;
	}
    }
  fprintf(out, "}\n");
  fflush(out);
  fclose(out);

  return;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////
//Modified
//////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * Function to write a graphviz representation (DOT language) of the
 * tree decomposition to the specified file.  The spline parameter
 * allows the user to specify whether or not spline curves are used in
 * edge layout.  This is (strongly) not recommended for large graphs
 * (> 100 vertices).  The final parameter takes one of the options
 * GV_BAG_LABELS, GV_SCALE_BAGS, GV_TREE_ONLY GV_BAG_LABELS lists the
 * entire bag of vertices inside each node's label.  GV_SCALE_BAGS
 * scales the bags to indicate bag size, and labels them with the
 * treenode ID. Note this may create dimacs files that neato cannot
 * handle on very large graphs!  GV_TREE_ONLY writes a file that just
 * draws the underlying tree.  GV_COLORS creates a colored version of
 * the decomposition showing the underlying tree with nodes colored
 * red if they have size k+1, blue if size k, and gradations of purple
 * for smaller nodes, getting lighter as you get smaller.
 **/
void TDTree::write_graphviz_file(bool spline, const char *GVIZ_file, const vector<double>& color_vector, const double max_color, const double min_color)
{
  int i, j; 
  FILE *out; 
  Graph::WeightedMutableGraph *G = this->G;
  list<int>::iterator L;
  char rgb[8];
  int k=0, rcolor; 

  if( (out = fopen(GVIZ_file, "w")) == NULL)
    fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

  fprintf(out, "Graph TD{\n"); 
  /* You can choose one of the two options below if you want to eliminate node overlap 
   * (the scale option makes a larger, more symmetric drawing than the false option). 
   * The problem is that the scale factors for these graphs are enormous: output files are in the range of 3MB (png) for a 50 node tree, 
   * and PDF rendering fails due to canvas size restrictions. 
   * Current config allows node overlap, which also makes diagrams difficult to read. We should think on how to fix this.
   */
	

  fprintf(out, "overlap=false;\n");
  if(spline==true)
    fprintf(out, "splines=true;\n");

  print_message(1, "k : %d\n", k);

  // Print out the bags
  // Feb 24 2011 - changed back to tree_nodes.size() and check for NULL.
  int size=this->tree_nodes.size();
  Graph::Node *n1;
  double val;
  for(i=0;i<size;i++)
    {
      fprintf(out,"%d [label=<",i+1);
      print_message(10, "Bag %d:\n", i+1);
      // Sort the bags before outputting
      this->tree_nodes[i]->bag.sort();
      for(L=this->tree_nodes[i]->bag.begin();
	  L!=this->tree_nodes[i]->bag.end();++L )
	{
	  // Get the actual node label if necessary (we have a pointer to G)
	  j=*L;
	  print_message(10, "Looking for %d\n", j);
	  n1=this->G->get_node(j);
	  j=n1->get_label();
	  //j=G->nodes[j].label;
	  print_message(10, "Found %d\n", j);
		  
	  if(max_color==min_color)
	    val = .5;
	  else
	    val = (color_vector[j]-min_color)/(max_color-min_color);
		  
	  get_rgb_value(rgb,val);
	  fprintf(out,"<FONT COLOR=\"%s\">",rgb);
	  fprintf(out,"%d ",j);
	  fprintf(out,"</FONT>\n");
	  // NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
	  // are disconnected nodes, then the union of all the bags
	  // will NOT be equal to (1,2,...,num_total_nodes)
	}
      fprintf(out,">];\n");

    }


  // Print out the edges
  //Feb 24 2011 - updated to allow for NULL entries in tree_nodes array
  for(i=0;i<size;i++)
    {   
      if(this->tree_nodes[i] != NULL)
	{
	  for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
	    {
	      j=*L;
	      if(i < j)
		// Make this 1-based in true DIMACS spirit
		fprintf(out,"%d -- %d;\n",i+1,j+1);
	    }
	}
    }

  fprintf(out, "}\n");

  fflush(out);  
  fclose(out); 
}

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
  else if(val>0.0)
    temp_color_value = 4.0*val + .5;
  
  blue = 255*temp_color_value;

  sprintf(rgb,"#%02x%02x%02x",red,green,blue);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//End Modified
/////////////////////////////////////////////////////////////////////////////////////////////////////////








/*
 * Highlights the subtree of the tree decomposition which contains the vertex v. Writes a
 * gviz file to the specified output file. 
 * The final parameter takes one of the options GV_BAG_LABELS, GV_SCALE_BAGS, GV_TREE_ONLY
 * GV_BAG_LABELS lists the entire bag of vertices inside each node's label. 
 * GV_SCALE_BAGS scales the bags to indicate bag size, and labels them with the treenode ID. Note this may 
 * create dimacs files that neato cannot handle on very large graphs! 
 * GV_TREE_ONLY creates small circles for all nodes in the tree decomposition.
 */
void TDTree::highlight_subtree_gviz(int v, const char *GVIZ_file, int style)
{
  int i, j; 
  FILE *out; 
  Graph::WeightedMutableGraph *G = this->G;
  list<int>::iterator L;
  char rgb[8];
  //You can modify this to make the uniform circle bags bigger or smaller.
  double fixed_width = .5;
  //for now we'll make the highlights a light blue                
  sprintf(rgb, "#00ced1");

  //this should be a valid test for whether or not the locations have been reocrded.
  //it does not guarantee they are up to date, necessarily.
  if(this->node_locations[0].size() == 0)
    {
      this->record_node_locations();
    }
  int size=(int) this->tree_nodes.size();
  vector<bool> containsv(size,false);
  for(L = (this->node_locations[v]).begin(); L != (this->node_locations[v]).end(); ++L)
    {
      //	    printf("Node %d is in bag %d \n", v, *L);
      containsv[*L] = true;
    }

  if( (out = fopen(GVIZ_file, "w")) == NULL)
    fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

  fprintf(out, "Graph TD{\n"); 
  //fprintf(out, "overlap=false;\n");
  Graph::Node *n1;
  for(i=0;i<size;i++)
    {
      if(this->tree_nodes[i] != NULL)
	{
	  if(style == GV_BAG_LABELS)
	    {
	      if(containsv[i])
		fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", label=\"",i+1, rgb);
	      else
		fprintf(out,"%d [label=\"",i+1);

	      print_message(10, "Bag %d:\n", i+1);
	      // Sort the bags before outputting
	      this->tree_nodes[i]->bag.sort();
	      for(L=this->tree_nodes[i]->bag.begin(); L!=this->tree_nodes[i]->bag.end();++L )
		{
		  // Get the actual node label if necessary (we have a pointer to G)
		  j=*L;
		  print_message(10, "Looking for %d\n", j);
		  n1=this->G->get_node(j);
		  j=n1->get_label();
		  //j=G->nodes[j].label;
		  print_message(10, "Found %d\n", j);
		  fprintf(out,"%d ",j);
		  // NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
		  // are disconnected nodes, then the union of all the bags
		  // will NOT be equal to (1,2,...,num_total_nodes)
		}
	      fprintf(out,"\"];\n");
	    }
	  else if(style == GV_SCALE_BAGS)
	    {
	      if(containsv[i])
		fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", width=%f, shape=ellipse];\n",i+1, rgb,
			(double)(this->tree_nodes[i]->bag.size())/sqrt((double)G->get_capacity()));
	      else
		fprintf(out,"%d [width=%f, shape=ellipse];\n",i+1, 
			(double)(this->tree_nodes[i]->bag.size())/sqrt((double)G->get_capacity()));
	    }
	  else if(style == GV_TREE_ONLY)
	    {
	      if(containsv[i])
		fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", width=%f, shape=ellipse];\n",i+1, rgb, fixed_width);
	      else
		fprintf(out,"%d [width=%f, shape=ellipse];\n",i+1, fixed_width);


	    }
	  else
	    fatal_error("%s: style parameter was not one of GV_BAG_LABELS, GV_TREE_ONLY, or GV_SCALE_BAGS.\n", __FUNCTION__);

	}
    }	   

  // Print out the edges
  for(i=0;i<size;i++)
    {   
      if(this->tree_nodes[i] != NULL)
	{
	  for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
	    {
	      j=*L;
	      // Make this 1-based in true DIMACS spirit
	      if(i < j)
		fprintf(out,"%d -- %d;\n",i+1,j+1);
	    }
	}
    }

  fprintf(out, "}\n");

  fflush(out);  
  fclose(out); 
}
