/*
Main file for parsing amazon-meta.txt obtained from the Stanford
Network Analysis Platform website: snap.stanford.edu.  Currently this
only parses out products which are books.

By Aaron Adcock
Stanford University, Ph.D. (Obtained Aug. 2014!)
Aug. 2014

Executable options:

-i input file...REQUIRED
-o output file prefix

The output is several files: 

1) output_graph.txt, a file containing an
edge list of product-product edges 

2) output_asin.txt, a file listing the ASIN's for each id in
output_graph.txt

3) output_product_categories.txt, a file listing the leaf category
numbers for each product.

4) output_category_tree.txt, a file listing the children of each
category.

5) output_category_name.txt, a file listing the name of each category.
Unclear if it will be full category tree...


*/

#include "graph_lib_boost.hpp"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <map>
#include <vector>
#include <fstream>
#include <queue>

vector<string> split(string & line);
vector<string> split(string & line, string chars_to_split);

int main(int argc, char * argv[])
{
  string input_file = "";
  string output_prefix = "output";

  if (argc < 2)
    {
      cout<<"-i <input_file> \n-o <output_prefix> \n";
      return (-1);
    }

  for (int i = 1; i < argc; i++)
    {
      if (strcmp(argv[i], "-i") == 0)
	input_file = argv[i + 1];
      else if (strcmp(argv[i], "-o") == 0)
	output_prefix = argv[i + 1];
    }

  if (input_file.length() == 0)
    {
      cerr<<"You must provide an input metadate file.\n";
      exit(EXIT_FAILURE);
    }

  ifstream meta_file;
  meta_file.open(input_file.c_str());

  map<string, int> asin_to_id;
  map<string, string> asin_to_title;
  map<int, vector<string> > category_to_books;
  map<string, vector<int> > book_to_categories;
  vector< pair<string, string> > edge_list;
  map<int, vector<int> > category_children;
  map<int, string> category_num_to_string;

  const int cat_root = 283155;

  if (meta_file.is_open())
    {
      string line;
      int curr_id = -1;
      string ASIN;
      string title;
       string group;

      while (meta_file.good())
	{
	  getline(meta_file, line);
	  stringstream ss(line);
	  
	  string sub;
	  ss>>sub;
	  if (strcmp(sub.c_str(), "Id:") == 0)
	    curr_id++;
	  else if (strcmp(sub.c_str(), "ASIN:") == 0) 
	    ss>>ASIN; 
	  else if (strcmp(sub.c_str(), "title:") == 0)
	    ss>>title;
	  else if (strcmp(sub.c_str(), "group:") == 0)
	    {
	      ss>>group;
	      if (curr_id >= 0 && strcmp(group.c_str(), "Book") == 0)
		{
		  asin_to_title[ASIN] = title;
		  asin_to_id[ASIN] = curr_id;
		}
	    }
	  else if (strcmp(sub.c_str(), "similar:") == 0)
	    {
	      size_t num_edges;
	      ss >> num_edges;
	      while (ss.good())
		{
		  string other_ASIN;
		  ss >> other_ASIN;
		  edge_list.push_back(make_pair(ASIN, other_ASIN));
		}
	    }
	  else if (strcmp(sub.c_str(), "categories:") == 0 && strcmp(group.c_str(), "Book") == 0)
	    {
	      size_t num_cat;
	      ss >> num_cat;

	      for (size_t i = 0; i < num_cat; ++i)
		{
		  getline(meta_file, line);
		  
		  // Now go through each line and pick out each category
		  vector<string> tokens = split(line, string(" |"));
		  int prev_num = cat_root;
		  for (size_t j = 0; j < tokens.size(); ++j)
		    {
		      vector<string> cat_num = split(tokens[j], string("[]"));
		      
		      if (cat_num.size() == 2)
			{
			  string category = cat_num[0];
			  stringstream ss(cat_num[1]);
			  int num;
			  ss>>num;
			  
			  // Set up category to books map
			  if (category_to_books.find(num) != category_to_books.end())
			    category_to_books[num].push_back(ASIN);
			  else
			    {
			      vector<string> dummy;
			      dummy.push_back(ASIN);
			      category_to_books[num] = dummy;
			    }

			  // Set up book to categories map
			  if (book_to_categories.find(ASIN) != book_to_categories.end())
			    book_to_categories[ASIN].push_back(num);
			  else
			    {
			      vector<int> dummy;
			      dummy.push_back(num);
			      book_to_categories[ASIN] = dummy;
			    }
			  
			  // Set up category num to string map
			  category_num_to_string[num] = category;

			  // Make category tree
			  if (num != cat_root)
			    {
			      if (category_children.find(prev_num) != category_children.end())
				category_children[prev_num].push_back(num);
			      else
				{
				  vector<int> dummy;
				  dummy.push_back(num);
				  category_children[prev_num] = dummy;
				}
			    }

			  prev_num = num;
			}
		      else
			cerr<<"Category did not match CAT[NUM] pattern.\n";
		    }
		}
	    }
	}
    }
  else
    {
      cerr<<"Error opening file.\n";
      exit(EXIT_FAILURE);
    }

  // Need to fix ASIN list with id list / write out files
  
  ofstream ofile;

  string network_file = output_prefix;
  network_file.append("_edgelist.txt");
  ofile.open(network_file.c_str());
  
  ofile<<"# A network file generated from "<<input_file<<", which is presumed to contain amazon meta data.\n";

  for (size_t i = 0; i < edge_list.size(); ++i)
    {
      pair<string, string> edge = edge_list[i];
      
      ofile<<asin_to_id[edge.first]<<" "<<asin_to_id[edge.second]<<"\n";
    }
  ofile.close();

  string asin_id_file = output_prefix;
  asin_id_file.append("_asin.txt");

  ofile.open(asin_id_file.c_str());
  ofile<<"# ASIN's associated with each product id.\n";
  
  for (map<string, int>::iterator it = asin_to_id.begin(); it != asin_to_id.end(); ++it)
    ofile<<it->second<<" "<<it->first<<"\n";

  ofile.close();

  string asin_title_file = output_prefix;
  asin_title_file.append("_titles.txt");

  ofile.open(asin_title_file.c_str());
  ofile<<"# Title of product for each id.\n";

  for (map<string, string>::iterator it = asin_to_title.begin(); it != asin_to_title.end(); ++it)
    ofile<<asin_to_id[it->first]<<" "<<it->second<<"\n";

  ofile.close();

  string category_books_file = output_prefix;
  category_books_file.append("_categories.txt");
  
  ofile.open(category_books_file.c_str());
  ofile<<"# Gives title of category followed by the list of books in each category.\n";

  for (map<int, vector<string> >::iterator it = category_to_books.begin(); it != category_to_books.end(); ++it)
    {
      ofile<<it->first<<" ";
      for (size_t i = 0; i < it->second.size(); ++i)
	ofile<<asin_to_id[it->second[i]]<<" ";
      
      ofile<<"\n";
    }

  ofile.close();

  string book_categories_file = output_prefix;
  book_categories_file.append("_book_categories.txt");

  ofile.open(book_categories_file.c_str());
  ofile<<"# This file contains a book id followed the categories the book is in.\n";

  for (map<string, vector<int> >::iterator it = book_to_categories.begin(); it != book_to_categories.end(); ++it)
    {
      ofile<<asin_to_id[it->first]<<" ";
      for (size_t i = 0; i < it->second.size(); ++it)
	ofile<<it->second[i]<<" ";

      ofile<<"\n";
    }

  ofile.close();

  string category_tree_file = output_prefix;
  category_tree_file.append("_category_tree.txt");
  
  ofile.open(category_tree_file.c_str());
  ofile<<"# Contains the category tree: parent <child list>, in the appropriate order.\n";

  queue<int> Q;
  Q.push(cat_root);

  while (!Q.empty())
    {
      int parent = Q.front();
      vector<int> children = category_children[parent];

      ofile<<parent<<" ";
      for (size_t i = 0; i < children.size(); ++i)
	{
	  ofile<<children[i]<<" ";
	  Q.push(children[i]);
	}
      ofile<<"\n";
      Q.pop();
    }
  
  ofile.close();

  string category_num_string_file = output_prefix;
  category_num_string_file.append("_category_names.txt");

  ofile.open(category_num_string_file.c_str());
  ofile<<"# File contains the names of each category.\n";

  for (map<int, string>::iterator it = category_num_to_string.begin(); it != category_num_to_string.end(); ++it)
    ofile<<it->first<<" "<<it->second<<"\n";
    
  ofile.close();

  return 0;
}

vector<string> split(string & line)
{
  string empty = "";
  return split(line, empty);
}

vector<string> split(string & line, string chars_to_split)
{
  // First get length of string
  size_t line_length = line.length();
  vector<string> tokens;

  if (line_length == 0)
    {
      cerr<<"String passed of zero length.\n";
      return tokens;
    }

  // If chars_to_split is empty, use whitespace
  if (chars_to_split.length() == 0)
    chars_to_split = " \t\n";

  // Cycle through characters
  size_t word_start = 0;
  for (size_t i = 0; i < line_length; ++i)
    {
      bool flag_split_char_found = false;
      for (size_t j = 0; j < chars_to_split.length(); ++j)
	if (line[i] == chars_to_split[j])
	  flag_split_char_found = true;

      if (flag_split_char_found)
	{
	  if (i > word_start)
	    tokens.push_back(line.substr(word_start, i - word_start));

	  word_start = i + 1;
	}
    }
  
  if (word_start != line_length)
    tokens.push_back(line.substr(word_start, line_length - word_start));

  return tokens;
}
