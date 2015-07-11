/*
This file writes out a color file for use in conjunction with Blair
Sullivan's tree decomposition code (INDDGO).  Using the correct
compiler flags with the executable produced by her code, the tree
decomposition can be visualized with labels colored according to the
colors produced by this file.

NOTE: This function really just prints out a line of numbers with some
fancy formatting.

NOTE: This assumes the colors are passed to the writer in ascending
nodeID order.

By Aaron Adcock, PhD Candidate, Stanford University, 2012

INPUT: 

string output_file Name of the output color file 
string input_file Name of the input network (for color file comments) 
vector<double> color Vector, element i is the color for the ith node

 */

#include "../graph_lib_boost.hpp"

void write_color_file(string output_file, string input_file, vector<double> color)
{

  ofstream out_file;
  out_file.open(output_file.c_str());
  
  out_file<<"#This file contains the color generated from the file "<<input_file<<" the colors are provided in ascending nodeID order\n";
  
  
  double color_max = color[0];
  double color_min = color[0];
  double color_avg = 0;

  for(int i = 0; i < color.size(); ++i)
    {
      if(color[i] > color_max)
	color_max = color[i];

      if(color[i] < color_min)
	color_min = color[i];

      color_avg += color[i];
    }

  color_avg = color_avg / color.size();
  out_file<<"#min "<<color_min<<"\n";
  out_file<<"#max "<<color_max<<"\n";
  out_file<<"#avg "<<color_avg<<"\n";

  for(int i = 0; i < color.size(); ++i)
    out_file<< color[i] <<" ";

  out_file<<"\n";

  out_file.close();

}

/*
This version of the function writes integer colors to file.

INPUT: 

string output_file Name of the output color file 
string input_file Name of the input network (for color file comments) 
vector<int> color Vector, element i is the color for the ith node

*/

void write_color_file(string output_file, string input_file, vector<int> color)
{

  ofstream out_file;
  out_file.open(output_file.c_str());
  
  out_file<<"#This file contains the color generated from the file "<<input_file<<" the colors are provided in ascending nodeID order\n";
  
  
  int color_max = color[0];
  int color_min = color[0];
  double color_avg = 0;

  for(int i=0;i<color.size();++i)
    {
      if(color[i] > color_max)
	color_max = color[i];

      if(color[i] < color_min)
	color_min = color[i];

      color_avg += color[i];
    }

  color_avg = color_avg/color.size();
  out_file<<"#min "<<color_min<<"\n";
  out_file<<"#max "<<color_max<<"\n";
  out_file<<"#avg "<<color_avg<<"\n";

  for(int i=0;i<color.size();++i)
    out_file<<color[i]<<" ";

  out_file<<"\n";

  out_file.close();

}
