/*
  This funtion uses a pipe to gnuplot to produce loglog plots.

  by Aaron Adcock, PhD Candidate, Stanford University
  Aug. 2011
*/

#include <string>
#include <stdio.h>
#include "../graph_lib_boost.hpp"

//using namespace std;

void produce_loglog_plot(string infile, string outfilename, string outplotname, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xtic)
{
  produce_loglog_plot(infile,outfilename,outplotname,columns_to_plot,label,title,xlabel,ylabel,xtic,true);
}
void produce_loglog_plot(string infile, string outfilename, string outplotname, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xtics, bool suppressOutput) 
{

  string script_end = outfilename;
  script_end.append(".plot");
  FILE *gp = fopen(script_end.c_str(),"w");
    
  string xticsGnuplot = "set xtics (";
  xticsGnuplot.append(xtics);
  xticsGnuplot.append(")\n");

  string xlabelGnuplot = "set xlabel \"";
  xlabelGnuplot.append(xlabel);
  xlabelGnuplot.append("\"\n");

  string ylabelGnuplot = "set ylabel \"";
  ylabelGnuplot.append(ylabel);
  ylabelGnuplot.append("\"\n");

  string titleGnuplot = "set title \"";
  titleGnuplot.append(title);
  titleGnuplot.append("\"\n");

  string fileGnuplot = "plot '";
  stringstream ts1,ts2;

  ts1<<(columns_to_plot[0]);
  ts2<<(label[0]);
  fileGnuplot.append(infile);
  fileGnuplot.append("' using 1:");
  fileGnuplot.append(ts1.str());
  fileGnuplot.append(" with lines title '");
  fileGnuplot.append(ts2.str());
  fileGnuplot.append("'\n");

 
  //std::cout<<xlabelGnuplot<<"\n";
  if(xtics!="")
    {
      fprintf(gp, xticsGnuplot.c_str());
    }  

  if(xlabel!="")
    {    
      fprintf(gp, xlabelGnuplot.c_str());
    }

  if(ylabel!="")
    {    
      fprintf(gp, ylabelGnuplot.c_str());
    }
  
  if(title!="")
    {
      fprintf(gp, titleGnuplot.c_str());
    }
  //cout<<"here?\n";
  fprintf(gp, "set log y\n");
  fprintf(gp, "set log x\n");
  fprintf(gp, "set key outside right\n");
  fprintf(gp, "set terminal wxt size 1280,480\n");

  if(1)
    {
      fprintf(gp, fileGnuplot.c_str());
      
      for(int i=1;i<columns_to_plot.size();i++)
	{
	  string nFile = "replot '";
	  stringstream s1,s2;
	  
	  s1<<(columns_to_plot[i]);
	  s2<<(label[i]);
	  nFile.append(infile);
	  nFile.append("' using 1:");
	  nFile.append(s1.str());
	  nFile.append(" with lines title '");
	  nFile.append(s2.str());
	  nFile.append("'\n");
	  
	  fprintf(gp, nFile.c_str());
	}
    }
 
  string filepng = "set output \"";
  filepng.append(outplotname);
  filepng.append(".png\"\n");
  fprintf(gp, "set term png size 1280,480\n");
  fprintf(gp, filepng.c_str());

  //fprintf(gp, fileGnuplot.c_str());
      
  // for(int i=1;i<columns_to_plot.size();i++)
  //   {
  //     string nFile = "replot '";
  //     stringstream s1,s2;
	  
  //     s1<<(columns_to_plot[i]);
  //     s2<<(label[i]);
  //     nFile.append(infile);
  //     nFile.append("' using 1:");
  //     nFile.append(s1.str());
  //     nFile.append(" with lines title '");
  //     nFile.append(s2.str());
  //     nFile.append("'\n");
      
  //     fprintf(gp, nFile.c_str());
  //   }

  fprintf(gp, "replot\n");
  fclose(gp);

}
