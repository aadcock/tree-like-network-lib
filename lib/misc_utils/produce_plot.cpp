/*
This funtion uses a pipe to gnuplot to produce plots for the
k_core_stats.  It also has more general functionality

by Aaron Adcock, PhD Candidate, Stanford University
Aug. 2011
 */

#include <string>
#include <stdio.h>
#include "../graph_lib_boost.hpp"

void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xtics)
{
  produce_plot(infile,outfilename,columns_to_plot,label,title,xlabel,ylabel,xtics,false);
}

void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xtics,bool suppressOutput)
{
  string script_end = outfilename;
  script_end.append(".plot");
  FILE *gp = fopen(script_end.c_str(),"w");
  
  string xticsGnuplot = "set xtics (";
  xticsGnuplot.append(xtics);
  xticsGnuplot.append(")\n");

  //string setxaxis = "set xzeroaxis lt -1\n set grid\n";

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
  fileGnuplot.append(" with lines lw 3 title '");
  fileGnuplot.append(ts2.str());
  fileGnuplot.append("'\n");

  // fprintf(gp, xticsGnuplot.c_str());
  // fprintf(gp, xlabelGnuplot.c_str());
  // fprintf(gp, ylabelGnuplot.c_str());
  // fprintf(gp, titleGnuplot.c_str());
  // fprintf(gp, fileGnuplot.c_str());

  if(xtics!="")
    fprintf(gp, xticsGnuplot.c_str());
  
  if(xlabel!="")
    fprintf(gp, xlabelGnuplot.c_str());

  if(ylabel!="")
    fprintf(gp, ylabelGnuplot.c_str());

  if(title!="")
    fprintf(gp, titleGnuplot.c_str());

  //fprintf(gp,setxaxis.c_str());
  //string log = "set log y\n";
  string key = "set key outside right\n";
  string termSize = "set terminal wxt size 1280,480\n";
  
  if(!suppressOutput)
    {
      //fprintf(gp, log.c_str());
      fprintf(gp, key.c_str());
      fprintf(gp, termSize.c_str());
      
      fprintf(gp, fileGnuplot.c_str());
            
      int count = 0;
      for(int i=1;i<columns_to_plot.size();i++)
	{
	  if(count<25)
	    {
	      string nFile = "replot '";
	      stringstream s1,s2;
	  
	      s1<<(columns_to_plot[i]);
	      s2<<(label[i]);
	      nFile.append(infile);
	      nFile.append("' using 1:");
	      nFile.append(s1.str());
	      nFile.append(" with lines lw 3 title '");
	      nFile.append(s2.str());
	      nFile.append("'\n");
	  
	      fprintf(gp, nFile.c_str());
	      count++;
	    }
	  else
	    {

	      string filepng = "set output \"";
	      stringstream ss;
	      ss<<i/25;
	      filepng.append(outfilename);
	      filepng.append(ss.str());
	      filepng.append(".png\"\n");

	      fprintf(gp, "set term png size 1000,480\n");
	      fprintf(gp, filepng.c_str());
  
	      fprintf(gp, "replot\n");

	      fprintf(gp, "set term x11\n");
	      //fprintf(gp, log.c_str());
	      fprintf(gp, key.c_str());
	      fprintf(gp, termSize.c_str());


	      string nFile = "plot '";
	      stringstream s1,s2;
	  
	      s1<<(columns_to_plot[i]);
	      s2<<(label[i]);
	      nFile.append(infile);
	      nFile.append("' using 1:");
	      nFile.append(s1.str());
	      nFile.append(" with lines lw 3 title '");
	      nFile.append(s2.str());
	      nFile.append("'\n");
	  
	      fprintf(gp, nFile.c_str());
	      count = 0;
	    }  
	}
    }

  
  //pclose(gp);  

  string filepng = "set output \"";
  stringstream ss;
  ss<<columns_to_plot.size()/25 + 1;
  filepng.append(outfilename);
  filepng.append(ss.str());
  filepng.append(".png\"\n");

  fprintf(gp, "set term png size 1000,480\n");
  fprintf(gp, filepng.c_str());
  
  fprintf(gp, "replot\n");      


  // fprintf(gp, fileGnuplot.c_str());


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
  // fprintf(gp, "replot\n");
  fclose(gp);
}



void produce_plot(string infile, string outfilename, int x_col, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xtics,bool suppressOutput)
{
  string script_end = outfilename;
  script_end.append(".plot");
  FILE *gp = fopen(script_end.c_str(),"w");
  
  string xticsGnuplot = "set xtics (";
  xticsGnuplot.append(xtics);
  xticsGnuplot.append(")\n");

  //string setxaxis = "set xzeroaxis lt -1\n set grid\n";

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
  stringstream xs;
  xs<<x_col;

  ts1<<(columns_to_plot[0]);
  ts2<<(label[0]);
  fileGnuplot.append(infile);
  fileGnuplot.append("' using ");
  fileGnuplot.append(xs.str());
  fileGnuplot.append(":");
  fileGnuplot.append(ts1.str());
  fileGnuplot.append(" with lines title '");
  fileGnuplot.append(ts2.str());
  fileGnuplot.append("'\n");

  // fprintf(gp, xticsGnuplot.c_str());
  // fprintf(gp, xlabelGnuplot.c_str());
  // fprintf(gp, ylabelGnuplot.c_str());
  // fprintf(gp, titleGnuplot.c_str());
  // fprintf(gp, fileGnuplot.c_str());

  if(xtics!="")
    fprintf(gp, xticsGnuplot.c_str());
  
  if(xlabel!="")
    fprintf(gp, xlabelGnuplot.c_str());

  if(ylabel!="")
    fprintf(gp, ylabelGnuplot.c_str());

  if(title!="")
    fprintf(gp, titleGnuplot.c_str());

  //fprintf(gp,setxaxis.c_str());
  //string log = "set log y\n";
  string key = "set key outside right\n";
  string termSize = "set terminal wxt size 1280,480\n";
  
  if(!suppressOutput)
    {
      //fprintf(gp, log.c_str());
      fprintf(gp, key.c_str());
      fprintf(gp, termSize.c_str());
      
      fprintf(gp, fileGnuplot.c_str());
            
      int count = 0;
      for(int i=1;i<columns_to_plot.size();i++)
	{
	  if(count<25)
	    {
	      string nFile = "replot '";
	      stringstream s1,s2;
	  
	      s1<<(columns_to_plot[i]);
	      s2<<(label[i]);
	      nFile.append(infile);
	      nFile.append("' using ");
	      nFile.append(xs.str());
	      nFile.append(":");
	      nFile.append(s1.str());
	      nFile.append(" with lines title '");
	      nFile.append(s2.str());
	      nFile.append("'\n");
	  
	      fprintf(gp, nFile.c_str());
	      count++;
	    }
	  else
	    {

	      string filepng = "set output \"";
	      stringstream ss;
	      ss<<i/25;
	      filepng.append(outfilename);
	      filepng.append(ss.str());
	      filepng.append(".png\"\n");

	      fprintf(gp, "set term png size 1000,480\n");
	      fprintf(gp, filepng.c_str());
  
	      fprintf(gp, "replot\n");

	      fprintf(gp, "set term x11\n");
	      //fprintf(gp, log.c_str());
	      fprintf(gp, key.c_str());
	      fprintf(gp, termSize.c_str());


	      string nFile = "plot '";
	      stringstream s1,s2;
	  
	      s1<<(columns_to_plot[i]);
	      s2<<(label[i]);
	      nFile.append(infile);
	      nFile.append("' using ");
	      nFile.append(xs.str());
	      nFile.append(":");
	      nFile.append(s1.str());
	      nFile.append(" with lines title '");
	      nFile.append(s2.str());
	      nFile.append("'\n");
	  
	      fprintf(gp, nFile.c_str());
	      count = 0;
	    }  
	}
    }

  
  //pclose(gp);  

  string filepng = "set output \"";
  stringstream ss;
  ss<<columns_to_plot.size()/25 + 1;
  filepng.append(outfilename);
  filepng.append(ss.str());
  filepng.append(".png\"\n");

  fprintf(gp, "set term png size 1000,480\n");
  fprintf(gp, filepng.c_str());
  
  fprintf(gp, "replot\n");      


  // fprintf(gp, fileGnuplot.c_str());


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
  // fprintf(gp, "replot\n");
  fclose(gp);
}



void produce_plot(string infile, string outfilename, vector<int> columns_to_plot, vector<string>label, string title, string xlabel, string ylabel, string xstart,string xtic,bool suppressOutput)
{
  string script_end = outfilename;
  script_end.append(".plot");
  FILE *gp = fopen(script_end.c_str(),"w");
  
  string xticsGnuplot = "set xtics ";
  xticsGnuplot.append(xstart);
  xticsGnuplot.append(",");
  xticsGnuplot.append(xtic);
  xticsGnuplot.append("\n");
  //  xticsGnuplot.append(xtics);
  //xticsGnuplot.append(")\n");

  string setxaxis = "set xzeroaxis lt -1\n set grid\n";

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

  // fprintf(gp, xticsGnuplot.c_str());
  // fprintf(gp, xlabelGnuplot.c_str());
  // fprintf(gp, ylabelGnuplot.c_str());
  // fprintf(gp, titleGnuplot.c_str());
  // fprintf(gp, fileGnuplot.c_str());

  if(xtic!="")
    fprintf(gp, xticsGnuplot.c_str());
  
  if(xlabel!="")
    fprintf(gp, xlabelGnuplot.c_str());

  if(ylabel!="")
    fprintf(gp, ylabelGnuplot.c_str());

  if(title!="")
    fprintf(gp, titleGnuplot.c_str());

  fprintf(gp,setxaxis.c_str());
  //string log = "set log y\n";
  string key = "set key outside right\n";
  string termSize = "set terminal wxt size 1280,480\n";
  
  if(!suppressOutput)
    {
      //fprintf(gp, log.c_str());
      fprintf(gp, key.c_str());
      fprintf(gp, termSize.c_str());
      
      fprintf(gp, fileGnuplot.c_str());
            
      int count = 0;
      for(int i=1;i<columns_to_plot.size();i++)
	{
	  if(count<25)
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
	      count++;
	    }
	  else
	    {

	      string filepng = "set output \"";
	      stringstream ss;
	      ss<<i/25;
	      filepng.append(outfilename);
	      filepng.append(ss.str());
	      filepng.append(".png\"\n");

	      fprintf(gp, "set term png size 1000,480\n");
	      fprintf(gp, filepng.c_str());
  
	      fprintf(gp, "replot\n");

	      fprintf(gp, "set term x11\n");
	      //fprintf(gp, log.c_str());
	      fprintf(gp, key.c_str());
	      fprintf(gp, termSize.c_str());


	      string nFile = "plot '";
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
	      count = 0;
	    }  
	}
    }

  
  //pclose(gp);  

  string filepng = "set output \"";
  stringstream ss;
  ss<<columns_to_plot.size()/25 + 1;
  filepng.append(outfilename);
  filepng.append(ss.str());
  filepng.append(".png\"\n");

  fprintf(gp, "set term png size 1000,480\n");
  fprintf(gp, filepng.c_str());
  
  fprintf(gp, "replot\n");      


  // fprintf(gp, fileGnuplot.c_str());


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
  // fprintf(gp, "replot\n");
  fclose(gp);
}

//   gp = popen("gnuplot -persist","w");
//   fprintf(gp, "set xtics (\"0\" 0, \"1-5\" 1, \"5-10\" 2, \"10-20\" 3, \">20\" 4)\n");
//   fprintf(gp,"set xlabel \"Jump Range\"\n");
//   fprintf(gp,"set ylabel \"Percent of Edges\"\n");

//   tempFile = "plot '";
//   ts1.flush();
//   ts2.flush();

//   ts1<<(cores_to_plot[0]+2);
//   ts2<<(cores_to_plot[0]);
//   tempFile.append(filein);
//   tempFile.append("' using 1:");
//   tempFile.append(ts1.str());
//   tempFile.append(" with lines title 'Core ");
//   tempFile.append(ts2.str());
//   tempFile.append("'\n");

//   fprintf(gp, "set title \"Out Edges\"\n");
//   fprintf(gp, tempFile.c_str());
  
      
//   for(int i=1;i<cores_to_plot.size();i++)
//     {
//       string nFile = "replot '";
//       stringstream s1,s2;
	  
//       s1<<(cores_to_plot[i]+2);
//       s2<<(cores_to_plot[i]);
//       nFile.append(fileout);
//       nFile.append("' using 1:");
//       nFile.append(s1.str());
//       nFile.append(" with lines title 'Core ");
//       nFile.append(s2.str());
//       nFile.append("'\n");
	  
//       fprintf(gp, nFile.c_str());
//     }

//   fprintf(gp, "set term png\n");
//   filepng = "set output \"";
//   filepng.append(prefix);
//   filepng.append("_jump_o.png\"\n");
//   fprintf(gp, filepng.c_str());
//   fprintf(gp, "replot\n");
//   fprintf(gp, "set term x11\n");
//   pclose(gp);
// }
