#!/usr/bin/env python

import numpy as np
import math

this_fxn = "datfile2_cfxn.py"

#make a c function to obtain one column of datfile;
outfile1 = "RTparam_hardcoded.c" #main code
outfile2 = "RTparam_hardcoded.h" #header file with function declarations
of1      = open(outfile1,'w')
of2      = open(outfile2,'w')

ss = "//Hard-coded parts of attenuation coefficents\n"
of1.write(ss)
of2.write(ss)
of2.write("\n")

##Header files;
"""
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
"""
ss = "#include <math.h>\n"
of1.write(ss)
ss = "#include <stdio.h>\n"
of1.write(ss)
ss = "#include <stdlib.h>\n"
of1.write(ss)
of1.write("\n")

##function needs to be like this;
"""
/***********************************************/
int Amn_fxn_L1(double *chebys,int ncol) {

  double Amn[Nl][Nc] = {
     {1.1,2.1,3.1},
     {4.1,5.1,6.1},
     {7.1,8.1,9.1},
     {10.1,11.1,12.1}
  };

  int M  = Nl;
  int N  = Nc;
  int i,j,r;

  for(i=0;i<M;i++) {
     r         = i+M*j;
     chebys[r] = Amn[i][ncol];
  }
}
/***********************************************/
"""


#filename eg ../datfiles/RT_Hparam_coeffsTn_region1.dat;
LH = ['L','H']
fdir="../datfiles/"
for lh in range(0,2):
   for regno in range(1,6):
      fname0   = "RT_"+LH[lh] +"param_coeffsTn_region"+str(regno)+".dat" 
      fname    = fdir+fname0

      #open file for reading;
      f2    = open(fname,'r')
      lines = f2.readlines()
      Nl    = len(lines)
      f2.close()
      #
      i     = 0
      vals  = lines[i].split()
      Nc    = len(vals)

      #write 1st lines of function;
      ss = "/***********************************************/\n"
      of1.write(ss)
      of2.write(ss)

      #specific to L/H and regno;
      funcname = "Amn_fxn_"+LH[lh]+str(regno)
      ss = "int "+funcname+"(double *chebys,int ncol) {\n"
      of1.write(ss)
      ss = "//Made by "+this_fxn+",\n"
      of1.write(ss)
      ss = "//from "+fname0+",\n"
      of1.write(ss)
      ss = "//which was made by RTparam.m"";\n"
      of1.write(ss)
      of1.write("\n")

      ##############################################################
      ##finish function declaration for header file;
      ss = "int "+funcname+"(double *,int);\n"
      of2.write(ss)
      ss = "/***********************************************/\n"
      of2.write(ss)
      of2.write("\n")
      ##############################################################

      
      ##############################################################
      ##1st line of matrix;
      ss  = "  double Amn["+str(Nl)+"]["+str(Nc)+"] = {\n"
      of1.write(ss)
      ss  = "     {"+vals[0]
      for j in range(1,Nc):
         ss  = ss+","+vals[j]

      ss  = ss+"},\n"
      of1.write(ss) 

      ## rest of the lines of matrix;
      for i in range(1,Nl):
         vals  = lines[i].split()
         ss  = "     {"+vals[0]
         for j in range(1,Nc):
            ss  = ss+","+vals[j]

         ss  = ss+"}"
         if i<Nl-1:
            ss  = ss+",\n"
         else:
            ss  = ss+"\n"

         of1.write(ss)

      ss  = "  };\n"
      of1.write(ss) 
      of1.write("\n") 
      ##############################################################

      ## rest of function;
      ss = "  int M  = "+str(Nl)+";\n"
      of1.write(ss) 
      ss = "  int N  = "+str(Nc)+";\n"
      of1.write(ss) 
      ss = "  int i,r;\n"
      of1.write(ss) 
      ss = "\n"
      of1.write(ss) 
      ss = "  for(i=0;i<M;i++) {\n"
      ## of1.write(ss) 
      ## ss = "     r         = i+M*ncol;\n"
      of1.write(ss) 
      ss = "     chebys[i] = Amn[i][ncol];\n"
      of1.write(ss) 
      ss = "  }\n"
      of1.write(ss) 
      ss = "}\n"
      of1.write(ss) 
      ss = "/***********************************************/\n"
      of1.write(ss) 
      of1.write("\n")

##finished - now close files
of1.close()
of2.close()
