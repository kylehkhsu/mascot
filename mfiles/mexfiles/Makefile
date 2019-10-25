#
# matlab root 
#
MATLABROOT	= /Applications/MATLAB_R2017b.app/

#
# compiler
#
#CC        = g++
CC       	  = $(MATLABROOT)/bin/mex
CXXFLAGS 		= -g -Wall -Wextra -std=c++11 -O3 

#
# scots
#
SCOTSROOT   = ../..
SCOTSINC		= -I$(SCOTSROOT)/bdd

#
# cudd 
#
CUDDPATH		=  /opt/local
CUDDINCa 		= -I$(CUDDPATH)/include
CUDDINCb		= -I$(CUDDPATH)/include/cudd
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib


TARGET 			= mexSymbolicSet

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) $(CUDDINCa) $(CUDDINCb) $(CUDDLPATH) $(CUDDLIBS) $(TARGET).cc

