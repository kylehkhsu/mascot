#
# matlab root 
#
MATLABROOT	= /usr/lib/matlab-9.4/

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
SCOTSINC		= -I$(SCOTSROOT)/bdd -I$(SCOTSROOT)/utils

#
# cudd 
#
CUDDPATH		=  /opt/local
CUDDINC 		= -I$(CUDDPATH)/include -I$(CUDDPATH)/include/cudd
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib


TARGET 			= mexSymbolicSet

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) $(CUDDINC) $(CUDDLPATH) $(CUDDLIBS) $(TARGET).cc

