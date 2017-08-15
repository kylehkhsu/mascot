#
# matlab root 
#
MATLABROOT	= /usr/lib/matlab-9.1/

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
CUDDPATH		=  /home/kylehsu/local
CUDDINC 		= -I$(CUDDPATH)/include
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib


TARGET 			= mexSymbolicSet

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) $(CUDDINC) $(CUDDLPATH) $(CUDDLIBS) $(TARGET).cc

