#
# matlab root 
#
MATLABROOT	= /usr/lib/matlab-9.4/
#
# compiler
#
CC       	  = $(MATLABROOT)/bin/mex
CXXFLAGS 		= -g -Wall -Wextra -std=c++11 -O3 -Wfatal-errors 
#
# scots
#
SCOTSROOT   = ../..
SCOTSINC		= -I$(SCOTSROOT)/utils -I$(SCOTSROOT)/spmat
#
# cudd 
#
CUDDPATH		=  /opt/local/
CUDDINC_a 		= -I$(CUDDPATH)/include
CUDDINC_b		= -I$(CUDDPATH)/include/cudd
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib



all: sparse bdd

sparse: GridPoints mexStaticController mexTransitions mexGoal

bdd: mexSymbolicSet


mexStaticController: mexStaticController.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexStaticController.cc

mexTransitions: mexTransitions.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexTransitions.cc

GridPoints: GridPoints.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) GridPoints.cc

mexGoal: mexGoal.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexGoal.cc

mexSymbolicSet: mexSymbolicSet.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -DSCOTS_BDD -largeArrayDims $(SCOTSINC) $(CUDDINC_a) $(CUDDINC_b) $(CUDDLPATH) $(CUDDLIBS) mexSymbolicSet.cc
	
