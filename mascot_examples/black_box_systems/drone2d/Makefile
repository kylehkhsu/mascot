#
# compiler
#
CC        = g++
#CC       	  = clang++
#CXXFLAGS		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -fopenmp -Wfatal-errors -g
#CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
#CXXFLAGS 		= -Wall -Wextra -std=c++11 -Wfatal-errors -g -fopenmp
CXXFLAGS                = -Wall -Wextra -std=c++11 -g -Wfatal-errors

#
# scots 
#
SCOTSROOT		= ../../../
SCOTSINC		= -I$(SCOTSROOT)/utils -I$(SCOTSROOT)/bdd

#
# cudd 
#
CUDDPATH		=  /opt/local
#CUDDPATH		= ../../../../../libraries
CUDDINC 		= -I$(CUDDPATH)/include -I$(CUDDPATH)/include/cudd
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib

TARGET = drone

all: $(TARGET)

%.o:%.cc
	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o

cleanAll:
	rm -r -f ./$(TARGET)  ./$(TARGET).o C E G O Z plotting InterimDom T *.log Zs* D* validZs* 
