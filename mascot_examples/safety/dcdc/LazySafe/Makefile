#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS		= -Wall -Wextra -std=c++11 -O3 -DNCEBUG -Wfatal-errors
#CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
#CXXFLAGS 		= -Wall -Wextra -std=c++11 -g

#
# scots 
#
SCOTSROOT		= ../../../../
SCOTSINC		= -I$(SCOTSROOT)/utils -I$(SCOTSROOT)/bdd 

#
# cudd 
#
CUDDPATH		=  /opt/local
CUDDINCa 		= -I$(CUDDPATH)/include
CUDDINCb		= -I$(CUDDPATH)/include/cudd
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib

TARGET = dcdc

all: $(TARGET)

%.o:%.cc
	$(CC) -c -g -O0  $(CXXFLAGS) $(CUDDINCa) $(CUDDINCb) $(SCOTSINC) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o

clean_all:
	rm -r -f ./$(TARGET)  ./$(TARGET).o C Z SafeSets *.txt *.fig plotting T *.log Zs* D* validZs* 

