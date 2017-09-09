#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CXXFLAGS 		= -Wall -Wextra -std=c++11 

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/utils -I$(SCOTSROOT)/bdd 

#
# cudd 
#
CUDDPATH		=  /opt/local
CUDDINC 		= -I$(CUDDPATH)/include
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib

TARGET = unicycle

all: $(TARGET)

%.o:%.cc
	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o
