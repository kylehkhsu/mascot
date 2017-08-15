#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -g -std=c++11 -O0 -DNDEBUG

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/bdd -I$(SCOTSROOT)/utils 

#
# cudd 
#
CUDDPATH		=  /home/kylehsu/local
CUDDINC 		= -I$(CUDDPATH)/include
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib

TARGET = simple

all: $(TARGET)

%.o:%.cc
	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm  ./$(TARGET)  ./$(TARGET).o
