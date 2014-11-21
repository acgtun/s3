####################################################
#                                                  #
#  RMP (Reads Mapping Protein Database)            #
#  Haifeng Chen (haifengc at usc dot edu)          #
#  University of Southern California               #
#  April 4, 2014                                   #
#                                                  #
####################################################


CC := gcc
CXX := g++
CXXFLAGS := -O3 -Wall -fmessage-length=0

TARGET := aligner

SRCS := iofile.cpp option.cpp RMP.cpp orf.cpp KlongestPathDP.cpp buildindex.cpp dbsearch.cpp

OBJS := $(SRCS:.cpp=.o)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	

                
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

all: $(TARGET)

clean:
	rm -rf $(TARGET) *.o *.exe
                                        
