####################################################
#                                                  #
#  Haifeng Chen (haifengc at usc dot edu)          #
#  University of Southern California               #
#  Jan 2015                                        #
#                                                  #
####################################################

CXX := g++
CXXFLAGS := -O3 -Wall -fmessage-length=0

SRCS := index.cc nearest_kmer.cc option.cc query_search.cc evalue.cc local_alignment.cc evaluation.cc
        
OBJS := $(SRCS:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

makedb: build_main.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
s3: aligner_main.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
evaluate: evaluate_main.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

all: makedb s3 evaluate

clean:
	rm -rf makedb s3 evaluate *.exe *.o
