####################################################
#                                                  #
#  Haifeng Chen (haifengc at usc dot edu)          #
#  University of Southern California               #
#  August 24, 2014                                 #
#                                                  #
####################################################

CXX := g++
CXXFLAGS := -O3 -Wall -fmessage-length=0

SRCS := query_search.cc \
		hits_matrix.cc \
		maximal_path.cc \
        ./seqalign/global_alignment.cc \
        ./seqalign/local_alignment.cc \
        ./index/index.cc \
        ./index/k_longest_path.cc \
        ./index/k_mer_position.cc \
        ./index/k_mer_nearest.cc \
        ./evaluation/evaluation.cc \
        ./util/option.cc \
        ./util/evalue.cc \
        ./util/bio_util.cc \
        ./util/fasta_file.cc
        
OBJS := $(SRCS:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
aligner: aligner_main.o $(OBJS)
	 $(CXX) $(CXXFLAGS) -o $@ $^
                
build: ./index/build_main.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
local_align: ./seqalign/local_align_main.o  $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
evaluate: ./evaluation/evaluation_main.o  $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
all: aligner build local_align evaluate

clean:
	rm -rf aligner build local_align evaluate *.exe *.o ./index/*.o ./util/*.o ./seqalign/*.o ./evaluation/*.o


                                        
