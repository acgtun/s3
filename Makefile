CXX := g++
CXXFLAGS := -O3 -Wall -fmessage-length=0

SRCS := index.cc nearest_kmer.cc option.cc
        
OBJS := $(SRCS:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

makedb: build_main.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

all: makedb

clean:
	rm -rf makedb *.exe *.o
