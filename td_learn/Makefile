CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/tree_decomp_utils/loadTreeDecomposition.cpp ../lib/tree_decomp_utils/td_learn.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/connected.cpp ../lib/graph_utils/subset.cpp ../lib/graph_utils/BFS.cpp
objects = $(sources:.cpp=.o)

all: td_learn

td_learn: td_learn_main.cpp ../lib/graph_lib_boost.hpp ../lib/tree_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) td_learn_main.cpp -o td_learn $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f td_learn