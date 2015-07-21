CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/k_core_utils/k_core.cpp ../lib/graph_utils/BFS.cpp ../lib/graph_utils/BFS_distance.cpp ../lib/graph_utils/machine_learning.cpp
objects = $(sources:.cpp=.o)

all: label_propagation

label_propagation: label_propagation_main.cpp ../lib/graph_lib_boost.hpp ../lib/tree_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) label_propagation_main.cpp -o label_propagation $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f label_propagation