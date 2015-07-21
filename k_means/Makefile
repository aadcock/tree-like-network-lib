CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/graph_utils/BFS.cpp ../lib/graph_utils/BFS_distance.cpp ../lib/k_core_utils/k_core.cpp ../lib/graph_utils/approx_Page_Rank.cpp ../lib/misc_utils/write_graphviz.cpp ../lib/graph_utils/snapshot.cpp
objects = $(sources:.cpp=.o)

all: k_means

k_means: k_means_main.cpp ../lib/graph_lib_boost.hpp ../lib/tree_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) k_means_main.cpp -o k_means $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f k_means