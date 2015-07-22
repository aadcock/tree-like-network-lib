CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/tree_decomp_utils/loadTreeDecomposition.cpp ../lib/tree_decomp_utils/td_learn.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/connected.cpp ../lib/graph_utils/subset.cpp ../lib/graph_utils/BFS.cpp ../lib/graph_utils/BFS_distance.cpp ../lib/misc_utils/write_graphviz.cpp ../lib/k_core_utils/k_core.cpp ../lib/misc_utils/write_color_file.cpp
objects = $(sources:.cpp=.o)

all: tree_decomp_bag_stats

tree_decomp_bag_stats: tree_decomp_bag_stats_main.cpp ../lib/graph_lib_boost.hpp ../lib/tree_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) tree_decomp_bag_stats_main.cpp -o tree_decomp_bag_stats $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f tree_decomp_bag_stats