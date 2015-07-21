CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/graph_utils/BFS.cpp ../lib/graph_utils/BFS_distance.cpp ../lib/k_core_utils/k_core.cpp ../lib/hyperbolicity_utils/calc_gromov.cpp ../lib/misc_utils/write_dimacs.cpp ../lib/graph_utils/approx_Page_Rank.cpp ../lib/k_core_utils/get_shells.cpp ../lib/k_core_utils/collapse_core.cpp ../lib/misc_utils/write_color_file.cpp 
objects = $(sources:.cpp=.o)

all: dissect_core

dissect_core: dissect_core_main.cpp ../lib/graph_lib_boost.hpp ../lib/tree_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) dissect_core_main.cpp -o dissect_core $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f dissect_core