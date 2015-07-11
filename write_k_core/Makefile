

CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/misc_utils/write_dimacs.cpp ../lib/graph_utils/BFS.cpp ../lib/k_core_utils/k_core.cpp ../lib/misc_utils/write_color_file.cpp
objects = $(sources:.cpp=.o)

all: write_k_core

write_k_core: write_k_core_main.cpp ../lib/graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_k_core_main.cpp -o write_k_core $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f write_k_core