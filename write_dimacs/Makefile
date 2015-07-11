

CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/misc_utils/write_dimacs.cpp ../lib/graph_utils/BFS.cpp
objects = $(sources:.cpp=.o)

all: write_dimacs

write_dimacs: write_dimacs_main.cpp ../lib/graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_dimacs_main.cpp -o write_dimacs $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f write_dimacs 