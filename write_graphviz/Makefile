

CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/misc_utils/write_graphviz.cpp ../lib/graph_utils/loadGraph.cpp 
objects = $(sources:.cpp=.o)

all: write_graphviz

write_graphviz: write_graphviz_main.cpp ../lib/graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_graphviz_main.cpp -o write_graphviz $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f write_graphviz 