

CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

sources = ../lib/graph_utils/connected.cpp ../lib/graph_utils/loadGraph.cpp ../lib/graph_utils/subset.cpp ../lib/graph_utils/BFS.cpp ../lib/misc_utils/write_color_file.cpp ../lib/graph_utils/approx_Page_Rank.cpp
objects = $(sources:.cpp=.o)

all: color_by_diffusion

color_by_diffusion: color_by_diffusion_main.cpp ../lib/graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) color_by_diffusion_main.cpp -o color_by_diffusion $(objects)

.cpp.o: ../lib/graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f write_k_core