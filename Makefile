#Makefile for testgraph

#Jaguar options
#CC = CC
#CFLAGS = -Ofast -Wall -fopenmp -march=native -m64
#LIB_DIR = -lrt 
#INC_DIR = -I ~/boost_1_48_0

#Smoky/Lens options
#CC = g++
#CFLAGS = -Ofast -Wall -fopenmp -march=native -m64
#LIB_DIR = -lrt 
#INC_DIR = -I ~/boost_1_48_0

#Standard Linux options
CC = g++
CFLAGS = -Ofast -Wall -fopenmp -m64 -g
LIB_DIR = -lrt 
INC_DIR = -I /usr/local/boost_1_48_0/

objects = approx_Page_Rank.o ncp_calc.o connected.o loadGraph.o subset.o k_core.o k_core_stat.o k_core_ncp.o k_shell_snapshot.o produce_plot.o produce_loglog_plot.o BFS_distance.o calc_gromov.o BFS.o get_shells.o write_dimacs.o collapse_core.o write_color_file.o

gromov_objects = connected.o loadGraph.o subset.o BFS_distance.o calc_gromov.o BFS.o write_dimacs.o write_color_file.o k_core.o

connected_graph_objects = connected.o subset.o BFS_distance.o BFS.o

snapshot_objects = loadGraph.o connected.o subset.o BFS_distance.o BFS.o k_core.o snapshot.o write_graphviz.o approx_Page_Rank.o

plot_objects = produce_plot.o produce_loglog_plot.o produce_3d_plot.o

td_objects = loadTreeDecomposition.o td_learn.o

all: testGraph ncp k_core_ncp k_shell_snapshot test_distance

k_means: k_means_main.cpp $(snapshot_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) k_means_main.cpp -o k_means $(snapshot_objects)

parse_amazon: parse_amazon_meta_main.cpp
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) parse_amazon_meta_main.cpp -o parse_amazon 

#these compile the specific executables 
distance_oracle: BFS_oracle.cpp graph_lib_boost.hpp loadGraph.o $(connected_graph_objects) k_core.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) BFS_oracle.cpp -o oracle loadGraph.o $(connected_graph_objects) k_core.o

label_prop: label_propagation_main.cpp graph_lib_boost.hpp loadGraph.o machine_learning.o $(connected_graph_objects) k_core.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) label_propagation_main.cpp -o label_prop loadGraph.o machine_learning.o $(connected_graph_objects) k_core.o

td_learn: td_learn_main.cpp graph_lib_boost.hpp tree_lib_boost.hpp $(td_objects) loadGraph.o $(connected_graph_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) td_learn_main.cpp -o td_learn $(td_objects) loadGraph.o $(connected_graph_objects)

tree_corr: tree_correlate_main.cpp graph_lib_boost.hpp loadGraph.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) tree_correlate_main.cpp -o tree_corr loadGraph.o

tree_stats: tree_decomp_bag_stats_main.cpp graph_lib_boost.hpp tree_lib_boost.hpp loadGraph.o write_graphviz.o $(connected_graph_objects) loadTreeDecomposition.o k_core.o write_color_file.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) tree_decomp_bag_stats_main.cpp -o tree_stats loadGraph.o write_graphviz.o $(connected_graph_objects) loadTreeDecomposition.o k_core.o write_color_file.o

community_in_tree: community_in_tree_main.cpp graph_lib_boost.hpp tree_lib_boost.hpp loadGraph.o $(connected_graph_objects) loadTreeDecomposition.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) community_in_tree_main.cpp -o community_in_tree loadGraph.o $(connected_graph_objects) loadTreeDecomposition.o

write_gviz: write_graphviz_main.cpp graph_lib_boost.hpp loadGraph.o write_graphviz.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_graphviz_main.cpp -o write_gviz loadGraph.o write_graphviz.o

testGraph: main.cpp graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) main.cpp -o testGraph $(objects)

k_core_stat: k_core_stat_main.cpp graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) k_core_stat_main.cpp -o k_stat $(objects)

dissect_core: dissect_core_main.cpp graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) dissect_core_main.cpp -o dissect_core $(objects)

graphSize: connected_size.cpp graph_lib_boost.hpp connected.o subset.o loadGraph.o BFS.o BFS_distance.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) connected_size.cpp -o graphSize connected.o subset.o loadGraph.o BFS.o BFS_distance.o

plot_gromov: plot_gromov_median_main.cpp graph_lib_boost.hpp produce_plot.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) plot_gromov_median_main.cpp -o plotGromovScale produce_plot.o

write_dimacs: write_dimacs_main.cpp graph_lib_boost.hpp connected.o loadGraph.o subset.o write_dimacs.o BFS.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_dimacs_main.cpp -o write_dimacs connected.o loadGraph.o subset.o write_dimacs.o BFS.o

ncp: ncp_main.cpp graph_lib_boost.hpp tree_lib_boost.hpp approx_Page_Rank.o $(connected_graph_objects) write_graphviz.o write_color_file.o k_core.o loadGraph.o ncp_calc.o produce_plot.o produce_loglog_plot.o TD_communities.o loadTreeDecomposition.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) ncp_main.cpp -o ncp approx_Page_Rank.o $(connected_graph_objects) write_graphviz.o write_color_file.o k_core.o loadGraph.o ncp_calc.o produce_plot.o produce_loglog_plot.o TD_communities.o loadTreeDecomposition.o

k_core_ncp: kcore_ncp_main.cpp graph_lib_boost.hpp $(objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) kcore_ncp_main.cpp -o k_core_ncp $(objects)

diffusion_label: diffusion_v_labelprop_main.cpp graph_lib_boost.hpp $(snapshot_objects) $(td_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) diffusion_v_labelprop_main.cpp -o diffusion_label $(snapshot_objects) $(td_objects)

snapshot: snapshot_main.cpp graph_lib_boost.hpp $(snapshot_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) snapshot_main.cpp -o k_shell_snapshot $(snapshot_objects)

test_distance: test_main.cpp graph_lib_boost.hpp $(gromov_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) test_main.cpp -o test_distance $(gromov_objects)

gromov: test_main.cpp graph_lib_boost.hpp $(gromov_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) test_main.cpp -o gromov $(gromov_objects)

write_k_core: write_k_core_main.cpp graph_lib_boost.hpp k_core.o loadGraph.o write_color_file.o connected.o subset.o BFS.o write_dimacs.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) write_k_core_main.cpp -o write_k_core k_core.o loadGraph.o write_color_file.o connected.o subset.o BFS.o write_dimacs.o

gen_rand_tree: gen_rand_tree_main.cpp graph_lib_boost.hpp write_dimacs.o 
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) gen_rand_tree_main.cpp -o gen_rand_tree write_dimacs.o

gen_rand_grid: gen_rand_grid_main.cpp graph_lib_boost.hpp write_dimacs.o 
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) gen_rand_grid_main.cpp -o gen_rand_grid write_dimacs.o

gen_rand_ring: gen_rand_ring_main.cpp graph_lib_boost.hpp write_dimacs.o 
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) gen_rand_ring_main.cpp -o gen_rand_ring write_dimacs.o

gen_rand_line: gen_rand_line_graph_main.cpp graph_lib_boost.hpp write_dimacs.o 
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) gen_rand_line_graph_main.cpp -o gen_rand_line write_dimacs.o

color_by_diffusion: color_by_diffusion_main.cpp graph_lib_boost.hpp loadGraph.o write_color_file.o connected.o subset.o BFS.o approx_Page_Rank.o
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) color_by_diffusion_main.cpp -o color_by_diffusion loadGraph.o write_color_file.o connected.o subset.o BFS.o approx_Page_Rank.o

plot_tw: plot_tw_data_main.cpp graph_lib_boost.hpp $(plot_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) plot_tw_data_main.cpp -o plot_tw $(plot_objects)

plot_3d_tw: plot_tw_data_3d_main.cpp graph_lib_boost.hpp $(plot_objects)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) plot_tw_data_3d_main.cpp -o plot_3d_tw $(plot_objects)
#This compiles our object files
.cpp.o: graph_lib_boost.hpp
	$(CC) $(CFLAGS) -c $(INC_DIR) $< -o $@

clean: 
	rm -f testGraph ncp k_core_ncp k_shell_snapshot test_distance graphSize plotGromovMedian *.o *~ 
############################
