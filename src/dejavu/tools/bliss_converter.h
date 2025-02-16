#ifndef DEJAVU_BLISS_CONVERTER_H
#define DEJAVU_BLISS_CONVERTER_H

// requires "bliss/graph.hh" and "dejavu/graph.h" to be loaded

static void convert_dejavu_to_bliss(dejavu::sgraph* dejavu_graph, int* dejavu_col, bliss::Graph* bliss_graph) {
    std::vector<int> sassy_v_to_bliss_v;
    sassy_v_to_bliss_v.reserve(dejavu_graph->v_size);
    for(int i = 0; i < dejavu_graph->v_size; ++i) {
        sassy_v_to_bliss_v.push_back(bliss_graph->add_vertex(dejavu_col[i]));
    }
    for(int i = 0; i < dejavu_graph->v_size; ++i) {
        const int ept = dejavu_graph->v[i];
        const int vd  = dejavu_graph->d[i];
        for(int j = ept; j < ept + vd; ++j) {
            const int to = dejavu_graph->e[j];
            if(i < to) {
                bliss_graph->add_edge(i, to);
            }
        }
    }
}

static void convert_dejavu_to_bliss(dejavu::static_graph* dejavu_graph, bliss::Graph* bliss_graph) {
    convert_dejavu_to_bliss(dejavu_graph->get_sgraph(), dejavu_graph->get_coloring(), bliss_graph);
}


#endif //DEJAVU_BLISS_CONVERTER_H
