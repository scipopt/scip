#ifndef SASSY_SAUCY_CONVERTER_H
#define SASSY_SAUCY_CONVERTER_H

#include <algorithm>

struct colorComparatorSaucy {
    colorComparatorSaucy(const int* colmap) : colmap(colmap) {}
    const int* colmap;

    bool operator()(const int & v1, const int & v2) {
        return (colmap[v1] < colmap[v2]);
    }
};

static void make_saucy_colmap(int* lab, int* smallcolmap, int* colmap, int colmap_sz) {
    for(int i = 0; i < colmap_sz; ++i) {
        lab[i] = i;
    }

    // int last_new_cell = 0;
    std::sort(lab, lab + colmap_sz, colorComparatorSaucy(colmap));

    int col_small = -1;
    int last_col_orig  = INT32_MAX;

    for(int i = 0; i < colmap_sz; i++) {
        if (colmap[lab[i]] != last_col_orig) {
            last_col_orig = colmap[lab[i]];
            ++col_small;
        }
        smallcolmap[lab[i]] = col_small;
    }
}

static void convert_sassy_to_saucy(sassy::sgraph *g, int* colmap, saucy_graph* _saucy_graph, int** colors) {
    *colors  = new int[g->v_size];
    int* lab = new int[g->v_size];
    const bool create_colmap = (colmap == nullptr);
    if(create_colmap) {
        colmap = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            colmap[i] = 0;
        }
    }
    make_saucy_colmap(lab, *colors, colmap, g->v_size);

    _saucy_graph->n = g->v_size;
    _saucy_graph->e = g->e_size;
    _saucy_graph->edg = new int[g->e_size];
    _saucy_graph->adj = new int[_saucy_graph->n + 1];

    int epos = 0;
    for(int i = 0; i < g->v_size; ++i) {
        const int npt = g->v[i];
        const int nd  = g->d[i];
        _saucy_graph->adj[i] = epos;
        for(int j = 0; j < nd; ++j) {
            const int adj_n = g->e[npt + j];
            _saucy_graph->edg[epos] = adj_n;
            ++epos;
        }
    }

    _saucy_graph->adj[g->v_size] = epos;
    delete[] lab;
    if(create_colmap) {
        delete[] colmap;
    }
}

static void convert_sassy_to_saucy(sassy::static_graph *sassy_graph, saucy_graph* _saucy_graph, int** colors){
    convert_sassy_to_saucy(sassy_graph->get_sgraph(), sassy_graph->get_coloring(), _saucy_graph, colors);
}


#endif //SASSY_SAUCY_CONVERTER_H
