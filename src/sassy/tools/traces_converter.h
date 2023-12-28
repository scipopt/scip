#ifndef SASSY_TRACES_CONVERTER_H
#define SASSY_TRACES_CONVERTER_H

#include "nauty_converter.h"

static void convert_sassy_to_traces(sassy::sgraph *g, int* colmap, sparsegraph* sg, int** lab, size_t* lab_sz, int** ptn, size_t* ptn_sz) {
    SG_INIT(*sg);
    SG_ALLOC(*sg, g->v_size, g->e_size, "malloc");
    sg->nv = g->v_size;
    sg->nde = g->e_size;

    DYNALLOC1(int, *lab, (*lab_sz), (g->v_size), "malloc");
    DYNALLOC1(int, *ptn, (*ptn_sz), (g->v_size), "malloc");

    for (int i = 0; i < g->v_size; ++i) {
        (*lab)[i] = i;
        (*ptn)[i] = 1;
        sg->v[i] = g->v[i];
        sg->d[i] = g->d[i];
    }
    if(g->v_size > 0)
        (*ptn)[g->v_size - 1] = 0;
    if(colmap != nullptr) {
        make_lab_ptn_from_colmap(*lab, *ptn, colmap, g->v_size);
    }

    for (int i = 0; i < g->e_size; ++i) {
        sg->e[i] = g->e[i];
    }
}


static void convert_sassy_to_traces(sassy::static_graph *sassy_graph, sparsegraph* sg, int** lab, size_t* lab_sz, int** ptn, size_t* ptn_sz) {
    convert_sassy_to_traces(sassy_graph->get_sgraph(), sassy_graph->get_coloring(), sg, lab, lab_sz, ptn, ptn_sz);
}


#endif //SASSY_TRACES_CONVERTER_H
