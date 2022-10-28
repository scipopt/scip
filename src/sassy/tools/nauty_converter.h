#ifndef SASSY_NAUTY_CONVERTER_H
#define SASSY_NAUTY_CONVERTER_H

#include <algorithm>

struct colorComparator {
    colorComparator(const int* colmap) : colmap(colmap) {}
    const int* colmap;

    bool operator()(const int & v1, const int & v2) {
        return (colmap[v1] < colmap[v2]);
    }
};

static void make_lab_ptn_from_colmap(int* lab, int* ptn, const int* colmap, int colmap_sz) {
    for(int i = 0; i < colmap_sz; ++i) {
        lab[i] = i;
        ptn[i] = 1;
    }
    int last_new_cell = 0;
    std::sort(lab, lab + colmap_sz, colorComparator(colmap));
    for(int i = 0; i < colmap_sz; i++) {
        if (i + 1 == colmap_sz) {
            ptn[last_new_cell] = (i - last_new_cell) > 0;
            ptn[i] = 0;
            break;
        }
        if (colmap[lab[i]] != colmap[lab[i + 1]]) {
            ptn[i] = 0;
            ptn[last_new_cell] = (i - last_new_cell) > 0;
            last_new_cell = i + 1;
            continue;
        }
    }
}

static void convert_sassy_to_nauty(sassy::sgraph *g, int* colmap, sparsegraph* sg, int** lab, size_t* lab_sz, int** ptn, size_t* ptn_sz) {
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

static void convert_sassy_to_nauty(sassy::static_graph *sassy_graph, sparsegraph* sg, int** lab, size_t* lab_sz, int** ptn, size_t* ptn_sz) {
    convert_sassy_to_nauty(sassy_graph->get_sgraph(), sassy_graph->get_coloring(), sg, lab, lab_sz, ptn, ptn_sz);
}


#endif //SASSY_NAUTY_CONVERTER_H
