#ifndef SASSY_GRAPH_H
#define SASSY_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <set>
#include <cstring>

namespace sassy {
    template<class vertex_t, class degree_t, class edge_t>
    class sgraph_t {
        struct vertexComparator {
            vertexComparator(const sgraph_t<vertex_t, degree_t, edge_t> &g) : g(g) {}

            const sgraph_t &g;

            bool operator()(const vertex_t &v1, const vertex_t &v2) {
                return g.d[v1] < g.d[v2];
            }
        };

        struct vertexComparatorColor {
            vertexComparatorColor(const sgraph_t<vertex_t, degree_t, edge_t> &g, const vertex_t *vertex_to_col) : g(g),
                                                                                                                  vertex_to_col(
                                                                                                                          vertex_to_col) {}

            const sgraph_t &g;
            const vertex_t *vertex_to_col;

            bool operator()(const vertex_t &v1, const vertex_t &v2) {
                //return (g.d[v1] < g.d[v2]) || ((g.d[v1] == g.d[v2]) && (vertex_to_col[v1] < vertex_to_col[v2]));
                return (vertex_to_col[v1] < vertex_to_col[v2]);
            }
        };

    public:
        bool initialized = false;
        edge_t *v;
        degree_t *d;
        vertex_t *e;

        int v_size;
        int d_size;
        int e_size;

        bool dense = false;

        int max_degree;

        void initialize(int nv, int ne) {
            initialized = true;
            v = new edge_t[nv];
            d = new degree_t[nv];
            e = new vertex_t[ne];
        }

        // initialize a coloring of this sgraph, partitioning degrees of vertices
        void initialize_coloring(coloring *c, vertex_t *vertex_to_col) {
            c->alloc(this->v_size);
            std::memset(c->ptn, 1, sizeof(int) * v_size);

            if (this->v_size < c->lab_sz) {
                c->lab_sz = this->v_size;
                c->ptn_sz = this->v_size;
            }

            if (v_size == 0)
                return;

            int cells = 0;
            int last_new_cell = 0;

            if (vertex_to_col == nullptr) {
                for (int i = 0; i < v_size; i++) {
                    c->lab[i] = i;
                }
                std::sort(c->lab, c->lab + c->lab_sz, vertexComparator(*this));
                for (int i = 0; i < c->lab_sz; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->lab_sz) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
                    if (this->d[c->lab[i]] != this->d[c->lab[i + 1]]) {
                        c->ptn[i] = 0;
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        last_new_cell = i + 1;
                        continue;
                    }
                }
            } else {
                int col = 0;

                int min_col = INT32_MAX;
                int max_col = INT32_MIN;
                for (int i = 0; i < v_size; i++) {
                    if (vertex_to_col[i] < min_col)
                        min_col = vertex_to_col[i];
                    if (vertex_to_col[i] > max_col)
                        max_col = vertex_to_col[i];
                }

                std::vector<int> colsize;
                colsize.reserve(std::min(this->v_size, (max_col - min_col) + 1));

                if (min_col < 0 || max_col > 4 * this->v_size) {
                    std::unordered_map<int, int> colors; // TODO: should not use unordered_map!
                    colors.reserve(this->v_size);
                    for (int i = 0; i < v_size; i++) {
                        auto it = colors.find(vertex_to_col[i]);
                        if (it == colors.end()) {
                            colors.insert(std::pair<int, int>(vertex_to_col[i], col));
                            colsize.push_back(1);
                            assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = it->second;
                            assert(found_col < this->v_size);
                            vertex_to_col[i] = found_col;
                            ++colsize[found_col];
                        }
                    }
                } else {
                    std::vector<int> colors;
                    colors.reserve(max_col + 1);
                    for (int i = 0; i < max_col + 1; ++i)
                        colors.push_back(-1);
                    for (int i = 0; i < v_size; i++) {
                        if (colors[vertex_to_col[i]] == -1) {
                            colors[vertex_to_col[i]] = col;
                            colsize.push_back(1);
                            assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = colors[vertex_to_col[i]];
                            assert(found_col < this->v_size);
                            vertex_to_col[i] = found_col;
                            ++colsize[found_col];
                        }
                    }
                }

                int increment = 0;
                for (size_t i = 0; i < colsize.size(); i++) {
                    const int col_sz = colsize[i];
                    colsize[i] += increment;
                    increment += col_sz;
                }
                assert(increment == v_size);

                // cache inefficiency probably starts here... try to make these procedures more sequential
                for (int i = 0; i < v_size; i++) {
                    const int v_col = vertex_to_col[i];
                    --colsize[v_col];
                    const int v_lab_pos = colsize[v_col];
                    c->lab[v_lab_pos] = i;
                }

                for (int i = 0; i < c->lab_sz; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->lab_sz) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    if (vertex_to_col[c->lab[i]] != vertex_to_col[c->lab[i + 1]]) {
                        c->ptn[i] = 0;
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        last_new_cell = i + 1;
                        continue;
                    }
                }
            }

            c->cells = cells;
        }

        void initialize_coloring_raw(coloring *c) {
            c->alloc(this->v_size);

            for (int i = 0; i < v_size; i++) {
                c->lab[i] = i;
                c->vertex_to_lab[i] = i;
            }

            std::memset(c->vertex_to_col, 0, sizeof(int) * v_size);
            c->ptn[0] = v_size - 1;
            c->ptn[v_size - 1] = 0;
            c->cells = 1;
        }

        void sanity_check() {
#ifndef NDEBUG
            for (int i = 0; i < v_size; ++i) {
                assert(d[i] > 0 ? v[i] < e_size : true);
                assert(d[i] > 0 ? v[i] >= 0 : true);
                assert(d[i] >= 0);
                assert(d[i] < v_size);
            }
            for (int i = 0; i < e_size; ++i) {
                assert(e[i] < v_size);
                assert(e[i] >= 0);
            }

            // multiedge test
            mark_set multiedge_test;
            multiedge_test.initialize(v_size);
            for (int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for (int j = 0; j < d[i]; ++j) {
                    const int neigh = e[v[i] + j];
                    assert(!multiedge_test.get(neigh));
                    multiedge_test.set(neigh);
                }
            }

            /*// fwd - bwd test
            multiedge_test.initialize(v_size);
            for (int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for (int j = 0; j < d[i]; ++j) {
                    const int neigh = e[v[i] + j];
                    bool found = false;
                    for (int k = 0; k < d[neigh]; ++k) {
                        const int neigh_neigh = e[v[neigh] + k];
                        if (neigh_neigh == i) {
                            found = true;
                            break;
                        }
                    }
                    assert(found);
                }
            }*/
#endif
        }

        void permute_graph(sgraph_t<vertex_t, degree_t, edge_t>* ng, int* p) {
            ng->initialize(v_size, e_size);
            ng->v_size = v_size;
            ng->d_size = d_size;
            ng->e_size = e_size;
            ng->max_degree = max_degree;

            bijection<vertex_t> p_inv;
            p_inv.initialize_empty(ng->v_size);
            //p_inv.copy(p);
            //p_inv.initialize_empty(ng->v_size);
            p_inv.read_from_array(p, ng->v_size);
            p_inv.inverse();

            int epos = 0;
            for(int i = 0; i < v_size; ++i) {
                int mapped_v = p[i];
                assert(p_inv.map_vertex(mapped_v) == i);
                assert(mapped_v < v_size);
                ng->d[i] = d[mapped_v];
                ng->v[i] = epos;
                for(int j = v[mapped_v]; j < v[mapped_v] + d[mapped_v]; j++) {
                    assert(j < e_size);
                    ng->e[epos] = p_inv.map_vertex(e[j]);
                    epos += 1;
                }
                //epos += ng->d[i];
            }

            assert(ng->v_size == v_size);
            assert(ng->e_size == e_size);
            assert(ng->d_size == d_size);
            assert(epos == ng->e_size);

            return;
        }

        void copy_graph(sgraph_t<vertex_t, degree_t, edge_t> *g) {
            if (initialized) {
                delete[] v;
                delete[] d;
                delete[] e;
            }
            initialize(g->v_size, g->e_size);

            memcpy(v, g->v, g->v_size * sizeof(edge_t));
            memcpy(d, g->d, g->d_size * sizeof(degree_t));
            memcpy(e, g->e, g->e_size * sizeof(vertex_t));
            v_size = g->v_size;
            d_size = g->d_size;
            e_size = g->e_size;
            max_degree = g->max_degree;
        }



        void print() {
            std::cout << "v_size: " << v_size << std::endl;
            std::cout << "d_size: " << d_size << std::endl;
            std::cout << "e_size: " << e_size << std::endl;
        }

        ~sgraph_t() {
            if (initialized) {
                delete[] v;
                delete[] d;
                delete[] e;
            }
        }
    };


    template<class vertex_type_src, class degree_type_src, class edge_type_src,
            class vertex_type_tgt, class degree_type_tgt, class edge_type_tgt>
    static void copy_graph(sgraph_t<vertex_type_src, degree_type_src, edge_type_src> *g1,
                           sgraph_t<vertex_type_tgt, degree_type_tgt, edge_type_tgt> *g2) {
        g2->v_size = g1->v_size;
        g2->d_size = g1->d_size;
        g2->e_size = g1->e_size;
        g2->max_degree = g1->max_degree;

        g2->initialize(g2->v_size, g2->e_size);

        for (int i = 0; i < g1->v_size; ++i) {
            g2->v[i] = static_cast<edge_type_tgt>(g1->v[i]);
        }
        for (int i = 0; i < g1->d_size; ++i) {
            g2->d[i] = static_cast<degree_type_tgt>(g1->d[i]);
        }
        for (int i = 0; i < g1->e_size; ++i) {
            g2->e[i] = static_cast<vertex_type_tgt>(g1->e[i]);
        }
    }

    enum sgraph_type {
        DSG_INT_INT_INT, DSG_SHORT_SHORT_INT, DSG_SHORT_SHORT_SHORT, DSG_CHAR_CHAR_SHORT, DSG_CHAR_CHAR_CHAR
    };

    template<class vertex_t, class degree_t, class edge_t>
    sgraph_t<vertex_t, degree_t, edge_t> *disjoint_union(sgraph_t<vertex_t, degree_t, edge_t> *g1,
                                                         sgraph_t<vertex_t, degree_t, edge_t> *g2) {
        sgraph_t<vertex_t, degree_t, edge_t> *union_g = new sgraph_t<vertex_t, degree_t, edge_t>();
        union_g->v_size = g1->v_size + g2->v_size;
        union_g->e_size = g1->e_size + g2->e_size;
        union_g->d_size = g1->d_size + g2->d_size;

        int g2_vshift = g1->v_size;
        int g2_eshift = g1->e_size;
        int g2_dshift = g1->d_size;

        union_g->initialize(union_g->v_size, union_g->e_size);

        for (int i = 0; i < g1->v_size; ++i)
            union_g->v[i] = g1->v[i];
        for (int i = 0; i < g2->v_size; ++i)
            union_g->v[i + g2_vshift] = g2->v[i] + g2_eshift;

        for (int i = 0; i < g1->d_size; ++i)
            union_g->d[i] = g1->d[i];
        for (int i = 0; i < g2->d_size; ++i)
            union_g->d[i + g2_dshift] = g2->d[i];

        for (int i = 0; i < g1->e_size; ++i)
            union_g->e[i] = g1->e[i];
        for (int i = 0; i < g2->e_size; ++i)
            union_g->e[i + g2_eshift] = g2->e[i] + g2_vshift;

        union_g->max_degree = std::max(g1->max_degree, g2->max_degree);

        return union_g;
    }

    typedef sgraph_t<int, int, int> sgraph;
}

#endif //SASSY_GRAPH_H
