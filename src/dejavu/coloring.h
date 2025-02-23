// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_COLORING_H
#define DEJAVU_COLORING_H

#include <vector>
#include <cassert>
#include <cstring>
#include "utility.h"

namespace dejavu { namespace ds {
    // (need to nest namespaces due to C++ 14)

    /**
     * \brief Vertex coloring for a graph
     *
     * Stores a vertex coloring for a graph with \a domain_size many vertices. Has mappings to and from colorings. The
     * datastructure mostly follows the design used by Traces, with some slight adjustments. The class mainly exposes
     * raw arrays, to be used by the solver internally.
     */
    class coloring {
    public:
        int *lab = nullptr;            /**< vertices of color `i` are listed in `lab[i]...lab[ptn[i]]`                */
        int *ptn = nullptr;            /**< color `i` has size `ptn[i] + 1`                                           */
        int *vertex_to_col = nullptr;  /**< vertex `v` has color `vertex_to_col[v]`                                   */
        int *vertex_to_lab = nullptr;  /**< vertex `v` is stored in lab at position `vertex_to_lab[v]`, i.e.,
                                        `lab[vertex_to_lab[v]] = v`                                                   */

        int cells       = 1;           /**< number of colors (i.e., cells) in the coloring                            */
        int domain_size = 0;           /**< domain size of the coloring (i.e., number of vertices of the graph)       */

    private:
        int *alloc_pt = nullptr;

        /**
         * Allocates this coloring to \a domain_size of \p sz. Does not initialize any content and runs in O(1).
         *
         * @param sz Domain size.
         */
        void alloc(int sz) {
            dej_assert(sz >= 0);

            if (alloc_pt) dealloc();
            alloc_pt = (int *) malloc(sz * 4 * sizeof(int));
            dej_assert(alloc_pt != nullptr);

            lab = alloc_pt;
            ptn = lab + sz;
            vertex_to_col = lab + sz * 2;
            vertex_to_lab = lab + sz * 3;

            domain_size = sz;
        }

        /**
         * Frees up memory used by this coloring.
         */
        void dealloc() {
            if (alloc_pt) free(alloc_pt);
            lab = nullptr;
            ptn = nullptr;
            vertex_to_col = nullptr;
            vertex_to_lab = nullptr;

            domain_size = 0;
        };

    public:
        ~coloring() {
            if (alloc_pt) {
                dealloc();
            }
        }

        /**
         * This function only copies the `ptn` array from \p c to this coloring. Assumes \p c and this coloring have the
         * same \a domain_size and are both allocated.
         *
         * @param c The coloring to copy from.
         */
        void copy_ptn(coloring *c) const {
            dej_assert(alloc_pt);
            dej_assert(c->alloc_pt);
            dej_assert(domain_size == c->domain_size);
            memcpy(ptn, c->ptn, c->domain_size * sizeof(int));
        }

        /**
         * Copies the given coloring into this coloring, but is only guaranteed to produce a correct coloring if \p c
         * stems from an ancestor in the IR tree of the currently stored coloring.
         *
         * Essentially, this function only copies the `vertex_to_col` array.
         *
         * @param c The coloring to copy from.
         */
        void copy_from_ir_ancestor(coloring *c) {
            if (alloc_pt) {
                if (domain_size != c->domain_size) {
                    dealloc();
                } else {
                    cells = c->cells;
                    memcpy(vertex_to_col, c->vertex_to_col, c->domain_size * sizeof(int));
                    return;
                }
            }

            if (!alloc_pt) alloc(c->domain_size);

            memcpy(ptn, c->ptn, c->domain_size * sizeof(int));
            memcpy(lab, c->lab, c->domain_size * sizeof(int));
            memcpy(vertex_to_col, c->vertex_to_col, c->domain_size * sizeof(int));
            memcpy(vertex_to_lab, c->vertex_to_lab, c->domain_size * sizeof(int));

            domain_size = c->domain_size;
            cells = c->cells;
        }

        /**
         * Copies the given coloring into this coloring.
         *
         * @param c The coloring to copy from.
         */
        void copy_any(coloring *c) {
            if (alloc_pt) {
                if (domain_size != c->domain_size) dealloc();
            }

            if (!alloc_pt) {
                alloc(c->domain_size);
            }

            memcpy(ptn, c->ptn, c->domain_size * sizeof(int));
            memcpy(lab, c->lab, c->domain_size * sizeof(int));
            memcpy(vertex_to_col, c->vertex_to_col, c->domain_size * sizeof(int));
            memcpy(vertex_to_lab, c->vertex_to_lab, c->domain_size * sizeof(int));

            domain_size = c->domain_size;
            cells = c->cells;
        }

        /**
         * Initialize this coloring with a given domain size. Allocates memory for the coloring. Does not initialize any
         * content and runs in O(1).
         *
         * @param new_domain_size domain size of the coloring
         */
        void initialize(int new_domain_size) {
            alloc(new_domain_size);
        }

        void check() const {
            bool comp = true;

            for (int i = 0; i < domain_size; ++i) {
                comp = comp && (lab[i] >= 0 && lab[i] < domain_size);
                comp = comp && (lab[vertex_to_lab[i]] == i);
            }

            int counter = 1;
            for (int i = 0; i < domain_size; ++i) {
                --counter;
                if (counter == 0) {
                    counter = ptn[i] + 1;
                    dej_assert(ptn[i] >= 0 && ptn[i] < domain_size);
                }
            }

            for (int i = 0; i < domain_size;) {
                dej_assert(vertex_to_col[lab[i]] == i);
                i += ptn[i] + 1;
            }
            dej_assert(comp);
        }
    };
}}

#endif //DEJAVU_COLORING_H
