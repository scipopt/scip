#ifndef SASSY_COLORING_H
#define SASSY_COLORING_H


#include <vector>
#include <cassert>
#include "utility.h"

namespace sassy {
    class coloring {
    public:
        int *lab;
        int *ptn;
        int lab_sz;
        int ptn_sz;

        bool init = false;
        int *vertex_to_col;
        int *vertex_to_lab;

        int cells = 1;

        ~coloring() {
            if (init) {
                dealloc();
            }
        }

        void alloc(int sz) {
            if (!init) {
                lab = new int[sz];
                ptn = new int[sz];
                vertex_to_col = new int[sz];
                vertex_to_lab = new int[sz];
                init = true;

                lab_sz = sz;
                ptn_sz = sz;
            }
        }

        void dealloc() {
            delete[] ptn;
            delete[] lab;
            delete[] vertex_to_lab;
            delete[] vertex_to_col;
        };

        void copy_ptn(coloring *c) {
            assert(init);
            assert(c->init);
            memcpy(ptn, c->ptn, c->ptn_sz * sizeof(int));
        }

        void copy(coloring *c) {
            if (init) {
                if (lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
                    dealloc();
                    init = false;
                } else {
                    cells = c->cells;
                    for (int i = 0; i < c->ptn_sz;) {
                        const int rd = c->ptn[i];
                        ptn[i] = rd;
                        i += rd + 1;
                    }
                    memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz * sizeof(int));
                    return;
                }
            }

            if (!init) {
                alloc(c->lab_sz);
            }

            if (c->cells > c->ptn_sz / 4) {
                memcpy(ptn, c->ptn, c->ptn_sz * sizeof(int));
            } else {
                for (int i = 0; i < c->ptn_sz;) {
                    const int rd = c->ptn[i];
                    ptn[i] = rd;
                    i += rd + 1;
                }
            }
            memcpy(lab, c->lab, c->lab_sz * sizeof(int));
            memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz * sizeof(int));
            memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz * sizeof(int));

            lab_sz = c->lab_sz;
            ptn_sz = c->ptn_sz;

            cells = c->cells;
            init = true;
        }

        void copy_force(coloring *c) {
            if (init) {
                if (lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
                    dealloc();
                    init = false;
                }
            }

            if (!init) {
                alloc(c->lab_sz);
            }

            if (c->cells > c->ptn_sz / 4) {
                memcpy(ptn, c->ptn, c->ptn_sz * sizeof(int));
            } else {
                for (int i = 0; i < c->ptn_sz;) {
                    const int rd = c->ptn[i];
                    ptn[i] = rd;
                    i += rd + 1;
                }
            }
            memcpy(lab, c->lab, c->lab_sz * sizeof(int));
            memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz * sizeof(int));
            memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz * sizeof(int));

            lab_sz = c->lab_sz;
            ptn_sz = c->ptn_sz;

            cells = c->cells;
            init = true;
        }

        void initialize(int domain_size) {
            alloc(domain_size);
        }

        bool check() {
            bool comp = true;

            for (int i = 0; i < lab_sz; ++i) {
                comp = comp && (lab[i] >= 0 && lab[i] < lab_sz);
                comp = comp && (lab[vertex_to_lab[i]] == i);
            }

            int last_col = -1;
            int counter = 1;
            for (int i = 0; i < ptn_sz; ++i) {
                --counter;
                if (counter == 0) {
                    counter = ptn[i] + 1;
                    assert(ptn[i] >= 0 && ptn[i] < lab_sz);
                } else {
                    assert(vertex_to_col[lab[i]] == last_col);
                }
            }

            for (int i = 0; i < ptn_sz;) {
                assert(vertex_to_col[lab[i]] == i);
                i += ptn[i] + 1;
            }
            return comp;
        }
    };
}

#endif //SASSY_COLORING_H
