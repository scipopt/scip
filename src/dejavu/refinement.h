// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_REFINEMENT_H
#define DEJAVU_REFINEMENT_H

#include "ds.h"
#include "coloring.h"
#include "graph.h"
#include "utility.h"

namespace dejavu {

    using namespace ds;

    namespace ir {

        /**
         * \brief Certifies automorphisms
         *
         * Contains methods to certify whether a given vertex mapping is a symmetry of a given graph.
         */
        class certification {
        private:
            static bool bijection_check(markset& scratch_set, int n, const int *p) {
                scratch_set.reset();
                bool comp = true;
                for(int i = 0; i < n && comp; ++i) {
                    const int v = p[i];
                    comp = !scratch_set.get(v);
                    scratch_set.set(v);
                }
                return comp;
            }

            static bool cycle_check(markset& scratch_set, const int *p, int supp, const int *supp_arr) {
                scratch_set.reset();
                bool comp = true;
                for(int i = 0; i < supp && comp; ++i) {
                    const int v = supp_arr[i];
                    if(scratch_set.get(v)) continue;
                    int v_next = p[v];
                    while(v_next != v && comp) {
                        comp = !scratch_set.get(v_next);
                        scratch_set.set(v_next);
                        v_next = p[v_next];
                    }
                }
                return comp;
            }

        public:
            // certify an automorphism on a graph
            static bool certify_automorphism(markset& scratch_set, sgraph *g, const int *p) {
                if(!bijection_check(scratch_set, g->v_size, p)) return false;

                for (int i = 0; i < g->v_size; ++i) {
                    const int image_i = p[i];
                    if (image_i == i) continue;

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    int found = 0;
                    const int start_pt = g->v[i];
                    const int end_pt   = g->v[i] + g->d[i];
                    for (int j = start_pt; j < end_pt; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    const int image_start_pt = g->v[image_i];
                    const int image_end_pt   = g->v[image_i] + g->d[image_i];
                    for (int j = image_start_pt; j < image_end_pt; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) return false;
                        scratch_set.unset(vertex_j);
                        found -= 1;
                    }
                    if (found != 0) return false;
                }

                return true;
            }

            // certify an automorphism on a graph, sparse
            static bool certify_automorphism_sparse(markset& scratch_set, const sgraph *g, const int *p, int supp,
                                                    const int *supp_arr) {
                int i, found;
                if(!cycle_check(scratch_set, p, supp, supp_arr)) return false;

                for (int f = 0; f < supp; ++f) {
                    i = supp_arr[f];
                    const int image_i = p[i];
                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    const int start_pt = g->v[i];
                    const int end_pt   = g->v[i] + g->d[i];
                    for (int j = start_pt; j < end_pt; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    const int image_start_pt = g->v[image_i];
                    const int image_end_pt   = g->v[image_i] + g->d[image_i];
                    for (int j = image_start_pt; j < image_end_pt; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            scratch_set.reset();
                            return false;
                        }
                        found -= 1;
                    }
                    if (found != 0) {
                        scratch_set.reset();
                        return false;
                    }
                }
                scratch_set.reset();
                return true;
            }
        };


        // return whether to continue color refinement
        // bool split_color_hook(int color_initial, int new_color, int new_color_sz);
        typedef bool type_split_color_hook(const int, const int, const int);

        // return whether to continue splitting the respective cell, or skip it
        // bool worklist_color_hook(int color, int color_sz);
        typedef bool type_worklist_color_hook(const int, const int);

        /**
         * \brief Color refinement and related algorithms
         *
         * Class that is used to preserve a workspace for color refinement, automorphism certification, and other miscellaneous
         * tasks. Once initialized to a certain size, the workspace can not be enlargened and methods can only be used for
         * graphs of the initial size, or smaller.
        */
        class refinement {
            bool g_early_out = false;
            const std::function<type_split_color_hook>* g_split_hook;

        public:
            /**
             * The color refinement algorithm. Refines a given coloring with respect to a given graph.
             * @param g The graph.
             * @param c The coloring to be refined.
             * @param init_color Initialize the worklist with a single color class (e.g., after individualization). The
             * default value -1 denotes that the worklist is initialized with all color classes of the coloring.
             * @param color_limit Integer which is used to stop refinement whenever the refined coloring reaches this number
             * of color classes. The default value -1 denotes that refinement is performed exhaustively.
             * @param split_hook Function pointer that is called whenever a color class is split. Return value can be used
                 * to stop refinement early.
             * @param worklist_hook Function pointer that is called whenever a color class is considered for refinement. Return value
             * can be used to skip refinement of that color class.
             */
            void refine_coloring(sgraph *g, coloring *c, int init_color = -1, int color_limit = -1,
                                 const std::function<type_split_color_hook>* split_hook = nullptr,
                                 const std::function<type_worklist_color_hook>* worklist_hook = nullptr) {
                assure_initialized(g);
                cell_todo.reset(queue_pointer);

                if (init_color < 0) {
                    // initialize queue with all classes (except for largest one)
                    for (int i = 0; i < c->domain_size;) {
                        cell_todo.add_cell(queue_pointer, i);
                        const int col_sz = c->ptn[i];
                        i += col_sz + 1;
                    }
                } else {
                    dej_assert(c->vertex_to_col[c->lab[init_color]] == init_color);
                    cell_todo.add_cell(queue_pointer, init_color);
                }

                g_early_out  = false;
                g_split_hook = split_hook;

                while (!cell_todo.empty()) {
                    const int next_color_class    = cell_todo.next_cell(queue_pointer, c);
                    const int next_color_class_sz = c->ptn[next_color_class] + 1;

                    if (worklist_hook && !(*worklist_hook)(next_color_class, next_color_class_sz)) continue;

                    // this scheme is reverse-engineered from the color refinement in Traces by Adolfo Piperno
                    // we choose a separate algorithm depending on the size and density of the graph and/or color class
                    const int  test_deg   = g->d[c->lab[next_color_class]];
                    const bool very_dense = test_deg > (g->v_size / (next_color_class_sz + 1));
                    //const bool cell_dense = test_deg > (c->cells);

                    if (next_color_class_sz == 1 && !(g->dense && very_dense)) { // singleton
                        refine_color_class_singleton(g, c, next_color_class);
                    } else if (g->dense) {
                        if (very_dense) { // dense-dense
                            refine_color_class_dense_dense(g, c, next_color_class,next_color_class_sz);
                        } //else if(cell_dense) { // dense-cell
                            //refine_color_class_dense_cell(g, c, next_color_class, next_color_class_sz);
                        //}
                        else { // dense-sparse
                            refine_color_class_dense(g, c, next_color_class, next_color_class_sz);
                        }
                    } else { // sparse
                        refine_color_class_sparse(g, c, next_color_class, next_color_class_sz);
                    }

                    // early out possible?
                    if (g_early_out || c->cells == g->v_size || c->cells == color_limit) {
                        cell_todo.reset(queue_pointer);
                        break;
                    }
                }
            }
        private:
            void report_split_color_class(coloring* c, const int old_class, const int new_class, const int new_class_sz,
                                          const bool is_largest) {
                c->cells += (old_class != new_class);
                dej_assert(c->ptn[new_class] + 1 == new_class_sz);

                if ((g_split_hook != nullptr) && !(*g_split_hook)(old_class, new_class, new_class_sz)) {
                    g_early_out = true;
                }

                if (!is_largest && old_class != new_class) {
                    cell_todo.add_cell(queue_pointer, new_class);
                } else if(is_largest) {
                    // since old color class is skipped above, this should be safe
                    int i = queue_pointer.get(old_class);
                    if (i >= 0)                cell_todo.replace_cell(queue_pointer, old_class, new_class);
                    if(old_class != new_class) cell_todo.add_cell(queue_pointer, old_class);
                }
            }

        public:
            /**
             * Individualizes a vertex in a coloring.
             * @param c Coloring in which the vertex is individualized.
             * @param v Vertex to be individualized.
             * @param split_hook Function to be called whenever a color class is split. Return value is not used.
             * @return
             */
            static int
            individualize_vertex(coloring *c, int v, const std::function<type_split_color_hook> &split_hook = nullptr) {
                const int color = c->vertex_to_col[v];
                const int pos = c->vertex_to_lab[v];

                int color_class_size = c->ptn[color];

                dej_assert(color_class_size > 0);

                const int vertex_at_pos = c->lab[color + color_class_size];
                c->lab[pos] = vertex_at_pos;
                c->vertex_to_lab[vertex_at_pos] = pos;

                c->lab[color + color_class_size] = v;
                c->vertex_to_lab[v] = color + color_class_size;
                c->vertex_to_col[v] = color + color_class_size;

                c->ptn[color] -= 1;
                c->ptn[color + color_class_size] = 0;
                c->ptn[color + color_class_size - 1] = 0;
                c->cells += 1;

                if (split_hook) {
                    split_hook(color, color + color_class_size, 1);
                    split_hook(color, color, c->ptn[color] + 1);
                }

                return color + color_class_size;
            }

            /**
             * The color refinement algorithm. Refines a given coloring with respect to a given graph. This color
             * refinement does not produce an isomorphism-invariant ordered partitioning, but uses more optimization
             * techniques -- meant to be used as the first refinement in automorphism computation
             *
             * @param g The graph.
             * @param c The coloring to be refined.
             * @param init_color Initialize the worklist with a single color class (e.g., after individualization). The
             * default value -1 denotes that the worklist is initialized with all color classes of the coloring.
             */
            void refine_coloring_first(sgraph *g, coloring *c, int init_color_class = -1) {
                assure_initialized(g);
                singleton_hint.reset();

                cell_todo.reset(queue_pointer);

                if (init_color_class < 0) {
                    for (int i = 0; i < c->domain_size;) {
                        cell_todo.add_cell(queue_pointer, i);
                        const int col_sz = c->ptn[i];
                        if (col_sz == 0) {
                            singleton_hint.push_back(i);
                        }
                        i += col_sz + 1;
                    }
                } else {
                    cell_todo.add_cell(queue_pointer, init_color_class);
                }

                while (!cell_todo.empty()) {
                    const int next_color_class = cell_todo.next_cell(queue_pointer, c, singleton_hint);
                    const int next_color_class_sz = c->ptn[next_color_class] + 1;
                    const bool very_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

                    if (next_color_class_sz == 1 && !(g->dense && very_dense)) {
                        // singleton
                        refine_color_class_singleton_first(g, c, next_color_class);
                    } else if (g->dense) {
                        if (very_dense) { // dense-dense
                            refine_color_class_dense_dense_first(g, c, next_color_class, next_color_class_sz);
                        } else { // dense-sparse
                            refine_color_class_dense_first(g, c, next_color_class, next_color_class_sz);
                        }
                    } else { // sparse
                        refine_color_class_sparse_first(g, c, next_color_class, next_color_class_sz);
                    }

                    if (c->cells == g->v_size) {
                        cell_todo.reset(queue_pointer);
                        return;
                    }
                }
            }

        private:
            void report_split_color_class_first(coloring* c, int old_class, int new_class, int new_class_sz,
                                                bool is_largest) {
                c->cells += (old_class != new_class);
                dej_assert(c->ptn[new_class] + 1 == new_class_sz);

                if (!is_largest && old_class != new_class) {
                    cell_todo.add_cell(queue_pointer, new_class);
                    if (new_class_sz == 1) singleton_hint.push_back(new_class);
                } else if(is_largest) {
                    // since old color class is skipped above, this should be safe
                    int i = queue_pointer.get(old_class);
                    if (i >= 0) {
                        cell_todo.replace_cell(queue_pointer, old_class, new_class);
                        if(new_class_sz == 1) singleton_hint.push_back(new_class);
                    }
                    if(old_class != new_class) {
                        cell_todo.add_cell(queue_pointer, old_class);
                        if(c->ptn[old_class] + 1 == 1) singleton_hint.push_back(old_class);
                    }
                }
            }

        public:
            /**
             * Certifies automorphism of a graph. Uses the `certification` methods, but provides the auxiliary data
             * structure.
             *
             * \sa certification
             *
             * @param g the graph
             * @param p the mapping to be certified
             * @return whether \a p is an automorphism of g
             */
            bool  certify_automorphism(sgraph *g, const int *p) {
                assure_initialized(g);
                return certification::certify_automorphism(scratch_set, g, p);
            }

            /**
             * Certifies (sparse) automorphism of a graph. Uses the `certification` methods, but provides the
             * auxiliary data structure. Does not certify points outside the given support.
             *
             * \sa certification
             *
             * @param g the graph
             * @param p the mapping to be certified
             * @param supp the number of points in the support of p
             * @param supp_arr the points in the support of p
             * @return whether \a p is an automorphism of g (assuming the support is provided correctly)
             */
            bool certify_automorphism_sparse(const sgraph *g, const int *p, int supp, const int *supp_arr) {
                assure_initialized(g);
                return certification::certify_automorphism_sparse(scratch_set, g, p, supp, supp_arr);
            }

        private:
            // worklist implementation for color refinement
            class cell_worklist {
            public:
                void initialize(int _domain_size) {
                    arr = new int[_domain_size];
                    arr_sz = _domain_size;
                    init = true;
                    cur_pos = 0;
                }

                ~cell_worklist() {
                    delete[] arr;
                }

                int add_cell(work_set_int& _queue_pointer, int col) {
                    dej_assert(init);
                    dej_assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
                    _queue_pointer.set(col, cur_pos);
                    arr[cur_pos] = col;
                    cur_pos++;
                    return 0;
                }

                int next_cell(work_set_int& _queue_pointer, coloring *c) {
                    // look at first 12 positions and pick the (first) smallest cell within these entries
                    int sm_j = cur_pos - 1; /*< AKA "smallest j" */
                    for (int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
                        if (c->ptn[arr[j]] < c->ptn[arr[sm_j]]) {
                            sm_j = j;
                            if (c->ptn[arr[sm_j]] == 0)
                                break;
                        }
                    }

                    // swap sm_j and j
                    const int sm_col = arr[sm_j];
                    arr[sm_j] = arr[cur_pos - 1];
                    _queue_pointer.set(arr[sm_j], sm_j);

                    cur_pos--;
                    _queue_pointer.set(sm_col, -1);
                    return sm_col;
                }

                int next_cell(work_set_int& _queue_pointer, coloring *c, worklist_t<int>& _singleton_hint) {
                    // use singleton_hint
                    int sm_j = -1; /*< AKA "smallest j" */
                    while (!_singleton_hint.empty() && sm_j == -1) {
                        const int next_hint = _singleton_hint.pop_back();
                        sm_j = _queue_pointer.get(next_hint);
                    }

                    // look at first 12 positions and pick the (first) smallest cell within these entries
                    if (sm_j == -1) {
                        sm_j = cur_pos - 1;
                        for (int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
                            const int size_sm_j = c->ptn[arr[sm_j]];
                            const bool smaller = (c->ptn[arr[j]] < size_sm_j);
                            sm_j = smaller ? j : sm_j;
                            if (size_sm_j - smaller <= 0) break;
                        }
                    }

                    // swap sm_j and j
                    const int sm_col = arr[sm_j];
                    arr[sm_j] = arr[cur_pos - 1];
                    _queue_pointer.set(arr[sm_j], sm_j);

                    cur_pos--;
                    _queue_pointer.set(sm_col, -1);
                    return sm_col;
                }

                void replace_cell(work_set_int& _queue_pointer, int col_old, int col) {
                    const int pos = _queue_pointer.get(col_old);
                    arr[pos] = col;
                    dej_assert(_queue_pointer.get(col_old) != -1);
                    _queue_pointer.set(col_old, -1);
                    _queue_pointer.set(col, pos);
                }

                void reset(work_set_int& _queue_pointer) {
                    while (cur_pos > 0) {
                        cur_pos--;
                        _queue_pointer.set(arr[cur_pos], -1);
                    }
                }

                dej_nodiscard bool empty() const {
                    return (cur_pos == 0);
                }

                dej_nodiscard int size() const {
                    return cur_pos;
                }

            private:
                int *arr = nullptr;
                int arr_sz = -1;
                int cur_pos = -1;
                bool init = false;
            };

            int domain_size = -1;

            // worklist of color refinement algorithm
            work_set_int     queue_pointer;
            cell_worklist    cell_todo;
            worklist_t<int>  singleton_hint;

            // helper data structures for color refinement
            markset        scratch_set;
            worklist_t<int> vertex_worklist;
            workset_t<int>  color_vertices_considered; // todo should use different datastructure, with n space not 2n
            workset_t<int>  neighbours;
            workset_t<int>  neighbour_sizes;
            worklist_t<int> old_color_classes;
            workspace       scratch;

            void assure_initialized(const sgraph *g) {
                if (g->v_size > domain_size) {
                    const int n = g->v_size;

                    vertex_worklist.allocate(n * 2);
                    singleton_hint.allocate(n);
                    old_color_classes.allocate(n);
                    neighbours.initialize(n);
                    neighbour_sizes.initialize(n);
                    queue_pointer.initialize(n);
                    color_vertices_considered.initialize(n);
                    scratch.resize(n);
                    scratch_set.initialize(n);
                    cell_todo.initialize(n * 2);
                    domain_size = n;
                }
            }

            void refine_color_class_sparse(sgraph *g, coloring *c, int color_class,
                                           int class_size) {
                // for all vertices of the color class...
                int i, j, cc, end_cc, largest_color_class_size, acc;
                int *vertex_to_lab = c->vertex_to_lab;
                int *lab = c->lab;
                int *ptn = c->ptn;
                int *vertex_to_col = c->vertex_to_col;

                cc = color_class; // iterate over color class

                old_color_classes.reset();
                neighbours.reset();
                color_vertices_considered.reset();

                end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = vertex_to_col[v];
                        if (ptn[col] == 0) {
                            continue;
                        }
                        neighbours.inc_nr(v);
                        if (neighbours.get(v) == 0) {
                            color_vertices_considered.inc_nr(col);
                            dej_assert(col + color_vertices_considered.get(col) < g->v_size);
                            scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                            if (color_vertices_considered.get(col) == 0) {
                                old_color_classes.push_back(col);
                            }
                        }
                    }
                    cc += 1;
                }

                // remove color classes from list if they don't split
                for(j = 0; j < old_color_classes.size(); ++j) {
                    const int _col = old_color_classes[j];
                    const int _col_sz = ptn[_col] + 1;
                    const int vcount = color_vertices_considered.get(_col) + 1;

                    if(vcount == _col_sz) {
                        const int v_first     = scratch[_col];
                        const int index_first = neighbours.get(v_first) + 1;
                        bool splits = false;

                        for (i = 1; i < vcount; ++i) {
                            const int v     = scratch[_col + i];
                            const int index = neighbours.get(v) + 1;
                            splits = index != index_first;
                            if(splits) break;
                        }

                        if(splits) continue;

                        i = 0;
                        while (i < vcount) {
                            const int v = scratch[_col + i];
                            neighbours.set(v, -1);
                            ++i;
                        }

                        color_vertices_considered.set(_col, -1);

                        old_color_classes[j] = old_color_classes[old_color_classes.size() - 1];
                        old_color_classes.pop_back();
                        continue;
                    }
                }

                // sort split color classes
                old_color_classes.sort();

                // now split the color classes according to neighbour count
                while (!old_color_classes.empty()) {
                    const int _col    = old_color_classes.pop_back();
                    const int _col_sz = ptn[_col] + 1;
                    const int vcount  = color_vertices_considered.get(_col) + 1;
                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    for (i = 0; i < vcount; ++i) {
                        const int v     = scratch[_col + i];
                        const int index = neighbours.get(v) + 1;
                        if (neighbour_sizes.inc(index) == 0) vertex_worklist.push_back(index);
                    }

                    if(vcount == _col_sz && vertex_worklist.cur_pos == 1) {
                        vertex_worklist.reset();
                        j = 0;
                        while (j < vcount) {
                            const int v = scratch[_col + j];
                            neighbours.set(v, -1);
                            ++j;
                        }
                        color_vertices_considered.set(_col, -1);
                        continue;
                    }

                    vertex_worklist.sort();

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _ncol = _col + _col_sz - (neighbour_sizes.get(k));
                            if (_ncol != _col)
                                ptn[_ncol] = -1;
                        }
                    }


                    vertex_worklist.reset();
                    j = 0;
                    color_vertices_considered.set(_col, -1);

                    // determine colors and rearrange vertices
                    while (j < vcount) {
                        const int v = scratch[_col + j];
                        ++j;
                        if ((neighbours.get(v) == -1))
                            continue;
                        const int v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
                        neighbours.set(v, -1);
                        if (v_new_color == _col)
                            continue;

                        const int lab_pt = v_new_color + ptn[v_new_color] + 1;
                        ptn[v_new_color] += 1;
                        ptn[_col] -= (_col != v_new_color);

                        const int vertex_old_pos = vertex_to_lab[v];
                        const int vertex_at_pos = lab[lab_pt];
                        lab[vertex_old_pos] = vertex_at_pos;
                        vertex_to_lab[vertex_at_pos] = vertex_old_pos;
                        lab[lab_pt] = v;
                        vertex_to_col[v] = v_new_color;
                        vertex_to_lab[v] = lab_pt;
                    }

                    // add new colors to worklist
                    largest_color_class_size = -1;
                    int largest_color_class = -1;
                    for (i = _col; i < _col + _col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (ptn[i] + 1);
                        largest_color_class_size = new_largest? (ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        i += ptn[i] + 1;
                    }
                    dej_assert(largest_color_class >= _col);
                    dej_assert(largest_color_class < _col + _col_sz);

                    for (i = _col; i < _col + _col_sz;) {
                        const int i_sz = ptn[i] + 1;
                        const bool is_largest = i == largest_color_class;
                        report_split_color_class(c, _col, i, i_sz, is_largest);
                        i += i_sz;
                    }
                }

                neighbour_sizes.reset();
                vertex_worklist.reset();
            }

            void refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size) {
                int i, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                scratch_set.reset();
                old_color_classes.reset();
                vertex_worklist.reset();

                const int end_cc = color_class + class_size;

                // for all vertices of the color class...
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v   = g->e[i];
                        const int col = c->vertex_to_col[v];
                        if (c->ptn[col] > 0) {
                            neighbours.inc(v); // want to use reset variant?
                            if (!scratch_set.get(col)) {
                                scratch_set.set(col);
                                old_color_classes.push_back(col);
                            }
                        }
                    }
                    cc += 1;
                }

                // remove color classes from list if they don't split
                for(int j = 0; j < old_color_classes.size(); ++j) {
                    const int _col = old_color_classes[j];
                    const int _col_sz = c->ptn[_col] + 1;

                    const int v_first     = c->lab[_col];
                    const int index_first = neighbours.get(v_first);
                    bool splits = false;

                    for (i = 1; i < _col_sz && !splits; ++i) {
                        const int v     = c->lab[_col + i];
                        const int index = neighbours.get(v);
                        splits = index != index_first;
                    }

                    if(splits) continue;

                    old_color_classes[j] = old_color_classes[old_color_classes.size() - 1];
                    old_color_classes.pop_back();
                }

                old_color_classes.sort();

                // for every cell to be split...
                while (!old_color_classes.empty()) {
                    const int col = old_color_classes.pop_back();
                    const int col_sz = c->ptn[col] + 1;
                    if (col_sz == 1) continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        int index = neighbours.get(v) + 1;
                        if (index == first_index) continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) continue;

                    //if(vertex_worklist.cur_pos < 8)
                    vertex_worklist.sort();
                    //else bucket_sort(vertex_worklist, scratch_set, vertex_worklist.max());

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int j = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(j) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(j, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(j));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    vertex_worklist.reset();

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;
                    int largest_color_class  = 0;
                    // determine largest class to throw away
                    for (i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                    dej_assert(largest_color_class >= col);
                    dej_assert(largest_color_class < col + col_sz);

                    // report splits
                    for (i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset();
            }

            void refine_color_class_dense_cell(sgraph *g, coloring *c, int color_class, int class_size) {
                int i, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                //neighbours.reset_hard();
                neighbours.reset();
                scratch_set.reset();
                old_color_classes.reset();
                vertex_worklist.reset();

                const int end_cc = color_class + class_size;

                // for all vertices of the color class...
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];
                        neighbours.inc(v);
                        scratch_set.set(col);
                    }
                    cc += 1;
                }

                // for every cell to be split...
                for(int _col = 0; _col < g->v_size;) {
                    const int col_sz = c->ptn[_col] + 1;
                    const int col = _col;
                    _col += col_sz;
                    if (col_sz == 1 || !scratch_set.get(col)) continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        int index = neighbours.get(v) + 1;
                        if (index == first_index) continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) continue;

                    vertex_worklist.sort();
                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int j = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(j) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(j, val + acc);
                            acc += val;
                            const int rcol = col + col_sz - (neighbour_sizes.get(j));
                            c->ptn[rcol] = -1; // this is val - 1, actually...
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    //vertex_worklist.cur_pos = col_sz;
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;
                    int largest_color_class  = 0;
                    // determine largest class to throw away
                    for (i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                    dej_assert(largest_color_class >= col);
                    dej_assert(largest_color_class < col + col_sz);

                    // report splits
                    for (i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset();
            }

            void refine_color_class_dense_dense(sgraph *g, coloring *c, int color_class, int class_size) {
                int i, j, acc, cc, largest_color_class_size;
                int deg = -1;
                cc = color_class; // iterate over color class

                neighbours.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    deg = g->d[vc];
                    const int end_i = pe + deg - (deg==g->v_size-1?deg:0); // special code for universal vertices

                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        neighbours.inc_nr(v);
                    }
                    cc += 1;
                }

                if(class_size == 1 && deg == g->v_size-1) return; // special code for universal singletons

                // in dense-dense we dont sort, instead we just iterate over all cells (which are already sorted)
                // for every cell...
                for (j = 0; j < g->v_size;) {
                    const int col    = j;
                    const int col_sz = c->ptn[col] + 1;
                    j += col_sz;
                    if (col_sz == 1) continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index) continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) continue; // no split

                    vertex_worklist.sort();
                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    vertex_worklist.reset();

                    // copy cell for rearranging
                    memcpy(vertex_worklist.get_array(), c->lab + col, col_sz * sizeof(int));
                    vertex_worklist.cur_pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (!vertex_worklist.empty()) {
                        const int v = vertex_worklist.pop_back();
                    //for(i = 0; i < col_sz; ++i) {
                        //const int v = vertex_worklist[i];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;
                    int largest_color_class  =  0;

                    // determine largest class to throw away and finish
                    for (i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                    dej_assert(largest_color_class >= col);
                    dej_assert(largest_color_class < col + col_sz);

                    // report splits
                    for (i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset_hard();
            }

            void refine_color_class_singleton(sgraph *g, coloring *c, int color_class) {
                int i, cc, deg1_write_pos, deg1_read_pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                vertex_worklist.reset();
                old_color_classes.reset();

                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int deg= g->d[vc];
                const int end_i = pe + deg;

                for (i = pe; i < end_i; i++) {
                    const int v   = g->e[i];
                    const int col = c->vertex_to_col[v];

                    // singletons can not be split, so let's just ignore them
                    if (c->ptn[col] == 0) continue;

                    if(neighbours.get(col) == -1) {
                        old_color_classes.push_back(col);
                        // neighbours acts as the write position for degree 1 vertices of col in the scratchpad
                        neighbours.set(col, col);
                    }
                    // write vertex to scratchpad for later use
                    scratch[neighbours.get(col)] = v;
                    neighbours.inc_nr(col); // we reset neighbours later, use old_color_classes for reset
                }

                /*
                // remove color classes from list if they don't split
                const int pre_size = old_color_classes.size();
                for(int j = 0; j < old_color_classes.size();) {
                    const int deg0_col    = old_color_classes[j];
                    const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                    const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                    const int deg1_col = deg0_col + deg0_col_sz;

                    if (deg0_col == deg1_col) {
                        neighbours.set(deg1_col, -1);
                        old_color_classes[j] = old_color_classes[old_color_classes.size() - 1];
                        old_color_classes.pop_back();
                        continue;
                    }

                    ++j;
                }*/

                old_color_classes.sort();

                for(int j = 0; j < old_color_classes.cur_pos; ++j) {
                    const int deg0_col = old_color_classes[j];
                    const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                    const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                    const int deg1_col = deg0_col + deg0_col_sz;

                    dej_assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);

                    // no split? done...
                    if (deg0_col == deg1_col) {
                        neighbours.set(deg1_col, -1);
                        dej_assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                        continue;
                    }

                    dej_assert(deg0_col_sz + deg1_col_sz - 1 == c->ptn[deg0_col]);

                    // set ptn
                    c->ptn[deg0_col] = deg0_col_sz - 1;
                    c->ptn[deg1_col] = deg1_col_sz - 1;
                    c->ptn[deg1_col - 1] = 0;

                    deg1_write_pos = deg1_col;
                    deg1_read_pos = neighbours.get(deg0_col) - 1;

                    //c->vertex_to_col[c->lab[deg1_col]] = deg1_col;

                    // rearrange vertices of deg1 to the back of deg0 color
                    dej_assert(deg1_read_pos >= deg0_col);

                    while (deg1_read_pos >= deg0_col) {
                        const int v = scratch[deg1_read_pos];
                        const int vertex_at_pos = c->lab[deg1_write_pos];
                        const int lab_pos = c->vertex_to_lab[v];

                        c->lab[deg1_write_pos] = v;
                        c->vertex_to_lab[v] = deg1_write_pos;
                        c->vertex_to_col[v] = deg1_col;

                        c->lab[lab_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = lab_pos;

                        deg1_write_pos++;
                        deg1_read_pos--;
                    }

                    dej_assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                    dej_assert(c->vertex_to_col[c->lab[deg1_col]] == deg1_col);

                    // add new classes to report_splits
                    const bool leq = deg1_col_sz > deg0_col_sz;

                    report_split_color_class(c, deg0_col, deg0_col, deg0_col_sz, !leq);
                    report_split_color_class(c, deg0_col, deg1_col, deg1_col_sz, leq);

                    // reset neighbours count to -1
                    neighbours.set(deg0_col, -1);
                }

                old_color_classes.reset();
            }

            void refine_color_class_singleton_first(sgraph *g, coloring *c, int color_class) {
                int i, cc, deg1_write_pos, deg1_read_pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                scratch_set.reset();
                vertex_worklist.reset();
                old_color_classes.reset();

                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int end_i = pe + g->d[vc];
                for (i = pe; i < end_i; i++) {
                    const int v = g->e[i];
                    const int col = c->vertex_to_col[v];

                    if (c->ptn[col] == 0)
                        continue;

                    if (!scratch_set.get(col)) {
                        scratch_set.set(col);
                        old_color_classes.push_back(col);
                        // neighbours acts as the write position for degree 1 vertices of col
                        // in the singleton scratchpad
                        neighbours.set(col, col);
                    }
                    // write vertex to singleton scratchpad for later use
                    scratch[neighbours.get(col)] = v;
                    neighbours.inc_nr(col); // we reset neighbours later...
                    // old_color_classes can be used as reset information
                }

                for(int j = 0; j < old_color_classes.cur_pos; ++j) {
                    const int deg0_col = old_color_classes[j];
                    const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                    const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                    const int deg1_col = deg0_col + deg0_col_sz;

                    // no split? done...
                    if (deg0_col == deg1_col) {
                        neighbours.set(deg1_col, -1);
                        continue;
                    }

                    // set ptn
                    c->ptn[deg0_col] = deg0_col_sz - 1;
                    c->ptn[deg1_col] = deg1_col_sz - 1;
                    c->ptn[deg1_col - 1] = 0;

                    deg1_write_pos = deg1_col;
                    deg1_read_pos = neighbours.get(deg0_col) - 1;

                    // rearrange vertices of deg1 to the back of deg0 color
                    while (deg1_read_pos >= deg0_col) {
                        const int v = scratch[deg1_read_pos];
                        const int vertex_at_pos = c->lab[deg1_write_pos];
                        const int lab_pos = c->vertex_to_lab[v];

                        c->lab[deg1_write_pos] = v;
                        c->vertex_to_lab[v] = deg1_write_pos;
                        c->vertex_to_col[v] = deg1_col;
                        c->lab[lab_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = lab_pos;

                        deg1_write_pos++;
                        deg1_read_pos--;
                    }

                    // add new classes to report_splits
                    const bool leq = deg1_col_sz > deg0_col_sz;

                    report_split_color_class_first(c, deg0_col, deg0_col, deg0_col_sz, !leq);
                    report_split_color_class_first(c, deg0_col, deg1_col, deg1_col_sz, leq);

                    // reset neighbours count to -1
                    neighbours.set(deg0_col, -1);
                }

                old_color_classes.reset();
            }

            void refine_color_class_dense_first(sgraph *g, coloring *c, int color_class, int class_size) {
                int i, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                scratch_set.reset();
                old_color_classes.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];

                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];
                        if (c->ptn[col] == 0)
                            continue;

                        neighbours.inc(v);
                        if (!scratch_set.get(col)) {
                            scratch_set.set(col);
                            old_color_classes.push_back(col);
                        }
                    }
                    cc += 1;
                }

                while (!old_color_classes.empty()) {
                    const int col = old_color_classes.pop_back();
                    const int col_sz = c->ptn[col] + 1;

                    if (col_sz == 1) {
                        const int v_degree = neighbours.get(c->lab[col]) + 1;
                        if (v_degree == -1)
                            continue;
                    }

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) {
                        continue;
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;
                    int largest_color_class      = -1;

                    // determine largest class to throw away
                    for (i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }

                    // report splits
                    for (i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class_first(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }
            }

            void refine_color_class_dense_dense_first(sgraph *g, coloring *c, int color_class, int class_size) {
                // for all vertices of the color class...
                int i, j, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                old_color_classes.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        neighbours.inc_nr(v);
                    }
                    cc += 1;
                }

                // dont sort, just iterate over all cells
                for (j = 0; j < g->v_size;) {
                    const int col = j;
                    const int c_ptn = c->ptn[j];
                    j += c_ptn + 1;
                    const int col_sz = c_ptn + 1;

                    if (col_sz == 1)
                        continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) {
                        continue;
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1;
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;
                    int largest_color_class  =  0;

                    // determine largest class to throw away and finish
                    for (i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }

                    // report splits
                    for (i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class_first(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset_hard();
            }

            // distance heuristic inspired by the heuristic of nauty
            /*void refine_color_class_singleton_distance_first(sgraph *g, coloring *c, int color_class) {
                // for all vertices of the color class...

                // worklist arrays
                vertex_worklist.reset();
                color_vertices_considered.reset();
                neighbours.reset(); // storing distances here
                neighbour_sizes.reset();

                // origin of breadth-first search
                const int v = c->lab[color_class];
                vertex_worklist.push_back(v);
                neighbours.set(v, 0);
                neighbour_sizes.set(v, 0);

                // breadth-first search through graph to determine distances
                for(int j = 0; j < vertex_worklist.size(); ++j) {
                    const int v_current     = vertex_worklist[j];
                    const int v_current_col = c->vertex_to_col[v_current];
                    color_vertices_considered.set(v_current, 1);
                    const int hash_current = neighbours.get(v_current);
                    const int dis_current  = neighbour_sizes.get(v_current);
                    const int pe = g->v[v_current];
                    const int end_i = pe + g->d[v_current];

                    for (int i = pe; i < end_i; i++) {
                        const int v_next = g->e[i];
                        const int state  = color_vertices_considered.get(v_next);
                        if(state == -1) {
                            color_vertices_considered.set(v_next, 0);
                            vertex_worklist.push_back(v_next);
                            neighbours.set(v_next, hash_current + 1);
                            neighbour_sizes.set(v_next, dis_current + 1);
                        }

                        const int dis_next = neighbour_sizes.get(v_next);
                        if(dis_current != dis_next) {
                            neighbours.set(v_next, neighbours.get(v_next) + v_current_col + hash_current + dis_current);
                        }
                    }
                }

                neighbour_sizes.reset_hard();

                for(int j = 0; j < g->v_size; j += 1) {
                    neighbours.set(j, abs(neighbours.get(j)) % g->v_size);
                    dej_assert(neighbours.get(j) >= 0);
                }

                // iterate all cells
                for (int j = 0; j < g->v_size;) {
                    const int col = j;
                    const int c_ptn = c->ptn[j];
                    j += c_ptn + 1;
                    const int col_sz = c_ptn + 1;

                    if (col_sz == 1) continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]);
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (int i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v);
                        if (index == first_index) continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0) vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) continue;

                    // enrich neighbour_sizes to accumulative counting array
                    int acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1;
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    int pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - neighbour_sizes.get(neighbours.get(v));
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    int largest_color_class_size = -1;
                    int largest_color_class  =  0;

                    // determine largest class to throw away and finish
                    for (int i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }

                    // report splits
                    for (int i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class_first(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset_hard();
                color_vertices_considered.reset_hard();
            }*/

            // distance heuristic inspired by the heuristic of nauty
            /*void refine_color_class_singleton_distance(sgraph *g, coloring *c, int color_class) {
                // for all vertices of the color class...

                // worklist arrays
                vertex_worklist.reset();
                color_vertices_considered.reset();
                neighbours.reset(); // storing distances here
                neighbour_sizes.reset();

                // origin of breadth-first search
                const int v = c->lab[color_class];
                vertex_worklist.push_back(v);
                neighbours.set(v, 0);
                neighbour_sizes.set(v, 0);

                // breadth-first search through graph to determine distances
                for(int j = 0; j < vertex_worklist.size(); ++j) {
                    const int v_current     = vertex_worklist[j];
                    const int v_current_col = c->vertex_to_col[v_current];
                    color_vertices_considered.set(v_current, 1);
                    const int hash_current = neighbours.get(v_current);
                    const int dis_current  = neighbour_sizes.get(v_current);
                    const int pe = g->v[v_current];
                    const int end_i = pe + g->d[v_current];

                    for (int i = pe; i < end_i; i++) {
                        const int v_next = g->e[i];
                        const int state  = color_vertices_considered.get(v_next);
                        if(state == -1) {
                            color_vertices_considered.set(v_next, 0);
                            vertex_worklist.push_back(v_next);
                            neighbours.set(v_next, hash_current + 1);
                            neighbour_sizes.set(v_next, dis_current + 1);
                        }

                        const int dis_next = neighbour_sizes.get(v_next);
                        if(dis_current != dis_next) {
                            neighbours.set(v_next, neighbours.get(v_next) + v_current_col + hash_current + dis_current);
                        }
                    }
                }

                neighbour_sizes.reset_hard();

                for(int j = 0; j < g->v_size; j += 1) {
                    neighbours.set(j, abs(neighbours.get(j)) % g->v_size);
                    dej_assert(neighbours.get(j) >= 0);
                }

                // iterate all cells
                for (int j = 0; j < g->v_size;) {
                    const int col = j;
                    const int c_ptn = c->ptn[j];
                    j += c_ptn + 1;
                    const int col_sz = c_ptn + 1;

                    if (col_sz == 1) continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]);
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (int i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v);
                        if (index == first_index) continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0) vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) continue;

                    // enrich neighbour_sizes to accumulative counting array
                    int acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1;
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch.get_array(), c->lab + col, col_sz * sizeof(int));
                    int pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - neighbour_sizes.get(neighbours.get(v));
                        dej_assert(v_new_color >= col);
                        dej_assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    int largest_color_class_size = -1;
                    int largest_color_class  =  0;

                    // determine largest class to throw away and finish
                    for (int i = col; i < col + col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }

                    // report splits
                    for (int i = col; i < col + col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class(c, col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbours.reset_hard();
                color_vertices_considered.reset_hard();
            }*/

            void refine_color_class_sparse_first(sgraph *g, coloring *c, int color_class, int class_size) {
                int v_new_color, cc, largest_color_class_size, acc;
                cc = color_class; // iterate over color class

                scratch_set.reset();
                old_color_classes.reset();
                neighbours.reset();
                color_vertices_considered.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (int i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];

                        if (c->ptn[col] == 0)
                            continue;

                        neighbours.inc_nr(v);
                        if (neighbours.get(v) == 0) {
                            color_vertices_considered.inc_nr(col);
                            dej_assert(col + color_vertices_considered.get(col) < g->v_size);
                            scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                            if (!scratch_set.get(col)) {
                                old_color_classes.push_back(col);
                                scratch_set.set(col);
                            }
                        }
                    }
                    cc += 1;
                }

                // split color classes according to neighbour count
                while (!old_color_classes.empty()) {
                    const int _col = old_color_classes.pop_back();
                    const int _col_sz = c->ptn[_col] + 1;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    for (int i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                        const int v = scratch[_col + i];
                        int index = neighbours.get(v) + 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int i = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(i) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(i, val + acc);
                            acc += val;
                        }
                    }

                    // determine colors
                    vertex_worklist.reset();

                    for (int i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                        const int v = scratch[_col + i];
                        if (neighbours.get(v) == -1) {
                            v_new_color = _col;
                            dej_assert(false);
                        } else {
                            dej_assert((neighbours.get(v) + 1) > 0);
                            v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
                        }
                        if (v_new_color != _col) {
                            vertex_worklist.push_back(v_new_color);
                            vertex_worklist.push_back(v);
                            c->ptn[v_new_color] = -1;
                        }
                        neighbours.set(v, -1);
                    }

                    color_vertices_considered.set(_col, -1);

                    // rearrange vertices
                    while (!vertex_worklist.empty()) {
                        const int v = vertex_worklist.pop_back();
                        const int v_new_color2 = vertex_worklist.pop_back();

                        const int vertex_old_pos = c->vertex_to_lab[v];
                        const int vertex_at_pos = c->lab[v_new_color2 + c->ptn[v_new_color2] + 1];
                        c->lab[vertex_old_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

                        c->lab[v_new_color2 + c->ptn[v_new_color2] + 1] = v;
                        c->vertex_to_col[v] = v_new_color2;
                        c->vertex_to_lab[v] = v_new_color2 + c->ptn[v_new_color2] + 1;
                        c->ptn[v_new_color2] += 1;

                        if (_col != v_new_color2) {
                            dej_assert(v_new_color2 > _col);
                            c->ptn[_col] -= 1;
                        } else
                            dej_assert(false);
                    }

                    // add new colors to worklist
                    largest_color_class_size = -1;
                    int largest_color_class      = -1;
                    int i;
                    for (i = _col; i < _col + _col_sz;) {
                        dej_assert(i >= 0 && i < c->domain_size);
                        dej_assert(c->ptn[i] + 1 > 0);
                        const bool new_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = new_largest? (c->ptn[i] + 1) : largest_color_class_size;
                        largest_color_class = new_largest? i : largest_color_class;
                        if (i != 0) c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                    if (i != 0) c->ptn[i - 1] = 0;


                    for (i = _col; i < _col + _col_sz;) {
                        const int i_sz = c->ptn[i] + 1;
                        report_split_color_class_first(c, _col, i, i_sz, i == largest_color_class);
                        i += i_sz;
                    }
                }

                neighbour_sizes.reset();
                vertex_worklist.reset();
            }
        };
    }
}

#endif //DEJAVU_REFINEMENT_H
