#ifndef SASSY_REFINEMENT_H
#define SASSY_REFINEMENT_H

#include "utility.h"
#include "coloring.h"
#include "sgraph.h"
#include "invariant.h"
#include "configuration.h"
#include <list>
#include <iostream>
#include <cstring>

namespace sassy {
// sorting utilizing minimal sorting networks for n <= 6
    template<class T>
    void sort_t(T *arr, int sz) {
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x, y) { const T a = min(arr[x], arr[y]); \
                    const T b = max(arr[x], arr[y]); \
                    arr[x] = a; arr[y] = b; }
        switch (sz) {
            case 0:
            case 1:
                break;
            case 2: SWAP(0, 1);
                break;
            case 3: SWAP(0, 1);
                SWAP(0, 2);
                SWAP(1, 2);
                break;
            case 4: SWAP(0, 1);
                SWAP(2, 3);
                SWAP(0, 2);
                SWAP(1, 3);
                SWAP(1, 2);
                break;
            case 5: SWAP(0, 1);
                SWAP(2, 3);
                SWAP(0, 2);
                SWAP(1, 4);
                SWAP(0, 1);
                SWAP(2, 3);
                SWAP(1, 2);
                SWAP(3, 4);
                SWAP(2, 3);
                break;
            case 6: SWAP(1, 2);
                SWAP(4, 5);
                SWAP(0, 2);
                SWAP(3, 5);
                SWAP(0, 1);
                SWAP(3, 4);
                SWAP(1, 4);
                SWAP(0, 3);
                SWAP(2, 5);
                SWAP(1, 3);
                SWAP(2, 4);
                SWAP(2, 3);
                break;
            default:
                std::sort(arr, arr + sz);
        }
#undef SWAP
#undef min
#undef max
    }

// work list / stack with fixed size limitation
    template<class T>
    class work_list_t {
    public:
        work_list_t() {};
        work_list_t(int size) {
            initialize(size);
        };

        void initialize(int size) {
            if(init)
                delete[] arr;
            arr = new T[size];
            arr_sz = size;
            init = true;
            cur_pos = 0;
        }

        void push_back(T value) {
            assert(cur_pos >= 0 && cur_pos < arr_sz);
            arr[cur_pos] = value;
            cur_pos += 1;
        }

        T pop_back() {
            return arr[--cur_pos];
        }

        T *last() {
            return &arr[cur_pos - 1];
        }

        bool empty() {
            return cur_pos == 0;
        }

        void reset() {
            cur_pos = 0;
        }

        ~work_list_t() {
            if (init) {
                delete[] arr;
            }
        }

        void sort() {
            sort_t<T>(arr, cur_pos);
        }

        T *get_array() {
            return arr;
        }

        T &operator[](int index) {
            assert(index >= 0);
            assert(index < arr_sz);
            return arr[index];
        }

        void sort_after_map(T *map) {
            struct comparator_map {
                T *map;

                comparator_map(T *map) {
                    this->map = map;
                }

                bool operator()(const T &a, const T &b) {
                    return map[a] < map[b];
                }
            };
            comparator_map c = comparator_map(map);
            std::sort(arr, arr + cur_pos, c);
        }

        int cur_pos = 0;
        bool init = false;
    private:
        T *arr = nullptr;
        int arr_sz = -1;
    };

    typedef work_list_t<std::pair<std::pair<int, int>, bool>> work_list_pair_bool;
    typedef work_list_t<int> work_list;

    class tiny_orbit {
        int sz;
        mark_set touched;
        work_list_t<int> reset_arr;
        work_list_t<int> map_arr;
    public:
        int find_and_cut_orbit(const int v1) {
            assert(v1 >= 0);
            assert(v1 < sz);
            int orbit1 = map_arr[v1];
            while (orbit1 != map_arr[orbit1])
                orbit1 = map_arr[orbit1];
            map_arr[v1] = orbit1;
            return orbit1;
        }

        void combine_orbits(const int v1, const int v2) {
            assert(v1 >= 0);
            assert(v2 >= 0);
            assert(v1 < sz);
            assert(v2 < sz);
            if (v1 != v2) {
                if (!touched.get(v1))
                    reset_arr.push_back(v1);
                if (!touched.get(v2))
                    reset_arr.push_back(v2);
                touched.set(v1);
                touched.set(v2);
                int orbit1 = find_and_cut_orbit(v1);
                int orbit2 = find_and_cut_orbit(v2);
                if (orbit1 < orbit2) {
                    map_arr[orbit2] = orbit1;
                } else {
                    map_arr[orbit1] = orbit2;
                }
            }
        }

        bool are_in_same_orbit(const int v1, const int v2) {
            assert(v1 >= 0);
            assert(v2 >= 0);
            assert(v1 < sz);
            assert(v2 < sz);
            if (v1 == v2)
                return true;
            int orbit1 = find_and_cut_orbit(v1);
            int orbit2 = find_and_cut_orbit(v2);
            return (orbit1 == orbit2);
        }

        void reset() {
            while (!reset_arr.empty()) {
                const int v = reset_arr.pop_back();
                map_arr[v] = v;
            }
            touched.reset();
        }

        void initialize(int domain_size) {
            sz = domain_size;
            touched.initialize(domain_size);
            reset_arr.initialize(domain_size);
            map_arr.initialize(domain_size);
            for (int i = 0; i < domain_size; ++i) {
                map_arr.push_back(i);
            }
        }
    };

// queue with fixed size limitation
    class work_queue {
    public:
        void initialize(int size) {
            assert(!init);
            sz = size;
            pos = 0;
            queue = new int[size];
            init = true;
        }

        void push(int val) {
            //assert(init);
            assert(pos != sz);
            queue[pos] = val;
            pos++;
        }

        int pop() {
            //assert(init);
            assert(pos > 0);
            pos--;
            return queue[pos];
        }

        bool empty() {
            //assert(init);
            return (pos == 0);
        }

        ~work_queue() {
            if (init) {
                delete[] queue;
            }
        }

        void reset() {
            pos = 0;
        }

        int *queue;
        int pos;
    private:
        bool init = false;
        int sz;
    };

// work set with arbitrary type
    template<class T>
    class work_set_t {
    public:
        void initialize(int size) {
            s = new T[size];
            reset_queue.initialize(size);

            memset(s, -1, size * sizeof(T));

            init = true;
            sz = size;
        }

        void set(int index, T value) {
            assert(index >= 0);
            assert(index < sz);
            s[index] = value;
        }

        T get(int index) {
            assert(index >= 0);
            assert(index < sz);
            return s[index];
        }

        void reset() {
            while (!reset_queue.empty())
                s[reset_queue.pop()] = -1;
        }

        void reset_hard() {
            memset(s, -1, sz * sizeof(T));
            reset_queue.pos = 0;
        }

        T inc(int index) {
            assert(index >= 0);
            assert(index < sz);
            if (s[index]++ == -1)
                reset_queue.push(index);
            return s[index];
        }

        void inc_nr(int index) {
            assert(index >= 0 && index < sz);
            ++s[index];
        }

        ~work_set_t() {
            if (init)
                delete[] s;
        }

        work_queue reset_queue;
    private:
        bool init = false;
        T *s;
        int sz;
    };

    typedef work_set_t<int> work_set_int;


// worklist implementation for color refinement
    class cell_worklist {
    public:
        void initialize(int domain_size) {
            arr = new int[domain_size];
            arr_sz = domain_size;
            init = true;
            cur_pos = 0;
        }

        ~cell_worklist() {
            delete[] arr;
        }

        int add_cell(work_set_int *queue_pointer, int col) {
            assert(init);
            assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
            queue_pointer->set(col, cur_pos);
            arr[cur_pos] = col;
            cur_pos++;
            return 0;
        }

        int next_cell(work_set_int *queue_pointer, coloring *c) {
            // look at first 12 positions and pick the (first) smallest cell within these entries
            int sm_j = cur_pos - 1;
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
            queue_pointer->set(arr[sm_j], sm_j);

            cur_pos--;
            queue_pointer->set(sm_col, -1);
            return sm_col;
        }

        int next_cell(work_set_int *queue_pointer, coloring *c, work_list_t<int> *singleton_hint) {
            // use singleton_hint
            int sm_j = -1;
            while (!singleton_hint->empty()) {
                const int next_hint = singleton_hint->pop_back();
                sm_j = queue_pointer->get(next_hint);
                if (sm_j == -1)
                    continue;
                else
                    break;
            }

            // look at first 4 positions and pick the (first) smallest cell within these entries
            if (sm_j == -1) {
                sm_j = cur_pos - 1;
                for (int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 4); --j) {
                    //bool smaller_d = g->d[arr[j]] > g->d[arr[sm_j]];
                    const int size_sm_j = c->ptn[arr[sm_j]];
                    if (c->ptn[size_sm_j] == 0)
                        break;
                    const bool smaller = (c->ptn[arr[j]] < size_sm_j);
                    //bool eq        = (c->ptn[arr[j]] == c->ptn[arr[sm_j]]);
                    sm_j = smaller ? j : sm_j;
                }
            }

            // swap sm_j and j
            const int sm_col = arr[sm_j];
            arr[sm_j] = arr[cur_pos - 1];
            queue_pointer->set(arr[sm_j], sm_j);

            cur_pos--;
            queue_pointer->set(sm_col, -1);
            return sm_col;
        }

        void replace_cell(work_set_int *queue_pointer, int col_old, int col) {
            const int pos = queue_pointer->get(col_old);
            arr[pos] = col;
            assert(queue_pointer->get(col_old) != -1);
            queue_pointer->set(col_old, -1);
            queue_pointer->set(col, pos);
        }

        void reset(work_set_int *queue_pointer) {
            while (cur_pos > 0) {
                cur_pos--;
                queue_pointer->set(arr[cur_pos], -1);
            }
        }

        bool empty() {
            return (cur_pos == 0);
        }

        int size() {
            return cur_pos;
        }

    private:
        int *arr = nullptr;
        int arr_sz = -1;
        int cur_pos = -1;
        bool init = false;
    };

// refinement manager, preserving the workspace between refinements
    class refinement {
    private:
        configstruct* config = nullptr;
    public:
        refinement()=delete;
        refinement(configstruct* config) {
            this->config = config;
        }

        // color refinement
        // includes several options for using invariants, blueprints and k-deviation
        bool refine_coloring(sgraph *g, coloring *c, invariant *I,
                             int init_color_class, strategy_metrics *m, int cell_early, int individualize_early,
                             std::vector<int> *early_individualized, mark_set *touched_color,
                             work_list *touched_color_list) {
            bool comp = true;
            singleton_hint.reset();
            assure_initialized(g);
            int deviation_expander = 0;

            cell_todo.reset(&queue_pointer);

            if (init_color_class < 0) {
                // initialize queue with all classes (except for largest one)
                for (int i = 0; i < c->ptn_sz;) {
                    cell_todo.add_cell(&queue_pointer, i);
                    const int col_sz = c->ptn[i];
                    if (col_sz == 0) {
                        //singleton_hint.push_back(i);
                    }
                    i += col_sz + 1;

                }
            } else {
                const int col_sz = c->ptn[init_color_class];
                assert(c->vertex_to_col[c->lab[init_color_class]] == init_color_class);
                cell_todo.add_cell(&queue_pointer, init_color_class);
                if (col_sz == 0) {
                    //singleton_hint.push_back(init_color_class);
                }
            }
            int its = 0;

            while (!cell_todo.empty()) {
                its += 1;
                color_class_splits.reset();
                const int next_color_class = cell_todo.next_cell(&queue_pointer, c);
                const int next_color_class_sz = c->ptn[next_color_class] + 1;
                if (m)
                    m->color_refinement_cost += next_color_class_sz;
                comp = I->write_top_and_compare(INV_MARK_STARTCELL, true) && comp;

                // if cell did not split anything in the target invariant, skip refinement until the end of this cell
                if (I->no_write && !I->never_fail && comp) {
                    const bool skip = !I->protocol_read(next_color_class);
                    if (skip && config->CONFIG_IR_IDLE_SKIP) {
                        I->fast_forward(INV_MARK_ENDCELL);
                        continue;
                    }
                }

                bool dense_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

                const bool pre_comp = comp;

                if (next_color_class_sz == 1 && !(g->dense && dense_dense)) {
                    // singleton
                    comp = refine_color_class_singleton(g, c, next_color_class, next_color_class_sz,
                                                        &color_class_splits, I);
                } else if (g->dense) {
                    if (dense_dense) { // dense-dense
                        comp = refine_color_class_dense_dense(g, c, next_color_class, next_color_class_sz,
                                                              &color_class_splits, I);
                    } else { // dense-sparse
                        comp = refine_color_class_dense(g, c, next_color_class, next_color_class_sz,
                                                        &color_class_splits, I);
                    }
                } else { // sparse
                    comp = refine_color_class_sparse(g, c, next_color_class, next_color_class_sz,
                                                     &color_class_splits, I);
                }

                comp = comp && pre_comp;

                deviation_expander -= (!comp);
                if (!comp && deviation_expander <= 0) {
                    break;
                }

                // add all new classes except for the first, largest one
                int skip = 0;

                int latest_old_class = -1;
                bool skipped_largest = false;
                const int pre_cells = c->cells;

                // color class splits are sorted in reverse
                // the old color class will always come last
                while (!color_class_splits.empty()) {
                    int old_class = color_class_splits.last()->first.first;
                    int new_class = color_class_splits.last()->first.second;
                    bool is_largest = color_class_splits.last()->second;


                    // record colors that were changed
                    if (touched_color) {
                        if (!touched_color->get(new_class)) {
                            touched_color->set(new_class);
                            touched_color_list->push_back(new_class);
                        }
                    }

                    c->cells += (old_class != new_class);

#ifndef NDEBUG // debug code
                    if (color_class_splits.empty()) {
                        int actual_cells = 0;
                        for (int i = 0; i < c->ptn_sz;) {
                            actual_cells += 1;
                            i += c->ptn[i] + 1;
                        }

                        assert(c->cells == actual_cells);
                    }
#endif

                    // detection if coloring is discrete
                    if (c->cells == g->v_size) {
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        I->write_cells(c->cells);
                        I->protocol_write(true, next_color_class);
                        I->protocol_mark();
                        if (I->no_write)
                            I->protocol_skip_to_mark();
                        comp = I->write_top_and_compare(INV_MARK_ENDCELL, true) && comp;
                        comp = I->write_top_and_compare(INV_MARK_ENDREF, true) && comp;
                        return comp;
                    }

                    // partition is as large as the one of target invariant, can skip to the end of the entire refinement
                    if (c->cells == cell_early && comp && !config->CONFIG_IR_REFINE_EARLYOUT_LATE) {
                        if (!I->only_acc) {
                            I->fast_forward(INV_MARK_ENDREF);
                            if (I->no_write)
                                I->protocol_skip_to_mark();
                        }
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        return comp;
                    }

                    // random blueprint individualization during color refinement
                    if (config->CONFIG_IR_INDIVIDUALIZE_EARLY) {
                    }

                    if (latest_old_class != old_class) {
                        latest_old_class = old_class;
                        skipped_largest = false;
                    }

                    // management code for skipping largest class resulting from old_class
                    color_class_splits.pop_back();
                    // int new_class_sz = c->ptn[new_class] + 1;

                    if (skipped_largest || !is_largest) {
                        cell_todo.add_cell(&queue_pointer, new_class);
                        //if(new_class_sz == 1)
                        //singleton_hint.push_back(new_class);
                    } else {
                        skipped_largest = true;
                        skip += 1;

                        // since old color class will always appear last, the queue pointer of old color class is still valid!
                        int i = queue_pointer.get(old_class);
                        if (i >= 0) {
                            cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                            //if(new_class_sz == 1)
                            //    singleton_hint.push_back(new_class);
                        }
                    }
                }
                const int new_cells = c->cells - pre_cells;

                // mark end of cell and denote whether this cell was splitting or non-splitting
                I->protocol_write(new_cells > 0, next_color_class);
                comp = I->write_top_and_compare(INV_MARK_ENDCELL, true) && comp;
                if (!comp && deviation_expander <= 0) break;
            }

            if (comp) {
                if (I->no_write)
                    I->protocol_skip_to_mark();
                I->protocol_mark();
                I->write_cells(c->cells);
                comp = I->write_top_and_compare(INV_MARK_ENDREF, true) && comp;
            }

            return comp;
        }

        // individualize a vertex in a coloring
        int individualize_vertex(coloring *c, int v) {
            const int color = c->vertex_to_col[v];
            const int pos = c->vertex_to_lab[v];

            int color_class_size = c->ptn[color];

            assert(color_class_size > 0);

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
            return color + color_class_size;
        }

        // color refinement that does not produce an isomorphism-invariant partitioning, but uses more optimization
        // techniques -- meant to be used as the first refinement in automorphism computation
        bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class = -1) {
            assure_initialized(g);
            singleton_hint.reset();
            cell_todo.reset(&queue_pointer);

            if (init_color_class < 0) {
                for (int i = 0; i < c->ptn_sz;) {
                    cell_todo.add_cell(&queue_pointer, i);
                    const int col_sz = c->ptn[i];
                    if (col_sz == 0) {
                        singleton_hint.push_back(i);
                    }
                    i += col_sz + 1;
                }
            } else {
                cell_todo.add_cell(&queue_pointer, init_color_class);
            }

            while (!cell_todo.empty()) {
                color_class_splits.reset();
                const int next_color_class    = cell_todo.next_cell(&queue_pointer, c, &singleton_hint);
                //const int next_color_class    = cell_todo.next_cell(&queue_pointer, c);
                const int next_color_class_sz = c->ptn[next_color_class] + 1;
                bool dense_dense = false;
                //if(g->d[c->lab[next_color_class]] == 0) {
                //   continue;
                //}
                if (g->dense && g->d[c->lab[next_color_class]] > 5) {
                    dense_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));
                }

                if (next_color_class_sz == 1 && !(g->dense && dense_dense)) {
                    // singleton
                    refine_color_class_singleton_first(g, c, next_color_class, &color_class_splits);
                } else if (g->dense) {
                    if (dense_dense) { // dense-dense
                        refine_color_class_dense_dense_first(g, c, next_color_class, next_color_class_sz,
                                                             &color_class_splits);
                    } else { // dense-sparse
                        refine_color_class_dense_first(g, c, next_color_class, next_color_class_sz,
                                                       &color_class_splits);
                    }
                } else { // sparse
                    refine_color_class_sparse_first(g, c, next_color_class, next_color_class_sz,
                                                    &color_class_splits);
                }

                // add all new classes except for the first, largest one
                int latest_old_class = -1;
                bool skipped_largest = false;

                // color class splits are sorted in reverse
                // the old color class will always come last
                while (!color_class_splits.empty()) {
                    const int old_class   = color_class_splits.last()->first.first;
                    const int new_class   = color_class_splits.last()->first.second;
                    const bool is_largest = color_class_splits.last()->second;

                    c->cells += (old_class != new_class);

#ifndef NDEBUG // debug code
                    if (color_class_splits.empty()) {
                        int actual_cells = 0;
                        for (int i = 0; i < c->ptn_sz;) {
                            actual_cells += 1;
                            i += c->ptn[i] + 1;
                        }

                        assert(c->cells == actual_cells);
                    }
#endif

                    if (c->cells == g->v_size) {
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        return true;
                    }

                    if (latest_old_class != old_class) {
                        latest_old_class = old_class;
                        skipped_largest = false;
                    }

                    color_class_splits.pop_back();


                    if (skipped_largest || !is_largest) {
                        cell_todo.add_cell(&queue_pointer, new_class);
                        const int new_class_sz = c->ptn[new_class] + 1;
                        if (new_class_sz == 1) {
                            singleton_hint.push_back(new_class);
                        }
                    } else {
                        skipped_largest = true;

                        // since old color class will always appear last, the queue pointer of old color class is still valid!
                        int i = queue_pointer.get(old_class);
                        if (i >= 0) {
                            cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                            const int new_class_sz = c->ptn[new_class] + 1;
                            if (new_class_sz == 1)
                                singleton_hint.push_back(new_class);
                        }
                    }
                }
            }

            // assert(assert_is_equitable(g, c));
            return true;
        }

        // certify an automorphism on a graph
        bool certify_automorphism(sgraph *g, const int *p) {
            int i, found;

            assure_initialized(g);

            for (i = 0; i < g->v_size; ++i) {
                const int image_i = p[i];
                if (image_i == i)
                    continue;
                if (g->d[i] != g->d[image_i]) // degrees must be equal
                    return false;

                scratch_set.reset();
                // automorphism must preserve neighbours
                found = 0;
                for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                    const int vertex_j = g->e[j];
                    const int image_j = p[vertex_j];
                    scratch_set.set(image_j);
                    found += 1;
                }
                for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                    const int vertex_j = g->e[j];
                    if (!scratch_set.get(vertex_j)) {
                        return false;
                    }
                    scratch_set.unset(vertex_j);
                    found -= 1;
                }
                if (found != 0) {
                    return false;
                }
            }

            return true;
        }

        bool is_adjacent(sgraph*g, int vert1, int vert2) {
            if(g->d[vert1] < g->d[vert2]) {
                for(int i = 0; i < g->d[vert1]; ++i) {
                    if(g->e[g->v[vert1] + i] == vert2)
                        return true;
                }
                return false;
            } else {
                for(int i = 0; i < g->d[vert2]; ++i) {
                    if(g->e[g->v[vert2] + i] == vert1)
                        return true;
                }
                return false;
            }
        }

        // certify an automorphism on a graph
        bool certify_automorphism(sgraph *g, const int *colmap, const int *p) {
            int i, found;

            assure_initialized(g);

            for (i = 0; i < g->v_size; ++i) {
                const int image_i = p[i];
                if (image_i == i)
                    continue;
                if (g->d[i] != g->d[image_i]) { // degrees must be equal
                    return false;
                }
                if (colmap[i] != colmap[image_i]) { // colors must be equal
                    return false;
                }

                scratch_set.reset();
                // automorphism must preserve neighbours
                found = 0;
                for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                    const int vertex_j = g->e[j];
                    const int image_j = p[vertex_j];
                    scratch_set.set(image_j);
                    found += 1;
                }
                for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                    const int vertex_j = g->e[j];
                    if (!scratch_set.get(vertex_j)) {
                        return false;
                    }
                    scratch_set.unset(vertex_j);
                    found -= 1;
                }
                if (found != 0) {
                    return false;
                }
            }

            return true;
        }

        // certify an automorphism on a graph, sparse
        bool certify_automorphism_sparse(const sgraph *g, const int *colmap, const int *p,
                                         int supp, const int *supp_arr) {
            int i, found;

            assure_initialized(g);

            //for(i = 0; i < g->v_size; ++i) {
            for (int f = 0; f < supp; ++f) {
                i = supp_arr[f];
                const int image_i = p[i];
                if (image_i == i)
                    continue;
                if (g->d[i] != g->d[image_i]) // degrees must be equal
                    return false;
                if (colmap[i] != colmap[image_i]) // colors must be equal
                    return false;

                scratch_set.reset();
                // automorphism must preserve neighbours
                found = 0;
                for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                    const int vertex_j = g->e[j];
                    const int image_j = p[vertex_j];
                    if (colmap[vertex_j] != colmap[image_j])
                        return false;
                    scratch_set.set(image_j);
                    //scratch[image_j] = vertex_j;
                    found += 1;
                }
                for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                    const int vertex_j = g->e[j];
                    if (!scratch_set.get(vertex_j)) {
                        return false;
                    }
                    //if(colmap[scratch[vertex_j]] != colmap[vertex_j])
                    //    return false;
                    scratch_set.unset(vertex_j);
                    found -= 1;
                }
                if (found != 0) {
                    return false;
                }
            }

            return true;
        }

        // certify an automorphism on a graph, sparse, report on which vertex failed
        std::pair<bool, int>
        certify_automorphism_sparse_report_fail(const sgraph *g, const int *colmap,
                                                const int *p, int supp, const int *supp_arr) {
            int i, found;

            assure_initialized(g);

            //for(i = 0; i < g->v_size; ++i) {
            for (int f = 0; f < supp; ++f) {
                i = supp_arr[f];
                const int image_i = p[i];
                if (image_i == i)
                    continue;
                if (g->d[i] != g->d[image_i]) // degrees must be equal
                    return std::pair<bool, int>(false, -1);
                if (colmap[i] != colmap[image_i]) // colors must be equal
                    return std::pair<bool, int>(false, -1);

                scratch_set.reset();
                // automorphism must preserve neighbours
                found = 0;
                for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                    const int vertex_j = g->e[j];
                    const int image_j = p[vertex_j];
                    if (colmap[vertex_j] != colmap[image_j])
                        return std::pair<bool, int>(false, i);
                    scratch_set.set(image_j);
                    //scratch[image_j] = vertex_j;
                    found += 1;
                }
                for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                    const int vertex_j = g->e[j];
                    if (!scratch_set.get(vertex_j)) {
                        return std::pair<bool, int>(false, i);
                    }
                    scratch_set.unset(vertex_j);
                    found -= 1;
                }
                if (found != 0) {
                    return std::pair<bool, int>(false, i);
                }
            }

            return std::pair<bool, int>(true, -1);
        }

        // certify an automorphism, for a single vertex
        bool check_single_failure(const sgraph *g, const int *colmap, const int *p,
                                  int failure) {
            int i, found;

            assure_initialized(g);

            i = failure;
            const int image_i = p[i];
            if (image_i == i)
                return true;
            if (g->d[i] != g->d[image_i]) // degrees must be equal
                return false;
            if (colmap[i] != colmap[image_i]) // colors must be equal
                return false;

            scratch_set.reset();
            // automorphism must preserve neighbours
            found = 0;
            for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                const int vertex_j = g->e[j];
                const int image_j = p[vertex_j];
                if (colmap[vertex_j] != colmap[image_j])
                    return false;
                scratch_set.set(image_j);
                //scratch[image_j] = vertex_j;
                found += 1;
            }
            for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                const int vertex_j = g->e[j];
                if (!scratch_set.get(vertex_j)) {
                    return false;
                }
                scratch_set.unset(vertex_j);
                found -= 1;
            }
            if (found != 0) {
                return false;
            }

            return true;
        }

        ~refinement() {
            if (initialized)
                delete[] workspace_int;
        }

    private:
        bool initialized = false;
        work_set_int queue_pointer;
        cell_worklist cell_todo;
        mark_set scratch_set;
        work_list_t<int> vertex_worklist;
        work_set_t<int> color_vertices_considered;
        work_set_t<int> neighbours; // degree type instead?
        work_set_t<int> neighbour_sizes;
        //work_list_t<int>  singletons;
        work_list_t<int> singleton_hint;
        work_list_t<int> old_color_classes;
        work_list_pair_bool color_class_splits;

        int *scratch;
        int *workspace_int;

        void assure_initialized(const sgraph *g) {
            if (!initialized) {
                const int n = g->v_size;

                // reducing contention on heap allocator through bulk allocation...
                workspace_int = new int[n * 2];

                vertex_worklist.initialize(n * 2);
                //singletons.initialize(n);
                singleton_hint.initialize(n);
                old_color_classes.initialize(n);
                neighbours.initialize(n);
                neighbour_sizes.initialize(n);
                queue_pointer.initialize(n);
                color_vertices_considered.initialize(n);

                scratch = (int *) workspace_int;
                scratch_set.initialize_from_array(workspace_int + n, n);

                color_class_splits.initialize(n);
                cell_todo.initialize(n * 2);

                memset(scratch, 0, n * sizeof(int));
                initialized = true;
            }
        }

        bool refine_color_class_sparse(sgraph *g, coloring *c,
                                       int color_class, int class_size,
                                       work_list_pair_bool *color_class_split_worklist, invariant *I) {
            // for all vertices of the color class...
            bool comp, mark_as_largest;
            int i, j, cc, end_cc, largest_color_class_size, acc_in, singleton_inv1, singleton_inv2, acc;
            int *vertex_to_lab = c->vertex_to_lab;
            int *lab = c->lab;
            int *ptn = c->ptn;
            int *vertex_to_col = c->vertex_to_col;

            cc = color_class; // iterate over color class
            comp = true;

            //singletons.reset();
            scratch_set.reset();
            old_color_classes.reset();
            neighbours.reset();
            color_vertices_considered.reset();

            end_cc = color_class + class_size;
            acc_in = 0;
            singleton_inv1 = 0;
            singleton_inv2 = 0;
            while (cc < end_cc) { // increment value of neighbours of vc by 1
                const int vc = lab[cc];
                const int pe = g->v[vc];
                const int end_i = pe + g->d[vc];
                for (i = pe; i < end_i; i++) {
                    const int v = g->e[i];
                    const int col = vertex_to_col[v];
                    if (ptn[col] == 0) {
                        singleton_inv1 += MASH2(col);
                        //singleton_inv2 += (col + 3) * (723732 - (col + 2));
                        continue;
                    }
                    neighbours.inc_nr(v);
                    if (neighbours.get(v) == 0) {
                        color_vertices_considered.inc_nr(col);
                        assert(col + color_vertices_considered.get(col) < g->v_size);
                        scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                        if (!scratch_set.get(col)) {
                            old_color_classes.push_back(col);
                            acc_in += MASH3(col);
                            scratch_set.set(col);
                        }
                    }
                }
                cc += 1;
            }

            // write invariant for singleton color classes
            comp = I->write_top_and_compare(singleton_inv1) && comp;
            comp = I->write_top_and_compare(singleton_inv2) && comp;
            comp = I->write_top_and_compare(-acc_in) && comp;

            // early out before sorting color classes
            /*if(!comp) {
                while(!old_color_classes.empty()) {
                    const int _col = old_color_classes.pop_back();
                    for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i)
                        neighbours.set(scratch[_col + i], -1);
                    color_vertices_considered.set(_col, -1);
                }
                return comp;
            }*/

            // sort split color classes
            old_color_classes.sort();

            // split color classes according to neighbour count
            while (!old_color_classes.empty()) {
                const int _col = old_color_classes.pop_back();
                const int _col_sz = ptn[_col] + 1;
                neighbour_sizes.reset();
                vertex_worklist.reset();

                for (i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                    const int v = scratch[_col + i];
                    const int index = neighbours.get(v) + 1;
                    if (neighbour_sizes.inc(index) == 0)
                        vertex_worklist.push_back(index);
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
                        const int v_degree = k;
                        comp = I->write_top_and_compare(_ncol + v_degree * g->v_size) && comp;
                        comp = I->write_top_and_compare(g->v_size * 7 + val + 1) && comp;
                        if (_ncol != _col)
                            ptn[_ncol] = -1;
                    }
                }

                const int vcount = color_vertices_considered.get(_col);

                // early out
                /*if(!comp) {
                    bool hard_reset = false;
                    if(2 * vcount > g->v_size) {
                        neighbours.reset_hard();
                        hard_reset = true;
                    } else {
                        j = 0;
                        while (j < vcount + 1) {
                            const int v = scratch[_col + j];
                            neighbours.set(v, -1);
                            ++j;
                        }
                    }
                    color_vertices_considered.set(_col, -1);
                    while (!old_color_classes.empty()) {
                        const int __col = old_color_classes.pop_back();
                        if(!hard_reset) {
                            for (i = 0; i < color_vertices_considered.get(__col) + 1; ++i)
                                neighbours.set(scratch[__col + i], -1);
                        }
                        color_vertices_considered.set(__col, -1);
                    }
                    neighbour_sizes.reset();
                    vertex_worklist.reset();
                    color_vertices_considered.reset();
                    return comp;
                }*/

                vertex_worklist.reset();
                j = 0;
                color_vertices_considered.set(_col, -1);

                // determine colors and rearrange vertices
                while (j < vcount + 1) {
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
                for (i = _col; i < _col + _col_sz;) {
                    //assert(i >= 0 && i < ptn_sz);
                    assert(ptn[i] + 1 > 0);
                    mark_as_largest = largest_color_class_size < (ptn[i] + 1);
                    largest_color_class_size = mark_as_largest ? (ptn[i] + 1) : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(_col, i), mark_as_largest));

                    i += ptn[i] + 1;
                }
            }

            neighbour_sizes.reset();
            vertex_worklist.reset();

            return comp;
        }

        bool refine_color_class_dense(sgraph *g, coloring *c,
                                      int color_class, int class_size,
                                      work_list_pair_bool *color_class_split_worklist, invariant *I) {
            bool comp;
            int i, cc, acc, largest_color_class_size, singleton_inv, pos;
            cc = color_class; // iterate over color class
            comp = true;

            neighbours.reset_hard();
            scratch_set.reset();
            old_color_classes.reset();
            vertex_worklist.reset();

            singleton_inv = 0;
            const int end_cc = color_class + class_size;

            // for all vertices of the color class...
            while (cc < end_cc) { // increment value of neighbours of vc by 1
                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int end_i = pe + g->d[vc];
                for (i = pe; i < end_i; i++) {
                    const int v = g->e[i];
                    const int col = c->vertex_to_col[v];
                    if (c->ptn[col] > 0) {
                        neighbours.inc(v); // want to use reset variant?
                        if (!scratch_set.get(col)) {
                            scratch_set.set(col);
                            old_color_classes.push_back(col);
                        }
                    } else {
                        if (config->CONFIG_IR_FULL_INVARIANT)
                            vertex_worklist.push_back(col);
                        else {
                            singleton_inv += MASH4(col);
                        }
                    }
                }
                cc += 1;
            }

            // write singletons
            comp = I->write_top_and_compare(g->v_size * 3 + old_color_classes.cur_pos) && comp;

            if (config->CONFIG_IR_FULL_INVARIANT) {
                vertex_worklist.sort();
                while (!vertex_worklist.empty()) {
                    const int col = vertex_worklist.pop_back();
                    comp = I->write_top_and_compare(g->v_size * 9 + col) && comp;
                }
            } else {
                comp = I->write_top_and_compare(singleton_inv) && comp;
            }

            old_color_classes.sort();

            // for every cell to be split...
            while (!old_color_classes.empty()) {
                const int col = old_color_classes.pop_back();
                const int col_sz = c->ptn[col] + 1;

                if (col_sz == 1)
                    continue;

                neighbour_sizes.reset();
                vertex_worklist.reset();

                const int first_index = neighbours.get(c->lab[col]) + 1;
                vertex_worklist.push_back(first_index);
                int total = 0;
                for (i = col; i < col + col_sz; ++i) {
                    const int v = c->lab[i];
                    int index = neighbours.get(v) + 1;
                    if (index == first_index)
                        continue;
                    total += 1;
                    if (neighbour_sizes.inc(index) == 0)
                        vertex_worklist.push_back(index);
                }

                neighbour_sizes.inc(first_index);
                neighbour_sizes.set(first_index, col_sz - total - 1);

                comp = I->write_top_and_compare(vertex_worklist.cur_pos) && comp;

                if (vertex_worklist.cur_pos == 1) {
                    // no split
                    const int v_degree = neighbours.get(c->lab[col]);
                    comp = I->write_top_and_compare(col + g->v_size * v_degree) && comp;
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
                        const int _col = col + col_sz - (neighbour_sizes.get(k));
                        c->ptn[_col] = -1; // this is val - 1, actually...

                        const int v_degree = k;
                        comp = I->write_top_and_compare(_col + g->v_size * v_degree) && comp;
                        comp = I->write_top_and_compare(_col + val + 1) && comp;
                    }
                }

                // copy cell for rearranging
                memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                //vertex_worklist.cur_pos = col_sz;
                pos = col_sz;

                // determine colors and rearrange
                // split color classes according to count in counting array
                while (pos > 0) {
                    const int v = scratch[--pos];
                    const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                    // i could immediately determine size of v_new_color here and write invariant before rearrange?
                    assert(v_new_color >= col);
                    assert(v_new_color < col + col_sz);

                    c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                    c->vertex_to_col[v] = v_new_color;
                    c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                    c->ptn[v_new_color] += 1;
                }

                largest_color_class_size = -1;

                // determine largest class to throw away and finish (fourth iteration)
                for (i = col; i < col + col_sz;) {
                    assert(i >= 0 && i < c->ptn_sz);
                    assert(c->ptn[i] + 1 > 0);
                    const int v_color = i;
                    const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                    largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(col, i), mark_as_largest));
                    if (v_color != 0)
                        c->ptn[v_color - 1] = 0;
                    i += c->ptn[i] + 1;
                }
            }

            return comp;
        }

        bool refine_color_class_dense_dense(sgraph *g, coloring *c,
                                            int color_class, int class_size,
                                            work_list_pair_bool *color_class_split_worklist, invariant *I) {
            bool comp;
            int i, j, acc, cc, largest_color_class_size;
            cc = color_class; // iterate over color class
            comp = true;

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

            // in dense-dense we dont sort, instead we just iterate over all cells (which are already sorted)
            old_color_classes.sort();
            // for every cell...
            for (j = 0; j < g->v_size;) {
                const int col = j;
                j += c->ptn[j] + 1;
                const int col_sz = c->ptn[col] + 1;
                if (col_sz == 1) {
                    const int v_degree = neighbours.get(c->lab[col]) + 1;
                    if (v_degree == -1) // not connected! skip...
                        continue;
                    comp = I->write_top_and_compare(col + v_degree * g->v_size) && comp;
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

                comp = I->write_top_and_compare(g->v_size * 12 + vertex_worklist.cur_pos) && comp;
                // if(!comp) {neighbours.reset_hard(); return comp;}

                if (vertex_worklist.cur_pos == 1) {
                    // no split
                    const int v_degree = neighbours.get(c->lab[col]);
                    //comp = comp && I->write_top_and_compare(-g->v_size * 10 - col);
                    comp = I->write_top_and_compare(col + v_degree * g->v_size) && comp;
                    comp = I->write_top_and_compare(col + c->ptn[col] + 1) && comp;
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
                        const int _col = col + col_sz - (neighbour_sizes.get(k));
                        c->ptn[_col] = -1; // this is val - 1, actually...

                        const int v_degree = k;
                        comp = I->write_top_and_compare(-g->v_size * 5 - _col) && comp;
                        comp = I->write_top_and_compare(_col + v_degree) && comp;
                        comp = I->write_top_and_compare(_col + val + 1) && comp;
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
                    const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                    // i could immediately determine size of v_new_color here and write invariant before rearrange?
                    assert(v_new_color >= col);
                    assert(v_new_color < col + col_sz);

                    c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                    c->vertex_to_col[v] = v_new_color;
                    c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                    c->ptn[v_new_color] += 1;
                }

                largest_color_class_size = -1;

                // determine largest class to throw away and finish (fourth iteration)
                for (i = col; i < col + col_sz;) {
                    assert(i >= 0 && i < c->ptn_sz);
                    assert(c->ptn[i] + 1 > 0);
                    const int v_color = i;
                    const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                    largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(col, i), mark_as_largest));
                    if (v_color != 0)
                        c->ptn[v_color - 1] = 0;
                    i += c->ptn[i] + 1;
                }
            }

            neighbours.reset_hard();
            return comp;
        }

        bool refine_color_class_singleton(sgraph *g, coloring *c,
                                          int color_class, int class_size,
                                          work_list_pair_bool *color_class_split_worklist, invariant *I) {
            bool comp;
            int i, cc, deg1_write_pos, deg1_read_pos, singleton_inv;
            cc = color_class; // iterate over color class
            comp = true;

            neighbours.reset();
            scratch_set.reset();
            vertex_worklist.reset();
            old_color_classes.reset();

            singleton_inv = 0;
            const int vc = c->lab[cc];
            const int pe = g->v[vc];
            const int end_i = pe + g->d[vc];
            for (i = pe; i < end_i; i++) {
                const int v = g->e[i];
                const int col = c->vertex_to_col[v];

                if (c->ptn[col] == 0) {
                    // if full invariant, sort -- else use a hash value
                    if (config->CONFIG_IR_FULL_INVARIANT)
                        vertex_worklist.push_back(col); // treat singletons in separate list (more efficient sorting)
                    else
                        singleton_inv += MASH5(col);
                    continue;
                }

                if (!scratch_set.get(col)) {
                    scratch_set.set(col);
                    old_color_classes.push_back(col);
                    // neighbours acts as the write position for degree 1 vertices of col in the scratchpad
                    neighbours.set(col, col);
                }
                // write vertex to scratchpad for later use
                scratch[neighbours.get(col)] = v;
                neighbours.inc_nr(col); // we reset neighbours later, since old_color_classes for reset
            }

            singleton_inv += old_color_classes.cur_pos;
            if (!config->CONFIG_IR_FULL_INVARIANT)
                comp = I->write_top_and_compare(g->v_size * 3 + singleton_inv) && comp;

            old_color_classes.sort();

            // write invariant first...
            for (i = 0; i < old_color_classes.cur_pos && comp; ++i) {
                //comp = comp && I->write_top_and_compare(old_color_classes.arr[i]); // color class
                comp = I->write_top_and_compare(g->v_size * 14 + neighbours.get(old_color_classes[i])) && comp;
                // contains information about color degree (= 1)
            }

            // sort and write down singletons in invariant
            if (config->CONFIG_IR_FULL_INVARIANT)
                vertex_worklist.sort();

            for (i = 0; i < vertex_worklist.cur_pos; ++i) {
                comp = I->write_top_and_compare(g->v_size * 11 + vertex_worklist[i]) && comp; // size
                // should contain information about color degree
            }

            while (!old_color_classes.empty()) {
                const int deg0_col = old_color_classes.pop_back();
                const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                const int deg1_col = deg0_col + deg0_col_sz;

                assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);

                // no split? done...
                if (deg0_col == deg1_col) {
                    neighbours.set(deg1_col, -1);
                    assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                    continue;
                }

                assert(deg0_col_sz + deg1_col_sz - 1 == c->ptn[deg0_col]);

                // set ptn
                c->ptn[deg0_col] = deg0_col_sz - 1;
                c->ptn[deg1_col] = deg1_col_sz - 1;
                c->ptn[deg1_col - 1] = 0;

                deg1_write_pos = deg1_col;
                deg1_read_pos = neighbours.get(deg0_col) - 1;

                //c->vertex_to_col[c->lab[deg1_col]] = deg1_col;

                // rearrange vertices of deg1 to the back of deg0 color
                assert(deg1_read_pos >= deg0_col);

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

                assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                assert(c->vertex_to_col[c->lab[deg1_col]] == deg1_col);

                // add new classes to color_class_split_worklist
                const bool leq = deg1_col_sz > deg0_col_sz;
                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(deg0_col, deg0_col), !leq));
                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(deg0_col, deg1_col), leq));

                // reset neighbours count to -1
                neighbours.set(deg0_col, -1);
            }

            return comp;
        }

        bool refine_color_class_singleton_first(sgraph *g, coloring *c,
                                                const int color_class, work_list_pair_bool *color_class_split_worklist) {

            const int vc = c->lab[color_class];
            const int pe = g->v[vc];
            const int deg = g->d[vc];

            int deg1_write_pos, deg1_read_pos;
            int* it;
            neighbours.reset();
            scratch_set.reset();
            old_color_classes.reset();
            int  first_col = -1; // buffer strategy

            // arbitrary degree
            const int* end_pt = g->e + pe + deg;
            for (it = g->e + pe; it < end_pt; it++) {
                const int v = *it;
                const int col    = c->vertex_to_col[v];
                const int col_sz = c->ptn[col];
                if (col_sz == 0) {
                    continue;
                }

                if (first_col == -1 || (col != first_col && !scratch_set.get(col))) {
                    if(first_col == -1) {
                        first_col = col;
                    } else {
                        scratch_set.set(col);
                    }
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

            while (!old_color_classes.empty()) {
                const int deg0_col    = old_color_classes.pop_back();
                const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                const int deg1_col    = deg0_col + deg0_col_sz;

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

                // add new classes to color_class_split_worklist
                const bool leq = deg1_col_sz > deg0_col_sz;
                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(deg0_col, deg0_col), !leq));
                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(deg0_col, deg1_col), leq));

                // reset neighbours count to -1
                neighbours.set(deg0_col, -1);
            }

            return true;
        }

        bool refine_color_class_dense_first(sgraph *g, coloring *c,
                                            int color_class, int class_size,
                                            work_list_pair_bool *color_class_split_worklist) {
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
                memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                pos = col_sz;

                // determine colors and rearrange
                // split color classes according to count in counting array
                while (pos > 0) {
                    const int v = scratch[--pos];
                    const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                    // i could immediately determine size of v_new_color here and write invariant before rearrange?
                    assert(v_new_color >= col);
                    assert(v_new_color < col + col_sz);

                    c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                    c->vertex_to_col[v] = v_new_color;
                    c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                    c->ptn[v_new_color] += 1;
                }

                largest_color_class_size = -1;

                // determine largest class to throw away and finish (fourth iteration)
                for (i = col; i < col + col_sz;) {
                    assert(i >= 0 && i < c->ptn_sz);
                    assert(c->ptn[i] + 1 > 0);
                    const int v_color = i;
                    const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                    largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(col, i), mark_as_largest));
                    if (v_color != 0)
                        c->ptn[v_color - 1] = 0;
                    i += c->ptn[i] + 1;
                }
            }

            return true;
        }

        bool refine_color_class_dense_dense_first(sgraph *g, coloring *c,
                                                  int color_class, int class_size,
                                                  work_list_pair_bool *color_class_split_worklist) {
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
                memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                pos = col_sz;

                // determine colors and rearrange
                // split color classes according to count in counting array
                while (pos > 0) {
                    const int v = scratch[--pos];
                    const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                    assert(v_new_color >= col);
                    assert(v_new_color < col + col_sz);

                    c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                    c->vertex_to_col[v] = v_new_color;
                    c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                    c->ptn[v_new_color] += 1;
                }

                largest_color_class_size = -1;

                // determine largest class to throw away and finish (fourth iteration)
                for (i = col; i < col + col_sz;) {
                    assert(i >= 0 && i < c->ptn_sz);
                    assert(c->ptn[i] + 1 > 0);
                    const int v_color = i;
                    const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                    largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(col, i), mark_as_largest));
                    if (v_color != 0)
                        c->ptn[v_color - 1] = 0;
                    i += c->ptn[i] + 1;
                }
            }

            neighbours.reset_hard();
            return true;
        }

        bool refine_color_class_sparse_first(sgraph *g, coloring *c,
                                             int color_class, int class_size,
                                             work_list_pair_bool *color_class_split_worklist) {
            bool comp;
            int v_new_color, cc, largest_color_class_size, acc;
            int* it;

            cc = color_class; // iterate over color class
            comp = true;

            scratch_set.reset();
            old_color_classes.reset();
            neighbours.reset();
            color_vertices_considered.reset();

            const int end_cc = color_class + class_size;
            while (cc < end_cc) { // increment value of neighbours of vc by 1
                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int* end_pt = g->e + pe + g->d[vc];
                for (it = g->e + pe; it < end_pt; it++) {
                    const int v = *it;
                    const int col = c->vertex_to_col[v];
                    if (c->ptn[col] == 0) {
                        continue;
                    }

                    neighbours.inc_nr(v);
                    if (neighbours.get(v) == 0) {
                        color_vertices_considered.inc_nr(col);
                        assert(col + color_vertices_considered.get(col) < g->v_size);
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

                int total = 0;
                for (int i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                    const int v = scratch[_col + i];
                    int index = neighbours.get(v) + 1;
                    total += 1;
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
                        assert(false);
                        v_new_color = _col;
                    } else {
                        assert((neighbours.get(v) + 1) > 0);
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
                        assert(v_new_color2 > _col);
                        c->ptn[_col] -= 1;
                    } else
                        assert(false);
                }

                // add new colors to worklist
                largest_color_class_size = -1;
                int i;
                for (i = _col; i < _col + _col_sz;) {
                    assert(i >= 0 && i < c->ptn_sz);
                    assert(c->ptn[i] + 1 > 0);
                    const bool mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
                    largest_color_class_size = mark_as_largest ? (c->ptn[i] + 1) : largest_color_class_size;

                    color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(_col, i), mark_as_largest));
                    if (i != 0)
                        c->ptn[i - 1] = 0;
                    i += c->ptn[i] + 1;
                }

                if (i != 0)
                    c->ptn[i - 1] = 0;
            }

            neighbour_sizes.reset();
            vertex_worklist.reset();

            return comp;
        }
    };

// ring queue for pairs of integers
    class ring_pair {
        std::pair<int, int> *arr = 0;
        bool init = false;
        int arr_sz = -1;
        int front_pos = -1;
        int back_pos = -1;

    public:
        void initialize(int size) {
            if(init)
                delete[] arr;
            arr = new std::pair<int, int>[size];
            arr_sz = size;
            back_pos = 0;
            front_pos = 0;
            init = true;
        }

        void push_back(std::pair<int, int> value) {
            arr[back_pos] = value;
            back_pos = (back_pos + 1) % arr_sz;
        }

        std::pair<int, int> *front() {
            return &arr[front_pos];
        }

        void pop() {
            front_pos = (front_pos + 1) % arr_sz;
        }

        bool empty() {
            return (front_pos == back_pos);
        }

        ~ring_pair() {
            if (init)
                delete[] arr;
        }

        void reset() {
            front_pos = back_pos;
        }
    };
}

#endif //SASSY_REFINEMENT_H
