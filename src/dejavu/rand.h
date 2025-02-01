// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_RAND_H
#define DEJAVU_RAND_H

#include "ir.h"
#include "groups.h"

namespace dejavu {namespace search_strategy {
    // (need to nest namespaces due to C++ 14)

    /**
     * \brief IR search using random walks.
     *
     * Performs random walks of the IR shared_tree, sifting resulting automorphisms into the given Schreier structure. If the
     * Schreier structure is complete with respect to the base, or the probabilistic abort criterion satisfied, the
     * process terminates. The algorithm guarantees to find all automorphisms up to the specified error bound.
     *
     * Alternatively, a limit for the amount of discovered differing leaves can be set.
     *
     * Objects contain a local workspace which contain datastructures only used for random_ir, as well as settings and
     * for statistics of the search.
     */
    class random_ir {
        random_source& rng;
        std::vector<int> heuristic_reroll;

        timed_print& gl_printer;
        groups::schreier_workspace&     gl_schreierw;
        groups::automorphism_workspace& gl_automorphism;

        /**
         * Loads the leaf into the \p local_state. Only works if base of leaf was stored, i.e., whenever
         * `full_save = false` was used in the constructor.
         *
         * @param leaf The leaf.
         * @param g The graph.
         * @param local_state Local state in which the leaf will be stored.
         * @param start_from State from which the walk to the leaf is performed.
         */
        static void load_state_from_leaf(sgraph *g, ir::controller &local_state, ir::limited_save &start_from,
                                  ir::stored_leaf *leaf) {
            dej_assert(leaf->get_store_type() == ir::stored_leaf::STORE_BASE);
            std::vector<int> base;
            base.reserve(leaf->get_lab_or_base_size());
            for(int i = 0; i < leaf->get_lab_or_base_size(); ++i) base.push_back(leaf->get_lab_or_base()[i]);
            local_state.walk(g, start_from, base);
        }

        /**
         * Co-routine which adds a leaf to leaf_storage, and sifts resulting automorphism into a given group.
         *
         * @returns whether an automorphism was found and it sifted successfully
         */
        bool add_leaf_to_storage_and_group(sgraph *g, dejavu_hook *hook, groups::compressed_schreier &group,
                                           ir::shared_leaves &leaf_storage, ir::controller &local_state,
                                           ir::controller &other_state, ir::limited_save &root_save, bool uniform) {
            dej_assert(g->v_size == local_state.c->cells);

            // outer loop to handle hash collisions -- at most h_hash_col_limit will be checked and stored
            for(int hash_offset = 0; hash_offset < h_hash_col_limit; ++hash_offset) {
                // after first hash collision, we write a stronger invariant
                if(hash_offset == 1) local_state.write_strong_invariant_quarter(g);

                // after second hash collision, an even stronger one
                if(hash_offset == 2) local_state.write_strong_invariant(g);

                // first, test whether leaf with same hash has already been stored
                const unsigned long hash_c = local_state.T->get_hash() + hash_offset; // '+hash_offset' is for hash
                                                                                      // collisions
                auto other_leaf = leaf_storage.lookup_leaf(hash_c);

                // if not, add leaf to leaf_storage
                if (other_leaf == nullptr) {
                    ++s_leaves;
                    ++s_paths_failany;
                    s_rolling_success = (9.0 * s_rolling_success + 0.0) / 10.0;
                    leaf_storage.add_leaf(hash_c, *local_state.c, local_state.base_vertex);
                    break;
                }

                // if there is a leaf with the same hash, load the leaf and test automorphism
                gl_automorphism.reset();

                if(other_leaf->get_store_type() == ir::stored_leaf::STORE_LAB) {
                    // sometimes, a lab array is already stored for the leaf
                    const int* lab = other_leaf->get_lab_or_base();
                    gl_automorphism.write_color_diff(local_state.c->vertex_to_col, lab);
                } else {
                    // other times, it is not and we need to recompute the walk to this leaf
                    load_state_from_leaf(g, other_state, root_save, other_leaf);
                    gl_automorphism.write_color_diff(local_state.c->vertex_to_col, other_state.c->lab);
                }

                // certify whether we actually found an automorphism
                const bool cert = local_state.certify(g, gl_automorphism);

                if (cert) {
                    // We found an automorphism!
                    s_rolling_success = (9.0 * s_rolling_success + 1.0) / 10.0;
                    ++s_succeed;

                    // update and check support limits
                    s_support_limit_reached = (h_schreier_support_limit >= 0) &&
                                              (group.get_support() >h_schreier_support_limit);

                    // output the automorphism
                    if(hook) (*hook)(g->v_size, gl_automorphism.p(), gl_automorphism.nsupp(),
                                     gl_automorphism.supp());

                    // sift automorphism into Schreier structure
                    bool sift = group.sift(gl_schreierw, gl_automorphism, uniform);
                    gl_automorphism.reset();

                    // if sifting changed the Schreier structure, consider sifting some more random elements to fill
                    // up the Schreier structure
                    if(sift && group.s_densegen() + group.s_sparsegen() > 1 && s_random_sift_success > -5) {
                        int fail = 3; // 3
                        bool any_changed = false;
                        while(fail >= 0) {
                            const bool sift_changed = group.sift_random(gl_schreierw, gl_automorphism, rng);
                            any_changed = sift_changed || any_changed;
                            fail -= !sift_changed;
                        }

                        s_random_sift_success += any_changed?1:-1;
                        s_random_sift_success  = std::max(std::min(s_random_sift_success, 5), -5);
                    }

                    // reset automorphism workspace to identity, return result of sift
                    gl_automorphism.reset();
                    return sift;
                }

                // actually, this wasn't an automorphism: there was a collision for the value of the invariant, even
                // though leaves are not equivalent
            }
            gl_automorphism.reset();
            return false;
        }

    public:
        // stats
        double    s_rolling_success = 0;                  /**< rolling probability how many random paths succeed     */
        double    s_rolling_first_level_success  = 1.0;   /**< rolling probability how many random paths succeed on the
                                                            *  first level*/

        long      s_trace_cost1   = 0;                    /**< total cost incurred on first level         */
        int       s_paths         = 0;                    /**< how many total paths have been computed    */
        int       s_paths_fail1   = 0;                    /**< how many total paths failed on first level */
        int       s_paths_failany = 0;                    /**< how many total paths have failed           */
        int       s_succeed       = 0;                    /**< how many total paths have succeeded        */
        int       s_leaves        = 0;                    /**< how many leaves were added                 */
        int       s_min_split_number = 0;
        bool      s_support_limit_reached = false;        /**< Schreier support limit was reached         */

        int s_random_sift_success = 0;

        // settings for heuristics
        bool      h_look_close      = false;              /**< whether to use trace early out on first level   */
        const int h_hash_col_limit  = 32;                 /**< limit for how many hash collisions are allowed  */
        bool      h_sift_random     = true;               /**< sift random elements into Schreier structure    */
        int       h_sift_random_lim = 8;                  /**< after how many paths random elements are sifted */
        int       h_randomize_up_to = INT32_MAX;          /**< randomize vertex selection up to this level */
        long      h_schreier_support_limit = -1;          /**< impose a limit on support of Schreier structure */

        void use_look_close(bool look_close = false) {
            h_look_close = look_close;
        }

        /**
         * Links this object to local workspaces.
         *
         * @param schreier Schreier workspace
         * @param automorphism workspace to store automorphisms
         */
        random_ir(timed_print& printer, groups::schreier_workspace& schreier,
                  groups::automorphism_workspace& automorphism, random_source& rgenerator) :
                  rng(rgenerator), gl_printer(printer), gl_schreierw(schreier), gl_automorphism(automorphism)
                  {}

        /**
         * Resets all recorded statistics.
         */
        void reset_statistics() {
            s_paths            = 0;
            s_paths_fail1      = 0;
            s_trace_cost1      = 0;
            s_paths_failany    = 0;
            s_leaves           = 0;
            s_succeed          = 0;
            s_rolling_success  = 0;
            s_rolling_first_level_success  = 1.0;
            s_min_split_number = INT32_MAX;
            s_support_limit_reached    = false;
        }

        static bool h_almost_done(groups::compressed_schreier &group) {
            return group.get_consecutive_success() >= 2;
        }

        static void
        specific_walk(sgraph *g, ir::shared_tree &ir_tree, ir::controller &local_state, std::vector<int> &base_vertex) {
            local_state.walk(g, *ir_tree.pick_node_from_level(0,0)->get_save(), base_vertex);
            auto other_leaf = ir_tree.stored_leaves.lookup_leaf(local_state.T->get_hash());
            if(other_leaf == nullptr) {
                ir_tree.stored_leaves.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
            }
        }

        /**
         * Performs Monte Carlo IR search.
         *
         * Returns either when probabilistic or deterministic abort criterion is satisfied, or whenever the
         * \p fail_limit is reached.
         *
         * May use non-uniform base-alignment to fill Schreier tables faster.
         *
         * @param g graph
         * @param hook hook to return automorphisms
         * @param selector cell selector
         * @param ir_tree ir tree computed so far
         * @param group Schreier structure to sift found automorphisms into, used to calculate probabilistic abort
         *              criterion
         * @param local_state Local workspace used to perform random walks of IR tree
         * @param fail_limit Limits the number of random IR walks not leading to a new automorphisms
         */
        void random_walks(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                          ir::shared_tree &ir_tree, groups::compressed_schreier &group, ir::controller &local_state,
                          ir::controller& other_state, int fail_limit) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            other_state.use_reversible(false);
            other_state.use_trace_early_out(false);

            // start from root of tree initially
            ir::limited_save my_own_save;
            ir::limited_save* root_save  = ir_tree.pick_node_from_level(0, 0)->get_save();
            ir::limited_save* start_from = root_save;

            // other_state is used to retrieve leaves that are stored as walks of the tree
            other_state.link_compare(&local_state);

            // we want to record how expensive walks are
            const int s_cell_initial   = start_from->get_coloring()->cells;
            double    progress_initial = 1.0 * s_cell_initial / g->v_size;

            const int target_level = ir_tree.get_finished_up_to();

            int  s_sifting_success = 0;

            // continue until budget exhausted, search successful, or support limits reached
            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion() &&
                    s_paths_failany < fail_limit && !s_support_limit_reached) {
                local_state.load_reduced_state(*start_from);

                int could_start_from = group.finished_up_to_level();

                // print progress, sometimes
                if(s_paths_failany > 8 && (s_paths & 0x00000FFF) == 0x000000FE)
                    gl_printer.progress_current_method("random", "leaves", ir_tree.stat_leaves(), "f1", s_paths_fail1,
                                                       "compress", group.s_compression_ratio);

                // can start from below the root if we finished Schreier table at the current root
                if(local_state.s_base_pos <= could_start_from) {
                    while (local_state.s_base_pos <= could_start_from) {
                        local_state.move_to_child(g, group.base_point(local_state.s_base_pos));
                    }
                    dej_assert(local_state.T->trace_equal());
                    local_state.save_reduced_state(my_own_save); // from now on, we start from this save!
                    start_from = &my_own_save;

                    const int    s_cells_now  = start_from->get_coloring()->cells;
                    const double progress_now = 1.0 * s_cells_now / g->v_size;

                    if(progress_now - progress_initial > 0.1) {
                        gl_printer.progress_current_method("random", "root_cells", 1.0 * s_cells_now / g->v_size,
                                                           "base_pos", could_start_from,
                                                           "sift", s_sifting_success, "rsift", s_random_sift_success);
                        progress_initial = progress_now;
                    }
                }

                const int start_from_base_pos = local_state.s_base_pos;
                int base_pos                  = local_state.s_base_pos;

                // track whether current walk is base-aware and/or uniform
                //bool base_aligned = true;
                bool uniform      = true;

                //walk down the tree as long as we are not in a leaf
                while (g->v_size != local_state.c->cells) {
                    const int col = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = rng() % col_sz;
                    int choose_pos = col + rand;

                    // if we are beyond where we need to sift, and we don't want to sample uniformly at random, we
                    // stop picking random elements whatsoever
                    /*if(base_pos >= h_randomize_up_to && !h_sift_random && ir_tree.stored_leaves.s_leaves <= 1 &&
                       base_pos < static_cast<int>(local_state.compare_base_vertex->size())) {
                        uniform    = false; // sampled leaf not uniform anymore now
                        choose_pos = col;   // just pick first vertex of color

                        // or even better: let's choose the base vertex, if it's in the correct color
                        const int v_base = (*local_state.compare_base_vertex)[base_pos];
                        const int v_base_col = local_state.c->vertex_to_col[v_base];
                        if(col == v_base_col) choose_pos = local_state.c->vertex_to_lab[v_base];
                    }*/

                    // in the case where we are not succeeding in sifting randomly generated elements in the group
                    // itself, let's try to create sparse generators by trying to stick to the base after a few
                    // individualizations
                    if(base_pos > start_from_base_pos + 1 && g->v_size > 5000 && s_sifting_success >= 0 &&
                       s_random_sift_success < 0 &&
                       base_pos < static_cast<int>(local_state.compare_base_vertex->size())) {
                        // or even better: let's choose the base vertex, if it's in the correct color
                        const int v_base     = (*local_state.compare_base_vertex)[base_pos];
                        const int v_base_col = local_state.c->vertex_to_col[v_base];

                        if(col == v_base_col) {
                            uniform = false;
                            choose_pos = local_state.c->vertex_to_lab[v_base];
                        }
                    }

                    int v = local_state.c->lab[choose_pos];

                    // base-aware search: if we are still walking along the base, and the vertex we picked is in the
                    // same orbit as the base -- we might as well keep walking on the base, or choose a different vertex
                    if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)
                         && ir_tree.stored_leaves.s_leaves <= 1) {
                        heuristic_reroll.clear();
                        for(int i = 0; i < col_sz; ++i) {
                            heuristic_reroll.push_back(local_state.c->lab[col + i]);
                        }
                        group.reduce_to_unfinished(gl_schreierw, heuristic_reroll, base_pos);
                        if(!heuristic_reroll.empty()) {
                            const int rand_reroll = rng() % static_cast<int>(heuristic_reroll.size());
                            v = heuristic_reroll[rand_reroll];
                        }
                    }

                    dej_assert(local_state.c->vertex_to_col[v] == col);
                    const int trace_pos_pre = local_state.T->get_position();
                    local_state.use_trace_early_out((base_pos == target_level) && !h_look_close);
                    local_state.move_to_child(g, v);

                    // keep track of some statistics for the first individualization (these statistics are used for
                    // decisions concerning breadth-first search)
                    if(base_pos == target_level) {
                        s_min_split_number = std::min(local_state.get_number_of_splits(), s_min_split_number);
                        s_trace_cost1 += local_state.T->get_position() - trace_pos_pre;
                        s_paths_fail1 += !local_state.T->trace_equal();
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        if(!h_look_close && !local_state.T->trace_equal()) break;
                    }

                    ++base_pos;
                    dej_assert(base_pos == local_state.s_base_pos);
                }

                ++s_paths;
                if(base_pos == target_level) { // did not arrive in a leaf
                    ++s_paths_failany;
                    continue;
                }

                // we arrived at a leaf... let's check whether we already have an equivalent leaf to form an
                // automorphism -- and otherwise we just add this leaf to the storage
                const bool sift = add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state,
                                                                other_state, *root_save, uniform);
                s_sifting_success += sift?1:-1;
                s_sifting_success += sift && !uniform?1:0;
                s_sifting_success = std::max(std::min(s_sifting_success, 10), -10);
            }
        }

        /**
         * Performs Monte Carlo IR search. Uses a given IR tree to initiate the search further down the tree (i.e., if
         * BFS has been performed, we can start from the furthest BFS level).
         *
         * Returns either when probabilistic or deterministic abort criterion is satisfied, or whenever the
         * \p fail_limit is reached.
         *
         * @param g graph
         * @param hook hook to return automorphisms
         * @param selector cell selector
         * @param ir_tree ir tree computed so far
         * @param group Schreier structure to sift found automorphisms into, used to calculate probabilistic abort
         *              criterion
         * @param local_state Local workspace used to perform random walks of IR tree
         * @param fail_limit Limits the number of random IR walks not leading to a new automorphisms
         */
        void random_walks_from_tree(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                                    ir::shared_tree &ir_tree, groups::compressed_schreier &group,
                                    ir::controller &local_state, ir::controller& other_state, int fail_limit) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            s_rolling_first_level_success = 1;
            const int pick_from_level = ir_tree.get_finished_up_to();

            // other_state is used to retrieve leaves that are stored as walks of the tree
            other_state.link_compare(&local_state);

            // continue until budget exhausted, or search successful
            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion() &&
                    s_paths_failany < fail_limit && !s_support_limit_reached) {

                // print progress, sometimes
                if((s_paths & 0x000000FF) == 0x000000FE)
                    gl_printer.progress_current_method("random", "leaves", ir_tree.stat_leaves(), "f1", s_paths_fail1,
                                                       "compress", group.s_compression_ratio);

                // start walk from random node at current breadth-first level
                auto node = ir_tree.pick_node_from_level(pick_from_level, rng());
                local_state.load_reduced_state(*node->get_save());

                int base_pos                  = local_state.s_base_pos;
                const int start_from_base_pos = base_pos;

                // now perform the random walk
                while (g->v_size != local_state.c->cells) {
                    const int col    = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = rng() % col_sz;
                    int v = local_state.c->lab[col + rand];
                    const int trace_pos_pre = local_state.T->get_position();
                    local_state.use_trace_early_out((base_pos == start_from_base_pos) && !h_look_close);
                    local_state.move_to_child(g, v);

                    if(base_pos == start_from_base_pos) {
                        s_min_split_number = std::min(local_state.get_number_of_splits(), s_min_split_number);
                        s_trace_cost1 += local_state.T->get_position() - trace_pos_pre;
                        s_paths_fail1 += !local_state.T->trace_equal();
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        if(!h_look_close && !local_state.T->trace_equal()) break;
                    }
                    ++base_pos;
                }

                ++s_paths;
                if(base_pos == start_from_base_pos) {
                    ++s_paths_failany;
                    continue;
                }

                // we arrived at a leaf... let's check whether we already have an equivalent leaf to form an
                // automorphism -- and otherwise we just add this leaf to the storage
                add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state, other_state,
                                              *ir_tree.pick_node_from_level(0, 0)->get_save(), true);
            }
        }
    };
}}

#endif //DEJAVU_RAND_H
