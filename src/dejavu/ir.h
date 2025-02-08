// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_IR_H
#define DEJAVU_IR_H

#include <cstdint>
#include <unordered_map>
#include "refinement.h"
#include "coloring.h"
#include "graph.h"
#include "trace.h"
#include "groups.h"

namespace dejavu {

    /**
     * \brief IR fundamentals.
     *
     * Contains fundamental algorithms and data structures to implement individualization-refinement algorithms. This
     * includes graphs, colorings, color refinement, cell selectors as well as higher level control mechanisms.
     */
    namespace ir {
        /**
         * \brief Mode of trace for IR search
         *
         * The `ir_mode` determines in which mode the trace is used: whether a new trace is recorded, or whether the
         * current computation is compared to a stored trace.
         *
         */
        enum ir_mode {
            IR_MODE_RECORD_TRACE, IR_MODE_COMPARE_TRACE_REVERSIBLE, IR_MODE_COMPARE_TRACE_IRREVERSIBLE
        };

        /**
         * int type_selector_hook(coloring* c, const int s_base_pos);
         */
        typedef int type_selector_hook(const coloring *, const int);

        /**
         * \brief Reduced IR save state
         *
         * Using this class, partial information of a state of an IR computation can be stored. Using this information,
         * IR computations can be resumed from this state either using BFS or random walks. The state in particular does not
         * keep enough information to resume using DFS.
         */
        class limited_save {
            // TODO enable a "compressed state" only consisting of base (but code outside should be oblivious to it)
            // TODO this is only supposed to be an "incomplete" state -- should there be complete states?

            std::vector<int> base_vertex; /**< base of vertices of this IR node  */
            coloring c;                   /**< vertex coloring of this IR node   */
            uint64_t invariant = 0;       /**< hash of invariant of this IR node */
            int trace_position = 0;       /**< position of trace of this IR node */
            int base_position  = 0;       /**< length of base of this IR node    */
        public:
            void save(std::vector<int>& s_base_vertex, coloring &s_c, uint64_t s_invariant, int s_trace_position,
                      int s_base_position) {
                this->base_vertex = s_base_vertex;
                this->c.copy_any(&s_c);
                this->invariant = s_invariant;
                this->trace_position = s_trace_position;
                this->base_position = s_base_position;
            }

            /**
             * @return coloring of this IR node
             */
            coloring *get_coloring() {
                return &c;
            }

            /**
             * @return hash of invariant of this IR node
             */
            dej_nodiscard uint64_t get_invariant_hash() const {
                return invariant;
            }

            /**
             * @return position of trace of this IR node
             */
            dej_nodiscard int get_trace_position() const {
                return trace_position;
            }

            /**
             * @return length of base of this IR node
             */
            dej_nodiscard int get_base_position() const {
                return base_position;
            }

            /**
             * @return base of this IR node
             */
            std::vector<int>& get_base() {
                return base_vertex;
            }
        };

        /**
         * \brief Tracks information for a base point
         */
        struct base_info {
            int color; /**< color of the base point */
            int color_sz; /**< color size of the base point */
            int cells; /**< number of cells of the coloring*/
            int touched_color_list_pt; /**< position of the touched color list */
            int singleton_pt; /**< position of the singleton list */

            int  trace_pos; /**< position of the trace */
            uint64_t trace_hash; /**< hash of the trace */

            base_info(int selColor, int colorSz, int cellNum, int touchedColorListPt, int singletonPt, int tracePos,
                      uint64_t traceHash) :
                      color(selColor), color_sz(colorSz), cells(cellNum), touched_color_list_pt(touchedColorListPt),
                      singleton_pt(singletonPt), trace_pos(tracePos), trace_hash(traceHash) {}
        };

        /**
         * \brief Controls movement in IR tree
         *
         * Keeps a state of an IR node. Enables the movement to a child of the node, or to the parent of the node.
         * The controller manages data structures and functions, which facilitate the trace \a T as well as the reversal of
         * color refinement.
         *
         * Has different modes (managed by \a mode) depending on whether color refinement should be reversible or not.
         */
        struct controller {
        public:
            coloring      *c  = nullptr; /**< coloring to be moved around in the tree */
            trace         *T  = nullptr; /**< trace of current root-to-node path */

            std::vector<int>       singletons;         /** singletons of the current trace */
            std::vector<int>*      compare_singletons = nullptr; /** singletons of the comparison trace */

            std::vector<int>       base_vertex; /** current base (i.e., root-to-node path) */
            std::vector<base_info> base;        /** additional info for the current base   */

            std::vector<int>*      compare_base_vertex = nullptr; /** comparison base */
            std::vector<base_info>*compare_base = nullptr;        /** additional info of the comparison base */

            coloring leaf_color; /** comparison leaf coloring */

            markset  touched_color; /**< were changes in this color already tracked? */
            worklist touched_color_list; /**< color touched_color_list[i] was changed... */
            worklist prev_color_list; /**< ...and its vertices were previously of color prev_color_list[i]  */

            // statistics
            int   s_base_pos = 0; /**< how large the base of the current IR node is*/
        private:
            refinement* R; /**< refinement workspace to use */
            trace *cT = nullptr; /**< comparison trace */
            trace internal_T1; /**< trace 1 -- trace 1 and trace 2 can be flipped to switch from reading to writing */
            trace internal_T2; /**< trace 2 -- trace 1 and trace 2 can be flipped to switch from reading to writing */

            std::vector<int>       internal_compare_singletons;  /** singletons of the comparison trace */
            std::vector<int>       internal_compare_base_vertex; /** comparison base */
            std::vector<base_info> internal_compare_base;        /** additional info of the comparison base */

            // the following workspaces are only used for the paired color dfs (AKA the saucy-style dfs) -- not needed
            // for normal bfs, dfs, or random search operation
            markset  diff_tester;        /**< used for algorithms for 'difference testing' */
            markset  diff_vertices;      /**< vertices that differ */
            markset  diff_is_singleton;  /**< vertices are singleton */
            worklist  diff_vertices_list; /**< vertices that differ, but in a list */
            worklist  diff_vertices_list_pt; /**< vertex v is stored at position diff_vertices_list_pt[v] in the list
                                               *  above */
            markset       diff_updated_color;
            workspace     diff_previous_color;
            workspace     diff_original_size;
            bool diff_diverge = false; /**< paths of difference are diverging, can not properly track it anymore, and
                                         *there is also no point in tracking it anymore since they are not isomorphic */
            int singleton_pt_start = 0; /**< a pointer were we need to start reading singletons */

            // hooks for color refinement -- this controller hooks into the color refinement algorithm to gather
            // necessary information
            std::function<type_split_color_hook>     my_split_hook; /**< split hook for color refinement */
            std::function<type_worklist_color_hook>  my_worklist_hook; /**< worklist hook for color refinement */

            // settings
            ir_mode mode = IR_MODE_RECORD_TRACE; /**< which mode are we operating in currently? */

            // internal flags for heuristics
            bool s_cell_active     = false; /**< is a cell supposed to be active right now? */
            bool s_individualize   = false; /**< are we individualizing right now? */

            bool h_trace_early_out = false;       /**< use trace early out                */
            bool h_deviation_inc_active  = false; /**< use increased trace deviation      */
            int  s_deviation_inc_current = 0;     /**< number of current trace deviations */
            int  h_deviation_inc = 48; /**< how many additional splits to wait before terminating whenever
                                         * use_increase_deviation is used  */

            bool h_use_split_limit = false;
            int  h_split_limit     = 0;

            int  s_splits = 0;

            /**
             * Marks all colors of the current coloring as "touched". Used on the initial coloring, since we never want
             * to "revert" these colors.
             */
            void touch_initial_colors() {
                int i = 0;
                while (i < c->domain_size) {
                    touched_color.set(i);
                    i += c->ptn[i] + 1;
                }
            }

            /**
             * Resets all information of touched colors (and touches the initial colors again).
             */
            void reset_touched() {
                touched_color.reset();
                touched_color_list.reset();
                prev_color_list.reset();
                touch_initial_colors();
            }

            /**
             * Makes the current trace the comparison trace (and vice versa).
             */
            void flip_trace() {
                trace* t_flip = T;
                T  = cT;
                cT = t_flip;
            }

            /**
             * The split hook function used for color refinement. Tracks the trace invariant, touched colors and
             * singletons. Provides early out functionality using the trace and/or number of cells.
             */
            bool split_hook(const int old_color, const int new_color, const int new_color_sz) {
                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE) {
                    // write singletons to singleton list
                    if (new_color_sz == 1) singletons.push_back(c->lab[new_color]);

                    // record colors that were changed
                    if (!touched_color.get(new_color)) {
                        touched_color.set(new_color);
                        prev_color_list.push_back(old_color);
                        touched_color_list.push_back(new_color);
                    }
                }

                // record split into trace invariant, unless we are individualizing
                if (!s_individualize && old_color != new_color) T->op_refine_cell_record(new_color);

                const bool cont = !h_trace_early_out || T->trace_equal();
                s_deviation_inc_current += (!cont);
                const bool deviation_override = h_deviation_inc_active && (s_deviation_inc_current <= h_deviation_inc);

                const bool continue_cell_limit = !((mode != IR_MODE_RECORD_TRACE) && T->trace_equal()
                                                   && s_base_pos - 1 < static_cast<int>(compare_base->size())
                                                   && (*compare_base)[s_base_pos - 1].cells == c->cells);
                ++s_splits;
                const bool continue_split_limit = !h_use_split_limit || (s_splits < h_split_limit);
                return continue_split_limit && continue_cell_limit && (cont || deviation_override);
            }

            std::function<type_split_color_hook> self_split_hook() {
                #ifndef dej_nolambda
                return [this](auto && PH1, auto && PH2, auto && PH3) { return
                            this->split_hook(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
                                std::forward<decltype(PH3)>(PH3));
                };
                #else
                return std::bind(&controller::split_hook, this, std::placeholders::_1, std::placeholders::_2,
                                 std::placeholders::_3);
                #endif
            }


            /**
             * The worklist hook function used for color refinement. Handles blueprints.
             */
            bool worklist_hook(const int color, const int color_sz) {
                if (s_cell_active) {
                    T->op_refine_cell_end();
                    s_cell_active = false;
                }

                T->op_refine_cell_start(color);
                if (!s_individualize) T->op_additional_info(color_sz);

                // blueprints
                /*if (T->trace_equal() && !T->blueprint_is_next_cell_active()) {
                    s_cell_active = false;
                    T->blueprint_skip_to_next_cell();
                    return false;
                }*/

                s_cell_active = true;
                return true;
            }

            std::function<type_worklist_color_hook> self_worklist_hook() {
                #ifndef dej_nolambda
                return [this](auto && PH1, auto && PH2) {
                    return this->worklist_hook(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2));
                };
                #else
                return std::bind(&controller::worklist_hook, this, std::placeholders::_1, std::placeholders::_2);
                #endif
            }

            /**
             * Vertex v is now differing.
             * @param v the vertex
             */
            void add_diff_vertex(int v) {
                //assert(trace || (!diff_vertices.get(v) && !diff_is_singleton.get(v)));
                if(!diff_vertices.get(v) && !diff_is_singleton.get(v)) {
                    diff_vertices.set(v);
                    diff_vertices_list.push_back(v);
                    diff_vertices_list_pt[v] = diff_vertices_list.cur_pos - 1;
                }
            }

            /**
             * Vertex v is now not differing anymore
             * @param v the vertex
             */
            void remove_diff_vertex(int v) {
                //assert(trace || (diff_vertices.get(v)));
                if(diff_vertices.get(v)) {
                    const int pt      = diff_vertices_list_pt[v];
                    dej_assert(diff_vertices_list[pt] == v);

                    // which element is in the back?
                    const int back_pt  = diff_vertices_list.cur_pos-1;
                    const int back_v = diff_vertices_list[back_pt];

                    // now swap back_col to col
                    diff_vertices_list[pt] = back_v;
                    dej_assert(diff_vertices.get(back_v));
                    dej_assert(diff_vertices_list_pt[back_v] == back_pt);
                    diff_vertices_list_pt[back_v] = pt;
                    diff_vertices_list_pt[v] = -1; // I don't trust myself

                    diff_vertices.unset(v);
                    --diff_vertices_list.cur_pos;
                }
            }

            /**
             * Adds and removes differences of the given color \p col between this state and the given \p state.
             * @param state the other state
             * @param col the coloring to check differences on
             */
            void color_fix_difference(const controller& state, int col) {
                const int col_sz = c->ptn[col] + 1;

                if(col_sz == 1) {
                    add_diff_vertex(state.c->lab[col]);
                    const int v = c->lab[col];
                    remove_diff_vertex(v);
                    if(!diff_is_singleton.get(v)) {
                        diff_is_singleton.set(v);
                    }
                    return;
                }

                diff_tester.reset();
                for(int i = 0; i < col_sz; ++i) diff_tester.set(state.c->lab[col + i]);
                for(int i = 0; i < col_sz; ++i)
                    if(!diff_tester.get(c->lab[col + i])) add_diff_vertex(c->lab[col + i]);

                diff_tester.reset();
                for(int i = 0; i < col_sz; ++i) diff_tester.set(c->lab[col + i]);
                for(int i = 0; i < col_sz; ++i)
                    if(!diff_tester.get(state.c->lab[col + i])) add_diff_vertex(state.c->lab[col + i]);
            }

        public:
            /**
             * Initialize this controller using a refinement and a graph coloring for the initial state.
             *
             * @param ref The refinement workspace to use.
             * @param col The initial coloring.
             */
            controller(refinement* ref, coloring *col) {
                this->c = col;
                this->R = ref;

                T  = &internal_T1;
                cT = &internal_T2;

                touched_color.initialize(col->domain_size);
                touched_color_list.allocate(col->domain_size);
                prev_color_list.allocate(col->domain_size);

                touch_initial_colors();

                my_split_hook = self_split_hook();
                my_worklist_hook = self_worklist_hook();

                diff_tester.initialize(col->domain_size);

                diff_vertices.initialize(col->domain_size);
                diff_vertices_list.resize(col->domain_size);
                diff_vertices_list_pt.resize(col->domain_size);
                diff_is_singleton.initialize(col->domain_size);

                diff_updated_color.initialize(col->domain_size);
                diff_previous_color.resize(col->domain_size);
                diff_original_size.resize(col->domain_size);

                singleton_pt_start = 0;
                singletons.reserve(col->domain_size);
            }

            int get_number_of_splits() {
                return s_splits;
            }

            /**
             * Sets internal trace into recording mode. We write a trace which we might want to compare to later.
             *
             * Always reversible (see also \a use_reversible).
             */
            void mode_write_base() {
                T->reset();
                cT->reset();

                reset_touched();
                mode = IR_MODE_RECORD_TRACE;
                T->set_compare(false);
                T->set_record(true);

                internal_compare_base.clear();
                internal_compare_base_vertex.clear();
                internal_compare_singletons.clear();
            }

            /**
             * We compare to a pre-existing trace, recorded earlier using mode_write_base. Must use \a compare_to_this or
             * provide a comparison trace in another manner before being able to use this mode in the intended manner.
             */
            void mode_compare_base() {
                mode = ir::IR_MODE_COMPARE_TRACE_REVERSIBLE;
            }

            /**
             * Checks whether the current coloring matches the coloring on the base (AKA in saucy terms, whether this
             * ordered partition and the one on the base are matched).
             *
             * @return whether the colorings match
             */
            bool there_is_difference_to_base() {
                for(int i = 0; i < c->domain_size;) {
                    const int col    = i;
                    const int col_sz = c->ptn[col] + 1;
                    i += col_sz;
                    if(col_sz == 1) continue;

                    diff_tester.reset();
                    for(int j = col; j < col + col_sz; ++j) {
                        const int v = c->lab[j];
                        diff_tester.set(v);
                    }
                    for(int j = col; j < col + col_sz; ++j) {
                        const int v = leaf_color.lab[j];
                        if(!diff_tester.get(v)) return true;
                    }
                }
                return false;
            }

            /**
             * Make an automorphism with the current singletons and the singletons of the base, i.e., useful when
             * \a there_is_difference_to_base returns `false`.
             *
             * @param automorphism workspace to write the automorphism to
             */
            void singleton_automorphism_base(groups::automorphism_workspace*  automorphism) {
                automorphism->reset();

                for(int i = 0; i < c->domain_size;) {
                    const int col_sz =  c->ptn[i] + 1;

                    if(col_sz == 1 && c->lab[i] != leaf_color.lab[i]) {
                        automorphism->write_single_map(c->lab[i], leaf_color.lab[i]);
                    }
                    i += col_sz;
                }
            }

            void color_diff_automorphism_base(groups::automorphism_workspace*  automorphism) {
                automorphism->reset();

                for(int i = 0; i < c->domain_size; ++i) {

                    const int v1 = i;
                    const int col = c->vertex_to_lab[i];
                    const int v2  = leaf_color.lab[col];

                    if(v1 != v2) {
                        automorphism->write_single_map(v1, v2);
                    }
                }
            }


            /**
             * Checks whether the current coloring matches the coloring on the base, including singletons (AKA, whether
             * the colorings are equal). Useful for debugging.
             *
             * @return whether the colorings are equal
             */
            bool there_is_difference_to_base_including_singles(int domain_size) {
                for(int i = 0; i < domain_size;) {
                    const int col    = i;
                    const int col_sz = c->ptn[col] + 1;
                    i += col_sz;

                    diff_tester.reset();
                    for(int j = col; j < col + col_sz; ++j) {
                        const int v = c->lab[j];
                        diff_tester.set(v);
                    }
                    for(int j = col; j < col + col_sz; ++j) {
                        const int v = leaf_color.lab[j];
                        if(!diff_tester.get(v)) return true;
                    }
                }
                return false;
            }

            /**
             * Need to use difference-checking before. Returns a differing pair of vertices of the same non-singleton
             * color.
             *
             * @param state the other state (for which we've checked differences before)
             * @return
             */
            std::pair<int, int> diff_pair(const controller& state) {
                dej_assert(diff_vertices_list.cur_pos > 0);

                // pick a vertex from the difference list
                const int right = diff_vertices_list[0];
                const int right_col = c->vertex_to_col[right];

                // find a corresponding counterpart in the other coloring
                //const int left      = state.c->lab[right_col];
                const int right_col_sz = c->ptn[right_col] + 1;
                int left = -1;
                for (int i = 0; i < right_col_sz; ++i) {
                    const int candidate = state.c->lab[right_col + i];
                    if(diff_vertices.get(candidate)) {
                        left = candidate;
                        break;
                    }
                }
                if(left == -1) left = state.c->lab[right_col];


                return {left, right};
            }

            /**
             * Resets all difference-checking information.
             */
            void reset_diff() {
                diff_vertices.reset();
                diff_vertices_list.reset();
                diff_is_singleton.reset();

                diff_tester.reset();
                singleton_pt_start = singletons.size();

                diff_diverge = false;
            }

            /**
             * Need to use difference-checking before. Makes a singleton automorphism with the singleton difference
             * between this state and the given \p state.
             *
             * @param state the other state
             * @param automorphism workspace to write the automorphism to
             */
            void singleton_automorphism(controller& state, groups::automorphism_workspace& automorphism) {
                automorphism.reset();
                for(int i = singleton_pt_start; i < static_cast<int>(singletons.size()); ++i) {
                    automorphism.write_single_map(singletons[i], state.singletons[i]);
                }
            }

            /**
             * Need to use difference-checking before. Returns whether the difference is diverging, i.e., there is no
             * point in further difference checking since branches are non-isomorphic.
             *
             * @return whether states are diverging
             */
            dej_nodiscard int get_diff_diverge() const {
                return diff_diverge;
            }


            /**
             * Difference-checking. Records differences between this state and the given \p other_state. Before using
             * this, this state and \p other_state should start from the same node of the IR tree, and the difference
             * should be reset at that point. Then, after each individualization this update function should be called.
             * The function keeps track of all colors that are non-matching.
             *
             * @param other_state the other state
             * @return whether there is a difference
             */
            bool update_diff_vertices_last_individualization(const controller& other_state) {
                int i;
                if(touched_color_list.cur_pos != other_state.touched_color_list.cur_pos) {
                    diff_diverge = true;
                    return true;
                }

                diff_updated_color.reset();

                // reconstruct color prior to individualization for each new color
                for(i = base.back().touched_color_list_pt; i < touched_color_list.cur_pos; ++i) {
                    const int old_color = prev_color_list[i];
                    const int new_color = touched_color_list[i];
                    diff_original_size[old_color] = 0;
                    if(!diff_updated_color.get(old_color)) {
                        diff_updated_color.set(old_color);
                        diff_previous_color[old_color] = old_color;
                    }
                    while(diff_previous_color[old_color] != diff_previous_color[diff_previous_color[old_color]])
                        diff_previous_color[old_color] = diff_previous_color[diff_previous_color[old_color]];
                    diff_updated_color.set(new_color);
                    diff_previous_color[new_color] = diff_previous_color[old_color];
                }

                // determine size of the original color classes
                for(i = base.back().touched_color_list_pt; i < touched_color_list.cur_pos; ++i) {
                    const int new_color = touched_color_list[i];
                    const int old_color = diff_previous_color[new_color];
                    diff_original_size[old_color] = std::max(diff_original_size[old_color],
                                                             new_color + c->ptn[new_color] + 1);
                }

                diff_updated_color.reset();

                // now fix colors according to fragments of original color classes
                for(i = base.back().touched_color_list_pt; i < touched_color_list.cur_pos; ++i) {
                    const int old_color = prev_color_list[i];
                    if(diff_updated_color.get(old_color)) {
                        dej_assert(diff_updated_color.get(touched_color_list[i]));
                        continue;
                    }
                    int old_color_end_pt = std::max(old_color + c->ptn[old_color] + 1, diff_original_size[old_color]);


                    // let's update diff for all the fragments of old_color
                    const int old_color_sz       = c->ptn[old_color] + 1;
                    const int old_color_sz_other = other_state.c->ptn[old_color] + 1;
                    diff_diverge = diff_diverge || (old_color_sz != old_color_sz_other);
                    if(diff_diverge) break;

                    // determine the largest fragment
                    int largest_fragment    = -1;
                    int largest_fragment_sz = -1;

                    for(int j = old_color; j < old_color_end_pt;) {
                        const int next_color    = j;
                        const int next_color_sz = c->ptn[next_color] + 1;
                        const int next_color_sz_other = other_state.c->ptn[next_color] + 1;
                        diff_diverge = diff_diverge || (next_color_sz != next_color_sz_other);
                        if(next_color_sz > largest_fragment_sz) {
                            largest_fragment    = next_color;
                            largest_fragment_sz = next_color_sz;
                        }

                        diff_updated_color.set(next_color);
                        j += next_color_sz;
                    }

                    if(diff_diverge) break;

                    // update all but the largest of those fragments
                    for(int j = old_color; j < old_color_end_pt;) {
                        const int next_color    = j;
                        const int next_color_sz = c->ptn[next_color] + 1;
                        if(j != largest_fragment) color_fix_difference(other_state, next_color);
                        j += next_color_sz;
                    }

                    // unless the largest fragment is a singleton, in that case we want to record the singleton
                    if(largest_fragment_sz == 1) color_fix_difference(other_state, largest_fragment);
                }
                return (diff_vertices_list.cur_pos != 0) || diff_diverge; // there is a diff!
            }


            /**
             * Copy the state of another controller. Does not copy the trace, but instead links the trace of this
             * controller to the comparison state of the other (such that both states compare to the same single
             * trace).
             *
             * @param state the other state which is copied
             */
            void link_compare(controller* state) {
                this->c->copy_any(state->c);

                T->set_compare(true);
                T->set_record(false);
                T->set_compare_trace(state->T->get_compare_trace());
                T->set_position(state->T->get_position());
                T->set_hash(state->T->get_hash());

                base_vertex = state->base_vertex;
                base        = state->base;

                compare_base_vertex = state->compare_base_vertex;
                compare_base        = state->compare_base;

                leaf_color.copy_any(&state->leaf_color);
                mode = state->mode;

                touched_color.copy(&state->touched_color);
                prev_color_list.copy(&state->prev_color_list);
                touched_color_list.copy(&state->touched_color_list);

                //singletons         = state->singletons;
                compare_singletons = state->compare_singletons;

                s_base_pos = state->s_base_pos;

                this->R = state->R;
            }

            /**
             * Copy the state of another controller. Does not copy the trace, but instead links the trace of this
             * controller to the comparison state of the other (such that both states compare to the same single
             * trace).
             *
             * @param state the other state which is copied
             */
            void link(controller* state) {
                this->c->copy_any(state->c);

                T->set_compare(true);
                T->set_record(false);
                T->set_compare_trace(state->T);
                T->set_position(state->T->get_position());
                T->set_hash(state->T->get_hash());

                base_vertex = state->base_vertex;
                base        = state->base;

                compare_base_vertex = &state->base_vertex;
                compare_base        = &state->base;

                leaf_color.copy_any(&state->leaf_color);
                mode = state->mode;

                touched_color.copy(&state->touched_color);
                prev_color_list.copy(&state->prev_color_list);
                touched_color_list.copy(&state->touched_color_list);
                singletons         = state->singletons;
                compare_singletons = &state->singletons;

                s_base_pos = state->s_base_pos;

                this->R = state->R;
            }

            int diff_num() {
                return diff_vertices_list.cur_pos;
            };

            /**
             * Compare all following computations to this IR node. The \a mode must be set to `IR_MODE_RECORD_TRACE`
             * (using \a mode_write_base) to call this function. Changes \a mode to `IR_MODE_COMPARE_TRACE_REVERSIBLE`
             * (i.e., such as calling \a mode_compare_base).
             */
            void compare_to_this() {
                dej_assert(mode == ir::IR_MODE_RECORD_TRACE);

                cT->set_compare(true);
                cT->set_record(false);
                cT->set_compare_trace(T);
                cT->set_position(T->get_position());
                flip_trace();

                internal_compare_base        = base;
                internal_compare_base_vertex = base_vertex;
                internal_compare_singletons  = singletons;

                compare_base        = &internal_compare_base;
                compare_base_vertex = &internal_compare_base_vertex;
                compare_singletons  = &internal_compare_singletons;


                leaf_color.copy_any(c);
                mode = ir::IR_MODE_COMPARE_TRACE_REVERSIBLE;
            }

            void reserve() {
                const int domain_size = sqrt(c->domain_size);
                int sqrt_domain = sqrt(domain_size);
                base.reserve(sqrt_domain);
                base_vertex.reserve(sqrt_domain);
                T->reserve(2*domain_size);
            }

            /**
             * Enables or disables whether following \a move_to_child calls will be reversible using \a move_to_parent,
             * or not. Using non-reversible calls to \a move_to_child is faster.
             *
             * When recording a trace using \a mode_write_base, computations are always reversible and calling this
             * function will have no effect.
             * @param reversible
             */
            void use_reversible(const bool reversible) {
                if(mode == IR_MODE_RECORD_TRACE) return;
                if(reversible) {
                    reset_touched();
                    mode = IR_MODE_COMPARE_TRACE_REVERSIBLE;
                } else {
                    mode = IR_MODE_COMPARE_TRACE_IRREVERSIBLE;
                }
            }

            /**
             * Whether to terminate color refinement whenever a deviation to its comparison trace is found.
             *
             * @param trace_early_out Flag that determines whether the early out is used.
             */
            void use_trace_early_out(bool trace_early_out) {
                this->h_trace_early_out = trace_early_out;
                if(!trace_early_out) use_increase_deviation(false);
            }

            /**
             * Whether to record additional deviation information once a deviation from the comparison trace is found.
             * Only applicable when \a use_trace_early_out is set. Essentially delays the termination of color refinement
             * to record more information into the trace invariant.
             *
             * @param deviation_inc_active Whether increased deviation is recorded or not.
             */
            void use_increase_deviation(bool deviation_inc_active) {
                h_deviation_inc_active = deviation_inc_active;
            }

            /**
             * Whether to record additional deviation information once a deviation from the comparison trace is found.
             * Only applicable when \a use_trace_early_out is set. Essentially delays the termination of color refinement
             * to record more information into the trace invariant.
             *
             * @param deviation_inc how many deviations to add before terminating (default = 48)
             */
            void set_increase_deviation(int deviation_inc = 48) {
                h_deviation_inc = deviation_inc;
            }

            /**
             * Whether to record additional deviation information once a deviation from the comparison trace is found.
             * Only applicable when \a use_trace_early_out is set. Essentially delays the termination of color
             * refinement to record more information into the trace invariant.
             *
             * @param deviation_inc_active Whether increased deviation is recorded or not.
             */
            void use_split_limit(bool use_split_limit, int limit = 0) {
                h_use_split_limit = use_split_limit;
                h_split_limit = limit;
            }

            /**
             * Resets whether the trace is deemed equal to its comparison trace.
             */
            void reset_trace_equal() {
                T->reset_trace_equal();
                s_deviation_inc_current = 0;
            }

            /**
             * Save a partial state of this controller.
             *
             * @param state A reference to the limited_save in which the state will be stored.
             */
            void save_reduced_state(limited_save &state) {
                state.save(base_vertex, *c, T->get_hash(), T->get_position(),
                           s_base_pos);
            }

            /**
             * Load a partial state into this controller.
             *
             * @param state A reference to the limited_save from which the state will be loaded.
             */
            void load_reduced_state(limited_save &state) {
                c->copy_any(state.get_coloring());

                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                s_base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless, so clear them out
                base.clear();

                // if reversible, need to reset touched colors
                if(mode == IR_MODE_COMPARE_TRACE_REVERSIBLE) reset_touched();
            }

            void load_reduced_state_without_coloring(limited_save &state) {
                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                s_base_pos  = state.get_base_position();
                base_vertex = state.get_base();

                // becomes meaningless
                base.clear();

                // deactivate reversability
                mode = IR_MODE_COMPARE_TRACE_IRREVERSIBLE;
            }

            dej_nodiscard coloring *get_coloring() const {
                return c;
            }

            dej_nodiscard int get_base_pos() const {
                return s_base_pos;
            }

            /**
             * Write additional information into the internal trace. Care must be taken that the written data is
             * isomorphism-invariant.
             *
             * @param d Data to be written to the trace.
             */
            void write_to_trace(const int d) const {
                T->op_additional_info(d);
            }

            /**
             * Writes a stronger invariant using the internal coloring and graph to the trace.
             *
             * @param g The graph.
             */
            void write_strong_invariant(const sgraph* g) const {
                for(int l = 0; l < g->v_size; ++l) {
                    const int v = c->lab[l];
                    unsigned int inv1 = 0;
                    const int start_pt = g->v[v];
                    const int end_pt   = start_pt + g->d[v];
                    for(int pt = start_pt; pt < end_pt; ++pt) {
                        const int other_v = g->e[pt];
                        inv1 += hash((unsigned int) c->vertex_to_col[other_v]);
                    }
                    T->op_additional_info((int) inv1);
                }
            }

            /**
             * Writes a stronger invariant using the internal coloring and graph to the trace, but less strong than the
             * one written by \a write_strong_invariant.
             *
             * @param g The graph.
             */
            void write_strong_invariant_quarter(const sgraph* g) const {
                for(int l = 0; l < g->v_size; l += 4) {
                    const int v = c->lab[l];
                    unsigned int inv1 = 0;
                    const int start_pt = g->v[v];
                    const int end_pt   = start_pt + g->d[v];
                    for(int pt = start_pt; pt < end_pt; ++pt) {
                        const int other_v = g->e[pt];
                        inv1 += hash((unsigned int) c->vertex_to_col[other_v]);
                    }
                    T->op_additional_info((int) inv1);
                }
            }

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child(sgraph *g, int v) {
                // always keep track of vertex base
                ++s_base_pos;
                s_splits = 0;
                base_vertex.push_back(v);

                dej_assert(!s_cell_active);

                // some info maybe needed later
                const int singleton_pt     = (int) singletons.size();
                const int touched_color_pt = touched_color_list.cur_pos;
                const int trace_pos        = T->get_position();
                const uint64_t trace_hash  = T->get_hash();

                // determine color
                const int prev_col    = c->vertex_to_col[v];
                const int prev_col_sz = c->ptn[prev_col] + 1;

                // individualize vertex
                s_individualize = true;
                const int init_color_class = refinement::individualize_vertex(c, v, my_split_hook);
                T->op_individualize(c->vertex_to_col[v]);
                s_individualize = false;

                // refine coloring
                T->op_refine_start();
                if (mode == IR_MODE_RECORD_TRACE) {
                    R->refine_coloring(g, c, init_color_class, -1, &my_split_hook, &my_worklist_hook);
                    if (s_cell_active) T->op_refine_cell_end();
                    T->op_refine_end();
                } else {
                    R->refine_coloring(g, c, init_color_class, -1, &my_split_hook, &my_worklist_hook);
                    if (T->trace_equal() && c->cells==(*compare_base)[s_base_pos - 1].cells && !h_use_split_limit) {
                        T->skip_to_individualization();
                    }
                }

                s_cell_active = false;

                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE)
                    base.emplace_back(prev_col, prev_col_sz, c->cells, touched_color_pt, singleton_pt, trace_pos,
                                      trace_hash);
            }

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child_no_trace(sgraph *g, int v) {
                const int init_color_class = R->individualize_vertex(c, v);
                R->refine_coloring_first(g, c, init_color_class);
            }

            /**
             * Perform color refinement on the internal coloring based on the given graph.
             *
             * @param g The graph.
             */
            void refine(sgraph *g) {
                R->refine_coloring_first(g, c);
            }

            /**
             * Move IR node kept in this controller back to its parent.
             */
            void move_to_parent() {
                dej_assert(mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE);
                dej_assert(base.size() > 0);

                // unwind invariant
                T->set_position(base.back().trace_pos);
                T->set_hash(base.back().trace_hash);

                // unwind colors introduced on this level of the tree
                --s_base_pos;
                while (prev_color_list.cur_pos > base.back().touched_color_list_pt) {
                    const int old_color = prev_color_list.pop_back();
                    const int new_color = touched_color_list.pop_back();

                    touched_color.unset(new_color);

                    const int new_color_sz = c->ptn[new_color] + 1;
                    c->ptn[old_color] += new_color_sz;
                    c->ptn[new_color] = 1;

                    for (int j = 0; j < new_color_sz; ++j) {
                        const int v = c->lab[new_color + j];
                        c->vertex_to_col[v] = old_color;
                        dej_assert(c->vertex_to_lab[v] == new_color + j);
                    }

                    --c->cells;
                }
                // unwind singletons
                int const new_singleton_pos = base.back().singleton_pt;
                singletons.resize(new_singleton_pos);

                // remove base information
                base_vertex.pop_back();
                base.pop_back();

                T->reset_trace_equal();
            }

            /**
             * Perform the given walk from the given IR node.
             *
             * @param g The graph.
             * @param start_from A limited_save describing the IR node from which the walk should be performed.
             * @param vertices A vector of vertices that describes the walk, i.e., these vertices will be individualized
             *                 in the given order.
             */
            void walk(sgraph *g, ir::limited_save &start_from, std::vector<int>& vertices) {
                load_reduced_state(start_from);

                while (s_base_pos < (int) vertices.size()) {
                    int v = vertices[s_base_pos];
                    move_to_child(g, v);
                }
            }

            /**
             * Certifies whether the provided \p automorphism is indeed an automorphism of the graph \p g. Runs in
             * time roughly in the support of \p automorphism.
             *
             * Uses internal workspace to perform this operation efficiently.
             *
             * @param g The graph.
             * @param automorphism The automorphism.
             * @return Whether \p automorphism is an automorphism of \p g.
             */
            bool certify(sgraph* g, groups::automorphism_workspace& automorphism) {
                if(automorphism.nsupp() == 0) return true;
                if(automorphism.nsupp() > g->v_size / 4) {
                    return R->certify_automorphism(g, automorphism.p());
                } else {
                    return R->certify_automorphism_sparse(g, automorphism.p(), automorphism.nsupp(),
                                                          automorphism.supp());
                }
            }
        };

        /**
         * \brief Creates cell selectors.
         *
         * Heuristics which enable the creation of different cell selectors, as well as moving an \ref controller to a
         * leaf of the IR tree.
         */
        class cell_selector_factory {
            const int locked_lim = 512;

            std::vector<int> saved_color_base;
            std::vector<int> saved_color_sizes;
            std::function<type_selector_hook> dynamic_seletor;
            markset test_set;
            std::vector<int> candidates;
            big_number ir_tree_size_estimate;

            int color_score(sgraph *g, controller *state, int color) {
                test_set.reset();
                const int v = state->get_coloring()->lab[color];
                const int d = g->d[v];
                const int ept_st  = g->v[v];
                const int ept_end = ept_st + d;
                int non_triv_col_d = 1;
                for (int i = ept_st; i < ept_end; ++i) {
                    const int test_col = state->get_coloring()->vertex_to_col[g->e[i]];
                    if (!test_set.get(test_col)) {
                        non_triv_col_d += 1;
                        test_set.set(test_col);
                    }
                }
                return state->get_coloring()->ptn[color] * non_triv_col_d;
            }

            static int color_score_size(controller *state, int color) {
                return state->get_coloring()->ptn[color];
            }

            static int color_score_anti_size(controller *state, int color) {
                return INT32_MAX - state->get_coloring()->ptn[color];
            }

        public:
            enum selector_style {SELECT_SPARSE, SELECT_COMB, SELECT_SMALL, SELECT_DRY_LAND};

            /**
             * Cell selector, chooses first non-trivial color class, unless color class from stored base is
             * applicable.
             *
             * @param c Coloring from which a color class shall be selected.
             * @param base_pos Current position in base.
             */
            int cell_selector(const coloring *c, const int base_pos) {
                if (base_pos >= 0 && base_pos < (int) saved_color_base.size() && c->ptn[saved_color_base[base_pos]] > 0 &&
                    c->vertex_to_col[c->lab[saved_color_base[base_pos]]] == saved_color_base[base_pos] &&
                        c->ptn[saved_color_base[base_pos]] + 1 == saved_color_sizes[base_pos]) {
                    return saved_color_base[base_pos];
                }
                int best_score = -1;
                int best_color = -1;
                for (int i = 0; i < c->domain_size;) {
                    if (c->ptn[i] > best_score && c->ptn[i] > 0) {
                        best_score = c->ptn[i];
                        best_color = i;
                    }
                    if(best_color >= 0) break;
                    i += c->ptn[i] + 1;
                }
                return best_color;
            }

            /**
             * @return A selector hook based on the saved base and configured dynamic selector.
             */
            std::function<type_selector_hook> *get_selector_hook() {
                dynamic_seletor = std::bind(&cell_selector_factory::cell_selector, this, std::placeholders::_1,
                                            std::placeholders::_2);
                return &dynamic_seletor;
            }

            /**
             * @return estimate for how large the IR tree of the last computed selector is
             */
            big_number get_ir_size_estimate() {
                return ir_tree_size_estimate;
            }

            /**
             * Creates and stores a new cell selector.
             *
             * @param g the graph
             * @param state ir controller that is navigated to the new target leaf
             * @param state_probe ir controller used for auxiliary probing
             * @param h_seed chooses strategy of the new cell selector
             * @param h_budget available budget
             */
            void make_cell_selector(sgraph *g, controller *state, controller *state_probe, const selector_style style,
                                    const int h_seed, const int h_budget) {
                int perturbe = 0;
                if(h_seed > 6) {
                    perturbe = (int) (hash(h_seed) % 256);
                }
                switch(style) {
                    case 0:
                        //find_combinatorial_optimized_base_recurse(g, state, perturbe);
                        find_sparse_optimized_base(g, state);
                        //find_first_base(g, state);
                        break;
                    case 1:
                        find_combinatorial_optimized_base_recurse(g, state, perturbe);
                        break;
                    case 2:
                        find_small_optimized_base(g, state, perturbe);
                        break;
                    case 3:
                        if(h_budget > 32) {
                            // hail mary selector... seems to be very good for some classes, but is very expensive and
                            // its use must be limited somehow... probing is only used for a fixed number of levels and
                            // then switches to find_combinatorial_optimized_base
                            // idea is to create trace deviations "higher up" in the tree, such that bfs may succeed
                            // earlier
                            find_early_trace_deviation_base(g, state, state_probe, perturbe);
                        } else find_combinatorial_optimized_base(g, state);
                        break;
                }

                // now save the selector and make some estimates
                ir_tree_size_estimate.mantissa = 1.0;
                ir_tree_size_estimate.exponent = 0;

                saved_color_base.clear();
                saved_color_base.reserve(state->base.size());
                for(auto& bi : state->base) {
                    saved_color_base.push_back(bi.color);
                }

                saved_color_sizes.clear();
                saved_color_sizes.reserve(state->base.size());
                for(auto& bi : state->base) {
                    ir_tree_size_estimate.multiply(bi.color_sz);
                    saved_color_sizes.push_back(bi.color_sz);
                }
            }
            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for simple,
             * sparse graphs.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_sparse_optimized_base(sgraph *g, controller *state) {
                // some settings for heuristics
                constexpr int base_lim = 1000;
                constexpr int test_lim_pre  = 512;
                constexpr int test_lim_post = 1;

                state->mode_write_base();

                int prev_color = -1;
                dej_assert(state->s_base_pos == 0);
                test_set.initialize(g->v_size);

                int start_test_from = 0;
                int start_test_from_inside = 0;

                int buffered_col       = -1;
                int buffered_col_score = -1;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = state->s_base_pos <= base_lim?test_lim_pre:test_lim_post;

                    // pick previous color if possible
                    if (prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                        best_color = prev_color;
                    } else if (prev_color >= 0) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        int i =  g->v[test_vertex] + start_test_from_inside;
                        const int end_pt = g->v[test_vertex] + g->d[test_vertex];
                        for (; i < end_pt; ++i) {
                            const int other_vertex = g->e[i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0) {
                                best_color = other_color;
                                break;
                            } else {
                                ++start_test_from_inside;
                            }
                        }
                    }

                    start_test_from_inside = 0;

                    // use buffered color if possible
                    if(best_color == -1 && buffered_col >= 0 && state->get_coloring()->ptn[buffered_col] > 0 &&
                       color_score(g, state, buffered_col) >= buffered_col_score) {
                        best_color = buffered_col;
                        buffered_col = -1;
                        buffered_col_score = -1;
                    }

                    // same color, neighbour of color and buffer failed, so now we try to pick a "good" color
                    if (best_color == -1) {
                        int best_score = -1;

                        for (int i = start_test_from; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                            const int col_sz = state->get_coloring()->ptn[i];

                            if (col_sz > 0) {
                                --test_lim;
                                const int test_score = color_score(g, state, i);
                                if (test_score > best_score) {
                                    best_color = i;
                                    best_score = test_score;
                                    buffered_col       = -1;
                                    buffered_col_score = -1;
                                } else if (test_score == best_score) {
                                    buffered_col       = i;
                                    buffered_col_score = test_score;
                                }
                            }

                            if(col_sz == 0 && start_test_from == i) start_test_from += col_sz + 1;
                            i += col_sz + 1;
                        }
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }
            }

            void find_first_base(sgraph *g, controller *state) {
                state->mode_write_base();
                dej_assert(state->s_base_pos == 0);

                test_set.initialize(g->v_size);

                int start_test_from = 0;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = 512;


                    for (int i = start_test_from; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                        const int col_sz = state->get_coloring()->ptn[i];

                        if (col_sz > 0) {
                                best_color = i;
                                break;
                        }

                        if(col_sz == 0 && start_test_from == i) start_test_from += col_sz + 1;
                        i += col_sz + 1;
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }
            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for simple,
             * sparse graphs.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_sparse_optimized_base_recurse(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                test_set.initialize(g->v_size);
                int prev_color    = -1;
                int prev_color_sz = 0;

                dej_assert(state->s_base_pos == 0);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = 512;

                    // recurse into previous color if possible
                    if(prev_color >= 0) {
                        for (int i = 0; i < prev_color + prev_color_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                best_color = i;
                                break;
                            }
                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    if (prev_color >= 0 && best_color == -1) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        for (int i = 0; i < g->d[test_vertex]; ++i) {
                            const int other_vertex = g->e[g->v[test_vertex] + i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0) {
                                best_color = other_color;
                            }
                        }
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        int best_score = -1;

                        for (int i = 0; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                int test_score = color_score(g, state, i);
                                if (test_score > best_score) {
                                    best_color = i;
                                    best_score = test_score;
                                }
                                --test_lim;
                            }

                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    prev_color    = best_color;
                    prev_color_sz = state->get_coloring()->ptn[best_color] + 1;
                    state->move_to_child(g, state->get_coloring()->lab[(best_color +
                                                                        (perturbe%state->get_coloring()->ptn[best_color]))]);
                }
            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attempts to find a with small color
             * classes.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_small_optimized_base(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim = 256;

                    // heuristic, try to pick "good" color
                    int best_score = INT32_MIN;
                    for (int i = 0; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            int test_score = color_score_anti_size(state, i);
                            if (test_score > best_score || best_color == -1) {
                                best_color = i;
                                best_score = test_score;
                                --test_lim;
                                if(best_score == INT32_MAX - 1) break;
                            }
                        }

                        i += state->get_coloring()->ptn[i] + 1;
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    state->move_to_child(g, state->get_coloring()->lab[(best_color +
                                                                    (perturbe%state->get_coloring()->ptn[best_color]))]);
                }
            }

            void find_combinatorial_optimized_base(sgraph *g, controller *state) {
                state->mode_write_base();

                test_set.initialize(g->v_size);
                candidates.clear();
                candidates.reserve(locked_lim);

                markset neighbour_color;
                neighbour_color.initialize(g->v_size);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // heuristic, try to pick "good" color
                    candidates.clear();
                    int best_score = -1;
                    for (int i = 0; i < state->get_coloring()->domain_size;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            candidates.push_back(i);
                        }
                        i += state->get_coloring()->ptn[i] + 1;
                    }
                    while (!candidates.empty()) {
                        const int test_color = candidates.back();
                        candidates.pop_back();

                        int test_score = color_score_size(state, test_color);
                        if (neighbour_color.get(test_color)) {
                            test_score *= 10;
                        }
                        if (test_score >= best_score) {
                            best_color = test_color;
                            best_score = test_score;
                        }
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }
            }



            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for
             * combinatorial graphs solved by bfs/random walks (i.e., by choosing large colors).
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_combinatorial_optimized_base_recurse(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                int prev_color    = -1;
                int prev_color_sz = 0;

                const bool perturbe_flip = perturbe % 2;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int best_score = -1;

                    // recurse into previous color if possible
                    if(prev_color >= 0) {
                        for (int i = 0; i < prev_color + prev_color_sz;) {
                            const int col_sz = state->get_coloring()->ptn[i];
                            int test_score = color_score_size(state, i);
                            if (((test_score > best_score) || (perturbe_flip && (test_score >= best_score))) && col_sz > 0) {
                                best_color = i;
                                best_score = test_score;
                            }
                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    if(best_color == -1) {
                        // heuristic, pick first largest color
                        for (int i = 0; i < state->get_coloring()->domain_size;) {
                            const int col_sz = state->get_coloring()->ptn[i];
                            int test_score = color_score_size(state, i);
                            if (test_score > best_score && col_sz > 0) {
                                best_color = i;
                                best_score = test_score;
                            }

                            i += col_sz + 1;
                        }
                    }

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    prev_color    = best_color;
                    prev_color_sz = state->get_coloring()->ptn[best_color] + 1;
                    state->move_to_child(g, state->get_coloring()->lab[(best_color + (perturbe% prev_color_sz))]);
                }

            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for
             * combinatorial graphs solved by bfs/random walks (i.e., by choosing large colors).
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_early_trace_deviation_base(sgraph *g, controller* state, controller* state_probe, int perturbe) {
                state->mode_write_base();
                constexpr int probe_limit = 5;

                state_probe->link(state);
                state_probe->mode_compare_base();
                state_probe->use_reversible(true);
                candidates.clear();

                int probe_limit_candidates = 8;

                bool use_probing = true;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int best_score_deviate = -1;
                    int best_score_cells   = -1;
                    int best_score_size    = -1;
                    int best_test_v        = -1;

                    probe_limit_candidates = std::max(probe_limit_candidates, 2);

                    if(use_probing) {
                        candidates.clear();
                        for (int i = 0; i < state->c->domain_size;) {
                            const int col_sz = state->c->ptn[i] + 1;
                            if (col_sz >= 2) {
                                candidates.push_back(i);
                            }
                            if (static_cast<int>(candidates.size()) >= probe_limit_candidates) break;
                            i += col_sz;
                        }

                        // now, probe and score the candidates... we are trying to find the color with most trace deviations
                        for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
                            const int col = candidates[i];
                            const int col_sz = state->c->ptn[col] + 1;
                            dej_assert(col_sz >= 2);
                            const int probe_lim_col = std::min(probe_limit, col_sz);

                            const int v_base = state->get_coloring()->lab[(col + (perturbe % col_sz))];
                            dej_assert(state->s_base_pos == state_probe->s_base_pos);
                            dej_assert(state->T->get_position() == state_probe->T->get_position());
                            dej_assert(state->c->cells == state_probe->c->cells);
                            #if defined(DEJDEBUG) &&  !defined(NDEBUG)
                            const int cells_pre = state->c->cells;
                            const int previous_pos = state->T->get_position();
                            #endif

                            state->move_to_child(g, v_base);
                            const int cells = state->c->cells;

                            int deviated = 0;
                            for (int j = 0; j < probe_lim_col; ++j) {
                                // pick random vertex
                                const int v = state_probe->c->lab[col + ((j + perturbe) % col_sz)];
                                // individualize in state_probe
                                state_probe->reset_trace_equal();
                                state_probe->use_trace_early_out(true);
                                dej_assert(state_probe->c->vertex_to_col[v] == state_probe->c->vertex_to_col[v_base]);
                                state_probe->move_to_child(g, v);
                                deviated += !state_probe->T->trace_equal();
                                dej_assert(v == v_base ? state_probe->T->trace_equal() : true);
                                state_probe->move_to_parent();
                                dej_assert(state_probe->c->cells == cells_pre);
                                dej_assert(state_probe->T->get_position() == previous_pos);
                            }

                            state->move_to_parent();
                            dej_assert(state->T->get_position() == previous_pos);
                            dej_assert(cells_pre == state->c->cells);
                            dej_assert(state->s_base_pos == state_probe->s_base_pos);
                            dej_assert(state->T->get_position() == state_probe->T->get_position());
                            dej_assert(state->c->cells == state_probe->c->cells);

                            if (deviated > best_score_deviate ||
                               ((deviated == best_score_deviate) && (cells > best_score_cells)) ||
                               ((deviated == best_score_deviate) && (cells == best_score_cells) &&
                                (col_sz > best_score_size))) {
                                best_color = col;
                                best_score_deviate = deviated;
                                best_score_cells   = cells;
                                best_score_size    = col_sz;
                                best_test_v        = v_base;
                                if (deviated > 0) break;
                            }
                        }

                        if (best_score_deviate > 0) --probe_limit_candidates;
                    } else {
                        // heuristic, try to pick "good" color
                        candidates.clear();
                        int best_score = -1;
                        for (int i = 0; i < state->get_coloring()->domain_size;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                candidates.push_back(i);
                            }
                            i += state->get_coloring()->ptn[i] + 1;
                        }
                        while (!candidates.empty()) {
                            const int test_color = candidates.back();
                            candidates.pop_back();

                            int test_score = color_score_size(state, test_color);
                            if (test_score >= best_score) {
                                best_color = test_color;
                                best_score = test_score;
                            }
                        }
                    }

                    if(probe_limit_candidates < 2 || state->s_base_pos > 5) use_probing = false;

                    dej_assert(best_color >= 0);
                    dej_assert(best_color < g->v_size);
                    dej_assert(!use_probing || state->s_base_pos == state_probe->s_base_pos);
                    dej_assert(!use_probing || state->T->get_position() == state_probe->T->get_position());
                    //assert(state->T->get_hash() == state_probe->T->get_hash());
                    dej_assert(!use_probing || state->c->cells == state_probe->c->cells);
                    const int col_sz = state->get_coloring()->ptn[best_color] + 1;
                    int v = best_test_v;
                    if(v == -1) v = state->get_coloring()->lab[(best_color + (perturbe% col_sz))];

                    dej_assert(!use_probing || state->c->vertex_to_col[v] == state_probe->c->vertex_to_col[v]);

                    state_probe->reset_trace_equal();
                    state->move_to_child(g, v);
                    if(use_probing) state_probe->move_to_child(g, v);

                    dej_assert(!use_probing || state->s_base_pos        == state_probe->s_base_pos);
                    dej_assert(!use_probing || state_probe->T->trace_equal());
                    dej_assert(!use_probing || state->c->cells == state_probe->c->cells);
                    dej_assert(!use_probing || state->T->get_position() == state_probe->T->get_position());
                }
            }
        };

        /**
         *  \brief Store deviations for a BFS level
         */
        class deviation_map {
        private:
            std::unordered_set<uint64_t> map;
            int computed_for_base = 0;
            int expected_for_base = 0;
            bool deviation_done = false;

            void check_finished() {
                if(computed_for_base == expected_for_base) deviation_done = true;
                dej_assert(computed_for_base <= expected_for_base);
            }

        public:
            void start(const int h_expected_for_base) {
                computed_for_base = 0;
                expected_for_base = h_expected_for_base;

                map.clear();
                deviation_done = false;
            }

            void record_deviation(uint64_t deviation) {
                map.insert(deviation);
                ++computed_for_base;
                dej_assert(computed_for_base <= expected_for_base);
                check_finished();
            }

            void record_no_deviation() {
                ++computed_for_base;
                dej_assert(computed_for_base <= expected_for_base);
                check_finished();
            }

            bool check_deviation(uint64_t deviation) {
                return !deviation_done || map.find(deviation) != map.end();
            }
        };

        /**
         * \brief IR leaf
         *
         * A stored leaf of an IR tree. The leaf can be stored in a dense manner (coloring of the leaf), or a sparse
         * manner (base of the walk that leads to this leaf).
         */
        class stored_leaf {
        public:
            enum stored_leaf_type { STORE_LAB, ///< stores coloring of the leaf
                                    STORE_BASE ///< stores only base of the leaf
            };

            stored_leaf(int* arr, int arr_sz, stored_leaf_type storetype) : store_type(storetype) {
                lab_or_base.allocate(arr_sz);
                memcpy(lab_or_base.get_array(), arr, arr_sz * sizeof(int));
                lab_or_base.set_size(arr_sz);
            }

            stored_leaf(std::vector<int>& arr, stored_leaf_type storetype) : store_type(storetype) {
                lab_or_base.allocate((int) arr.size());
                std::copy(arr.begin(), arr.end(), lab_or_base.get_array());
                lab_or_base.set_size((int) arr.size());
            }

            const int* get_lab_or_base() {
                return lab_or_base.get_array();
            }

            int get_lab_or_base_size() {
                return lab_or_base.size();
            }

            dej_nodiscard stored_leaf_type get_store_type() const {
                return store_type;
            }

        private:
            worklist lab_or_base;
            stored_leaf_type store_type;
        };

        /**
         * \brief Collection of leaves
         *
         * Can be used across multiple threads.
         *
         */
        class shared_leaves {
            std::unordered_multimap<uint64_t, stored_leaf*> leaf_store;
            std::vector<stored_leaf*> garbage_collector;

        public:
            int s_leaves          = 0;   /**< number of leaves stored */
            int h_full_save_limit = 256; /**< number of leaves which will be stored fully */

            shared_leaves() {
                leaf_store.reserve(20);
            }

            /**
             * Free up all memory.
             */
            ~shared_leaves() {
                for(auto & l : garbage_collector) {
                    delete l;
                }
            }

            /**
             * Lookup whether a leaf with the given hash already exists.
             *
             * @param hash
             * @return
             */
            stored_leaf* lookup_leaf(uint64_t hash) {
                auto find = leaf_store.find(hash);
                if(find != leaf_store.end()) {
                    return find->second;
                } else {
                    return nullptr;
                }
            }

            /**
             * Add leaf with the given hash. Does not add the leaf, if a leaf with the given hash already exists.
             *
             * @param hash
             * @param ptr
             */
            void add_leaf(uint64_t hash, coloring& c, std::vector<int>& base) {
                // check whether hash already exists
                //if(leaf_store.contains(hash)) return;
                if(leaf_store.find(hash) != leaf_store.end()) return;

                // if not, add the leaf
                const bool full_save = s_leaves < h_full_save_limit;
                auto type
                       = full_save?stored_leaf::stored_leaf_type::STORE_LAB:stored_leaf::stored_leaf_type::STORE_BASE;
                auto new_leaf
                       = full_save?new stored_leaf(c.lab,c.domain_size, type):new stored_leaf(base, type);
                leaf_store.insert(std::pair<uint64_t, stored_leaf*>(hash, new_leaf));
                garbage_collector.push_back(new_leaf);
                ++s_leaves;
            }

            /**
             * Empty this leaf container.
             */
            void clear() {
                s_leaves = 0;
                leaf_store.clear();
                // TODO clear garbage collector?
            }
        };

        /**
         * \brief A node of an IR tree
         */
        class tree_node {
            limited_save* data = nullptr;
            bool          owns_data = true;
            tree_node*    next;
            tree_node*    parent;
            bool          is_base = false;
            bool          is_pruned = false;
            uint64_t hash = 0;
        public:

            int           nodes_below  = 0;
            int           pruned_below = 0;

            tree_node(limited_save* _data, tree_node* _next, tree_node* _parent, bool ownsdata) {
                this->data = _data;
                this->next = _next;
                this->owns_data = ownsdata;
                this->parent = _parent;
            }

            tree_node* get_next() {
                return next;
            }
            tree_node* get_parent() {
                return parent;
            }
            void set_next(tree_node* new_next) {
                this->next = new_next;
            }
            limited_save* get_save() {
                return data;
            }

            void prune() {
                is_pruned = true;
            }
            dej_nodiscard bool get_prune() const {
                return is_pruned;
            }
            void add_hash(uint64_t add) {
                this->hash += add;
            }

            dej_nodiscard uint64_t get_hash() const {
                return hash;
            }

            void base() {
                is_base = true;
            }

            dej_nodiscard bool get_base() const {
                return is_base;
            }

            ~tree_node() {
                if(owns_data) delete data;
            }
        };

        typedef std::pair<ir::tree_node*, int> missing_node;

        /**
         * \brief IR tree structure
         *
         * Datastructure to explicitly store parts of an IR tree, such as a level-wise store, leaf store, as well as
         * further information used for pruning in BFS.
         *
         * Can be used across multiple threads.
         */
        class shared_tree {
            stack_t<missing_node>         missing_nodes;
            std::vector<tree_node*>              tree_data;
            std::vector<std::vector<tree_node*>> tree_data_jump_map;
            std::vector<int>        tree_level_size;
            std::vector<tree_node*> garbage_collector;
            int                     finished_up_to = 0;

            std::vector<int> current_base;

            std::vector<uint64_t> node_invariant; // TODO: move this to inprocessor?

            bool init = false;
        public:
            groups::orbit h_bfs_top_level_orbit;
            int h_bfs_automorphism_pw = 0;
            shared_leaves stored_leaves;    /**< stores leaves of the IR tree */
            deviation_map stored_deviation; /**< stores trace deviations of a BFS level*/

            explicit shared_tree(int domain_size) {
                h_bfs_top_level_orbit.initialize(domain_size);
            };

            // TODO: move this to inprocessor?
            void make_node_invariant() {
                if(finished_up_to > 1) {
                    for(int j = finished_up_to; j >= 1; --j) {
                        finish_level(j);
                        for(auto node : tree_data_jump_map[j]) {
                            const int v = node->get_save()->get_base()[0];
                            uint64_t add_hash = 1;
                            add_hash = add_to_hash(add_hash, node->get_hash());
                            add_hash = add_to_hash(add_hash, hash(j));

                            if(j == finished_up_to) {
                                int cnt = 0;
                                for (auto v_pre: node->get_save()->get_base()) {
                                    node_invariant[v_pre] += add_hash * hash(cnt);
                                    ++cnt;
                                }
                            }

                            node_invariant[node->get_save()->get_base()[j-1]] += add_hash;
                            node_invariant[v] += add_hash + 1;
                        }
                    }
                }
            }

            void introspection() {
                /*for(int j = finished_up_to; j >= 1; --j) {
                    finish_level(j);
                    for(auto node : tree_data_jump_map[j]) {
                        node->nodes_below = 0;
                    }
                }

                for(int j = finished_up_to; j >= 1; --j) {
                    for(auto node : tree_data_jump_map[j]) {
                        node->get_parent()->nodes_below += node->nodes_below + 1;
                    }
                }

                for(int j = 1; j <= finished_up_to; ++j) {
                    int nodes_below_compare = 0;
                    for(auto node : tree_data_jump_map[j]) {
                        if(node->get_base()) {
                            nodes_below_compare = node->nodes_below;
                            break;
                        }
                    }
                    int could_prune = 0;
                    for(auto node : tree_data_jump_map[j]) {
                        if(nodes_below_compare != node->nodes_below) {
                            ++could_prune;
                        }
                    }
                    std::cout << could_prune << "@" << j << std::endl;
                }*/
            }

            std::vector<uint64_t>* get_node_invariant() {
                return &node_invariant;
            }

            void initialize(std::vector<int> &base, ir::limited_save* root) {
                tree_data.resize(base.size() + 1);
                tree_level_size.resize(base.size() + 1);
                tree_data_jump_map.resize(base.size() + 1);
                add_node(0, root, nullptr, true);
                node_invariant.resize(root->get_coloring()->domain_size);
                current_base = base;
                init = true;
            }

            void clear_leaves() {
                stored_leaves.clear();
            }

            /**
             * @return How many leaves were stored during random search.
             */
            int stat_leaves() const {
                return stored_leaves.s_leaves;
            }

            bool reset(std::vector<int> &new_base, ir::limited_save* root, bool keep_old) {
                if(!init) {
                    initialize(new_base, root);
                    return false;
                }

                for(int i = 0; i < root->get_coloring()->domain_size; ++i) node_invariant[i] = 0;

                const int old_size = (int) current_base.size();
                const int new_size = (int) new_base.size();

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                if(keep_old) {
                    for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                        if (current_base[keep_until] != new_base[keep_until]) break;
                    }
                }

                if(keep_until == 0) {
                    for(auto t : garbage_collector) delete t;
                    garbage_collector.clear();
                }


                if(keep_until == new_size && new_size == old_size) return false;

                finished_up_to = std::min(keep_until, finished_up_to);

                h_bfs_top_level_orbit.reset();
                tree_data.resize(new_size + 1);
                tree_level_size.resize(new_size + 1);
                tree_data_jump_map.resize(new_size + 1);

                for (int i = keep_until+1; i < new_size+1; ++i) {
                    tree_data[i] = nullptr;
                    tree_level_size[i] = 0;
                    tree_data_jump_map[i].clear();
                }

                if(keep_until == 0) {
                    tree_level_size[0] = 0;
                    tree_data_jump_map[0].clear();
                    tree_data[0] = nullptr;
                    add_node(0, root, nullptr, true);
                }
                dej_assert(missing_nodes.empty());

                current_base = new_base;

                return true;
            }


            void queue_reserve(const int n) {
                missing_nodes.reserve(n);
            }

            void queue_missing_node(missing_node node) {
                missing_nodes.add(node);
            }

            bool queue_missing_node_empty() {
                return missing_nodes.empty();
            }

            missing_node queue_missing_node_pop() {
                return missing_nodes.pop();
            }

            void mark_first_level(markset& marks) {
                if(tree_data[1] == nullptr) return;

                marks.reset();

                tree_node * first = tree_data[1];
                tree_node * next = first;
                do {
                    marks.set(next->get_save()->get_base()[0]);
                    next = next->get_next();
                } while (next != first);
            }

            void record_invariant(int v, uint64_t inv) {
                node_invariant[v] = inv;
            }

            void record_add_invariant(int v, uint64_t inv) {
                node_invariant[v] += inv;
            }

            void add_node(int level, limited_save* data, tree_node* parent, bool is_base = false) {
                dej_assert(data != nullptr);
                if(tree_data[level] == nullptr) {
                    tree_level_size[level] = 0;
                    tree_data[level] = new tree_node(data, nullptr, parent, level != 0);
                    tree_data[level]->set_next(tree_data[level]);

                    garbage_collector.push_back( tree_data[level]);
                    if(is_base) tree_data[level]->base();
                } else {
                    tree_node* a_node    = tree_data[level];
                    tree_node* next_node = a_node->get_next();
                    auto       new_node  = new tree_node(data, next_node, parent, level != 0);
                    garbage_collector.push_back(new_node);
                    if(is_base) new_node->base();
                    a_node->set_next(new_node);
                    tree_data[level] = new_node;
                }
                ++tree_level_size[level];
            }

            void finish_level(int level) {
                if(tree_data_jump_map[level].empty()) {
                    tree_node * first = tree_data[level];
                    tree_data_jump_map[level].reserve(tree_level_size[level]);
                    tree_node * next = first;
                    do {
                        if(next->get_parent() == nullptr || !next->get_parent()->get_prune()) {
                            tree_data_jump_map[level].push_back(next);
                        }
                        next = next->get_next();
                    } while (next != first);
                    tree_level_size[level] = tree_data_jump_map[level].size();
                    dej_assert(static_cast<int>(tree_data_jump_map[level].size()) == tree_level_size[level]);
                }
            }

            ir::tree_node* pick_node_from_level(const int level, int num) {
                finish_level(level);
                num = num % tree_level_size[level];
                return tree_data_jump_map[level][num];
            }

            dej_nodiscard int get_finished_up_to() const {
                return finished_up_to;
            }

            void set_finished_up_to(const int new_finished_up_to) {
                this->finished_up_to = new_finished_up_to;
            }

            tree_node* get_level(int level) {
                return tree_data[level];
            }

            int get_current_level_size() {
                return tree_level_size[finished_up_to];
            }

            int get_current_level_tracepos() {
                return tree_data[finished_up_to]->get_save()->get_trace_position();
            }

            int get_level_size(int level) {
                return tree_level_size[level];
            }

            ~shared_tree() {
                for(auto & i : garbage_collector) delete i;
            };
        };
    }
}

#endif //DEJAVU_IR_H
