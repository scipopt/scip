// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_GROUPS_H
#define DEJAVU_GROUPS_H

#include "ds.h"
#include "utility.h"
#include "coloring.h"
#include "graph.h"
#include "trace.h"

namespace dejavu {

    /**
     * \brief Data structures and algorithms to deal with groups.
     *
     * Contains basic data structures to construct and deal with automorphisms, as well as a Schreier structure.
     */
    namespace groups {

        using namespace dejavu::ds;

        /**
         * Reset an automorphism given in \p automorphism with support \p support to the identity.
         *
         * @param automorphism Dense notation of the given automorphism, i.e., vertex `i` is mapped to `automorphism[i]`.
         * @param support Support of the automorphism, contains all vertices where `i != automorphism[i]`.
         */
        static void reset_automorphism(int *automorphism, worklist *support) {
            for (int i = 0; i < support->cur_pos; ++i) {
                automorphism[(*support)[i]] = (*support)[i];
            }
            support->reset();
        }

        /**
         * Create an automorphism from two discrete vertex colorings.
         *
         * @param domain_size Size of the underlying domain (i.e., number of vertices of the graph).
         * @param vertex_to_col Vertex-to-color mapping of the first coloring, i.e., vertex `i` is mapped to color
         * `vertex_to_col[i]`.
         * @param col_to_vertex Color-to-vertex mapping of the second coloring, i.e., the color i contains the vertex
         * `col_to_vertex[i]` (since the colorings are assumed to be discrete, every `i` must be a distinct color).
         * @param automorphism Dense notation for the automorphism to be written. Assumed to be the identity upon
         * calling the method.
         * @param support Support for the automorphism \p automorphism.
         */
        static void color_diff_automorphism(int domain_size, const int *vertex_to_col, const int *col_to_vertex,
                                            int *automorphism, worklist *support) {
            support->reset();
            for (int v1 = 0; v1 < domain_size; ++v1) {
                const int col = vertex_to_col[v1];
                const int v2  = col_to_vertex[col];
                if (v1 != v2) {
                    automorphism[v1] = v2;
                    support->push_back(v1);
                }
            }
        }

        /**
         * \brief Workspace for sparse automorphisms
         *
         * Enables O(1) lookup on a sparse automorphism by using an O(n) workspace.
         */
        class automorphism_workspace {
            worklist automorphism;
            worklist automorphism_supp;
            int domain_size = 0;
            bool support01 = false;

            workspace inverse_automorphism;
            bool have_inverse_automorphism = false;

            /**
             * Reconstruct inverse automorphism
             */
            void update_inverse_automorphism() {
                if(!have_inverse_automorphism) {
                    for (int i = 0; i < domain_size; ++i) {
                        const int j = automorphism[i];
                        inverse_automorphism[j] = i;
                    }
                    have_inverse_automorphism = true;
                }
            }

            /**
             * Invalidate inverse automorphism
             */
            void invalidate_inverse_automorphism() {
                have_inverse_automorphism = false;
            }

            /**
             * Updates the support using the internal dense notation.
             */
            void update_support() {
                // rewrite support
                if (!support01) {
                    automorphism_supp.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (i != automorphism[i])
                            automorphism_supp.push_back(i);
                    }
                } else {
                    automorphism_supp.reset();
                    int i;
                    for (i = 0; i < domain_size; ++i) {
                        if (i != automorphism[i]) break;
                    }
                    automorphism_supp.cur_pos = (i != domain_size);
                }
            }

        public:
            /**
             * Initializes the stored automorphism to the identity.
             *
             * @param domain_sz Size of the domain on which automorphisms operate
             */
            explicit automorphism_workspace(int domain_sz = 0) {
                automorphism.allocate(domain_sz);
                for (int i = 0; i < domain_sz; ++i) {
                    automorphism[i] = i;
                }
                automorphism_supp.allocate(domain_sz);
                inverse_automorphism.allocate(domain_sz);
                invalidate_inverse_automorphism();
                this->domain_size = domain_sz;
            }

            /**
             * Resizes the datastructure.
             *
             * @param new_domain_size allocate and initialize space for a domain of this size
             */
            void resize(int new_domain_size) {
                automorphism.resize(new_domain_size);
                for (int i = 0; i < new_domain_size; ++i) {
                    automorphism[i] = i;
                }
                automorphism_supp.resize(new_domain_size);
                automorphism_supp.reset();
                inverse_automorphism.resize(new_domain_size);
                invalidate_inverse_automorphism();
                this->domain_size = new_domain_size;
            }

            /**
             * Returns where `point` is mapped to under currently stored automorphism.
             * @param point
             * @return
             */
            inline int &operator[](int point) const {
                dej_assert(point >= 0);
                dej_assert(point < domain_size);
                return automorphism[point];
            }

            /**
             * If 0/1 support is activated, this automorphism_workspace will only track whether the support of the
             * underlying automorphism is trivial (=0) or non-trivial (>=1). Some operations are more efficient, if
             * 0/1 support is used (e.g., \a apply)
             *
             * @param new_support01 Flag of whether to use 0/1 support, or not.
             */
            void set_support01(const bool new_support01) {
                bool update_support_necessary = support01 && !new_support01;
                this->support01 = new_support01;
                if(update_support_necessary) update_support();
            }

            bool get_support01() {
                return support01;
            }

            /**
             * Create automorphism from two discrete vertex colorings.
             *
             * @param vertex_to_col Vertex-to-color mapping of the first coloring, i.e., vertex `i` is mapped to color
             * `vertex_to_col[i]`.
             * @param col_to_vertex Color-to-vertex mapping of the second coloring, i.e., the color i contains the vertex
             * `col_to_vertex[i]` (since the colorings are assumed to be discrete, every `i` must be a distinct color).
             */
            void write_color_diff(const int *vertex_to_col, const int *col_to_vertex) {
                color_diff_automorphism(domain_size, vertex_to_col, col_to_vertex, automorphism.get_array(),
                                        &automorphism_supp);
                invalidate_inverse_automorphism();
            }

            /**
             * @return Size of the support, independent of support01. Runs in O(n) if support01 is set, and in O(1)
             * otherwise.
             */
            int compute_support() {
                if(support01) {
                    int supp = 0;
                    for (int i = 0; i < domain_size; ++i) supp += (i != automorphism[i]);
                    return supp;
                } else {
                    return nsupp();
                }
            }

            /**
             * Apply another automorphism to the stored automorphism. To be more precise, `other^pwr` is applied to
             * the automorphism stored in this object.
             *
             * Closely follows the implementation in nauty / Traces.
             *
             * @param other Automorphism in sparse notation that is applied to this automorphism in.
             * @param pwr Power with which \p other is applied to this automorphism.
             */
            void apply(automorphism_workspace& other, int pwr = 1) {
                #ifndef dej_nothreadlocal
                thread_local worklist scratch_apply1;
                thread_local worklist scratch_apply2;
                thread_local markset  scratch_apply3;
                #else
                static worklist scratch_apply1;
                static worklist scratch_apply2;
                static markset  scratch_apply3;
                #endif
                scratch_apply1.allocate(domain_size);
                scratch_apply2.allocate(domain_size);
                scratch_apply3.initialize(domain_size);

                apply(scratch_apply1, scratch_apply2, scratch_apply3, other, pwr);
            }

            /**
             * Apply another automorphism to the stored automorphism. To be more precise, `other^pwr` is applied to
             * the automorphism stored in this object.
             *
             * Closely follows the implementation in nauty / Traces.
             *
             * @param scratch_apply1 Auxiliary workspace used for the operation.
             * @param scratch_apply2 Auxiliary workspace used for the operation.
             * @param scratch_apply3 Auxiliary workspace used for the operation.
             * @param other Automorphism in sparse notation that is applied to this automorphism in.
             * @param pwr Power with which \p other is applied to this automorphism.
             */
            void apply(worklist &scratch_apply1, worklist &scratch_apply2, markset &scratch_apply3,
                                        automorphism_workspace& other, int pwr = 1) {
                if(!other.support01 && other.nsupp() <= domain_size / 4) {
                    apply_sparse(scratch_apply1, scratch_apply2, scratch_apply3, other.p(), other.supp(),
                                 other.nsupp(), pwr);
                } else {
                    apply(scratch_apply1, scratch_apply2, scratch_apply3, other.p(), pwr);

                }
            }

            /**
             * Apply another automorphism to the stored automorphism. To be more precise, `other^pwr` is applied to
             * the automorphism stored in this object.
             *
             * Closely follows the implementation in nauty / Traces.
             *
             * @param scratch_apply1 Auxiliary workspace used for the operation.
             * @param scratch_apply2 Auxiliary workspace used for the operation.
             * @param scratch_apply3 Auxiliary workspace used for the operation.
             * @param other Automorphism in dense notation that is applied to this automorphism.
             * @param pwr Power with which \p other is applied to this automorphism.
             */
            void apply(worklist &scratch_apply1, worklist &scratch_apply2, markset &scratch_apply3,
                       const int *p, int pwr = 1) {
                if (pwr == 0) return;
                if (pwr <= 5) {
                    if (pwr == 1)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                    else if (pwr == 2)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                    else if (pwr == 3)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[automorphism[i]]]];
                    else if (pwr == 4)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[automorphism[i]]]]];
                    else if (pwr == 5)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[p[automorphism[i]]]]]];
                } else if (pwr <= 19) {
                    // apply other automorphism
                    for (int j = 0; j < domain_size; ++j) {
                        scratch_apply1[j] = p[p[p[j]]];
                    }
                    for (; pwr >= 6; pwr -= 6)
                        for (int j = 0; j < domain_size; ++j)
                            automorphism[j] = scratch_apply1[scratch_apply1[automorphism[j]]];

                    if (pwr == 1)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                    else if (pwr == 2)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                    else if (pwr == 3)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply1[automorphism[i]];
                    else if (pwr == 4)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[scratch_apply1[automorphism[i]]];
                    else if (pwr == 5)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[scratch_apply1[automorphism[i]]]];
                } else {
                    // 1 cycle at a time

                    scratch_apply3.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (scratch_apply3.get(i)) continue;
                        if (p[i] == i)
                            scratch_apply2[i] = i;
                        else {
                            int cyclen = 1;
                            scratch_apply1[0] = i;
                            for (int j = p[i]; j != i; j = p[j]) {
                                scratch_apply1[cyclen++] = j;
                                scratch_apply3.set(j);
                            }
                            int kk = pwr % cyclen;
                            for (int j = 0; j < cyclen; ++j) {
                                scratch_apply2[scratch_apply1[j]] = scratch_apply1[kk];
                                if (++kk == cyclen) kk = 0;
                            }
                        }
                    }
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply2[automorphism[i]];
                    scratch_apply3.reset();
                }

                // rewrite support
                update_support();
                invalidate_inverse_automorphism();
            }

            void apply_sparse(worklist &scratch_apply1, worklist &scratch_apply2, markset &scratch_apply3,
                              const int *p, const int* support, const int nsupport, int pwr = 1) {
                if (pwr == 0) return;
                // sparse version only implemented for pwr <= 6 right now
                if(pwr <= 6) {
                    // we need the inverse automorphism for the algorithm
                    update_inverse_automorphism();
                    scratch_apply1.reset();
                    scratch_apply2.reset();

                    switch(pwr) {
                        case 1:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[i];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                        break;
                        case 2:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[p[i]];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                            break;
                        case 3:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[p[p[i]]];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                            break;
                        case 4:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[p[p[p[i]]]];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                            break;
                        case 5:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[p[p[p[p[i]]]]];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                            break;
                        case 6:
                            for (int k = 0; k < nsupport; ++k) {
                                const int i = support[k];
                                dej_assert(automorphism[inverse_automorphism[i]] == i);
                                const int inv_i = inverse_automorphism[i];
                                const int p_i = p[p[p[p[p[p[i]]]]]];
                                automorphism[inv_i] = p_i;
                                scratch_apply1.push_back(inv_i);
                                scratch_apply2.push_back(p_i);
                            }
                            break;
                        default:
                            dej_assert(false); // unreachable
                            break;
                    }

                    // fix inverse automorphism
                    for (int k = 0; k < scratch_apply1.size(); ++k)
                        inverse_automorphism[scratch_apply2[k]] = scratch_apply1[k];
                    scratch_apply1.reset();
                    scratch_apply2.reset();
                } else {
                    // otherwise just use dense version...
                    apply(scratch_apply1, scratch_apply2, scratch_apply3, p, pwr);
                }

                // rewrite support
                update_support();
            }


            /**
             * Create mapping from two canonically-ordered vectors of singletons. The resulting automorphism maps
             * `singletons1[i]` to `singletons2[i]` for `i` in `pos_start, ..., pos_end`.
             *
             * May result in a mapping that is not a bijection.
             *
             * @param singletons1 first vector of singletons
             * @param singletons2 second vector of singletons
             * @param pos_start start reading the vectors at this position
             * @param pos_end stop reading the vecvtors at this position.
             */
            void write_singleton(const std::vector<int> *singletons1, const std::vector<int> *singletons2,
                                 const int pos_start, const int pos_end) {
                for (int i = pos_start; i < pos_end; ++i) {
                    const int from = (*singletons1)[i];
                    const int to   = (*singletons2)[i];
                    dej_assert(automorphism[from] == from);
                    if (from != to) {
                        automorphism_supp.push_back(from);
                        automorphism[from] = to;
                    }
                }
                invalidate_inverse_automorphism();
            }

            /**
             * Heuristic that turns the internal mapping into a bijection (in some way).
             *
             * @param scratch_set some auxiliary workspace
             */
            void cycle_completion(markset& scratch_set) {
                scratch_set.reset();
                for(int i = 0; i < automorphism_supp.size(); ++i) {
                    const int v = automorphism_supp[i];
                    if(scratch_set.get(v)) continue;
                    int v_next = automorphism[v];
                    while(v_next != v) {
                        if(scratch_set.get(v_next)) break;
                        scratch_set.set(v_next);
                        if(v_next == automorphism[v_next]) {
                            automorphism_supp.push_back(v_next);
                            automorphism[v_next] = v;
                        }
                        v_next = automorphism[v_next];
                    }
                }
                invalidate_inverse_automorphism();
            }

            /**
             * Write
             * @param from
             * @param to
             */
            void write_single_map(const int from, const int to) {
                dej_assert(automorphism[from] == from);
                if (from != to) {
                    automorphism_supp.push_back(from);
                    automorphism[from] = to;
                }
                invalidate_inverse_automorphism();
            }

            /**
             * Reset the contained automorphism back to the identity.
             */
            void reset() {
                invalidate_inverse_automorphism();
                reset_automorphism(automorphism.get_array(), &automorphism_supp);
            }

            /**
             * @return Integer array \p p describing the stored automorphism, where point v is mapped to \p p[v].
             */
            dej_nodiscard int *p() const {
                return automorphism.get_array();
            }

            /**
             * @return Integer array which contains all vertices in the support of the contained automorphism.
             */
            dej_nodiscard int* supp() const {
                return automorphism_supp.get_array();
            }

            /**
             *
             * @return Size of the support. If support01 is set, only returns 0 or 1 depending on whether support is
             * trivial.
             */
            dej_nodiscard int nsupp() const {
                return automorphism_supp.cur_pos;
            }
        };

        /**
         * \brief Orbit partition
         *
         * Keeps track of an orbit partition, and provides tools to manipulate the orbits within.
         */
        class orbit {
            int        sz = 0;
            worklist  map_arr;
            worklist  orb_sz;
        public:

            /**
             * Retrieve the orbit of the given vertex.
             *
             * @param v Vertex of the specified domain.
             * @return The orbit of \p v.
             */
            int find_orbit(const int v) {
                dej_assert(v >= 0);
                dej_assert(v < sz);
                int last_map;
                int next_map = v;
                do {
                    last_map = next_map;
                    next_map = map_arr[last_map];
                } while(next_map != last_map);
                map_arr[v] = next_map;
                return next_map;
            }

            /**
             * Returns the size of an orbit.
             *
             * @param v Vertex of the specified domain.
             * @return Size of the orbit of \p v.
             */
            int orbit_size(const int v) {
                dej_assert(v >= 0);
                dej_assert(v < sz);
                return orb_sz[find_orbit(v)];
            }

            /**
             * Every orbit has precisely one representative. This function enables to test this.
             *
             * @param v Vertex of the specified domain.
             * @return Whether \p v is the representative of the orbit of \p v.
             */
            dej_nodiscard bool represents_orbit(const int v) const {
                return v == map_arr[v];
            }

            /**
             * Combines the orbits of two given vertices.
             *
             * @param v1 The first vertex.
             * @param v2  The second vertex.
             */
            void combine_orbits(const int v1, const int v2) {
                dej_assert(v1 >= 0);
                dej_assert(v2 >= 0);
                dej_assert(v1 < sz);
                dej_assert(v2 < sz);
                if(v1 == v2) return;
                const int orbit1 = find_orbit(v1);
                const int orbit2 = find_orbit(v2);
                if(orbit1 == orbit2) return;
                if(orbit1 < orbit2) {
                    map_arr[orbit2] = orbit1;
                    orb_sz[orbit1] += orb_sz[orbit2];
                } else {
                    map_arr[orbit1] = orbit2;
                    orb_sz[orbit2] += orb_sz[orbit1];
                }
            }


            /**
             * Checks whether two given vertices are in the same orbit
             *
             * @param v1 The first vertex.
             * @param v2  The second vertex.
             * @return Whether \p v1 and \p v2 are in the same orbit.
             */
            bool are_in_same_orbit(const int v1, const int v2) {
                dej_assert(v1 >= 0);
                dej_assert(v2 >= 0);
                dej_assert(v1 < sz);
                dej_assert(v2 < sz);
                if(v1 == v2) return true;
                const int orbit1 = find_orbit(v1);
                const int orbit2 = find_orbit(v2);
                return (orbit1 == orbit2);
            }

            /**
             * Resets the datastructure to the trivial orbit partition. Takes time O(n) where n is the current domain
             * size.
             */
            void reset() {
                for(int v = 0; v < sz; ++v) {
                    map_arr[v] = v;
                    orb_sz[v]  = 1;
                }
            }

            /**
             * Applies an automorphism to the orbit structure.
             *
             * @param aut Automorphism workspace which is applied.
             */
            void add_automorphism_to_orbit(const int* p, const int nsupp, const int* supp) {
                for (int i = 0; i < nsupp; ++i) {
                    combine_orbits(p[supp[i]], supp[i]);
                }
            }

            /**
             * Applies an automorphism to the orbit structure.
             *
             * @param aut Automorphism workspace which is applied.
             */
            void add_automorphism_to_orbit(groups::automorphism_workspace& aut) {
                const int  nsupp = aut.nsupp();
                const int* supp  = aut.supp();
                const int* p     = aut.p();
                for (int i = 0; i < nsupp; ++i) {
                    combine_orbits(p[supp[i]], supp[i]);
                }
            }

            /**
             * Initializes the orbit structure with the given size.
             *
             * @param domain_size size of the underlying domain
             */
            void initialize(int domain_size) {
                if(sz == domain_size) {
                    reset();
                    return;
                }
                sz = domain_size;
                map_arr.allocate(domain_size);
                orb_sz.allocate(domain_size);
                for(int i = 0; i < domain_size; ++i) {
                    map_arr.push_back(i);
                    orb_sz.push_back(1);
                }
            }

            /**
             * Initializes orbit with the given domain size.
             *
             * @param domain_size size of the underlying domain
             */
            explicit orbit(int domain_size) {
                initialize(domain_size);
            }

            orbit() = default;
            orbit(const orbit& other)  {
                sz = other.sz;
                map_arr = other.map_arr;
                orb_sz = other.orb_sz;
            }



            bool operator==(orbit& other_orbit) {
                bool comp = (other_orbit.sz == sz) ;
                for(int i = 0; i < sz && comp; ++i)
                    comp = comp && (find_orbit(i) == other_orbit.find_orbit(i));
                return comp;
            }
        };

        /**
         * \brief Stores a link to an automorphism.
         *
         * Used to implement indiscriminate loading of dense and sparse automorphisms.
         */
        class dense_sparse_arbiter {
            int *loaded_automorphism = nullptr;
        public:
            void load(int *automorphism) {
                loaded_automorphism = automorphism;
            }

            int *p() {
                return loaded_automorphism;
            }
        };

        /**
         * \brief Stores an automorphism in a dense or sparse manner, dynamically.
         */
        class stored_automorphism {
        public:
            enum stored_automorphism_type { STORE_DENSE, ///< stored densely, in size O(\a domain_size)
                                            STORE_SPARSE ///< stored in minimal encoding size of automorphism
                                          };

        private:
            worklist data;
            int domain_size = 0;
            stored_automorphism_type store_type = STORE_SPARSE;

        public:
            /**
             * @return whether the automorphism is stored in a dense or sparse manner
             */
            stored_automorphism_type get_store_type() {
                return store_type;
            }

            /**
             * Load this stored automorphism into workspace.
             *
             * @param loader Store a pointer to permutation of automorphism in the arbiter.
             * @param space Auxiliary space that may or may not be used, depending on whether the loaded automorphism is
             *              sparse. Should only be reset once \p loader is not used anymore.
             */
            void load(dense_sparse_arbiter &loader, automorphism_workspace &space) {
                if (store_type == STORE_SPARSE) {
                    space.reset();
                    int first_of_cycle = 0;
                    for (int i = 0; i < data.size(); ++i) {
                        const int j = abs(data[i]) - 1;
                        const bool is_last = data[i] < 0;
                        dej_assert(i == data.size() - 1 ? is_last : true);
                        if (is_last) {
                            space.write_single_map(j, abs(data[first_of_cycle]) - 1);
                            first_of_cycle = i + 1;
                        } else {
                            dej_assert(i + 1 < data.size());
                            space.write_single_map(j, abs(data[i + 1]) - 1);
                        }
                    }
                    loader.load(space.p());
                } else {
                    // store_type == STORE_DENSE
                    loader.load(data.get_array());
                }
            }

            void load(automorphism_workspace &automorphism) {
                if (store_type == STORE_SPARSE) {
                    automorphism.reset();
                    int first_of_cycle = 0;
                    for (int i = 0; i < data.size(); ++i) {
                        const int j = abs(data[i]) - 1;
                        const bool is_last = data[i] < 0;
                        dej_assert(i == data.size() - 1 ? is_last : true);
                        if (is_last) {
                            automorphism.write_single_map(j, abs(data[first_of_cycle]) - 1);
                            first_of_cycle = i + 1;
                        } else {
                            dej_assert(i + 1 < data.size());
                            automorphism.write_single_map(j, abs(data[i + 1]) - 1);
                        }
                    }
                } else {
                    automorphism.reset();
                    for(int i = 0; i < domain_size; ++i) {
                        const int v_from = i;
                        const int v_to   = data.get_array()[i];
                        if(v_from != v_to) automorphism.write_single_map(v_from, v_to);
                    }
                }
            }

            /*void store(int new_domain_size, automorphism_workspace &automorphism, markset &helper) {
                domain_size = new_domain_size;
                dej_assert(data.empty());

                int support = 0;
                for (int i = 0; i < domain_size; ++i) support += (automorphism.p()[i] != i);

                // decide whether to store dense or sparse representation
                if (support < domain_size / 4) {
                    store_type = STORE_SPARSE;
                    helper.reset();

                    data.allocate(support);
                    for (int i = 0; i < domain_size; ++i) {
                        if (automorphism.p()[i] == i) continue;
                        const int j = i;
                        if (helper.get(j)) continue;
                        helper.set(j);
                        int map_j = automorphism.p()[j];
                        dej_assert(map_j != j);
                        while (!helper.get(map_j)) {
                            data.push_back(map_j + 1);
                            helper.set(map_j);
                            map_j = automorphism.p()[map_j];
                        }
                        dej_assert(map_j == j);
                        data.push_back(-(j + 1));
                    }
                    helper.reset();
                    dej_assert(data.size() == support);
                } else {
                    store_type = STORE_DENSE;
                    data.allocate(domain_size);
                    data.set_size(domain_size);
                    memcpy(data.get_array(), automorphism.p(), domain_size * sizeof(int));
                    dej_assert(data.size() == domain_size);
                }
            }*/

            /**
             * Store the given automorphism workspace. Depending on the support of the automorphism, it is either stored
             * using a dense encoding, or a sparse encoding based on cycle notation. Overall, uses O(|support|) space.
             *
             * @param new_domain_size The domain size of the automorphism to be stored.
             * @param automorphism The automorphism to be stored.
             * @param helper A `markset` which is used as auxiliary space for the algorithm.
             */
            void store(int new_domain_size, automorphism_workspace &automorphism, markset &helper) {
                domain_size = new_domain_size;
                dej_assert(data.empty());

                // automorphism may be encoded with support01, hence we compute the support if necessary
                int support = automorphism.compute_support();

                // decide whether to store dense or sparse representation
                if (support < domain_size / 4) {
                    store_type = STORE_SPARSE;
                    helper.reset();

                    // for sparse encoding, we need access to the support, so we deactivate support01 encoding
                    automorphism.set_support01(false);

                    // store cycle notation of the automorphism
                    data.allocate(support);
                    for (int i = 0; i < automorphism.nsupp(); ++i) {
                        const int j = automorphism.supp()[i];
                        if (automorphism.p()[j] == j) continue; // no need to consider trivially mapped vertices
                        if (helper.get(j)) continue; // we have already considered cycle of this vertex
                        helper.set(j); // mark that we have already considered the vertex

                        // move along the cycle of j, until we come back to j
                        int map_j = automorphism.p()[j];
                        dej_assert(map_j != j);
                        while (!helper.get(map_j)) {
                            data.push_back(map_j + 1); // offset +1 because we use negation in the encoding, see below
                            helper.set(map_j); // mark that we have already considered the vertex
                            map_j = automorphism.p()[map_j];
                        }

                        // finally we store j, the vertex we started with
                        dej_assert(map_j == j);
                        data.push_back(-(j + 1)); // `-` marks the end of the cycle
                    }
                    helper.reset();
                    dej_assert(data.size() == support);
                } else {
                    // we simply store the dense portion of the automorphism_workspace
                    store_type = STORE_DENSE;
                    data.allocate(domain_size);
                    data.set_size(domain_size);
                    memcpy(data.get_array(), automorphism.p(), domain_size * sizeof(int));
                    dej_assert(data.size() == domain_size);
                }
            }
        };

        /**
         * \brief Auxiliary workspace used for Schreier computations
         *
         * A global (thread local) state used for computations in Schreier structures. Used such that auxiliary space
         * does not have to be allocated or re-allocated for every single operation, but only once. Also, the same space
         * is used across different operations.
         */
        class schreier_workspace {
        public:
            /**
             * Initialize this workspace.
             *
             * @param domain_size Size of the underlying domain (i.e., number of vertices of the graph).
             */
            explicit schreier_workspace(int new_domain_size) : scratch_auto(new_domain_size) {
                scratch1.initialize(new_domain_size);
                scratch2.initialize(new_domain_size);
                scratch_apply1.allocate(new_domain_size);
                scratch_apply2.allocate(new_domain_size);
                scratch_apply3.initialize(new_domain_size);
            }

            dense_sparse_arbiter loader; /**< used for indiscriminate loading of dense and sparse automorphisms */

            markset scratch1;        /**< auxiliary space */
            markset scratch2;        /**< auxiliary space */
            worklist scratch_apply1; /**< auxiliary space used for `apply` operations */
            worklist scratch_apply2; /**< auxiliary space used for `apply` operations */
            markset  scratch_apply3; /**< auxiliary space used for `apply` operations */
            automorphism_workspace scratch_auto; /**< used to store a sparse automorphism*/
        };

        /**
         * \brief Stores a generating set.
         *
         */
        class generating_set {
            std::vector<stored_automorphism *> generators; /** list of generators */
            int  domain_size  = 0;
            long support_size = 0;
        public:
            int s_stored_sparse = 0; /**< how many generators are stored in a sparse manner */
            int s_stored_dense  = 0; /**< how many generators are stored in a dense manner  */

            generating_set() = default;
            generating_set(const generating_set&) = delete;
            generating_set& operator=(const generating_set&) = delete;

            /**
             * Set up this generating set.
             * @param domain_size Size of the domain of the stored generators.
             */
            void initialize(int new_domain_size) {
                this->domain_size = new_domain_size;
            }

            /**
             * Add a generator to this generating set.
             *
             * @param w The Schreier workspace.
             * @param automorphism The automorphism to be stored as a generator
             * @return Identifier of the new generator in the generating set.
             */
            int add_generator(schreier_workspace &w, automorphism_workspace &automorphism) {
                // store the automorphism
                generators.emplace_back(new stored_automorphism);
                const auto num = generators.size() - 1;
                generators[num]->store(domain_size, automorphism, w.scratch2);

                // update statistic on sparse and dense generators
                s_stored_sparse += (generators[num]->get_store_type() ==
                                    stored_automorphism::stored_automorphism_type::STORE_SPARSE);
                s_stored_dense += (generators[num]->get_store_type() ==
                                   stored_automorphism::stored_automorphism_type::STORE_DENSE);

                // update support size
                if (generators[num]->get_store_type() == stored_automorphism::stored_automorphism_type::STORE_DENSE)
                    support_size += domain_size;
                else
                    support_size += automorphism.nsupp();

                return static_cast<int>(num);
            }

            void remove_generator(size_t num) {
                dej_assert(num >= 0);
                dej_assert(num < generators.size());
                delete generators[num];
                generators[num] = nullptr;
            }

            /**
             * Remove generators from generating set which do not fix points in \p fix_points.
             *
             * @param w A Schreier workspace.
             * @param fix_points Points to be fixed by generating set.
             */
            void filter(schreier_workspace& w, std::vector<int> &fix_points) {
                for(int test_pt : fix_points) {
                    for(size_t j = 0; j < generators.size(); ++j) {
                        const auto gen = generators[j];
                        if(gen == nullptr) continue;
                        gen->load(w.loader, w.scratch_auto);
                        if(w.loader.p()[test_pt] != test_pt) remove_generator(j);
                        w.scratch_auto.reset();
                    }
                }
                compact_generators();
            }

            /**
             * Remove generators marked for deletion.
             */
            void compact_generators() {
                std::vector<stored_automorphism *> new_gens;
                for(auto gen : generators) {
                    if(gen) new_gens.push_back(gen);
                }
                generators.swap(new_gens);
            }

            /**
             * Retrieve a generator from the stored generating set.
             *
             * @param num Identifier of a generator.
             * @return A pointer to the generator which corresponds to the identifier.
             */
            dej_nodiscard stored_automorphism* get_generator(const size_t num) const {
                return generators[num];
            }

            stored_automorphism* operator[](const size_t num) const {
                return get_generator(num);
            }

            /**
             * @return number of stored generators
             */
            dej_nodiscard int size() const {
                return static_cast<int>(generators.size());
            }

            /**
             * @return memory size of stored generators
             */
            dej_nodiscard long get_support() const {
                return support_size;
            }

            void clear() {
                for(auto & generator : generators) {
                    delete generator;
                }
                generators.clear();
            }

            ~generating_set() {
                clear();
            }
        };

        /**
         * \brief A transversal in a Schreier structure.
         *
         * Can be used across multiple threads.
         */
        class shared_transversal {
            /*enum stored_transversal_type {
                STORE_DENSE, STORE_SPARSE
            };*/
            int fixed = -1;                       /**< vertex fixed by this transversal */
            int sz_upb = INT32_MAX;               /**< upper bound for size of the transversal (e.g. color class size) */
            int level = -1;
            bool finished = false;

            std::vector<int> fixed_orbit;         /**< contains vertices of orbit at this schreier level */
            std::vector<int> fixed_orbit_to_perm; /**< maps fixed_orbit[i] to generators[i] in class \ref schreier. */
            std::vector<int> fixed_orbit_to_pwr;  /**< power that should to be applied to generators[i] in class \ref schreier. */

            //stored_transversal_type store_type = STORE_SPARSE; /**< whether above structures are stored dense or sparse */

            /**
             * We load the stored orbit to a given schreier_workspace to unlock O(1) access.
             *
             * @param w The schreier_workspace to which the orbit is loaded.
             */
            void load_orbit_to_scratch(schreier_workspace &w) const {
                w.scratch1.reset();
                for (int p : fixed_orbit) {
                    w.scratch1.set(p);
                }
            }

            /**
             * Add \p vertex to the orbit of vertex \a fixed. Store in the transversal that we begin mapping \p vertex
             * to \a fixed by applying `perm^pwr` (further applications of permutations might be necessary, the
             * transversal stores a Schreier vector).
             *
             * @param vertex Vertex added to orbit of \a fixed.
             * @param perm Corresponding permutation used in Schreier vector stored in the transversal.
             * @param pwr Corresponding power used in Schreier vector stored in the transversal.
             */
            void add_to_fixed_orbit(const int vertex, const int perm, const int pwr) {
                fixed_orbit.push_back(vertex);
                fixed_orbit_to_perm.push_back(perm);
                fixed_orbit_to_pwr.push_back(pwr);
            }

            /**
             * The method performs two operations: first, it loads the generator \p gen_num of the generating set
             * \p generators. We denote the generator with \p gen.
             *
             * Second, it applies \p gen with power \p pwr (i.e., `gen^pwr`) to the permutation stored in \p automorphism.
             *
             *
             * @param w Schreier workspace used as auxiliary space for computations.
             * @param automorphism The automorphism to which the permutation is applied.
             * @param generators A generating set.
             * @param gen_num The number of the generator in \p generators, which is applied to \p automorphism.
             * @param pwr A power that is applied to the generator first.
             */
            static void apply_perm(schreier_workspace &w, automorphism_workspace &automorphism,
                                   generating_set &generators, const int gen_num, const int pwr) {
                // load perm into workspace
                auto generator = generators[gen_num];

                // apply generator
                dej_assert(pwr >= 0);
                if (pwr > 0) { // use generator
                    generator->load(w.loader, w.scratch_auto);
                    // multiply
                    if(generator->get_store_type() == stored_automorphism::STORE_DENSE) {
                        automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), pwr);
                    } else {
                        automorphism.apply_sparse(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(),
                                                  w.scratch_auto.supp(), w.scratch_auto.nsupp(), pwr);
                    }
                    w.scratch_auto.reset();
                }
            }

        public:

            /**
             * @return Size of the transversal.
             */
            dej_nodiscard int size() const {
                return (int) fixed_orbit.size();
            }

            void set_size_upper_bound(const int new_sz_upb) {
                sz_upb = new_sz_upb;
            }

            dej_nodiscard int get_size_upper_bound() const {
                return sz_upb;
            }
            /**
             * Check whether a point \p p is contained in transversal.
             *
             * @param p Point to check.
             * @return Position of point \p in transversal, or -1 if not contained.
             */
            dej_nodiscard int find_point(const int p) const {
                for (size_t i = 0; i < fixed_orbit.size(); ++i) {
                    if (p == fixed_orbit[i]) {
                        return (int) i;
                    }
                }
                return -1;
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in this transversal
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection) {
                load_orbit_to_scratch(w);
                int back_swap = (int) selection.size() - 1;
                int front_pt;
                for (front_pt = 0; front_pt <= back_swap;) {
                    if (!w.scratch1.get(selection[front_pt])) {
                        ++front_pt;
                    } else {
                        selection[front_pt] = selection[back_swap];
                        --back_swap;
                    }
                }
                selection.resize(front_pt);
                w.scratch1.reset();
            }

            /**
             * @return Whether size of this transversal matches the given upper bound.
             */
            dej_nodiscard bool is_finished() const {
                return finished;
            }

            /**
             * @return Point fixed by this transversal.
             */
            dej_nodiscard int get_fixed_point() const {
                return fixed;
            }

            dej_nodiscard const std::vector<int>& get_fixed_orbit() {
                return fixed_orbit;
            }

            dej_nodiscard const std::vector<int>& get_generators() {
                return fixed_orbit_to_perm;
            }

            /**
             * Initialize this transversal.
             *
             * @param fixed_vertex Vertex fixed by the transversal.
             * @param level Position of the transversal in the base of Schreier structure.
             * @param sz_upb Upper bound for the size of transversal (e.g., color class size in combinatorial base).
             */
            void initialize(const int fixed_vertex, const int new_level, const int new_sz_upb) {
                dej_assert(fixed_vertex >= 0);
                dej_assert(new_level >= 0);
                fixed = fixed_vertex;
                this->level  = new_level;
                this->sz_upb = new_sz_upb;
                add_to_fixed_orbit(fixed_vertex, -1, 0);
            }

            /**
             * Extend this transversal using the given \p automorphism.
             *
             * @param generators Generating set of this Schreier structure.
             * @param automorphism Extend the transversal using this automorphism.
             * @return whether transversal was extended
             */
            bool extend_with_automorphism(schreier_workspace &w, generating_set &generators,
                                          automorphism_workspace &automorphism) {
                if (finished) return false; // Already finished transversal? Nothing to do, then!

                // load orbit of this transversal to our workspace, such that we have random access to points in O(1)
                load_orbit_to_scratch(w);
                bool changed = false; /*< we changed this transversal? */
                int  gen_num = -1;    /*< the generator we added to extend this transversal, -1 means no generator */

                // watch out, we may enlarge fixed_orbit in the loop below
                for (int i = 0; i < static_cast<int>(fixed_orbit.size()); ++i) {
                    const int p = fixed_orbit[i];
                    int j = automorphism.p()[p];
                    dej_assert(j >= 0);
                    if (j == p || w.scratch1.get(j)) continue;

                    int pwr = 0; // power we save in Schreier vector
                    for (int jj = j; !w.scratch1.get(jj); jj = automorphism.p()[jj]) {
                        ++pwr;
                    }

                    while (!w.scratch1.get(j)) {
                        // we change this transversal
                        changed = true;

                        // add generator to generating set (once)
                        if (gen_num == -1) gen_num = generators.add_generator(w, automorphism);

                        // add entry to transversal
                        add_to_fixed_orbit(j, gen_num, pwr);
                        w.scratch1.set(j);

                        // we check out the entire cycle of j now
                        j = automorphism.p()[j];
                        --pwr;
                    }
                }

                // We reached upper bound for the size of this transversal? Then mark transversal as "finished"!
                if (sz_upb == (int) fixed_orbit.size() && !finished) {
                    finished = true;
                }

                dej_assert(sz_upb >= (int) fixed_orbit.size());

                // reset our workspace
                w.scratch1.reset();
                return changed;
            }

            /**
             * Fix \a fixed in \p automorphism.
             *
             * @param generators The underlying generating set.
             * @param automorphism Automorphism where \a fixed will be fixed.
             * @return whether \p automorphism is now the identity
             */
            bool fix_automorphism(schreier_workspace &w, generating_set &generators,
                                  automorphism_workspace &automorphism) const {
                int fixed_map = automorphism.p()[fixed]; // Where is fixed mapped to?

                // as long as `fixed` is not yet fixed, we apply automorphisms from the transversal as prescribed by
                // the Schreier vector
                while (fixed != fixed_map) {
                    const int pos = find_point(fixed_map); // Where are we storing the information for `fixed_map`?
                    dej_assert(pos >= 0);
                    const int perm = fixed_orbit_to_perm[pos]; // generator to apply for `fixed_map`
                    const int pwr  = fixed_orbit_to_pwr[pos];  // power to use for `fixed_map`
                    apply_perm(w, automorphism, generators, perm, pwr);
                    fixed_map = automorphism.p()[fixed]; // Fixed now? Or we need to go again?
                }
                dej_assert(automorphism.p()[fixed] == fixed);
                return automorphism.nsupp() == 0;
            }
        };

        /**
         * \brief Schreier structure
         *
         * Enables sifting of automorphisms into a Schreier structure with given base. Intended for internal use in
         * dejavu.
         *
         */
        class random_schreier_internal {
            int  domain_size    = -1;
            int  finished_up_to = -1;

            generating_set generators;
            std::vector<shared_transversal> transversals;
            std::vector<int> stabilized_generators;

            bool init = false;

        public:
            int s_consecutive_success = 0;  /**< track repeated sifts for probabilistic abort criterion */
            int h_error_bound         = 10; /**< determines error probability                           */
            big_number s_grp_sz; /**< size of the automorphism group computed */

            /**
             * @return Number of stored generators using a sparse data structure.
             */
            dej_nodiscard int s_sparsegen() const {
                return generators.s_stored_sparse;
            }

            /**
             * @return Number of stored generators using a dense data structure.
             */
            dej_nodiscard int s_densegen() const {
                return generators.s_stored_dense;
            }

            /**
             * Set up this Schreier structure using the given base. The base is then fixed and can not be adjusted
             * later on.
             *
             * @param base the base
             * @param base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             */
            void initialize(const int new_domain_size, std::vector<int> &base, std::vector<int> &base_sizes,
                            int stop = INT32_MAX) {
                stop = std::min(static_cast<int>(base.size()), stop);
                dej_assert(static_cast<int>(base.size()) >= stop);
                domain_size = new_domain_size;
                dej_assert(this->domain_size > 0);
                generators.initialize(domain_size);
                generators.clear();
                transversals.reserve(stop);
                transversals.clear();
                finished_up_to = -1;

                for (int i = 0; i < stop && i < static_cast<int>(base.size()); ++i) {
                    transversals.emplace_back();
                    transversals[i].initialize(base[i], i, base_sizes[i]);
                }
                init = true;
            }

            void set_base(schreier_workspace &w, automorphism_workspace& automorphism, random_source& rng,
                          std::vector<int> &new_base, int err = 10, bool resift_generators = false) {
                dej_assert(init);
                const int old_size = static_cast<int>(transversals.size());
                const int new_size = static_cast<int>(new_base.size());

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                    if (transversals[keep_until].get_fixed_point() != new_base[keep_until]) break;
                }

                if (keep_until == new_size) return;

                finished_up_to = -1;
                transversals.resize(new_size);


                for (int i = 0; i < keep_until; ++i) { transversals[i].set_size_upper_bound(INT32_MAX); }
                for (int i = keep_until; i < new_size; ++i) {
                    transversals[i] = shared_transversal();
                    dej_assert(new_base[i] >= 0);
                    transversals[i].initialize(new_base[i], i, INT32_MAX);
                }

                std::vector<int> stabilized_generators_copy;
                stabilized_generators_copy.swap(stabilized_generators);
                stabilized_generators.clear();

                if(resift_generators) {
                    for (int gen_id : stabilized_generators_copy) sift_generator(w, automorphism, gen_id, true);
                }
                sift_random(w, automorphism, rng, err);
            }

            void sift_random(schreier_workspace &w, automorphism_workspace& automorphism,
                             random_source& rng, int err) {
                int fail = err;
                bool any_changed = false;
                while(fail >= 0) {
                    const bool sift_changed = sift_random(w, automorphism, rng, true);
                    any_changed = sift_changed || any_changed;
                    fail -= !sift_changed;
                }
            }

            /**
             * Reset up this Schreier structure with a new base.
             *
             * @param new_base the new base
             * @param new_base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             * @param keep_old If true, attempt to keep parts of the base that is already stored.
             */
            bool reset(int new_domain_size, schreier_workspace& w, std::vector<int> &new_base,
                       std::vector<int> &new_base_sizes, const int stop, bool keep_old, bool remove_generators,
                       std::vector<int> &global_fixed_points) {
                const int old_size = static_cast<int>(transversals.size());
                if(!init || new_domain_size != domain_size) {
                    initialize(new_domain_size, new_base, new_base_sizes, stop);
                    return false;
                }
                const int new_size = stop;

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                if(remove_generators) generators.clear();
                if(keep_old) {
                    for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                        if (transversals[keep_until].get_fixed_point() != new_base[keep_until]) break;
                    }
                } else {
                    //generators.clear();
                    generators.filter(w, global_fixed_points);
                }

                if(keep_until == new_size && new_size == old_size) return false;

                finished_up_to = -1;

                transversals.resize(new_size);
                //transversals.set_size(new_size);


                for(int i = 0; i < keep_until; ++i) {
                    transversals[i].set_size_upper_bound(new_base_sizes[i]);
                }
                for (int i = keep_until; i < stop; ++i) {
                    //if(i < old_size) delete transversals[i];
                    transversals[i] = shared_transversal();
                    dej_assert(new_base[i] >= 0);
                    transversals[i].initialize(new_base[i], i, new_base_sizes[i]);
                }

                stabilized_generators.clear();

                return true;
            }

            /**
             * Returns a vertex to individualize for each color of \p root_coloring that matches in size a corresponding
             * transversal.
             *
             * @param save_to_individualize Vector in which vertices deemed save to individualize are pushed.
             * @param root_coloring The coloring with which the stored transversals are compared.
             */
            void determine_potential_individualization(std::vector<std::pair<int, int>>* save_to_individualize,
                                                       coloring* root_coloring) {
                for (int i = base_size()-1; i >= 0; --i) {
                    const int corresponding_root_color_sz =
                            root_coloring->ptn[root_coloring->vertex_to_col[transversals[i].get_fixed_point()]] + 1;
                    if(transversals[i].size() >= corresponding_root_color_sz && corresponding_root_color_sz > 1) {
                        save_to_individualize->emplace_back(transversals[i].get_fixed_point(), corresponding_root_color_sz);
                    }
                }
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            dej_nodiscard int base_point(int pos) const {
                return transversals[pos].get_fixed_point();
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            dej_nodiscard int base_size() const {
                return static_cast<int>(transversals.size());
            }

            /**
             * Checks whether a vertex \p v is contained in the transversal at position \p s_base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the transversal at position \p s_base_pos.
             */
            bool is_in_fixed_orbit(const int base_pos, const int v) {
                if (base_pos >= base_size()) return false;
                if(v < 0) return false;
                dej_assert(base_pos >= 0);
                dej_assert(base_pos < base_size());
                const int search = transversals[base_pos].find_point(v);
                return search != -1;
            }

            const std::vector<int>& get_fixed_orbit(const int base_pos) {
                dej_assert(base_pos >= 0);
                dej_assert(base_pos < base_size());
                return transversals[base_pos].get_fixed_orbit();
            }

            const std::vector<int>& get_stabilizer_generators(const int base_pos) {
                dej_assert(base_pos >= 0);
                dej_assert(base_pos < base_size());
                return transversals[base_pos].get_generators();
            }

            const std::vector<int>& get_stabilized_generators() {
                return stabilized_generators;
            }

            dej_nodiscard int get_fixed_orbit_size(const int base_pos) {
                if (base_pos >= static_cast<int>(transversals.size())) return 0;
                return transversals[base_pos].size();
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in transversal at position
             * \p s_base_pos in Schreier structure.
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             * @param base_pos Position in base.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection, int base_pos) {
                transversals[base_pos].reduce_to_unfinished(w, selection);
            }

            /**
             * Checks whether the transversal at position \p s_base_pos matches its size upper bound.
             *
             * @param base_pos Position in base.
             * @return Bool that indicates whether the transversal at position \p s_base_pos matches its size upper bound.
             */
            dej_nodiscard bool is_finished(const int base_pos) const {
                return transversals[base_pos].is_finished();
            }

            /**
             * @return memory size of stored generators
             */
            dej_nodiscard long get_support() const {
                return generators.get_support();
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace &w, automorphism_workspace &automorphism, bool uniform = false,
                      bool keep_at_end = false) {
                bool changed = false; /*< keeps track of whether we changed the Schreier structure while sifting */

                automorphism.set_support01(true); // we don't need to track full support
                for (int level = 0; level < static_cast<int>(transversals.size()); ++level) { // sift level-by-level
                    // first, we try to extend the transversal using the new automorphism
                    changed = transversals[level].extend_with_automorphism(w, generators, automorphism)
                              || changed;

                    // secondly, we fix the point of this transversal in automorphism
                    const bool identity = transversals[level].fix_automorphism(w, generators, automorphism);
                    if (identity) break; // if automorphism is the identity now, no need to sift further
                }

                if(automorphism.nsupp() != 0 && keep_at_end) {
                    const int gen_id = generators.add_generator(w, automorphism);
                    stabilized_generators.push_back(gen_id);
                }

                // keep track of how far this Schreier structure is finished (matches given upper bounds)
                for (int level = finished_up_to + 1; level < static_cast<int>(transversals.size()); ++level) {
                    if (finished_up_to == level - 1 && transversals[level].is_finished()) {
                        ++finished_up_to;
                    } else {
                        break;
                    }
                }

                    // reset the automorphism_workspace
                automorphism.set_support01(false);
                automorphism.reset();

                if(uniform) record_sift_result(changed); // uniform automorphisms count towards abort criterion

                return changed;
            }

            /**
             * Generate a (semi-)random element from the generators.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Random element is saved in this automorphism_workspace.
             * @param rn_generator Random number generator used for the generation.
             */
            void generate_random(schreier_workspace& w, automorphism_workspace& automorphism,
                                 random_source& rng) {
                automorphism.reset();

                const int random_mult = static_cast<int>(rng() % 7); // 7
                const int num_mult = 1 + (random_mult);
                for(int i = 0; i < num_mult; ++i) {
                    // load generator
                    const int next_gen_num = static_cast<int>(rng() % generators.size());
                    dej_assert(next_gen_num >= 0);
                    dej_assert(next_gen_num < generators.size());
                    auto next_gen = generators[next_gen_num];
                    dej_assert(next_gen != nullptr);
                    next_gen->load(w.loader, w.scratch_auto);

                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), 1);
                    w.scratch_auto.reset();
                }
            }

            void load_generator(automorphism_workspace& automorphism, int generator) {
                auto next_gen = generators[generator];
                next_gen->load(automorphism);
            }

            /**
             * Sift a (semi-)randomly generated element into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism An automorphism_workspace used to store the randomly generated element.
             * @return Whether the generated automorphism was added to the Schreier structure or not.
             */
            bool sift_random(schreier_workspace &w, automorphism_workspace& automorphism,
                             random_source& rng, bool keep_at_end = false) {
                if(generators.size() <= 0) return false;
                automorphism.reset();
                automorphism.set_support01(true);
                generate_random(w, automorphism, rng);
                const bool added_generator = sift(w, automorphism, false, keep_at_end);
                return added_generator;
            }

            bool sift_generator(schreier_workspace &w, automorphism_workspace& automorphism, int generator,
                                bool keep_at_end = false) {
                if(generators.size() <= 0) return false;
                automorphism.reset();
                automorphism.set_support01(true);
                load_generator(automorphism, generator);
                const bool added_generator = sift(w, automorphism, false, keep_at_end);
                return added_generator;
            }

            /**
             * Records a sift result for the probabilistic abort criterion.
             *
             * @param changed Whether the sift was successful or not.
             */
            void record_sift_result(const bool changed) {
                if(!changed) {
                    ++s_consecutive_success;
                } else {
                    if(s_consecutive_success > 0) {
                        ++h_error_bound;
                        s_consecutive_success = 0;
                    }
                }
            }

            /**
             * Reset the probabilistic abort criterion.
             */
            void reset_probabilistic_criterion() {
                s_consecutive_success = 0;
            }

            /**
             * @return Whether the probabilistic abort criterion allows termination or not.
             */
            dej_nodiscard bool probabilistic_abort_criterion() const {
                return (s_consecutive_success > h_error_bound);
            }

            /**
             * @return Whether the deterministic abort criterion allows termination or not.
             */
            dej_nodiscard bool deterministic_abort_criterion() const {
                return (finished_up_to_level() + 1 == base_size());
            }

            /**
             * @return Level up to which Schreier structure is guaranteed to be complete according to given upper bounds.
             *         -1 indicates no level has been finished.
             */
            dej_nodiscard int finished_up_to_level() const {
                return finished_up_to;
            }

            /**
             * @return Size of group described by this Schreier structure.
             */
            void compute_group_size() {
                s_grp_sz.mantissa = 1.0;
                s_grp_sz.exponent = 0;
                // multiply the sizes of the individual levels in the Schreier table
                for (auto & transversal : transversals) s_grp_sz.multiply(transversal.size());
            }
        };

        /**
         * \brief API for the dejavu Schreier structure
         *
         * Enables sifting of automorphisms into a Schreier structure with given base. The Schreier structure does not
         * compute proper random elements of the group, hence no guarantees regarding the error bound are given. The
         * group described by the Schreier structure is always guaranteed to be a subgroup of the actual group generated
         * by the contained elements.
         *
         */
        class random_schreier {
        private:
            int h_domain_size    = -1;          /**< size of the underlying domain */
            random_schreier_internal schreier;  /**< Schreier structure */

            automorphism_workspace ws_auto;     /**< a workspace to store automorphisms */
            schreier_workspace     ws_schreier; /**< a workspace for Schreier-Sims computations */
            random_source rng;                  /**< source of randomness */

            markset auxiliary_set;              /**< auxiliary workspace used by some methods */

            int h_error_bound = 10; /**< higher value reduces likelihood of error (AKA missing generators) */
            bool init = false;      /**< whether the structure has been initialized */

        public:
            /**
             * @return Number of stored generators using a sparse data structure.
             */
            dej_nodiscard int get_number_of_generators() const {
                return schreier.s_sparsegen() + schreier.s_densegen();
            }

            /**
             * Loads the `i-th` generator into the given automorphism_workspace. Valid values for `i` are between `0`
             * and `get_number_of_generators()-1`.
             *
             * @param i will load the `i-th` generator
             * @param automorphism workspace to load the generator into
             */
            void get_generator(int i, automorphism_workspace& automorphism) {
                schreier.load_generator(automorphism, i);
            }

            /**
             * Initializes a random Schreier structure with the given parameters.
             *
             * @param domain_size the size of the domain (e.g., number of vertices of the graph)
             * @param error higher values reduce the likelihood of missing generators (default is 10)
             * @param true_random whether to use true random number generation
             * @param seed the seed if pseudo random numbers are being used
             */
            explicit random_schreier(int domain_size, int error = 10, bool true_random = false, int seed = 0) :
                                     ws_auto(domain_size), ws_schreier(domain_size), rng(true_random, seed),
                                     auxiliary_set(domain_size) {
                h_error_bound = error;
                h_domain_size = domain_size;
            }

            /**
             * Sets the base of the Schreier structure.
             *
             * @param new_base the base
             */
            void set_base(std::vector<int> &new_base, bool resift_generators = false) {
                if(!init) {
                    std::vector<int> base_sizes;
                    base_sizes.reserve(new_base.size());
                    for(int i = 0; i < static_cast<int>(new_base.size()); ++i) base_sizes.push_back(INT32_MAX);
                    schreier.initialize(h_domain_size, new_base, base_sizes);
                    init = true;
                } else schreier.set_base(ws_schreier, ws_auto, rng, new_base, h_error_bound, resift_generators);
                sift_random();
            }

            /**
             * Generate random elements to complete the Schreier structure.
             */
            void sift_random() {
                schreier.sift_random(ws_schreier, ws_auto, rng, h_error_bound);
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            dej_nodiscard int base_size() const {
                return schreier.base_size();
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            dej_nodiscard int get_fixed_point(int pos) const {
                return schreier.base_point(pos);
            }

            /**
             * Checks whether a vertex \p v is contained in the transversal at position \p s_base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the transversal at position \p s_base_pos.
             */
            bool is_in_fixed_orbit(const int base_pos, const int v) {
                return schreier.is_in_fixed_orbit(base_pos, v);
            }

            /**
             * Returns the orbit size of the fixed point at \p base_pos.
             * @param base_pos position of base to look at
             * @return the orbit size
             */
            dej_nodiscard int get_fixed_orbit_size(const int base_pos) {
                return schreier.get_fixed_orbit_size(base_pos);
            }

            /**
             * Returns the orbit of the fixed point at \p base_pos.
             *
             * @param base_pos position of base to look at
             * @return the fixed orbit
             */
            const std::vector<int>& get_fixed_orbit(const int base_pos) {
                return schreier.get_fixed_orbit(base_pos);
            }

            /**
             * Computes the entire orbit partition at \p base_pos.
             *
             * @param base_pos position of base to look at
             * @param orbit_partition orbits will be read into this orbit structure
             */
            void get_stabilizer_orbit(int base_pos, orbit& orbit_partition) {
                auxiliary_set.initialize(get_number_of_generators());
                auxiliary_set.reset();
                orbit_partition.reset();

                // TODO this can be done much more efficiently...
                for(int j = base_pos; j < schreier.base_size(); ++j) {
                    const std::vector<int> &generators = schreier.get_stabilizer_generators(j);
                    for (auto i: generators) {
                        if (i < 0) continue;
                        dej_assert(i < get_number_of_generators());
                        if (auxiliary_set.get(i)) continue;
                        auxiliary_set.set(i);
                        schreier.load_generator(ws_auto, i);
                        orbit_partition.add_automorphism_to_orbit(ws_auto);
                    }
                }

                for(auto i : schreier.get_stabilized_generators()) {
                    schreier.load_generator(ws_auto, i);
                    orbit_partition.add_automorphism_to_orbit(ws_auto);
                }
            }

            /**
             * Computes a list of generators for the pointwise stabilizer fixing the first \p base_pos points of the
             * current base. The returned vector contains a list of generator ID's as stored in the Schreier structure.
             * The generators can then be loaded using `get_generator`.
             *
             * @param base_pos position of base to consider
             * @returns a list of generator numbers
             */
            std::vector<int> get_stabilizer_generators(int base_pos) {
                auxiliary_set.initialize(get_number_of_generators());
                auxiliary_set.reset();
                std::vector<int> all_generators;

                // TODO this can be done much more efficiently...
                for(int j = base_pos; j < schreier.base_size(); ++j) {
                    const std::vector<int> &generators = schreier.get_stabilizer_generators(j);
                    for (auto i: generators) {
                        if (i < 0) continue;
                        dej_assert(i < get_number_of_generators());
                        if (auxiliary_set.get(i)) continue;
                        auxiliary_set.set(i);
                        all_generators.push_back(i);
                    }
                }

                for(auto i : schreier.get_stabilized_generators()) {
                    if (auxiliary_set.get(i)) continue;
                    auxiliary_set.set(i);
                    all_generators.push_back(i);
                }

                return all_generators;
            }

            /**
             * Sift automorphism into the Schreier structure. If the automorphism is not yet represented in the Schreier
             * structure, the structure is extended using the automorphism. I.e., the routine will add the given
             * automorphism to the structure, if necessary.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(automorphism_workspace &automorphism, bool known_in_group = false) {
                return sift(h_domain_size, automorphism.p(), automorphism.nsupp(), automorphism.supp(), known_in_group);
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param p complete bijection of the automorphism
             * @param nsupp number of points in the support of automorphism
             * @param supp support of automorphism
             * @param known_in_group whether we know that the given automorphism is already in the group or not
             * @return whether a transversal changed
             */
            bool sift(int, const int *p, int nsupp, const int *supp, bool known_in_group = false) {
                ws_auto.reset();
                for(int i = 0; i < h_domain_size; ++i) dej_assert(ws_auto.p()[i] == i);

                for(int i = 0; i < nsupp; ++i) {
                    const int v_from = supp[i];
                    const int v_to   = p[v_from];
                    if(v_from != v_to) ws_auto.write_single_map(v_from, v_to);
                }

                const bool result = schreier.sift(ws_schreier, ws_auto, false, !known_in_group);
                ws_auto.reset();
                return result;
            }

            /**
             * @return Order of group described by this Schreier structure. Note that if the base is incomplete, the
             * group order will not match the actual order the group.
             */
            big_number group_size() {
                sift_random();
                schreier.compute_group_size();
                return schreier.s_grp_sz;
            }
        };

        /**
         * \brief Compresses vertex set to smaller window
         */
        class domain_compressor {
            std::vector<int> map_to_small;
            std::vector<int> map_to_big;
            bool do_compress = false;
            int s_vertices_active = 0;
            int domain_size = 0;

            void reset() {
                s_compression_ratio = 1.0;
                do_compress = false;
            }
        public:
            double s_compression_ratio = 1.0;

            /**
             * Determines which colors of given vertex coloring are not needed for the given base
             * @param c the vertex coloring
             * @param base the base
             * @param stop place to stop reading the \a base
             */
            void determine_compression_map(coloring& c, std::vector<int>& base, int stop) {
                do_compress = true;
                domain_size = c.domain_size;

                markset color_was_added(c.domain_size);
                markset test_vertices_active(c.domain_size);

                color_was_added.reset();
                test_vertices_active.reset();

                s_vertices_active = 0;

                for(int i = 0; i < stop; ++i) {
                    const int base_vertex = base[i];
                    const int col = c.vertex_to_col[base_vertex];
                    if(!color_was_added.get(col)) {
                        color_was_added.set(col);
                        const int col_sz = c.ptn[col] + 1;
                        s_vertices_active += col_sz;
                        for(int j = 0; j < col_sz; ++j) {
                            const int add_v = c.lab[col + j];
                            test_vertices_active.set(add_v);
                        }
                    }
                }

                map_to_small.resize(c.domain_size);
                map_to_big.resize(s_vertices_active);
                int next_active_v_maps_to = 0;
                for(int v = 0; v < c.domain_size; ++v) {
                    if(test_vertices_active.get(v)) {
                        map_to_small[v] = next_active_v_maps_to;
                        map_to_big[next_active_v_maps_to] = v;
                        ++next_active_v_maps_to;
                    } else {
                        map_to_small[v] = -1;
                    }
                }
                s_compression_ratio = 1.0;
                if(domain_size > 0) s_compression_ratio = 1.0 * s_vertices_active / domain_size;
            }

            /**
             * @return how large the compressed domain is
             */
            dej_nodiscard int compressed_domain_size() const {
                return s_vertices_active;
            }

            /**
             * @return how large the original domain is
             */
            dej_nodiscard int decompressed_domain_size() const {
                return domain_size;
            }

            /**
             * Mapping from original domain to compressed domain.
             * @param v vertex of the original domain
             * @return \a v in the compressed domain, or `-1` if not contained in compressed domain
             */
            int compress(int v) {
                if(!do_compress) return v;
                return map_to_small[v];
            }

            /**
             * Mapping from compressed domain to original domain
             * @param v vertex of compressed domain
             * @return \a v in the original domain
             */
            int decompress(int v) {
                if(!do_compress) return v;
                return map_to_big[v];
            }
        };

        /**
         * \brief Compressed Schreier structure.
         *
         * Internally, stores a `random_schreier_internal` structure in a compressed form using `domain_compressor`.
         */
        class compressed_schreier {
            random_schreier_internal internal_schreier;
            domain_compressor*       compressor;
            automorphism_workspace   compressed_automorphism;
            std::vector<int> original_base;
            std::vector<int> original_base_sizes;

        public:
            double h_min_compression_ratio = 0.4; /*< only use compression if at least this ratio can be achieved */
            double s_compression_ratio     = 1.0;

            /**
             * @return Number of stored generators using a sparse data structure.
             */
            dej_nodiscard int s_sparsegen() const {
                return internal_schreier.s_sparsegen();
            }

            /**
             * @return Number of stored generators using a dense data structure.
             */
            dej_nodiscard int s_densegen() const {
                return internal_schreier.s_densegen();
            }

            /**
             * Reset up this Schreier structure with a new base.
             *
             * @param new_base the new base
             * @param new_base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             * @param keep_old If true, attempt to keep parts of the base that is already stored.
             */
            bool reset(domain_compressor* new_compressor, int new_domain_size, schreier_workspace& w,
                       std::vector<int> &new_base, std::vector<int> &new_base_sizes, const int stop, bool keep_old,
                       std::vector<int> &global_fixed_points) {
                bool remove_generators = (compressor != nullptr);

                compressor = new_compressor;
                // do not use compression at all if ration too bad (to save on translation cost, and such that keep_old
                // works for difficult graphs)
                if(compressor != nullptr &&
                    compressor->s_compression_ratio > h_min_compression_ratio) compressor = nullptr;

                remove_generators = remove_generators || (compressor != nullptr);

                if(compressor != nullptr) {
                    // need more management if we are using compression
                    original_base       = new_base;
                    original_base_sizes = new_base_sizes;
                    compressed_automorphism.resize(compressor->compressed_domain_size());
                    keep_old = false;
                    std::vector<int> new_basec;
                    new_basec.reserve(new_base.size());
                    for (int i = 0; i < static_cast<int>(new_base.size()); ++i)
                        new_basec.push_back(compressor->compress(new_base[i]));

                    s_compression_ratio = compressor->s_compression_ratio;
                    return internal_schreier.reset(new_compressor->compressed_domain_size(), w, new_basec,
                                                    new_base_sizes, stop, keep_old, remove_generators,
                                                    global_fixed_points);
                } else {
                    // we are not using compression, so simply use the Schreier structure normally
                    s_compression_ratio = 1.0;
                    return internal_schreier.reset(new_domain_size, w, new_base,
                                                    new_base_sizes, stop, keep_old, remove_generators,
                                                    global_fixed_points);
                }
            }

            /**
             * Returns a vertex to individualize for each color of \p root_coloring that matches in size a corresponding
             * transversal.
             *
             * @param save_to_individualize Vector in which vertices deemed save to individualize are pushed.
             * @param root_coloring The coloring with which the stored transversals are compared.
             */
            void determine_potential_individualization(std::vector<std::pair<int, int>>* save_to_individualize,
                                                       coloring* root_coloring) {
                if(compressor == nullptr) {
                    internal_schreier.determine_potential_individualization(save_to_individualize, root_coloring);
                    return;
                }

                for (int i = base_size()-1; i >= 0; --i) {
                    const int corresponding_root_color_sz = root_coloring->ptn[root_coloring->vertex_to_col[original_base[i]]] + 1;
                    if(internal_schreier.get_fixed_orbit_size(i) >= corresponding_root_color_sz && corresponding_root_color_sz > 1) {
                        save_to_individualize->emplace_back(original_base[i], corresponding_root_color_sz);
                    }
                }
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            dej_nodiscard int base_point(int pos) const {
                if(compressor != nullptr) {
                    return original_base[pos];
                } else {
                    return internal_schreier.base_point(pos);
                }
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            dej_nodiscard int base_size() const {
                return internal_schreier.base_size();
            }

            /**
             * Checks whether a vertex \p v is contained in the transversal at position \p s_base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the transversal at position \p s_base_pos.
             */
            bool is_in_base_orbit(const int base_pos, const int v) {
                int v_translate = compressor != nullptr? compressor->compress(v) : v;
                return internal_schreier.is_in_fixed_orbit(base_pos, v_translate);
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in transversal at position
             * \p s_base_pos in Schreier structure.
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             * @param base_pos Position in base.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection, int base_pos) {
                if(compressor != nullptr) {
                    for(int i = 0; i < static_cast<int>(selection.size()); ++i) {
                        dej_assert(compressor->compress(selection[i]) >= 0);
                        selection[i] = compressor->compress(selection[i]);
                    }
                }
                internal_schreier.reduce_to_unfinished(w, selection, base_pos);
                if(compressor != nullptr) {
                    for(int i = 0; i < static_cast<int>(selection.size()); ++i) {
                        selection[i] = compressor->decompress(selection[i]);
                        dej_assert(selection[i] >= 0);
                    }
                }

            }

            void compress_automorphism(automorphism_workspace &automorphism, automorphism_workspace &automorphism_compress) {
                automorphism_compress.reset();

                const int support = automorphism.nsupp();
                for(int i = 0; i < support; ++i) {
                    const int vsupport = automorphism.supp()[i];
                    const int vsupportc = compressor->compress(vsupport);
                    const int vmapsto = automorphism.p()[vsupport];
                    const int vmapstoc = compressor->compress(vmapsto);
                    if(vsupportc >= 0) {
                        dej_assert(vmapstoc >= 0);
                        automorphism_compress.write_single_map(vsupportc, vmapstoc);
                    }
                }
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace &w, automorphism_workspace &automorphism, bool uniform = false) {
                if(compressor != nullptr) {
                    compress_automorphism(automorphism, compressed_automorphism);
                    return internal_schreier.sift(w, compressed_automorphism, uniform);
                } else {
                    return internal_schreier.sift(w, automorphism, uniform);
                }
            }

            /**
             * Sift a (semi-)randomly generated element into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism An automorphism_workspace used to store the randomly generated element.
             * @return Whether the generated automorphism was added to the Schreier structure or not.
             */
            bool sift_random(schreier_workspace &w, automorphism_workspace& automorphism,
                             random_source& rng) {
                if(compressor != nullptr) {
                    return internal_schreier.sift_random(w, compressed_automorphism, rng);
                } else {
                    return internal_schreier.sift_random(w, automorphism, rng);
                }
            }

            /**
             * Reset the probabilistic abort criterion.
             */
            void reset_probabilistic_criterion() {
                internal_schreier.reset_probabilistic_criterion();
            }

            /**
             * @return memory size of stored generators
             */
            dej_nodiscard long get_support() const {
                return internal_schreier.get_support();
            }

            /**
             * @return Whether the probabilistic abort criterion allows termination or not.
             */
            dej_nodiscard bool probabilistic_abort_criterion() const {
                return internal_schreier.probabilistic_abort_criterion();
            }

            /**
             * @return Whether the deterministic abort criterion allows termination or not.
             */
            dej_nodiscard bool deterministic_abort_criterion() const {
                return internal_schreier.deterministic_abort_criterion();
            }

            /**
             * @return Whether any abort criterion allows termination or not.
             */
            dej_nodiscard bool any_abort_criterion() const {
                return internal_schreier.probabilistic_abort_criterion() ||
                       internal_schreier.deterministic_abort_criterion();
            }

            void set_error_bound(int error_bound) {
                internal_schreier.h_error_bound = error_bound;
            }

            dej_nodiscard int get_consecutive_success() const {
                return internal_schreier.s_consecutive_success;
            }

            dej_nodiscard big_number get_group_size() const {
                return internal_schreier.s_grp_sz;
            }

            /**
             * @return Level up to which Schreier structure is guaranteed to be complete according to given upper bounds.
             *         -1 indicates no level has been finished.
             */
            dej_nodiscard int finished_up_to_level() const {
                return internal_schreier.finished_up_to_level();
            }

            /**
             * @return Size of group described by this Schreier structure.
             */
            void compute_group_size() {
                internal_schreier.compute_group_size();
            }
        };
    }
}

#endif //DEJAVU_GROUPS_H
