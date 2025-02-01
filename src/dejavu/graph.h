// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef SASSY_GRAPH_BUILDER_H
#define SASSY_GRAPH_BUILDER_H

#include<unordered_map>
#include <fstream>
#include "utility.h"
#include "coloring.h"

namespace dejavu {
    /**
     * \brief Internal graph data structure
     *
     * Graph data structure as used internally by the dejavu solver.
     */
    class sgraph {
        struct vertexComparator {
            explicit vertexComparator(const sgraph &g) : graph(g) {}

            const sgraph &graph;

            bool operator()(const int &v1, const int &v2) const {
                return graph.d[v1] < graph.d[v2];
            }
        };

        /*struct vertexComparatorColor {
            vertexComparatorColor(const sgraph &g, const int *vtocol) : graph(g), vertex_to_col(vtocol) {}

            const sgraph &graph;
            const int *vertex_to_col;

            bool operator()(const int &v1, const int &v2) const {
                //return (g.d[v1] < g.d[v2]) || ((g.d[v1] == g.d[v2]) && (vertex_to_col[v1] < vertex_to_col[v2]));
                return (vertex_to_col[v1] < vertex_to_col[v2]);
            }
        };*/

    public:
        bool initialized = false;
        int *v = nullptr; /**< maps vertices to positions in edge array `e`*/
        int *d = nullptr; /**< maps vertices to their degree */
        int *e = nullptr; /**< vertices adjacent to given vertex w are stored in `e[v[w]] ... e[v[w] + d[w]]`*/

        int v_size = 0; /**< number of vertices of the graph */
        int e_size = 0; /**< number of (directed) edges of the graph, i.e., the size of the array `e` */

        bool dense = false; /**< flag whether the graph is considered dense, set internally by the solver */

        /**
         * Allocate arrays according to the given sizes.
         * @param nv number of vertices
         * @param ne number of (directed) edges
         */
        void initialize(int nv, int ne) {
            initialized = true;
            v = new int[nv];
            d = new int[nv];
            e = new int[ne];
        }

        // initialize a vertex coloring of this sgraph, that already partitions according to the degrees of vertices
        void initialize_coloring(ds::coloring *c, int *vertex_to_col) {
            c->initialize(this->v_size);
            std::memset(c->ptn, 1, sizeof(int) * v_size);

            if (this->v_size < c->domain_size) {
                c->domain_size = this->v_size;
            }

            if (v_size == 0)
                return;

            int cells = 0;
            int last_new_cell = 0;

            if (vertex_to_col == nullptr) {
                for (int i = 0; i < v_size; i++) {
                    c->lab[i] = i;
                }
                std::sort(c->lab, c->lab + c->domain_size, vertexComparator(*this));
                for (int i = 0; i < c->domain_size; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->domain_size) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    dej_assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
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
                            dej_assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = it->second;
                            dej_assert(found_col < this->v_size);
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
                            dej_assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = colors[vertex_to_col[i]];
                            dej_assert(found_col < this->v_size);
                            vertex_to_col[i] = found_col;
                            ++colsize[found_col];
                        }
                    }
                }

                int increment = 0;
                for (int &i: colsize) {
                    const int col_sz = i;
                    i += increment;
                    increment += col_sz;
                }
                dej_assert(increment == v_size);

                for (int i = 0; i < v_size; i++) {
                    const int v_col = vertex_to_col[i];
                    --colsize[v_col];
                    const int v_lab_pos = colsize[v_col];
                    c->lab[v_lab_pos] = i;
                }

                /*for(int i = 0; i < v_size; i++) {
                    c->lab[i] = i;
                }
                std::sort(c->lab, c->lab + c->lab_sz, vertexComparatorColor(*this, vertex_to_col));*/
                for (int i = 0; i < c->domain_size; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->domain_size) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    //assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
                    //if(this->d[c->lab[i]] < this->d[c->lab[i + 1]]  || (this->d[c->lab[i]] == this->d[c->lab[i + 1]]
                    //&& (vertex_to_col[c->lab[i]] < vertex_to_col[c->lab[i + 1]]))) {
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

        void sanity_check() {
#if defined(DEJDEBUG) &&  !defined(NDEBUG)
            for(int i = 0; i < v_size; ++i) {
                dej_assert(d[i]>0?v[i] < e_size:true);
                dej_assert(d[i]>0?v[i] >= 0:true);
                dej_assert(d[i] >= 0);
                dej_assert(d[i] < v_size);
            }
            for(int i = 0; i < e_size; ++i) {
                dej_assert(e[i] < v_size);
                dej_assert(e[i] >= 0);
            }

            // multiedge test
            dejavu::ds::markset multiedge_test;
            multiedge_test.initialize(v_size);
            for(int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for(int j = 0; j < d[i]; ++j) {
                    const int neigh = e[v[i] + j];
                    dej_assert(!multiedge_test.get(neigh));
                    multiedge_test.set(neigh);
                }
            }

            // fwd - bwd test
            multiedge_test.initialize(v_size);
            for(int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for(int j = 0; j < d[i]; ++j) {
                    const int neigh = e[v[i] + j];
                    dej_assert(neigh >= 0);
                    dej_assert(neigh < v_size);
                    bool found = false;
                    for(int k = 0; k < d[neigh]; ++k) {
                        const int neigh_neigh = e[v[neigh] + k];
                        if(neigh_neigh == i) {
                            found = true;
                            break;
                        }
                    }
                    dej_assert(found);
                }
            }
#endif
        }

        void copy_graph(const sgraph *g) {
            if (initialized) {
                delete[] v;
                delete[] d;
                delete[] e;
            }
            initialize(g->v_size, g->e_size);

            memcpy(v, g->v, g->v_size * sizeof(int));
            memcpy(d, g->d, g->v_size * sizeof(int));
            memcpy(e, g->e, g->e_size * sizeof(int));
            v_size = g->v_size;
            e_size = g->e_size;
        }

        void sort_edgelist() const {
            for (int i = 0; i < v_size; ++i) {
                const int estart = v[i];
                const int eend = estart + d[i];
                std::sort(e + estart, e + eend);
            }
        }

        ~sgraph() {
            if (initialized) {
                delete[] v;
                delete[] d;
                delete[] e;
            }
        }
    };

    /**
     * \brief Graph with static number of vertices and edges
     *
     * Graph format based on the internal format of dejavu (`sgraph`), but with additional sanity checks and easy access
     * to the construction. Essentially, this class provides a more convenient interface to construct `sgraph`s.
     *
     * The graph must first be initialized (either using the respective constructor or using initialize_graph). For the
     * initialization, the final number of vertices or edges must be given. The number of vertices or edges can not be
     * changed. Then, using add_vertex and add_edge, the precise number of defined vertices and edges must be added.
     * The `add_vertex(color, deg)` function requires a color and a degree. Both can not be changed later.
     *
     * The `add_edge(v1, v2)` function adds an undirected edge from `v1` to `v2`. It is always required that `v1 < v2`
     * holds, to prevent the accidental addition of hyper-edges.
     *
     * After the graph is built, the internal graph (sgraph) can be accessed either by the user, or the provided
     * functions. Once the graph construction is finished, the internal sgraph can be changed arbitrarily.
     */
    class static_graph {
    private:
        sgraph   g;
        int*     c        = nullptr;
        int*     edge_cnt = nullptr;
        unsigned int num_vertices_defined  = 0;
        unsigned int num_edges_defined     = 0;
        unsigned int num_deg_edges_defined = 0;
        bool initialized = false;
        bool finalized   = false;

        void finalize() {
            if(!finalized) {
                if (!initialized)
                    throw std::logic_error("uninitialized graph");
                if (num_vertices_defined != (unsigned int) g.v_size)
                    throw std::logic_error("did not add the number of vertices requested by constructor");
                if (num_edges_defined != (unsigned int) g.e_size) {
                    std::cout << num_edges_defined << " vs. " << g.e_size << std::endl;
                    throw std::logic_error("did not add the number of edges requested by constructor");
                }
                sanity_check();
                finalized = true;
            }
        }

        void alloc(const int nv) {
            c        = new int[nv];
            edge_cnt = new int[nv];
            initialized = true;
        }

        void dealloc() {
            if(initialized && c != nullptr)
                delete[] c;
            if(initialized && edge_cnt != nullptr)
                delete[] edge_cnt;
            initialized = false;
        }
    public:

        /**
         * @param nv number of vertices
         * @param ne number of undirected edges
         */
        static_graph(const int nv, const int ne) {
            initialize_graph(nv, ne);
        };

        static_graph() = default;

        static_graph(const static_graph& other)  {
            copy(&other);
        }

        ~static_graph() {
            dealloc();
        }

        /**
         * @param nv number of vertices
         * @param ne number of undirected edges
         */
        void initialize_graph(const int nv, const int ne) {
            if(nv < 0) throw std::out_of_range("number of vertices must be positive");
            if(ne < 0) throw std::out_of_range("number of edges must be positive");
            if(initialized || finalized)
                throw std::logic_error("can not initialize a graph that is already initialized");
            alloc(nv);
            g.initialize((int) nv, (int) (2*ne));
            g.v_size = (int) nv;
            g.e_size = (int) (2*ne);
            for(int i = 0; i < nv; ++i)
                edge_cnt[i] = 0;
        };

        /**
         * Add a vertex to the graph with given color and degree. Returns the index of the vertex, which is simply the
         * number of vertices previously added.
         *
         * @param color the vertex color
         * @param deg the number of outgoing edges of the vertex
         * @return index of the newly added vertex (= number of previously added vertices)
         */
        unsigned int add_vertex(const int color, const int deg) {
            if(!initialized)
                throw std::logic_error("uninitialized graph");
            if(finalized)
                throw std::logic_error("can not change finalized graph");
            const unsigned int vertex = num_vertices_defined;
            ++num_vertices_defined;
            if(num_vertices_defined > (unsigned int) g.v_size)
                throw std::out_of_range("vertices out-of-range, define more vertices initially");
            c[vertex]   = color;
            g.d[vertex] = deg;
            g.v[vertex] = (int) num_deg_edges_defined;
            num_deg_edges_defined += deg;
            if(num_deg_edges_defined > (unsigned int) g.e_size)
                throw std::out_of_range("edges out-of-range, define more edges initially");
            return vertex;
        };

        /**
         * Add an undirected edge between two given vertices.
         *
         * @param v1 first vertex incident to the edge
         * @param v2 second vertex incident to the edge
         */
        void add_edge(const unsigned int v1, const unsigned int v2) {
            if(!initialized)
                throw std::logic_error("uninitialized graph");
            if(finalized)
                throw std::logic_error("can not change finalized graph");
            if(v1 > v2 || v1 == v2)
                throw std::invalid_argument("invalid edge: v1 < v2 must hold");
            if(v1 >= num_vertices_defined)
                throw std::out_of_range("v1 is not a defined vertex, use add_vertex to add vertices");
            if(v2 >= num_vertices_defined)
                throw std::out_of_range("v2 is not a defined vertex, use add_vertex to add vertices");
            if(static_cast<int>(num_edges_defined + 2) > g.e_size)
                throw std::out_of_range("too many edges");
            if(v1 > INT32_MAX)
                throw std::out_of_range("v1 too large, must be < INT32_MAX");
            if(v2 > INT32_MAX)
                throw std::out_of_range("v2 too large, must be < INT32_MAX");
            ++edge_cnt[v1];
            const int edg_cnt1 = edge_cnt[v1];
            if(edg_cnt1 > g.d[v1])
                throw std::out_of_range("too many edges incident to v1");
            g.e[g.v[v1] + edg_cnt1 - 1] = static_cast<int>(v2);

            ++edge_cnt[v2];
            const int edg_cnt2 = edge_cnt[v2];
            if(edg_cnt2 > g.d[v2])
                throw std::out_of_range("too many edges incident to v2");
            g.e[g.v[v2] + edg_cnt2 - 1] = static_cast<int>(v1);

            num_edges_defined += 2;
        };

        void sanity_check() {
            g.sanity_check();
        }

        /**
         * Export the graph to a file. Outputs the graph in DIMACS format.
         *
         * @param filename the file
         */
        void dump_dimacs(const std::string& filename) {
            finalize();
            std::ofstream dumpfile;
            dumpfile.open (filename, std::ios::out);

            dumpfile << "p edge " << g.v_size << " " << g.e_size/2 << std::endl;

            for(int i = 0; i < g.v_size; ++i) {
                dumpfile << "n " << i+1 << " " << c[i] << std::endl;
            }

            for(int i = 0; i < g.v_size; ++i) {
                for(int j = g.v[i]; j < g.v[i]+g.d[i]; ++j) {
                    const int neighbour = g.e[j];
                    if(neighbour < i) {
                        dumpfile << "e " << neighbour+1 << " " << i+1 << std::endl;
                    }
                }
            }
        }

        /**
         * @return internal graph representation as used by dejavu
         */
        dejavu::sgraph* get_sgraph() {
            finalize();
            return &g;
        };

        /**
         * @return internal vertex coloring representation as used by dejavu
         */
        int* get_coloring() {
            finalize();
            return c;
        }

        void copy(const static_graph* other) {
            finalized   = false;
            dealloc();
            if(!other->initialized) return;
            g.copy_graph(&other->g);
            alloc(other->g.v_size);
            memcpy(c, other->c, other->g.v_size  * sizeof(int));
            memcpy(edge_cnt, other->edge_cnt, other->g.v_size  * sizeof(int));
            num_vertices_defined  = other->num_vertices_defined;
            num_edges_defined     = other->num_edges_defined;
            num_deg_edges_defined = other->num_deg_edges_defined;
            initialized           = other->initialized;
            finalized             = other->finalized;
        }

        static_graph& operator=(const static_graph& other) {
            if(&other == this) return *this;
            copy(&other);
            return *this;
        }
    };
}

static inline bool parse_dimacs(const std::string& filename, dejavu::sgraph* g, int** colmap, bool silent=true,
                         int seed_permute=0) {
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    std::ifstream infile(filename);

    std::vector<int> reshuffle;

    std::vector<std::vector<int>> incidence_list;
    std::set<int> degrees;
    std::set<int> colors;
    std::string line;
    std::string nv_str, ne_str;
    std::string nv1_string, nv2_string;
    int nv1, nv2;
    size_t i;
    int nv = 0;
    int ne = 0;
    int line_number = 0;

    bool fail = false;
    bool initialized = false;
    while (std::getline(infile, line)) {
        ++line_number;
        const char m = line[0];
        int average_d;
        switch (m) {
            case 'p':
                for(i = 7; i < line.size() && line[i] != ' '; ++i) nv_str += line[i];

                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) ne_str += line[i];

                try {
                    nv = std::stoi(nv_str);
                    ne = std::stoi(ne_str);
                } catch(...) {
                    fail = true;
                    break;
                }

                if(line[2] != 'e' || line[3] != 'd'  || line[4] != 'g' ||  line[5] != 'e' || nv < 0 || ne < 0) {
                    fail = true;
                    break;
                }

                initialized = true;
                average_d = (ne / nv) + 3;
                g->initialize(nv, ne * 2);

                reshuffle.reserve(nv);
                for(int j = 0; j < nv; ++j) reshuffle.push_back(j + 1);
                if(seed_permute != 0) {
                    std::mt19937 eng(seed_permute);
                    std::shuffle(std::begin(reshuffle), std::end(reshuffle), eng);
                }

                incidence_list.reserve(nv);
                for(int j = 0; j < nv; ++j) {
                    incidence_list.emplace_back();
                    incidence_list[incidence_list.size() - 1].reserve(average_d);
                }
                break;
            case 'e':
                if(!initialized) {
                    fail = true;
                    break;
                }

                nv1_string = "";
                nv2_string = "";
                for(i = 2; i < line.size() && line[i] != ' '; ++i) nv1_string += line[i];

                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) nv2_string += line[i];

                nv1 = reshuffle[std::stoi(nv1_string)-1];
                nv2 = reshuffle[std::stoi(nv2_string)-1];

                if(nv1 < 0 || nv2 < 0 || nv1 > nv || nv2 > nv) {
                    fail = true;
                    break;
                }

                incidence_list[nv1 - 1].push_back(nv2 - 1);
                incidence_list[nv2 - 1].push_back(nv1 - 1);
                break;
            case 'n':
                if(!initialized) {
                    fail = true;
                    break;
                }

                if(*colmap == nullptr) *colmap = (int *) calloc(nv, sizeof(int));

                nv1_string = "";
                nv2_string = "";
                for(i = 2; i < line.size() && line[i] != ' '; ++i) nv1_string += line[i];

                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) nv2_string += line[i];

                nv1 = reshuffle[std::stoi(nv1_string)-1];
                nv2 = std::stoi(nv2_string);

                if(nv1 < 0 || nv1 > nv) {
                    fail = true;
                    break;
                }

                (*colmap)[nv1 - 1] = nv2;
                break;
            default:
                break;
        }
        if(fail) break;
    }

    if(fail) {
        std::cout << "parsing failed in line " << line_number << std::endl;
        std::cout << "> \'" << line << "\'" << std::endl;
        return false;
    }

    if(!initialized) {
        std::cout << "file does not contain a graph " << std::endl;
        return false;
    }

    int epos = 0;
    int vpos = 0;

    for(auto & incident : incidence_list) {
        g->v[vpos] = epos;
        g->d[vpos] = (int) incident.size();
        vpos += 1;
        for(int j : incident) {
            g->e[epos] = j;
            epos += 1;
        }
    }

    g->v_size = nv;
    g->e_size = 2 * ne;

    dej_assert(nv == g->v_size);
    dej_assert(2 * ne == g->e_size);
    const double parse_time = (double) (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
    if(!silent) std::cout << std::setprecision(2) << "parse_time=" << parse_time / 1000000.0 << "ms";

    return true;
}

#endif //SASSY_GRAPH_BUILDER_H
