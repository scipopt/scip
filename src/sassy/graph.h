#ifndef SASSY_GRAPH_BUILDER_H
#define SASSY_GRAPH_BUILDER_H

#include "sgraph.h"

namespace sassy {
    // Graph format based on the internal format of sassy, but adding sanity checks and easy access to the construction.
    // The graph must first be initialized (either using the respective constructor or using initialize_graph).
    // For the initialization, the final number of vertices or edges must be given. The number of vertices or edges
    // can not be changed. Then, using add_vertex and add_edge, the precise number of defined vertices and edges must be
    // added. The add_vertex(color, deg) function requests a color and a degree. Both can not be changed later.
    // The add_edge(v1, v2) function adds an undirected edge from v1 to v2. It is always required that v1 < v2 holds, to
    // prevent the accidental addition of hyper-edges.
    //
    // After the graph was built, the internal sassy graph (sgraph) can be accessed either by the user, or the provided
    // functions. Once the graph construction is finished, the internal sgraph can be changed arbitrarily.
    class static_graph {
    private:
        sgraph   g;
        int*     c = nullptr;
        int*     edge_cnt;
        unsigned int num_vertices_defined  = 0;
        unsigned int num_edges_defined     = 0;
        unsigned int num_deg_edges_defined = 0;
        bool initialized;
        bool finalized = false;

    private:
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
    public:
        static_graph(const unsigned int nv, const unsigned int ne) {
            g.initialize(nv, 2*ne);
            g.v_size = nv;
            g.d_size = nv;
            g.e_size = 2*ne;
            c = new int[nv];
            edge_cnt = new int[nv];
            for(unsigned int i = 0; i < nv; ++i)
                edge_cnt[i] = 0;
            initialized = true;
        };

        static_graph() {
            initialized = false;
        }

        ~static_graph() {
            if(initialized && c != nullptr)
                delete[] c;
            if(initialized && edge_cnt != nullptr)
                delete[] edge_cnt;
        }

        void initialize_graph(const unsigned int nv, const unsigned int ne) {
            if(initialized)
                throw std::logic_error("can not initialize a graph that is already initialized");
            if(finalized)
                throw std::logic_error("can not change finalized graph");
            initialized = true;
            g.initialize(nv, 2*ne);
            g.v_size = nv;
            g.d_size = nv;
            g.e_size = 2*ne;
            c = new int[nv];
            edge_cnt = new int[nv];
            for(unsigned int i = 0; i < nv; ++i)
                edge_cnt[i] = 0;
        };

        int add_vertex(const int color, const int deg) {
            if(!initialized)
                throw std::logic_error("uninitialized graph");
            if(finalized)
                throw std::logic_error("can not change finalized graph");
            const int vertex = num_vertices_defined;
            ++num_vertices_defined;
            if(num_vertices_defined > (unsigned int) g.v_size)
                throw std::out_of_range("vertices out-of-range, define more vertices initially");
            c[vertex]   = color;
            g.d[vertex] = deg;
            g.v[vertex] = num_deg_edges_defined;
            num_deg_edges_defined += deg;
            return vertex;
        };

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
            ++edge_cnt[v1];
            const int edg_cnt1 = edge_cnt[v1];
            if(edg_cnt1 > g.d[v1])
                throw std::out_of_range("too many edges incident to v1");
            g.e[g.v[v1] + edg_cnt1 - 1] = v2;

            ++edge_cnt[v2];
            const int edg_cnt2 = edge_cnt[v2];
            if(edg_cnt2 > g.d[v2])
                throw std::out_of_range("too many edges incident to v2");
            g.e[g.v[v2] + edg_cnt2 - 1] = v1;

            num_edges_defined += 2;
        };

        void sanity_check() {
            g.sanity_check();
        }

        void dump_dimacs(std::string filename) {
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

        void shuffle() {
            std::cout << "graph shuffle ";
            finalize();
            std::vector<int> vec;
            std::vector<int> old_v;
            std::vector<int> old_d;
            std::vector<int> old_e;
            std::vector<int> old_col;
            for(int i = 0; i < g.v_size; ++i) {
                vec.push_back(i);
                old_v.push_back(g.v[i]);
                old_d.push_back(g.d[i]);
                old_col.push_back(c[i]);
            }

            std::random_device rd;
            auto rng = std::default_random_engine { rd() };
            std::shuffle(std::begin(vec), std::end(vec), rng);

            for(int i = 0; i < g.v_size; ++i) {
                std::cout << vec[i] << " ";
            }
            std::cout << std::endl;

            for(int i = 0; i < g.v_size; ++i) {
                g.v[vec[i]] = old_v[i];
                g.d[vec[i]] = old_d[i];
                c[vec[i]]   = old_col[i];
            }

            for(int i = 0; i < g.e_size; ++i) {
                g.e[i] = vec[g.e[i]];
            }
        }

        sgraph* get_sgraph() {
            finalize();
            return &g;
        };

        int*    get_coloring() {
            finalize();
            return c;
        };
    };
}

#endif //SASSY_GRAPH_BUILDER_H
