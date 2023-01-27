#ifndef SASSY_BIJECTION_H
#define SASSY_BIJECTION_H

#include "coloring.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>
#include <assert.h>

namespace sassy {
    template<class vertex_t>
    class bijection {
        bool init = false;

        void alloc(int n) {
            dealloc();
            map = new vertex_t[n];
            init = true;
        }

        void dealloc() {
            if (init) {
                delete[] map;
                init = false;
            }
        };
    public:
        vertex_t *map = nullptr;
        bool mark = false;
        int res_sz;
        int map_sz;
        bool non_uniform = false;
        bool foreign_base = false;
        bool certified = false;

        int map_vertex(int v) {
            return map[v];
        }

        void initialize_empty(int reserve) {
            alloc(reserve);
            res_sz = reserve;
            map_sz = 0;
            init = true;
        }

        void append(int v) {
            assert(init);
            assert(map_sz < res_sz);
            map[map_sz] = v;
            ++map_sz;
        }

        void copy(bijection *p) {
            assert(p->init);
            alloc(p->map_sz);
            res_sz = p->map_sz;
            init = p->init;
            mark = p->mark;
            map_sz = p->map_sz;
            non_uniform = p->non_uniform;
            foreign_base = p->foreign_base;
            certified = p->certified;
            memcpy(map, p->map, p->map_sz * sizeof(vertex_t));
        }

        void swap(bijection *p) {
            bool s_init = p->init;
            vertex_t *s_map = p->map;
            int s_map_sz = p->map_sz;

            p->init = init;
            p->map = map;
            p->map_sz = map_sz;

            init = s_init;
            map = s_map;
            map_sz = s_map_sz;
        }

        int *extract_map() {
            int *r_map = map;
            map = nullptr;
            init = false;
            return r_map;
        }

        void copy_map(int *b) {
            assert(init);
            assert(map != nullptr);
            assert(b != nullptr);
            memcpy(b, map, map_sz * sizeof(vertex_t));
        }

        void print() {
            for (int i = 0; i < map_sz; ++i)
                std::cout << map[i] << " ";
            std::cout << std::endl;
        }

        void read_from_array(const vertex_t *_map, int _map_sz) {
            alloc(_map_sz);
            res_sz = _map_sz;
            init = true;
            map_sz = _map_sz;
            for (int i = 0; i < _map_sz; ++i) {
                map[i] = _map[i];
            }
        }

        void read_from_coloring(coloring *c) {
            alloc(c->lab_sz);
            init = true;
            map_sz = c->lab_sz;
            res_sz = c->lab_sz;
            for (int i = 0; i < c->lab_sz; ++i) {
                map[i] = c->lab[i];
            }
        }

        void inverse() {
            assert(init);
            // ToDo: buffer this map
            // thread local with unique_ptr?
            vertex_t *switch_map;
            bool switch_map_init = false;

            if (!switch_map_init) {
                switch_map_init = true;
                switch_map = new vertex_t[map_sz];

            }

            vertex_t *swap = map;
            map = switch_map;
            switch_map = swap;

            for (int i = 0; i < map_sz; ++i) {
                map[switch_map[i]] = i;
            }

            if (switch_map_init) {
                switch_map_init = false;
                delete[] switch_map;
            }
        }

        void compose(bijection<vertex_t> *p) {
            assert(p->init);
            assert(init);
            for (int i = 0; i < map_sz; ++i) {
                map[i] = p->map[map[i]];
            }
        }

        bijection() {
            init = false;
        }

        ~bijection() {
            dealloc();
        }

        static void random_bijection(bijection<vertex_t> *p, int n, unsigned seed) {
            if (p->init)
                delete[] p->map;

            p->map = new vertex_t[n];
            p->init = true;
            p->map_sz = n;
            for (int i = 0; i < n; ++i) {
                p->map[i] = i;
            }
            std::default_random_engine re = std::default_random_engine(seed);
            std::shuffle(p->map, p->map + p->map_sz, re);
        }
    };
}

#endif //SASSY_BIJECTION_H
