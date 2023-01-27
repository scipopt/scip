#include <atomic>
#include <iostream>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <set>
#include <cstring>
#include <queue>
#include <memory>

#ifndef SASSY_UTILITY_H
#define SASSY_UTILITY_H

namespace sassy {

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define MASH0(i) ((unsigned long) i * (35235237 - (unsigned long) i * 5))
#define MASH1(i) ((unsigned long) i * (352355 - (unsigned long) i * 3))
#define MASH2(i) (((unsigned long) i + 1) * (423733 - ((unsigned long) i + 1)))
#define MASH3(i) (((unsigned long) i + 1) * (423233 - ((unsigned long) i + 1)))
#define MASH4(i) (((unsigned long) i + 1) * (23524361 - (unsigned long) i * 3))
#define MASH5(i) (((unsigned long) i + 1) * (23524361 - (unsigned long) i * 3))

//#define PRINT(str) {if(config->CONFIG_PRINT) {std::cout << str << std::endl;}}
#define PRINT(str) {}

// metrics used to compare strategies
    struct strategy_metrics {
        int color_refinement_cost = 0;
    };

// set specialized for quick resets
    class mark_set {
        int mark = 0;
        int *s = nullptr;
        int sz = -1;
        bool init = false;
    public:
        mark_set() {};
        mark_set(int size) {
            initialize(size);
        };
        void initialize(int size) {
            if(init)
                delete[] s;

            s = new int[size];
            sz = size;
            init = true;
            memset(s, mark, sz * sizeof(int));
            reset();
        }

        void initialize_from_array(int *arr, int size) {
            if(init)
                delete[] s;
            s = arr;
            sz = size;
            init = false;
            memset(s, mark, sz * sizeof(int));
            reset();
        }

        bool get(int pos) {
            return s[pos] == mark;
        }

        void set(int pos) {
            s[pos] = mark;
        }

        void unset(int pos) {
            s[pos] = mark - 1;
        }

        void reset() {
            if (mark == -1) {
                memset(s, mark, sz * sizeof(int));
            }
            ++mark;
        }

        ~mark_set() {
            if (init)
                delete[] s;
        }
    };

    // set specialized for quick resets with (badly implemented) bloom filter
    class bloom_mark_set {
        int mark = 0;
        int *s;
        int *bloom;
        const int bloom_sz = 8;
        int sz;
        bool init = false;
    public:
        void initialize(int size) {
            if(init) {
                delete[] s;
            }

            s     = new int[size + bloom_sz];
            bloom = s + size;
            sz = size;
            init = true;
            memset(s, mark, (sz + bloom_sz) * sizeof(int));
            reset();
        }

        void initialize_from_array(int *arr, int size) {
            if(init)
                delete[] s;
            s = arr;
            sz = size;
            init = false;
            memset(s, mark, sz * sizeof(int));
            reset();
        }

        bool get(int pos) {
            return s[pos] == mark;
        }

        void set(int pos) {
            s[pos] = mark;
        }

        void unset(int pos) {
            s[pos] = mark - 1;
        }

        void reset() {
            if (mark == -1) {
                memset(s, mark, (sz + bloom_sz) * sizeof(int));
            }
            memset(bloom, 0, bloom_sz * sizeof(int));
            ++mark;
        }

        ~bloom_mark_set() {
            if (init) {
                delete[] s;
            }
        }
    };
}

#endif //SASSY_UTILITY_H
