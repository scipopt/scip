// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DS_H
#define DEJAVU_DS_H

#include <list>
#include <iostream>
#include <cstring>
#include <functional>
#include <algorithm>
#include <cassert>
#include "utility.h"
#include "coloring.h"

namespace dejavu {

    /**
     * \brief General-purpose datastructures.
     *
     */
    namespace ds {

        /**
         * Sorting.
         *
         * @tparam T Template parameter for the type of array elements.
         * @param arr Array of elements of type \p T.
         * @param sz Length of the array \p arr.
         */
        template<class T>
        void inline sort_t(T *arr, int sz) {
            std::sort(arr, arr + sz);
        }

        /**
         * \brief Stack datastructure
         *
         * @tparam T Type of elements stored on the stack.
         */
        template<class T>
        class stack_t {
        private:
            std::vector<T> queue;
        public:
            /**
             * Add an element \p t to the queue.
             * @param t Element to be added.
             */
            void add(T &t) {
                queue.emplace_back(t);
            }

            /**
             * Reserve space for at least \n elements.
             * @param n Space to be reserved.
             */
            void reserve(int n) {
                queue.reserve(n);
            }

            /**
             * @return Whether the queue is empty.
             */
            bool empty() {
                return queue.empty();
            }

            /**
             * Pop an element from the queue and return it.
             *
             * @return The element popped from the queue.
             */
            T pop() {
                auto element = queue.back();
                queue.pop_back();
                return element;
            }
        };

        /**
         * \brief Fixed-size array, uninitialized
         *
         * An array of fixed size, with some further convenience functions.
         *
         * @tparam T The type of array elements.
         */
        template<class T>
        class worklist_t {
        private:
            /**
             * Allocate an array of size \p size.
             * @param size Space to allocate.
             */
            void alloc(const int size) {
                dealloc();
                arr = new T[size];
                arr_sz = size;
            }

            void dealloc() {
                //if(arr) free(arr);
                if(arr) delete[] arr;
                arr = nullptr;
            }

        public:
            /**
             * Default constructor, does not allocate any memory.
             */
            worklist_t() = default;

            worklist_t(const worklist_t<T>& other) {
                copy(&other);
            }

            worklist_t(const worklist_t<T>&& other) {
                swap(other);
            }

            worklist_t<T>& operator=(const worklist_t<T>& other) {
                if(other.arr == arr) return *this;
                copy(&other);
                return *this;
            }

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit worklist_t(int size) {
                dej_assert(size >= 0);
                allocate(size);
            }

            void swap(worklist_t<T>& other) {
                auto swap_arr = arr;
                auto swap_arr_sz = arr_sz;
                auto swap_cur_pos = cur_pos;
                arr = other.arr;
                arr_sz = other.arr_sz;
                cur_pos = other.cur_pos;
                other.arr = swap_arr;
                other.arr_sz = swap_arr_sz;
                other.cur_pos = swap_cur_pos;
            }

            void copy(const worklist_t<T>* other) {
                alloc(other->arr_sz);
                for(int i = 0; i < other->arr_sz; ++i) {
                    arr[i] = other->arr[i];
                }
                arr_sz  = other->arr_sz;
                cur_pos = other->cur_pos;
            }

            /**
             * Allocates the internal array with size \p size. The allocated memory is not
             * initialized. Initializes the internal position \a cur_pos of the array at 0.
             *
             * @param size Size to allocate.
             */
            void allocate(int size) {
                dej_assert(size >= 0);
                alloc(size);
                cur_pos = 0;
            }

            /**
             * Push back an element at position \a cur_pos. Increments the internal position \a cur_pos.
             *
             * @param value Element to push back.
             */
            inline void push_back(T value) {
                dej_assert(cur_pos >= 0 && cur_pos < arr_sz);
                arr[cur_pos] = value;
                cur_pos += 1;
            }

            /**
             * Pop an element at position \a cur_pos. Decreases the internal position \a cur_pos. There is no safeguard
             * in place if `cur_pos <= 0`.
             *
             * @return Element popped at \a cur_pos.
             *
             * \sa The function empty() tests whether `cur_pos == 0`.
             */
            inline T pop_back() {
                dej_assert(cur_pos > 0);
                return arr[--cur_pos];
            }

            /**
             * @return Element at \a cur_pos.
             */
            inline T *last() const {
                return &arr[cur_pos - 1];
            }

            /**
             * @return Whether `cur_pos == 0`.
             */
            dej_nodiscard bool empty() const {
                return cur_pos == 0;
            }

            /**
             * Sets \a cur_pos to \p size.
             *
             * @param size Value to set \a cur_pos to.
             */
            void set_size(const int size) {
                cur_pos = size;
            }

            /**
             * @return The current position \a cur_pos.
             */
            dej_nodiscard int size() const {
                return cur_pos;
            }

            /**
             * Sets \a cur_pos to `0`.
             */
            inline void reset() {
                cur_pos = 0;
            }

            int max() {
                int max_val = arr[0];
                for(int i = 1; i < cur_pos; ++i) {
                    const int val = arr[i];
                    if(val > max_val) max_val = val;
                }
                return max_val;
            }

            /**
             * Resizes the internal array to `size`. Copies the contents of the old array into the new one. If the new
             * size is larger than the old one, the new space is only allocated and not initialized.
             *
             * @param size New size to allocate the array to.
             */
            void resize(const int size) {
                dej_assert(size >= 0);
                if (arr && size <= arr_sz) return;
                T *old_arr = nullptr;
                int old_arr_sz = arr_sz;
                if (arr) old_arr = arr;
                arr = nullptr;
                alloc(size);
                if (old_arr != nullptr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(T));
                    delete[] old_arr;
                }
            }

            /**
             * Deallocates the internal array.
             */
            ~worklist_t() {
                dealloc();
            }

            /**
             * Sort the internal array up to position \a cur_pos.
             */
            void sort() {
                sort_t<T>(arr, cur_pos);
            }

            /**
             * @return A pointer to the internal memory.
             */
            inline T *get_array() const {
                return arr;
            }

            /**
             * Access element \p index in the internal array \p arr.
             *
             * @param index Index of the internal array.
             * @return The element `arr[index]`.
             */
            inline T &operator[](int index) const {
                dej_assert(index >= 0);
                dej_assert(index < arr_sz);
                return arr[index];
            }

            void sort_after_map(T *map) {
                struct comparator_map {
                    T *map;

                    explicit comparator_map(T *m) {
                        this->map = m;
                    }

                    bool operator()(const T &a, const T &b) {
                        return map[a] < map[b];
                    }
                };
                comparator_map c = comparator_map(map);
                std::sort(arr, arr + cur_pos, c);
            }

            int cur_pos = 0; /**< current position */
        private:
            int arr_sz = 0;       /**< size to which \a arr is currently allocated*/
            T *arr     = nullptr; /**< internal array */
        };

        typedef worklist_t<int> worklist;

        /**
         * \brief Fixed-size integer array, 0-initialized
         *
         * An array of fixed size, with some further convenience functions.
         *
         */
        class workspace {
        private:
            /**
             * Allocate an array of size \p size.
             * @param size Space to allocate.
             */
            void alloc(const int size) {
                dej_assert(size >= 0);
                dealloc();
                arr = (int*) calloc(size, sizeof(int));
                arr_sz = size;
            }

            void dealloc() {
                if(arr) free(arr);
                arr = nullptr;
            }

        public:
            /**
             * Default constructor, does not allocate any memory.
             */
            workspace() = default;

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit workspace(int size) {
                allocate(size);
            }

            workspace(const workspace& other)  {
                copy(other);
            }

            workspace& operator=(const workspace& other) {
                if(other.arr == arr) return *this;
                copy(other);
                return *this;
            }

            workspace& operator=(workspace&& other) {
                swap(other);
                return *this;
            }

            void swap(workspace& other) {
                auto swap_arr = arr;
                auto swap_arr_sz = arr_sz;
                arr = other.arr;
                arr_sz = other.arr_sz;
                other.arr = swap_arr;
                other.arr_sz = swap_arr_sz;
            }

            void copy(const workspace& other) {
                alloc(other.arr_sz);
                memcpy(arr, other.arr, other.arr_sz  * sizeof(int));
                arr_sz  = other.arr_sz;
            }

            /**
             * Allocates the internal array with size \p size. The allocated memory is not
             * initialized. Initializes the internal position \a cur_pos of the array at 0.
             *
             * @param size Size to allocate.
             */
            void allocate(int size) {
                alloc(size);
            }

            /**
             * Sets the entire array to 0.
             */
            inline void reset() {
                memset(arr, 0, arr_sz * sizeof(int));
            }

            /**
             * Resizes the internal array to `size`. Copies the contents of the old array into the new one. If the new
             * size is larger than the old one, the new space is only allocated and not initialized.
             *
             * @param size New size to allocate the array to.
             */
            void resize(const int size) {
                dej_assert(size >= 0);
                if (arr && size <= arr_sz) return;
                int *old_arr = arr;
                arr = nullptr;
                int old_arr_sz = arr_sz;
                alloc(size);
                if (old_arr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(int));
                    free(old_arr);
                }
            }

            /**
             * Deallocates the internal array.
             */
            ~workspace() {
                dealloc();
            }

            /**
             * @return A pointer to the internal memory.
             */
            dej_nodiscard inline int* get_array() const {
                return arr;
            }

            /**
             * Access element \p index in the internal array \p arr.
             *
             * @param index Index of the internal array.
             * @return The element `arr[index]`.
             */
            inline int &operator[](int index) const {
                dej_assert(index >= 0);
                dej_assert(index < arr_sz);
                return arr[index];
            }

        private:
            int arr_sz = 0;     /**< size to which \a arr is currently allocated*/
            int *arr = nullptr; /**< internal array */
        };


        /**
         * \brief Set with counting
         *
         * Set on a statically specified domain of elements 1, ..., size, where each element holds an additional
         * counter. Time complexity is O(1) for \a set, \a inc, \a inc_nr, \a get. 
         * Amortized O(1) for \a reset if only \a inc and \a get were used.
         */
        template<class T>
        class workset_t {
        public:
            /**
             * Default constructor, does not allocate any memory.
             */
            workset_t() = default;

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit workset_t(int size) {
                initialize(size);
            }

            workset_t(const workset_t<T>& other) {
                copy(other);
            }

            workset_t(const workset_t<T>&& other) {
                swap(other);
            }

            workset_t<T>& operator=(const workset_t<T>& other) {
                if(other.s == s) return *this;
                copy(other);
                return *this;
            }

            workset_t<T>& operator=(workset_t<T>&& other) {
                swap(other);
                return *this;
            }

            void swap(workset_t<T>& other) {
                auto swap_arr = s;
                auto swap_arr_sz = sz;
                s = other.s;
                sz = other.sz;
                other.s = swap_arr;
                other.sz = swap_arr_sz;
                reset_queue.swap(other.reset_queue);
            }

            void copy(const workset_t<T>& other) {
                initialize(other.sz);
                memcpy(s, other.s, other.sz  * sizeof(T));
                sz  = other.sz;
                reset_queue.copy(&other.reset_queue);
            }

            void initialize(int size) {
                if(init) delete[] s;
                s = new T[size];
                reset_queue.allocate(size);

                memset(s, -1, size * sizeof(T));

                init = true;
                sz = size;
            }

            /** 
             * Set entry at position \p index to given \p value. Will not be reset using \a reset method.
             *
             * @param index the position of entry 
             * @param value value to set the entry to
             *
             */
            void set(int index, T value) {
                dej_assert(index >= 0);
                dej_assert(index < sz);
                s[index] = value;
            }

            /** 
             * Get value at position \p index.
             *
             * @param index the position of entry 
             *
             */
            T get(int index) {
                dej_assert(index >= 0);
                dej_assert(index < sz);
                return s[index];
            }

            /**
             * Reset all entries that were manipulated using \a inc method.
             */
            void reset() {
                while (!reset_queue.empty()) s[reset_queue.pop_back()] = -1;
            }

            /**
             * Reset all entries.
             */
            void reset_hard() {
                memset(s, -1, sz * sizeof(T));
                reset_queue.reset();
            }

            /**
             * Increment entry at position \p index by 1.
             *
             * @param index the position of entry that will be manipulated
             *
             */
            T inc(int index) {
                dej_assert(index >= 0);
                dej_assert(index < sz);
                if (s[index]++ == -1)
                    reset_queue.push_back(index);
                return s[index];
            }

            /**
             * Increment entry at position \p index by 1, but do not keep track of reset information.
             *
             * @param index the position of entry that will be manipulated
             *
             */
            void inline inc_nr(int index) {
                dej_assert(index >= 0 && index < sz);
                ++s[index];
            }

            ~workset_t() {
                if (init)
                    delete[] s;
            }

        private:
            worklist reset_queue;
            bool init = false;
            T   *s = nullptr;
            int sz = 0;
        };

        typedef workset_t<int> work_set_int;

        /**
         * \brief Set specialized for quick resets
         *
         * Set on a statically specified domain of elements 1, ..., size, with O(1) \a set and \a get.
         */
        class markset {
            int *s   = nullptr;
            int mark = 0;
            int sz = 0;

            void full_reset() {
                memset(s, mark, sz * sizeof(int));
            }

        public:
            /**
             * Initializes a set of size 0
             */
            markset() = default;

            markset(const markset& other)  {
                copy(&other);
            }

            markset(markset&& other)  {
                swap(other);
            }

            markset& operator=(const markset& other) {
                if(other.s == s) return *this;
                copy(&other);
                return *this;
            }

            /**
             * Initialize this set with the given \p size.
             * @param size size to initialize this set to
             */
            explicit markset(int size) {
                initialize(size);
            }

            /**
             * Resizes this set to size \p size, resets the set
             *
             * @param size new size of this set
             */
            void initialize(int size) {
                dej_assert(size >= 0);
                if(s && sz == size) {
                    reset();
                    return;
                }
                if(s) free(s);
                s = (int*) calloc((unsigned int) size, sizeof(int));
                sz   = size;
                mark = 0;
                reset();
            }

            /**
             * @param pos element to check
             * @return Is element \p pos in set?
             */
            inline bool get(int pos) {
                dej_assert(pos >= 0);
                dej_assert(pos < sz);
                return s[pos] == mark;
            }

            /**
             * Adds element \p pos to set
             *
             * @param pos element to set
             */
            inline void set(int pos) {
                dej_assert(pos >= 0);
                dej_assert(pos < sz);
                s[pos] = mark;
            }

            /**
             * Removes element \p pos from set
             * @param pos element to remove
             */
            inline void unset(int pos) {
                dej_assert(pos >= 0);
                dej_assert(pos < sz);
                s[pos] = mark - 1;
            }
            /**
             * Resets this set to the empty set
             */
            void reset() {
                if(mark == -1) full_reset();
                ++mark;
            }

            void swap(markset& other) {
                auto swap_s = s;
                auto swap_mark = mark;
                auto swap_sz = sz;
                s = other.s;
                mark = other.mark;
                sz = other.sz;
                other.s = swap_s;
                other.mark = swap_mark;
                other.sz = swap_sz;
            }

            void copy(const markset* other) {
                initialize(other->sz);
                for(int i = 0; i < other->sz; ++i) {
                    s[i] = other->s[i];
                }
                mark = other->mark;
                sz   = other->sz;
            }

            ~markset() {
                if(s) free(s);
            }
        };

//        static void bucket_sort(worklist& list,  markset& buckets, int limit) {
//            buckets.reset();
//            for(int i = 0; i < list.size(); ++i) buckets.set(list[i]);
//            list.reset();
//            for(int i = 0; i <= limit; ++i) if(buckets.get(i)) list.push_back(i);
//            buckets.reset();
//        }
    }
}

#endif //DEJAVU_DS_H
