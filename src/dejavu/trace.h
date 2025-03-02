// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_TRACE_H
#define DEJAVU_TRACE_H

#include <vector>
#include <cstdint>
#include "utility.h"

namespace dejavu {
    namespace ir {
#define TRACE_MARKER_INDIVIDUALIZE     (INT32_MAX-7)
#define TRACE_MARKER_REFINE_START      (INT32_MAX-2)
#define TRACE_MARKER_REFINE_END        (INT32_MAX-3)
#define TRACE_MARKER_REFINE_CELL_START (INT32_MAX-4)
#define TRACE_MARKER_REFINE_CELL_END   (INT32_MAX-5)

        /**
         * \brief The trace invariant.
         *
         * Class that serves to store and compare the trace of a walk in an individualization-refinement shared_tree.
         * The class provides several different modes in which information is recorded and/or compared.
         *
         * Specifically, it is possible to (1) record a full trace, (2) compare to a full trace, or (3) compare to a
         * full trace and recording a hash invariant as soon as the new computation deviates from the stored trace (see
         * also \ref ir_mode ).
         *
         * While comparing to a stored trace (2, 3), the class facilitates the use of the blueprint heuristic, which
         * enables skipping of non-splitting cells in the stored trace.
        */
        class trace {
        private:
            // the trace
            std::vector<int> data; /**< keeps all the data of the trace */

            trace *compare_trace = nullptr; /**< link to a stored trace to compare to */
            uint64_t hash = 0; /**< hash value to summarize all operations performed on this trace */

            // mode
            bool compare = false; /**< whether to compare operations to a stored trace*/
            bool record  = false; /**< whether to record a trace */

            // housekeeping
            // int cell_act_spot      = -1;
            // int cell_old_color     = -1;
            bool assert_cell_act   = false;
            bool assert_refine_act = false;

            // comparison variables
            int position = 0;
            bool comp = true;

            void inline add_to_hash(int d) {
                uint64_t ho = hash & 0xff00000000000000; // extract high-order 8 bits from hash
                hash    = hash << 8;                    // shift hash left by 5 bits
                hash    = hash ^ (ho >> 56);            // move the highorder 5 bits to the low-order
                hash    = hash ^ d;                     // XOR into hash
            }

            void write_compare(int d) {
                d = std::min(INT32_MAX-10,d);
                dej_assert(d != TRACE_MARKER_INDIVIDUALIZE && d != TRACE_MARKER_REFINE_START);
                write_compare_no_limit(d);
                dej_assert(record?static_cast<int>(data.size())==position:true);
            }

            void write_compare_no_limit(int d) {
                add_to_hash(d);
                if (record)  data.push_back(d);
                if (compare) comp = comp && (position < ((int) compare_trace->data.size()))
                                    && (compare_trace->data[position] == d);
                ++position;
            }

            /*void write_skip_compare(int d) {
                d = std::min(INT32_MAX-10,d);
                dej_assert(d != TRACE_MARKER_INDIVIDUALIZE && d != TRACE_MARKER_REFINE_START);
                if (record) data.push_back(d);
                ++position;
            }*/

        public:

            trace* get_compare_trace() {
                return compare_trace;
            }

            /**
             * Records an individualization.
             * @param color The color being individualized.
             */
            void op_individualize(const int ind_color) {
                dej_assert(ind_color >= 0);
                write_compare_no_limit(TRACE_MARKER_INDIVIDUALIZE);
                write_compare(ind_color);
            }

            /**
             * Records the start of a refinement.
             */
            void op_refine_start() {
                dej_assert(!comp || !assert_refine_act);
                write_compare_no_limit(TRACE_MARKER_REFINE_START);
                assert_refine_act = true;
            }

            /**
             * Records the start of a refinement with respect to a color.
             * @param color The color in respect to which the coloring is refined.
             */
            void op_refine_cell_start(int) {
                dej_assert(!comp || !assert_cell_act);
                write_compare_no_limit(TRACE_MARKER_REFINE_CELL_START);
                //write_compare(color);
                // cell_old_color = color;
                //cell_act_spot = (int) data.size();
                //write_skip_compare(false);
                assert_cell_act = true;
            }

            /**
             * Records a that a new color appeared while refining with respect to a color.
             * @param new_color The new color that was refined.
             */
            void op_refine_cell_record(int new_color) {
                dej_assert(!comp || assert_cell_act);
                write_compare(new_color);
                //if (new_color != cell_old_color && record) data[cell_act_spot] = true;
            }

            void op_additional_info(int d) {
                write_compare(d);
            }

            /**
             * Records the end of a refinement with respect to a color.
             */
            void op_refine_cell_end() {
                dej_assert(!comp || assert_cell_act);
                assert_cell_act = false;
                //
                write_compare_no_limit(TRACE_MARKER_REFINE_CELL_END);
            }

            /**
             * Records the end of a refinement.
             */
            void op_refine_end() {
                dej_assert(!comp || !assert_cell_act);
                dej_assert(assert_refine_act);
                assert_refine_act = false;
                write_compare_no_limit(TRACE_MARKER_REFINE_END);
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * @return Determines whether in the stored trace, the next color in respect to which was refined created new colors
             * (i.e., whether the next color is splitting).
             */
            dej_nodiscard bool blueprint_is_next_cell_active() {
                if (!compare || !comp || position > static_cast<int>(compare_trace->data.size())) return true;

                dej_assert(compare_trace);
                size_t read_pt = position;
                dej_assert(compare_trace->data.size() > read_pt);
                for(; read_pt > 0 && compare_trace->data[read_pt] != TRACE_MARKER_REFINE_CELL_START; --read_pt);
                dej_assert(compare_trace->data[read_pt] == TRACE_MARKER_REFINE_CELL_START);
                ++read_pt;

                dej_assert(compare_trace->data.size() > read_pt);
                dej_assert((compare_trace->data[read_pt] == false) || (compare_trace->data[read_pt] == true));
                return compare_trace->data[read_pt];
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * Skips the \a position to the start of the next refinement with respect to a color. To be used after
             * \a blueprint_is_next_cell_active() determined the current color to be non-splitting.
             */
            void blueprint_skip_to_next_cell() {
                while (position < static_cast<int>(compare_trace->data.size()) &&
                       compare_trace->data[position] != TRACE_MARKER_REFINE_CELL_END) {
                    dej_assert(compare_trace->data.size() > (size_t) position);
                    ++position;
                }
                ++position;
                assert_cell_act = false;
            }

            /**
             * Rewinds the \a position to the previous individualization.
             */
            void rewind_to_individualization() {
                assert_cell_act = false;
                assert_refine_act = false;
                if (record) {
                    int read_pt = std::max((int) data.size() - 1, 0);
                    while (read_pt > 0 && data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        --read_pt;
                    }
                    data.resize(read_pt);
                    position = read_pt;
                }
                if (compare) {
                    int read_pt = std::max(position - 1, 0);
                    while (read_pt > 0 && compare_trace->data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        --read_pt;
                    }
                    position = read_pt;
                }
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * Skips the \a position to the next individualization in the \a compare_trace.
             */
            void skip_to_individualization() {
                assert_cell_act = false;
                assert_refine_act = false;
                if (compare) {
                    int read_pt = position - 1;
                    while ((size_t) read_pt < compare_trace->data.size() &&
                            compare_trace->data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        ++read_pt;
                    }
                    position = read_pt;
                }
            }

            /**
             * Stores a new trace to compare to.
             *
             * @param compare_trace
             */
            void set_compare_trace(trace *new_compare_trace) {
                this->compare_trace = new_compare_trace;
            }

            /**
             * Determines whether we want to compare the following operations with those of the trace stored in
             * \a compare_trace.
             *
             * @param compare
             */
            void set_compare(bool new_compare) {
                this->compare = new_compare;
            }

            /**
             * @return A hash value summarizing the operations recorded in this trace.
             */
            dej_nodiscard uint64_t get_hash() const {
                return hash;
            }

            /**
             * Sets the hash value to a pre-determined value.
             * @param hash The hash value.
             */
            void set_hash(uint64_t new_hash) {
                this->hash = new_hash;
            }

            /**
             * @return Whether the recorded operations deviated from the stored trace in \a compare_trace.
             */
            dej_nodiscard bool trace_equal() const {
                return comp;
            }

            /**
             * Resets the trace to find new deviations from the stored trace.
             */
            void reset_trace_equal() {
                comp = true;
            }
            /**
 * Resets the trace to find new deviations from the stored trace.
 */
            void reset_trace_unequal() {
                comp = false;
            }

            void reset() {
                data.clear();
                compare_trace = nullptr;
                position = 0;
                hash = 0;
                dej_assert(record?static_cast<int>(data.size())==position:true);
                reset_trace_equal();
            }

            // record to trace
            void set_record(bool new_record) {
                this->record = new_record;
            }

            void reserve(const int n) {
                data.reserve(n);
            }

            // position the trace
            void set_position(int new_position) {
                assert_cell_act   = false;
                assert_refine_act = false;
                this->position = new_position;
                if(record) data.resize(position);
                dej_assert(record?static_cast<int>(data.size())==position:true);
            }

            dej_nodiscard int get_position() const {
                return position;
            }
        };
    }
}


#endif //DEJAVU_TRACE_H
