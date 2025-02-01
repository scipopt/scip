// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include <iostream>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <set>
#include <cstring>
#include <queue>
#include <memory>
#include <chrono>
#include <iomanip>
#include <functional>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
    #define OS_WINDOWS
#endif

#if __APPLE__
    #define OS_MAC
#endif

#if __linux__
    #define OS_LINUX
#endif

#define PRINT_NO_NEWLINE(str) std::cout << str << std::flush;
#define PRINT(str) std::cout << str << std::endl;

#ifdef DEJDEBUG
#define dej_assert(expr) (assert(expr))
#else
#define dej_assert(expr) (void)0
#endif

#if ((defined(_MSVC_LANG) && _MSVC_LANG > 201402L) || __cplusplus > 201402L)
#define dej_nodiscard [[nodiscard]]
#else
#define dej_nodiscard
#endif

// use these options to prevent dejavu from using certain C++ language features
// #define dej_nolambda
// #define dej_nothreadlocal


/**
 * Hash function for unsigned integers.
 *
 * @param x the unsigned integer
 * @return hashed integer
 */
static inline unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

/**
 * Accumulate a hash, for example to be used to hash strings of integers.
 *
 * @param hash hash computed so far
 * @param d integer to accumulate to \p hash
 * @return the new hash
 */
static inline unsigned long add_to_hash(unsigned long hash, const int d) {
    const unsigned long ho = hash & 0xff00000000000000; // extract high-order 8 bits from hash
    hash    = hash << 8;                    // shift hash left by 5 bits
    hash    = hash ^ (ho >> 56);            // move the highorder 5 bits to the low-order
    hash    = hash ^ d;                     // XOR into hash

    return hash;
}

/**
 * Does the file with filename \p name exist?
 *
 * @param name filename to look for
 * @return whether file \p name exists
 */
static inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

typedef void type_dejavu_hook(int, const int*, int, const int*);
typedef std::function<void(int, const int*, int, const int*)> dejavu_hook;

namespace dejavu {

    /**
     * \brief Random number generation
     *
     * Either uses a pseudo random number generator, or random device of OS.
     */
    class random_source {
        bool true_random = false;
        std::mt19937 pseudo_random_device;
        std::random_device true_random_device;
    public:
        /**
         * Creates random source with given parameters.
         *
         * @param set_true_random use random device of OS if true, or pseudo random if false
         * @param set_seed sets seed of pseudo random number generator
         */
        random_source(bool set_true_random, int set_seed) {
            true_random = set_true_random;
            pseudo_random_device.seed(set_seed);
        }

        /**
         * Returns a random number
         * @return
         */
        int operator()() {
            return true_random?static_cast<int>(true_random_device()&INT32_MAX):
                               static_cast<int>(pseudo_random_device()&INT32_MAX);
        }
    };

    /**
     * \brief Stores big numbers
     *
     * A simple class to store big, positive numbers. Consists of a \a mantissa and a \a exponent, where the value of
     * the number is `mantissa^exponent`.
     *
     * Used to store automorphism group sizes.
     */
    class big_number {
    public:
        long double mantissa = 1.0; /**< mantissa, number is `mantissa^exponent`
                                          * \sa exponent */
        int         exponent = 0;   /**< exponent, number is `mantissa^exponent`
                                          * \sa mantissa  */

        friend bool operator<(const big_number& l, const big_number& r)
        {
            return (l.exponent < r.exponent) || (l.exponent == r.exponent && l.mantissa+0.01 < r.mantissa);
        }

        friend bool operator==(const big_number& l, const big_number& r)
        {
            return (l.exponent == r.exponent) && (l.mantissa > r.mantissa-0.01) && (l.mantissa < r.mantissa+0.01);
        }

        void set(long double set_mantissa, int set_exponent) {
            mantissa = set_mantissa;
            exponent = set_exponent;
        }

        /**
         * Multiply a \p number to this big_number.
         *
         * @param number The number to multiply.
         */
        void multiply(int number) {
            multiply(number, 0);
        }

        /**
         * Multiply a \p number to this big_number.
         *
         * @param number The number to multiply.
         */
        void multiply(big_number number) {
            multiply(number.mantissa, number.exponent);
        }


        /**
         * Multiply a number consisting of a mantissa (\p other_mantissa) and exponent (\p other_exponent) to this
         * big_number.
         *
         * @param other_mantissa Mantissa of number to multiply.
         * @param other_exponent Exponent of number to multiply.
         */
        void multiply(long double other_mantissa, int other_exponent) {
            if(std::fpclassify(other_mantissa) == FP_INFINITE ||  std::fpclassify(other_mantissa) == FP_NAN) {
                return;
            }
            while (other_mantissa >= 10.0) {
                exponent += 1;
                other_mantissa = other_mantissa / 10;
            }
            exponent += other_exponent;
            mantissa *= other_mantissa;
            man_to_exp();
        }

    private:
        void man_to_exp() {
            if(std::fpclassify(mantissa) == FP_INFINITE ||  std::fpclassify(mantissa) == FP_NAN) {
                return;
            }
            while(mantissa >= 10.0) {
                exponent += 1;
                mantissa = mantissa / 10;
            }
        }
    };

    inline std::ostream& operator<<(std::ostream& out, big_number number) {
        return out << number.mantissa << "*10^" << number.exponent;
    }

    /**
     * Used to make the output look a bit more structured.
     */
    static void progress_print_split() {
        PRINT("\r______________________________________________________________");
    }

    /**
     * Prints heading of table in the output.
     */
    static void progress_print_header() {
        progress_print_split();
        PRINT(std::setw(11) << std::left <<"T (ms)" << std::setw(11) << "delta(ms)" << std::setw(12) << "proc"
              << std::setw(16) << "p1"        << std::setw(16)        << "p2");
        progress_print_split();
    }

    /**
     * \brief Prints information to the console.
     *
     * Contains additional facilities to measure elapsed time in-between prints.
     */
    class timed_print {
        std::chrono::high_resolution_clock::time_point first;
        std::chrono::high_resolution_clock::time_point previous;
    public:

        bool h_silent = false;

        timed_print() {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
        }

        void print_header() const {
            if(h_silent) return;
            progress_print_header();
        }

        void print_split() const {
            if(h_silent) return;
            progress_print_split();
        }

        void print(const std::string& str) const {
            if(h_silent) return;
            PRINT("\r" << str);
        }

        void timer_print(const std::string& proc, const std::string& p1, const std::string& p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void timer_split() {
            previous = std::chrono::high_resolution_clock::now();
        }

        void timer_print(const std::string& proc, const int p1, const int p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void timer_print(const std::string& proc, const int p1, const double p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void progress_current_method(const std::string& print) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << print);
        }
        void progress_current_method(const std::string& method_name, const std::string& var1, double var1_val,
                                            const std::string& var2, double var2_val) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " " << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val);
        }
        void progress_current_method(const std::string& method_name, const std::string& var1, int var1_val,
                                            const std::string& var2, int var2_val,
                                            const std::string& var3, double var3_val) const {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " "  << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val
                                   << ", " << var3 << "=" << var3_val);
        }

        void progress_current_method(const std::string& method_name, const std::string& var1, double var1_val,
                                            const std::string& var2, int var2_val,
                                            const std::string& var3, int var3_val,
                                            const std::string& var4, int var4_val) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " "  << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val
                                   << ", " << var3 << "=" << var3_val << ", " << var4 << "=" << var4_val);
        }
    };
}

#endif //DEJAVU_UTILITY_H
