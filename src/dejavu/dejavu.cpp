// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include <iostream>
#include "dejavu.h"
#include <chrono>
#include <string>

typedef std::chrono::high_resolution_clock Clock;

dejavu::ir::refinement test_r;
dejavu::sgraph dej_test_graph;

void empty_hook(int, const int*, int, const int *) {}

int commandline_mode(int argc, char **argv) {
    std::string filename;
    bool entered_file = false;
    bool permute_graph = false;
    bool permute_graph_have_seed  = false;
    int  permute_graph_given_seed = 0;
    bool print = true;

    bool true_random = false;
    bool true_random_seed = false;

    int error_bound = 10;

    bool write_grp_sz = false;
    bool write_benchmark_lines = false;
    bool write_auto_stdout = false;
    bool        write_auto_file      = false;
    std::string write_auto_file_name;

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
        std::replace(arg.begin(), arg.end(), '-', '_');

        if (arg == "__HELP" || arg == "_H") {
            std::cout << "Usage: dejavu [file] [options]" << std::endl;
            std::cout << "Computes the automorphism group of undirected graph described in FILE." << std::endl;
            std::cout << "FILE is expected to be in DIMACS format." << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--err [n]" << std::setw(16) <<
            "Sets the error to be bounded by 1/2^N, assuming uniform random numbers" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--silent" << std::setw(16) <<
            "Does not print progress of the solver" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--gens" << std::setw(16) <<
            "Prints found generators line-by-line to console" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--gens-file [f]" << std::setw(16) <<
           "Writes found generators line-by-line to file F" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--grp-sz" << std::setw(16) <<
            "Prints group size to console (even if --silent)" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--pseudo-random" << std::setw(16) <<
            "Uses pseudo random numbers (default)" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--true-random" << std::setw(16) <<
            "Uses random device of OS" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--true-random-seed" << std::setw(16) <<
            "Seeds pseudo random with random device of OS" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--permute" << std::setw(16) <<
            "Randomly permutes the given graph" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--permute-seed [n]" << std::setw(16) <<
            "Seed for the previous option with N" << std::endl;
            return 0;
        } else if (arg == "__VERSION" || arg == "_V") {
            std::cout << DEJAVU_VERSION_MAJOR << "." << DEJAVU_VERSION_MINOR <<
                        (DEJAVU_VERSION_IS_PREVIEW?"preview":"") << std::endl;
            return 0;
        } else if (arg == "__FILE") {
            if (i + 1 < argc) {
                i++;
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "--file option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__ERR") {
            if (i + 1 < argc) {
                i++;
                error_bound = atoi(argv[i]);
            } else {
                std::cerr << "--err option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__GRP_SZ") {
            write_grp_sz = true;
        }  else if (arg == "__BENCHMARK_LINES") {
            write_benchmark_lines = true;
        } else if (arg == "__GENS") {
            write_auto_stdout = true;
        }  else if (arg == "__GENS_FILE") {
            if (i + 1 < argc) {
                i++;
                write_auto_file = true;
                write_auto_file_name = argv[i];
            } else {
                std::cerr << "--write-gens-file option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__TRUE_RANDOM") {
            if(true_random_seed) {
                std::cerr << "--true-random and --true-random-seed can not be activated at the same time:" <<
                          "--true-random-seed seeds a pseudo random number generator with a random number" << std::endl;
                return 1;
            }
            true_random = true;
        } else if (arg == "__PSEUDO_RANDOM") {
            true_random = false;
        }  else if (arg == "__TRUE_RANDOM_SEED") {
            if(true_random) {
                std::cerr << "--true-random and --true-random-seed can not be activated at the same time:" <<
                "--true-random-seed seeds a pseudo random number generator with a random number" << std::endl;
                return 1;
            }
            true_random_seed = true;
        } else if (arg == "__PERMUTE") {
            permute_graph = true;
        }  else if (arg == "__PERMUTE_SEED") {
            if (i + 1 < argc) {
                i++;
                permute_graph_have_seed  = true;
                permute_graph_given_seed = atoi(argv[i]);
            } else {
                std::cerr << "--permute_seed option requires one argument." << std::endl;
                return 1;
            }
        }  else if (arg == "__SILENT") {
            print = false;
        }  else if (argv[i][0] != '-') {
            if(!entered_file) {
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "Extraneous file '" << argv[i] << "'. Only 1 file required." << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Invalid commandline option '" << argv[i] << "'." << std::endl;
            return 1;
        }
    }

    if (!entered_file) {
        std::cerr << "no file was specified, usage: dejavu [file] [options], use --help for options" << std::endl;
        return 1;
    }

    if(!file_exists(filename)) {
        std::cerr << "File '" << filename << "' does not exist." << std::endl;
        return 1;
    }

    if(print) std::cout << "dejavu version=" << DEJAVU_VERSION_MAJOR << "." << DEJAVU_VERSION_MINOR <<
                        (DEJAVU_VERSION_IS_PREVIEW?"preview":"") << std::endl;
    if(print) std::cout << "------------------------------------------------------------------" << std::endl;

    dejavu::sgraph g;
    if(print) std::cout << "parsing '" << filename << "'" << std::endl;
    int* colmap = nullptr;

    int permute_seed = 0;
    if(permute_graph) {
        permute_seed = permute_graph_given_seed;
        if(!permute_graph_have_seed) {
            std::random_device r;
            permute_seed = static_cast<int>(r());
        }
        if(print) std::cout << (true_random?"true_random=true, ":"") << (true_random_seed?"true_random_seed=true":"");
        if(print) std::cout << "permutation_seed=" << permute_seed << ", ";
    }
    const bool parse_success = parse_dimacs(filename, &g, &colmap, !print, permute_seed);
    if(!parse_success) return 1;
    if(print) std::cout << ", n=" << g.v_size << ", " << "m=" << g.e_size/2 << std::endl << std::endl;

    // manage hooks
    auto empty_hook_func = dejavu_hook(empty_hook);
    dejavu::hooks::multi_hook hooks;
    std::ofstream output_file;
    dejavu::hooks::ostream_hook file_hook(output_file);
    dejavu::hooks::ostream_hook cout_hook(std::cout);
    dejavu_hook* hook;

    // write automorphism to file or cout
    if(write_auto_stdout) hooks.add_hook(cout_hook.get_hook());
    if(write_auto_file) {
        output_file.open(write_auto_file_name);
        hooks.add_hook(file_hook.get_hook());
    }

    // debug hook
#ifndef NDEBUG
#ifdef DEJDEBUG
    auto test_hook_func = dejavu_hook(dejavu::test_hook);
    dej_test_graph.copy_graph(&g);
    hooks.add_hook(&test_hook_func);
#endif
#endif

    // use multi-hook or empty hook
    if(hooks.size() == 0) hook = &empty_hook_func; // using empty_hook_func for fair benchmarks, 'nullptr' is faster
    else hook = hooks.get_hook();

    // no coloring given? let's insert the trivial coloring
    if (colmap == nullptr) colmap = (int *) calloc(g.v_size, sizeof(int));

    // now run the solver with the given options...
    Clock::time_point timer = Clock::now();

    dejavu::solver d;
    d.set_error_bound(error_bound);
    d.set_print(print);
    if (true_random_seed) d.randomize_seed();
    d.set_true_random(true_random);
    d.automorphisms(&g, colmap, hook);

    long dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    dejavu::big_number grp_sz = d.get_automorphism_group_size();

    if (print) std::cout << "------------------------------------------------------------------" << std::endl;
    if (print || write_benchmark_lines)
        std::cout << std::setprecision(4) << "symmetries=" << grp_sz
                  << ", deterministic=" << (d.get_deterministic_termination() ? "true" : "false")
                  << ", error=1/2^" << d.get_error_bound() << "," << std::endl;

    if(print || write_benchmark_lines) std::cout << "solve_time=" <<
                                       static_cast<double>(dejavu_solve_time) / 1000000.0 << "ms" << std::endl;
    if(!print && write_grp_sz) std::cout << grp_sz << std::endl;

    free(colmap);
    return 0;
}

int main(int argc, char *argv[]) {
    return commandline_mode(argc, argv);
}
