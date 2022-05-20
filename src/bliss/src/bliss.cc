/*
  Copyright (c) 2003-2021 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <functional>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include "bliss/defs.hh"
#include "bliss/timer.hh"
#include "bliss/utils.hh"
#include "bliss/graph.hh"
#include "bliss/digraph.hh"

/**
 * \page executable The bliss executable
 * \include bliss.cc
 */

/* Input file name */
static const char* infilename = 0;

static bool opt_directed = false;
static bool opt_canonize = false;
static const char* opt_output_can_file = 0;
static const char* opt_splitting_heuristics = "fsm";
static bool opt_use_failure_recording = true;
static bool opt_use_component_recursion = true;


/* Verbosity level and target stream */
static unsigned int verbose_level = 1;
static FILE* verbstr = stdout;



static void
usage(FILE* const fp, const char* argv0)
{
  const char* program_name = strrchr(argv0, '/');

  if(program_name) program_name++;
  else program_name = argv0;
  if(!program_name or *program_name == 0) program_name = "bliss";

  fprintf(fp, "bliss version %s (compiled %s)\n", bliss::version, __DATE__);
  fprintf(fp, "Copyright 2003-2021 Tommi Junttila\n");
  fprintf(fp,
"\n"
"Usage: %s [options] [<graph file>]\n"
"\n"
"  -directed   the input graph is directed\n"
"  -can        compute canonical form\n"
"  -ocan=f     compute canonical form and output it in file f\n"
"  -v=N        set verbose level to N [N >= 0, default: 1]\n"
"  -sh=X       select splitting heuristics, where X is\n"
"                f    first non-singleton cell\n"
"                fl   first largest non-singleton cell\n"
"                fs   first smallest non-singleton cell\n"
"                fm   first maximally non-trivially connected\n"
"                     non-singleton cell\n"
"                flm  first largest maximally non-trivially connected\n"
"                     non-singleton cell\n"
"                fsm  first smallest maximally non-trivially connected\n"
"                     non-singleton cell [default]\n"
"  -fr=X       use failure recording? [X=y/n, default: y]\n"
"  -cr=X       use component recursion? [X=y/n, default: y]\n"
"  -version    print the version number and exit\n"
"  -help       print this help and exit\n"
          ,program_name
          );
}



static void
parse_options(const int argc, const char** argv)
{
  unsigned int tmp;
  for(int i = 1; i < argc; i++)
    {
      if(strcmp(argv[i], "-can") == 0)
        opt_canonize = true;
      else if((strncmp(argv[i], "-ocan=", 6) == 0) and (strlen(argv[i]) > 6))
        {
          opt_canonize = true;
          opt_output_can_file = argv[i]+6;
        }
      else if(sscanf(argv[i], "-v=%u", &tmp) == 1)
        verbose_level = tmp;
      else if(strcmp(argv[i], "-directed") == 0)
        opt_directed = true;
      else if(strcmp(argv[i], "-fr=n") == 0)
        opt_use_failure_recording = false;
      else if(strcmp(argv[i], "-fr=y") == 0)
        opt_use_failure_recording = true;
      else if(strcmp(argv[i], "-cr=n") == 0)
        opt_use_component_recursion = false;
      else if(strcmp(argv[i], "-cr=y") == 0)
        opt_use_component_recursion = true;
      else if((strncmp(argv[i], "-sh=", 4) == 0) and (strlen(argv[i]) > 4))
        {
          opt_splitting_heuristics = argv[i]+4;
        }
      else if(strcmp(argv[i], "-version") == 0)
        {
          fprintf(stdout, "bliss version %s\n", bliss::version);
          exit(0);
        }
      else if(strcmp(argv[i], "-help") == 0)
        {
          usage(stdout, argv[0]);
          exit(0);
        }
      else if(argv[i][0] == '-')
        {
          fprintf(stderr, "Unknown command line argument `%s'\n", argv[i]);
          usage(stderr, argv[0]);
          exit(1);
        }
      else
        {
          if(infilename)
            {
              fprintf(stderr, "Too many file arguments\n");
              usage(stderr, argv[0]);
              exit(1);
            }
          else
            {
              infilename = argv[i];
            }
        }
    }
}



/* Output an error message and exit the whole program with the exit value 1. */
static void
_fatal(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap); fprintf(stderr, "\n");
  va_end(ap);
  exit(1);
}



int
main(const int argc, const char** argv)
{
#ifndef _WIN32
  bliss::Timer timer;
#endif
  bliss::AbstractGraph* g = 0;

  parse_options(argc, argv);

  /* Parse splitting heuristics */
  bliss::Digraph::SplittingHeuristic shs_directed = bliss::Digraph::shs_fsm;
  bliss::Graph::SplittingHeuristic shs_undirected = bliss::Graph::shs_fsm;
  if(opt_directed)
    {
      if(strcmp(opt_splitting_heuristics, "f") == 0)
        shs_directed = bliss::Digraph::shs_f;
      else if(strcmp(opt_splitting_heuristics, "fs") == 0)
        shs_directed = bliss::Digraph::shs_fs;
      else if(strcmp(opt_splitting_heuristics, "fl") == 0)
        shs_directed = bliss::Digraph::shs_fl;
      else if(strcmp(opt_splitting_heuristics, "fm") == 0)
        shs_directed = bliss::Digraph::shs_fm;
      else if(strcmp(opt_splitting_heuristics, "fsm") == 0)
        shs_directed = bliss::Digraph::shs_fsm;
      else if(strcmp(opt_splitting_heuristics, "flm") == 0)
        shs_directed = bliss::Digraph::shs_flm;
      else
        _fatal("Illegal option -sh=%s, aborting", opt_splitting_heuristics);
    }
  else
    {
      if(strcmp(opt_splitting_heuristics, "f") == 0)
        shs_undirected = bliss::Graph::shs_f;
      else if(strcmp(opt_splitting_heuristics, "fs") == 0)
        shs_undirected = bliss::Graph::shs_fs;
      else if(strcmp(opt_splitting_heuristics, "fl") == 0)
        shs_undirected = bliss::Graph::shs_fl;
      else if(strcmp(opt_splitting_heuristics, "fm") == 0)
        shs_undirected = bliss::Graph::shs_fm;
      else if(strcmp(opt_splitting_heuristics, "fsm") == 0)
        shs_undirected = bliss::Graph::shs_fsm;
      else if(strcmp(opt_splitting_heuristics, "flm") == 0)
        shs_undirected = bliss::Graph::shs_flm;
      else
        _fatal("Illegal option -sh=%s, aborting", opt_splitting_heuristics);
    }

  /* Open the input file */
  FILE* infile = stdin;
  if(infilename)
    {
      infile = fopen(infilename, "r");
      if(!infile)
        _fatal("Cannot not open `%s' for input, aborting", infilename);
    }

  /* Read the graph from the file */
  if(opt_directed)
    {
      /* Read directed graph in the DIMACS format */
      g = bliss::Digraph::read_dimacs(infile);
    }
  else
    {
      /* Read undirected graph in the DIMACS format */
      g = bliss::Graph::read_dimacs(infile);
    }

  if(infile != stdin)
    fclose(infile);

  if(!g)
  {
    _fatal("Failed to read the graph, aborting");
    exit(1); /* this is here because scan-build does not recognize _fatal exits. */
  }

#ifndef _WIN32
  if(verbose_level >= 2)
    {
      fprintf(verbstr, "Graph read in %.2f seconds\n", timer.get_duration());
      fflush(verbstr);
    }
#endif


  bliss::Stats stats;

  /* Set splitting heuristics and verbose level */
  if(opt_directed)
    ((bliss::Digraph*)g)->set_splitting_heuristic(shs_directed);
  else
    ((bliss::Graph*)g)->set_splitting_heuristic(shs_undirected);
  g->set_verbose_level(verbose_level);
  g->set_verbose_file(verbstr);
  g->set_failure_recording(opt_use_failure_recording);
  g->set_component_recursion(opt_use_component_recursion);


  auto report_aut = [&](const unsigned int n, const unsigned int* aut) -> void {
    fprintf(stdout, "Generator: ");
    bliss::print_permutation(stdout, n, aut, 1);
    fprintf(stdout, "\n");
  };

  if(opt_canonize == false)
    {
      /* No canonical labeling, only automorphism group */
      g->find_automorphisms(stats, report_aut);
    }
  else
    {
      /* Canonical labeling and automorphism group */
      const unsigned int* cl = g->canonical_form(stats, report_aut);

      fprintf(stdout, "Canonical labeling: ");
      bliss::print_permutation(stdout, g->get_nof_vertices(), cl, 1);
      fprintf(stdout, "\n");

      if(opt_output_can_file)
        {
          bliss::AbstractGraph* cf = g->permute(cl);
          FILE* const fp = fopen(opt_output_can_file, "w");
          if(!fp)
            _fatal("Cannot open '%s' for outputting the canonical form, aborting", opt_output_can_file);
          cf->write_dimacs(fp);
          fclose(fp);
          delete cf;
        }
    }

  /* Output search statistics */
  if(verbose_level > 0 and verbstr)
    stats.print(verbstr);

#ifndef _WIN32
  if(verbose_level > 0)
    {
      fprintf(verbstr, "Total time:\t%.2f seconds\n", timer.get_duration());
      fflush(verbstr);
    }
#endif


  delete g; g = 0;

  return 0;
}
