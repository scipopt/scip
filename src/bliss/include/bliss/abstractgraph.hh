#pragma once

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

/**
 * \namespace bliss
 * The namespace bliss contains all the classes and functions of the bliss
 * tool except for the C programming language API.
 */
namespace bliss {
  class AbstractGraph;
}

#include <cstdio>
#include <functional>
#include <vector>
#include "bliss/stats.hh"
#include "bliss/kqueue.hh"
#include "bliss/heap.hh"
#include "bliss/orbit.hh"
#include "bliss/partition.hh"
#include "bliss/uintseqhash.hh"

namespace bliss {


/**
 * \brief An abstract base class for different types of graphs.
 */
class AbstractGraph
{
  friend class Partition;

public:
  AbstractGraph();
  virtual ~AbstractGraph();

  /**
   * Set the verbose output level for the algorithms.
   * \param level  the level of verbose output, 0 means no verbose output
   */
  void set_verbose_level(const unsigned int level);

  /**
   * Set the file stream for the verbose output.
   * \param fp  the file stream; if null, no verbose output is written
   */
  void set_verbose_file(FILE * const fp);

  /**
   * Add a new vertex with color \a color in the graph and return its index.
   */
  virtual unsigned int add_vertex(const unsigned int color = 0) = 0;

  /**
   * Add an edge between vertices \a source and \a target.
   * Duplicate edges between vertices are ignored but try to avoid introducing
   * them in the first place as they are not ignored immediately but will
   * consume memory and computation resources for a while.
   */
  virtual void add_edge(const unsigned int source, const unsigned int target) = 0;

  /**
   * Get the color of a vertex.
   */
  virtual unsigned int get_color(const unsigned int vertex) const = 0;

  /**
   * Change the color of the vertex \a vertex to \a color.
   */
  virtual void change_color(const unsigned int vertex, const unsigned int color) = 0;


  /** Activate/deactivate failure recording.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate failure recording, deactivate otherwise
   */
  void set_failure_recording(const bool active) {
    assert(not in_search);
    opt_use_failure_recording = active;
  }

  /** Activate/deactivate component recursion.
   * The choice affects the computed canonical labelings;
   * therefore, if you want to compare whether two graphs are isomorphic by
   * computing and comparing (for equality) their canonical versions,
   * be sure to use the same choice for both graphs.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate component recursion, deactivate otherwise
   */
  void set_component_recursion(const bool active) {
    assert(not in_search);
    opt_use_comprec = active;
  }



  /**
   * Return the number of vertices in the graph.
   */
  virtual unsigned int get_nof_vertices() const = 0;

  /**
   * Return a new graph that is the result of applying the permutation \a perm
   * to this graph. This graph is not modified.
   * \a perm must contain N=this.get_nof_vertices() elements and be a bijection
   * on {0,1,...,N-1}, otherwise the result is undefined or a segfault.
   */
  virtual AbstractGraph* permute(const unsigned int* const perm) const = 0;

  /**
   * Return a new graph that is the result of applying the permutation \a perm
   * to this graph. This graph is not modified.
   * \a perm must contain N=this.get_nof_vertices() elements and be a bijection
   * on {0,1,...,N-1}, otherwise the result is undefined or a segfault.
   */
  virtual AbstractGraph* permute(const std::vector<unsigned int>& perm) const = 0;

  /**
   * Check whether \a perm is an automorphism of this graph.
   * Unoptimized, mainly for debugging purposes.
   */
  virtual bool is_automorphism(unsigned int* const perm) const = 0;

  /**
   * Check whether \a perm is an automorphism of this graph.
   * Unoptimized, mainly for debugging purposes.
   */
  virtual bool is_automorphism(const std::vector<unsigned int>& perm) const = 0;

  /**
   * Find a set of generators for the automorphism group of the graph.
   * The function \a report (if non-null) is called each time a new generator
   * for the automorphism group is found.
   * The first argument \a n for the function
   * is the length of the automorphism (equal to get_nof_vertices()), and
   * the second argument \a aut is the automorphism
   * (a bijection on {0,...,get_nof_vertices()-1}).
   * The memory for the automorphism \a aut will be invalidated immediately
   * after the return from the \a report function;
   * if you want to use the automorphism later, you have to take a copy of it.
   * Do not call any member functions from the \a report function.
   *
   * The search statistics are copied in \a stats.
   *
   * If the \a terminate function argument is given,
   * it is called in each search tree node: if the function returns true,
   * then the search is terminated and thus not all the automorphisms
   * may have been generated.
   * The \a terminate function may be used to limit the time spent in bliss
   * in case the graph is too difficult under the available time constraints.
   * If used, keep the function simple to evaluate so that
   * it does not consume too much time.
   */
  void find_automorphisms(Stats& stats,
                          const std::function<void(unsigned int n, const unsigned int* aut)>& report = nullptr,
                          const std::function<bool()>& terminate = nullptr);

  /**
   * Otherwise the same as find_automorphisms() except that
   * a canonical labeling of the graph (a bijection on
   * {0,...,get_nof_vertices()-1}) is returned.
   * The memory allocated for the returned canonical labeling will remain
   * valid only until the next call to a member function with the exception
   * that constant member functions (for example, bliss::Graph::permute()) can
   * be called without invalidating the labeling.
   * To compute the canonical version of an undirected graph, call this
   * function and then bliss::Graph::permute() with the returned canonical
   * labeling.
   * Note that the computed canonical version may depend on the applied version
   * of bliss as well as on some other options (for instance, the splitting
   * heuristic selected with bliss::Graph::set_splitting_heuristic()).
   *
   * If the \a terminate function argument is given,
   * it is called in each search tree node: if the function returns true,
   * then the search is terminated and thus (i) not all the automorphisms
   * may have been generated and (ii) the returned labeling may not
   * be canonical.
   * The \a terminate function may be used to limit the time spent in bliss
   * in case the graph is too difficult under the available time constraints.
   * If used, keep the function simple to evaluate so that
   * it does not consume too much time.
   */
  const unsigned int* canonical_form(Stats& stats,
                                     const std::function<void(unsigned int n, const unsigned int* aut)>& report = nullptr,
                                     const std::function<bool()>& terminate = nullptr);

  /**
   * Write the graph to a file in a variant of the DIMACS format.
   * See the <a href="https://users.aalto.fi/tjunttil/bliss">bliss website</a>
   * for the definition of the file format.
   * Note that in the DIMACS file the vertices are numbered from 1 to N while
   * in this C++ API they are from 0 to N-1.
   * Thus the vertex n in the file corresponds to the vertex n-1 in the API.
   * \param fp  the file stream where the graph is written
   */
  virtual void write_dimacs(FILE * const fp) = 0;

  /**
   * Write the graph to a file in the graphviz dotty format.
   * \param fp  the file stream where the graph is written
   */
  virtual void write_dot(FILE * const fp) = 0;

  /**
   * Write the graph in a file in the graphviz dotty format.
   * Do nothing if the file cannot be written.
   * \param file_name  the name of the file to which the graph is written
   */
  virtual void write_dot(const char * const file_name) = 0;

  /**
   * Get a hash value for the graph.
   * \return  the hash value
   */
  virtual unsigned int get_hash() = 0;

  /**
   * Disable/enable the "long prune" method.
   * The choice affects the computed canonical labelings;
   * therefore, if you want to compare whether two graphs are isomorphic by
   * computing and comparing (for equality) their canonical versions,
   * be sure to use the same choice for both graphs.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate "long prune", deactivate otherwise
   */
  void set_long_prune_activity(const bool active) {
    assert(not in_search);
    opt_use_long_prune = active;
  }



protected:
  /** \internal
   * How much verbose output is produced (0 means none) */
  unsigned int verbose_level;

  /** \internal
   * The output stream for verbose output. */
  FILE *verbstr;

  /** \internal
   * The ordered partition used in the search algorithm. */
  Partition p;

  /** \internal
   * Whether the search for automorphisms and a canonical labeling is
   * in progress.
   */
  bool in_search;

  /** \internal
   * Is failure recording in use?
   */
  bool opt_use_failure_recording;

  /* The "tree-specific" invariant value for the point when current path
   * got different from the first path */
  unsigned int failure_recording_fp_deviation;

  /** \internal
   * Is component recursion in use?
   */
  bool opt_use_comprec;


  unsigned int refine_current_path_certificate_index;
  bool refine_compare_certificate;
  bool refine_equal_to_first;
  unsigned int refine_first_path_subcertificate_end;
  int refine_cmp_to_best;
  unsigned int refine_best_path_subcertificate_end;

  static const unsigned int CERT_SPLIT = 0; //UINT_MAX;
  static const unsigned int CERT_EDGE  = 1; //UINT_MAX-1;
  /** \internal
   * Add a triple (v1,v2,v3) in the certificate.
   * May modify refine_equal_to_first and refine_cmp_to_best.
   * May also update eqref_hash and failure_recording_fp_deviation. */
  void cert_add(const unsigned int v1,
                const unsigned int v2,
                const unsigned int v3);

  /** \internal
   * Add a redundant triple (v1,v2,v3) in the certificate.
   * Can also just dicard the triple.
   * May modify refine_equal_to_first and refine_cmp_to_best.
   * May also update eqref_hash and failure_recording_fp_deviation. */
  void cert_add_redundant(const unsigned int x,
                          const unsigned int y,
                          const unsigned int z);

  /**\internal
   * Is the long prune method in use?
   */
  bool opt_use_long_prune;
  /**\internal
   * Maximum amount of memory (in megabytes) available for
   * the long prune method
   */
  static const unsigned int long_prune_options_max_mem = 50;
  /**\internal
   * Maximum amount of automorphisms stored for the long prune method;
   * less than this is stored if the memory limit above is reached first
   */
  static const unsigned int long_prune_options_max_stored_auts = 100;

  unsigned int long_prune_max_stored_autss;
  std::vector<std::vector<bool> *> long_prune_fixed;
  std::vector<std::vector<bool> *> long_prune_mcrs;
  std::vector<bool> long_prune_temp;
  unsigned int long_prune_begin;
  unsigned int long_prune_end;
  /** \internal
   * Initialize the "long prune" data structures.
   */
  void long_prune_init();
  /** \internal
   * Release the memory allocated for "long prune" data structures.
   */
  void long_prune_deallocate();
  void long_prune_add_automorphism(const unsigned int *aut);
  std::vector<bool>& long_prune_get_fixed(const unsigned int index);
  std::vector<bool>& long_prune_allocget_fixed(const unsigned int index);
  std::vector<bool>& long_prune_get_mcrs(const unsigned int index);
  std::vector<bool>& long_prune_allocget_mcrs(const unsigned int index);
  /** \internal
   * Swap the i:th and j:th stored automorphism information;
   * i and j must be "in window, i.e. in [long_prune_begin,long_prune_end[
   */
  void long_prune_swap(const unsigned int i, const unsigned int j);

  /*
   * Data structures and routines for refining the partition p into equitable
   */
  Heap neighbour_heap;
  virtual bool split_neighbourhood_of_unit_cell(Partition::Cell * const) = 0;
  virtual bool split_neighbourhood_of_cell(Partition::Cell * const) = 0;
  void refine_to_equitable();
  void refine_to_equitable(Partition::Cell * const unit_cell);
  void refine_to_equitable(Partition::Cell * const unit_cell1,
                           Partition::Cell * const unit_cell2);


  /** \internal
   * \return false if it was detected that the current certificate
   *         is different from the first and/or best (whether this is checked
   *         depends on in_search and refine_compare_certificate flags.
   */
  bool do_refine_to_equitable();

  unsigned int eqref_max_certificate_index;
  /** \internal
   * Whether eqref_hash is updated during equitable refinement process.
   */
  bool compute_eqref_hash;
  UintSeqHash eqref_hash;


  /** \internal
   * Check whether the current partition p is equitable.
   * Performance: very slow, use only for debugging purposes.
   */
  virtual bool is_equitable() const = 0;

  unsigned int *first_path_labeling;
  unsigned int *first_path_labeling_inv;
  Orbit         first_path_orbits;
  unsigned int *first_path_automorphism;

  unsigned int *best_path_labeling;
  unsigned int *best_path_labeling_inv;
  Orbit         best_path_orbits;
  unsigned int *best_path_automorphism;

  void update_labeling(unsigned int * const lab);
  void update_labeling_and_its_inverse(unsigned int * const lab,
                                       unsigned int * const lab_inv);
  void update_orbit_information(Orbit &o, const unsigned int *perm);

  void reset_permutation(unsigned int *perm);

  std::vector<unsigned int> certificate_current_path;
  std::vector<unsigned int> certificate_first_path;
  std::vector<unsigned int> certificate_best_path;

  unsigned int certificate_index;
  virtual void initialize_certificate() = 0;

  /* Remove duplicates from seq.
   * If m is the largest element in seq, them m < tmp.size() must hold.
   * All entries in tmp must be false when called.
   * Under that condition, all entries in tmp are false on exit as well.
   */
  static void remove_duplicates(std::vector<unsigned int>& seq, std::vector<bool>& tmp);

  virtual void remove_duplicate_edges() = 0;
  virtual void make_initial_equitable_partition() = 0;
  virtual Partition::Cell* find_next_cell_to_be_splitted(Partition::Cell *cell) = 0;


  /** \struct PathInfo
   *
   * A structure for holding first, current, and best path information.
   */
  typedef struct {
    unsigned int splitting_element;
    unsigned int certificate_index;
    unsigned int subcertificate_length;
    UintSeqHash eqref_hash;
  } PathInfo;

  void search(const bool canonical, Stats &stats,
              const std::function<void(unsigned int n, const unsigned int* aut)>& report_function = nullptr,
              const std::function<bool()>& terminate = nullptr);


  /*
   *
   * Nonuniform component recursion (NUCR)
   *
   */

  /* The currently traversed component */
  unsigned int cr_level;

  /** @internal @class CR_CEP
   * The "Component End Point" data structure
   */
  class CR_CEP {
  public:
    /** At which level in the search was this CEP created */
    unsigned int creation_level;
    /** The current component has been fully traversed when the partition has
     * this many discrete cells left */
    unsigned int discrete_cell_limit;
    /** The component to be traversed after the current one */
    unsigned int next_cr_level;
    /** The next component end point */
    unsigned int next_cep_index;
    bool first_checked;
    bool best_checked;
  };
  /** \internal
   * A stack for storing Component End Points
   */
  std::vector<CR_CEP> cr_cep_stack;

  /** \internal
   * Find the first non-uniformity component at the component recursion
   * level \a level.
   * The component is stored in \a cr_component.
   * If no component is found, \a cr_component is empty.
   * Returns false if all the cells in the component recursion level \a level
   * were discrete.
   * Modifies the max_ival and max_ival_count fields of Partition:Cell
   * (assumes that they are 0 when called and
   *  quarantees that they are 0 when returned).
   */
  virtual bool nucr_find_first_component(const unsigned int level) = 0;
  virtual bool nucr_find_first_component(const unsigned int level,
                                         std::vector<unsigned int>& component,
                                         unsigned int& component_elements,
                                         Partition::Cell*& sh_return) = 0;
  /** \internal
   * The non-uniformity component found by nucr_find_first_component()
   * is stored here.
   */
  std::vector<unsigned int> cr_component;
  /** \internal
   * The number of vertices in the component \a cr_component
   */
  unsigned int cr_component_elements;






};

} // namespace bliss
