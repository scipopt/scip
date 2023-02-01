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

#include "bliss/abstractgraph.hh"

namespace bliss {

/**
 * \brief The class for undirected, vertex colored graphs.
 *
 * Multiple edges between vertices are not allowed (i.e., are ignored).
 */
class Graph : public AbstractGraph
{
public:
  /**
   * The possible splitting heuristics.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  typedef enum {
    /** First non-unit cell.
     * Very fast but may result in large search spaces on difficult graphs.
     * Use for large but easy graphs. */
    shs_f = 0,
    /** First smallest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fs,
    /** First largest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fl,
    /** First maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fm,
    /** First smallest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fsm,
    /** First largest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_flm
  } SplittingHeuristic;

protected:
  class Vertex {
  public:
    Vertex();
    ~Vertex();
    void add_edge(const unsigned int other_vertex);
    void remove_duplicate_edges(std::vector<bool>& tmp);
    void sort_edges();

    unsigned int color;
    std::vector<unsigned int> edges;
    unsigned int nof_edges() const {return edges.size(); }
  };
  std::vector<Vertex> vertices;
  void sort_edges();
  void remove_duplicate_edges();

  /** \internal
   * Partition independent invariant.
   * Returns the color of the vertex.
   * Time complexity: O(1).
   */
  static unsigned int vertex_color_invariant(const Graph* const g,
                                             const unsigned int v);

  /** \internal
   * Partition independent invariant.
   * Returns the degree of the vertex.
   * DUPLICATE EDGES MUST HAVE BEEN REMOVED BEFORE.
   * Time complexity: O(1).
   */
  static unsigned int degree_invariant(const Graph* const g,
                                       const unsigned int v);

  /** \internal
   * Partition independent invariant.
   * Returns 1 if there is an edge from the vertex to itself, 0 if not.
   * Time complexity: O(k), where k is the number of edges leaving the vertex.
   */
  static unsigned int selfloop_invariant(const Graph* const g,
                                         const unsigned int v);


  bool refine_according_to_invariant(unsigned int (*inv)(const Graph* const g,
                                                         const unsigned int v));

  /*
   * Routines needed when refining the partition p into equitable
   */
  bool split_neighbourhood_of_unit_cell(Partition::Cell * const);
  bool split_neighbourhood_of_cell(Partition::Cell * const);

  /** \internal
   * \copydoc AbstractGraph::is_equitable() const
   */
  bool is_equitable() const;

  /* Splitting heuristics, documented in more detail in graph.cc */
  SplittingHeuristic sh;
  Partition::Cell* find_next_cell_to_be_splitted(Partition::Cell *cell);
  Partition::Cell* sh_first();
  Partition::Cell* sh_first_smallest();
  Partition::Cell* sh_first_largest();
  Partition::Cell* sh_first_max_neighbours();
  Partition::Cell* sh_first_smallest_max_neighbours();
  Partition::Cell* sh_first_largest_max_neighbours();


  /* A data structure used in many functions.
   * Allocated only once to reduce allocation overhead,
   * may be used only in one function at a time.
   */
  std::vector<Partition::Cell*> _neighbour_cells;

  void make_initial_equitable_partition();

  void initialize_certificate();



  bool nucr_find_first_component(const unsigned int level);
  bool nucr_find_first_component(const unsigned int level,
                                 std::vector<unsigned int>& component,
                                 unsigned int& component_elements,
                                 Partition::Cell*& sh_return);




public:
  /**
   * Create a new graph with \a N vertices and no edges.
   */
  Graph(const unsigned int N = 0);

  /**
   * Destroy the graph.
   */
  ~Graph();

  /**
   * Read the graph from the file \a fp in a variant of the DIMACS format.
   * See the <a href="https://users.aalto.fi/tjunttil/bliss">bliss website</a>
   * for the definition of the file format.
   * Note that in the DIMACS file the vertices are numbered from 1 to N while
   * in this C++ API they are from 0 to N-1.
   * Thus the vertex n in the file corresponds to the vertex n-1 in the API.
   *
   * \param fp      the file stream for the graph file
   * \param errstr  if non-null, the possible error messages are printed
   *                in this file stream
   * \return        a new Graph object or 0 if reading failed for some
   *                reason
   */
  static Graph* read_dimacs(FILE* const fp, FILE* const errstr = stderr);

  /**
   * Write the graph to a file in a variant of the DIMACS format.
   * See the <a href="https://users.aalto.fi/tjunttil/bliss">bliss website</a>
   * for the definition of the file format.
   */
  void write_dimacs(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(FILE * const fp)
   */
  void write_dot(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(const char * const file_name)
   */
  void write_dot(const char* const file_name);



  /**
   * \copydoc AbstractGraph::get_hash()
   */
  virtual unsigned int get_hash();

  /**
   * Return the number of vertices in the graph.
   */
  unsigned int get_nof_vertices() const {return vertices.size(); }

  /**
   * \copydoc AbstractGraph::permute(const unsigned int* const perm) const
   */
  Graph* permute(const unsigned int* const perm) const;
  /**
   * \copydoc AbstractGraph::permute(const std::vector<unsigned int>& perm) const
   */
  Graph* permute(const std::vector<unsigned int>& perm) const;

  /**
   * \copydoc AbstractGraph::is_automorphism(unsigned int* const perm) const
   */
  bool is_automorphism(unsigned int* const perm) const;

  /**
   * \copydoc AbstractGraph::is_automorphism(const std::vector<unsigned int>& perm) const
   */
  bool is_automorphism(const std::vector<unsigned int>& perm) const;

  /**
   * Add a new vertex with color \a color in the graph and return its index.
   */
  unsigned int add_vertex(const unsigned int color = 0);

  /**
   * Add an edge between vertices \a v1 and \a v2.
   * Duplicate edges between vertices are ignored but try to avoid introducing
   * them in the first place as they are not ignored immediately but will
   * consume memory and computation resources for a while.
   */
  void add_edge(const unsigned int v1, const unsigned int v2);

  /**
   * \copydoc AbstractGraph::get_color(const unsigned int vertex) const
   */
  unsigned int get_color(const unsigned int vertex) const;

  /**
   * Change the color of the vertex \a vertex to \a color.
   */
  void change_color(const unsigned int vertex, const unsigned int color);

  /**
   * Get a copy of the graph.
   */
  Graph* copy() const;

  /**
   * Compare this graph to the \a other graph in a total orger on graphs.
   * \return 0 if the graphs are equal,
   *         -1 if this graph is "smaller than" the other, and
   *         1 if this graph is "greater than" the other.
   */
  int cmp(Graph& other);

  /**
   * Set the splitting heuristic used by the automorphism and canonical
   * labeling algorithm.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  void set_splitting_heuristic(const SplittingHeuristic shs) {sh = shs; }


};

} // namespace bliss
