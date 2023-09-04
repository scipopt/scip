#include <unordered_set>
#include <string>
#include "graph.hh"  // from bliss

/** \file
 * A small example that enumerates all non-isomorphic graphs of n vertices.
 * Uses the "folklore" method.
 */

/* The equals operator needed in std::unordered_set of BlissGraph*. */
struct equals {
  bool operator()(bliss::Graph* g1, bliss::Graph* g2) const {
    return g1->cmp(*g2) == 0;
  }
};

/* The hash operator needed in std::unordered_set of BlissGraph*. */
struct hash {
  unsigned int operator()(bliss::Graph* g) const {
    return g->get_hash();
  }
};

void
traverse(const int n, std::unordered_set<bliss::Graph*, hash, equals>& seen, bliss::Graph* g, unsigned int& nof_graphs) {
  // Build the canonical form of g
  bliss::Stats stats;
  bliss::Graph* g_canonical = g->permute(g->canonical_form(stats));
  // Do not expand this node if an isomorphic graph have already been expanded
  if(seen.find(g_canonical) != seen.end()) {
    delete g_canonical;
    return;
  }
  seen.insert(g_canonical);

  const unsigned int k = g->get_nof_vertices();
  if(k == n) {
    // A graph with n vertices has been generated
    nof_graphs++;
    return;
  }
  // Expand children (graphs with one more vertex)
  for(unsigned long i = 0; i < (1L << k); i++) {
    bliss::Graph* child = g->copy();
    int v = child->add_vertex();
    for(int j = 0; j < n; j++)
      if((i >> j) & 0x01)
	child->add_edge(j, v);
    traverse(n, seen, child, nof_graphs);
    delete child;
  }
}


int main(const int argc, const char** argv)
{
  int n = 7;
  if(argc >= 2)
    n = std::stoi(std::string(argv[1]));
  bliss::Graph* root = new bliss::Graph();
  std::unordered_set<bliss::Graph*, hash, equals> seen;
  unsigned int nof_graphs = 0;
  traverse(n, seen, root, nof_graphs);
  printf("There are %u non-isomorphic simple graphs of %d vertices\n", nof_graphs, n);
  // Clear everything
  for(auto g: seen)
    delete g;
  seen.clear();
  delete root;
  return 0;
}
