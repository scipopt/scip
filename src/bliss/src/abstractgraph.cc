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

#include <new>
#include <set>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cassert>
#include <climits>

#include "bliss/defs.hh"
#include "bliss/timer.hh"
#include "bliss/abstractgraph.hh"
#include "bliss/partition.hh"
#include "bliss/utils.hh"

namespace bliss {


/*-------------------------------------------------------------------------
 *
 * Constructor and destructor routines for the abstract graph class
 *
 *-------------------------------------------------------------------------*/


AbstractGraph::AbstractGraph()
{
  /* Initialize stuff */
  first_path_labeling = nullptr;
  first_path_labeling_inv = nullptr;
  best_path_labeling = nullptr;
  best_path_labeling_inv = nullptr;
  first_path_automorphism = nullptr;
  best_path_automorphism = nullptr;
  in_search = false;

  /* Default value for using "long prune" */
  opt_use_long_prune = true;
  /* Default value for using failure recording */
  opt_use_failure_recording = true;
  /* Default value for using component recursion */
  opt_use_comprec = true;


  verbose_level = 0;
  verbstr = stdout;
}


AbstractGraph::~AbstractGraph()
{
  delete[] first_path_labeling; first_path_labeling = nullptr;
  delete[] first_path_labeling_inv; first_path_labeling_inv = nullptr;
  delete[] first_path_automorphism; first_path_automorphism = nullptr;

  delete[] best_path_labeling; best_path_labeling = nullptr;
  delete[] best_path_labeling_inv; best_path_labeling_inv = nullptr;
  delete[] best_path_automorphism; best_path_automorphism = nullptr;
}



/*-------------------------------------------------------------------------
 *
 * Verbose output management routines
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::set_verbose_level(const unsigned int level)
{
  verbose_level = level;
}

void
AbstractGraph::set_verbose_file(FILE* const fp)
{
  verbstr = fp;
}


/*-------------------------------------------------------------------------
 *
 * Helpers
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::remove_duplicates(std::vector<unsigned int>& seq, std::vector<bool>& tmp) {
#if defined(BLISS_CONSISTENCY_CHECKS)
  /* Pre-conditions  */
  for(unsigned int i = 0; i < tmp.size(); i++) assert(tmp[i] == false);
#endif
  const unsigned int m = seq.size();
  unsigned int j = 0;
  for(unsigned int i = 0; i < m; i++) {
    if(tmp[seq[i]] == false) {
      // Not yet seen
      tmp[seq[i]] = true;
      seq[j++] = seq[i];
    }
  }
  seq.resize(j);
  /* Clear tmp */
  for(int e: seq) tmp[e] = false;
#if defined(BLISS_CONSISTENCY_CHECKS)
  /* Post-conditions  */
  for(unsigned int i = 0; i < tmp.size(); i++) assert(tmp[i] == false);
#endif
}

/*-------------------------------------------------------------------------
 *
 * Routines for refinement to equitable partition
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::refine_to_equitable()
{

  /* Start refinement from all cells -> push 'em all in the splitting queue */
  for(Partition::Cell* cell = p.first_cell; cell; cell = cell->next)
    p.splitting_queue_add(cell);

  do_refine_to_equitable();

}

void
AbstractGraph::refine_to_equitable(Partition::Cell* const unit_cell)
{

  p.splitting_queue_add(unit_cell);

  do_refine_to_equitable();
}



void
AbstractGraph::refine_to_equitable(Partition::Cell* const unit_cell1,
                                   Partition::Cell* const unit_cell2)
{

  p.splitting_queue_add(unit_cell1);
  p.splitting_queue_add(unit_cell2);

  do_refine_to_equitable();
}



bool
AbstractGraph::do_refine_to_equitable()
{

  eqref_hash.reset();

  while(!p.splitting_queue_is_empty())
    {
      Partition::Cell* const cell = p.splitting_queue_pop();

      if(cell->is_unit())
        {
          if(in_search) {
            const unsigned int index = cell->first;
            if(first_path_automorphism)
              {
                /* Build the (potential) automorphism on-the-fly */
                first_path_automorphism[first_path_labeling_inv[index]] =
                  p.elements[index];
              }
            if(best_path_automorphism)
              {
                /* Build the (potential) automorphism on-the-fly */
                best_path_automorphism[best_path_labeling_inv[index]] =
                  p.elements[index];
              }
          }
          const bool worse = split_neighbourhood_of_unit_cell(cell);
          if(in_search and worse)
            goto worse_exit;
        }
      else
        {
          const bool worse = split_neighbourhood_of_cell(cell);
          if(in_search and worse)
            goto worse_exit;
        }
    }

  return true;

 worse_exit:
  /* Clear splitting_queue */
  p.splitting_queue_clear();
  return false;
}
















/*-------------------------------------------------------------------------
 *
 * Routines for handling the canonical labeling
 *
 *-------------------------------------------------------------------------*/

/** \internal
 * Assign the labeling induced by the current partition 'this.p' to
 * \a labeling.
 * That is, if the partition is [[2,0],[1]],
 * then \a labeling will map 0 to 1, 1 to 2, and 2 to 0.
 */
void
AbstractGraph::update_labeling(unsigned int* const labeling)
{
  const unsigned int N = get_nof_vertices();
  unsigned int* ep = p.elements;
  for(unsigned int i = 0; i < N; i++, ep++)
    labeling[*ep] = i;
}

/** \internal
 * The same as update_labeling() except that the inverse of the labeling
 * is also produced and assigned to \a labeling_inv.
 */
void
AbstractGraph::update_labeling_and_its_inverse(unsigned int* const labeling,
                                               unsigned int* const labeling_inv)
{
  const unsigned int N = get_nof_vertices();
  unsigned int* ep = p.elements;
  unsigned int* clip = labeling_inv;

  for(unsigned int i = 0; i < N; ) {
    labeling[*ep] = i;
    i++;
    *clip = *ep;
    ep++;
    clip++;
  }
}





/*-------------------------------------------------------------------------
 *
 * Routines for handling automorphisms
 *
 *-------------------------------------------------------------------------*/


/** \internal
 * Reset the permutation \a perm to the identity permutation.
 */
void
AbstractGraph::reset_permutation(unsigned int* perm)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++, perm++)
    *perm = i;
}



/*-------------------------------------------------------------------------
 *
 * Certificate building
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::cert_add(const unsigned int v1,
                        const unsigned int v2,
                        const unsigned int v3)
{
  if(refine_compare_certificate)
    {
      if(refine_equal_to_first)
        {
          /* So far equivalent to the first path... */
          unsigned int index = certificate_current_path.size();
          if(index >= refine_first_path_subcertificate_end)
            {
              refine_equal_to_first = false;
            }
          else if(certificate_first_path[index] != v1)
            {
              refine_equal_to_first = false;
            }
          else if(certificate_first_path[++index] != v2)
            {
              refine_equal_to_first = false;
            }
          else if(certificate_first_path[++index] != v3)
            {
              refine_equal_to_first = false;
            }
          if(opt_use_failure_recording and !refine_equal_to_first)
            {
              /* We just became different from the first path,
               * remember the deviation point tree-specific invariant
               * for the use of failure recording */
              UintSeqHash h;
              h.update(v1);
              h.update(v2);
              h.update(v3);
              h.update(index);
              h.update(eqref_hash.get_value());
              failure_recording_fp_deviation = h.get_value();
            }
        }
      if(refine_cmp_to_best == 0)
        {
          /* So far equivalent to the current best path... */
          unsigned int index = certificate_current_path.size();
          if(index >= refine_best_path_subcertificate_end)
            {
              refine_cmp_to_best = 1;
            }
          else if(v1 > certificate_best_path[index])
            {
              refine_cmp_to_best = 1;
            }
          else if(v1 < certificate_best_path[index])
            {
              refine_cmp_to_best = -1;
            }
          else if(v2 > certificate_best_path[++index])
            {
              refine_cmp_to_best = 1;
            }
          else if(v2 < certificate_best_path[index])
            {
              refine_cmp_to_best = -1;
            }
          else if(v3 > certificate_best_path[++index])
            {
              refine_cmp_to_best = 1;
            }
          else if(v3 < certificate_best_path[index])
            {
              refine_cmp_to_best = -1;
            }
        }
      if((refine_equal_to_first == false) and
         (refine_cmp_to_best < 0))
        return;
    }
  /* Update the current path certificate */
  certificate_current_path.push_back(v1);
  certificate_current_path.push_back(v2);
  certificate_current_path.push_back(v3);
}


void
AbstractGraph::cert_add_redundant(const unsigned int v1,
                                  const unsigned int v2,
                                  const unsigned int v3)
{
  return cert_add(v1, v2, v3);
}











/*-------------------------------------------------------------------------
 *
 * Long prune code
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::long_prune_init()
{
  const unsigned int N = get_nof_vertices();
  long_prune_temp.clear();
  long_prune_temp.resize(N);
  /* Of how many automorphisms we can store information in
     the predefined, fixed amount of memory? */
  const unsigned int nof_fitting_in_max_mem =
    (long_prune_options_max_mem * 1024 * 1024) / (((N * 2) / 8)+1);
  long_prune_max_stored_autss = long_prune_options_max_stored_auts;
  /* Had some problems with g++ in using (a<b)?a:b when constants involved,
     so had to make this in a stupid way... */
  if(nof_fitting_in_max_mem < long_prune_options_max_stored_auts)
    long_prune_max_stored_autss = nof_fitting_in_max_mem;

  long_prune_deallocate();
  long_prune_fixed.resize(N, 0);
  long_prune_mcrs.resize(N, 0);
  long_prune_begin = 0;
  long_prune_end = 0;
}

void
AbstractGraph::long_prune_deallocate()
{
  while(!long_prune_fixed.empty())
    {
      delete long_prune_fixed.back();
      long_prune_fixed.pop_back();
    }
  while(!long_prune_mcrs.empty())
    {
      delete long_prune_mcrs.back();
      long_prune_mcrs.pop_back();
    }
}

void
AbstractGraph::long_prune_swap(const unsigned int i, const unsigned int j)
{
  const unsigned int real_i = i % long_prune_max_stored_autss;
  const unsigned int real_j = j % long_prune_max_stored_autss;
  std::vector<bool>* tmp = long_prune_fixed[real_i];
  long_prune_fixed[real_i] = long_prune_fixed[real_j];
  long_prune_fixed[real_j] = tmp;
  tmp = long_prune_mcrs[real_i];
  long_prune_mcrs[real_i] = long_prune_mcrs[real_j];
  long_prune_mcrs[real_j] = tmp;
}

std::vector<bool>&
AbstractGraph::long_prune_allocget_fixed(const unsigned int index)
{
  const unsigned int i = index % long_prune_max_stored_autss;
  if(!long_prune_fixed[i])
    long_prune_fixed[i] = new std::vector<bool>(get_nof_vertices());
  return *long_prune_fixed[i];
}

std::vector<bool>&
AbstractGraph::long_prune_get_fixed(const unsigned int index)
{
  return *long_prune_fixed[index % long_prune_max_stored_autss];
}

std::vector<bool>&
AbstractGraph::long_prune_allocget_mcrs(const unsigned int index)
{
  const unsigned int i = index % long_prune_max_stored_autss;
  if(!long_prune_mcrs[i])
    long_prune_mcrs[i] = new std::vector<bool>(get_nof_vertices());
  return *long_prune_mcrs[i];
}

std::vector<bool>&
AbstractGraph::long_prune_get_mcrs(const unsigned int index)
{
  return *long_prune_mcrs[index % long_prune_max_stored_autss];
}

void
AbstractGraph::long_prune_add_automorphism(const unsigned int* aut)
{
  if(long_prune_max_stored_autss == 0)
    return;

  const unsigned int N = get_nof_vertices();


  /* If the buffer of stored auts is full, remove the oldest aut */
  if(long_prune_end - long_prune_begin == long_prune_max_stored_autss)
    {
      long_prune_begin++;
    }
  long_prune_end++;
  std::vector<bool>& fixed = long_prune_allocget_fixed(long_prune_end-1);
  std::vector<bool>& mcrs = long_prune_allocget_mcrs(long_prune_end-1);
  /* Mark nodes that are (i) fixed or (ii) minimal orbit representatives
   * under the automorphism 'aut' */
  for(unsigned int i = 0; i < N; i++)
    {
      fixed[i] = (aut[i] == i);
      if(long_prune_temp[i] == false)
        {
          mcrs[i] = true;
          unsigned int j = aut[i];
          while(j != i)
            {
              long_prune_temp[j] = true;
              j = aut[j];
            }
        }
      else
        {
          mcrs[i] = false;
        }
      /* Clear the temp array on-the-fly... */
      long_prune_temp[i] = false;
    }


}



/*-------------------------------------------------------------------------
 *
 * Routines for handling orbit information
 *
 *-------------------------------------------------------------------------*/

void
AbstractGraph::update_orbit_information(Orbit& o, const unsigned int* perm)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++)
    if(perm[i] != i)
      o.merge_orbits(i, perm[i]);
}








/*-------------------------------------------------------------------------
 *
 * The actual backtracking search
 *
 *-------------------------------------------------------------------------*/

/** \internal \brief Search tree node information.
 */
class TreeNode
{
  //friend class AbstractGraph;
public:
  unsigned int split_cell_first;

  int split_element;
  static const int SPLIT_START = -1;
  static const int SPLIT_END   = -2;

  Partition::BacktrackPoint partition_bt_point;

  unsigned int certificate_index;

  static const char NO = -1;
  static const char MAYBE = 0;
  static const char YES = 1;

  /* First path stuff */
  bool fp_on;
  bool fp_cert_equal;
  char fp_extendable;

  /* Best path stuff */
  bool in_best_path;
  int cmp_to_best_path;

  unsigned int failure_recording_ival;

  /* Component recursion related data */
  unsigned int cr_cep_stack_size;
  unsigned int cr_cep_index;
  unsigned int cr_level;

  bool needs_long_prune;
  unsigned int long_prune_begin;
  std::set<unsigned int, std::less<unsigned int> > long_prune_redundant;

  UintSeqHash eqref_hash;
  unsigned int subcertificate_length;
};





void
AbstractGraph::search(const bool canonical,
                      Stats& stats,
                      const std::function<void(unsigned int n, const unsigned int* aut)>& report,
                      const std::function<bool()>& terminate)
{
  const unsigned int N = get_nof_vertices();

  unsigned int all_same_level = UINT_MAX;

  p.graph = this;

  /*
   * Must be done!
   */
  remove_duplicate_edges();

  /*
   * Reset search statistics
   */
  stats.reset();
  stats.nof_nodes = 1;
  stats.nof_leaf_nodes = 1;

  /* Free old first path data structures */
  delete[] first_path_labeling; first_path_labeling = nullptr;
  delete[] first_path_labeling_inv; first_path_labeling_inv = nullptr;
  delete[] first_path_automorphism; first_path_automorphism = nullptr;

  /* Free old best path data structures */
  delete[] best_path_labeling; best_path_labeling = nullptr;
  delete[] best_path_labeling_inv; best_path_labeling_inv = nullptr;
  delete[] best_path_automorphism; best_path_automorphism = nullptr;

  if(N == 0)
    {
      /* Nothing to do, return... */
      return;
    }

  /* Initialize the partition ... */
  p.init(N);
  /* ... and the component recursion data structures in the partition */
  if(opt_use_comprec)
    p.cr_init();

  // REMOVE
  // neighbour_heap.init(N);

  in_search = false;
  /* Do not compute certificate when building the initial partition */
  refine_compare_certificate = false;
  refine_equal_to_first = false;
  /* The 'eqref_hash' hash value is not computed when building
   * the initial partition as it is not used for anything at the moment.
   * This saves some cycles. */
  compute_eqref_hash = false;

  Timer timer1;

  make_initial_equitable_partition();

#if defined(BLISS_VERIFY_EQUITABLEDNESS)
  if(!is_equitable())
    throw std::logic_error("consistency check failed - initial partition not equitable");
#endif

  if(verbstr and verbose_level >= 2)
    {
      fprintf(verbstr, "Initial partition computed in %.2f seconds\n",
              timer1.get_duration());
      fflush(verbstr);
    }



  /*
   * Allocate space for the "first path" and "best path" labelings
   */
  delete[] first_path_labeling;
  first_path_labeling = new unsigned int[N];

  delete[] best_path_labeling;
  best_path_labeling = new unsigned int[N];
  for(unsigned int i = 0; i < N; i++) best_path_labeling[i] = i;

  /*
   * Is the initial partition discrete?
   */
  if(p.is_discrete())
    {
      /* Make the best path labeling i.e. the canonical labeling */
      update_labeling(best_path_labeling);
      /* Update statistics */
      stats.nof_leaf_nodes = 1;
      /* Release component recursion data in partition */
      if(opt_use_comprec)
        p.cr_free();
      return;
    }

  /*
   * Allocate the inverses of the "first path" and "best path" labelings
   */
  delete[] first_path_labeling_inv;
  first_path_labeling_inv = new unsigned int[N];
  std::fill_n(first_path_labeling_inv, N, 0);
  delete[] best_path_labeling_inv;
  best_path_labeling_inv = new unsigned int[N];
  std::fill_n(best_path_labeling_inv, N, 0);

  /*
   * Allocate space for the automorphisms
   */
  delete[] first_path_automorphism;
  first_path_automorphism = new unsigned int[N];
  delete[] best_path_automorphism;
  best_path_automorphism = new unsigned int[N];

  /*
   * Initialize orbit information so that all vertices are in their own orbits
   */
  first_path_orbits.init(N);
  best_path_orbits.init(N);

  /*
   * Initialize certificate memory
   */
  initialize_certificate();

  std::vector<TreeNode> search_stack;
  std::vector<PathInfo> first_path_info;
  std::vector<PathInfo> best_path_info;

  search_stack.clear();

  /* Initialize "long prune" data structures */
  if(opt_use_long_prune)
    long_prune_init();

  /*
   * Initialize failure recording data structures
   */
  typedef std::set<unsigned int, std::less<unsigned int> > FailureRecordingSet;
  std::vector<FailureRecordingSet> failure_recording_hashes;

  /*
   * Initialize component recursion data structures
   */
  cr_cep_stack.clear();
  unsigned int cr_cep_index = 0;
  {
    /* Inset a sentinel "component end point" */
    CR_CEP sentinel;
    sentinel.creation_level = 0;
    sentinel.discrete_cell_limit = get_nof_vertices();
    sentinel.next_cr_level = 0;
    sentinel.next_cep_index = 0;
    sentinel.first_checked = false;
    sentinel.best_checked = false;
    cr_cep_index = 0;
    cr_cep_stack.push_back(sentinel);
  }
  cr_level = 0;
  if(opt_use_comprec and
     nucr_find_first_component(cr_level) == true and
     p.nof_discrete_cells() + cr_component_elements <
     cr_cep_stack[cr_cep_index].discrete_cell_limit)
    {
      cr_level = p.cr_split_level(0, cr_component);
      CR_CEP cep;
      cep.creation_level = 0;
      cep.discrete_cell_limit = p.nof_discrete_cells() + cr_component_elements;
      cep.next_cr_level = 0;
      cep.next_cep_index = cr_cep_index;
      cep.first_checked = false;
      cep.best_checked = false;
      cr_cep_index = cr_cep_stack.size();
      cr_cep_stack.push_back(cep);
    }

  /*
   * Build the root node of the search tree
   */
  {
    TreeNode root;
    Partition::Cell* split_cell = find_next_cell_to_be_splitted(p.first_cell);
    root.split_cell_first = split_cell->first;
    root.split_element = TreeNode::SPLIT_START;
    root.partition_bt_point = p.set_backtrack_point();
    root.certificate_index = 0;
    root.fp_on = true;
    root.fp_cert_equal = true;
    root.fp_extendable = TreeNode::MAYBE;
    root.in_best_path = false;
    root.cmp_to_best_path = 0;
    root.long_prune_begin = 0;
    root.needs_long_prune = true;

    root.failure_recording_ival = 0;

    /* Save component recursion info for backtracking */
    root.cr_level = cr_level;
    root.cr_cep_stack_size = cr_cep_stack.size();
    root.cr_cep_index = cr_cep_index;
    search_stack.push_back(root);
  }

  /*
   * Set status and global flags for search related procedures
   */
  in_search = true;
  /* Do not compare certificates during refinement until the first path has been traversed to the leaf */
  refine_compare_certificate = false;




  /*
   * The actual backtracking search
   */
  while(!search_stack.empty())
    {
      if(terminate and terminate()) {
        break;
      }
      TreeNode&          current_node  = search_stack.back();
      const unsigned int current_level = (unsigned int)search_stack.size()-1;


      if(opt_use_comprec)
        {
          CR_CEP& cep = cr_cep_stack[current_node.cr_cep_index];
          if(cep.first_checked == true and
             current_node.fp_extendable == TreeNode::MAYBE and
             !search_stack[cep.creation_level].fp_on)
            {
              current_node.fp_extendable = TreeNode::NO;
            }
        }

      if(current_node.fp_on)
        {
          if(current_node.split_element == TreeNode::SPLIT_END)
            {
              search_stack.pop_back();
              continue;
            }
        }
      else
        {
          if(current_node.fp_extendable == TreeNode::YES)
            {
              search_stack.pop_back();
              continue;
            }
          if(current_node.split_element == TreeNode::SPLIT_END)
            {
              if(opt_use_failure_recording)
                {
                  TreeNode& parent_node = search_stack[current_level-1];
                  if(parent_node.fp_on)
                    failure_recording_hashes[current_level-1].insert(current_node.failure_recording_ival);
                }
              search_stack.pop_back();
              continue;
            }
          if(current_node.fp_extendable == TreeNode::NO and
             (!canonical or current_node.cmp_to_best_path < 0))
            {
              if(opt_use_failure_recording)
                {
                  TreeNode& parent_node = search_stack[current_level-1];
                  if(parent_node.fp_on)
                    failure_recording_hashes[current_level-1].insert(current_node.failure_recording_ival);
                }
              search_stack.pop_back();
              continue;
            }
        }

      /* Restore partition ... */
      p.goto_backtrack_point(current_node.partition_bt_point);
      /* ... and re-remember backtracking point */
      current_node.partition_bt_point = p.set_backtrack_point();

      /* Restore current path certificate */
      certificate_index = current_node.certificate_index;
      refine_current_path_certificate_index = current_node.certificate_index;
      certificate_current_path.resize(certificate_index);

      /* Fetch split cell information */
      Partition::Cell * const cell =
        p.get_cell(p.elements[current_node.split_cell_first]);

      /* Restore component recursion information */
      cr_level = current_node.cr_level;
      cr_cep_stack.resize(current_node.cr_cep_stack_size);
      cr_cep_index = current_node.cr_cep_index;


      /*
       * Update long prune redundancy sets
       */
      if(opt_use_long_prune and current_level >= 1 and !current_node.fp_on)
        {
          unsigned int begin = (current_node.long_prune_begin>long_prune_begin)?current_node.long_prune_begin:long_prune_begin;
          for(unsigned int i = begin; i < long_prune_end; i++)
            {
              const std::vector<bool>& fixed = long_prune_get_fixed(i);
#if defined(BLISS_CONSISTENCY_CHECKS)
              for(unsigned int l = 0; l < search_stack.size()-2; l++)
                assert(fixed[search_stack[l].split_element]);
#endif
              if(fixed[search_stack[search_stack.size()-1-1].split_element] ==
                 false)
                {
                  long_prune_swap(begin, i);
                  begin++;
                  current_node.long_prune_begin = begin;
                  continue;
                }
            }

          if(current_node.split_element == TreeNode::SPLIT_START)
            {
              current_node.needs_long_prune = true;
            }
          else if(current_node.needs_long_prune)
            {
              current_node.needs_long_prune = false;
              begin = (current_node.long_prune_begin>long_prune_begin)?current_node.long_prune_begin:long_prune_begin;
              for(unsigned int i = begin; i < long_prune_end; i++)
                {
                  const std::vector<bool>& fixed = long_prune_get_fixed(i);
#if defined(BLISS_CONSISTENCY_CHECKS)
                  for(unsigned int l = 0; l < search_stack.size()-2; l++)
                    assert(fixed[search_stack[l].split_element]);
#endif
                  assert(fixed[search_stack[current_level-1].split_element] == true);
                  if(fixed[search_stack[current_level-1].split_element] == false)
                    {
                      long_prune_swap(begin, i);
                      begin++;
                      current_node.long_prune_begin = begin;
                      continue;
                    }
                  const std::vector<bool>& mcrs = long_prune_get_mcrs(i);
                  unsigned int* ep = p.elements + cell->first;
                  for(unsigned int j = cell->length; j > 0; j--, ep++) {
                    if(mcrs[*ep] == false)
                      current_node.long_prune_redundant.insert(*ep);
                  }
                }
            }
        }


      /*
       * Find the next smallest, non-isomorphic element in the cell and
       * store it in current_node.split_element
       */
      {
        unsigned int  next_split_element = UINT_MAX;
        //unsigned int* next_split_element_pos = 0;
        unsigned int* ep = p.elements + cell->first;
        if(current_node.fp_on)
          {
            /* Find the next larger splitting element that is
             * a minimal orbit representative w.r.t. first_path_orbits */
            for(unsigned int i = cell->length; i > 0; i--, ep++) {
              if((int)(*ep) > current_node.split_element and
                 *ep < next_split_element and
                 first_path_orbits.is_minimal_representative(*ep)) {
                next_split_element = *ep;
                //next_split_element_pos = ep;
              }
            }
          }
        else if(current_node.in_best_path)
          {
            /* Find the next larger splitting element that is
             * a minimal orbit representative w.r.t. best_path_orbits */
            for(unsigned int i = cell->length; i > 0; i--, ep++) {
              if((int)(*ep) > current_node.split_element and
                 *ep < next_split_element and
                 best_path_orbits.is_minimal_representative(*ep) and
                 (!opt_use_long_prune or
                  current_node.long_prune_redundant.find(*ep) ==
                  current_node.long_prune_redundant.end())) {
                next_split_element = *ep;
                //next_split_element_pos = ep;
              }
            }
          }
        else
          {
            /* Find the next larger splitting element */
            for(unsigned int i = cell->length; i > 0; i--, ep++) {
              if((int)(*ep) > current_node.split_element and
                 *ep < next_split_element and
                 (!opt_use_long_prune or
                  current_node.long_prune_redundant.find(*ep) ==
                  current_node.long_prune_redundant.end())) {
                next_split_element = *ep;
                //next_split_element_pos = ep;
              }
            }
          }
        if(next_split_element == UINT_MAX)
          {
            /* No more (unexplored children) in the cell */
            current_node.split_element = TreeNode::SPLIT_END;
            if(current_node.fp_on)
              {
                /* Update group size */
                const unsigned int index = first_path_orbits.orbit_size(first_path_info[search_stack.size()-1].splitting_element);
                stats.group_size.multiply(index);
                stats.group_size_approx *= (long double)index;
                /*
                 * Update all_same_level
                 */
                if(index == cell->length and all_same_level == current_level+1)
                  all_same_level = current_level;
                if(verbstr and verbose_level >= 2) {
                  fprintf(verbstr,
                          "Level %u: orbits=%u, index=%u/%u, all_same_level=%u\n",
                          current_level,
                          first_path_orbits.nof_orbits(),
                          index, cell->length,
                          all_same_level);
                  fflush(verbstr);
                }
              }
            continue;
          }

        /* Split on smallest */
        current_node.split_element = next_split_element;
      }

      const unsigned int child_level = current_level+1;
      /* Update some statistics */
      stats.nof_nodes++;
      if(search_stack.size() > stats.max_level)
        stats.max_level = search_stack.size();



      /* Set flags and indices for the refiner certificate builder */
      refine_equal_to_first = current_node.fp_cert_equal;
      refine_cmp_to_best = current_node.cmp_to_best_path;
      if(!first_path_info.empty())
        {
          if(refine_equal_to_first)
            refine_first_path_subcertificate_end =
              first_path_info[search_stack.size()-1].certificate_index +
              first_path_info[search_stack.size()-1].subcertificate_length;
          if(canonical)
            {
              if(refine_cmp_to_best == 0)
                refine_best_path_subcertificate_end =
                  best_path_info[search_stack.size()-1].certificate_index +
                  best_path_info[search_stack.size()-1].subcertificate_length;
            }
          else
            refine_cmp_to_best = -1;
        }

      const bool was_fp_cert_equal = current_node.fp_cert_equal;

      /* Individualize, i.e. split the cell in two, the latter new cell
       * will be a unit one containing info.split_element */
      Partition::Cell* const new_cell =
        p.individualize(cell, current_node.split_element);

      /*
       * Refine the new partition to equitable
       */
      if(cell->is_unit())
        refine_to_equitable(cell, new_cell);
      else
        refine_to_equitable(new_cell);




      /* Update statistics */
      if(p.is_discrete())
        stats.nof_leaf_nodes++;


      if(!first_path_info.empty())
        {
          /* We are no longer on the first path */
          const unsigned int subcertificate_length =
            certificate_current_path.size() - certificate_index;
          if(refine_equal_to_first)
            {
              /* Was equal to the first path so far */
              PathInfo& first_pinfo = first_path_info[current_level];
              assert(first_pinfo.certificate_index == certificate_index);
              if(subcertificate_length != first_pinfo.subcertificate_length)
                {
                  refine_equal_to_first = false;
                  if(opt_use_failure_recording)
                    failure_recording_fp_deviation = subcertificate_length;
                }
              else if(first_pinfo.eqref_hash.cmp(eqref_hash) != 0)
                {
                  refine_equal_to_first = false;
                  if(opt_use_failure_recording)
                    failure_recording_fp_deviation = eqref_hash.get_value();
                }
            }
          if(canonical and (refine_cmp_to_best == 0))
            {
              /* Was equal to the best path so far */
              PathInfo& bestp_info = best_path_info[current_level];
              assert(bestp_info.certificate_index == certificate_index);
              if(subcertificate_length < bestp_info.subcertificate_length)
                {
                  refine_cmp_to_best = -1;
                }
              else if(subcertificate_length > bestp_info.subcertificate_length)
                {
                  refine_cmp_to_best = 1;
                }
              else if(bestp_info.eqref_hash.cmp(eqref_hash) > 0)
                {
                  refine_cmp_to_best = -1;
                }
              else if(bestp_info.eqref_hash.cmp(eqref_hash) < 0)
                {
                  refine_cmp_to_best = 1;
                }
            }

          if(opt_use_failure_recording and
             was_fp_cert_equal and
             !refine_equal_to_first)
            {
              UintSeqHash k;
              k.update(failure_recording_fp_deviation);
              k.update(eqref_hash.get_value());
              failure_recording_fp_deviation = k.get_value();

              if(current_node.fp_on)
                failure_recording_hashes[current_level].insert(failure_recording_fp_deviation);
              else
                {
                  for(unsigned int i = current_level; i > 0; i--)
                    {
                      if(search_stack[i].fp_on)
                        break;
                      const FailureRecordingSet& s = failure_recording_hashes[i];
                      if(i == current_level and
                         s.find(failure_recording_fp_deviation) != s.end())
                        break;
                      if(s.find(0) != s.end())
                        break;
                      search_stack[i].fp_extendable = TreeNode::NO;
                    }
                }
            }


          /* Check if no longer equal to the first path and,
           * if canonical labeling is desired, also worse than the
           * current best path */
          if(refine_equal_to_first == false and
             (!canonical or (refine_cmp_to_best < 0)))
            {
              /* Yes, backtrack */
              stats.nof_bad_nodes++;
              if(current_node.fp_cert_equal == true and
                 current_level+1 > all_same_level)
                {
                  assert(all_same_level >= 1);
                  for(unsigned int i = all_same_level;
                      i < search_stack.size();
                      i++)
                    {
                      search_stack[i].fp_extendable = TreeNode::NO;
                    }
                }

              continue;
            }
        }


      /*
       * Next level search tree node info
       */
      TreeNode child_node;

      /* No more in the first path */
      child_node.fp_on = false;
      /* No more in the best path */
      child_node.in_best_path = false;

      child_node.fp_cert_equal = refine_equal_to_first;
      if(current_node.fp_extendable == TreeNode::NO or
         (current_node.fp_extendable == TreeNode::MAYBE and
          child_node.fp_cert_equal == false))
        child_node.fp_extendable = TreeNode::NO;
      else
        child_node.fp_extendable = TreeNode::MAYBE;
      child_node.cmp_to_best_path = refine_cmp_to_best;

      child_node.failure_recording_ival = 0;
      child_node.cr_cep_stack_size = current_node.cr_cep_stack_size;
      child_node.cr_cep_index = current_node.cr_cep_index;
      child_node.cr_level = current_node.cr_level;

      certificate_index = certificate_current_path.size();

      current_node.eqref_hash = eqref_hash;
      current_node.subcertificate_length =
        certificate_index - current_node.certificate_index;


      /*
       * The first encountered leaf node at the end of the "first path"?
       */
      if(p.is_discrete() and first_path_info.empty())
        {
          //fprintf(stdout, "Level %u: FIRST\n", child_level); fflush(stdout);
          stats.nof_canupdates++;
          /*
           * Update labelings and their inverses
           */
          update_labeling_and_its_inverse(first_path_labeling,
                                          first_path_labeling_inv);
          update_labeling_and_its_inverse(best_path_labeling,
                                          best_path_labeling_inv);
          /*
           * Reset automorphism array
           */
          reset_permutation(first_path_automorphism);
          reset_permutation(best_path_automorphism);
          /*
           * Reset orbit information
           */
          first_path_orbits.reset();
          best_path_orbits.reset();
          /*
           * Reset group size
           */
          stats.group_size.assign(1);
          stats.group_size_approx = 1.0;
          /*
           * Reset all_same_level
           */
          all_same_level = child_level;
          /*
           * Mark the current path to be the first and best one and save it
           */
          const unsigned int base_size = search_stack.size();
          best_path_info.clear();
          //fprintf(stdout, " New base is: ");
          for(unsigned int i = 0; i < base_size; i++) {
            search_stack[i].fp_on = true;
            search_stack[i].fp_cert_equal = true;
            search_stack[i].fp_extendable = TreeNode::YES;
            search_stack[i].in_best_path = true;
            search_stack[i].cmp_to_best_path = 0;
            PathInfo path_info;
            path_info.splitting_element = search_stack[i].split_element;
            path_info.certificate_index = search_stack[i].certificate_index;
            path_info.eqref_hash = search_stack[i].eqref_hash;
            path_info.subcertificate_length = search_stack[i].subcertificate_length;
            first_path_info.push_back(path_info);
            best_path_info.push_back(path_info);
            //fprintf(stdout, "%u ", search_stack[i].split_element);
          }
          //fprintf(stdout, "\n"); fflush(stdout);
          /* Copy certificates */
          certificate_first_path = certificate_current_path;
          certificate_best_path = certificate_current_path;

          /* From now on, compare certificates when refining */
          refine_compare_certificate = true;

          if(opt_use_failure_recording)
            failure_recording_hashes.resize(base_size);

          /*
          for(unsigned int j = 0; j < search_stack.size(); j++)
            fprintf(stderr, "%u ", search_stack[j].split_element);
          fprintf(stderr, "\n");
          p.print(stderr); fprintf(stderr, "\n");
          */

          /*
           * Backtrack to the previous level
           */
          continue;
        }


      if(p.is_discrete() and child_node.fp_cert_equal)
        {
          /*
           * A leaf node that is equal to the first one.
           * An automorphism found: aut[i] = elements[first_path_labeling[i]]
           */
          goto handle_first_path_automorphism;
        }


      if(!p.is_discrete())
        {
          Partition::Cell* next_split_cell = 0;
          /*
           * An internal, non-leaf node
           */
          if(opt_use_comprec)
            {
              assert(p.nof_discrete_cells() <=
                     cr_cep_stack[cr_cep_index].discrete_cell_limit);
              assert(cr_level == child_node.cr_level);


              if(p.nof_discrete_cells() ==
                 cr_cep_stack[cr_cep_index].discrete_cell_limit)
                {
                  /* We have reached the end of a component */
                  assert(cr_cep_index != 0);
                  CR_CEP& cep = cr_cep_stack[cr_cep_index];

                  /* First, compare with respect to the first path */
                  if(first_path_info.empty() or child_node.fp_cert_equal) {
                    if(cep.first_checked == false)
                      {
                        /* First time, go to the next component */
                        cep.first_checked = true;
                      }
                    else
                      {
                        assert(!first_path_info.empty());
                        assert(cep.creation_level < search_stack.size());
                        TreeNode& old_info = search_stack[cep.creation_level];
                        /* If the component was found when on the first path,
                         * handle the found automorphism as the other
                         * first path automorphisms */
                        if(old_info.fp_on)
                          goto handle_first_path_automorphism;
                      }
                  }

                  if(canonical and
                     !first_path_info.empty() and
                     child_node.cmp_to_best_path >= 0) {
                    if(cep.best_checked == false)
                      {
                        /* First time, go to the next component */
                        cep.best_checked = true;
                      }
                    else
                      {
                        assert(cep.creation_level < search_stack.size());
                        TreeNode& old_info = search_stack[cep.creation_level];
                        if(child_node.cmp_to_best_path == 0) {
                          /* If the component was found when on the best path,
                           * handle the found automorphism as the other
                           * best path automorphisms */
                          if(old_info.in_best_path)
                            goto handle_best_path_automorphism;
                          /* Otherwise, we do not remember the automorhism as
                           * we didn't memorize the path that was invariant
                           * equal to the best one and passed through the
                           * component.
                           * Thus we can only backtrack to the previous level */
                          child_node.cmp_to_best_path = -1;
                          if(!child_node.fp_cert_equal)
                            {
                              continue;
                            }
                        }
                        else {
                          assert(child_node.cmp_to_best_path > 0);
                          if(old_info.in_best_path)
                            {
                              stats.nof_canupdates++;
                              /*
                               * Update canonical labeling and its inverse
                               */
                              for(unsigned int i = 0; i < N; i++) {
                                if(p.get_cell(p.elements[i])->is_unit()) {
                                  best_path_labeling[p.elements[i]] = i;
                                  best_path_labeling_inv[i] = p.elements[i];
                                }
                              }
                              //update_labeling_and_its_inverse(best_path_labeling, best_path_labeling_inv);
                              /* Reset best path automorphism */
                              reset_permutation(best_path_automorphism);
                              /* Reset best path orbit structure */
                              best_path_orbits.reset();
                              /* Mark to be the best one and save prefix */
                              unsigned int postfix_start = cep.creation_level;
                              assert(postfix_start < best_path_info.size());
                              while(p.get_cell(best_path_info[postfix_start].splitting_element)->is_unit()) {
                                postfix_start++;
                                assert(postfix_start < best_path_info.size());
                              }
                              unsigned int postfix_start_cert = best_path_info[postfix_start].certificate_index;
                              std::vector<PathInfo> best_path_temp = best_path_info;
                              best_path_info.clear();
                              for(unsigned int i = 0; i < search_stack.size(); i++) {
                                TreeNode& ss_info = search_stack[i];
                                PathInfo  bp_info;
                                ss_info.cmp_to_best_path = 0;
                                ss_info.in_best_path = true;
                                bp_info.splitting_element = ss_info.split_element;
                                bp_info.certificate_index = ss_info.certificate_index;
                                bp_info.subcertificate_length = ss_info.subcertificate_length;
                                bp_info.eqref_hash = ss_info.eqref_hash;
                                best_path_info.push_back(bp_info);
                              }
                              /* Copy the postfix of the previous best path */
                              for(unsigned int i = postfix_start;
                                  i < best_path_temp.size();
                                  i++)
                                {
                                  best_path_info.push_back(best_path_temp[i]);
                                  best_path_info[best_path_info.size()-1].certificate_index =
                                    best_path_info[best_path_info.size()-2].certificate_index +
                                    best_path_info[best_path_info.size()-2].subcertificate_length;
                                }
                              std::vector<unsigned int> certificate_best_path_old = certificate_best_path;
                              certificate_best_path = certificate_current_path;
                              for(unsigned int i = postfix_start_cert;  i < certificate_best_path_old.size(); i++)
                                certificate_best_path.push_back(certificate_best_path_old[i]);
                              assert(certificate_best_path.size() == best_path_info.back().certificate_index + best_path_info.back().subcertificate_length);
                              /* Backtrack to the previous level */
                              continue;
                            }
                        }
                      }
                  }

                  /* No backtracking performed, go to next componenet */
                  cr_level = cep.next_cr_level;
                  cr_cep_index = cep.next_cep_index;
                }

              /* Check if the current component has been split into
               * new non-uniformity subcomponents */
              //if(nucr_find_first_component(cr_level) == true and
              // p.nof_discrete_cells() + cr_component_elements <
              // cr_cep_stack[cr_cep_index].discrete_cell_limit)
              if(nucr_find_first_component(cr_level, cr_component,
                                           cr_component_elements,
                                           next_split_cell) == true and
                 p.nof_discrete_cells() + cr_component_elements <
                 cr_cep_stack[cr_cep_index].discrete_cell_limit)
                {
                  const unsigned int next_cr_level =
                    p.cr_split_level(cr_level, cr_component);
                  CR_CEP cep;
                  cep.creation_level = search_stack.size();
                  cep.discrete_cell_limit =
                    p.nof_discrete_cells() + cr_component_elements;
                  cep.next_cr_level = cr_level;
                  cep.next_cep_index = cr_cep_index;
                  cep.first_checked = false;
                  cep.best_checked = false;
                  cr_cep_index = cr_cep_stack.size();
                  cr_cep_stack.push_back(cep);
                  cr_level = next_cr_level;
                }
            }


          /*
           * Build the next node info
           */
          /* Find the next cell to be splitted */
          if(!next_split_cell)
            next_split_cell = find_next_cell_to_be_splitted(p.get_cell(p.elements[current_node.split_cell_first]));
          //Partition::Cell * const next_split_cell = find_next_cell_to_be_splitted(p.get_cell(p.elements[current_node.split_cell_first]));
          child_node.split_cell_first = next_split_cell->first;
          child_node.split_element = TreeNode::SPLIT_START;
          child_node.certificate_index = certificate_index;
          child_node.partition_bt_point = p.set_backtrack_point();
          child_node.long_prune_redundant.clear();
          child_node.long_prune_begin = current_node.long_prune_begin;

          /* Save component recursion info for backtracking */
          child_node.cr_level = cr_level;
          child_node.cr_cep_stack_size = cr_cep_stack.size();
          child_node.cr_cep_index = cr_cep_index;

          search_stack.push_back(child_node);
          continue;
        }

      /*
       * A leaf node not in the first path or equivalent to the first path
       */



      if(child_node.cmp_to_best_path > 0)
        {
          /*
           * A new, better representative found
           */
          //fprintf(stdout, "Level %u: NEW BEST\n", child_level); fflush(stdout);
          stats.nof_canupdates++;
          /*
           * Update canonical labeling and its inverse
           */
          update_labeling_and_its_inverse(best_path_labeling,
                                          best_path_labeling_inv);
          /* Reset best path automorphism */
          reset_permutation(best_path_automorphism);
          /* Reset best path orbit structure */
          best_path_orbits.reset();
          /*
           * Mark the current path to be the best one and save it
           */
          const unsigned int base_size = search_stack.size();
          assert(current_level+1 == base_size);
          best_path_info.clear();
          for(unsigned int i = 0; i < base_size; i++) {
            search_stack[i].cmp_to_best_path = 0;
            search_stack[i].in_best_path = true;
            PathInfo path_info;
            path_info.splitting_element = search_stack[i].split_element;
            path_info.certificate_index = search_stack[i].certificate_index;
            path_info.subcertificate_length = search_stack[i].subcertificate_length;
            path_info.eqref_hash = search_stack[i].eqref_hash;
            best_path_info.push_back(path_info);
          }
          certificate_best_path = certificate_current_path;
          /*
           * Backtrack to the previous level
           */
          continue;
        }


    handle_best_path_automorphism:
      {
        /*
         *
         * Best path automorphism handling
         *
         */

        /*
         * Equal to the previous best path
         */
        if(p.is_discrete())
          {
#if defined(BLISS_CONSISTENCY_CHECKS)
            /* Verify that the automorphism is correctly built */
            for(unsigned int i = 0; i < N; i++)
              assert(best_path_automorphism[i] ==
                     p.elements[best_path_labeling[i]]);
#endif
          }
        else
          {
            /* An automorphism that was found before the partition was discrete.
             * Set the image of all elements in non-disrete cells accordingly */
            for(Partition::Cell* c = p.first_nonsingleton_cell; c;
                c = c->next_nonsingleton) {
              for(unsigned int i = c->first; i < c->first+c->length; i++)
                if(p.get_cell(p.elements[best_path_labeling[p.elements[i]]])->is_unit())
                  best_path_automorphism[p.elements[best_path_labeling[p.elements[i]]]] = p.elements[i];
                else
                  best_path_automorphism[p.elements[i]] = p.elements[i];
            }
          }

#if defined(BLISS_VERIFY_AUTOMORPHISMS)
        /* Verify that it really is an automorphism */
        if(!is_automorphism(best_path_automorphism))
          throw std::logic_error("Best path automorhism validation check failed");
#endif

        unsigned int gca_level_with_first = 0;
        for(unsigned int i = search_stack.size(); i > 0; i--) {
          if((int)first_path_info[gca_level_with_first].splitting_element !=
             search_stack[gca_level_with_first].split_element)
            break;
          gca_level_with_first++;
        }

        unsigned int gca_level_with_best = 0;
        for(unsigned int i = search_stack.size(); i > 0; i--) {
          if((int)best_path_info[gca_level_with_best].splitting_element !=
             search_stack[gca_level_with_best].split_element)
            break;
          gca_level_with_best++;
        }

        if(opt_use_long_prune)
          {
            /* Record automorphism */
            long_prune_add_automorphism(best_path_automorphism);
          }

        /*
         * Update orbit information
         */
        update_orbit_information(best_path_orbits, best_path_automorphism);

        /*
         * Update orbit information
         */
        const unsigned int nof_old_orbits = first_path_orbits.nof_orbits();
        update_orbit_information(first_path_orbits, best_path_automorphism);
        if(nof_old_orbits != first_path_orbits.nof_orbits())
          {
            /* Some orbits were merged */
            /* Report automorphism */
            if(report)
              report(get_nof_vertices(), best_path_automorphism);
            /* Update statistics */
            stats.nof_generators++;
          }

        /*
         * Compute backjumping level
         */
        unsigned int backjumping_level;
        if(!first_path_orbits.is_minimal_representative(search_stack[gca_level_with_first].split_element))
          {
            backjumping_level = gca_level_with_first;
          }
        else
          {
            assert(!best_path_orbits.is_minimal_representative(search_stack[gca_level_with_best].split_element));
            backjumping_level = gca_level_with_best;
          }
        /* Backtrack */
        search_stack.resize(backjumping_level + 1);
        continue;
      }

    handle_first_path_automorphism:
      /*
       *
       * A first-path automorphism: aut[i] = elements[first_path_labeling[i]]
       *
       */


      if(p.is_discrete())
        {
#if defined(BLISS_CONSISTENCY_CHECKS)
          /* Verify that the complete automorphism is correctly built */
          for(unsigned int i = 0; i < N; i++)
            assert(first_path_automorphism[i] ==
                   p.elements[first_path_labeling[i]]);
#endif
        }
      else
        {
          /* An automorphism that was found before the partition was discrete.
           * Set the image of all elements in non-disrete cells accordingly */
          for(Partition::Cell* c = p.first_nonsingleton_cell; c;
              c = c->next_nonsingleton) {
            for(unsigned int i = c->first; i < c->first+c->length; i++)
              if(p.get_cell(p.elements[first_path_labeling[p.elements[i]]])->is_unit())
                first_path_automorphism[p.elements[first_path_labeling[p.elements[i]]]] = p.elements[i];
              else
                first_path_automorphism[p.elements[i]] = p.elements[i];
          }
        }

#if defined(BLISS_VERIFY_AUTOMORPHISMS)
      /* Verify that it really is an automorphism */
      if(!is_automorphism(first_path_automorphism))
        throw std::logic_error("First path automorhism validation check failed");
#endif

      if(opt_use_long_prune)
        {
          long_prune_add_automorphism(first_path_automorphism);
        }

      /*
       * Update orbit information
       */
      update_orbit_information(first_path_orbits, first_path_automorphism);

      /*
       * Compute backjumping level
       */
      for(unsigned int i = 0; i < search_stack.size(); i++) {
        TreeNode& n = search_stack[i];
        if(n.fp_on) {
          ;
        } else {
          n.fp_extendable = TreeNode::YES;
        }
      }

      /* Report automorphism */
      if(report)
        report(get_nof_vertices(), first_path_automorphism);
      /* Update statistics */
      stats.nof_generators++;
      continue;

    } /* while(!search_stack.empty()) */




  /* Free "long prune" technique memory */
  if(opt_use_long_prune)
    long_prune_deallocate();

  /* Release component recursion data in partition */
  if(opt_use_comprec)
    p.cr_free();
}




void
AbstractGraph::find_automorphisms(Stats& stats,
                                  const std::function<void(unsigned int n, const unsigned int* aut)>& report,
                                  const std::function<bool()>& terminate)
{
  search(false, stats, report, terminate);

  delete[] first_path_labeling; first_path_labeling = nullptr;
  delete[] best_path_labeling; best_path_labeling = nullptr;
}


const unsigned int *
AbstractGraph::canonical_form(Stats& stats,
                              const std::function<void(unsigned int n, const unsigned int* aut)>& report,
                              const std::function<bool()>& terminate)
{
  search(true, stats, report, terminate);

  return best_path_labeling;
}

} // namespace bliss
