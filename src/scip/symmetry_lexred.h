#ifndef __SCIP_SYMMETRY_LEXRED_H__
#define __SCIP_SYMMETRY_LEXRED_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_event.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Data structures
 */

/** data for dynamic lexicographic reduction propagator */
struct SCIP_LexicographicReductionData;
typedef struct SCIP_LexicographicReductionData SCIP_LEXICOGRAPHICREDUCTIONDATA;

/*
 * Interface methods
 */

/** apply lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred,               /**< number of domain reductions */
   SCIP_Bool*            didrun              /**< whether propagator actually ran */
   );


/** adds permutation for lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionAddPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int*                  perm                /**< permutation */
   );


/** resets lexicographic reduction propagation (removes all permutations) */
SCIP_RETCODE SCIPlexicographicReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata/**< pointer to global data for lexicogrpahic order propagator */
   );


/** free lexicographic reduction data */
SCIP_RETCODE SCIPlexicographicReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA** masterdata/**< pointer to global data for lexicogrpahic order propagator */
   );


/** initializes structures needed for lexicographic reduction propagation
 * This is only done exactly once.
 */
SCIP_RETCODE SCIPlexicographicReductionInclude(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA** masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   );

#ifdef __cplusplus
}
#endif

#endif
