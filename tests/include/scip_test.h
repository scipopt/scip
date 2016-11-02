#include "scip/scip.h"

/* the method TESTsetSCIPStage(scip, stage) can be called in SCIP_STAGE_PROBLEM and can get to
 *  SCIP_STAGE_TRANSFORMING
 *  SCIP_STAGE_TRANSFORMED
 *  SCIP_STAGE_INITPRESOLVE
 *  SCIP_STAGE_PRESOLVING
 *  SCIP_STAGE_EXITPRESOLVE
 *  SCIP_STAGE_PRESOLVED
 *  SCIP_STAGE_SOLVING
 *  SCIP_STAGE_SOLVED
 *  SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE TESTscipSetStage(SCIP* scip, SCIP_STAGE stage);
#include "scip_test.c"

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include <criterion/parameterized.h>
#include <criterion/theories.h>

#undef SCIP_CALL
#define SCIP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             cr_assert(FALSE, "Error <%d> in function call\n", _restat_);                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )
