#include "scip/scip.h"

/* the method TESTsetSCIPStage(scip, stage) can be called in SCIP_STAGE_PROBLEM and can get to
 *  SCIP_STAGE_TRANSFORMED
 *  SCIP_STAGE_PRESOLVING
 *  SCIP_STAGE_PRESOLVED
 *  SCIP_STAGE_SOLVING
 *  SCIP_STAGE_SOLVED
 *
 *  If stage == SCIP_STAGE_SOLVING and enableNLP is true, then SCIP will build its NLP
 */
SCIP_RETCODE TESTscipSetStage(SCIP* scip, SCIP_STAGE stage, SCIP_Bool enableNLP);

/** assembles path of testfile using another files directory as directory name */
void TESTsetTestfilename(
   char*                 filename,           /**< buffer to write to, assumed to have length at least SCIP_MAXSTRLEN */
   const char*           file,               /**< name of file, usually including full path, from which to take directory name */
   const char*           testfile            /**< name of file to append, assumed to be in same directory as file */
);

/* Include the .c file here because the plugins implemented in scip_test.c should
 * be available every test.
 * */
#include "scip_test.c"
#include "locale.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#pragma GCC diagnostic ignored "-Wdeclaration-after-statement"
#endif

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include <criterion/parameterized.h>
#include <criterion/theories.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wredundant-decls"
#pragma GCC diagnostic warning "-Wstrict-prototypes"
#endif

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

CR_API int main(int argc, char *argv[]) {
    struct criterion_test_set *tests = criterion_initialize();

    int result = 0;
    if (criterion_handle_args(argc, argv, true))
    {
       setlocale(LC_ALL, "C");
       result = !criterion_run_all_tests(tests);
    }

    criterion_finalize(tests);
    return result;
}
