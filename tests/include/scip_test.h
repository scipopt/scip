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

/** sets filename to full path of testfile, assuming that testfile is in same directory as file */
static
void setfilename(
   char*                 filename,
   const char*           file,
   const char*           testfile
)
{
   char* pathsep;

   assert(filename != NULL);
   assert(testfile != NULL);

   /* get file to read: testfile that lives in the same directory as this file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", file);
   /* find last path separator */
#ifdef _WIN32
   pathsep = strrchr(filename, '\\');
#else
   pathsep = strrchr(filename, '/');
#endif
   /* overwrite filename */
   if( pathsep != NULL )
      (void)SCIPsnprintf(pathsep+1, SCIP_MAXSTRLEN - (pathsep+1 - filename), testfile);
   else
      (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, testfile);
}

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
