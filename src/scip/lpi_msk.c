/*  
  MOSEK SCIP interface. 
  
  Heavily revised in feb. 2007 by Bo Jensen bo.jensen@mosek.com
*/

#undef NDEBUG
#include <assert.h>

#define MSKCONST const
#include "mosek.h"

#include "scip/lpi.h"
#include "scip/bitencode.h"
#include <string.h>

#if MSK_VERSION_MAJOR >= 5
#define scipmskobjsen MSKobjsensee
#else
#define scipmskobjsen int
#endif

#if MSK_VERSION_MAJOR >= 5
#define SENSE2MOSEK(objsen) (((objsen)==SCIP_OBJSEN_MINIMIZE)?(MSK_OBJECTIVE_SENSE_MINIMIZE):(MSK_OBJECTIVE_SENSE_MAXIMIZE))
#else
#define SENSE2MOSEK(objsen) (((objsen)==SCIP_OBJSEN_MINIMIZE)?(MSK_OBJECTIVE_SENSE_MIN):(MSK_OBJECTIVE_SENSE_MAX))
#endif

#ifdef DEBUG_PRINT_MSK_CALLS
#define MOSEK_CALL(x)                                                               \
{ MSKrescodee _restat_ = (x);                                                             \
  printf("Calling MOSEK[%d]: %s\n",__LINE__,#x);                                    \
  if( (_restat_ ) != MSK_RES_OK && (_restat_ ) != MSK_RES_TRM_MAX_NUM_SETBACKS && (_restat_ ) != MSK_RES_ERR_PARAM_IS_TOO_SMALL  && (_restat_) != MSK_RES_ERR_PARAM_IS_TOO_LARGE)            \
  {                                                                                 \
    printf("optimizecount is %d retval %d MSK_RES_OK %d\n",optimizecount,(x),MSK_RES_OK);                                  \
    fflush(stdout);                                                                 \
    SCIPerrorMessage("LP Error: MOSEK returned %d in line %d\n",_restat_,__LINE__); \
    r = SCIP_LPERROR; goto CLEANUP;                                                 \
  }                                                                                 \
}
#else
#define MOSEK_CALL(x)                                                               \
{ MSKrescodee _restat_ = (x);                                                             \
  if( (_restat_ ) != MSK_RES_OK && (_restat_ ) != MSK_RES_TRM_MAX_NUM_SETBACKS && (_restat_ ) != MSK_RES_ERR_PARAM_IS_TOO_SMALL  && (_restat_) != MSK_RES_ERR_PARAM_IS_TOO_LARGE)            \
  {                                                                                 \
    printf("optimizecount is %d retval %d MSK_RES_OK %d\n",optimizecount,(_restat_),MSK_RES_OK);                                  \
    fflush(stdout);                                                                 \
    SCIPerrorMessage("LP Error: MOSEK returned %d in line %d\n",_restat_,__LINE__); \
    r = SCIP_LPERROR; goto CLEANUP;                                                 \
  }                                                                                 \
}
#endif

#define CHECK_ALLOC(x)                                                              \
{                                                                                   \
  if ((x) == NULL)                                                                  \
  {                                                                                 \
    fflush(stdout);                                                                 \
    SCIPerrorMessage("No memory in function call, line %d\n",__LINE__);             \
    r = SCIP_NOMEMORY; goto CLEANUP;                                                \
  }                                                                                 \
}

#define SCIP_CALL_W_R(x)                                                            \
{                                                                                   \
  r = (x);                                                                          \
  if (r != SCIP_OKAY)                                                               \
  {                                                                                 \
    fflush(stdout);                                                                 \
    SCIPerrorMessage("LP Error: SCIP return code %d\n",r);                          \
    goto CLEANUP;                                                                   \
  }                                                                                 \
}

#define RAISE_SCIP_ERROR(e) \
{                           \
  r = (e);                  \
  goto CLEANUP;             \
}

#define FREE_AND_NULL(ptr) { free(ptr); ptr = NULL; }
#define FREE_IF(ptr) { if (ptr != NULL) { FREE_AND_NULL(ptr); } }

#define IS_POSINF(x) ((x) >= SCIP_DEFAULT_INFINITY)
#define IS_NEGINF(x) ((x) <= -SCIP_DEFAULT_INFINITY)

#define STD_ASSERT { assert(MosekEnv != NULL); assert(lpi != NULL); assert(lpi->task != NULL); }

#define true TRUE
#define false FALSE


static MSKenv_t MosekEnv =           NULL;
static int numlp         =           0;

static int optimizecount            =  0;
static int forcescaling             =  0;
static int nextlpid                 =  1;
static int numstrongbranchmaxiterup =  0;
static int numstrongbranchmaxiterdo =  0;
static int numprimalmaxiter         =  0;
static int numdualmaxiter           =  0;
static int numstrongbranchobjup     =  0;
static int numstrongbranchobjdo     =  0;
static int numprimalobj             =  0;
static int numdualobj               =  0;

#define DEBUG_PRINT_STAT             0
#define DEBUG_PRINT_CALLS            0
#define DEBUG_PRINT_OPT_CALLS        0
#define DEBUG_CHECK_DATA             0
#define DEBUG_PRINT_MSK_CALLS        0
#define FORCE_MOSEK_LOG              0
#define FORCE_MOSEK_SUMMARY          0
#define FORCE_FREE_SIMPLEX           0
#define FORCE_NO_MAXITER             0
#define FORCE_SILENCE                1
#define FORCE_NO_PRESOLVE            1
#define FORCE_NO_SING                0
#define FORCE_NO_SAVE_LU             0
#define FORCE_CONTROL_SCALING        1
#define ASSERT_ON_NUMERICAL_TROUBLES 0
#define DEBUG_DO_INTPNT_FEAS_CHECK   0
#define ASSERT_ON_WARNING            0
#define DEBUG_CHECK_STATE            0
#define DEBUG_CHECK_STATE_TOL        1e-5
#define PACK_LPISTATE                1
#define SETBACK_LIMIT                250
#define SCIP_CONTROLS_PRICING        1
#define SCIP_CONTROLS_TOLERANCES     1
#define STRONGBRANCH_PRICING         MSK_SIM_SELECTION_FREE
#define USE_NAMES                    0
#define SUPRESS_NAME_ERROR           1
#define SHOW_WARNINGS                0
#define SHOW_ERRORS                  0
#define WRITE_DUAL                   0
#define WRITE_PRIMAL                 0
#define WRITE_INTPNT                 0
#define WRITE_ABOVE                 -1

#if MSK_VERSION_MAJOR >= 5
#define DEGEN_LEVEL                  MSK_SIM_DEGEN_MODERATE
#define ALWAYS_SOLVE_PRIMAL          1
#define AUTO_SWITCH_OPTIMIZER        1
#define CONTROL_DUALIZATION          0
#else
#define ALWAYS_SOLVE_PRIMAL          1
#define AUTO_SWITCH_OPTIMIZER        1
#define CONTROL_DUALIZATION          0
#endif


/**********************************************/

struct SCIP_LPi
{
    MSKtask_t task;
    MSKrescodee termcode;
    int itercount;
    int pricing;
    int lpid;
    int skxsize;
    int skcsize;
    MSKstakeye *skx;
    MSKstakeye *skc;
};

#if PACK_LPISTATE
typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE
/** returns the number of packets needed to store column packet information */
#define COLPACKET_NUM(ncols) ((ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET)
/** returns the number of packets needed to store row packet information */
#define ROWPACKET_NUM(nrows) ((nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET)
#else
#define COLPACKET MSKstakeye
#define ROWPACKET MSKstakeye
#define COLPACKET_NUM(ncols) (ncols)
#define ROWPACKET_NUM(nrows) (nrows)
#endif

struct SCIP_LPiState
{
    int ncols;
    int nrows;
    COLPACKET* skx;
    ROWPACKET* skc;
};

void MSKAPI printstr (void                 *handle, 
                      MSKCONST char        *str)
{
  #if SUPRESS_NAME_ERROR && !FORCE_SILENCE
  char errstr[32];
  snprintf(errstr,32,"MOSEK Error %d",MSK_RES_ERR_DUP_NAME);
  if (0 == strncmp(errstr,str,strlen(errstr)))
      return;
  #endif

  #if !FORCE_SILENCE
  printf("MOSEK: %s",str); fflush (stdout);
  #endif
}

#if DEBUG_CHECK_DATA > 0   
static void scip_checkdata(SCIP_LPI* lpi,
                           char      where[])
{
  MSKrescodee  r = MSK_RES_OK;
  int          i,numcon,numvar,gotbasicsol;
  MSKboundkeye *tbkc=NULL,*tbkx=NULL;
  MSKstakeye *tskc=NULL,*tskx=NULL;
  double       *tblc=NULL,*tbuc=NULL,*tblx=NULL,*tbux=NULL;

  MOSEK_CALL( MSK_solutiondef(lpi->task,
                              MSK_SOL_BAS,
                              &gotbasicsol) );


  MOSEK_CALL( MSK_getnumvar(lpi->task,&numvar) );
  MOSEK_CALL( MSK_getnumcon(lpi->task,&numcon) );

  /* Check bounds */
  CHECK_ALLOC( tbkc = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * numcon) );
  CHECK_ALLOC( tskc = (MSKstakeye*) malloc(sizeof(MSKstakeye) * numcon) );
  CHECK_ALLOC( tblc = (double*) malloc(sizeof(double) * numcon) );
  CHECK_ALLOC( tbuc = (double*) malloc(sizeof(double) * numcon) );

  CHECK_ALLOC( tbkx = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * numvar) );
  CHECK_ALLOC( tskx = (MSKstakeye*) malloc(sizeof(MSKstakeye) * numvar) );
  CHECK_ALLOC( tblx = (double*) malloc(sizeof(double) * numvar) );
  CHECK_ALLOC( tbux = (double*) malloc(sizeof(double) * numvar) );

  if( gotbasicsol )
    MOSEK_CALL( MSK_getsolution(lpi->task,
                                MSK_SOL_BAS,
                                NULL,
                                NULL,
                                tskc,
                                tskx,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL) );


  for( i = 0; i < numvar; i++ )
    MOSEK_CALL( MSK_getbound(lpi->task,MSK_ACC_VAR,i,&tbkx[i],&tblx[i],&tbux[i]) );
   
  for( i = 0; i < numcon; i++ )
    MOSEK_CALL( MSK_getbound(lpi->task,MSK_ACC_CON,i,&tbkc[i],&tblc[i],&tbuc[i]) );


  for( i = 0; i < numcon; ++i )
  {
    if( gotbasicsol )
    {
      if( ( tskc[i] == MSK_SK_FIX && tbkc[i] != MSK_BK_FX ) ||
          ( tskc[i] == MSK_SK_LOW && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_RA ) ) ||
          ( tskc[i] == MSK_SK_UPR && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_RA ) ) )
      {
        printf("STATUS KEY ERROR i %d bkc %d skc %d %s\n",i,tbkc[i],tskc[i],where);
        getchar();
      }
    }

    if( tbkc[i] == MSK_BK_LO ||
        tbkc[i] == MSK_BK_FX ||
        tbkc[i] == MSK_BK_RA ) 
    {
      if( isnan(tblc[i]) )
      {
        printf("nan in blc : %s\n",where);
        getchar();
      }
    }  

    if( tbkc[i] == MSK_BK_UP ||
        tbkc[i] == MSK_BK_FX ||
        tbkc[i] == MSK_BK_RA ) 
    {
      if( isnan(tbuc[i]) )
      {
        printf("nan in bux : %s\n",where);
        getchar();
      }
    }  
  }

  for( i = 0; i < numvar; ++i )
  {
    if( tbkx[i] == MSK_BK_LO ||
        tbkx[i] == MSK_BK_FX ||
        tbkx[i] == MSK_BK_RA ) 
    {
      if( isnan(tblx[i]) )
      {
        printf("nan in blx : %s\n",where);
        getchar();
      }
    }  

    if( tbkx[i] == MSK_BK_UP ||
        tbkx[i] == MSK_BK_FX ||
        tbkx[i] == MSK_BK_RA ) 
    {
      if( isnan(tbux[i]) )
      {
        printf("nan in bux : %s\n",where);
        getchar();
      }
    }  
  }

  CLEANUP:
    FREE_IF(tbkc);
    FREE_IF(tskc);
    FREE_IF(tblc);
    FREE_IF(tbuc);
    FREE_IF(tbkx);
    FREE_IF(tskx);
    FREE_IF(tblx);
    FREE_IF(tbux);
}
#endif


/*
 * Local functions
 */

static void generate_mskbounds(int          n, 
                               const double *lb, 
                               const double *ub, 
                               MSKboundkeye *bk, 
                               double       *msklb, 
                               double       *mskub)
{
  int i;
  assert(lb);
  assert(ub);
  assert(bk);
  assert(msklb);
  assert(mskub);

  for (i=0; i<n; i++)
  {
    msklb[i] = lb[i];
    mskub[i] = ub[i];
    if (IS_NEGINF(lb[i]))
    {
      msklb[i] = -MSK_INFINITY;
      if (IS_POSINF(ub[i]))
      {
          mskub[i] = MSK_INFINITY;
          bk[i] = MSK_BK_FR;
      }
      else
      {
          assert(!IS_NEGINF(ub[i]));
          bk[i] = MSK_BK_UP;
      }
    }
    else
    {
      assert(!IS_POSINF(lb[i]));
      if (IS_POSINF(ub[i]))
      {
          mskub[i] = MSK_INFINITY;
          bk[i] = MSK_BK_LO;
      }
      else if (lb[i] == ub[i])
      {
          assert(lb[i]-ub[i]==0);
          assert(ub[i]-lb[i]==0);
          bk[i] = MSK_BK_FX;
      }
      else
      {
          assert(lb[i] < ub[i]);
          bk[i] = MSK_BK_RA;
      }
    }
  }
}

static int* get_endptrs(int        n, 
                        const int* beg, 
                        int        nnonz)
{
  int i;
  int *aptre;

  assert(beg != NULL || nnonz == 0);
  aptre = (int *) malloc(sizeof(int) * n);

  if (aptre == NULL) 
    return NULL;

  if (nnonz > 0)
  {
    for(i=0; i < n-1; i++)
    {
      aptre[i] = beg[i+1];
      assert(aptre[i] >= beg[i]);
    }

    aptre[n-1] = nnonz;
    assert(aptre[n-1] >= beg[n-1]);
  }
  else
  {
      for (i=0; i<n; i++)
          aptre[i] = 0;
  }

  return aptre;
}

static int* get_indices_range(int first, 
                              int last)
{
  int i;
  int *sub;

  assert(first <= last);

  sub = (int*) malloc(sizeof(int) * (last-first+1));

  if (sub == NULL) 
    return NULL;

  for (i=first; i<=last; i++)
  {
    sub[i-first] = i;
  }

  return sub;
}

static int* get_indices_from_dense(int *dstat, 
                                   int n, 
                                   int *count)
{
  int i, j;
  int *sub;
  assert(dstat != NULL);

  *count = 0;
  for (i=0; i<n; i++) 
  {
    if (dstat[i] == 1)
    {
      (*count)++;
    }
  }

  sub = (int*) malloc(sizeof(int) * (*count));

  if (sub == NULL) 
    return NULL;

  j = 0;
  for (i=0; i<n; i++)
  { 
    if (dstat[i] == 1)
    {
      sub[j++] = i;
    }
  }
  
  return sub;
}

static void scale_vec(int    len, 
                      double *vec, 
                      double s)
{
  int i;
  for (i=0; i<len; i++)
  {
    vec[i] *= s;
  }
}

static void scale_bound(MSKboundkeye *bk, 
                        double       *bl, 
                        double       *bu, 
                        double       s)
{
  switch(*bk)
  {
    case MSK_BK_LO:
      *bl *= s;
      if (s < 0) *bk = MSK_BK_UP;
      break;
    case MSK_BK_UP:
      *bu *= s;
      if (s < 0) *bk = MSK_BK_LO;
      break;
    case MSK_BK_FX:
    case MSK_BK_RA:
      *bl *= s;
      *bu *= s;
      break;
    case MSK_BK_FR:
      break;
    default:
      assert(false);
      break;
  }

  if (s < 0)
  {
    double tmp;
    tmp = *bl;
    *bl = *bu;
    *bu = tmp;
  }
}

static SCIP_RETCODE ensureStateMem(SCIP_LPI* lpi, 
                                   int       ncols, 
                                   int       nrows)
{
  int r = SCIP_OKAY;
  if (lpi->skxsize < ncols)
  {
    int newsize;
    newsize = MAX(2*lpi->skxsize, ncols);

    CHECK_ALLOC( lpi->skx = (MSKstakeye*) realloc(lpi->skx,sizeof(MSKstakeye) * newsize) );

    lpi->skxsize = newsize;
  }

  if (lpi->skcsize < nrows)
  {
    int newsize;
    newsize = MAX(2*lpi->skcsize, nrows);

    CHECK_ALLOC( lpi->skc = (MSKstakeye*) realloc(lpi->skc,sizeof(MSKstakeye) * newsize) );

    lpi->skcsize = newsize;
  }

  CLEANUP:

  return r;
}

static SCIP_RETCODE getbase(SCIP_LPI* lpi, 
                            int       ncols, 
                            int       nrows)
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling getbase (%d)\n",lpi->lpid);
  #endif

  SCIP_CALL_W_R( ensureStateMem(lpi,ncols,nrows) );

  MOSEK_CALL( MSK_getsolution(lpi->task,
                              MSK_SOL_BAS,
                              NULL,
                              NULL,
                              lpi->skc,
                              lpi->skx,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL) );

  CLEANUP:

  return r;
}

static SCIP_RETCODE setbase(SCIP_LPI* lpi)
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling setbase (%d)\n",lpi->lpid);
  #endif

  MOSEK_CALL( MSK_putsolution(lpi->task,
                              MSK_SOL_BAS,
                              lpi->skc,
                              lpi->skx,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL) );

  #if MSK_VERSION_MAJOR >= 5
  /* We only have status keys (recalc dual solution with out dual superbasics) */
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_HOTSTART,
                              MSK_SIM_HOTSTART_STATUS_KEYS) );
  #endif


  CLEANUP:

  return r;
}



/*
 * Miscellaneous Methods
 */

static char mosekname[] = "MOSEK";

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(void)
{
    return mosekname;
}


/** gets pointer for LP solver - use only with great care */
void* SCIPgetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->task;
}


/*
 * LPI Creation and Destruction Methods
 */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
                           const char*           name,               /**< problem name */
                           SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiCreate\n");
  #endif

  assert(lpi != NULL);

  if (!MosekEnv)
  {
    MOSEK_CALL( MSK_makeenv(&MosekEnv,
                            NULL,
                            NULL,
                            NULL,
                            NULL) );

    MOSEK_CALL( MSK_linkfunctoenvstream(MosekEnv,
                                        MSK_STREAM_LOG,
                                        NULL,
                                        printstr) );

    MOSEK_CALL( MSK_initenv(MosekEnv) );
  }

  numlp++;

  CHECK_ALLOC( *lpi = (SCIP_LPI*) malloc(sizeof(SCIP_LPI)) );

  MOSEK_CALL( MSK_makeemptytask(MosekEnv,
                                &((*lpi)->task)) );

  MOSEK_CALL( MSK_linkfunctotaskstream((*lpi)->task,
                                       MSK_STREAM_LOG,
                                       NULL,
                                       printstr) );

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putobjsense((*lpi)->task,
                               SENSE2MOSEK(objsen)) );
  #else
  MOSEK_CALL( MSK_putintparam((*lpi)->task,
                              MSK_IPAR_OBJECTIVE_SENSE,
                              SENSE2MOSEK(objsen)) );
  #endif

  MOSEK_CALL( MSK_putintparam((*lpi)->task,
                              MSK_IPAR_SIM_MAX_NUM_SETBACKS,
                              SETBACK_LIMIT) );

  MOSEK_CALL( MSK_putintparam((*lpi)->task,
                              MSK_IPAR_OPTIMIZER,
                              MSK_OPTIMIZER_FREE_SIMPLEX) );

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putintparam((*lpi)->task,
                              MSK_IPAR_SIM_DEGEN,
                              DEGEN_LEVEL) );
  #endif


  #if AUTO_SWITCH_OPTIMIZER
  MOSEK_CALL( MSK_putintparam((*lpi)->task,
                               MSK_IPAR_SIM_SWITCH_OPTIMIZER,
                               MSK_ON) );
  #endif

  MOSEK_CALL( MSK_puttaskname((*lpi)->task,
                              (char*) name) );

  (*lpi)->termcode   = MSK_RES_OK;
  (*lpi)->itercount  = 0;
  (*lpi)->pricing    = SCIP_PRICING_AUTO;
  (*lpi)->lpid       = nextlpid++;
  (*lpi)->skxsize    = 0;
  (*lpi)->skcsize    = 0;
  (*lpi)->skx        = NULL;
  (*lpi)->skc        = NULL;

  CLEANUP:

  return r;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
  int r = SCIP_OKAY;

  assert(lpi != NULL);
  assert(*lpi != NULL);
  assert(numlp > 0);

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiFree (%d)\n",(*lpi)->lpid);
  #endif

  MOSEK_CALL( MSK_deletetask(&(*lpi)->task) );

  FREE_IF((*lpi)->skx);
  FREE_IF((*lpi)->skc);

  FREE_AND_NULL(*lpi);

  numlp--;
  if (numlp == 0)
  {
    MOSEK_CALL( MSK_deleteenv(&MosekEnv) );
    MosekEnv = NULL;
  }

  CLEANUP:
  return r;
}

/*
 * Modification Methods
 */

#if USE_NAMES
static MSKrescodee ignore_dup_names(MSKrescodee res)
{
  if (res == MSK_RES_ERR_DUP_NAME)
  {
    #if SHOW_WARNINGS
    printf("Ignoring MSK_RES_ERR_DUP_NAME\n");
    #endif

    return MSK_RES_OK;
  }
  else
    return res;
}
#endif

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(SCIP_LPI*             lpi,                /**< LP interface structure */
                              SCIP_OBJSEN           objsen,             /**< objective sense */
                              int                   ncols,              /**< number of columns */
                              const SCIP_Real*      obj,                /**< objective function values of columns */
                              const SCIP_Real*      lb,                 /**< lower bounds of columns */
                              const SCIP_Real*      ub,                 /**< upper bounds of columns */
                              char**                colnames,           /**< column names, or NULL */
                              int                   nrows,              /**< number of rows */
                              const SCIP_Real*      lhs,                /**< left hand sides of rows */
                              const SCIP_Real*      rhs,                /**< right hand sides of rows */
                              char**                rownames,           /**< row names, or NULL */
                              int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
                              const int*            beg,                /**< start index of each column in ind- and val-array */
                              const int*            ind,                /**< row indices of constraint matrix entries */
                              const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
  #if DEBUG_CHECK_DATA > 0   
  char where[] = "SCIPlpiLoadColLP";
  #endif

  int r = SCIP_OKAY;

  int *aptre        = NULL;
  MSKboundkeye *bkc = NULL;
  MSKboundkeye *bkx = NULL;
  double *blc       = NULL;
  double *buc       = NULL;
  double *blx       = NULL;
  double *bux       = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiLoadColLP (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  if (nrows > 0)
  {
      CHECK_ALLOC( bkc = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * nrows) );
      CHECK_ALLOC( blc = (double*) malloc(sizeof(double) * nrows) );
      CHECK_ALLOC( buc = (double*) malloc(sizeof(double) * nrows) );

      generate_mskbounds(nrows,lhs,rhs,bkc,blc,buc);
  }

  if (ncols > 0)
  {
      CHECK_ALLOC( bkx = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * ncols) );
      CHECK_ALLOC( blx = (double*) malloc(sizeof(double) * ncols) );
      CHECK_ALLOC( bux = (double*) malloc(sizeof(double) * ncols) );

      generate_mskbounds(ncols,lb,ub,bkx,blx,bux);

      CHECK_ALLOC( aptre = get_endptrs(ncols,beg,nnonz) );
  }

  MOSEK_CALL( MSK_inputdata(lpi->task,
                            nrows,
                            ncols,
                            nrows,
                            ncols,
                            obj,
                            0.0,
                            beg,
                            aptre,
                            ind,
                            val,
                            bkc,
                            blc,
                            buc,
                            bkx,
                            blx,
                            bux) );

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putobjsense(lpi->task,
                              SENSE2MOSEK(objsen)) );
  #else
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_OBJECTIVE_SENSE,
                              SENSE2MOSEK(objsen)) );
  #endif

  #if USE_NAMES
  if (colnames)
  {
    int i;
    for (i=0; i<ncols; i++)
    {
      /*printf("Assigning name of var[%d] \"%s\"\n",i,colnames[i]);*/
      MOSEK_CALL( ignore_dup_names(MSK_putname(lpi->task,
                                               MSK_PI_VAR,
                                               i,
                                               colnames[i])) );
    }
  }

  if (rownames)
  {
    int i;
    for (i=0; i<nrows; i++)
    {
      /*printf("Assigning name of con[%d] \"%s\"\n",i,rownames[i]);*/
      MOSEK_CALL( ignore_dup_names(MSK_putname(lpi->task,
                                   MSK_PI_CON,
                                   i,
                                   rownames[i])) );
    }
  }
  #endif

  CLEANUP:

  FREE_IF(aptre);
  FREE_IF(bkc);
  FREE_IF(blc);
  FREE_IF(buc);
  FREE_IF(bkx);
  FREE_IF(blx);
  FREE_IF(bux);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   ncols,              /**< number of columns to be added */
                            const SCIP_Real*      obj,                /**< objective function values of new columns */
                            const SCIP_Real*      lb,                 /**< lower bounds of new columns */
                            const SCIP_Real*      ub,                 /**< upper bounds of new columns */
                            char**                colnames,           /**< column names, or NULL */
                            int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
                            const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
                            const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
                            const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiAddCols";
  #endif

  int r             = SCIP_OKAY;
  const int *aptrb  = NULL;
  int *aptre        = NULL;
  MSKboundkeye *bkx = NULL;
  double *blx       = NULL;
  double *bux       = NULL;
  int oldcols;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiAddCols (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  if (ncols == 0) 
    goto CLEANUP;

  CHECK_ALLOC( bkx = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * ncols) );
  CHECK_ALLOC( blx = (double*) malloc(sizeof(double) * ncols) );
  CHECK_ALLOC( bux = (double*) malloc(sizeof(double) * ncols) );
  generate_mskbounds(ncols,lb,ub,bkx,blx,bux);

  MOSEK_CALL( MSK_getnumvar(lpi->task,&oldcols) );

  CHECK_ALLOC( aptre = get_endptrs(ncols,beg,nnonz) );

  if (nnonz == 0)
    aptrb = aptre;
  else
    aptrb = beg;

  MOSEK_CALL( MSK_appendvars(lpi->task,
                             ncols,
                             obj,
                             aptrb,
                             aptre,
                             ind,
                             val,
                             bkx,
                             blx,
                             bux) );

  #if USE_NAMES
  if (colnames)
  {
    int i;
    for (i=0; i<ncols; i++)
    {
      /*printf("Assigning name of var[%d] \"%s\"\n",oldcols+i,colnames[i]);*/
      MOSEK_CALL( ignore_dup_names(MSK_putname(lpi->task,
                                               MSK_PI_VAR,
                                               oldcols+i,
                                               colnames[i])) );
    }
  }
  #endif

  CLEANUP:
  FREE_IF(aptre);
  FREE_IF(bkx);
  FREE_IF(blx);
  FREE_IF(bux);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   firstcol,           /**< first column to be deleted */
                            int                   lastcol             /**< last column to be deleted */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiDelCols";
  #endif
  int r    = SCIP_OKAY;
  int *sub = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiDelCols (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  CHECK_ALLOC( sub = get_indices_range(firstcol,lastcol) );

  /*printf("Deleting vars %d to %d\n",firstcol,lastcol);*/
  MOSEK_CALL( MSK_remove(lpi->task,MSK_ACC_VAR,
                         lastcol-firstcol+1,
                         sub) );

  CLEANUP:
  FREE_IF(sub);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(SCIP_LPI*             lpi,                /**< LP interface structure */
                              int*                  dstat               /**< deletion status of columns
                                               *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiDelColset";
  #endif
  int r = SCIP_OKAY;
  int *sub = NULL;
  int count, i, ncols, col;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiDelColset (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  CHECK_ALLOC( sub = get_indices_from_dense(dstat,
                                            ncols,
                                            &count) );

  col = 0;
  for (i=0; i<ncols; i++)
  {
    if (dstat[i] == 1)
    {
      dstat[i] = -1;
    }
    else
    {
      dstat[i] = col;
      col++;
    }
  }

  if (count > 0)
  {
    /*printf("Deleting %d vars %d,...\n",count,sub[0]);*/
    MOSEK_CALL( MSK_remove(lpi->task,
                           MSK_ACC_VAR,
                           count,
                           sub) );
  }

  CLEANUP:
  FREE_IF(sub);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   nrows,              /**< number of rows to be added */
                            const SCIP_Real*      lhs,                /**< left hand sides of new rows */
                            const SCIP_Real*      rhs,                /**< right hand sides of new rows */
                            char**                rownames,           /**< row names, or NULL */
                            int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
                            const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
                            const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
                            const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiAddRows";
  #endif
  int r = SCIP_OKAY;


  const int *aptrb  = NULL;
  int *aptre        = NULL;
  MSKboundkeye *bkc = NULL;
  double *blc       = NULL;
  double *buc       = NULL;
  int oldrows;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiAddRows (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  if (nrows == 0) goto CLEANUP;

  CHECK_ALLOC( bkc = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * nrows) );
  CHECK_ALLOC( blc = (double*) malloc(sizeof(double) * nrows) );
  CHECK_ALLOC( buc = (double*) malloc(sizeof(double) * nrows) );

  generate_mskbounds(nrows,
                     lhs,
                     rhs,
                     bkc,
                     blc,
                     buc);

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &oldrows) );

  CHECK_ALLOC( aptre = get_endptrs(nrows,
                                   beg,
                                   nnonz) );

  if (nnonz == 0)
    aptrb = aptre;
  else
    aptrb = beg;

  MOSEK_CALL( MSK_appendcons(lpi->task,
                             nrows,
                             aptrb,
                             aptre,
                             ind,
                             val,
                             bkc,
                             blc,
                             buc) );

  #if USE_NAMES
  if (rownames)
  {
    int i;
    for (i=0; i<nrows; i++)
    {
      /*printf("Assigning name of con[%d] \"%s\"\n",oldrows+i,rownames[i]);*/
      MOSEK_CALL( ignore_dup_names(MSK_putname(lpi->task,
                                               MSK_PI_CON,
                                               oldrows+i,
                                               rownames[i])) );
    }
  }
  #endif

  CLEANUP:
  FREE_IF(aptre);
  FREE_IF(bkc);
  FREE_IF(blc);
  FREE_IF(buc);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   firstrow,           /**< first row to be deleted */
                            int                   lastrow             /**< last row to be deleted */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiDelRows";
  #endif
  int r = SCIP_OKAY;
  int* sub = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiDelRows (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  CHECK_ALLOC( sub = get_indices_range(firstrow,lastrow) );

  /*printf("Deleting cons %d to %d\n",firstrow,lastrow);*/

  MOSEK_CALL( MSK_remove(lpi->task,
                         MSK_ACC_CON,
                         lastrow-firstrow+1,
                         sub) );

  CLEANUP:
  FREE_IF(sub);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(SCIP_LPI*             lpi,                /**< LP interface structure */
                              int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiDelRowset";
  char bwhere[] = "SCIPlpiDelRowset end";
  #endif
  int r = SCIP_OKAY;
  int *sub = NULL;
  int count, i, nrows, row;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiDelRowset (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  CHECK_ALLOC( sub = get_indices_from_dense(dstat, 
                                            nrows,
                                            &count) );

  row = 0;
  for (i=0; i<nrows; i++)
  {
    if (dstat[i] == 1)
    {
      dstat[i] = -1;
    }
    else
    {
      dstat[i] = row;
      row++;
    }
  }

  if (count > 0)
  {
    /*printf("Deleting %d cons %d,...\n",count,sub[0]);*/

    MOSEK_CALL( MSK_remove(lpi->task,
                           MSK_ACC_CON,
                           count,
                           sub) );
  }

  CLEANUP:
  FREE_IF(sub);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 bwhere);
  #endif

  return r;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  int nrows, ncols;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiClear (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  SCIP_CALL_W_R( SCIPlpiDelRows(lpi,
                                0,
                                nrows) );

  SCIP_CALL_W_R( SCIPlpiDelCols(lpi,0,ncols) );

  CLEANUP:

  return r;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(SCIP_LPI*             lpi,                /**< LP interface structure */
                              int                   ncols,              /**< number of columns to change bounds for */
                              const int*            ind,                /**< column indices */
                              const SCIP_Real*      lb,                 /**< values for the new lower bounds */
                              const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiChgBounds";
  #endif
  int r = SCIP_OKAY;
  MSKboundkeye *bkx = NULL;
  double *blx       = NULL;
  double *bux       = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiChgBounds (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  if (ncols == 0) goto CLEANUP;

  CHECK_ALLOC( bkx = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * ncols) );
  CHECK_ALLOC( blx = (double*) malloc(sizeof(double) * ncols) );
  CHECK_ALLOC( bux = (double*) malloc(sizeof(double) * ncols) );

  generate_mskbounds(ncols,
                     lb,
                     ub,
                     bkx,
                     blx,
                     bux);

  MOSEK_CALL( MSK_putboundlist(lpi->task,
                               MSK_ACC_VAR,
                               ncols,
                               ind, 
                               bkx,
                               blx,
                               bux) );


  CLEANUP:
  FREE_IF(bkx);
  FREE_IF(blx);
  FREE_IF(bux);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int                   nrows,              /**< number of rows to change sides for */
                             const int*            ind,                /**< row indices */
                             const SCIP_Real*      lhs,                /**< new values for left hand sides */
                             const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiChgSides";
    #endif
  int r = SCIP_OKAY;
  MSKboundkeye *bkc = NULL;
  double *blc       = NULL;
  double *buc       = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiChgSides (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  if (nrows == 0) goto CLEANUP;

  CHECK_ALLOC( bkc = (MSKboundkeye*) malloc(sizeof(MSKboundkeye) * nrows) );
  CHECK_ALLOC( blc = (double*) malloc(sizeof(double) * nrows) );
  CHECK_ALLOC( buc = (double*) malloc(sizeof(double) * nrows) );

  generate_mskbounds(nrows,
                     lhs,
                     rhs,
                     bkc,
                     blc,
                     buc);

  MOSEK_CALL( MSK_putboundlist(lpi->task,
                               MSK_ACC_CON,
                               nrows,
                               ind,
                               bkc,
                               blc,
                               buc) );

  CLEANUP:
  FREE_IF(bkc);
  FREE_IF(blc);
  FREE_IF(buc);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   row,                /**< row number of coefficient to change */
                            int                   col,                /**< column number of coefficient to change */
                            SCIP_Real             newval              /**< new value of coefficient */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiChgCoef";
    #endif
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiChgCoef (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  MOSEK_CALL( MSK_putaij(lpi->task,
                         row,
                         col,
                         newval) );

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(SCIP_LPI*             lpi,                /**< LP interface structure */
                              SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiChgObjsen (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putobjsense(lpi->task,
                              SENSE2MOSEK(objsen)) );
  #else
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_OBJECTIVE_SENSE,
                              SENSE2MOSEK(objsen)) );
  #endif

  CLEANUP:

  return r;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(SCIP_LPI*             lpi,                /**< LP interface structure */
                           int                   ncols,              /**< number of columns to change objective value for */
                           int*                  ind,                /**< column indices to change objective value for */
                           SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiChgObj";
    #endif
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiChgObj (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_putclist(lpi->task,
                           ncols,
                           ind,
                           obj) );

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int                   row,                /**< row number to scale */
                             SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiScaleRow";
    #endif
  int r = SCIP_OKAY;
  int nnonz;
  int *sub    = NULL;
  double *val = NULL;
  MSKboundkeye bkc;
  double blc, buc;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiScaleRow (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  assert(scaleval != 0);

  MOSEK_CALL( MSK_getavecnumnz(lpi->task,
                               MSK_ACC_CON,
                               row,
                               &nnonz) );

  if (nnonz != 0)
  {
    CHECK_ALLOC( sub = (int*) malloc(sizeof(int) * nnonz) );
    CHECK_ALLOC( val = (double*) malloc(sizeof(double) * nnonz) );

    MOSEK_CALL( MSK_getavec(lpi->task,
                            MSK_ACC_CON,
                            row,
                            &nnonz,
                            sub,
                            val) );

    scale_vec(nnonz,
              val, 
              scaleval);

    MOSEK_CALL( MSK_putavec(lpi->task,
                            MSK_ACC_CON,
                            row,
                            nnonz,
                            sub,
                            val) );
  }

  MOSEK_CALL( MSK_getbound(lpi->task,
                           MSK_ACC_CON,
                           row,
                           &bkc,
                           &blc,
                           &buc) );

  scale_bound(&bkc,
              &blc,
              &buc,
              scaleval);

  MOSEK_CALL( MSK_putbound(lpi->task,
                           MSK_ACC_CON,
                           row,
                           bkc,
                           blc,
                           buc) );

  CLEANUP:
  FREE_IF(sub);
  FREE_IF(val);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int                   col,                /**< column number to scale */
                             SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiScaleCol";
  #endif

  int r = SCIP_OKAY;
  int nnonz;
  int *sub    = NULL;
  double *val = NULL;
  MSKboundkeye bkx;
  double blx, bux, c;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiScaleCol (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  
  assert(scaleval != 0);

  MOSEK_CALL( MSK_getavecnumnz(lpi->task,
                               MSK_ACC_VAR,
                               col,
                               &nnonz) );

  if (nnonz != 0)
  {
    CHECK_ALLOC( sub = (int*) malloc(sizeof(int) * nnonz) );
    CHECK_ALLOC( val = (double*) malloc(sizeof(double) * nnonz) );

    MOSEK_CALL( MSK_getavec(lpi->task,
                            MSK_ACC_VAR,
                            col,
                            &nnonz,
                            sub,
                            val) );

    scale_vec(nnonz,
              val,
              scaleval);

    MOSEK_CALL( MSK_putavec(lpi->task,
                            MSK_ACC_VAR,
                            col,
                            nnonz,
                            sub,
                            val) );
  }

  MOSEK_CALL( MSK_getbound(lpi->task,
                           MSK_ACC_VAR,
                           col,
                           &bkx,
                           &blx,
                           &bux) );

  scale_bound(&bkx,
              &blx,
              &bux,
              1.0/scaleval);

  MOSEK_CALL( MSK_putbound(lpi->task,
                           MSK_ACC_VAR,
                           col,
                           bkx,
                           blx,
                           bux) );

  MOSEK_CALL( MSK_getcslice(lpi->task,
                            col,
                            col+1,
                            &c) );

  MOSEK_CALL( MSK_putcj(lpi->task,
                        col,
                        c*scaleval) );

  CLEANUP:
  FREE_IF(sub);
  FREE_IF(val);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}


/*
 * Data Accessing Methods
 */


/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int*                  nrows               /**< pointer to store the number of rows */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetNRows (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            nrows) );

  CLEANUP:

  return r;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int*                  ncols               /**< pointer to store the number of cols */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetNCols (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            ncols) );

  CLEANUP:

  return r;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetNNonz (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumanz(lpi->task,
                            nnonz) );

  CLEANUP:

  return r;
}

static SCIP_RETCODE getaslice(SCIP_LPI* lpi,
                              int       iscon,
                              int       first,
                              int       last,
                              int*      nnonz,
                              int*      beg,
                              int*      ind,
                              double*   val)
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "getaslice";
  #endif

  int r = SCIP_OKAY;
  int *aptre = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetNNonz (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;
  assert(first <= last);

  if (nnonz)
  {
    int surplus;

    assert(beg);
    assert(ind);
    assert(val);

    CHECK_ALLOC( aptre = (int*) malloc(sizeof(int) * (last - first + 1)) );

    MOSEK_CALL( MSK_getaslicenumnz(lpi->task,
                                   iscon,
                                   first,
                                   last+1,nnonz) );

    surplus = *nnonz;

    MOSEK_CALL( MSK_getaslice(lpi->task,
                              iscon,
                              first,
                              last+1,
                              *nnonz,
                              &surplus,
                              beg,
                              aptre,
                              ind,
                              val) );

    assert(surplus == 0);
  }

  CLEANUP:
  FREE_IF(aptre);

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   firstcol,           /**< first column to get from LP */
                            int                   lastcol,            /**< last column to get from LP */
                            SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
                            SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
                            int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
                            int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
                            int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
                            SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetCols (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  SCIP_CALL_W_R( SCIPlpiGetBounds(lpi,
                                  firstcol,
                                  lastcol,
                                  lb,
                                  ub) );

  SCIP_CALL_W_R( getaslice(lpi,
                           MSK_ACC_VAR,
                           firstcol,
                           lastcol,
                           nnonz,
                           beg,
                           ind,
                           val) );

  CLEANUP:

  return r;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   firstrow,           /**< first row to get from LP */
                            int                   lastrow,            /**< last row to get from LP */
                            SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
                            SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
                            int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
                            int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
                            int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
                            SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiGetRows";
  #endif

  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetRows (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  SCIP_CALL_W_R( SCIPlpiGetSides(lpi,
                                 firstrow,
                                 lastrow,
                                 lhs,
                                 rhs) );

  SCIP_CALL_W_R( getaslice(lpi,
                           MSK_ACC_CON,
                           firstrow,
                           lastrow,
                           nnonz,
                           beg,
                           ind,
                           val) );

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(SCIP_LPI*             lpi,                /**< LP interface structure */
                           int                   firstcol,           /**< first column to get objective coefficient for */
                           int                   lastcol,            /**< last column to get objective coefficient for */
                           SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetObj (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getcslice(lpi->task,
                            firstcol,
                            lastcol+1,  
                            vals) );

  CLEANUP:

  return r;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(SCIP_LPI*             lpi,                /**< LP interface structure */
                              int                   firstcol,           /**< first column to get bounds for */
                              int                   lastcol,            /**< last column to get bounds for */
                              SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
                              SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiGetBounds";
  #endif

  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBounds (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getboundslice(lpi->task,
                                MSK_ACC_VAR,
                                firstcol,
                                lastcol+1,
                                NULL,
                                lbs,
                                ubs) );

  CLEANUP:

  return r;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(SCIP_LPI*             lpi,                /**< LP interface structure */
                             int                   firstrow,           /**< first row to get sides for */
                             int                   lastrow,            /**< last row to get sides for */
                             SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
                             SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiGetSides";
  #endif

  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetSides (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getboundslice(lpi->task,
                                MSK_ACC_CON,
                                firstrow,
                                lastrow+1,
                                NULL,
                                lhss,
                                rhss) );

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int                   row,                /**< row number of coefficient */
                            int                   col,                /**< column number of coefficient */
                            SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiGetCoef";
  #endif

  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetCoef (%d)\n",lpi->lpid);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getaij(lpi->task,
                         row,
                         col,
                         val) );

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/*
 * Solving Methods
 */

static MSKrescodee filterTRMrescode(MSKrescodee* termcode, 
                                    MSKrescodee  res)
{
  if (   res == MSK_RES_TRM_MAX_ITERATIONS
      || res == MSK_RES_TRM_MAX_TIME
      || res == MSK_RES_TRM_OBJECTIVE_RANGE
    /*|| res == MSK_RES_TRM_USER_BREAK*/
    /*|| res == MSK_RES_TRM_STALL*/
      #if ASSERT_ON_NUMERICAL_TROUBLES > 0
      || res == MSK_RES_TRM_MAX_NUM_SETBACKS
      || res == MSK_RES_TRM_NUMERICAL_PROBLEM
      #endif
    /*|| res == MSK_RES_TRM_INTERNAL*/)
  {
    *termcode = res;
    if (res == MSK_RES_TRM_MAX_NUM_SETBACKS || res == MSK_RES_TRM_NUMERICAL_PROBLEM)
    {
        #if SHOW_WARNINGS && !FORCE_SILENCE
        printf("WARNING: Return code %d in [%d]\n",res,optimizecount);
        #endif
        
        #if ASSERT_ON_WARNING
        assert(0);
        #endif
    }

    return MSK_RES_OK;
  }
  else
  {
    *termcode = MSK_RES_OK;

    return res;
  }
}

static SCIP_RETCODE SolveWSimplex(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SolveWSimplex";
  #endif

  int r = SCIP_OKAY;
  int itercount_primal, itercount_dual,gotbasicsol;
  int maxiter;
  MSKprostae prosta;
  MSKsolstae solsta;
  double pobj,dobj;

  #if FORCE_MOSEK_LOG
  if( optimizecount > WRITE_ABOVE )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_ON) );

    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_LOG_FREQ,
                                1) );
  }
  else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_OFF) );
  }
  #else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_OFF) );
  }
  #endif

  MOSEK_CALL( MSK_solutiondef(lpi->task,
                              MSK_SOL_BAS,
                              &gotbasicsol) );

  #if MSK_VERSION_MAJOR >= 5 && !FORCE_NO_SING
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_NON_SINGULAR,
                              MSK_OFF) );
  #endif


  #if FORCE_NO_PRESOLVE 
  if( gotbasicsol )
  { 
    #if MSK_VERSION_MAJOR >= 5
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_PRESOLVE_USE,
                                MSK_PRESOLVE_MODE_OFF) );
    #else
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_PRESOLVE_USE,
                                MSK_OFF) );
    #endif
  }
  else
  { 
    #if MSK_VERSION_MAJOR >= 5
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_PRESOLVE_USE,
                                MSK_PRESOLVE_MODE_ON) );
    #else
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_PRESOLVE_USE,
                                MSK_ON) );
    #endif
  }  
  #endif

  #if CONTROL_DUALIZATION > 0
  int old_solveform;
  int nrows, ncols;
  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_SIM_SOLVE_FORM,
                              &old_solveform) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_SOLVE_FORM,
                              MSK_SOLVE_PRIMAL) );

  if(nrows > 2*ncols)
  {
    int optimizer, other_optimizer;

    MOSEK_CALL( MSK_getintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                &optimizer) );

    other_optimizer = (optimizer == MSK_OPTIMIZER_PRIMAL_SIMPLEX ?
                       MSK_OPTIMIZER_DUAL_SIMPLEX : MSK_OPTIMIZER_PRIMAL_SIMPLEX);

    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_SOLVE_FORM,
                                MSK_SOLVE_DUAL) );

    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                other_optimizer) );
  }
  #endif

  #if ALWAYS_SOLVE_PRIMAL > 0
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_SOLVE_FORM,
                              MSK_SOLVE_PRIMAL) );
  #endif

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_SIM_MAX_ITERATIONS,
                              &maxiter) );

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  #if !SCIP_CONTROLS_PRICING
  {
    int optimizer;

    MOSEK_CALL( MSK_getintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                &optimizer) );

    if( maxiter < 100000 && optimizer == MSK_OPTIMIZER_PRIMAL_SIMPLEX )
    {
      /* We control pricing, so raise max iter */
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  5*maxiter) );
    }
  }
  #elif MSK_VERSION_MAJOR >= 5
  if( gotbasicsol && maxiter < 2000000 )
  {
    int optimizer;

    /* Since max iter often is set, we switch off restricted pricing */
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION,
                                0) );  


    MOSEK_CALL( MSK_getintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                &optimizer) );

    if( optimizer == MSK_OPTIMIZER_PRIMAL_SIMPLEX )
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  2*maxiter) );
  }
  #endif

  if( FORCE_NO_MAXITER )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_MAX_ITERATIONS,
                                2000000000) );
  }


  #if DEBUG_CHECK_DATA > 0
  {
    char begin[] = "Begin optimize with simplex";
    scip_checkdata(lpi,
                   begin);
  }
  #endif

  #if FORCE_MOSEK_SUMMARY
  if( optimizecount > WRITE_ABOVE )
    MOSEK_CALL( MSK_solutionsummary(lpi->task,MSK_STREAM_LOG) );
  #endif

  #if FORCE_CONTROL_SCALING > 0
  if( gotbasicsol )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_SCALING,
                                MSK_SCALING_NONE) );
  }
  #endif

  if( forcescaling )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_SCALING,
                                MSK_SCALING_AGGRESSIVE) );
  }

  MOSEK_CALL( filterTRMrescode(&lpi->termcode, 
                               MSK_optimize(lpi->task)) );

  if( !forcescaling && lpi->termcode == MSK_RES_TRM_MAX_NUM_SETBACKS )
  {    
    forcescaling = 1;

    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_SCALING,
                                MSK_SCALING_AGGRESSIVE) );

    MOSEK_CALL( filterTRMrescode(&lpi->termcode, 
                                 MSK_optimize(lpi->task)) );
  }

  #if FORCE_MOSEK_SUMMARY
  if( optimizecount > WRITE_ABOVE )
    MOSEK_CALL( MSK_solutionsummary(lpi->task,MSK_STREAM_LOG) );
  #endif

  #if DEBUG_CHECK_DATA > 0
  {
    char end[] = "End optimize with simplex";
    scip_checkdata(lpi,
                   end);
  }
  #endif

  MOSEK_CALL( MSK_getintinf(lpi->task,
                            MSK_IINF_SIM_PRIMAL_ITER,
                            &itercount_primal) );

  MOSEK_CALL( MSK_getintinf(lpi->task,
                            MSK_IINF_SIM_DUAL_ITER,
                            &itercount_dual) );

  lpi->itercount = itercount_primal + itercount_dual;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 &pobj,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 &dobj,
                                 NULL,
                                 NULL,
                                 NULL) );

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS || DEBUG_PRINT_STAT
  printf("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, "
          "objval = %g : %g, iter = %d+%d\n",
          maxiter,lpi->termcode,prosta,solsta,
          pobj,dobj,itercount_primal,itercount_dual);
  fflush(stdout);
  #endif

  /*printf("Iter dual %d primal %d\n",itercount_dual,itercount_primal);*/

  switch (solsta)
  {
    case MSK_SOL_STA_OPTIMAL:
    case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
    case MSK_SOL_STA_PRIM_FEAS:
    case MSK_SOL_STA_DUAL_FEAS:
    case MSK_SOL_STA_PRIM_INFEAS_CER:
    case MSK_SOL_STA_DUAL_INFEAS_CER:
        break;
    case MSK_SOL_STA_UNKNOWN:
    case MSK_SOL_STA_NEAR_OPTIMAL:
    case MSK_SOL_STA_NEAR_PRIM_FEAS:
    case MSK_SOL_STA_NEAR_DUAL_FEAS:
    case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      #if SHOW_WARNINGS && !FORCE_SILENCE
      printf("WARNING Simplex[%d] returned solsta = %d\n",optimizecount,solsta);
      #endif

      if (lpi->termcode == MSK_RES_OK)
          lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      break;
    case MSK_SOL_STA_INTEGER_OPTIMAL:
    case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
    default:
      #if SHOW_ERRORS && !FORCE_SILENCE
      printf("ERROR! Simplex[%d] returned solsta = %d\n",optimizecount,solsta);
      #endif

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      RAISE_SCIP_ERROR( SCIP_LPERROR );
  }
  switch (prosta)
  {
    case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
    case MSK_PRO_STA_PRIM_FEAS:
    case MSK_PRO_STA_DUAL_FEAS:
    case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
    case MSK_PRO_STA_PRIM_INFEAS:
    case MSK_PRO_STA_DUAL_INFEAS:
        break;
    case MSK_PRO_STA_UNKNOWN:
    case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
    case MSK_PRO_STA_NEAR_PRIM_FEAS:
    case MSK_PRO_STA_NEAR_DUAL_FEAS:
    case MSK_PRO_STA_ILL_POSED:
    case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      #if SHOW_WARNINGS && !FORCE_SILENCE
      printf("WARNING Simplex[%d] returned prosta = %d\n",optimizecount,prosta);
      #endif

      if (lpi->termcode == MSK_RES_OK)
          lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      break;
    default:
      #if SHOW_ERRORS && !FORCE_SILENCE
      printf("ERROR! Simplex[%d] returned prosta = %d\n",optimizecount,prosta);
      #endif

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      RAISE_SCIP_ERROR( SCIP_LPERROR );
  }

  if (       solsta == MSK_SOL_STA_OPTIMAL
          && fabs(dobj)+fabs(dobj) > 1.0e-6
          && fabs(pobj-dobj)>0.0001*(fabs(pobj)+fabs(dobj)))
  {
    #if SHOW_ERRORS && !FORCE_SILENCE
    printf("ERROR! Simplex[%d] returned optimal solution with different objvals %g != %g reldiff %.2g%%\n",
            optimizecount,pobj,dobj,100*fabs(pobj-dobj)/MAX(fabs(pobj),fabs(dobj)));
    #endif

    #if ASSERT_ON_WARNING
    assert(0);
    #endif
  }

  if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
  {
    if (solsta != MSK_SOL_STA_DUAL_FEAS &&
        solsta != MSK_SOL_STA_OPTIMAL &&
        solsta != MSK_SOL_STA_PRIM_AND_DUAL_FEAS)
    {
      #if SHOW_ERRORS && !FORCE_SILENCE
      printf("ERROR! [%d] Terminated on obj range without dual feasible solsta.\n",
              optimizecount);
      #endif

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      SCIP_CALL_W_R( SCIPlpiSolveBarrier(lpi,true) );
    }
    else
    {
      scipmskobjsen objsen;
      double bound;

      #if MSK_VERSION_MAJOR >= 5
      MOSEK_CALL( MSK_getobjsense(lpi->task,
                                  &objsen) );
      #else
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_OBJECTIVE_SENSE,
                                  &objsen) );
      #endif


      #if MSK_VERSION_MAJOR >= 5
      if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
      #else
      if (objsen == MSK_OBJECTIVE_SENSE_MIN)
      #endif
      {
        MOSEK_CALL( MSK_getdouparam(lpi->task,
                                    MSK_DPAR_UPPER_OBJ_CUT,
                                    &bound) );

        if (0.0000001*(fabs(bound)+fabs(dobj)) < bound-dobj)
        {
          #if SHOW_ERRORS && !FORCE_SILENCE
          printf("ERROR! [%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount,dobj,bound);
          #endif
 
          #if ASSERT_ON_WARNING
          assert(0);
          #endif

          SCIP_CALL_W_R( SCIPlpiSolveBarrier(lpi,true) );
        }
      }
      else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
      {
        MOSEK_CALL( MSK_getdouparam(lpi->task,
                                    MSK_DPAR_LOWER_OBJ_CUT,
                                    &bound) );

        if (0.0000001*(fabs(bound)+fabs(dobj)) < dobj-bound)
        {
          #if SHOW_ERRORS && !FORCE_SILENCE
          printf("ERROR! [%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount,dobj,bound);
          #endif

          #if ASSERT_ON_WARNING
          assert(0);
          #endif

          SCIP_CALL_W_R( SCIPlpiSolveBarrier(lpi,true) );
        }
      }
    }
  }

  if (maxiter >= 2000000000)
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_SIM_MAX_ITERATIONS,
                                maxiter) );

    if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
    {
      #if SHOW_WARNINGS && !FORCE_SILENCE
      printf("WARNING: Simplex[%d] failed to terminate in 10000 iterations, "
              "switching to interior point\n",optimizecount);
      #endif

      #if ASSERT_ON_WARNING
      assert(0);
      #endif

      SCIP_CALL_W_R( SCIPlpiSolveBarrier(lpi,true) );
    }
  }

  #if DEBUG_DO_INTPNT_FEAS_CHECK
  if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
  {
    printf("Checking infeasibility[%d]... ",optimizecount);

    SCIP_CALL_W_R( SCIPlpiSolveBarrier(lpi,true) );

    MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                   MSK_SOL_BAS,
                                   &prosta,
                                   &solsta,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL) );

    if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
    {
      printf("ok\n");
    }
    else
    {
      printf("wrong [%d] prosta = %d, solsta = %d\n",optimizecount,prosta,solsta);
    }
  }
  #endif


  #if DEBUG_PRINT_STAT > 0
  printf("Max iter stat    : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
         optimizecount,numstrongbranchmaxiterup,numstrongbranchmaxiterdo,numprimalmaxiter,numdualmaxiter);
  printf("Objcut iter stat : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
         optimizecount,numstrongbranchobjup,numstrongbranchobjdo,numprimalobj,numdualobj);
  #endif

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiSolvePrimal";
  #endif

  int r = SCIP_OKAY,gotbasicsol;
  optimizecount++;

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  printf("Calling SCIPlpiSolvePrimal[%d] (%d) ",optimizecount,lpi->lpid);
  fflush(stdout);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  MOSEK_CALL( MSK_solutiondef(lpi->task,
                              MSK_SOL_BAS,
                              &gotbasicsol) );

  if( gotbasicsol || FORCE_FREE_SIMPLEX == 0 )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                MSK_OPTIMIZER_PRIMAL_SIMPLEX) );
  }
  else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                MSK_OPTIMIZER_FREE_SIMPLEX) );
  }

  #if WRITE_PRIMAL > 0
  if( optimizecount > WRITE_ABOVE )
  {
    char fname[40];
    sprintf(fname,"primal_%d.mbt",optimizecount);
    printf("\nWriting mbt %s\n",fname);
    /*MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
    MSK_writedata(lpi->task,fname);
  }
  #endif

  SCIP_CALL_W_R( SolveWSimplex(lpi) );

  if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
  {
    MSKprostae prosta;
    MSKsolstae solsta;

    MOSEK_CALL(MSK_getsolutioninf ( lpi->task, MSK_SOL_BAS, &prosta, &solsta,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                         NULL, NULL ));


    if( solsta != MSK_SOL_STA_PRIM_FEAS )
    {
      SCIP_CALL_W_R( SolveWSimplex(lpi) );
    }
  }

  if ( lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS )
  {
    MSKprostae prosta;
    MSKsolstae solsta;

    MOSEK_CALL(MSK_getsolutioninf ( lpi->task, MSK_SOL_BAS, &prosta, &solsta,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                         NULL, NULL ));


    if( solsta != MSK_SOL_STA_PRIM_FEAS )
    {
      int maxiter;

      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  &maxiter) );

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  maxiter+100) );

      SCIP_CALL_W_R( SolveWSimplex(lpi) );
    }
  }

  if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
   ++numprimalobj;

  if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
   ++numprimalmaxiter;


  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY,gotbasicsol;
  optimizecount++;

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  printf("Calling SCIPlpiSolveDual[%d] (%d) ",optimizecount,lpi->lpid);
  fflush(stdout);
  #endif


  MOSEK_CALL( MSK_solutiondef(lpi->task,
                              MSK_SOL_BAS,
                              &gotbasicsol) );

  if( gotbasicsol || FORCE_FREE_SIMPLEX == 0 )
  {
    SCIP_Real bound;
    scipmskobjsen objsen;

    #if MSK_VERSION_MAJOR >= 5
    MOSEK_CALL( MSK_getobjsense(lpi->task,
                                &objsen) );
    #else
    MOSEK_CALL( MSK_getintparam(lpi->task,
                                MSK_IPAR_OBJECTIVE_SENSE,
                                &objsen) );
    #endif


    #if MSK_VERSION_MAJOR >= 5 && !FORCE_NO_SAVE_LU
    if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
    {
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_UPPER_OBJ_CUT,
                                  &bound) );
    }
    else 
    {
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_LOWER_OBJ_CUT,
                                  &bound) );
    }

    if( fabs(bound) < 1.0e+20  )
    {
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_SAVE_LU,
                                  MSK_ON) );
    }
    #endif

    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                MSK_OPTIMIZER_DUAL_SIMPLEX) );
  }
  else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_OPTIMIZER,
                                MSK_OPTIMIZER_FREE_SIMPLEX) );
  }

  #if WRITE_DUAL > 0
  if( optimizecount > WRITE_ABOVE )
  {
    char fname[40];
    sprintf(fname,"dual_%d.mbt",optimizecount);
    printf("\nWriting mbt %s\n",fname);
    MSK_writedata(lpi->task,fname);
  }
  #endif

  SCIP_CALL_W_R( SolveWSimplex(lpi) );

  if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
  {
    MSKprostae prosta;
    MSKsolstae solsta;

    MOSEK_CALL(MSK_getsolutioninf ( lpi->task, MSK_SOL_BAS, &prosta, &solsta,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                         NULL, NULL ));


    if( solsta != MSK_SOL_STA_DUAL_FEAS )
    {
      SCIP_CALL_W_R( SolveWSimplex(lpi) );
    }
  }

  if ( lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS )
  {
    MSKprostae prosta;
    MSKsolstae solsta;

    MOSEK_CALL(MSK_getsolutioninf ( lpi->task, MSK_SOL_BAS, &prosta, &solsta,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                         NULL, NULL ));


    if( solsta != MSK_SOL_STA_DUAL_FEAS )
    {
      int maxiter;

      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  &maxiter) );

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  maxiter+100) );

      SCIP_CALL_W_R( SolveWSimplex(lpi) );
    }
  }

  if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
   ++numdualobj;

  if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
   ++numdualmaxiter;


  CLEANUP:

  return r;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(SCIP_LPI*             lpi,                 /**< LP interface structure */
                                 SCIP_Bool             crossover            /**< perform crossover */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiSolveBarrier";
  #endif

  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  MSKprostae prosta;
  MSKsolstae solsta;
  #endif

  optimizecount++;

  #if FORCE_MOSEK_LOG
  if( optimizecount > WRITE_ABOVE )
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_ON) );
  }
  else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_OFF) );
  }
  #else
  {
    MOSEK_CALL( MSK_putintparam(lpi->task,
                                MSK_IPAR_LOG,
                                MSK_OFF) );
  }
  #endif


  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  printf("Calling SCIPlpiSolveBarrier[%d] (%d) ",optimizecount,lpi->lpid);
  fflush(stdout);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_INTPNT_BASIS,
                              crossover ? MSK_BI_ALWAYS : MSK_BI_NEVER) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_OPTIMIZER,
                              MSK_OPTIMIZER_INTPNT) );


  #if WRITE_INTPNT > 0
  if( optimizecount > WRITE_ABOVE )
  {
    char fname[40];
    sprintf(fname,"intpnt_%d.mbt",optimizecount);
    printf("\nWriting mbt %s\n",fname);
    /*MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
    MSK_writedata(lpi->task,fname);
  }
  #endif

  MOSEK_CALL( filterTRMrescode(&lpi->termcode, 
                               MSK_optimize(lpi->task)) );

  if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
   ++numdualmaxiter;

  MOSEK_CALL( MSK_getintinf(lpi->task,
                            MSK_IINF_INTPNT_ITER,
                            &lpi->itercount) );

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  printf("termcode = %d, prosta = %d, solsta = %d, iter = %d\n",
          lpi->termcode,prosta,solsta,lpi->itercount);
  fflush(stdout);
  #endif

  CLEANUP:

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  return r;
}

/** performs strong branching iterations on all candidates */
SCIP_RETCODE SCIPlpiStrongbranch(SCIP_LPI*             lpi,                /**< LP interface structure */
                                 int                   col,                /**< column to apply strong branching on */
                                 SCIP_Real             psol,               /**< current primal solution value of column */
                                 int                   itlim,              /**< iteration limit for strong branchings */
                                 SCIP_Real*            down,               /**< stores dual bound after branching column down */
                                 SCIP_Real*            up,                 /**< stores dual bound after branching column up */
                                 SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                                                            *   otherwise, it can only be used as an estimate value */
                                 SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                                                            *   otherwise, it can only be used as an estimate value */
                                 int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
  #if DEBUG_CHECK_DATA > 0
  char where[] = "SCIPlpiStrongbranch";
  #endif

  int r = SCIP_OKAY;
  scipmskobjsen objsen;
  int olditerlim;
  int oldselection;

  #if MSK_VERSION_MAJOR >= 5
  int oldhotstart;
  #else
  int oldhotstart;
  #endif

  double bound;
  int ncols, nrows;
  MSKboundkeye bkx;
  double blx, bux, newub, newlb;

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  printf("Calling SCIPlpiStrongbranch (%d)\n",lpi->lpid);
  fflush(stdout);
  #endif

  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  STD_ASSERT;

  #if SHOW_WARNINGS && !FORCE_SILENCE
  if (lpi->termcode != MSK_RES_OK)
    printf("SB Warning: Previous termcode is %d\n",lpi->termcode);
  #endif

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  SCIP_CALL_W_R( getbase(lpi,
                         ncols,
                         nrows) );

  #if MSK_VERSION_MAJOR >= 5 && !FORCE_NO_SAVE_LU
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_SAVE_LU,
                              MSK_ON) );
  #endif

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_getobjsense(lpi->task,
                              &objsen) );
  #else
  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_OBJECTIVE_SENSE,
                              &objsen) );
  #endif

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_SIM_MAX_ITERATIONS,
                              &olditerlim) );

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_SIM_DUAL_SELECTION,
                              &oldselection) );

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_SIM_HOTSTART,
                              &oldhotstart) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_MAX_ITERATIONS,
                              itlim) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_DUAL_SELECTION,
                              STRONGBRANCH_PRICING) );

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_HOTSTART,
                              MSK_SIM_HOTSTART_STATUS_KEYS) );
  #endif

  #if MSK_VERSION_MAJOR >= 5
  if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
  #else
  if (objsen == MSK_OBJECTIVE_SENSE_MIN)
  #endif
  {
    MOSEK_CALL( MSK_getdouparam(lpi->task,
                                MSK_DPAR_UPPER_OBJ_CUT,
                                &bound) );
  }
  else // objsen == MSK_OBJECTIVE_SENSE_MAX
  {
    MOSEK_CALL( MSK_getdouparam(lpi->task,
                                MSK_DPAR_LOWER_OBJ_CUT,
                                &bound) );
  }

  MOSEK_CALL( MSK_getbound(lpi->task,
                           MSK_ACC_VAR,
                           col,
                           &bkx,
                           &blx,
                           &bux) );

  *iter = 0;

  newub = EPSCEIL(psol-1.0, 1e-06);

  if (newub < blx - 0.5) // infeasible
  {
    *down = bound;
    *downvalid = true;
  }
  else
  {
    int newbk;

    if (IS_NEGINF(blx))
      newbk = MSK_BK_UP;
    else if (EPSEQ(blx,newub,1.0e-6))
    {
      newbk = MSK_BK_FX;
      newub = blx;
    }
    else
      newbk = MSK_BK_RA;

    MOSEK_CALL( MSK_putbound(lpi->task,
                             MSK_ACC_VAR,
                             col,
                             newbk,
                             blx,
                             newub) );


    SCIP_CALL_W_R( SCIPlpiSolveDual(lpi) );


    *iter += lpi->itercount;

    if (SCIPlpiIsStable(lpi))
        *downvalid = true;
    else
        *downvalid = false;

    if (SCIPlpiExistsPrimalRay(lpi))
    {
      #if SHOW_ERRORS && !FORCE_SILENCE
      printf("SB ERROR: Lp [%d] is dual infeasible\n",optimizecount);
      #endif

      *down = -1e20;
      *downvalid = false;
    }
    else if (SCIPlpiExistsDualRay(lpi))
    {
      *down = bound;
    }
    else
    {
      SCIP_Bool dfeas;

      SCIP_CALL_W_R( SCIPlpiGetSolFeasibility(lpi,
                                              NULL,
                                              &dfeas) );

      if (!dfeas)
      {
        #if SHOW_ERRORS && !FORCE_SILENCE
        printf("SB ERROR: Lp [%d] is not dual feasible\n",optimizecount);
        getchar();
        #endif

        *down = -1e20;
        *downvalid = false;
      }
      else
      {
        MOSEK_CALL( MSK_getdualobj(lpi->task,
                                   MSK_SOL_BAS,
                                   down) );
      }
    }

    if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
     ++numstrongbranchobjup;

    if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
     ++numstrongbranchmaxiterup;
  }

  // Reset basis solution before doing the up branch

  MOSEK_CALL( MSK_putbound(lpi->task,
                           MSK_ACC_VAR,
                           col,
                           bkx,
                           blx,
                           bux) );

  SCIP_CALL_W_R( setbase(lpi) );

  newlb = EPSFLOOR(psol+1.0, 1e-06);
  if (newlb > bux + 0.5) // infeasible
  {
    *up = bound;
    *upvalid = true;
  }
  else
  {
    int newbk;

    if (IS_POSINF(bux))
        newbk = MSK_BK_LO;
    else if (EPSEQ(bux,newlb,1.0e-6))
    {
      newbk = MSK_BK_FX;
      newlb = bux;
    }
    else
      newbk = MSK_BK_RA;

    MOSEK_CALL( MSK_putbound(lpi->task,
                             MSK_ACC_VAR,
                             col,
                             newbk,
                             newlb,
                             bux) );

    SCIP_CALL_W_R( SCIPlpiSolveDual(lpi) );

    *iter += lpi->itercount;

    if (SCIPlpiIsStable(lpi))
      *upvalid = true;
    else
      *upvalid = false;

    if (SCIPlpiExistsPrimalRay(lpi))
    {
      *up = -1e20;
      *upvalid = false;
    }
    else if (SCIPlpiExistsDualRay(lpi))
    {
      *up = bound;
    }
    else
    {
      SCIP_Bool dfeas;

      SCIP_CALL_W_R( SCIPlpiGetSolFeasibility(lpi,
                                              NULL,
                                              &dfeas) );

      if (!dfeas)
      {
        #if SHOW_ERRORS && !FORCE_SILENCE
        printf("SB ERROR: Lp [%d] is not dual feasible\n",optimizecount);
        getchar();
        #endif

        *up = -1e20;
        *upvalid = false;
      }
      else
      {
        MOSEK_CALL( MSK_getdualobj(lpi->task,
                                   MSK_SOL_BAS,
                                   up) );
      }
    }

    if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
     ++numstrongbranchobjdo;

    if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
     ++numstrongbranchmaxiterdo;
  }

  MOSEK_CALL( MSK_putbound(lpi->task, 
                           MSK_ACC_VAR,
                           col,
                           bkx,
                           blx,
                           bux) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_MAX_ITERATIONS,
                              olditerlim) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_DUAL_SELECTION,
                              oldselection) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_HOTSTART,
                              oldhotstart) );

  #if MSK_VERSION_MAJOR >= 5
  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_SIM_SAVE_LU,
                              MSK_OFF) );
  #endif

  SCIP_CALL_W_R( setbase(lpi) );

  lpi->termcode = MSK_RES_OK;
  lpi->itercount = 0;

  CLEANUP:
  #if DEBUG_CHECK_DATA > 0
  scip_checkdata(lpi,
                 where);
  #endif

  #if DEBUG_PRINT_CALLS || DEBUG_PRINT_OPT_CALLS
  printf("End SCIPlpiStrongbranch (%d)\n",lpi->lpid);
  fflush(stdout);
  #endif

  return r;
}


/*
 * Solution Information Methods
 */


/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiWasSolved (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(SCIP_LPI*             lpi,                /**< LP interface structure */
                                      SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
                                      SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;
  SCIP_Bool pfeas, dfeas;
  pfeas = dfeas = false;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetSolFeasibility (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  switch (solsta)
  {
    case MSK_SOL_STA_OPTIMAL:
    case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      pfeas = dfeas = true;
      break;
    case MSK_SOL_STA_PRIM_FEAS:
      pfeas = true;
      break;
    case MSK_SOL_STA_DUAL_FEAS:
      dfeas = true;
      break;
    case MSK_SOL_STA_UNKNOWN:
    case MSK_SOL_STA_NEAR_OPTIMAL:
    case MSK_SOL_STA_NEAR_PRIM_FEAS:
    case MSK_SOL_STA_NEAR_DUAL_FEAS:
    case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
    case MSK_SOL_STA_PRIM_INFEAS_CER:
    case MSK_SOL_STA_DUAL_INFEAS_CER:
    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
    case MSK_SOL_STA_INTEGER_OPTIMAL:
    case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
      break;
    default:
      RAISE_SCIP_ERROR( SCIP_LPERROR );
  }

  if (primalfeasible != NULL) 
    *primalfeasible = pfeas;

  if (dualfeasible   != NULL) 
    *dualfeasible   = dfeas;

  CLEANUP:

  return r;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiExistsPrimalRay (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL, 
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (   solsta == MSK_SOL_STA_DUAL_INFEAS_CER
          || prosta == MSK_PRO_STA_DUAL_INFEAS
          || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiHasPrimalRay (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (solsta == MSK_SOL_STA_DUAL_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  return false;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  return SCIPlpiExistsDualRay(lpi);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiIsPrimalFeasible (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (prosta == MSK_PRO_STA_PRIM_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiExistsDualRay (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (   solsta == MSK_SOL_STA_PRIM_INFEAS_CER
          || prosta == MSK_PRO_STA_PRIM_INFEAS
          || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiHasDualRay (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (solsta == MSK_SOL_STA_PRIM_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  return false;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  return SCIPlpiExistsPrimalRay(lpi);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiIsDualFeasible (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (prosta == MSK_PRO_STA_DUAL_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiIsOptimal (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  STD_ASSERT;

  return (   lpi->termcode == MSK_RES_OK
          || lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS
          || lpi->termcode == MSK_RES_TRM_MAX_TIME
          || lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  STD_ASSERT;

  return lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  STD_ASSERT;

  return lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  STD_ASSERT;

  return lpi->termcode == MSK_RES_TRM_MAX_TIME;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  int r = SCIP_OKAY;
  MSKprostae prosta;
  MSKsolstae solsta;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetInternalStatus (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolutioninf(lpi->task,
                                 MSK_SOL_BAS,
                                 &prosta,
                                 &solsta,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL) );

  CLEANUP:

  assert(r == SCIP_OKAY);

  return solsta;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(SCIP_LPI*             lpi,                /**< LP interface structure */
                                      SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiIgnoreInstability (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  *success = false;

  return r;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetObjval (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getprimalobj(lpi->task,
                               MSK_SOL_BAS,
                               objval) );

  /* TODO: tjek lighed med dual objektiv i de fleste tilfaelde. */

  CLEANUP:

  return r;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiGetSol(SCIP_LPI*             lpi,                /**< LP interface structure */
                           SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
                           SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
                           SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
                           SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
                           SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
  int r = SCIP_OKAY;
  int ncols, i;
  double *sux = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetSol (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  if (redcost)
  {
    MOSEK_CALL( MSK_getnumvar(lpi->task,
                              &ncols) );

    CHECK_ALLOC( sux = (double *) malloc(sizeof(double) * ncols) );
  }

  MOSEK_CALL( MSK_getsolution(lpi->task,
                              MSK_SOL_BAS,
                              NULL,
                              NULL, 
                              NULL,
                              NULL,
                              NULL,
                              activity,
                              primsol,
                              dualsol,
                              NULL,
                              NULL,
                              redcost,
                              sux,
                              NULL) );

  if (redcost) 
    for (i=0; i<ncols; i++)
    {
      redcost[i] -= sux[i];
    }

  CLEANUP:
  FREE_IF(sux);

  return r;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(SCIP_LPI*             lpi,                /**< LP interface structure */
                                 SCIP_Real*            ray                 /**< primal ray */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetPrimalRay (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolution(lpi->task,
                              MSK_SOL_BAS,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              ray,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL) );

  CLEANUP:

  return r;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(SCIP_LPI*             lpi,                /**< LP interface structure */
                                  SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetDualfarkas (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getsolution(lpi->task,
                              MSK_SOL_BAS,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              dualfarkas,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL) );

  CLEANUP:
  return r;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(SCIP_LPI*             lpi,                /**< LP interface structure */
                                  int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
  int r = SCIP_OKAY;
  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetIterations (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  *iterations = lpi->itercount;

  return r;
}




/*
 * LP Basis Methods
 */

static
SCIP_RETCODE convertstat_mosek2scip(SCIP_LPI*   lpi, 
                                    MSKaccmodee acc, 
                                    MSKstakeye* sk, 
                                    int         n, 
                                    int*        stat)
{
  int r = SCIP_OKAY;
  int i;

  for (i=0; i<n; i++)
  {
    double sl, su;
    switch (sk[i])
    {
      case MSK_SK_BAS:
        stat[i] = SCIP_BASESTAT_BASIC;
        break;
      case MSK_SK_SUPBAS:
        stat[i] = SCIP_BASESTAT_ZERO;
        break;
      case MSK_SK_FIX:
        MOSEK_CALL( MSK_getsolutioni(lpi->task,
                                     acc,
                                     i,
                                     MSK_SOL_BAS,
                                     NULL,
                                     NULL,
                                     &sl,
                                     &su,
                                     NULL) );

        if (sl < su) /* Negative reduced cost */
            stat[i] = SCIP_BASESTAT_UPPER;
        else
            stat[i] = SCIP_BASESTAT_LOWER;
        break;
      case MSK_SK_UNK:
        stat[i] = SCIP_BASESTAT_LOWER;
        break;
      case MSK_SK_INF:
        stat[i] = SCIP_BASESTAT_LOWER;
        break;
      case MSK_SK_LOW:
        stat[i] = SCIP_BASESTAT_LOWER;
        break;
      case MSK_SK_UPR:
        stat[i] = SCIP_BASESTAT_UPPER;
        break;
    }
  }

  CLEANUP:

  return r;
}

static void convertstat_scip2mosek(int        *stat, 
                                   int        n, 
                                   MSKstakeye *resstat)
{
  int i;
  for (i=0; i<n; i++)
  {
    switch (stat[i])
    {
      case SCIP_BASESTAT_LOWER:
        resstat[i] = MSK_SK_LOW;
        break;
      case SCIP_BASESTAT_BASIC:
        resstat[i] = MSK_SK_BAS;
        break;
      case SCIP_BASESTAT_UPPER:
        resstat[i] = MSK_SK_UPR;
        break;
      case SCIP_BASESTAT_ZERO:
        resstat[i] = MSK_SK_SUPBAS;
        break;
    }
  }
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int*                  cstat,              /**< array to store column basis status, or NULL */
                            int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
  int r = SCIP_OKAY;
  int nrows, ncols;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBase (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  SCIP_CALL_W_R( getbase(lpi,
                         ncols,
                         nrows) );

  if (cstat) 
    SCIP_CALL_W_R( convertstat_mosek2scip(lpi,
                                          MSK_ACC_VAR,
                                          lpi->skx,
                                          ncols,
                                          cstat) );

  if (rstat) 
    SCIP_CALL_W_R( convertstat_mosek2scip(lpi,
                                          MSK_ACC_CON,
                                          lpi->skc,
                                          nrows,
                                          rstat) );

  CLEANUP:

  return r;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(SCIP_LPI*             lpi,                /**< LP interface structure */
                            int*                  cstat,              /**< array with column basis status */
                            int*                  rstat               /**< array with row basis status */
   )
{
  int r = SCIP_OKAY;
  int nrows, ncols;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiSetBase (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  SCIP_CALL_W_R( ensureStateMem(lpi,
                                ncols,
                                nrows) );

  convertstat_scip2mosek(cstat,
                         ncols,
                         lpi->skx);

  convertstat_scip2mosek(rstat,
                         nrows,
                         lpi->skc);

  SCIP_CALL_W_R( setbase(lpi) );

  CLEANUP:

  return r;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(SCIP_LPI*             lpi,                /**< LP interface structure */
                                int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
  int r = SCIP_OKAY;
  int i, nrows;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBasisInd (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  MOSEK_CALL( MSK_initbasissolve(lpi->task,
                                 bind) );

  for (i=0; i<nrows; i++)
  {
    if (bind[i] < nrows) // row bind[i] is basic
    {
      bind[i] = -1 - bind[i];
    }
    else // column bind[i]-nrows is basic
    {
      bind[i] = bind[i] - nrows;
    }
  }

  CLEANUP:

  return r;
}

SCIP_RETCODE SCIPlpiGetBInvCol(SCIP_LPI*             lpi,                /**< LP interface structure */
                               int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                                                          *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                                                          *   B^-1 column numbers to the row and column numbers of the LP!
                                                                          *   c must be between 0 and nrows-1, since the basis has the size
                                                                          *   nrows * nrows */
                               SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
  int r = SCIP_OKAY;
  int i;
  int nrows, numnz;
  int *sub = NULL;
  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBInvCol (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,&nrows) );
  CHECK_ALLOC( sub = (int *) malloc(sizeof(int) * nrows) );

  for (i=0; i<nrows; i++)
      coef[i] = 0;

  numnz   = 1; 
  sub[0]  = c; 
  coef[c] = 1; // Unit vector e_col

  MOSEK_CALL( MSK_solvewithbasis(lpi->task,
                                 0,
                                 &numnz,
                                 sub,
                                 coef) );

  CLEANUP:
  FREE_IF(sub);

  return r;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(SCIP_LPI*             lpi,                /**< LP interface structure */
                                int                   c,                  /**< column number */
                                SCIP_Real*            coef                /**< vector to return coefficients */
   )
{ 
  /*lint --e{715}*/
  int r = SCIP_OKAY;
  int i;
  int nrows, numnz;
  int *sub = NULL;
  SCIP_Real *val=NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBInvACol (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,&nrows) );
  MOSEK_CALL( MSK_getavecnumnz(lpi->task,MSK_ACC_VAR,c,&numnz) );
  CHECK_ALLOC( sub = (int *) malloc(sizeof(int) * (nrows)) );
  CHECK_ALLOC( val = (SCIP_Real *) malloc(sizeof(SCIP_Real) * (numnz+1)) );

  /* Avoid a bug in scip  */
  MOSEK_CALL( MSK_initbasissolve(lpi->task,
                                 sub) ); /*used as dummy*/

  for (i=0; i<nrows; i++)
      coef[i] = 0;

  MOSEK_CALL( MSK_getavec(lpi->task,
                          MSK_ACC_VAR,
                          c,
                          &numnz,
                          sub,
                          val) );

  for (i=0; i<numnz; i++)
      coef[sub[i]] = val[i];

  MOSEK_CALL( MSK_solvewithbasis(lpi->task,
                                 0,
                                 &numnz,
                                 sub,
                                 coef) );

  CLEANUP:
  FREE_IF(sub);
  FREE_IF(val);

  return r;
}


/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(SCIP_LPI*             lpi,                /**< LP interface structure */
                               int                   row,                /**< row number */
                               SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
  int r = SCIP_OKAY;
  int i;
  int nrows, numnz;
  int *sub = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBInvRow (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  CHECK_ALLOC( sub = (int *) malloc(sizeof(int) * nrows) );

  for (i=0; i<nrows; i++)
    coef[i] = 0;

  numnz     = 1; 
  sub[0]    = row; 
  coef[row] = 1; // Unit vector e_row

  MOSEK_CALL( MSK_solvewithbasis(lpi->task,
                                 1,
                                 &numnz,
                                 sub,
                                 coef) );

  CLEANUP:
  FREE_IF(sub);

  return r;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(SCIP_LPI*             lpi,                /**< LP interface structure */
                                int                   row,                /**< row number */
                                const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
                                SCIP_Real*            val                 /**< vector to return coefficients */
   )
{
  int r = SCIP_OKAY;
  int i,k;
  int nrows, ncols, numnz;
  int *csub = NULL,didalloc=0;
  double *cval = NULL,*binv = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetBInvARow (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  CHECK_ALLOC( csub = (int *) malloc(sizeof(int) * nrows) );

  /* Avoid a bug in scip  */
  MOSEK_CALL( MSK_initbasissolve(lpi->task,
                                 csub) ); /*used as dummy*/

  for(i=0; i<ncols; i++)
    val[i] = 0;

  if( binvrow == NULL )
  {
    didalloc = 1;

    CHECK_ALLOC( binv = (SCIP_Real *) malloc(sizeof(SCIP_Real) * nrows) );

    SCIP_CALL_W_R( SCIPlpiGetBInvRow(lpi,
                                     row,
                                     binv) );
  }
  else
    binv = (SCIP_Real*)binvrow;

  /* binvrow*A */
  for (i=0; i<ncols; i++)
  {
    val[i] = 0;

    MOSEK_CALL( MSK_getavec(lpi->task,
                            MSK_ACC_VAR,
                            i,
                            &numnz,
                            csub,
                            cval) );

 
    for( k = 0; k < numnz; ++k )
      val[i] += binv[csub[k]]*cval[k];
  }

  CLEANUP:
  FREE_IF(csub);

  if( didalloc )
    FREE_IF(binv);

  return r;
}

/*
 * LP State Methods
 */

static SCIP_LPISTATE* allocLpiState(BMS_BLKMEM* blkmem, int ncols, int nrows)
{
  SCIP_LPISTATE* lpistate = NULL;

  if (blkmem)
  {
    BMSallocBlockMemory(blkmem, &lpistate);

    if (lpistate == NULL) 
      return NULL;

    lpistate->ncols = ncols;
    lpistate->nrows = nrows;

    BMSallocBlockMemoryArray(blkmem, &lpistate->skx, COLPACKET_NUM(ncols));
    BMSallocBlockMemoryArray(blkmem, &lpistate->skc, ROWPACKET_NUM(nrows));

    if (lpistate->skc == NULL || lpistate->skx == NULL)
    {
        if (lpistate->skx != NULL) 
        { 
          BMSfreeBlockMemoryArray(blkmem, &lpistate->skx, COLPACKET_NUM(ncols)); 
        }

        if (lpistate->skc != NULL) 
        { 
          BMSfreeBlockMemoryArray(blkmem, &lpistate->skc, ROWPACKET_NUM(nrows)); 
        }

        BMSfreeBlockMemory(blkmem, &lpistate);

        return NULL;
    }
  }
  else
  {
    lpistate = malloc(sizeof(SCIP_LPISTATE));

    if (lpistate == NULL) 
      return NULL;

    lpistate->ncols = ncols;
    lpistate->nrows = nrows;

    lpistate->skx = (COLPACKET *) malloc(sizeof(COLPACKET) * COLPACKET_NUM(ncols));
    lpistate->skc = (ROWPACKET *) malloc(sizeof(ROWPACKET) * ROWPACKET_NUM(nrows));

    if (lpistate->skc == NULL || lpistate->skx == NULL)
    {
      if (lpistate->skx != NULL) 
      { 
        free(lpistate->skx); 
      }

      if (lpistate->skc != NULL) 
      { 
        free(lpistate->skc); 
      }

      free(lpistate);

      return NULL;
    }
  }

  return lpistate;
}

static void freeLpiState(BMS_BLKMEM* blkmem, SCIP_LPISTATE** lpistate)
{
    int nrows, ncols;
    assert(lpistate);
    if (*lpistate == NULL) return;
    ncols = (*lpistate)->ncols;
    nrows = (*lpistate)->nrows;
    if (blkmem)
    {
        BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skx, COLPACKET_NUM(ncols));
        BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skc, ROWPACKET_NUM(nrows));
        BMSfreeBlockMemory(blkmem, lpistate);
    }
    else
    {
        free((*lpistate)->skx);
        free((*lpistate)->skc);
        free(*lpistate);
    }
    *lpistate = NULL;
}

static SCIP_RETCODE handle_singular(SCIP_LPI*   lpi, 
                                    int*        basis, 
                                    MSKrescodee res)
{
  int r = SCIP_OKAY;
  if (res == MSK_RES_ERR_BASIS_SINGULAR)
  {
    SCIP_CALL_W_R( SCIPlpiSolvePrimal(lpi) );

    MOSEK_CALL( MSK_initbasissolve(lpi->task,
                                   basis) );
  }
  else
  {
    MOSEK_CALL( res );
  }

  CLEANUP:

  return r;
}

#if DEBUG_CHECK_STATE
static
SCIP_RETCODE checkState1(SCIP_LPI*   lpi, 
                         int         n, 
                         MSKstakeye* sk, 
                         MSKaccmodee accmode, 
                         char        xc)
{
  int r = SCIP_OKAY;
  int i;

  /* printout for all except LOW, UPR, FIX and BAS with sl[xc]==su[xc] */
  for (i=0; i<n; i++)
  {
    double sl, su;
    switch (sk[i])
    {
      case MSK_SK_UNK:
        printf("STATE[%d]: %c[%d] = unk\n",optimizecount,xc,i);
        break;
      case MSK_SK_BAS:
        MOSEK_CALL( MSK_getsolutioni(lpi->task,
                                     accmode,
                                     i,
                                     MSK_SOL_BAS,
                                     NULL,
                                     NULL,
                                     &sl,
                                     &su,
                                     NULL) );

        if (fabs(sl-su) > DEBUG_CHECK_STATE_TOL)
            printf("STATE[%d]: %c[%d] = bas, sl%c = %g, su%c = %g\n",
                    optimizecount,xc,i,xc,sl,xc,su);
        break;
      case MSK_SK_SUPBAS:
        printf("STATE[%d]: %c[%d] = supbas\n",optimizecount,xc,i);
        break;
      case MSK_SK_LOW:
      case MSK_SK_UPR:
      case MSK_SK_FIX:
        break;
      case MSK_SK_INF:
        printf("STATE[%d]: %c[%d] = inf\n",optimizecount,xc,i);
        break;
      default:
        printf("STATE[%d]: %c[%d] = ????????????????\n",optimizecount,xc,i);
        break;
    }
  }

  CLEANUP:

  return r;
}

static
SCIP_RETCODE checkState(SCIP_LPI* lpi, 
                        int       ncols, 
                        int       nrows)
{
  int r = SCIP_OKAY;

  SCIP_CALL_W_R( checkState1(lpi,
                             ncols,
                             lpi->skx,
                             MSK_ACC_VAR,
                             'x') );

  SCIP_CALL_W_R( checkState1(lpi,
                             nrows,
                             lpi->skc,
                             MSK_ACC_CON,
                             'c') );

  CLEANUP:

  return r;
}
#endif

static
SCIP_RETCODE lpistatePack(SCIP_LPI*      lpi, 
                          SCIP_LPISTATE* lpistate)
{
  int r = SCIP_OKAY;

  #if PACK_LPISTATE
  int *skxi = (int *) lpi->skx;
  int *skci = (int *) lpi->skc;

  assert(sizeof(int) == sizeof(MSKstakeye));

  SCIP_CALL_W_R( convertstat_mosek2scip(lpi,
                                        MSK_ACC_VAR,
                                        lpi->skx,
                                        lpistate->ncols,
                                        skxi) );

  SCIP_CALL_W_R( convertstat_mosek2scip(lpi,
                                        MSK_ACC_CON,
                                        lpi->skc,lpistate->nrows,skci) );

  SCIPencodeDualBit(skxi, lpistate->skx, lpistate->ncols);

  SCIPencodeDualBit(skci, lpistate->skc, lpistate->nrows);
  #else
  int i;
  for (i=0; i<lpistate->ncols; i++)
    lpistate->skx[i] = lpi->skx[i];

  for (i=0; i<lpistate->nrows; i++)
    lpistate->skc[i] = lpi->skc[i];
  #endif

  CLEANUP:

  return r;
}

static
void lpistateUnpack(SCIP_LPISTATE* lpistate, 
                    MSKstakeye*    skx, 
                    MSKstakeye*    skc)
{
  #if PACK_LPISTATE
  int *skxi = (int *) skx;
  int *skci = (int *) skc;

  assert(sizeof(int) == sizeof(MSKstakeye));

  SCIPdecodeDualBit(lpistate->skx, 
                    skxi, 
                    lpistate->ncols);

  SCIPdecodeDualBit(lpistate->skc, 
                    skci, 
                    lpistate->nrows);

  convertstat_scip2mosek(skxi,
                         lpistate->ncols,
                         skx);

  convertstat_scip2mosek(skci,
                         lpistate->nrows,
                         skc);
  #else
  int i;
  for (i=0; i<lpistate->ncols; i++)
    skx[i] = lpistate->skx[i];

  for (i=0; i<lpistate->nrows; i++)
    skc[i] = lpistate->skc[i];
  #endif
}

/** stores LP state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(SCIP_LPI*             lpi,                /**< LP interface structure */
                             BMS_BLKMEM*           blkmem,             /**< block memory */
                             SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
  int r = SCIP_OKAY;
  int gotbasicsol;
  int nrows, ncols;
  int *dummy = NULL;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetState (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  assert(lpistate);

  *lpistate = NULL;

  MOSEK_CALL( MSK_solutiondef(lpi->task,
                              MSK_SOL_BAS,
                              &gotbasicsol) );

  if (!gotbasicsol || SCIPlpiExistsDualRay(lpi)) 
    goto CLEANUP;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );

  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  CHECK_ALLOC( dummy = (int*) malloc(sizeof(int) * nrows) );

  SCIP_CALL_W_R( handle_singular(lpi,dummy,MSK_initbasissolve(lpi->task,
                                                              dummy)) );

  CHECK_ALLOC( *lpistate = allocLpiState(blkmem,
                                         ncols,
                                         nrows) );

  SCIP_CALL_W_R( getbase(lpi,
                         ncols,
                         nrows) );

  #if DEBUG_CHECK_STATE
  SCIP_CALL_W_R( checkState(lpi,
                            ncols,
                            nrows) );
  #endif

  SCIP_CALL_W_R( lpistatePack(lpi,
                              *lpistate) );

  CLEANUP:
  FREE_IF(dummy);

  if (r != SCIP_OKAY) 
    freeLpiState(blkmem, lpistate);

  return r;
}

/** loads LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpiSetState(SCIP_LPI*             lpi,                /**< LP interface structure */
                             BMS_BLKMEM*           blkmem,             /**< block memory */
                             SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
  int r = SCIP_OKAY;
  int nrows, ncols, i;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiSetState (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  if (!lpistate)
  {
    printf("Setting NULL state\n");
    goto CLEANUP;
  }

  if (lpistate->nrows == 0 || lpistate->ncols == 0)
    goto CLEANUP;

  MOSEK_CALL( MSK_getnumcon(lpi->task,
                            &nrows) );


  MOSEK_CALL( MSK_getnumvar(lpi->task,
                            &ncols) );

  assert(lpistate->nrows <= nrows);
  assert(lpistate->ncols <= ncols);

  SCIP_CALL_W_R( ensureStateMem(lpi,ncols,nrows) );

  lpistateUnpack(lpistate,lpi->skx,lpi->skc);

  for (i=lpistate->ncols; i<ncols; i++)
    lpi->skx[i] = MSK_SK_UNK;

  for (i=lpistate->nrows; i<nrows; i++)
    lpi->skc[i] = MSK_SK_BAS;

  SCIP_CALL_W_R( setbase(lpi) );

  CLEANUP:

  return r;
}

/** frees LP state information */
SCIP_RETCODE SCIPlpiFreeState(SCIP_LPI*             lpi,                /**< LP interface structure */
                              BMS_BLKMEM*           blkmem,             /**< block memory */
                              SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiFreeState (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  freeLpiState(blkmem, lpistate);

  return r;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(SCIP_LPI*             lpi,                /**< LP interface structure */
                               SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
  STD_ASSERT;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiHasStateBasis (%d)\n",lpi->lpid);
  #endif

  return lpistate != NULL;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(SCIP_LPI*             lpi,                /**< LP interface structure */
                              const char*           fname               /**< file name */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiReadState (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_readsolution(lpi->task,
                               MSK_SOL_BAS,
                               fname) );

  CLEANUP:

  return r;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(SCIP_LPI*             lpi,                /**< LP interface structure */
                               const char*           fname               /**< file name */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiWriteState (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_writesolution(lpi->task,
                                MSK_SOL_BAS,
                                fname) );

  CLEANUP:

  return r;
}

/*
 * Parameter Methods
 */

#if DEBUG_PRINT_CALLS
static const char* paramname[] = {
    "SCIP_LPPAR_FROMSCRATCH",
    "SCIP_LPPAR_FASTMIP",
    "SCIP_LPPAR_SCALING",
    "SCIP_LPPAR_PRESOLVING",
    "SCIP_LPPAR_PRICING",
    "SCIP_LPPAR_LPINFO",
    "SCIP_LPPAR_FEASTOL",
    "SCIP_LPPAR_DUALFEASTOL",
    "SCIP_LPPAR_BARRIERCONVTOL",
    "SCIP_LPPAR_LOBJLIM",
    "SCIP_LPPAR_UOBJLIM",
    "SCIP_LPPAR_LPITLIM",
    "SCIP_LPPAR_LPTILIM",
    "SCIP_LPPAR_MARKOWITZ", };

static const char* paramty2str(int type)
{
  return paramname[type];
}
#endif

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(SCIP_LPI*             lpi,                /**< LP interface structure */
                              SCIP_LPPARAM          type,               /**< parameter number */
                              int*                  ival                /**< buffer to store the parameter value */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetIntpar (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  switch (type)
  {
    case SCIP_LPPAR_FROMSCRATCH:               /**< solver should start from scratch at next call? */
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_SIM_HOTSTART,
                                  ival) );


      #if MSK_VERSION_MAJOR >= 5
      *ival = (*ival == MSK_SIM_HOTSTART_NONE);
      #else
      *ival = (*ival == MSK_OFF);
      #endif

      break;
    case SCIP_LPPAR_FASTMIP:                   /**< fast mip setting of LP solver */
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
    case SCIP_LPPAR_SCALING:                   /**< should LP solver use scaling? */
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_SIM_SCALING,
                                  ival) );

      *ival = (*ival != MSK_SCALING_NONE);

      break;
    case SCIP_LPPAR_PRESOLVING:                /**< should LP solver use presolving? */
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_PRESOLVE_USE,
                                  ival) );

      #if MSK_VERSION_MAJOR >= 5
      *ival = (*ival != MSK_PRESOLVE_MODE_OFF);
      #else
      *ival = (*ival == MSK_ON);
      #endif

      break;
    case SCIP_LPPAR_PRICING:                   /**< pricing strategy */
      *ival = lpi->pricing;

      break;
    case SCIP_LPPAR_LPINFO:                    /**< should LP solver output information to the screen? */
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_LOG,
                                  ival) );

      *ival = (*ival == MSK_ON);

      break;
    case SCIP_LPPAR_LPITLIM:                   /**< LP iteration limit */
      MOSEK_CALL( MSK_getintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  ival) );

      break;
    default:
        RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );
        break;
  }

  CLEANUP:

  return r;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
  int r = SCIP_OKAY;
  int scaling;

  #if MSK_VERSION_MAJOR >= 5
  #if SCIP_CONTROLS_PRICING
  static int pricing[6] = {
      MSK_SIM_SELECTION_FREE,
      MSK_SIM_SELECTION_FULL,
      MSK_SIM_SELECTION_PARTIAL,
      MSK_SIM_SELECTION_DEVEX,
      MSK_SIM_SELECTION_ASE,
      MSK_SIM_SELECTION_DEVEX,
  };
  #endif
  #else
  static int pricing[6] = {
      MSK_SIM_SELECTION_FREE,
      MSK_SIM_SELECTION_FULL,
      MSK_SIM_SELECTION_PARTIAL,
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_ASE,
      MSK_SIM_SELECTION_ASE,
  };
  #endif

  assert(SCIP_PRICING_AUTO        == 0);
  assert(SCIP_PRICING_FULL        == 1);
  assert(SCIP_PRICING_PARTIAL     == 2);
  assert(SCIP_PRICING_STEEP       == 3);
  assert(SCIP_PRICING_STEEPQSTART == 4);
  assert(SCIP_PRICING_DEVEX == 5);

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiSetIntpar (%d) %s = %d\n",lpi->lpid,paramty2str(type),ival);
  #endif

  STD_ASSERT;

  switch (type)
  {
    case SCIP_LPPAR_FROMSCRATCH:               /**< solver should start from scratch at next call? */
      #if MSK_VERSION_MAJOR >= 5
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_HOTSTART,
                                  ival ? MSK_SIM_HOTSTART_NONE : MSK_SIM_HOTSTART_FREE ) );
      #else
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_HOTSTART,
                                  ival ? MSK_OFF : MSK_ON ) );
      #endif

      break;
    case SCIP_LPPAR_FASTMIP:                   /**< fast mip setting of LP solver */
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
    case SCIP_LPPAR_SCALING:                   /**< should LP solver use scaling? */
      scaling = (ival ? MSK_SCALING_FREE : MSK_SCALING_NONE);

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_SCALING,
                                  scaling) );

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_INTPNT_SCALING,
                                  scaling) );

      break;
    case SCIP_LPPAR_PRESOLVING:                /**< should LP solver use presolving? */
      #if MSK_VERSION_MAJOR >= 5
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_PRESOLVE_USE, 
                                  ival ? MSK_PRESOLVE_MODE_FREE : MSK_PRESOLVE_MODE_OFF) );
      #else
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_PRESOLVE_USE, 
                                  ival ? MSK_ON : MSK_OFF) );
      #endif

      break;
    case SCIP_LPPAR_PRICING:                   /**< pricing strategy */
      assert(ival>=0 && ival<=SCIP_PRICING_DEVEX);

      lpi->pricing = ival;

      #if SCIP_CONTROLS_PRICING
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_PRIMAL_SELECTION,
                                  pricing[ival]) );
      #endif

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_DUAL_SELECTION,
                                  MSK_SIM_SELECTION_FREE) );

      break;
    case SCIP_LPPAR_LPINFO:                    /**< should LP solver output information to the screen? */
      #if FORCE_MOSEK_LOG
      printf("Ignoring log setting!\n");
      #else
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_LOG, 
                                  ival ? MSK_ON : MSK_OFF) );
      #endif

      break;
    case SCIP_LPPAR_LPITLIM:                   /**< LP iteration limit */
      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_SIM_MAX_ITERATIONS,
                                  ival) );

      MOSEK_CALL( MSK_putintparam(lpi->task,
                                  MSK_IPAR_INTPNT_MAX_ITERATIONS,
                                  ival) );

      break;
    default:
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
  }

  CLEANUP:
  return r;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(SCIP_LPI*             lpi,                /**< LP interface structure */
                               SCIP_LPPARAM          type,               /**< parameter number */
                               SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiGetRealpar (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  switch (type)
  {
    case SCIP_LPPAR_FEASTOL:                   /**< feasibility tolerance for primal variables and slacks */
      #if SCIP_CONTROLS_TOLERANCES
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_BASIS_TOL_X,
                                  dval) );
      #endif

      break;
    case SCIP_LPPAR_DUALFEASTOL:               /**< feasibility tolerance for dual variables and reduced costs */
      #if SCIP_CONTROLS_TOLERANCES
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_BASIS_TOL_S,
                                  dval) );
      #endif

      break;
    case SCIP_LPPAR_BARRIERCONVTOL:            /**< convergence tolerance used in barrier algorithm */
      #if SCIP_CONTROLS_TOLERANCES
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_INTPNT_TOL_REL_GAP,
                                  dval) );
      #else
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );
      #endif

      break;
    case SCIP_LPPAR_LOBJLIM:                   /**< lower objective limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_LOWER_OBJ_CUT,
                                  dval) );

      break;
    case SCIP_LPPAR_UOBJLIM:                   /**< upper objective limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_UPPER_OBJ_CUT,
                                  dval) );

      break;
    case SCIP_LPPAR_LPTILIM:                   /**< LP time limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task,
                                  MSK_DPAR_OPTIMIZER_MAX_TIME,
                                  dval) );

      break;
    case SCIP_LPPAR_MARKOWITZ:                 /**< Markowitz tolerance */
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
    default:
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
  }

  CLEANUP:

  return r;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(SCIP_LPI*             lpi,                /**< LP interface structure */
                               SCIP_LPPARAM          type,               /**< parameter number */
                               SCIP_Real             dval                /**< parameter value */
   )
{
  int r = SCIP_OKAY;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiSetRealpar (%d) %s = %g\n",lpi->lpid,paramty2str(type),dval);
  #endif

  STD_ASSERT;

  /* TODO: Limits shouldn't be hardcoded */

  switch (type)
  {
    case SCIP_LPPAR_FEASTOL:                   /**< feasibility tolerance for primal variables and slacks */
      #if SCIP_CONTROLS_TOLERANCES
      if (dval < 1e-10) dval = 1e-10;
        MOSEK_CALL( MSK_putdouparam(lpi->task,
                                    MSK_DPAR_BASIS_TOL_X,
                                    dval) );

      break;
      #endif
    case SCIP_LPPAR_DUALFEASTOL:               /**< feasibility tolerance for dual variables and reduced costs */
      #if SCIP_CONTROLS_TOLERANCES
      if (dval < 1e-10) dval = 1e-10;
        MOSEK_CALL( MSK_putdouparam(lpi->task,
                                    MSK_DPAR_BASIS_TOL_S,
                                    dval) );

      break;
      #endif
    case SCIP_LPPAR_BARRIERCONVTOL:            /**< convergence tolerance used in barrier algorithm */
      #if SCIP_CONTROLS_TOLERANCES
      MOSEK_CALL( MSK_putdouparam(lpi->task,
                                  MSK_DPAR_INTPNT_TOL_REL_GAP,
                                  dval) );
      #else
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );
      #endif

      break;
    case SCIP_LPPAR_LOBJLIM:                   /**< lower objective limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task,
                                  MSK_DPAR_LOWER_OBJ_CUT,
                                  dval) );

      break;
    case SCIP_LPPAR_UOBJLIM:                   /**< upper objective limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task,
                                  MSK_DPAR_UPPER_OBJ_CUT,
                                  dval) );

      break;
    case SCIP_LPPAR_LPTILIM:                   /**< LP time limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task,
                                  MSK_DPAR_OPTIMIZER_MAX_TIME,
                                  dval) );

      break;
    case SCIP_LPPAR_MARKOWITZ:                 /**< Markowitz tolerance */
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
    default:
      RAISE_SCIP_ERROR( SCIP_PARAMETERUNKNOWN );

      break;
  }

  CLEANUP:

  return r;
}


/*
 * Numerical Methods
 */


/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
  STD_ASSERT;

  return MSK_INFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(SCIP_LPI*             lpi,                /**< LP interface structure */
                            SCIP_Real             val
   )
{
  STD_ASSERT;

  return IS_POSINF(val);
}


/*
 * File Interface Methods
 */


/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(SCIP_LPI*             lpi,                /**< LP interface structure */
                           const char*           fname               /**< file name */
   )
{
  int r = SCIP_OKAY;
  int olddataformat;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiReadLP (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_READ_DATA_FORMAT,
                              &olddataformat) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_READ_DATA_FORMAT,
                              MSK_DATA_FORMAT_MBT) );

  MOSEK_CALL( MSK_readdata(lpi->task,
                           fname) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_READ_DATA_FORMAT,
                              olddataformat) );

  CLEANUP:

  return r;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(SCIP_LPI*             lpi,                /**< LP interface structure */
                            const char*           fname               /**< file name */
   )
{
  int r = SCIP_OKAY;
  int olddataformat;

  #if DEBUG_PRINT_CALLS
  printf("Calling SCIPlpiWriteLP (%d)\n",lpi->lpid);
  #endif

  STD_ASSERT;

  MOSEK_CALL( MSK_getintparam(lpi->task,
                              MSK_IPAR_WRITE_DATA_FORMAT,
                              &olddataformat) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_WRITE_DATA_FORMAT,
                              MSK_DATA_FORMAT_MBT) );

  MOSEK_CALL( MSK_writedata(lpi->task,
                            fname) );

  MOSEK_CALL( MSK_putintparam(lpi->task,
                              MSK_IPAR_WRITE_DATA_FORMAT,
                              olddataformat) );

  CLEANUP:

  return r;
}

