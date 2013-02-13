/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpi_msk.c
 * @ingroup LPIS
 * @brief  LP interface for MOSEK
 * @author Bo Jensen
 */

#include <assert.h>

#define MSKCONST const
#include "mosek.h"

#include "scip/lpi.h"
#include "scip/bitencode.h"
#include <string.h>

/* do defines for windows directly her to make the lpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#endif

#define scipmskobjsen MSKobjsensee
#define SENSE2MOSEK(objsen) (((objsen)==SCIP_OBJSEN_MINIMIZE)?(MSK_OBJECTIVE_SENSE_MINIMIZE):(MSK_OBJECTIVE_SENSE_MAXIMIZE))

#define MOSEK_CALL(x)  do                                                                                     \
                       {  /*lint --e{641}*/                                                                   \
                          MSKrescodee _restat_;                                                               \
                          _restat_ = (x);                                                                     \
                          if( (_restat_) != MSK_RES_OK && (_restat_ ) != MSK_RES_TRM_MAX_NUM_SETBACKS )       \
                          {                                                                                   \
                             SCIPerrorMessage("LP Error: MOSEK returned %d\n", (int)_restat_);                \
                             return SCIP_LPERROR;                                                             \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/* this macro is only called in functions returning SCIP_Bool; thus, we return FALSE if there is an error in optimized mode */
#define ABORT_FALSE(x) { int _restat_;                                  \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPerrorMessage("LP Error: MOSEK returned %d\n", (int)_restat_); \
         SCIPABORT();                                                   \
         return FALSE;                                                  \
      }                                                                 \
   }


#define IS_POSINF(x) ((x) >= SCIP_DEFAULT_INFINITY)
#define IS_NEGINF(x) ((x) <= -SCIP_DEFAULT_INFINITY)

static MSKenv_t MosekEnv =           NULL;
static int numlp         =           0;

static int optimizecount            =  0;
static int nextlpid                 =  1;
static int numstrongbranchmaxiterup =  0;
static int numstrongbranchmaxiterdo =  0;
static int numprimalmaxiter         =  0;
static int numdualmaxiter           =  0;
static int numstrongbranchobjup     =  0;
static int numstrongbranchobjdo     =  0;
static int numprimalobj             =  0;
static int numdualobj               =  0;

#define DEBUG_PARAM_SETTING          0
#define DEBUG_PRINT_STAT             0
#define DEBUG_CHECK_DATA             0
#define DEBUG_EASY_REPRODUCE         0
#define DEBUG_DO_INTPNT_FEAS_CHECK   0
#define DEBUG_CHECK_STATE_TOL        1e-5
#define SHOW_ERRORS                  0
#define ASSERT_ON_NUMERICAL_TROUBLES 0
#define ASSERT_ON_WARNING            0
#define FORCE_MOSEK_LOG              0
#define FORCE_MOSEK_SUMMARY          0
#define FORCE_NO_MAXITER             0
#define FORCE_SILENCE                1
#define SETBACK_LIMIT                250
#define SCIP_CONTROLS_PRICING        1
#define SCIP_CONTROLS_TOLERANCES     1
#define STRONGBRANCH_PRICING         MSK_SIM_SELECTION_SE
#define SUPRESS_NAME_ERROR           1
#define WRITE_DUAL                   0
#define WRITE_PRIMAL                 0
#define WRITE_INTPNT                 0
#define WRITE_ABOVE                  0
#define DEGEN_LEVEL                  MSK_SIM_DEGEN_FREE
#define ALWAYS_SOLVE_PRIMAL          1

#if MSK_VERSION_MAJOR >= 7
/** gives problem and solution status for a Mosek Task
 *
 * With Mosek 7.0, the routine MSK_getsolutionstatus was replaced by
 * MSK_getprosta and MSK_getsolsta.
 */
static
MSKrescodee MSK_getsolutionstatus(
   MSKtask_t             task,               /**< Mosek Task */
   MSKsoltypee           whichsol,           /**< for which type of solution a status is requested */
   MSKprostae*           prosta,             /**< buffer to store problem status, or NULL if not needed */
   MSKsolstae*           solsta              /**< buffer to store solution status, or NULL if not needed */
   )
{
   if( prosta != NULL )
   {
      MOSEK_CALL( MSK_getprosta(task, whichsol, prosta) );
   }
   if( solsta != NULL )
   {
      MOSEK_CALL( MSK_getsolsta(task, whichsol, solsta) );
   }

   return MSK_RES_OK;
}
#endif

/**********************************************/

struct SCIP_LPi
{
   MSKtask_t             task;
   MSKrescodee           termcode;
   int                   itercount;
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   int                   lpid;
   int                   skxsize;
   int                   skcsize;
   MSKstakeye*           skx;
   MSKstakeye*           skc;
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

struct SCIP_LPiState
{
   int                   num;
   MSKsolstae            solsta;
   int                   ncols;
   int                   nrows;
   COLPACKET*            skx;
   ROWPACKET*            skc;
};

/** returns the number of packets needed to store column packet information */
static 
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static 
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** create error string */
static
void MSKAPI printstr(
   void*                 handle,             /**< error handle */
   const char*           str                 /**< string that contains string on output */
   )
{  /*lint --e{715}*/
#if SUPRESS_NAME_ERROR && !FORCE_SILENCE
   char errstr[32];
   snprintf(errstr,32,"MOSEK Error %d",MSK_RES_ERR_DUP_NAME);
   if (0 == strncmp(errstr,str,strlen(errstr)))
      return;
#endif

   SCIPdebugMessage("MOSEK: %s",str);
}

#if DEBUG_CHECK_DATA > 0
/** check data */
static SCIP_RETCODE scip_checkdata(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   const char*           functionname        /**< function name */
   )
{
   int i;
   int numcon;
   int numvar;
   int gotbasicsol;
   MSKboundkeye* tbkc;
   MSKboundkeye* tbkx;
   MSKstakeye *tskc;
   MSKstakeye* tskx;
   double* tblc;
   double* tbuc;
   double* tblx;
   double* tbux;

   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   MOSEK_CALL( MSK_getnumvar(lpi->task,&numvar) );
   MOSEK_CALL( MSK_getnumcon(lpi->task,&numcon) );

   /* allocate memory */
   SCIP_ALLOC( BMSallocMemoryArray( &tbkc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tskc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tblc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tbuc, numcon) );

   SCIP_ALLOC( BMSallocMemoryArray( &tbkx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tskx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tblx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tbux, numvar) );

   /* Check bounds */
   if( gotbasicsol )
   {
      MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, tskc, tskx, 
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   for( i = 0; i < numvar; i++ )
   {
      MOSEK_CALL( MSK_getbound(lpi->task,MSK_ACC_VAR,i,&tbkx[i],&tblx[i],&tbux[i]) );
   }
   
   for( i = 0; i < numcon; i++ )
   {
      MOSEK_CALL( MSK_getbound(lpi->task,MSK_ACC_CON,i,&tbkc[i],&tblc[i],&tbuc[i]) );
   }
   
   for( i = 0; i < numcon; ++i )
   {
      if( gotbasicsol )
      {
         if( ( tskc[i] == MSK_SK_FIX && tbkc[i] != MSK_BK_FX ) ||
            ( tskc[i] == MSK_SK_LOW && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_RA ) ) ||
            ( tskc[i] == MSK_SK_UPR && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_RA ) ) )
         {
            SCIPerrorMessage("STATUS KEY ERROR i %d bkc %d skc %d %s\n", i, tbkc[i], tskc[i], functionname);
         }
      }
      
      if( tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA ) 
      {
         if( isnan(tblc[i]) )
         {
            SCIPdebugMessage("nan in blc : %s\n", functionname);
         }
      }  
      
      if( tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA ) 
      {
         if( isnan(tbuc[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
         }
      }  
   }
   
   for( i = 0; i < numvar; ++i )
   {
      if( tbkx[i] == MSK_BK_LO || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA ) 
      {
         if( isnan(tblx[i]) )
         {
            SCIPdebugMessage("nan in blx : %s\n",functionname);
         }
      }  

      if( tbkx[i] == MSK_BK_UP || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA ) 
      {
         if( isnan(tbux[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
            getchar();
         }
      }  
   }

   BMSfreeMemoryArray(&tbkc);
   BMSfreeMemoryArray(&tskc);
   BMSfreeMemoryArray(&tblc);
   BMSfreeMemoryArray(&tbuc);
   BMSfreeMemoryArray(&tbkx);
   BMSfreeMemoryArray(&tskx);
   BMSfreeMemoryArray(&tblx);
   BMSfreeMemoryArray(&tbux);

   return SCIP_OKAY;
}
#endif


/*
 * Local functions
 */

static 
void generateMskBounds(
   int                   n, 
   const double*         lb, 
   const double*         ub, 
   MSKboundkeye*         bk, 
   double*               msklb, 
   double*               mskub
   )
{
   int i;
 
   assert(lb != NULL);
   assert(ub != NULL);
   assert(bk != NULL);
   assert(msklb != NULL);
   assert(mskub != NULL);

   for( i = 0; i < n; i++ )
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
         else if (lb[i] == ub[i])  /**@todo is this good idea to compare the bound without any epsilontic? */
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

/** get end pointers of arrays */
static 
SCIP_RETCODE getEndptrs(
   int                   n,                  /**< array size */
   const int*            beg,                /**< array of beginning indices */
   int                   nnonz,              /**< number of nonzeros */
   int**                 aptre               /**< pointer to store the result */
   )
{
   int i;

   assert(beg != NULL || nnonz == 0);
   
   SCIP_ALLOC( BMSallocMemoryArray( aptre, n) );

   /*    if (aptre == NULL) 
         return NULL;
   */

   if (nnonz > 0)
   {
      assert(beg != NULL);
      for(i = 0; i < n-1; i++)
      {
         (*aptre)[i] = beg[i+1];
         assert((*aptre)[i] >= beg[i]);
      }

      (*aptre)[n-1] = nnonz;
      assert((*aptre)[n-1] >= beg[n-1]);
   }
   else
   {
      for( i = 0; i < n; i++ )
         (*aptre)[i] = 0;
   }

   return SCIP_OKAY;
}

/** compute indices from range */
static 
SCIP_RETCODE getIndicesRange(
   int                   first,              /**< first index */ 
   int                   last,               /**< last index */ 
   int**                 sub                 /**< pointer to store the indices ranges */
   )
{
   int i;
   
   assert(first <= last);
   
   SCIP_ALLOC( BMSallocMemoryArray( sub, (last-first+1)) );
   
   for( i = first; i <= last; i++ )
   {
      (*sub)[i-first] = i;
   }
   
   return SCIP_OKAY;
}

/** compute indices from dense array */
static 
SCIP_RETCODE getIndicesFromDense(
   int*                  dstat,              /**< array */
   int                   n,                  /**< size of array */
   int*                  count,              /**< array of counts (sizes) */
   int**                 sub                 /**< pointer to store array of indices */
   )
{
   int i;
   int j;

   assert(dstat != NULL);

   *count = 0;
   for( i = 0; i < n; i++ ) 
   {
      if (dstat[i] == 1)
      {
         (*count)++;
      }
   }

   if( (*count) > 0 )
   {
      SCIP_ALLOC( BMSallocMemoryArray( sub, (*count)) );
   }
   else 
      return SCIP_OKAY;

   j = 0;
   for( i = 0; i < n; i++ )
   { 
      if (dstat[i] == 1)
      {
         (*sub)[j++] = i;
      }
   }
  
   return SCIP_OKAY;
}

static 
void scale_vec(
   int                   len, 
   double*               vec, 
   double                s
   )
{
   int i;
   for( i = 0; i < len; i++ )
   {
      vec[i] *= s;
   }
}

static 
void scale_bound(
   MSKboundkeye*         bk, 
   double*               bl, 
   double*               bu, 
   double                s
   )
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
      assert(FALSE);
      break;
   }  /*lint !e788*/

   if (s < 0)
   {
      double tmp;
      tmp = *bl;
      *bl = *bu;
      *bu = tmp;
   }
}

static 
SCIP_RETCODE ensureStateMem(
   SCIP_LPI*             lpi, 
   int                   ncols, 
   int                   nrows
   )
{
   if (lpi->skxsize < ncols)
   {
      int newsize;
      newsize = MAX(2*lpi->skxsize, ncols);
      
      SCIP_ALLOC( BMSreallocMemoryArray( &(lpi->skx), newsize) );
      lpi->skxsize = newsize;
   }

   if (lpi->skcsize < nrows)
   {
      int newsize;
      newsize = MAX(2*lpi->skcsize, nrows);

      SCIP_ALLOC( BMSreallocMemoryArray( &(lpi->skc), newsize) );
      lpi->skcsize = newsize;
   }

   return SCIP_OKAY;
}

static 
SCIP_RETCODE getbase(
   SCIP_LPI*             lpi, 
   int                   ncols, 
   int                   nrows
   )
{
   SCIPdebugMessage("Calling getbase (%d)\n",lpi->lpid);

   SCIP_CALL( ensureStateMem(lpi,ncols,nrows) );
   MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, lpi->skc, lpi->skx, 
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

static 
SCIP_RETCODE setbase(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   SCIPdebugMessage("Calling setbase (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_putsolution(lpi->task, MSK_SOL_BAS, lpi->skc, lpi->skx, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   
   return SCIP_OKAY;
}



/*
 * Miscellaneous Methods
 */

static char mskname[100];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   sprintf(mskname, "MOSEK %.2f", (SCIP_Real)MSK_VERSION_MAJOR);
   return mskname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by MOSEK Optimization Software (www.mosek.com)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->task;
}


/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != NULL);
   assert(numlp >= 0);

   SCIPdebugMessage("Calling SCIPlpiCreate\n");

   if (!MosekEnv)
   {
#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_makeenv(&MosekEnv, NULL, NULL, NULL, NULL) );
#else
      MOSEK_CALL( MSK_makeenv(&MosekEnv, NULL) );
#endif
      MOSEK_CALL( MSK_linkfunctoenvstream(MosekEnv, MSK_STREAM_LOG, NULL, printstr) );
      MOSEK_CALL( MSK_initenv(MosekEnv) );
   }

   numlp++;

   SCIP_ALLOC( BMSallocMemory(lpi) );

   MOSEK_CALL( MSK_makeemptytask(MosekEnv, &((*lpi)->task)) );

   MOSEK_CALL( MSK_linkfunctotaskstream((*lpi)->task, MSK_STREAM_LOG, NULL, printstr) );

   MOSEK_CALL( MSK_putobjsense((*lpi)->task, SENSE2MOSEK(objsen)) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_MAX_NUM_SETBACKS, SETBACK_LIMIT) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE_SIMPLEX) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_DEGEN, DEGEN_LEVEL) );
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_SWITCH_OPTIMIZER, MSK_ON) );
   /* We only have status keys (recalculate dual solution without dual superbasics) */
   MOSEK_CALL( MSK_putintparam((*lpi)->task, MSK_IPAR_SIM_HOTSTART, MSK_SIM_HOTSTART_STATUS_KEYS) );
   MOSEK_CALL( MSK_puttaskname((*lpi)->task, (char*) name) );

   (*lpi)->termcode = MSK_RES_OK;
   (*lpi)->itercount = 0;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->lpid = nextlpid++;
   (*lpi)->skxsize = 0;
   (*lpi)->skcsize = 0;
   (*lpi)->skx = NULL;
   (*lpi)->skc = NULL;
   (*lpi)->messagehdlr = messagehdlr;

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert(numlp > 0);

   SCIPdebugMessage("Calling SCIPlpiFree (%d)\n",(*lpi)->lpid);

   MOSEK_CALL( MSK_deletetask(&(*lpi)->task) );

   BMSfreeMemoryArrayNull(&(*lpi)->skx);
   BMSfreeMemoryArrayNull(&(*lpi)->skc);
   BMSfreeMemory(lpi);

   numlp--;
   if (numlp == 0)
   {
      MOSEK_CALL( MSK_deleteenv(&MosekEnv) );
      MosekEnv = NULL;
   }

   return SCIP_OKAY;
}

/*
 * Modification Methods
 */


/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
{  /*lint --e{715}*/
   int* aptre;
   MSKboundkeye* bkc;
   MSKboundkeye* bkx;
   double* blc;
   double* buc;
   double* blx;
   double* bux;

   SCIPdebugMessage("Calling SCIPlpiLoadColLP (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   /* initialize all array with NULL */
   aptre = NULL;
   bkc = NULL;
   bkx = NULL;
   blc = NULL;
   buc = NULL;
   blx = NULL;
   bux = NULL;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiLoadColLP") );
#endif

   if (nrows > 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &bkc, nrows) );
      SCIP_ALLOC( BMSallocMemoryArray( &blc, nrows) );
      SCIP_ALLOC( BMSallocMemoryArray( &buc, nrows) );

      generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);
   }

   if (ncols > 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &bkx, ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &blx, ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &bux, ncols) );

      generateMskBounds(ncols, lb, ub, bkx, blx, bux);

      SCIP_CALL( getEndptrs(ncols, beg, nnonz, &aptre) );
   }

   MOSEK_CALL( MSK_inputdata(lpi->task, nrows, ncols, nrows, ncols, obj, 0.0, beg, aptre, ind, val, 
         bkc, blc, buc, bkx, blx, bux) );

   MOSEK_CALL( MSK_putobjsense(lpi->task, SENSE2MOSEK(objsen)) );


   if( ncols > 0 )
   {
      BMSfreeMemoryArray(&aptre);
      BMSfreeMemoryArray(&bux);
      BMSfreeMemoryArray(&blx);
      BMSfreeMemoryArray(&bkx);
   }

   if( nrows > 0 )
   {
      BMSfreeMemoryArray(&buc);
      BMSfreeMemoryArray(&blc);
      BMSfreeMemoryArray(&bkc);
   }
   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiLoadColLP") );
#endif

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
{  /*lint --e{715}*/
#if MSK_VERSION_MAJOR < 7
   const int* aptrb;
#endif
   int* aptre;
   MSKboundkeye* bkx;
   double* blx;
   double* bux;
   int oldcols;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiAddCols (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddCols") );
#endif

   if (ncols == 0) 
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&blx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&bux, ncols) );
   generateMskBounds(ncols, lb, ub, bkx, blx, bux);

   MOSEK_CALL( MSK_getnumvar(lpi->task, &oldcols) );
   
#if MSK_VERSION_MAJOR < 7
   SCIP_CALL( getEndptrs(ncols, beg, nnonz, &aptre) );
   
   if (nnonz == 0)
      aptrb = aptre;
   else
      aptrb = beg;

   MOSEK_CALL( MSK_appendvars(lpi->task, ncols, obj, aptrb, aptre, ind, val, bkx, blx, bux) );

   BMSfreeMemoryArray(&aptre);

#else
   MOSEK_CALL( MSK_appendvars(lpi->task, ncols) );
   MOSEK_CALL( MSK_putcslice(lpi->task, oldcols, oldcols+ncols, obj) );
   MOSEK_CALL( MSK_putvarboundslice(lpi->task, oldcols, oldcols+ncols, bkx, blx, bux) );

   if( nnonz > 0 )
   {
      SCIP_CALL( getEndptrs(ncols, beg, nnonz, &aptre) );
      MOSEK_CALL( MSK_putacolslice(lpi->task, oldcols, oldcols+ncols, beg, aptre, ind, val) );
      BMSfreeMemoryArray(&aptre);
   }
#endif

   BMSfreeMemoryArray(&bux);
   BMSfreeMemoryArray(&blx);
   BMSfreeMemoryArray(&bkx);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddCols") );
#endif
   
   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int* sub;
   
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelCols (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelCols") );
#endif

   SCIP_CALL( getIndicesRange(firstcol, lastcol, &sub) );

   /*printf("Deleting vars %d to %d\n",firstcol,lastcol);*/
#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_remove(lpi->task,MSK_ACC_VAR, lastcol-firstcol+1, sub) );
#else
   MOSEK_CALL( MSK_removevars(lpi->task, lastcol-firstcol+1, sub) );
#endif

   BMSfreeMemoryArray(&sub);
   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelCols") );
#endif
   
   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int* sub;
   int count;
   int ncols;
   int col;
   int i;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelColset (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelColset") );
#endif

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, ncols, &count, &sub) );

   col = 0;
   for( i = 0; i < ncols; i++)
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
      SCIPdebugMessage("Deleting %d vars %d,...\n", count, sub[0]);
#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_remove(lpi->task, MSK_ACC_VAR, count, sub) );
#else
      MOSEK_CALL( MSK_removevars(lpi->task, count, sub) );
#endif
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelColset") );
#endif

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{  /*lint --e{715}*/
#if MSK_VERSION_MAJOR < 7
   const int* aptrb;
#endif
   int* aptre;
   MSKboundkeye* bkc;
   double* blc;
   double* buc;
   int oldrows;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiAddRows (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddRows") );
#endif

   if (nrows == 0) 
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&blc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&buc, nrows) );

   generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &oldrows) );

#if MSK_VERSION_MAJOR < 7
   SCIP_CALL( getEndptrs(nrows, beg, nnonz, &aptre) );

   if (nnonz == 0)
      aptrb = aptre;
   else
      aptrb = beg;

   MOSEK_CALL( MSK_appendcons(lpi->task, nrows, aptrb, aptre, ind, val, bkc, blc, buc) );

   BMSfreeMemoryArray(&aptre);

#else
   MOSEK_CALL( MSK_appendcons(lpi->task, nrows) );
   MOSEK_CALL( MSK_putconboundslice(lpi->task, oldrows, oldrows+nrows, bkc, blc, buc) );

   if( nnonz > 0 )
   {
      SCIP_CALL( getEndptrs(nrows, beg, nnonz, &aptre) );
      MOSEK_CALL( MSK_putarowslice(lpi->task, oldrows, oldrows+nrows, beg, aptre, ind, val) );
      BMSfreeMemoryArray(&aptre);
   }
#endif

   BMSfreeMemoryArray(&buc);
   BMSfreeMemoryArray(&blc);
   BMSfreeMemoryArray(&bkc);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiAddRows") );
#endif

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int* sub;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelRows (%d)\n",lpi->lpid);
   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRows") );
#endif

   SCIP_CALL( getIndicesRange(firstrow, lastrow, &sub) );

   SCIPdebugMessage("Deleting cons %d to %d\n",firstrow,lastrow);

#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_remove(lpi->task, MSK_ACC_CON, lastrow-firstrow+1, sub) );
#else
   MOSEK_CALL( MSK_removecons(lpi->task, lastrow-firstrow+1, sub) );
#endif

   BMSfreeMemoryArray(&sub);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRows") );
#endif

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int* sub;
   int count;
   int nrows;
   int row;
   int i;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiDelRowset (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRowset") );
#endif

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, nrows, &count, &sub) );

   row = 0;
   for( i = 0; i < nrows; i++ )
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
      SCIPdebugMessage("Deleting %d cons %d,...\n",count,sub[0]);
#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_remove(lpi->task, MSK_ACC_CON, count, sub) );
#else
      MOSEK_CALL( MSK_removecons(lpi->task, count, sub) );
#endif
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiDelRowset end") );
#endif

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int nrows;
   int ncols;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiClear (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   
   SCIP_CALL( SCIPlpiDelRows(lpi, 0, nrows) );
   SCIP_CALL( SCIPlpiDelCols(lpi, 0, ncols) );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   MSKboundkeye* bkx;
   double* blx;
   double* bux;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgBounds (%d)\n",lpi->lpid);
   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgBounds") );
#endif
   
   if (ncols == 0) 
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&blx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&bux, ncols) );

   generateMskBounds(ncols, lb, ub, bkx, blx, bux);
   MOSEK_CALL( MSK_putboundlist(lpi->task, MSK_ACC_VAR, ncols, ind, bkx, blx, bux) );

   BMSfreeMemoryArray(&bux);
   BMSfreeMemoryArray(&blx);
   BMSfreeMemoryArray(&bkx);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgBounds") );
#endif

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   MSKboundkeye* bkc;
   double* blc;
   double* buc;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgSides (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgSides") );
#endif

   if (nrows == 0) 
      return SCIP_OKAY;
   
   SCIP_ALLOC( BMSallocMemoryArray(&bkc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&blc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&buc, nrows) );

   generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);
   MOSEK_CALL( MSK_putboundlist(lpi->task, MSK_ACC_CON, nrows, ind, bkc, blc, buc) );

   BMSfreeMemoryArray(&buc);
   BMSfreeMemoryArray(&blc);
   BMSfreeMemoryArray(&bkc);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgSides") );
#endif

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgCoef (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgCoef") );
#endif

   MOSEK_CALL( MSK_putaij(lpi->task, row, col, newval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgCoef") );
#endif

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgObjsen (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_putobjsense(lpi->task, SENSE2MOSEK(objsen)) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiChgObj (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiChgObj") );
#endif

   MOSEK_CALL( MSK_putclist(lpi->task, ncols, ind, obj) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi,"SCIPlpiChgObj") );
#endif

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   int nnonz;
   int* sub;
   double* val;
   MSKboundkeye bkc;
   double blc;
   double buc;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiScaleRow (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleRow") );
#endif

   assert(scaleval != 0);

#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_getavecnumnz(lpi->task, MSK_ACC_CON, row, &nnonz) );
#else
   MOSEK_CALL( MSK_getarownumnz(lpi->task, row, &nnonz) );
#endif

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_getavec(lpi->task, MSK_ACC_CON, row, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putavec(lpi->task, MSK_ACC_CON, row, nnonz, sub, val) );
#else
      MOSEK_CALL( MSK_getarow(lpi->task, row, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putarow(lpi->task, row, nnonz, sub, val) );
#endif

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_CON, row, &bkc, &blc, &buc) );
   scale_bound(&bkc, &blc, &buc, scaleval);
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_CON, row, bkc, blc, buc) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleRow") );
#endif
   
   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   int nnonz;
   int *sub    = NULL;
   double *val = NULL;
   MSKboundkeye bkx;
   double blx, bux, c;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiScaleCol (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleCol") );
#endif

   assert(scaleval != 0);
#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_getavecnumnz(lpi->task, MSK_ACC_VAR, col, &nnonz) );
#else
   MOSEK_CALL( MSK_getacolnumnz(lpi->task, col, &nnonz) );
#endif

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_getavec(lpi->task, MSK_ACC_VAR, col, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putavec(lpi->task, MSK_ACC_VAR, col, nnonz, sub, val) );
#else
      MOSEK_CALL( MSK_getacol(lpi->task, col, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putacol(lpi->task, col, nnonz, sub, val) );
#endif

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );
   scale_bound(&bkx, &blx, &bux, 1.0/scaleval);
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );

   MOSEK_CALL( MSK_getcslice(lpi->task, col, col+1, &c) );
   MOSEK_CALL( MSK_putcj(lpi->task, col, c*scaleval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiScaleCol") );
#endif

   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */


/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNRows (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNCols (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getnumvar(lpi->task, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetNNonz (%d)\n",lpi->lpid);
   
   MOSEK_CALL( MSK_getnumanz(lpi->task, nnonz) );

   return SCIP_OKAY;
}

static 
SCIP_RETCODE getASlice(
   SCIP_LPI*             lpi,
   MSKaccmodee           iscon,
   int                   first,
   int                   last,
   int*                  nnonz,
   int*                  beg,
   int*                  ind,
   double*               val
   )
{
   int* aptre;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   assert(first <= last);

   SCIPdebugMessage("Calling SCIPlpiGetNNonz (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "getASlice") );
#endif

   if( nnonz != 0 )
   {
      int surplus;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      SCIP_ALLOC( BMSallocMemoryArray(&aptre, last - first + 1) );

      MOSEK_CALL( MSK_getaslicenumnz(lpi->task, iscon, first, last+1,nnonz) );
      surplus = *nnonz;
      MOSEK_CALL( MSK_getaslice(lpi->task, iscon, first, last+1, *nnonz, &surplus, beg, aptre, ind, val) );

      assert(surplus == 0);

      BMSfreeMemoryArray(&aptre);
   }
   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "getASlice") );
#endif

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetCols (%d)\n",lpi->lpid);

   SCIP_CALL( SCIPlpiGetBounds(lpi, firstcol, lastcol, lb, ub) );
   SCIP_CALL( getASlice(lpi, MSK_ACC_VAR, firstcol, lastcol, nnonz, beg, ind, val) );

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetRows (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetRows") );
#endif


   SCIP_CALL( SCIPlpiGetSides(lpi, firstrow, lastrow, lhs, rhs) );
   SCIP_CALL( getASlice(lpi, MSK_ACC_CON, firstrow, lastrow, nnonz, beg, ind, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetRows") );
#endif

   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet.\n");
   return SCIP_LPERROR;
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");
   return SCIP_LPERROR;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPerrorMessage("SCIPlpiGetObjsen() has not been implemented yet.\n");
   return SCIP_LPERROR;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   SCIPdebugMessage("Calling SCIPlpiGetObj (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getcslice(lpi->task, firstcol, lastcol+1, vals) );

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBounds (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetBounds") );
#endif
   
   MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_VAR, firstcol, lastcol+1, NULL, lbs, ubs) );
   
   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetSides (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetSides") );
#endif

   MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_CON, firstrow, lastrow+1, NULL, lhss, rhss) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetSides") );
#endif

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetCoef (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetCoef") );
#endif

   MOSEK_CALL( MSK_getaij(lpi->task, row, col, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiGetCoef") );
#endif

   return SCIP_OKAY;
}

/*
 * Solving Methods
 */


/** gets the internal solution status of the solver */
static
SCIP_RETCODE getSolutionStatus(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   MSKprostae*           prosta,             /**< pointer to store the problem status */
   MSKsolstae*           solsta              /**< pointer to store the solution status */
   )
{
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   MOSEK_CALL( MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, prosta, solsta) );
   
   return SCIP_OKAY;
}


static
MSKrescodee filterTRMrescode(
   SCIP_MESSAGEHDLR*     messagehdlr,
   MSKrescodee*          termcode,
   MSKrescodee           res
   )
{
   if (   res == MSK_RES_TRM_MAX_ITERATIONS || res == MSK_RES_TRM_MAX_TIME
      || res == MSK_RES_TRM_OBJECTIVE_RANGE || res == MSK_RES_TRM_STALL
#if ASSERT_ON_NUMERICAL_TROUBLES > 0
      || res == MSK_RES_TRM_MAX_NUM_SETBACKS
      || res == MSK_RES_TRM_NUMERICAL_PROBLEM
#endif
          )
   {
      *termcode = res;
      if (res == MSK_RES_TRM_MAX_NUM_SETBACKS || res == MSK_RES_TRM_NUMERICAL_PROBLEM)
      {
         SCIPmessagePrintWarning(messagehdlr, "Return code %d in [%d]\n", res, optimizecount);

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

static 
SCIP_RETCODE SolveWSimplex(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int itercount_primal;
   int itercount_dual;
   int gotbasicsol;
   int presolve;
   int maxiter;
   MSKprostae prosta;
   MSKsolstae solsta;
   double pobj,dobj;

   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, &presolve) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &maxiter) );

#if DEBUG_EASY_REPRODUCE
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_AUTO_SORT_A_BEFORE_OPT, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_OFF) );
#else
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_AUTO_UPDATE_SOL_INFO, MSK_OFF) );

#if FORCE_MOSEK_LOG
   
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_ON) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM_FREQ, 1) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#endif
 
   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if( gotbasicsol )
   { 
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF) );
   }
   else
   { 
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_ON) );
   }

#if ALWAYS_SOLVE_PRIMAL > 0
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SOLVE_FORM, MSK_SOLVE_PRIMAL) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SolveWSimplex") );
#endif

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   if( gotbasicsol && maxiter < 20000 )
   {
      /* Since max iter often is set, we switch off restricted pricing */
      MOSEK_CALL( MSK_putintparam(lpi->task,  MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, 0) );  
   }

   if( FORCE_NO_MAXITER )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, 2000000000) );
   }

   
#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "Begin optimize with simplex") );
#endif
   
#if FORCE_MOSEK_SUMMARY > 1
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(lpi->task,MSK_STREAM_LOG) );
   }
#endif
   
#if !FORCE_SILENCE
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, 100) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM_FREQ, 100) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM, 100) );
#endif

   MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );

   if( lpi->termcode == MSK_RES_TRM_MAX_NUM_SETBACKS )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_AGGRESSIVE) );

      MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );
   }

#if FORCE_MOSEK_SUMMARY
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(lpi->task,MSK_STREAM_LOG) );
   }
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "End optimize with simplex") );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, presolve) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_SIM_PRIMAL_ITER, &itercount_primal) );
   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_SIM_DUAL_ITER, &itercount_dual) );

   lpi->itercount = itercount_primal + itercount_dual;

   MOSEK_CALL( MSK_getprimalobj(lpi->task, MSK_SOL_BAS, &pobj) );
   MOSEK_CALL( MSK_getdualobj(lpi->task, MSK_SOL_BAS, &dobj) );
   MOSEK_CALL( MSK_getsolutionstatus(lpi->task, MSK_SOL_BAS, &prosta, &solsta) );

#if DEBUG_PRINT_STAT
   SCIPdebugMessage("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, objval = %g : %g, iter = %d+%d\n",
      maxiter, lpi->termcode, prosta, solsta, pobj, dobj, itercount_primal, itercount_dual);
#endif

   SCIPdebugMessage("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, "
      "objval = %g : %g, iter = %d+%d\n",
      maxiter,lpi->termcode,prosta,solsta,
      pobj,dobj,itercount_primal,itercount_dual);
   
   /*  SCIPdebugMessage("Iter dual %d primal %d\n",itercount_dual,itercount_primal); */
   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
   case MSK_SOL_STA_UNKNOWN:
      break;
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] returned solsta = %d\n", optimizecount, solsta);
      
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
      SCIPerrorMessage("Simplex[%d] returned solsta = %d\n", optimizecount, solsta);
#endif
      
#if ASSERT_ON_WARNING
      assert(0);
#endif
      
      return SCIP_LPERROR;
   }  /*lint !e788*/

   switch (prosta)
   {
   case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_FEAS:
   case MSK_PRO_STA_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
   case MSK_PRO_STA_PRIM_INFEAS:
   case MSK_PRO_STA_DUAL_INFEAS:
   case MSK_PRO_STA_UNKNOWN:
      break;
   case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_NEAR_PRIM_FEAS:
   case MSK_PRO_STA_NEAR_DUAL_FEAS:
   case MSK_PRO_STA_ILL_POSED:
   case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] returned prosta = %d\n", optimizecount, prosta);

      if (lpi->termcode == MSK_RES_OK)
         lpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
   default:
#if SHOW_ERRORS && !FORCE_SILENCE
      SCIPerrorMessage("Simplex[%d] returned prosta = %d\n", optimizecount, prosta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_LPERROR;
   }  /*lint !e788*/

   if( solsta == MSK_SOL_STA_OPTIMAL && fabs(dobj)+fabs(dobj) > 1.0e-6 && fabs(pobj-dobj)>0.0001*(fabs(pobj)+fabs(dobj)))
   {
      SCIPerrorMessage("Simplex[%d] returned optimal solution with different objvals %g != %g reldiff %.2g%%\n",
         optimizecount, pobj, dobj, 100*fabs(pobj-dobj)/ MAX(fabs(pobj),fabs(dobj)));
   }
   
   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
   {
      if(solsta != MSK_SOL_STA_DUAL_FEAS && solsta != MSK_SOL_STA_OPTIMAL && solsta != MSK_SOL_STA_PRIM_AND_DUAL_FEAS)
      {
         SCIPerrorMessage("[%d] Terminated on objective range without dual feasible solsta.\n", optimizecount);

         SCIP_CALL( SCIPlpiSolveBarrier(lpi, TRUE) );
      }
      else
      {
         scipmskobjsen objsen;
         double bound;

         MOSEK_CALL( MSK_getobjsense(lpi->task, &objsen) );

         if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
         {
            MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, &bound) );

            if (1.0e-6*(fabs(bound)+fabs(dobj)) < bound-dobj)
            {
               SCIPerrorMessage("[%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount, dobj, bound);
               
               SCIP_CALL( SCIPlpiSolveBarrier(lpi, TRUE) );
            }
         }
         else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
         {
            MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, &bound) );

            if (1.0e-6*(fabs(bound)+fabs(dobj)) < dobj-bound)
            {
               SCIPerrorMessage("[%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount, dobj, bound);
               
               SCIP_CALL( SCIPlpiSolveBarrier(lpi, TRUE) );
            }
         }
      }
   }
   
   if (maxiter >= 2000000000)
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      {
         SCIPmessagePrintWarning(lpi->messagehdlr, "Simplex[%d] failed to terminate in 10000 iterations, switching to interior point\n",
            optimizecount);

         SCIP_CALL( SCIPlpiSolveBarrier(lpi,TRUE) );
      }
   }

#if DEBUG_DO_INTPNT_FEAS_CHECK
   if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
   {
      SCIPdebugMessage("Checking infeasibility[%d]... ",optimizecount);

      SCIP_CALL( SCIPlpiSolveBarrier(lpi,true) );

      MOSEK_CALL(MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, &prosta, &solsta));

      if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
      {
         SCIPdebugPrintf("ok\n");
      }
      else
      {
         SCIPdebugPrintf("wrong [%d] prosta = %d, solsta = %d\n",optimizecount,prosta,solsta);
      }
   }
#endif
   

#if DEBUG_PRINT_STAT > 0
   SCIPdebugMessage("Max iter stat    : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      optimizecount, numstrongbranchmaxiterup, numstrongbranchmaxiterdo, numprimalmaxiter, numdualmaxiter);
   SCIPdebugMessage("Objcut iter stat : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      optimizecount, numstrongbranchobjup, numstrongbranchobjdo, numprimalobj, numdualobj);
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SolveWSimplex") );
#endif

   return SCIP_OKAY;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   optimizecount++;

   SCIPdebugMessage("Calling SCIPlpiSolvePrimal[%d] (%d) ",optimizecount,lpi->lpid);

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolvePrimal") );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_PRIMAL_SIMPLEX) );
   
#if WRITE_PRIMAL > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"primal_%d.lp",optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n",fname);
      /*MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
      MSK_writedata(lpi->task,fname);
   }
#endif

   SCIP_CALL( SolveWSimplex(lpi) );

   if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, NULL, &solsta) );


      if( solsta != MSK_SOL_STA_PRIM_FEAS )
      {
         SCIP_CALL( SolveWSimplex(lpi) );
      }
   }

   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numprimalobj;

   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numprimalmaxiter;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolvePrimal") );
#endif

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   optimizecount++;

   SCIPdebugMessage("Calling SCIPlpiSolveDual[%d] (%d)\n",optimizecount,lpi->lpid);


   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_INTEGER, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_DUAL_SIMPLEX) );

#if WRITE_DUAL > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"dual_%d.lp",optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n",fname);
      MSK_writedata(lpi->task,fname);
   }
#endif

#if !FORCE_SILENCE
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG_SIM_FREQ, 1) );
#endif

   SCIP_CALL( SolveWSimplex(lpi) );

   if ( lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, NULL, &solsta) );

      if( solsta != MSK_SOL_STA_DUAL_FEAS )
      {
         SCIP_CALL( SolveWSimplex(lpi) );
      }
   }

   if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numdualobj;

   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   optimizecount++;

#if FORCE_MOSEK_LOG
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_ON) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#else
   {
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#endif


   SCIPdebugMessage("Calling SCIPlpiSolveBarrier[%d] (%d) ",optimizecount,lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolveBarrier") );
#endif

   MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_INTPNT_BASIS, crossover ? MSK_BI_ALWAYS : MSK_BI_NEVER) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT) );


#if WRITE_INTPNT > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"intpnt_%d.lp",optimizecount);
      SCIPdebugMessage("\nWriting lp %s\n",fname);
      /*MOSEK_CALL( MSK_putintparam(lpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
      MSK_writedata(lpi->task,fname);
   }
#endif

   MOSEK_CALL( filterTRMrescode(lpi->messagehdlr, &lpi->termcode, MSK_optimize(lpi->task)) );

   if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;

   MOSEK_CALL( MSK_getintinf(lpi->task, MSK_IINF_INTPNT_ITER, &lpi->itercount) );

#ifdef SCIP_DEBUG
   {
      MSKprostae prosta;
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, &prosta, &solsta) );
      SCIPdebugMessage("termcode = %d, prosta = %d, solsta = %d, iter = %d\n",
         lpi->termcode, prosta, solsta, lpi->itercount);
   }
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiSolveBarrier") );
#endif

   return SCIP_OKAY;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on all candidates */
static
SCIP_RETCODE SCIPlpiStrongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
   MSKobjsensee objsen;
   int olditerlim;
   int oldselection;
   int oldhotstart;

   double bound;
   int ncols;
   int nrows;
   MSKboundkeye bkx;
   double blx;
   double bux;
   double newub;
   double newlb;
   
   SCIPdebugMessage("Calling SCIPlpiStrongbranch (%d)\n",lpi->lpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiStrongbranch") );
#endif

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   if (lpi->termcode != MSK_RES_OK)
   {
      SCIPmessagePrintWarning(lpi->messagehdlr, "SB Warning: Previous termcode is %d\n",lpi->termcode);
   }
   
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( getbase(lpi, ncols, nrows) );

   MOSEK_CALL( MSK_getobjsense(lpi->task, &objsen) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &olditerlim) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, &oldselection) );
   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, &oldhotstart) );

   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, itlim) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, STRONGBRANCH_PRICING) );

   if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
   {
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, &bound) );
   }
   else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
   {
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, &bound) );
   }

   MOSEK_CALL( MSK_getbound(lpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );

   *iter = 0;

   newub = EPSCEIL(psol-1.0, 1e-06);

   if (newub < blx - 0.5) /* infeasible */
   {
      *down = bound;
      *downvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_NEGINF(blx))
         newbk = MSK_BK_UP;
      else if (EPSEQ(blx,newub,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newub = blx;
      }
      else
         newbk = MSK_BK_RA;
      
      MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, newbk, blx, newub) );

      SCIP_CALL( SCIPlpiSolveDual(lpi) );
      
      *iter += lpi->itercount;

      if (SCIPlpiIsStable(lpi))
         *downvalid = TRUE;
      else
         *downvalid = FALSE;

      if (SCIPlpiExistsPrimalRay(lpi))
      {
         SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is dual infeasible\n",optimizecount);
         
         *down = -1e20;
         *downvalid = FALSE;
      }
      else if (SCIPlpiExistsDualRay(lpi))
      {
         *down = bound;
      }
      else
      {
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPlpiGetSolFeasibility(lpi, NULL, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n", optimizecount);

            *down = -1e20;
            *downvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(lpi->task, MSK_SOL_BAS, down) );
         }
      }
      
      if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjup;

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterup;
   }

   /* Reset basis solution before doing the up branch */
   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
   SCIP_CALL( setbase(lpi) );

   newlb = EPSFLOOR(psol+1.0, 1e-06);
   if (newlb > bux + 0.5) /* infeasible */
   {
      *up = bound;
      *upvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_POSINF(bux))
         newbk = MSK_BK_LO;
      else if (EPSEQ(bux,newlb,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newlb = bux;
      }
      else
         newbk = MSK_BK_RA;

      MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, newbk, newlb, bux) );
      SCIP_CALL( SCIPlpiSolveDual(lpi) );

      *iter += lpi->itercount;

      if (SCIPlpiIsStable(lpi))
         *upvalid = TRUE;
      else
         *upvalid = FALSE;

      if (SCIPlpiExistsPrimalRay(lpi))
      {
         *up = -1e20;
         *upvalid = FALSE;
      }
      else if (SCIPlpiExistsDualRay(lpi))
      {
         *up = bound;
      }
      else
      {
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPlpiGetSolFeasibility(lpi, NULL, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(lpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n",optimizecount);
            
            *up = -1e20;
            *upvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(lpi->task, MSK_SOL_BAS, up) );
         }
      }

      if (lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjdo;

      if (lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterdo;
   }

   MOSEK_CALL( MSK_putbound(lpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, olditerlim) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, oldselection) );
   MOSEK_CALL( MSK_putintparam(lpi->task,   MSK_IPAR_SIM_HOTSTART, oldhotstart) );

   SCIP_CALL( setbase(lpi) );

   lpi->termcode = MSK_RES_OK;
   lpi->itercount = 0;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(lpi, "SCIPlpiStrongbranch") );
#endif

   SCIPdebugMessage("End SCIPlpiStrongbranch (%d)\n", lpi->lpid);

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
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
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( SCIPlpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( SCIPlpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
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
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( SCIPlpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );
   
   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( SCIPlpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}


/*
 * Solution Information Methods
 */


/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiWasSolved (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus (lpi, &prosta, &solsta) );

   return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   MSKsolstae solsta;
   SCIP_Bool pfeas;
   SCIP_Bool dfeas;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetSolFeasibility (%d)\n",lpi->lpid);

   pfeas = FALSE;
   dfeas = FALSE;

   MOSEK_CALL( MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, NULL, &solsta) );

   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      pfeas = TRUE;
      dfeas = TRUE;
      break;
   case MSK_SOL_STA_PRIM_FEAS:
      pfeas = TRUE;
      break;
   case MSK_SOL_STA_DUAL_FEAS:
      dfeas = TRUE;
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
      return SCIP_LPERROR;
   }  /*lint !e788*/

   if( primalfeasible != NULL ) 
      *primalfeasible = pfeas;

   if( dualfeasible != NULL ) 
      *dualfeasible = dfeas;
   
   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiExistsPrimalRay (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus (lpi, &prosta, &solsta));

   return (   solsta == MSK_SOL_STA_DUAL_INFEAS_CER
      || prosta == MSK_PRO_STA_DUAL_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
 
   SCIPdebugMessage("Calling SCIPlpiHasPrimalRay (%d)\n",lpi->lpid);
   
   SCIP_CALL_ABORT( getSolutionStatus (lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_DUAL_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return SCIPlpiExistsDualRay(lpi);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsPrimalFeasible (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus (lpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_PRIM_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiExistsDualRay (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus(lpi, &prosta, &solsta) );

   return (   solsta == MSK_SOL_STA_PRIM_INFEAS_CER
      || prosta == MSK_PRO_STA_PRIM_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiHasDualRay (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus (lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_PRIM_INFEAS_CER);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   return FALSE;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return SCIPlpiExistsPrimalRay(lpi);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKprostae prosta;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsDualFeasible (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus(lpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_DUAL_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}


/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiIsOptimal (%d)\n",lpi->lpid);

   SCIP_CALL_ABORT( getSolutionStatus(lpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   return (   lpi->termcode == MSK_RES_OK
      || lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS
      || lpi->termcode == MSK_RES_TRM_MAX_TIME
      || lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   return lpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   return lpi->termcode == MSK_RES_TRM_MAX_ITERATIONS;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   return lpi->termcode == MSK_RES_TRM_MAX_TIME;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   MSKsolstae solsta;
   SCIP_RETCODE retcode;

   assert(MosekEnv != NULL);
   assert(lpi != NULL);
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetInternalStatus (%d)\n", lpi->lpid);

   retcode = getSolutionStatus(lpi, NULL, &solsta);
   if ( retcode != SCIP_OKAY )
      return 0;

   return solsta; /*lint !e641*/
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   SCIPdebugMessage("Calling SCIPlpiIgnoreInstability (%d)\n",lpi->lpid);

   *success = FALSE;
   
   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
 
   SCIPdebugMessage("Calling SCIPlpiGetObjval (%d)\n",lpi->lpid);
 
   MOSEK_CALL( MSK_getprimalobj(lpi->task, MSK_SOL_BAS, objval) );

   /* TODO: tjek lighed med dual objektiv i de fleste tilfaelde. */

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{  /*lint --e{715}*/
   double* sux;
   int ncols;
   int i;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
 
   SCIPdebugMessage("Calling SCIPlpiGetSol (%d)\n",lpi->lpid);

   sux = NULL;
   ncols = 0;

   if( redcost )
   {
      MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &sux, ncols) );
   }

   MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, NULL, NULL, NULL, activity,
         primsol, dualsol, NULL, NULL, redcost, sux, NULL) );

   if( redcost ) 
   {
      for( i = 0; i < ncols; i++ )
      {
         assert(sux != NULL);
         redcost[i] -= sux[i];
      }
   }
  
   BMSfreeMemoryArray(&sux);

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetPrimalRay (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, NULL, NULL, NULL, NULL, ray,
         NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetDualfarkas (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getsolution(lpi->task, MSK_SOL_BAS, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dualfarkas,
         NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("Calling SCIPlpiGetIterations (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   *iterations = lpi->itercount;

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert(lpi != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/** handle singular basis */
static 
SCIP_RETCODE handle_singular(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  basis,              /**< array of basis indices */
   MSKrescodee           res                 /**< result */
   )
{
   if (res == MSK_RES_ERR_BASIS_SINGULAR)
   {
      SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

      MOSEK_CALL( MSK_initbasissolve(lpi->task, basis) );
   }
   else
   {
      MOSEK_CALL( res );
   }

   return SCIP_OKAY;
}


/*
 * LP Basis Methods
 */

/** convert Mosek status to SCIP status */
static
SCIP_RETCODE convertstat_mosek2scip(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   MSKaccmodee           acc,                /**< ??? */
   MSKstakeye*           sk,                 /**< ??? */
   int                   n,                  /**< size */
   int*                  stat                /**< status array */
   )
{
   int i;

   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;

      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
         MOSEK_CALL( MSK_getsolutioni(lpi->task, acc, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_INF:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_LOW:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UPR:
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
      case MSK_SK_END:
         break;
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert Mosek to SCIP status - slack variables */
static
SCIP_RETCODE convertstat_mosek2scip_slack(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   MSKaccmodee           acc,                /**< ??? */
   MSKstakeye*           sk,                 /**< ??? */
   int                   n,                  /**< size */
   int*                  stat                /**< status array */
   )
{
   int i;

   /* slacks are stored as -1 in Mosek, i.e., bounds are reversed compared to SCIP  */

   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
         MOSEK_CALL( MSK_getsolutioni(lpi->task, acc, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
      case MSK_SK_INF:
      case MSK_SK_UPR: /* Reversed */
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_LOW: /* Reversed */
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
      case MSK_SK_END:
         break;
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert SCIP to Mosek status */
static
void convertstat_scip2mosek(
   int*                  stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   int i;
   for( i = 0; i < n; i++ )
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
      default:
         SCIPABORT();
      }
   }
}

/** convert SCIP to Mosek status - slack variables */
static 
void convertstat_scip2mosek_slack(
   int*                  stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   /* slacks are stored as -1 in Mosek, i.e., bounds are reversed compared to SCIP  */
   int i;

   for( i = 0; i < n; i++ )
   {
      switch (stat[i])
      {
      case SCIP_BASESTAT_LOWER:
         resstat[i] = MSK_SK_UPR;/* Reversed */
         break;
      case SCIP_BASESTAT_BASIC:
         resstat[i] = MSK_SK_BAS;
         break;
      case SCIP_BASESTAT_UPPER:
         resstat[i] = MSK_SK_LOW; /* Reversed */
         break;
      case SCIP_BASESTAT_ZERO:
         resstat[i] = MSK_SK_SUPBAS;
         break;
      default:
         SCIPABORT();
      }
   }
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPlpiGetBase (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( getbase(lpi, ncols, nrows) );

   if (cstat) 
   {
      SCIP_CALL( convertstat_mosek2scip(lpi, MSK_ACC_VAR, lpi->skx, ncols, cstat) );
   }
   
   if (rstat) 
   {
      SCIP_CALL( convertstat_mosek2scip_slack(lpi, MSK_ACC_CON, lpi->skc, nrows, rstat) );
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPlpiSetBase (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

   SCIP_CALL( ensureStateMem(lpi, ncols, nrows) );

   convertstat_scip2mosek(cstat, ncols, lpi->skx);
   convertstat_scip2mosek_slack(rstat, nrows, lpi->skc);

   SCIP_CALL( setbase(lpi) );

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int nrows;
   int i;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBasisInd (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );

#if 0
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_OFF) );
#endif

   SCIP_CALL( handle_singular(lpi,bind,MSK_initbasissolve(lpi->task, bind)) );

#if 0
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
#endif

   for (i = 0; i < nrows; i++ )
   {
      if (bind[i] < nrows) /* row bind[i] is basic */
         bind[i] = -1 - bind[i];
      else                 /* column bind[i]-nrows is basic */
         bind[i] = bind[i] - nrows;
   }

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   int* sub;
   int nrows;
   int numnz;
   int i;

   SCIPdebugMessage("Calling SCIPlpiGetBInvCol (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getnumcon(lpi->task,&nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );

   for (i=0; i<nrows; i++)
      coef[i] = 0;
   
   numnz = 1; 
   sub[0]= c; 
   coef[c] = 1; /* Unit vector e_col */

   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_OFF) );
   MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   SCIP_Real* val;
   int* sub;
   int nrows;
   int numnz;
   int i;

   SCIPdebugMessage("Calling SCIPlpiGetBInvACol (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getnumcon(lpi->task,&nrows) );
#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_getavecnumnz(lpi->task,MSK_ACC_VAR,c,&numnz) );
#else
   MOSEK_CALL( MSK_getacolnumnz(lpi->task,c,&numnz) );
#endif
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &val, numnz+1) );
   
   for (i=0; i<nrows; i++)
      coef[i] = 0;

#if MSK_VERSION_MAJOR < 7
   MOSEK_CALL( MSK_getavec(lpi->task, MSK_ACC_VAR, c, &numnz, sub, val) );
#else
   MOSEK_CALL( MSK_getacol(lpi->task, c, &numnz, sub, val) );
#endif

   for (i=0; i<numnz; i++)
      coef[sub[i]] = val[i];

   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_OFF) );
   MOSEK_CALL( MSK_solvewithbasis(lpi->task, 0, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);
   BMSfreeMemoryArray(&val);
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   int* sub;
   int nrows;
   int numnz;
   int i;
 
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBInvRow (%d)\n",lpi->lpid);

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );
   
   for (i=0; i<nrows; i++)
      coef[i] = 0;

   numnz = 1; 
   sub[0] = row; 
   coef[row] = 1; /* Unit vector e_row */

   MOSEK_CALL( MSK_putnaintparam(lpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_ON) );
   MOSEK_CALL( MSK_solvewithbasis(lpi->task, 1, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);

   SCIPdebugMessage("End SCIPlpiGetBInvRow (%d)\n",lpi->lpid);
   
   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            val                 /**< vector to return coefficients */
   )
{
   int nrows;
   int ncols;
   int numnz;
   int* csub;
   int didalloc;
   double* cval;
   double* binv;
   int i;
   int k;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiGetBInvARow (%d)\n",lpi->lpid);
   
   didalloc = 0;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   SCIP_ALLOC( BMSallocMemoryArray(&csub, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&cval, nrows) );
   
   for( i = 0; i < ncols; i++ )
      val[i] = 0;

   if( binvrow == NULL )
   {
      didalloc = 1;

      SCIP_ALLOC( BMSallocMemoryArray( &binv, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, row, binv) );
   }
   else
      binv = (SCIP_Real*)binvrow;

   /* binvrow*A */
   for( i = 0; i < ncols; i++)
   {
      val[i] = 0;

#if MSK_VERSION_MAJOR < 7
      MOSEK_CALL( MSK_getavec(lpi->task, MSK_ACC_VAR, i, &numnz, csub, cval) );
#else
      MOSEK_CALL( MSK_getacol(lpi->task, i, &numnz, csub, cval) );
#endif

      for( k = 0; k < numnz; ++k )
         val[i] += binv[csub[k]] * cval[k];
   }
   
   /* free memory arrays */
   BMSfreeMemoryArray(&cval);
   BMSfreeMemoryArray(&csub);
   
   if( didalloc > 0 )
   {
      BMSfreeMemoryArray(&binv);
   }

   return SCIP_OKAY;
}

/*
 * LP State Methods
 */

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->skx, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->skc, rowpacketNum(nrows)) );

   if( lpistate[0] )
   {
      lpistate[0]->solsta = MSK_SOL_STA_UNKNOWN; 
      lpistate[0]->num    = -1;
      lpistate[0]->ncols  = ncols;
      lpistate[0]->nrows  = nrows;
   }

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skx, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->skc, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}

#ifndef NDEBUG
/** check state */
static
SCIP_RETCODE checkState1(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   n,                  /**< number of rows or columns */
   MSKstakeye*           sk,                 /**< ??? */
   MSKaccmodee           accmode,            /**< ??? */
   char                  xc                  /**< ??? */
   )
{
   int i;

   /* printout for all except LOW, UPR, FIX and BAS with sl[xc]==su[xc] */
   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_UNK:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unk\n", optimizecount, xc, i);
         break;
      case MSK_SK_BAS:
         MOSEK_CALL( MSK_getsolutioni(lpi->task, accmode, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );
         if (fabs(sl-su) > DEBUG_CHECK_STATE_TOL)
            SCIPdebugMessage("STATE[%d]: %c[%d] = bas, sl%c = %g, su%c = %g\n", optimizecount, xc, i, xc, sl, xc, su);
         break;
      case MSK_SK_SUPBAS:
         SCIPdebugMessage("STATE[%d]: %c[%d] = supbas\n", optimizecount, xc, i);
         break;
      case MSK_SK_LOW:
      case MSK_SK_UPR:
      case MSK_SK_FIX:
         break;
      case MSK_SK_INF:
         SCIPdebugMessage("STATE[%d]: %c[%d] = inf\n", optimizecount, xc, i);
         break;
      default:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unknown status <%d>\n", optimizecount, xc, i, sk[i]);
         break;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** check state */
static
SCIP_RETCODE checkState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns */
   int                   nrows               /**< number of rows */
   )
{
   SCIP_CALL( checkState1(lpi, ncols, lpi->skx, MSK_ACC_VAR, 'x') );
   SCIP_CALL( checkState1(lpi, nrows, lpi->skc, MSK_ACC_CON, 'c') );

   return SCIP_OKAY;
 }
#endif

/** store row and column basis status in a packed LPi state object */
static
SCIP_RETCODE lpistatePack(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< pointer to LPi state data */
   )
{
   int *skxi = (int *) lpi->skx; /* Used as temp. buffer */
   int *skci = (int *) lpi->skc; /* Used as temp. buffer */

   assert(sizeof(int) == sizeof(MSKstakeye));

   SCIP_CALL( convertstat_mosek2scip(lpi, MSK_ACC_VAR, lpi->skx, lpistate->ncols, skxi) );
   SCIP_CALL( convertstat_mosek2scip_slack(lpi, MSK_ACC_CON, lpi->skc, lpistate->nrows, skci) );

   SCIPencodeDualBit(skxi, lpistate->skx, lpistate->ncols);
   SCIPencodeDualBit(skci, lpistate->skc, lpistate->nrows);

   return SCIP_OKAY;
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   SCIP_LPISTATE*        lpistate,           /**< pointer to LPi state data */
   MSKstakeye*           skx,                /**< ??? */
   MSKstakeye*           skc                 /**< ??? */
   )
{
   assert(sizeof(int) == sizeof(MSKstakeye));

   SCIPdecodeDualBit(lpistate->skx, (int*) skx, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->skc, (int*) skc, lpistate->nrows);
   
   convertstat_scip2mosek((int*) skx, lpistate->ncols, skx);
   convertstat_scip2mosek_slack((int*) skc, lpistate->nrows, skc);
}

/** stores LP state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   int gotbasicsol;
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPlpiGetState (%d)\n",lpi->lpid);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   assert(lpistate != NULL);

   *lpistate = NULL;

   MOSEK_CALL( MSK_solutiondef(lpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if ( gotbasicsol == 0 || SCIPlpiExistsDualRay(lpi) ) 
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   lpistate[0]->num = optimizecount;

   MOSEK_CALL(MSK_getsolutionstatus ( lpi->task, MSK_SOL_BAS, NULL, &lpistate[0]->solsta));

   SCIP_CALL( getbase(lpi, ncols, nrows) );

#ifndef NDEBUG
   SCIP_CALL( checkState(lpi, ncols, nrows) );
#endif

   SCIP_CALL( lpistatePack(lpi, lpistate[0]) );

   SCIPdebugMessage("Store into state from iter : %d\n",optimizecount);
   
   /*    if (r != SCIP_OKAY) 
    *    lpistateFree(lpistate, blkmem );
    */
   
   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   int nrows;
   int ncols;
   int i;

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   if (lpistate == NULL)
   {
      SCIPdebugMessage("Setting NULL state\n");
      return SCIP_OKAY;
   }

   if (lpistate->nrows == 0 || lpistate->ncols == 0)
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(lpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(lpi->task, &ncols) );
   assert(lpistate->nrows <= nrows);
   assert(lpistate->ncols <= ncols);

   SCIP_CALL( ensureStateMem(lpi, ncols, nrows) );
   SCIP_CALL( getbase(lpi, ncols, nrows) );

   lpistateUnpack(lpistate, lpi->skx, lpi->skc);

   /* extend the basis to the current LP beyond the previously existing columns */
   for (i = lpistate->ncols; i < ncols; ++i)
   {
      SCIP_Real lb;
      SCIP_Real ub;
      MOSEK_CALL( MSK_getboundslice(lpi->task, MSK_ACC_VAR, i, i, NULL, &lb, &ub) );
      if ( SCIPlpiIsInfinity(lpi, REALABS(lb)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         if ( SCIPlpiIsInfinity(lpi, REALABS(ub)) )
            lpi->skx[i] = MSK_SK_SUPBAS;  /* variable is free (super basic) */
         else
            lpi->skx[i] = MSK_SK_UPR;     /* use finite upper bound */
      }
      else
         lpi->skx[i] = MSK_SK_LOW;        /* use finite lower bound */
   }
   for (i = lpistate->nrows; i < nrows; ++i)
      lpi->skc[i] = MSK_SK_BAS;

   /* load basis information into MOSEK */
   SCIP_CALL( setbase(lpi) );

   SCIPdebugMessage("Store from state into task iter : %d with solsta : %d\n", lpistate->num, lpistate->solsta);

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /**@todo implement SCIPlpiClearState() for MOSEK */
   SCIPmessagePrintWarning(lpi->messagehdlr, "MOSEK interface does not implement SCIPlpiClearState()\n");

   return SCIP_OKAY;
}

/** frees LP state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiFreeState (%d)\n",lpi->lpid);

   if( *lpistate != NULL )
   {
      lpistateFree(lpistate, blkmem);
   }

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("Calling SCIPlpiHasStateBasis (%d)\n",lpi->lpid);
   
   return ( lpistate != NULL && lpistate->num >= 0);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   MOSEK_CALL( MSK_readsolution(lpi->task, MSK_SOL_BAS, fname) );

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("writing LP state to file <%s>\n", fname);
   
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   /* set parameter to be able to write */
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_HEAD, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_VARIABLES, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_SOL_CONSTRAINTS, MSK_ON) );
   
   MOSEK_CALL( MSK_writesolution(lpi->task, MSK_SOL_BAS, fname) );

   return SCIP_OKAY;
}

/*
 * Parameter Methods
 */

/** constant array containing the parameter names */
static const char* paramname[] = {
   "SCIP_LPPAR_FROMSCRATCH",                 /**< solver should start from scratch at next call? */             
   "SCIP_LPPAR_FASTMIP",                     /**< fast mip setting of LP solver */                              
   "SCIP_LPPAR_SCALING",                     /**< should LP solver use scaling? */                              
   "SCIP_LPPAR_PRESOLVING",                  /**< should LP solver use presolving? */                           
   "SCIP_LPPAR_PRICING",                     /**< pricing strategy */                                           
   "SCIP_LPPAR_LPINFO",                      /**< should LP solver output information to the screen? */         
   "SCIP_LPPAR_FEASTOL",                     /**< feasibility tolerance for primal variables and slacks */      
   "SCIP_LPPAR_DUALFEASTOL",                 /**< feasibility tolerance for dual variables and reduced costs */ 
   "SCIP_LPPAR_BARRIERCONVTOL",              /**< convergence tolerance used in barrier algorithm */            
   "SCIP_LPPAR_LOBJLIM",                     /**< lower objective limit */                                      
   "SCIP_LPPAR_UOBJLIM",                     /**< upper objective limit */                                      
   "SCIP_LPPAR_LPITLIM",                     /**< LP iteration limit */                                         
   "SCIP_LPPAR_LPTILIM",                     /**< LP time limit */                                              
   "SCIP_LPPAR_MARKOWITZ",                   /**< Markowitz tolerance */
   "SCIP_LPPAR_ROWREPSWITCH",                /**< simplex algorithm shall use row representation of the basis
                                               *  if number of rows divided by number of columns exceeds this value */
   "SCIP_LPPAR_THREADS"                      /**< number of threads used to solve the LP */
};

/** method mapping parameter index to parameter name */
static 
const char* paramty2str(
   SCIP_LPPARAM          type
   )
{  /*lint --e{641}*/
   /* check if the parameters in this order */
   assert(SCIP_LPPAR_FROMSCRATCH == 0);      /**< solver should start from scratch at next call? */             
   assert(SCIP_LPPAR_FASTMIP == 1);          /**< fast mip setting of LP solver */                              
   assert(SCIP_LPPAR_SCALING == 2);          /**< should LP solver use scaling? */                              
   assert(SCIP_LPPAR_PRESOLVING == 3);       /**< should LP solver use presolving? */                           
   assert(SCIP_LPPAR_PRICING == 4);          /**< pricing strategy */                                           
   assert(SCIP_LPPAR_LPINFO == 5);           /**< should LP solver output information to the screen? */         
   assert(SCIP_LPPAR_FEASTOL == 6);          /**< feasibility tolerance for primal variables and slacks */      
   assert(SCIP_LPPAR_DUALFEASTOL == 7);      /**< feasibility tolerance for dual variables and reduced costs */ 
   assert(SCIP_LPPAR_BARRIERCONVTOL == 8);   /**< convergence tolerance used in barrier algorithm */            
   assert(SCIP_LPPAR_LOBJLIM == 9);          /**< lower objective limit */                                      
   assert(SCIP_LPPAR_UOBJLIM == 10);         /**< upper objective limit */                                      
   assert(SCIP_LPPAR_LPITLIM == 11);         /**< LP iteration limit */                                         
   assert(SCIP_LPPAR_LPTILIM == 12);         /**< LP time limit */                                              
   assert(SCIP_LPPAR_MARKOWITZ == 13);       /**< Markowitz tolerance */
   assert(SCIP_LPPAR_ROWREPSWITCH == 14);    /**< row representation switch */
   assert(SCIP_LPPAR_THREADS == 15);         /**< number of threads used to solve the LP */
                                           
   return paramname[type]; 
}

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
                              )
{  /*lint --e{641}*/
   SCIPdebugMessage("getting int parameter %s\n", paramty2str(type));

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   switch (type)
   {
   case SCIP_LPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_HOTSTART, ival) );
      *ival = (*ival == MSK_SIM_HOTSTART_NONE); 
      break;
   case SCIP_LPPAR_FASTMIP:                   /* fast mip setting of LP solver */
      return  SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:                   /* should LP solver use scaling? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_SCALING, ival) );
      *ival = (*ival != MSK_SCALING_NONE);
      break;
   case SCIP_LPPAR_PRESOLVING:                /* should LP solver use presolving? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, ival) );
      *ival = (*ival != MSK_PRESOLVE_MODE_OFF); 
      break;
   case SCIP_LPPAR_PRICING:                   /* pricing strategy */
      *ival = lpi->pricing; 
      break;
   case SCIP_LPPAR_LPINFO:                    /* should LP solver output information to the screen? */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_LOG, ival) );
      *ival = (*ival == MSK_ON);
      break;
   case SCIP_LPPAR_LPITLIM:                   /* LP iteration limit */
      MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   int scaling;

#if SCIP_CONTROLS_PRICING
   /*lint --e{641}*/
   static int pricing[7] = {
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_FULL,
      MSK_SIM_SELECTION_PARTIAL,
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_ASE,
      MSK_SIM_SELECTION_DEVEX,
   };
#endif

   SCIPdebugMessage("Calling SCIPlpiSetIntpar (%d) Parameter=<%s>  Value=<%d>\n", lpi->lpid, paramty2str(type), ival);

   assert(SCIP_PRICING_LPIDEFAULT == 0);
   assert(SCIP_PRICING_AUTO == 1);
   assert(SCIP_PRICING_FULL == 2);
   assert(SCIP_PRICING_PARTIAL == 3);
   assert(SCIP_PRICING_STEEP == 4);
   assert(SCIP_PRICING_STEEPQSTART == 5);
   assert(SCIP_PRICING_DEVEX == 6);

   SCIPdebugMessage("Calling SCIPlpiSetIntpar (%d) %s = %d\n",lpi->lpid,paramty2str(type),ival);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   switch (type)
   {
   case SCIP_LPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_HOTSTART,
            ival ? MSK_SIM_HOTSTART_NONE : MSK_SIM_HOTSTART_STATUS_KEYS ) );
      break;
   case SCIP_LPPAR_FASTMIP:                   /* fast mip setting of LP solver */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:                   /* should LP solver use scaling? */
      scaling = (ival ? MSK_SCALING_FREE : MSK_SCALING_NONE);
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_SCALING, scaling) );
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_INTPNT_SCALING, scaling) );
      break;
   case SCIP_LPPAR_PRESOLVING:                /* should LP solver use presolving? */
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_PRESOLVE_USE, 
            ival ? MSK_PRESOLVE_MODE_FREE : MSK_PRESOLVE_MODE_OFF) );

#ifdef SCIP_DEBUG
      if( ival )
      {
         SCIPdebugMessage("Setting presolve to on\n");
      }
#endif
      break;
   case SCIP_LPPAR_PRICING:                   /* pricing strategy */
      assert(ival >= 0 && ival <= SCIP_PRICING_DEVEX);
      lpi->pricing = (SCIP_PRICING)ival;

#ifdef SCIP_DEBUG
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_AUTO:
         SCIPdebugMessage("Setting pricing to auto\n");
         break;
      case SCIP_PRICING_FULL:
         SCIPdebugMessage("Setting pricing to full\n");
         break;
      case SCIP_PRICING_PARTIAL:
         SCIPdebugMessage("Setting pricing to partial\n");
         break;
      case SCIP_PRICING_LPIDEFAULT:
         SCIPdebugMessage("Setting pricing to lpi default\n");
         break;
      case SCIP_PRICING_STEEP:
         SCIPdebugMessage("Setting pricing to steep\n");
         break;
      case SCIP_PRICING_STEEPQSTART:
         SCIPdebugMessage("Setting pricing to steep quick start\n");
         break;
      case SCIP_PRICING_DEVEX:
         SCIPdebugMessage("Setting pricing to devex\n");
         break;
      }
#endif

#if SCIP_CONTROLS_PRICING
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_SELECTION, pricing[ival]) );

      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, pricing[ival]) );

      if( !(lpi->pricing == SCIP_PRICING_PARTIAL || lpi->pricing == SCIP_PRICING_AUTO ) )
      {
         /* No restrict */
         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION, 0) );  

         MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, 0) );  
      }
#else
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_PRIMAL_SELECTION, MSK_SIM_SELECTION_FREE) );

      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_DUAL_SELECTION, MSK_SIM_SELECTION_FREE) );
#endif
      break;
   case SCIP_LPPAR_LPINFO:                    
      /* should LP solver output information to the screen? */
#if FORCE_MOSEK_LOG
      SCIPdebugMessage("Ignoring log setting!\n");
#else
      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_LOG, ival ? MSK_ON : MSK_OFF) );
#endif
      break;
   case SCIP_LPPAR_LPITLIM:                   /* LP iteration limit */
#if DEBUG_PARAM_SETTING
      if( ival )
      {
         SCIPdebugMessage("Setting max iter to : %d\n",ival);
      }
#endif

      MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("getting real parameter %s\n", paramty2str(type));
   
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   switch (type)
   {
#if SCIP_CONTROLS_TOLERANCES
   case SCIP_LPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
#endif
   case SCIP_LPPAR_LOBJLIM:                   /* lower objective limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      break;
   case SCIP_LPPAR_UOBJLIM:                   /* upper objective limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:                   /* LP time limit */
      MOSEK_CALL( MSK_getdouparam(lpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:                 /* Markowitz tolerance */
   default:
      return SCIP_PARAMETERUNKNOWN;
   } /*lint !e788*/ 
   
   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("setting real parameter %s to %g\n", paramty2str(type), dval);
   
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   /**@todo Limits shouldn't be hardcoded */
   
   switch (type)
   {
#if SCIP_CONTROLS_TOLERANCES
   case SCIP_LPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      if (dval < 1e-9) 
         dval = 1e-9;
      
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      if (dval < 1e-9) 
         return SCIP_PARAMETERUNKNOWN;
      /*         dval = 1e-9; */

      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
#endif
   case SCIP_LPPAR_LOBJLIM:                   /* lower objective limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      break;
   case SCIP_LPPAR_UOBJLIM:                   /* upper objective limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:                   /* LP time limit */
      MOSEK_CALL( MSK_putdouparam(lpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:                 /* Markowitz tolerance */
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/
   
   return SCIP_OKAY;
}


/*
 * Numerical Methods
 */


/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   return MSK_INFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);
   
   return IS_POSINF(val);
}


/*
 * File Interface Methods
 */


/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int olddataformat;

   SCIPdebugMessage("Calling SCIPlpiReadLP (%d), filename <%s>\n",lpi->lpid, fname);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_readdata(lpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_READ_DATA_FORMAT, olddataformat) );
   
   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int olddataformat;

   SCIPdebugMessage("Calling SCIPlpiReadLP (%d), filename <%s>\n",lpi->lpid, fname);

   assert(MosekEnv != NULL); 
   assert(lpi != NULL); 
   assert(lpi->task != NULL);

   MOSEK_CALL( MSK_getintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_writedata(lpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(lpi->task, MSK_IPAR_WRITE_DATA_FORMAT, olddataformat) );
   
   return SCIP_OKAY;
}
 

