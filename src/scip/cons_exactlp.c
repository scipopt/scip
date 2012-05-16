/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_exactlp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for exactlp constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define BOUNDCHG_OUT   /** uncomment to get more debug msgs (SCIP_DEBUG) about processing bound changes of variables,
//                        *  including pseudo objective value update */
//#define FPRELAXVAL_OUT /** uncomment to get info about relaxing exact values to floating-point values */
//#define FPRELAXLP_OUT  /** uncomment to get info about creation of floating-point relaxation or approximation */
//#define DBAUTO_OUT     /** uncomment to get info about automatic selection of safe dual bounding method */
//#define LPEX_OUT       /** uncomment to get info about solving process of exact LP solver */
//#define BASISCHECK_OUT /** uncomment to get info about LP basis consistency check between exact and inexact LP solvers */
//#define PS_OUT         /** uncomment to get additional info about project-and-shift dual bounding method */

//#define TESTDBVERIFY   /** uncomment to test basis verification method by comparing it to exact LP solver behavior */
//#define USEOBJLIM      /** uncomment to pass objlimit to exact LP solver; same as in lpiex_qso.c
//                        *  warning: QSopt_ex allows objlimits but the support is buggy; if the limit is reached,
//                        *  QSopt_ex does not stop but increases the precision */
//#define USETIMELIM     /** uncomment to pass timelimit to exact LP solver
//                        *  warning: QSopt_ex allows timelimits but the support is buggy; if the limit is reached,
//                        *  QSopt_ex does not stop but increases the precision */

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "scip/cons_exactlp.h"
#include "scip/intervalarith.h"
#include "scip/lp.h"
#include "scip/lpiex.h"
#include "scip/misc.h"
#include "scip/primalex.h"
#include "scip/scip.h"
#include "scip/solex.h"
#include "scip/struct_scip.h"
#include "rectlu/rectlu.h"

#ifdef WITH_EXACTSOLVE

/** @todo exiptodo: mpfr library is needed to control the rounding mode in the gmp-double conversion. these methods might
 *   move to another file. delete library inclusion then.
 */
#include "mpfr.h" /* mpfr.h has to be included before gmp.h */
#include "gmp.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "exactlp"
#define CONSHDLR_DESC          "LP relaxation of a MIP that is given by rational data"
#define CONSHDLR_SEPAPRIORITY    950000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -400000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -400000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */

#define DEFAULT_PSOBJWEIGHT         0.0 /**< weight of the original objective function in lp to compute interior point */
#define DEFAULT_PSREDUCEAUXLP     FALSE /**< should the number of constraints in lp to compute interior point be reduced? */
#define DEFAULT_PSDUALCOLSELECTION  'A' /**< strategy to select which dual columns to use for lp to compute interior point
                                         *   ('n'o sel, 'a'ctive rows of exact primal LP, 'A'ctive rows of inexact primal LP,
                                         *   'b'asic rows of exact primal LP, 'B'asic rows of inexact primal LP)" */
#define DEFAULT_PSINTPOINTSELECTION 'o' /**< method to select interior point ('a'rbitrary interior point, 'o'ptimized interior
                                         *   point, 'A'rbitrary interior point in dual form, 't'wo stage optimized interior point */
#define DEFAULT_PSUSEINTPOINT      TRUE /**< should correction shift use an interior pt? (otherwise use interior ray of recession cone) */
#define DEFAULT_PSINFEASVERSION    TRUE /**< should infeasibility version of project-and-shift be used? */
#define DEFAULT_INTERLEAVEDBFREQ      0 /**< frequency at which dual bounding strategy is interleaved (-1: never, 0: if prommising, x: xth node) */
#define DEFAULT_COMPUTEDBDIFF     FALSE /**< should the quality of the safe dual bounds be computed (involves solving LPs exactly)? */
#define DEFAULT_SMALLDBDIFF        1e-9 /**< maximum value of small relativ difference of safe dual bound */
#define DEFAULT_MEDIUMDBDIFF       1e-3 /**< maximum value of medium relativ difference of safe dual bound */

#define EVENTHDLR_NAME        "exactlp"
#define EVENTHDLR_DESC        "bound change event handler for exactlp constraints"

#define PSPOSTPROCESSDUALSOL       TRUE

#define OBJSCALE_MAXFINALSCALE   1000.0 /**< maximal final value to apply as scaling */
#define PSBIGM                      100
#define PSWARMSTARTAUXPROB         TRUE
#define SPARSEMAXNNONZ            10000 /**< maximal number of nonzero coefficients in matrix to be considered sparse */
#define FEWLBOUNDSRATIO             0.2 /**< maximal percentage of vars with large bounds that is regarded to be small;
                                         *   same as in lp.c */
#define STORE_INTERA                    /**< should interval constraint matrix for provedBoundInterval() be stored? */


/*
 * Data structures
 */

/** @todo exiptodo: posinfinity and neginfinity should be moved to a more appropriate place (scipex.h, def.h or defex.h) */
/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   mpq_t                 posinfinity;        /**< value considered to be infinity */
   mpq_t                 neginfinity;        /**< value considered to be infinity */
   SCIP_LPIEX*           lpiex;              /**< Exact LP solver interface */
   SCIP_Bool             lpexconstructed;    /**< was the exact LP of some prior node already constructed (constraints)? */
   mpq_t                 lpiexuobjlim;       /**< current upper objective limit in LPIEX */
   SCIP_NODE*            lastenfopsnode;     /**< last node at which enfops was called */
   mpq_t                 pseudoobjval;       /**< pseudo solution value at node where all unprocessed bound changes were
                                              *   applied last with all variables set to their best bounds, ignoring
                                              *   variables with infinite best bound */
   int                   pseudoobjvalinf;    /**< number of variables with infinite best bound in pseudo solution at node
                                              *   where all unprocessed bound changes were applied last */
   mpq_t*                interiorpt;         /**< stores S-interior point for root node dual problem */
   mpq_t*                interiorray;        /**< stores S-interior ray for root node dual problem */
   int*                  includedcons;       /**< 1 if constraints dual variable is included in original S-interior point/ray */
   int                   nextendedconss;     /**< dimension of S-interior point/ray = 2*(nvars+nconss) */
   int*                  psbasis;            /**< mapping for basis used in factorization */
   int                   npsbasis;           /**< length of psbasis */
   qsnum_factor_work*    rectfactor;         /**< stores factorized matrix for project-and-shift */
   mpq_t                 commonslack;        /**< slack by which S-interior point/ray satisfies inequalities */
   SCIP_PRIMALEX*        primal;             /**< exact primal data and solution storage */
   SCIP_Bool             psdatacon;          /**< was project-and-shift data structure constructed? */
   SCIP_Bool             psdatafail;         /**< did the construction of the project-and-shift root node data fail? */
   SCIP_Bool             pshaspoint;         /**< has an S-interior point successfully been constructed? */
   SCIP_Bool             pshasray;           /**< has an S-interior ray successfully been constructed? */
   SCIP_Real             psobjweight;        /**< weight of the original objective function in lp to compute interior point */
   SCIP_Bool             psreduceauxlp;      /**< should the number of constraints in lp to compute interior point be reduced? */
   char                  psdualcolselection; /**< strategy to select which dual columns to use for lp to compute interior point
                                              *   ('n'o sel, 'a'ctive rows of exact primal LP, 'A'ctive rows of inexact primal LP) */
   char                  psintpointselection;/**< method to select interior point ('a'rbitrary interior point, 'o'ptimized interior point
                                              *   'A'rbitrary interior point in dual form, 't'wo stage optimized interior point) */
   SCIP_Bool             psuseintpoint;      /**< should correction shift use an interior pt? (otherwise use interior ray of recession cone) */
   SCIP_INTERVAL*        intervalval;        /**< stores interval representation of matrix values interval*/
   SCIP_Bool             intervalvalcon;     /**< was the interval representation of the matrix data constructed? */
   SCIP_Bool*            rootactivity;       /**< which primal constraints are active at the optimal root node LP solution */
   SCIP_Bool             rootactivitycon;    /**< has the rootactivity been stored? */
   SCIP_Bool             psinfeasversion;    /**< should infeasibility version of project-and-shift be used? */
   int                   interleavedbfreq;   /**< frequency at which dual bounding strategy is interleaved (-1: never, 0: if prommising, x: xth node) */
   SCIP_Longint          nprovedfeaslp;      /**< number of times, safe dual bound computation for feasible LPs was called */
   SCIP_Longint          nprovedinfeaslp;    /**< number of times, safe verification for infeasible LPs was called */
   SCIP_Longint          nfailprovedfeaslp;  /**< number of times, safe dual bound computation for feasible LPs failed */
   SCIP_Longint          nfailprovedinfeaslp;/**< number of times, safe verification for infeasible LPs failed */
   SCIP_Longint          nabortprovedinfeaslp;/**< number of times, safe verification for infeasible LPs was aborted
                                               *  because no approximate certificate was given */
   SCIP_CLOCK*           provedfeaslptime;   /**< time needed for safe dual bound computation for feasible LPs */
   SCIP_CLOCK*           provedinfeaslptime; /**< time needed for safe verification for infeasible LPs */
   SCIP_Longint          nexactfeaslp;       /**< number of times, exact integrality verification for feasible LPs was called */
   SCIP_Longint          nexactinfeaslp;     /**< number of times, exact verification for infeasible LPs was called */
   SCIP_Longint          nexactunsollp;      /**< number of times, exact LP solver is called for unsolved LP */
   SCIP_Longint          nwrongexactfeaslp;  /**< number of times, claimed integral LP solution was fractional */
   SCIP_Longint          nwronginfeasexactfeaslp; /**< number of times, claimed integral LP solution was infeasible and even the LP was infeasible */
   SCIP_Longint          nwrongexactinfeaslp;/**< number of times, claimed infeasible LP was feasible */
   SCIP_CLOCK*           exactfeaslptime;    /**< time needed for exact integrality verification for feasible LPs */
   SCIP_CLOCK*           exactinfeaslptime;  /**< time needed for exact verification for infeasible LPs */
   SCIP_CLOCK*           exactunsollptime;   /**< time needed for exact LP solve for unsolved LP */
   SCIP_Longint          nzerodbdiff;        /**< number of safe dual bounds with zero relativ difference */
   SCIP_Longint          nsmalldbdiff;       /**< number of safe dual bounds with small relativ difference */
   SCIP_Longint          nmediumdbdiff;      /**< number of safe dual bounds with medium relativ difference */
   SCIP_Longint          nlargedbdiff;       /**< number of safe dual bounds with large relativ difference */
   SCIP_Bool             computedbdiff;      /**< should the quality of the safe dual bounds be computed (involves solving LPs exactly)? */
   SCIP_Real             smalldbdiff;        /**< maximum value of small relativ difference of safe dual bound */
   SCIP_Real             mediumdbdiff;       /**< maximum value of medium relativ difference of safe dual bound */
};

/** constraint data for exactlp constraints */
struct SCIP_ConsData
{
   SCIP_OBJSENSE         objsense;           /**< objective sense */
   mpq_t                 objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   mpq_t                 objscale;           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   int                   nvars;              /**< number of variables */
   int                   ninfbounds;         /**< number of variables with infinite bound in safe dual bounding method */
   int                   ninfintbounds;      /**< number of integer variables with infinite bound in safe db method */
   int                   nlargebounds;       /**< number of variables with large bound in safe dual bounding method */
   mpq_t*                obj;                /**< objective function values of variables */
   SCIP_Bool             objneedscaling;     /**< do objective function values need to be scaled because some are not FP representable? */
   mpq_t*                lb;                 /**< lower bounds of variables */
   mpq_t*                ub;                 /**< upper bounds of variables */
   mpq_t*                lbloc;              /**< local lower bounds of variables at node where all unprocessed
                                              *   bound changes were applied last */
   mpq_t*                ubloc;              /**< local upper bounds of variables when exact LP was last updated */
   int*                  bndchglbpos;        /**< positions of variables in unprocessed lower bound change arrays, or -1 */
   int*                  bndchgubpos;        /**< positions of variables in unprocessed upper bound change arrays, or -1 */
   int*                  lockdown;           /**< rounding down locks of variables */
   int*                  lockup;             /**< rounding up locks of variables */
   int                   nconss;             /**< number of constraints in exactlp data structure */
   int                   nsplitconss;        /**< number of constraints that were split in exactlp data structure */
   int                   conssize;           /**< size of constraint specific array in exactlp data structure */
   mpq_t*                lhs;                /**< left hand sides of constraints */
   mpq_t*                rhs;                /**< right hand sides of constraints */
   int                   nnonz;              /**< number of nonzero elements in the constraint matrix */
   int                   nintegral;          /**< number of integral nonzero elements in the constraint matrix */
   int*                  beg;                /**< start index of each constraint in ind and val array */
   int*                  len;                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind;                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val;                /**< values of nonzero constraint matrix entries (and some zeros) */
   mpq_t                 minabsval;          /**< minimum absolute nonzero constraint matrix, lhs, or rhs entry */
   mpq_t                 maxabsval;          /**< maximum absolute nonzero constraint matrix, lhs, or rhs entry */
   SCIP_ROW**            rows;               /**< rows for LP relaxation/approximation (FP data) of exactlp constraint */
   SCIP_EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   int                   nbndchglb;          /**< number of variables with unprocessed lower bound changes */
   int                   nbndchgub;          /**< number of variables with unprocessed upper bound changes */
   mpq_t*                bndchglb;           /**< local lower bounds of variables with unprocessed bound changes */
   mpq_t*                bndchgub;           /**< local upper bounds of variables with unprocessed bound changes */
   int*                  bndchglbind;        /**< indices of variables with unprocessed lower bound changes */
   int*                  bndchgubind;        /**< indices of variables with unprocessed upper bound changes */
};

/** event data for bound change event */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< exactlp constraint data to process the bound change for */
   int                   varind;             /**< index of variable in exactlp consdata arrays (lb, ub, obj) */
   int                   filterpos;          /**< position of event in variable's event filter */
};




/*
 * Local methods
 */


/*
 * local methods for consistency checks
 */

#ifndef NDEBUG
/** checks for correct status of original variables */
static
void checkOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** origvars;
   int i;

   origvars = SCIPgetOrigVars(scip);

   for( i = 0; i < SCIPgetNOrigVars(scip); ++i )
   {
      assert(SCIPvarIsOriginal(origvars[i]));
   }
}
#endif

#ifndef NDEBUG
/** checks whether basis of exact and inexact LP solvers match after basis was loaded from one solver to the other */
static
SCIP_RETCODE checkLoadState(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata       /**< exactlp constraint handler data */
   )
{
   int* cstat;
   int* rstat;
   int* cstatex;
   int* rstatex;
   int i;
   int ncols;
   int nrows;
   int ncolsex;
   int nrowsex;

   ncols = SCIPgetNLPCols(scip);
   nrows = SCIPgetNLPRows(scip);
   SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncolsex);
   SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrowsex);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cstat, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rstat, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cstatex, ncolsex) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rstatex, nrowsex) );

   /** @todo exiptodo: using scip->lp is not SCIP conform.
    *  - implement method corresponding to SCIPlpGetBase() in scip.c
    *  - there is at least one other position where lp is accessed directly
    */
   /* get exact and inexact basis */
   SCIP_CALL( SCIPlpGetBase(scip->lp, cstat, rstat) );
   SCIP_CALL( SCIPlpiexGetBase(conshdlrdata->lpiex, cstatex, rstatex) );

#ifdef BASISCHECK_OUT
   printf("basis of exact and inexact LP solver should match after basis was loaded from one to the other:\n");
   printf("  inexact LP basis (ncols<%d>, nrows<%d>):\n", ncols, nrows);
   //printCurrentLP(scip); /* uncomment to display inexact LP */
   for( i = 0; i < ncols; ++i)
   {
      if( cstat[i] == SCIP_BASESTAT_LOWER)
         printf("     c=%d: %d (at lower)\n", i, cstat[i]);
      else if( cstat[i] == SCIP_BASESTAT_BASIC)
         printf("     c=%d: %d (basic)\n", i, cstat[i]);
      else if( cstat[i] ==  SCIP_BASESTAT_UPPER)
         printf("     c=%d: %d (at upper)\n", i, cstat[i]);
      else if( cstat[i] ==  SCIP_BASESTAT_ZERO)
         printf("     c=%d: %d (at zero)\n", i, cstat[i]);
      else
         printf("     c=%d: %d (error)\n", i, cstat[i]);
   }
   for( i = 0; i < nrows; ++i)
   {
      if( rstat[i] == SCIP_BASESTAT_LOWER)
         printf("     r=%d: %d (at lower)\n", i, rstat[i]);
      else if( rstat[i] == SCIP_BASESTAT_BASIC)
         printf("     r=%d: %d (basic)\n", i, rstat[i]);
      else if( rstat[i] ==  SCIP_BASESTAT_UPPER)
         printf("     r=%d: %d (at upper)\n", i, rstat[i]);
      else
         printf("     r=%d: %d (error)\n", i, rstat[i]);
   }
   printf("  exact LP basis (ncols<%d>, nrows<%d>):\n", ncolsex, nrowsex);
#endif

   /* check basis status of columns */
   assert(ncols == ncolsex);
   for( i = 0; i < ncolsex; ++i)
   {
#ifdef BASISCHECK_OUT
      if( cstatex[i] == SCIP_BASESTAT_LOWER)
         printf("     ex c=%d: %d (at lower)\n", i, cstatex[i]);
      else if( cstatex[i] == SCIP_BASESTAT_BASIC)
         printf("     ex c=%d: %d (basic)\n", i, cstatex[i]);
      else if( cstatex[i] ==  SCIP_BASESTAT_UPPER)
         printf("     ex c=%d: %d (at upper)\n", i, cstatex[i]);
      else if( cstatex[i] ==  SCIP_BASESTAT_ZERO)
         printf("     ex c=%d: %d (at zero)\n", i, cstatex[i]);
      else
         printf("     ex c=%d: %d (error)\n", i, cstatex[i]);
#endif

      assert(cstatex[i] == cstat[i]);
   }

   /* check basis status of row slacks */
   assert(nrows == nrowsex);
   for( i = 0; i < nrowsex; ++i)
   {
#ifdef BASISCHECK_OUT
      if( rstatex[i] == SCIP_BASESTAT_LOWER)
         printf("     ex r=%d: %d (at lower)\n", i, rstatex[i]);
      else if( rstatex[i] == SCIP_BASESTAT_BASIC)
         printf("     ex r=%d: %d (basic)\n", i, rstatex[i]);
      else if( rstatex[i] ==  SCIP_BASESTAT_UPPER)
         printf("     ex r=%d: %d (at upper)\n", i, rstatex[i]);
      else
         printf("     ex r=%d: %d (error)\n", i, rstatex[i]);
#endif

      assert(rstatex[i] == rstat[i]);
   }

#if 0 /* activate this for detailed debugging. writes both LP basis to files and stops afterwards */
   {
      SCIP_LPI* lpi;

      SCIPgetLPI(scip, &lpi);
      SCIP_CALL( SCIPlpiWriteState(lpi, "testset/state_lp.bas") );
      SCIP_CALL( SCIPlpiexWriteState(conshdlrdata->lpiex, "testset/state_lpex.bas") );
      SCIPABORT();
   }
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &rstatex);
   SCIPfreeBufferArray(scip, &cstatex);
   SCIPfreeBufferArray(scip, &rstat);
   SCIPfreeBufferArray(scip, &cstat);

   return SCIP_OKAY;
}
#endif

#ifndef NDEBUG
/** checks whether pseudo objective value has been updated correctly */
static
void checkPseudoobjval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   mpq_t pseudoobjval;
   int pseudoobjvalinf;
   mpq_t prod;
   int i;

   mpq_init(pseudoobjval);
   mpq_set_d(pseudoobjval, 0.0);
   pseudoobjvalinf = 0;
   mpq_init(prod);

   /* compute pseudo objective value manually */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( mpq_sgn(consdata->obj[i]) > 0 )
      {
         if( isNegInfinity(conshdlrdata, consdata->lbloc[i]) )
            pseudoobjvalinf++;
         else
         {
            /* pseudoobjval += lb * obj */
            mpq_mul(prod, consdata->lbloc[i], consdata->obj[i]);
            mpq_add(pseudoobjval, pseudoobjval, prod);
         }
      }
      else if( mpq_sgn(consdata->obj[i]) < 0 )
      {
         if( isPosInfinity(conshdlrdata, consdata->ubloc[i]) )
            pseudoobjvalinf++;
         else
         {
            /* pseudoobjval += ub * obj */
            mpq_mul(prod,  consdata->ubloc[i], consdata->obj[i]);
            mpq_add(pseudoobjval, pseudoobjval, prod);
         }
      }
   }

   /* check whether manually computed pseudo objective value is equal to updated one */
   assert(pseudoobjvalinf == conshdlrdata->pseudoobjvalinf);
   assert(mpq_equal(pseudoobjval, conshdlrdata->pseudoobjval) != 0);

   mpq_clear(pseudoobjval);
   mpq_clear(prod);
}
#endif


/*
 * local methods for working with rational numbers
 */

/** returns value treated as negative infinite in exactlp constraint handler */
const mpq_t* negInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   )
{  /*lint --e{715} */
   return (const mpq_t*) (&conshdlrdata->neginfinity);
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
const mpq_t* posInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   )
{  /*lint --e{715} */
   return (const mpq_t*) (&conshdlrdata->posinfinity);
}

/** checks if value is treated as negative infinite in exactlp constraint handler */
SCIP_Bool isNegInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   const mpq_t           val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsNegInfinity(conshdlrdata->lpiex, val);
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
SCIP_Bool isPosInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   const mpq_t           val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsPosInfinity(conshdlrdata->lpiex, val);
}

/** returns whether given rational number can be stored as FP number without roundinf errors */
SCIP_Bool mpqIsReal(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val                 /**< given rational number */
   )
{
   SCIP_Bool result;
   mpq_t tmp;

   mpq_init(tmp);

   mpq_set_d(tmp, mpqGetRealApprox(scip, val));

   if( mpq_equal(tmp, val) == 0 )
      result = FALSE; /* approx(a) =/= a */
   else
      result = TRUE;

#ifndef NDEBUG
   {
      SCIP_Bool resultnew;

      mpq_set_d(tmp, mpqGetRealRelax(scip, val, GMP_RNDU));
      if( mpq_equal(tmp, val) == 0 )
         resultnew = FALSE; /* approx(a) =/= a */
      else
         resultnew = TRUE;
      assert(result == resultnew);

      mpq_set_d(tmp, mpqGetRealRelax(scip, val, GMP_RNDD));
      if( mpq_equal(tmp, val) == 0 )
         resultnew = FALSE; /* approx(a) =/= a */
      else
         resultnew = TRUE;
      assert(result == resultnew);
   }
#endif

   mpq_clear(tmp);

   return result;
}

/** converts given rational number into an FP number; uses given rounding mode during conversion
 * (should be used to construct an FP relaxation of a constraint)
 */
SCIP_Real mpqGetRealRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const mpq_t           val,                /**< given rational number */
   mp_rnd_t              roundmode           /**< rounding mode to be used for the conversion */
   )
{
   SCIP_Real valrelax;

   mpfr_t valmpfr;

   mpfr_init_set_q(valmpfr, val, roundmode);
   valrelax = (SCIP_Real) mpfr_get_d(valmpfr, roundmode);
   mpfr_clear(valmpfr);

#ifndef NDEBUG
   {
      mpq_t result;

      mpq_init(result);
      mpq_set_d(result, valrelax);

      if( roundmode == GMP_RNDU )
         assert(mpq_cmp(result, val) >= 0);
      if( roundmode == GMP_RNDD )
         assert(mpq_cmp(result, val) <= 0);
#ifdef FPRELAXVAL_OUT
      if( roundmode == GMP_RNDU )
      {
         gmp_printf("   relaxing value: given<%Qd> <= relaxed<%Qd|", val, result);
         printf("%g>\n", mpq_get_d(result));
      }
      if( roundmode == GMP_RNDD )
      {
         gmp_printf("   relaxing value: given<%Qd> >= relaxed<%Qd|", val, result);
         printf("%g>\n", mpq_get_d(result));
      }
#endif
      mpq_clear(result);
   }
#endif

   /** @todo exiptodo: after reworking the handling of (different) infinity values, we should check the treatment of
    *  infinite values here, in particular, when constructing an FP relaxation
    */
   if( SCIPisInfinity(scip, valrelax) )
      valrelax = SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -valrelax) )
      valrelax = -SCIPinfinity(scip);

   return valrelax;
}

/** converts given rational number into an FP number; uses default rounding mode during conversion
 * (should be used to construct an FP approximation of a constraint)
 */
SCIP_Real mpqGetRealApprox(
   SCIP*                 scip,               /**< SCIP data structure */
   const mpq_t           val                 /**< given rational number */
   )
{
   SCIP_Real valapprox;

   valapprox = (SCIP_Real) mpq_get_d(val);

   /** @todo exiptodo: after reworking the handling of (different) infinity values, we should check the treatment of
    *  infinite values here
    */
   if( SCIPisInfinity(scip, valapprox) )
      valapprox = SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -valapprox) )
      valapprox = -SCIPinfinity(scip);

   return valapprox;
}

/** returns rational number rounded down */
static
void mpqFloor(
   mpq_t                 floorval,           /**< pointer to store rational number rounded down */
   const mpq_t           val                 /**< given rational number */
   )
{
   mpz_t floorint;

   mpz_init(floorint);
   mpz_fdiv_q(floorint, mpq_numref(val), mpq_denref(val));
   mpq_set_z(floorval, floorint);
   mpz_clear(floorint);
}

/** returns rational number rounded up */
static
void mpqCeil(
   mpq_t                 ceilval,             /**< pointer to store rational number rounded up */
   const mpq_t           val                  /**< given rational number */
   )
{
   mpz_t ceilint;

   mpz_init(ceilint);
   mpz_cdiv_q(ceilint, mpq_numref(val), mpq_denref(val));
   mpq_set_z(ceilval, ceilint);
   mpz_clear(ceilint);
}

/** checks, if value is integral */
SCIP_Bool mpqIsIntegral(
   const mpq_t           val                  /**< value to process */
   )
{
   return (mpz_cmp_ui(mpq_denref(val), 1) == 0);
}

/*
 * local methods for managing event data
 */

/** creates event data */
static
SCIP_RETCODE eventdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata,          /**< pointer to store event data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   ind                 /**< index of variable in consdata arrays (lb, ub, obj) to catch events for */
   )
{
   assert(eventdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, eventdata) );
   (*eventdata)->consdata = consdata;
   (*eventdata)->varind = ind;

   return SCIP_OKAY;
}

/** frees event data */
static
SCIP_RETCODE eventdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata           /**< pointer to event data */
   )
{
   assert(eventdata != NULL);

   SCIPfreeBlockMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** catches bound change events for variables in exactlp constraint */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   int i;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );

   /* get transformed problem variables */
   origvars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);
   SCIP_CALL( SCIPgetTransformedVars(scip, consdata->nvars, origvars, vars) );

   for( i = 0; i < consdata->nvars; i++)
   {
      assert(SCIPvarIsOriginal(origvars[i]));
      assert(SCIPvarIsTransformed(vars[i]));
      assert(consdata->eventdatas[i] == NULL);

      SCIP_CALL( eventdataCreate(scip, &consdata->eventdatas[i], consdata, i) );
      SCIP_CALL( SCIPcatchVarEvent(scip, vars[i], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, consdata->eventdatas[i],
            &consdata->eventdatas[i]->filterpos) );
   }

#ifndef NDEBUG
   /* check for correct status of original variables */
   checkOrigVars(scip);
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** drops bound change events for variables in exactlp constraint */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   int i;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );

   /* get transformed problem variables */
   origvars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);
   SCIP_CALL( SCIPgetTransformedVars(scip, consdata->nvars, origvars, vars) );

   for( i = 0; i < consdata->nvars; i++)
   {
      assert(SCIPvarIsOriginal(origvars[i]));
      assert(SCIPvarIsTransformed(vars[i]));
      assert(consdata->eventdatas[i] != NULL);

      SCIP_CALL( SCIPdropVarEvent(scip, vars[i], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, consdata->eventdatas[i],
            consdata->eventdatas[i]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdatas[i]) );
   }

#ifndef NDEBUG
   /* check for correct status of original variables */
   checkOrigVars(scip);
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * local methods for managing constraint handler data and constraint data
 */

/** creates constaint handler data for exactlp constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* initialize exact LP solver */
   SCIPlpiexStart();

   /* open exact LP Solver interface */
   SCIP_CALL( SCIPlpiexCreate(&(*conshdlrdata)->lpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

   /* open exact LP Solver interface for factorization */

   (*conshdlrdata)->rectfactor = (qsnum_factor_work*) NULL;
   (*conshdlrdata)->nextendedconss = 0;

   mpq_init((*conshdlrdata)->commonslack);

   mpq_init((*conshdlrdata)->posinfinity);
   mpq_init((*conshdlrdata)->neginfinity);
   mpq_init((*conshdlrdata)->lpiexuobjlim);

   SCIPlpiexPosInfinity((*conshdlrdata)->lpiex, &(*conshdlrdata)->posinfinity);
   SCIPlpiexNegInfinity((*conshdlrdata)->lpiex, &(*conshdlrdata)->neginfinity);
   SCIPlpiexPosInfinity((*conshdlrdata)->lpiex, &(*conshdlrdata)->lpiexuobjlim);

#ifdef USEOBJLIM
   SCIP_CALL( SCIPlpiexSetRealpar((*conshdlrdata)->lpiex, SCIP_LPPAR_UOBJLIM, (*conshdlrdata)->lpiexuobjlim) );
#endif

   (*conshdlrdata)->lpexconstructed = FALSE;
   (*conshdlrdata)->psdatacon = FALSE;
   (*conshdlrdata)->lastenfopsnode = NULL;
   (*conshdlrdata)->intervalvalcon = FALSE;

   mpq_init((*conshdlrdata)->pseudoobjval);
   mpq_set_d((*conshdlrdata)->pseudoobjval, 0.0);
   (*conshdlrdata)->pseudoobjvalinf = 0;

   (*conshdlrdata)->interiorpt = NULL;
   (*conshdlrdata)->interiorray = NULL;
   (*conshdlrdata)->includedcons = NULL;
   (*conshdlrdata)->psbasis = NULL;
   (*conshdlrdata)->rootactivity = NULL;
   (*conshdlrdata)->rootactivitycon = FALSE;

   /* get event handler for updating bounds of variables in the exact LP */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for exactlp constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   (*conshdlrdata)->primal = NULL;
   (*conshdlrdata)->psdatafail = FALSE;
   (*conshdlrdata)->pshaspoint = FALSE;
   (*conshdlrdata)->pshasray = FALSE;

   (*conshdlrdata)->nprovedfeaslp = 0;
   (*conshdlrdata)->nprovedinfeaslp = 0;
   (*conshdlrdata)->nprovedinfeaslp = 0;
   (*conshdlrdata)->nfailprovedfeaslp = 0;
   (*conshdlrdata)->nfailprovedinfeaslp = 0;
   (*conshdlrdata)->nabortprovedinfeaslp = 0;
   (*conshdlrdata)->provedfeaslptime = NULL;
   (*conshdlrdata)->provedinfeaslptime = NULL;
   (*conshdlrdata)->nexactfeaslp = 0;
   (*conshdlrdata)->nexactinfeaslp = 0;
   (*conshdlrdata)->nexactunsollp = 0;
   (*conshdlrdata)->nwrongexactfeaslp = 0;
   (*conshdlrdata)->nwronginfeasexactfeaslp = 0;
   (*conshdlrdata)->nwrongexactinfeaslp = 0;
   (*conshdlrdata)->exactfeaslptime = NULL;
   (*conshdlrdata)->exactinfeaslptime = NULL;
   (*conshdlrdata)->exactunsollptime = NULL;
   (*conshdlrdata)->nzerodbdiff = 0;
   (*conshdlrdata)->nsmalldbdiff = 0;
   (*conshdlrdata)->nmediumdbdiff = 0;
   (*conshdlrdata)->nlargedbdiff = 0;

   return SCIP_OKAY;
}

/** frees constraint handler data for exactlp constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   if( (*conshdlrdata)->nextendedconss > 0 )
   {
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->psbasis);
      if( (*conshdlrdata)->interiorpt != NULL )
      {
         for( i = 0; i < (*conshdlrdata)->nextendedconss; i++ )
            mpq_clear((*conshdlrdata)->interiorpt[i]);
         SCIPfreeMemoryArray(scip, &(*conshdlrdata)->interiorpt);
      }
      if( (*conshdlrdata)->interiorray != NULL )
      {
         for( i = 0; i < (*conshdlrdata)->nextendedconss; i++ )
            mpq_clear((*conshdlrdata)->interiorray[i]);
         SCIPfreeMemoryArray(scip, &(*conshdlrdata)->interiorray);
      }
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->includedcons);
   }

   if( (*conshdlrdata)->intervalvalcon != FALSE )
   {
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->intervalval);
   }

   if( (*conshdlrdata)->rootactivitycon != FALSE )
   {
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->rootactivity);
   }

   if( (*conshdlrdata)->rectfactor != NULL)
      RECTLUfreeFactorization((*conshdlrdata)->rectfactor);

   if( (*conshdlrdata)->lpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*conshdlrdata)->lpiex) );
   }
   assert((*conshdlrdata)->lpiex == NULL);

   mpq_clear((*conshdlrdata)->pseudoobjval);
   mpq_clear((*conshdlrdata)->commonslack);
   mpq_clear((*conshdlrdata)->lpiexuobjlim);
   mpq_clear((*conshdlrdata)->posinfinity);
   mpq_clear((*conshdlrdata)->neginfinity);

   SCIPfreeMemory(scip, conshdlrdata);

   /* deinitialize exact LP solver */
   SCIPlpiexEnd();

   return SCIP_OKAY;
}

/** gets number of LP rows needed for the LP relaxation of the exactlp constraint */
static
int consdataGetNRows(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   return consdata->nconss;
}

/** updates current pseudo objective values for a change in a variable's objective value or bounds */
static
SCIP_RETCODE updateVar(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   const mpq_t           oldobj,             /**< old objective value of variable */
   const mpq_t           oldlb,              /**< old objective value of variable */
   const mpq_t           oldub,              /**< old objective value of variable */
   mpq_t                 newobj,             /**< new objective value of variable */
   mpq_t                 newlb,              /**< new objective value of variable */
   mpq_t                 newub               /**< new objective value of variable */
   )
{
   mpq_t deltaval;
   mpq_t prod;
   int deltainf;

   assert(conshdlrdata->pseudoobjvalinf >= 0);

   mpq_init(deltaval);
   mpq_init(prod);
   mpq_set_d(deltaval, 0.0);
   deltainf = 0;

#ifdef BOUNDCHG_OUT
   gmp_printf("   change pseudoobjval for obj (%Qd->%Qd) lb (%Qd->%Qd) ub (%Qd->%Qd): pseudoobjval %Qd (inf<%d>) --> ",
      oldobj, newobj, oldlb, newlb, oldub, newub, conshdlrdata->pseudoobjval, conshdlrdata->pseudoobjvalinf);
#endif

   /* subtract old pseudo objective value */
   if( mpq_sgn(oldobj) > 0 )
   {
      if( isNegInfinity(conshdlrdata, oldlb) )
         deltainf--;
      else
      {
         /* deltaval -= oldlb * oldobj */
         mpq_mul(prod, oldlb, oldobj);
         mpq_sub(deltaval, deltaval, prod);
      }
   }
   else if( mpq_sgn(oldobj) < 0 )
   {
      if( isPosInfinity(conshdlrdata, oldub) )
         deltainf--;
      else
      {
         /* deltaval -= oldub * oldobj */
         mpq_mul(prod, oldub, oldobj);
         mpq_sub(deltaval, deltaval, prod);
      }
   }

   /* add new pseudo objective value */
   if( mpq_sgn(newobj) > 0 )
   {
      if( isNegInfinity(conshdlrdata, newlb) )
         deltainf++;
      else
      {
         /* deltaval += newlb * newobj */
         mpq_mul(prod, newlb, newobj);
         mpq_add(deltaval, deltaval, prod);
      }
   }
   else if( mpq_sgn(newobj) < 0 )
   {
      if( isPosInfinity(conshdlrdata, newub) )
         deltainf++;
      else
      {
         /* deltaval += newub * newobj */
         mpq_mul(prod, newub, newobj);
         mpq_add(deltaval, deltaval, prod);
      }
   }

   /* update the pseudo objective values */
   mpq_add(conshdlrdata->pseudoobjval, conshdlrdata->pseudoobjval, deltaval);
   conshdlrdata->pseudoobjvalinf += deltainf;

#ifdef BOUNDCHG_OUT
   gmp_printf("%Qd (inf<%d>)\n", conshdlrdata->pseudoobjval, conshdlrdata->pseudoobjvalinf);
#endif

   assert(conshdlrdata->pseudoobjvalinf >= 0);

   mpq_clear(deltaval);
   mpq_clear(prod);

   return SCIP_OKAY;
}

/** creates exactlp constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   int                   ninfbounds,         /**< number of variables with infinite bound in safe dual bounding method */
   int                   ninfintbounds,      /**< number of integer variables with infinite bound in safe db method */
   int                   nlargebounds,       /**< number of variables with large bound in safe dual bounding method */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,             /**< number of constraints in exactlp data structure */
   int                   nsplitconss,        /**< number of constraints that were split in exactlp data structure */
   int                   conssize,           /**< size of constraint specific array in exactlp data structure */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int                   nintegral,          /**< number of integral nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each constraint in ind and val array */
   int*                  len,                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries (and some zeros) */
   mpq_t                 minabsval,          /**< minimum absolute nonzero constraint matrix, lhs, or rhs entry */
   mpq_t                 maxabsval,          /**< maximum absolute nonzero constraint matrix, lhs, or rhs entry */
   SCIP_Bool             objneedscaling      /**< do objective function values need to be scaled because some are not FP representable? */
   )
{
   int j;

   assert(consdata != NULL);
   assert(nvars > 0 || (nconss == 0 && nconss == 0) );
   assert(nconss > 0 || (nvars >= 0 && nnonz == 0) );
   assert(nnonz > 0 || (nconss >= 0 && nvars >= 0) );
   assert(nconss <= conssize);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* store variable specific information */
   if( nvars > 0 )
   {
      mpq_t tmpzero;

      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->obj, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lb, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ub, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lbloc, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ubloc, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lockdown, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lockup, nvars) );
      mpq_init(tmpzero);
      mpq_set_d(tmpzero, 0.0);

      /* reset pseudo objective value */
      mpq_set_d(conshdlrdata->pseudoobjval, 0.0);
      conshdlrdata->pseudoobjvalinf = 0;
      for( j = 0; j < nvars; ++j )
      {
         mpq_init((*consdata)->obj[j]);
         mpq_init((*consdata)->lb[j]);
         mpq_init((*consdata)->ub[j]);
         mpq_init((*consdata)->lbloc[j]);
         mpq_init((*consdata)->ubloc[j]);

         mpq_set((*consdata)->obj[j], obj[j]);
         mpq_set((*consdata)->lb[j], lb[j]);
         mpq_set((*consdata)->ub[j], ub[j]);
         mpq_set((*consdata)->lbloc[j], lb[j]);
         mpq_set((*consdata)->ubloc[j], ub[j]);

         /* initialize pseudo objective value */
         SCIP_CALL( updateVar(conshdlrdata, tmpzero, tmpzero, tmpzero, (*consdata)->obj[j], (*consdata)->lbloc[j],
               (*consdata)->ubloc[j]) );
      }
      mpq_clear(tmpzero);
      BMSclearMemoryArray((*consdata)->lockdown, nvars);
      BMSclearMemoryArray((*consdata)->lockup, nvars);
   }
   else
   {
      (*consdata)->obj = NULL;
      (*consdata)->lb = NULL;
      (*consdata)->ub = NULL;
      (*consdata)->lbloc = NULL;
      (*consdata)->ubloc = NULL;
      (*consdata)->lockdown = NULL;
      (*consdata)->lockup = NULL;
   }

   /* store constraint specific information */
   if( conssize > 0 )
   {
      /* allocate and copy integer array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->beg, beg, conssize + 1) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->len, len, conssize) );

      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lhs, conssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->rhs, conssize) );

      for( j = 0; j < conssize; ++j )
      {
         mpq_init((*consdata)->lhs[j]);
         mpq_init((*consdata)->rhs[j]);

         if( j < nconss )
         {
            mpq_set((*consdata)->lhs[j], lhs[j]);
            mpq_set((*consdata)->rhs[j], rhs[j]);
         }
      }
   }
   else
   {
      (*consdata)->lhs = NULL;
      (*consdata)->rhs = NULL;
      (*consdata)->beg = NULL;
      (*consdata)->len = NULL;
   }

   /* store matrix specific information */
   if( nnonz > 0 )
   {
      /* allocate and copy integer array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->ind, ind, nnonz) );

      /* allocate, initialize, and copy rational array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->val, nnonz) );

      for( j = 0; j < nnonz; ++j )
      {
         mpq_init((*consdata)->val[j]);

         mpq_set((*consdata)->val[j], val[j]);
      }
   }
   else
   {
      (*consdata)->ind = NULL;
      (*consdata)->val = NULL;
   }

   (*consdata)->objsense = objsense;
   mpq_init((*consdata)->objscale);
   mpq_set_d((*consdata)->objscale, 1.0);
   mpq_init((*consdata)->objoffset);
   mpq_set_d((*consdata)->objoffset, 0.0);
   (*consdata)->nvars = nvars;
   (*consdata)->ninfbounds = ninfbounds;
   (*consdata)->ninfintbounds = ninfintbounds;
   (*consdata)->nlargebounds = nlargebounds;
   (*consdata)->nconss = nconss;
   (*consdata)->nsplitconss = nsplitconss;
   (*consdata)->conssize = conssize;
   (*consdata)->nnonz = nnonz;
   (*consdata)->nintegral = nintegral;
   (*consdata)->rows = NULL;
   (*consdata)->eventdatas = NULL;
   (*consdata)->nbndchglb = 0;
   (*consdata)->nbndchgub = 0;
   (*consdata)->objneedscaling = objneedscaling;
   mpq_init((*consdata)->minabsval);
   mpq_set((*consdata)->minabsval, minabsval);
   mpq_init((*consdata)->maxabsval);
   mpq_set((*consdata)->maxabsval, maxabsval);

   /* catch events for variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      /* allocate bound change event specific data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchglbpos, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchgubpos, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchglb, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchgub, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchglbind, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->bndchgubind, nvars) );
      assert((*consdata)->eventdatas != NULL);
      BMSclearMemoryArray((*consdata)->eventdatas, nvars);
      for( j = 0; j < nvars; ++j )
      {
         mpq_init((*consdata)->bndchglb[j]);
         mpq_init((*consdata)->bndchgub[j]);

         (*consdata)->bndchglbpos[j] = -1;
         (*consdata)->bndchgubpos[j] = -1;
      }

      /* catch events */
      SCIP_CALL( catchEvents(scip, *consdata, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->rows != NULL )
   {
      int nrows;
      int r;

      nrows = consdataGetNRows(consdata);

      for( r = nrows-1; r >= 0; --r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, nrows);
   }

   return SCIP_OKAY;
}

/** frees exactlp constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int j;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release and free the rows */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   /* drop bound change events for variables in exactlp constraints */
   if( (*consdata)->eventdatas != NULL )
   {
      SCIP_CALL( dropEvents(scip, *consdata, eventhdlr) );

      for( j = 0; j < (*consdata)->nvars; ++j )
      {
         mpq_clear((*consdata)->bndchglb[j]);
         mpq_clear((*consdata)->bndchgub[j]);
      }

      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgubind, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglbind, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgub, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglb, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgubpos, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglbpos, (*consdata)->nvars);
   }

   mpq_clear((*consdata)->objoffset);
   mpq_clear((*consdata)->objscale);

   mpq_clear((*consdata)->minabsval);
   mpq_clear((*consdata)->maxabsval);

   /* free matrix specific information */
   for( j = 0; j < (*consdata)->nnonz; ++j )
   {
      mpq_clear((*consdata)->val[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ind, (*consdata)->nnonz);

   /* free constraint specific information */
   for( j = 0; j < (*consdata)->conssize; ++j )
   {
      mpq_clear((*consdata)->rhs[j]);
      mpq_clear((*consdata)->lhs[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->rhs, (*consdata)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lhs, (*consdata)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->len, (*consdata)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->beg, (*consdata)->conssize+1);

   /* free variable specific information */
   for( j = 0; j < (*consdata)->nvars; ++j )
   {
      mpq_clear((*consdata)->obj[j]);
      mpq_clear((*consdata)->lb[j]);
      mpq_clear((*consdata)->ub[j]);
      mpq_clear((*consdata)->lbloc[j]);
      mpq_clear((*consdata)->ubloc[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lockup, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lockdown, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ubloc, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lbloc, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ub, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lb, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->obj, (*consdata)->nvars);

   /* free consdata */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** print single constraint of exactlp constraint in CIP format to file stream */
static
void printSingleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   i                   /**< index of constraint in exactlp data structures */
   )
{
   char s[SCIP_MAXSTRLEN];
   int v;

   assert(i >= 0 && i < consdata->nconss);

   /* print left hand side for ranged rows */
   if( !isNegInfinity(conshdlrdata, consdata->lhs[i])
      && !isPosInfinity(conshdlrdata, consdata->rhs[i])
      && mpq_equal(consdata->lhs[i], consdata->rhs[i]) == 0 )
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, "%Qd <= ", consdata->lhs[i]);
      SCIPinfoMessage(scip, file, s);
   }

   /* print coefficients and variables */
   if( consdata->len[i] == 0 )
      SCIPinfoMessage(scip, file, "0 ");
   else
   {
      for( v = consdata->beg[i]; v < consdata->beg[i] + consdata->len[i]; v++ )
      {
         assert(consdata->ind[v] >= 0 && consdata->ind[v] < consdata->nvars);
         assert(SCIPvarGetProbindex(vars[consdata->ind[v]]) == consdata->ind[v]);

         gmp_snprintf(s, SCIP_MAXSTRLEN, "%+Qd<%s> ", consdata->val[v], SCIPvarGetName(vars[consdata->ind[v]]));
         SCIPinfoMessage(scip, file, s);
      }
   }

   /* print right hand side */
   if( mpq_equal(consdata->lhs[i], consdata->rhs[i]) != 0 )
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, "== %Qd\n", consdata->rhs[i]);
      SCIPinfoMessage(scip, file, s);
   }
   else if( !isPosInfinity(conshdlrdata, consdata->rhs[i]) )
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, "<= %Qd\n", consdata->rhs[i]);
      SCIPinfoMessage(scip, file, s);
   }
   else if( !isNegInfinity(conshdlrdata, consdata->lhs[i]) )
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, ">= %Qd\n", consdata->lhs[i]);
      SCIPinfoMessage(scip, file, s);
   }
   else
      SCIPinfoMessage(scip, file, " [free]\n");
}


/** print variable of exactlp constraint to file stream */
static
void printVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR*             var                 /**< problem variable to be printed */
   )
{
   char s[SCIP_MAXSTRLEN];
   int i;

   i = SCIPvarGetProbindex(var);

   /* name */
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "<%s>:", SCIPvarGetName(var));

   /* objective value */
   gmp_snprintf(s, SCIP_MAXSTRLEN, " obj=%Qd", consdata->obj[i]);
   SCIPinfoMessage(scip, file, s);

   /* bounds (global bounds for transformed variables, original bounds for original variables) */
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ", bounds=");

   if( isPosInfinity(conshdlrdata, consdata->lb[i]) )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "[+inf,");
   else if( isNegInfinity(conshdlrdata, consdata->lb[i]) )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "[-inf,");
   else
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, "[%Qd,", consdata->lb[i]);
      SCIPinfoMessage(scip, file, s);
   }

   if( isPosInfinity(conshdlrdata, consdata->ub[i]) )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "+inf]");
   else if( isNegInfinity(conshdlrdata, consdata->ub[i]) )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "-inf]");
   else
   {
      gmp_snprintf(s, SCIP_MAXSTRLEN, "%Qd]", consdata->ub[i]);
      SCIPinfoMessage(scip, file, s);
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
}

/** prints exactlp constraint in CIP format to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   int i;

   assert(conshdlrdata != NULL);
   assert(consdata != NULL);

   /* get problem variables */
   vars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);

   SCIPinfoMessage(scip, file, "\n");

   /* print objective sense */
   SCIPinfoMessage(scip, file, "  OBJECTIVE EXACT\n");
   SCIPinfoMessage(scip, file, "    Sense            : %s\n", consdata->objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");

   /* print variable specific information */
   if( consdata->nvars > 0 )
   {
      SCIPinfoMessage(scip, file, "  VARIABLES EXACT\n");
      for( i = 0; i < consdata->nvars; ++i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "    ");
         printVar(scip, conshdlrdata, consdata, file, vars[i]);
      }
   }

   /* print constraint and matrix specific information */
   if( consdata->nconss > 0 )
   {
      SCIPinfoMessage(scip, file, "  CONSTRAINTS EXACT\n");

      for( i = 0; i < consdata->nconss; ++i )
      {
         SCIPinfoMessage(scip, file, "    ");

         printSingleCons(scip, conshdlrdata, consdata, file, vars, i);
      }
   }
}

/** checks exactlp constraint for feasibility of given solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONS*            cons,               /**< exactlp constraint */
   SCIP_SOL*             sol,                /**< solution to be checked; or NULL */
   SCIP_SOLEX*           solex,              /**< exact solution to be checked; or NULL */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool*            feasible            /**< pointer to store whether the constraint is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   char s[SCIP_MAXSTRLEN];
   mpq_t activity;
   mpq_t solval;
   mpq_t prod;
   mpq_t violation;
   int v;
   int c;
   int i;

   assert(feasible != NULL);
   assert((sol != NULL && solex == NULL) || (sol == NULL && solex != NULL));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( sol != NULL )
   {
      SCIPdebugMessage("checking exactlp constraint <%s> for feasibility of solution %p\n",
         SCIPconsGetName(cons), (void*)sol);
   }
   else
   {
      SCIPdebugMessage("checking exactlp constraint <%s> for feasibility of exact solution %p\n",
         SCIPconsGetName(cons), (void*)solex);
   }

   *feasible = TRUE;

   /* increase age of constraint; age is reset to zero, if a violation was found */
   SCIP_CALL( SCIPincConsAge(scip, cons) );

   /* get problem variables */
   vars = SCIPgetVars(scip);
   assert(SCIPgetNVars(scip) == consdata->nvars);

   /* initializes gmp data */
   mpq_init(activity);
   mpq_init(solval);
   mpq_init(prod);
   mpq_init(violation);

   /* check bounds of variables for exact feasibility */
   for( v = 0; v < consdata->nvars && *feasible; ++v )
   {
      assert(SCIPvarGetProbindex(vars[v]) == v);

      if( sol != NULL )
      {
         /** @todo exiptodo: presolving extension
          *  - this only works if presolving is disabled (solval may already be an approximation since
          *    solution values of aggregated variables are calculated in FP arithmetic in SCIPgetSolVal())
          */
         mpq_set_d(solval, SCIPgetSolVal(scip, sol, vars[v]));
      }
      else
         SCIPsolexGetVal(solex, vars[v], solval);

      if( mpq_cmp(solval, consdata->ub[v]) > 0 || mpq_cmp(solval, consdata->lb[v]) < 0 )
      {
         *feasible = FALSE;

         if( printreason )
         {
            printVar(scip, conshdlrdata, consdata, NULL, vars[v]);
            if( mpq_cmp(solval, consdata->ub[v]) > 0 )
            {
               /* upper bound of current variable is violated */
               mpq_sub(violation, solval, consdata->ub[v]);
               assert(mpq_sgn(violation) > 0);
               gmp_snprintf(s, SCIP_MAXSTRLEN, "violation: upper bound is violated by %Qd\n", violation);
            }
            else
            {
               assert(mpq_cmp(solval, consdata->lb[v]) < 0);

               /* lower bound of current variable is violated */
               mpq_sub(violation, consdata->lb[v], solval);
               assert(mpq_sgn(violation) > 0);
               gmp_snprintf(s, SCIP_MAXSTRLEN, "violation: lower bound is violated by %Qd\n", violation);
            }
            SCIPinfoMessage(scip, NULL, s);
         }
      }
   }

   /* check each linear constraint of exactlp constraint for exact feasibility */
   for( c = 0; c < consdata->nconss && *feasible; ++c )
   {
      mpq_set_d(activity, 0.0);

      /* calculate activity of current constraint */
      for( i = consdata->beg[c]; i < consdata->beg[c] + consdata->len[c]; ++i )
      {
         assert(consdata->ind[i] >= 0 && consdata->ind[i] < SCIPgetNVars(scip));

         if( sol != NULL )
         {
            /** @todo exiptodo: presolving extension
             *  - this only works if presolving is disabled (solval may already be an approximation since
             *    solution values of aggregated variables are calculated in FP arithmetic in SCIPgetSolVal())
             */
            mpq_set_d(solval, SCIPgetSolVal(scip, sol, vars[consdata->ind[i]]));
         }
         else
            SCIPsolexGetVal(solex, vars[consdata->ind[i]], solval);

         mpq_mul(prod, consdata->val[i], solval);
         mpq_add(activity, activity, prod);
      }

      /* current constraint is violated */
      if( mpq_cmp(activity, consdata->rhs[c]) > 0 || mpq_cmp(activity, consdata->lhs[c]) < 0 )
      {
         *feasible = FALSE;

         if( printreason )
         {
            printSingleCons(scip, conshdlrdata, consdata, NULL, vars, c);
            if( mpq_cmp(activity, consdata->rhs[c]) > 0 )
            {
               /* rhs of current constraint is violated */
               mpq_sub(violation, activity, consdata->rhs[c]);
               assert(mpq_sgn(violation) > 0);
               gmp_snprintf(s, SCIP_MAXSTRLEN, "violation: right hand side of constraint is violated by %Qd\n", violation);
            }
            else
            {
               assert(mpq_cmp(activity, consdata->lhs[c]) < 0);

               /* lhs of current constraint is violated */
               mpq_sub(violation, consdata->lhs[c], activity);
               assert(mpq_sgn(violation) > 0);
               gmp_snprintf(s, SCIP_MAXSTRLEN, "violation: left hand side of constraint is violated by %Qd\n", violation);
            }
            SCIPinfoMessage(scip, NULL, s);
         }
      }
   }

   /* frees gmp data */
   mpq_clear(activity);
   mpq_clear(prod);
   mpq_clear(solval);
   mpq_clear(violation);

   return SCIP_OKAY;
}

/** applies all unprocessed bounds changes of variables, i.e., updates exact local bound information in constraint data,
 *  exact LP, and in project-and-shift data structure
 */
static
SCIP_RETCODE processBoundchgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(consdata != NULL);

   /* apply and delete unprocessed lower bound changes */
   assert(consdata->eventdatas != NULL);
   if( consdata->nbndchglb > 0 )
   {
      /* update lower bounds in exact LP */
      if( conshdlrdata->lpexconstructed )
      {
#ifdef BOUNDCHG_OUT
         SCIPdebugMessage("apply lower bound changes to LPEX (nbndchglb=%d, nbndchgub=%d):\n", consdata->nbndchglb,
            consdata->nbndchgub);
#endif
         SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchglb, consdata->bndchglbind,
               consdata->bndchglb, NULL) );
      }

#ifdef BOUNDCHG_OUT
      SCIPdebugMessage("apply lower bound changes to pseudoobjval:\n");
#endif
      /* update pseudo objective value and local lower bounds */
      for( i = 0; i < consdata->nbndchglb; ++i )
      {
         int probidx;

         probidx = consdata->bndchglbind[i];
         SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
               consdata->obj[probidx], consdata->bndchglb[i], consdata->ubloc[probidx]) );
         mpq_set(consdata->lbloc[probidx], consdata->bndchglb[i]);
      }

      /* delete unprocessed lower bound change information */
      for( i = 0; i < consdata->nvars; ++i )
         consdata->bndchglbpos[i] = -1;
      consdata->nbndchglb = 0;
   }

   /* apply and delete unprocessed upper bound changes */
   if( consdata->nbndchgub > 0 )
   {
      /* update upper bounds in exact LP */
      if( conshdlrdata->lpexconstructed )
      {
#ifdef BOUNDCHG_OUT
         SCIPdebugMessage("apply upper bound changes to LPEX (nbndchglb=%d, nbndchgub=%d):\n", consdata->nbndchglb,
            consdata->nbndchgub);
#endif
         SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchgub, consdata->bndchgubind,
               NULL, consdata->bndchgub) );
      }

#ifdef BOUNDCHG_OUT
      SCIPdebugMessage("apply upper bound changes to pseudo objval:\n");
#endif
      /* update pseudo objective value and local upper bounds */
      for( i = 0; i <  consdata->nbndchgub; ++i )
      {
         int probidx;

         probidx = consdata->bndchgubind[i];
         SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
               consdata->obj[probidx], consdata->lbloc[probidx], consdata->bndchgub[i]) );
         mpq_set(consdata->ubloc[probidx], consdata->bndchgub[i]);
      }

      /* delete unprocessed upper bound change information */
      for( i = 0; i < consdata->nvars; ++i )
         consdata->bndchgubpos[i] = -1;
      consdata->nbndchgub = 0;
   }

#ifndef NDEBUG
   /* check whether pseudo objective value was updated correctly */
   checkPseudoobjval(scip, conshdlrdata, consdata);
#endif

   return SCIP_OKAY;
}

/*
 * local methods for managing exact primal solutions
 */

/** returns the exact external value of the given exact internal objective value */
static
void getExternObjvalex(
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   const mpq_t           objintern,          /**< exact internal objective value */
   mpq_t                 objextern           /**< pointer to store exact external objective value */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   mpq_t objsense;

   assert(cons != NULL);

   /* get exactlp constraint handler */
   conshdlr = SCIPconsGetHdlr(cons);

   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_init(objsense);
   mpq_set_d(objsense, (SCIP_Real) consdata->objsense);

   if( isPosInfinity(conshdlrdata, objintern) )
      mpq_mul(objextern, objsense, *posInfinity(conshdlrdata));
   if( isNegInfinity(conshdlrdata, objintern) )
      mpq_mul(objextern, objsense, *negInfinity(conshdlrdata));
   else
   {
      /* objextern = objsense * objscale * (objintern + objoffset) */
      mpq_add(objextern, objintern, consdata->objoffset);
      mpq_mul(objextern, consdata->objscale, objextern);
      mpq_mul(objextern, objsense, objextern);
   }

   mpq_clear(objsense);
}

/** gets current pseudo objective value */
static
const mpq_t* getPseudoObjval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data to process the bound change for */
   )
{
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->pseudoobjvalinf >= 0);

   /* apply all unprocessed bounds changes of variables */
   processBoundchgs(scip, conshdlrdata, consdata);
   assert(consdata->nbndchgub == 0 && consdata->nbndchglb == 0);

   if( conshdlrdata->pseudoobjvalinf > 0 || SCIPignorePseudosol(scip) )
      return negInfinity(conshdlrdata);
   else
      return (const mpq_t*) (&conshdlrdata->pseudoobjval);
}

/** sets integral objective value flag, if all variables with non-zero objective values are integral and have
 *  integral objective value
 */
static
SCIP_RETCODE checkObjIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
#endif
   int v;

   assert(consdata != NULL);

   SCIPdebugMessage("check whether objective is always integral:\n");

   /* if we know already, that the objective value is integral, nothing has to be done */
   if( SCIPisObjIntegral(scip) )
      return SCIP_OKAY;

   /* if there exist unknown variables, we cannot conclude that the objective value is always integral */
   if( SCIPgetNActivePricers(scip) != 0 )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   gmp_snprintf(s, SCIP_MAXSTRLEN, "objoffset=%Qd (integral=%d)\n", consdata->objoffset,
      mpqIsIntegral(consdata->objoffset));
   SCIPdebugMessage(s);
#endif

   /* if the objective value offset is fractional, the value itself is possibly fractional */
   if( !mpqIsIntegral(consdata->objoffset) )
     return SCIP_OKAY;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );

   /* get transformed problem variables */
   origvars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);
   SCIP_CALL( SCIPgetTransformedVars(scip, consdata->nvars, origvars, vars) );

   /* scan through the variables */
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(SCIPvarIsOriginal(origvars[v]));
      assert(SCIPvarIsTransformed(vars[v]));

      /* check, if objective value is non-zero */
      if( mpq_sgn(consdata->obj[v]) != 0 )
      {
         /* if variable with non-zero objective value is continuous, the problem's objective value may be fractional */
         if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, " -> v=%d [type=%d (3=cont)]: objval=%Qd --> continuous var\n",
               v, SCIPvarGetType(vars[v]), consdata->obj[v]);
            SCIPdebugMessage(s);
#endif
            break;
         }

         /* if variable's objective value is fractional, the problem's objective value may also be fractional */
         if( !mpqIsIntegral(consdata->obj[v]) )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, " -> v=%d [type=%d (3=cont)]: objval=%Qd --> nonintegral val for int var\n",
               v, SCIPvarGetType(vars[v]), consdata->obj[v]);
            SCIPdebugMessage(s);
#endif
            break;
         }

#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, " -> v=%d [type=%d (3=cont)]: objval=%Qd\n",
               v, SCIPvarGetType(vars[v]), consdata->obj[v]);
            SCIPdebugMessage(s);
#endif
      }
   }

   /* objective value is integral, if the variable loop scanned all variables */
   if( v == consdata->nvars )
      SCIPsetObjIntegral(scip);

   SCIPdebugMessage("--> obj is %s always integral\n", SCIPisObjIntegral(scip) ? "" : "NOT");

#ifndef NDEBUG
   /* check for correct status of original variables */
   checkOrigVars(scip);
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** if possible, scales objective function such that it is integral with gcd = 1 */
static
SCIP_RETCODE scaleObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_VAR** vars;
   int nints;
   int v;

   assert(consdata != NULL);

   SCIPdebugMessage("scale objective to be integral with gcd = 1:\n");

   /* if there exist unknown variables, we cannot conclude that the objective value is always integral */
   if( SCIPgetNActivePricers(scip) != 0 )
      return SCIP_OKAY;

   /* get problem variables */
   vars = SCIPgetVars(scip);
   assert(SCIPgetNVars(scip) == consdata->nvars);

   nints = consdata->nvars - SCIPgetNContVars(scip);

   /* scan through the continuous variables */
   for( v = nints; v < consdata->nvars; ++v )
   {
      /* get objective value of variable; it it is non-zero, no scaling can be applied */
      if( mpq_sgn(consdata->obj[v]) != 0 )
         break;
   }

   /* only continue if all continuous variables have obj = 0 */
   if( v == consdata->nvars )
   {
      mpq_t* objvals;
      SCIP_Bool success;
      mpq_t intscalar;
      mpq_t one;
#ifdef SCIP_DEBUG
      char s[SCIP_MAXSTRLEN];
#endif

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nints) );

      /* get objective values of integer variables */
      for( v = 0; v < nints; ++v )
      {
         mpq_init(objvals[v]);
         mpq_set(objvals[v], consdata->obj[v]);
      }

      mpq_init(intscalar);
      mpq_init(one);
      mpq_set_d(intscalar, 1.0);
      mpq_set_d(one, 1.0);
      success = TRUE;

      /* calculate integral scalar */
      SCIP_CALL( SCIPmpqCalcIntegralScalar((const mpq_t*) objvals, nints, OBJSCALE_MAXFINALSCALE, intscalar, &success) );

#ifdef SCIP_DEBUG
      gmp_snprintf(s, SCIP_MAXSTRLEN, "integral objective scalar: success=%u, intscalar=%Qd\n", success, intscalar);
      SCIPdebugMessage(s);
#endif

      /* apply scaling */
      if( success && !mpq_equal(intscalar, one) )
      {
         /* calculate scaled objective values */
         for( v = 0; v < nints; ++v )
         {
	    mpq_mul(objvals[v], objvals[v], intscalar);
            assert(mpqIsIntegral(objvals[v]));

	    /* large integral values might not be FP representable which is required when we work with an FP relaxation */
	    if( SCIPuseFPRelaxation(scip) && !mpqIsReal(scip, objvals[v]) )
               break;
         }

         /* change the variables' objective values and adjust objscale */
         if( v == nints )
         {
            for( v = 0; v < nints; ++v )
            {
               assert(mpqIsIntegral(objvals[v]));
               assert(mpqIsReal(scip, objvals[v]));

#ifdef SCIP_DEBUG
               gmp_snprintf(s, SCIP_MAXSTRLEN, " -> v=%d: old val=%Qd --> \n", v, consdata->obj[v]);
               SCIPdebugMessage(s);
#endif

               SCIP_CALL( updateVar(conshdlrdata, consdata->obj[v], consdata->lbloc[v], consdata->ubloc[v],
                     objvals[v], consdata->lbloc[v], consdata->ubloc[v]) );
               mpq_set(consdata->obj[v], objvals[v]);

               SCIP_CALL( SCIPchgVarObj(scip, vars[v], mpqGetRealApprox(scip, objvals[v])) );

#ifdef SCIP_DEBUG
               gmp_snprintf(s, SCIP_MAXSTRLEN, "                           scaled val=%Qd\n", consdata->obj[v]);
               SCIPdebugMessage(s);
#endif
            }

            /* update objscale of transformed constraint;
             * in contrast to the objscale of the transformed problem this value here is exact. therefore, it is used
             * to compute primal and dual bounds (external values), see getPrimalbound() and getDualbound() in scip.c
             */
            mpq_div(consdata->objscale, consdata->objscale, intscalar);

            /* update objscale of transforemd problem;
             * note that this value is not reliable and should never be used to compute dual and primal bounds
             */
            SCIPsetTransObjscale(scip, mpqGetRealApprox(scip, consdata->objscale));

            SCIPsetObjIntegral(scip);

#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "integral objective scalar: objscale=%Qd\n", consdata->objscale);
            SCIPdebugMessage(s);
#endif

            /* update upperbound and cutoffbound in primal data structure:
             * there could already be solutions for the original problem or the upper and cutoff bound could be set
             * according to an objective limit (see handling in SCIPprimalUpdateObjoffset());
             * but currently we do not support this here, i.e., cutoff and upper boudn do not need to be updated
             */
            assert(SCIPgetNSolexs(scip) == 0 && SCIPgetNSols(scip) == 0 );
            assert(SCIPisInfinity(scip, SCIPgetObjlimit(scip)));
            assert(SCIPisInfinity(scip, SCIPgetCutoffbound(scip)));
            assert(SCIPisInfinity(scip, SCIPgetUpperbound(scip)));
         }
      }

      /* free temporary memory */
      mpq_clear(intscalar);
      for( v = 0; v < nints; ++v )
         mpq_clear(objvals[v]);
      SCIPfreeBufferArray(scip, &objvals);

   }

   return SCIP_OKAY;
}


/*
 * local methods for managing the LP relaxation
 */

/** creates LP rows corresponding to exactlp constraint */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR** rowvars;
   SCIP_Real* rowvals;
   SCIP_Real rowlhs;
   SCIP_Real rowrhs;

   char rowname[SCIP_MAXSTRLEN];
   int nrows;
   int v;
   int c;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nrows = consdataGetNRows(consdata);
   assert(SCIPgetNVars(scip) == consdata->nvars);
   assert(nrows == consdata->nconss);

   /* allocate memory for all rows */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, nrows) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvals, SCIPgetNVars(scip)) );

#ifdef FPRELAXLP_OUT
   printf("constructing FP %s of exactlp constraint<%s>:\n", SCIPuseFPRelaxation(scip) ? "relaxation" : "approximation",
      SCIPconsGetName(cons));
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif

   if( SCIPuseFPRelaxation(scip) )
   {
      /* for each row of the exactlp constraint, create a row with FP data that defines a relaxation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         assert(consdata->len[c] >= 0 && consdata->len[c] <= SCIPgetNVars(scip));

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealRelax(scip, consdata->lhs[c], GMP_RNDD);
         rowrhs = mpqGetRealRelax(scip, consdata->rhs[c], GMP_RNDU);

         /* calculate coefficients of all variables in the row */
         for( i = consdata->beg[c], v = 0; v < consdata->len[c]; ++i, ++v )
         {
            int probidx;

            probidx = consdata->ind[i];

            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) >= SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) <= SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) >= SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) <= SCIPvarGetUbLocal(vars[probidx]));
            assert(mpq_cmp(consdata->lb[probidx], consdata->ub[probidx]) <= 0);

            /* x_j >= 0 holds always ==>
             *   cons with rhs: underestimate(a_j) * x_j <= a_j * x_j <= rhs
             *   cons with lhs:  overestimate(a_j) * x_j >= a_j * x_j >= lhs
             */
            if( mpq_sgn(consdata->lb[probidx]) >= 0 )
            {
               if( !SCIPisInfinity(scip, rowrhs) )
               {
                  assert(SCIPisInfinity(scip, -rowlhs) || mpqIsReal(scip, consdata->val[i]));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -rowlhs) && SCIPisInfinity(scip, rowrhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU);
               }
            }
            /* x_j <= 0 always holds ==>
             *   cons with rhs:  overestimate(a_j) * x_j <= a_j * x_j <= rhs
             *   cons with lhs: underestimate(a_j) * x_j >= a_j * x_j >= lhs
             */
            else if( mpq_sgn(consdata->ub[probidx]) <= 0 )
            {
               if( !SCIPisInfinity(scip, rowrhs) )
               {
                  assert(SCIPisInfinity(scip, -rowlhs) || mpqIsReal(scip, consdata->val[i]));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -rowlhs) && SCIPisInfinity(scip, rowrhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD);
               }
            }
            /* x_j <= 0 and x_j >= 0 may hold but a_j is FP representable */
            else if( mpqIsReal(scip, consdata->val[i]) )
            {
               rowvals[v] = mpqGetRealApprox(scip, consdata->val[i]);
            }
            /* x_j <= 0 and x_j >= 0 may hold and a_j is not FP representable ==>
             *   split x_j into negative and positive part
             */
            else
            {
               /** @todo exiptodo: reading extension
                *  - split variable into positive and negative part and modify consdata structures in order to later be
                *    able to transform a basis to the original problem
                */
               SCIPerrorMessage("FP-relaxation not supported yet (requires splitting <%s>[%g,%g] into pos/neg-part) --> switch to FP-approximation via \"set misc usefprelax FALSE\"\n",
                  SCIPvarGetName(vars[probidx]), mpqGetRealApprox(scip, consdata->lb[probidx]),
                  mpqGetRealApprox(scip, consdata->ub[probidx]));

               /* free temporary memory */
               SCIPfreeBufferArray(scip, &rowvals);
               SCIPfreeBufferArray(scip, &rowvars);

               return SCIP_ERROR;
            }

            rowvars[v] = vars[probidx];
         }

         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_relax_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], consdata->len[c], rowvars, rowvals) );

#ifdef FPRELAXLP_OUT
         SCIPprintRow(scip, consdata->rows[c], NULL);
#endif
      }
   }
   else
   {
      /* for each row of the exactlp constraint, create a row with FP data that defines an approximation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         assert(consdata->len[c] >= 0 && consdata->len[c] <= SCIPgetNVars(scip));

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealApprox(scip, consdata->lhs[c]);
         rowrhs = mpqGetRealApprox(scip, consdata->rhs[c]);

         /* add all variables to the row */
         for( i = consdata->beg[c], v = 0; v < consdata->len[c]; ++i, ++v )
         {
            int probidx;
            probidx = consdata->ind[i];

            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) >= SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) <= SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) >= SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) <= SCIPvarGetUbLocal(vars[probidx]));
            assert(mpq_cmp(consdata->lb[probidx], consdata->ub[probidx]) <= 0);

            rowvals[v] = mpqGetRealApprox(scip, consdata->val[i]);
            rowvars[v] = vars[probidx];
         }

         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_approx_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], consdata->len[c], rowvars, rowvals) );

#ifdef FPRELAXLP_OUT
         SCIPprintRow(scip, consdata->rows[c], NULL);
#endif
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &rowvals);
   SCIPfreeBufferArray(scip, &rowvars);

   return SCIP_OKAY;
}

/** adds linear relaxation of exactlp constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   int nrows;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert( consdata->rows != NULL );

   nrows = consdataGetNRows(consdata);

   for( r = 0; r < nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         /** @todo exiptodo: is it better e.g. for basis transfer between LP and LPEX to force the cut to enter the LP? */
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], TRUE) );
      }
   }

   return SCIP_OKAY;
}

#if 0 /* uncomment if method is needed */
#ifdef BASISCHECK_OUT
static
SCIP_RETCODE printCurrentLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_ROW** rows;
   int i;
   int nrows;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIPdebugMessage("current inexact LP rows:\n");

   for( i = 0; i < nrows; ++i )
   {
      SCIP_CALL( SCIPprintRow(scip, rows[i], NULL) );
   }

   return SCIP_OKAY;
}
#endif
#endif

/** loads LP state from exact LP into inexact LP solver */
static
SCIP_RETCODE loadLPEXState(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_Bool*            success             /**< pointer to store whether LP state was excepted by inexact LP solver */
   )
{
   SCIP_LPISTATE* lpistate;
   int ncolsex;
   int nrowsex;

   assert(success != NULL);

   *success = TRUE;

   SCIP_CALL( SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncolsex) );
   SCIP_CALL( SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrowsex) );

   if( ncolsex == SCIPgetNLPCols(scip) && nrowsex == SCIPgetNLPRows(scip) )
   {
      /* stores LP state (like basis information) of exact LP solver into LP state object */
      SCIP_CALL( SCIPlpiexGetState(conshdlrdata->lpiex, SCIPblkmem(scip), &lpistate) );

      /* loads LP state (like basis information) into exact LP solver */
      SCIP_CALL( SCIPsetLPState(scip, lpistate) );

#ifndef NDEBUG
      /* check whether basis of exact and inexact LP solvers match it was loaded from inexact solver to exact one */
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif

      /* frees LP state */
      SCIP_CALL( SCIPlpiexFreeState(conshdlrdata->lpiex, SCIPblkmem(scip), &lpistate) );

      assert(*success);
   }
   else
   {
      *success = FALSE;

#ifdef BASISCHECK_OUT
      SCIPdebugMessage("loading LP basis into LPEX failed because dimensions do not match (LP: %d x %d, LPEX: %d x %d)\n",
         SCIPgetNLPRows(scip), SCIPgetNLPCols(scip), nrowsex, ncolsex);
#endif
   }

   return SCIP_OKAY;
}


/*
 * local methods for managing the exact LP relaxation
 */

/** constructs the exact LP of the current node, but does not load the LP state and warmstart information  */
static
SCIP_RETCODE constructCurrentLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(consdata != NULL);

   /* exact LP is constructed when it is needed for the first time */
   if( !conshdlrdata->lpexconstructed )
   {
      SCIP_VAR** vars;
      char** colnames;

      SCIPdebugMessage("constructing initial exact LP\n");

      /* allocate and initialize temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, consdata->nvars) );

      /* get names of problem variables */
      vars = SCIPgetVars(scip);
      assert(SCIPgetNVars(scip) == consdata->nvars);
      for( i = 0; i < consdata->nvars; ++i )
      {
         /* allocate and initialize temporary memory */
         colnames[i] = (char*) (SCIPvarGetName(vars[i]));

         assert(SCIPvarGetProbindex(vars[i]) == i);
         assert(mpqGetRealApprox(scip, consdata->obj[i]) == SCIPvarGetObj(vars[i]));
         assert(SCIPisLE(scip, mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD), SCIPvarGetLbGlobal(vars[i])));
         assert(SCIPisGE(scip, mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU), SCIPvarGetUbGlobal(vars[i])));
         assert(SCIPisLE(scip, mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD), SCIPvarGetLbLocal(vars[i])));
         assert(SCIPisGE(scip, mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU), SCIPvarGetUbLocal(vars[i])));
      }

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(conshdlrdata->lpiex, consdata->nvars, consdata->obj, consdata->lbloc, consdata->ubloc,
            colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(conshdlrdata->lpiex, consdata->nconss, (const mpq_t*) consdata->lhs,
            (const mpq_t*) consdata->rhs, NULL, consdata->nnonz, consdata->beg, consdata->len, consdata->ind,
            consdata->val) );

#if 0 /* activate this for detailed debugging. writes LPEX to a file and stops afterwards */
      SCIP_CALL( SCIPlpiexWriteLP(conshdlrdata->lpiex, "testset/debug.lp") );
      SCIPABORT();
#endif
      conshdlrdata->lpexconstructed = TRUE;

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &colnames);
   }

   /* apply all unprocessed bounds changes of variables */
   processBoundchgs(scip, conshdlrdata, consdata);

   return SCIP_OKAY;
}

/** loads LP state from inexact LP into exact LP solver */
static
SCIP_RETCODE loadLPState(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   )
{
   SCIP_LPISTATE* lpistate;
   int ncolsex;
   int nrowsex;

   SCIP_CALL( SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncolsex) );
   SCIP_CALL( SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrowsex) );

   if( ncolsex == SCIPgetNLPCols(scip) && nrowsex == SCIPgetNLPRows(scip) )
   {
      /** @todo exiptodo:
       *  - for infeasible LPs with unbounded dual ray, check whether this always gives the basis corresponding to the
       *    node from which the unbounded ray was found
       */
      /* stores LP state (like basis information) of inexact LP solver into LP state object */
      SCIP_CALL( SCIPgetLPState(scip, &lpistate) );

      /* loads LP state (like basis information) into exact LP solver */
      SCIP_CALL( SCIPlpiexSetState(conshdlrdata->lpiex, SCIPblkmem(scip), lpistate) );

#ifndef NDEBUG
      /* check whether basis of exact and inexact LP solvers match it was loaded from inexact solver to exact one */
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif

      /* free LP state */
      SCIP_CALL( SCIPfreeLPState(scip, &lpistate) );
   }
#ifdef BASISCHECK_OUT
   else
   {
      SCIPdebugMessage("loading LP basis into LPEX failed because dimensions don't match (LP: %d x %d, LPEX: %d x %d)\n",
         SCIPgetNLPRows(scip), SCIPgetNLPCols(scip), nrowsex, ncolsex);
   }
#endif
   return SCIP_OKAY;
}

/** solves the exact LP with the given algorithm and evaluates return status */
static
SCIP_RETCODE solveLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool             usetimelimit,       /**< should timelimit be set for exact LP solver? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;

   int iterations;
   int ncols;
   int nrows;

   assert(scip != NULL);
   assert(lperror != NULL);

   SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncols);
   SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrows);

   *lperror = FALSE;

#if defined(LPEX_OUT) || defined(TESTDBVERIFY)
   SCIP_CALL( SCIPlpiexSetIntpar(conshdlrdata->lpiex, SCIP_LPPAR_LPINFO, TRUE) );
#endif

#ifdef USEOBJLIM
   {
      /* set objective limit for exact LP solver */
      mpq_t uobjlim;

      mpq_init(uobjlim);

      if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
         SCIPlpiexPosInfinity(conshdlrdata->lpiex, &uobjlim);
      else
         mpq_set_d(uobjlim, SCIPgetCutoffbound(scip));

      if( mpq_cmp(uobjlim, conshdlrdata->lpiexuobjlim) < 0 )
      {
         SCIP_CALL( SCIPlpiexSetRealpar(conshdlrdata->lpiex, SCIP_LPPAR_UOBJLIM, uobjlim) );
         mpq_set(conshdlrdata->lpiexuobjlim, uobjlim);
         SCIPdebugMessage("---> set uobjlim!\n");
      }
      mpq_clear(uobjlim);
   }
#endif

#ifdef USETIMELIM
   {
      if( usetimelimit )
      {
         mpq_t lptimelimit;
         SCIP_Real timelimit;

         mpq_init(lptimelimit);

         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

         mpq_set_d(lptimelimit, timelimit - SCIPgetSolvingTime(scip));
         if( mpq_sgn(lptimelimit) > 0 )
         {
            SCIP_CALL( SCIPlpiexSetRealpar(conshdlrdata->lpiex, SCIP_LPPAR_LPTILIM, lptimelimit) );
         }
         mpq_clear(lptimelimit);
      }
   }
#endif

   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      SCIPdebugMessage("   solving current primal exact LP (%d cols, %d rows):\n", ncols, nrows);

      /* call primal simplex */
      retcode = SCIPlpiexSolvePrimal(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("      (node %"SCIP_LONGINT_FORMAT") primal simplex solving error in current exact LP\n",
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }

      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("      solved primal exact LP in %d iterations\n", iterations);
      break;

   case SCIP_LPALGO_DUALSIMPLEX:
      SCIPdebugMessage("   solving current dual exact LP (%d cols, %d rows)\n", ncols, nrows);

      /* call dual simplex */
      retcode = SCIPlpiexSolveDual(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("      (node %"SCIP_LONGINT_FORMAT") dual simplex solving error in current exact LP\n",
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }

      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("      solved dual exact LP in %d iterations\n", iterations);
      break;

   default:
      SCIPerrorMessage("invalid exact LP algorithm\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** checks whether primal solution of exact LP solver satisfies all integrality restrictions */
static
SCIP_RETCODE checkIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_RESULT*          result              /**< pointer to store the result of the lp enforcement call */
   )
{
   SCIP_VAR** vars;
   mpq_t* primsol;
   mpq_t lpobjval;
   int nvars;
   int nbin;
   int nint;
   int v;
   int ncols;
   SCIP_Bool inrange;
   SCIP_Bool integral;
   int branchvar;

#ifdef SCIP_DEBUG
   SCIP_Real oldlb;
   oldlb = SCIPgetLocalLowerbound(scip);
#endif

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(result != NULL);

   SCIPdebugMessage("checking integrality of exact LP solution:\n");

   SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncols);
   assert(SCIPgetNVars(scip) == ncols);

   /* allocate and initialize temporary arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &primsol, ncols) );
   for( v = 0; v < ncols; ++v )
   {
      mpq_init(primsol[v]);
   }
   mpq_init(lpobjval);

   /* gets primal solution vector of exact LP */
   SCIP_CALL( SCIPlpiexGetSol(conshdlrdata->lpiex, &lpobjval, primsol, NULL, NULL, NULL) );

   /* try to improve current local lower bound by using the result of the exact LP solver */
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, lpobjval, GMP_RNDD)) );

#ifdef SCIP_DEBUG
   if( oldlb < SCIPgetLocalLowerbound(scip) )
   {
      SCIPdebugMessage("   local lower bound improved: %.50f --> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
   }
   else
   {
      SCIPdebugMessage("   local lower bound did not improve: %.50f -/-> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
   }
#endif

   integral = TRUE;
   inrange = TRUE;

   /* check for infeasible node by bounding */
   if( SCIPgetLocalLowerbound(scip) >= SCIPgetCutoffbound(scip) )
   {
      SCIPdebugMessage("   --> cut off node (lower=%g >= upper=%g)\n", SCIPgetLocalLowerbound(scip),
         SCIPgetCutoffbound(scip));

      *result = SCIP_CUTOFF;
   }
   else
   {
      /* get all problem variables and integer region in vars array */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbin, &nint, NULL, NULL) );

      /* check whether primal solution satisfies all integrality restrictions */
      branchvar = -1;
      for( v = 0; v < nbin + nint && integral; ++v )
      {
         assert(SCIPvarGetProbindex(vars[v]) == v);
         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER );

         if( !mpqIsIntegral(primsol[v]) )
         {
#ifdef SCIP_DEBUG
            char s[SCIP_MAXSTRLEN];
            gmp_snprintf(s, SCIP_MAXSTRLEN, "   exact LP value<%Qd> of intvar<%s> is not integral", primsol[v],
               SCIPvarGetName(vars[v]));
            SCIPdebugMessage(s);
#endif
            integral = FALSE;
            branchvar = v;
         }
      }

      /* store exact LP solution as feasible MIP solution and cut off current node */
      if( integral )
      {
         SCIP_SOL* sol;
         SCIP_SOLEX* solex;
         mpq_t tmp;
         SCIP_Bool stored;
         SCIP_Bool fpvalue;

         SCIPdebugMessage("   ---> exact LP solution is integral\n");

         assert(branchvar == -1);

         mpq_init(tmp);

         /* create exact and approximate primal solution */
         SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
         SCIP_CALL( SCIPsolexCreate(&solex, SCIPblkmem(scip), NULL) );

         fpvalue = TRUE;

         /* check whether it is likely that the lpobjval is treated as infinite in scip, i.e., is out of range in scip */
         if( SCIPisInfinity(scip, REALABS(mpqGetRealApprox(scip, lpobjval))) )
            inrange = FALSE;

         for( v = 0; v < nvars; ++v )
         {
            SCIP_Real scipsolval;

            scipsolval = mpqGetRealApprox(scip, primsol[v]);

            /* check whether we can store the primal solution without FP errors */
            mpq_set_d(tmp, scipsolval);
            if( SCIPisInfinity(scip, REALABS(mpqGetRealApprox(scip, primsol[v]))) )
               inrange = FALSE;
            if( !mpq_equal(primsol[v], tmp) )
               fpvalue = FALSE;

            /* set value of variabel in exact and approximate solution */
            /** @todo exiptodo:
             *  - in SCIPsetSolVal() the objval is computed using interval arithmetic which leads to a safe objval in
             *    case the varsolvals are FP representable. We needed that in an older version where we did not yet store
             *    the exact solution.
             *  - now, we could compute the objval unsafe in SCIPsetSolVal() (below it's set to a safe value anyway)
             *    and also delete the fpvalue flag here.
             */
            SCIP_CALL( SCIPsetSolVal(scip, sol, vars[v], scipsolval) );
            /** @todo exiptodo: using scip->set is not SCIP conform
             *  - think about a more general way to get information from cons_exactlp.c.
             *  - maybe, implement all methods that operate/use exact problem data (like solex methods) in cons_exactlp.c
             */
            SCIP_CALL( SCIPsolexSetVal(solex, scip->set, vars[v], consdata->obj[v], consdata->lb[v], primsol[v]) );
         }
         mpq_set_d(tmp, SCIPgetSolTransObj(scip, sol));
         assert(!inrange || !fpvalue || !SCIPuseFPRelaxation(scip) || mpq_cmp(lpobjval, tmp) <= 0);

         /* add approximate primal solution */
         if( !inrange )
            SCIPfreeSol(scip, &sol);
         else
         {
            /** @todo exiptodo:
             *  - setting the upperbound in SCIP via the FP-solution using SCIPsetSolTransObj() is more a workaround
             *  - maybe it is better to set the upperbound directly via the exact solution. this way the code also works
             *    when we do not store a FP-solution for every exact solution.
             *  - remove SCIPsetSolTransObj() if it is not needed anymore
             */
            /* store safe objective value for approximate primal solution */
            SCIP_CALL( SCIPsetSolTransObj(scip, sol, mpqGetRealRelax(scip, lpobjval, GMP_RNDU)) );
            SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         }

         /** @todo exiptodo: using scip->set is not SCIP conform
          *  - think about a more general way to get information from cons_exactlp.c.
          *  - maybe, implement all methods that operate/use exact problem data (like solex methods) in cons_exactlp.c
          */
         /* add exact primal solution */
         SCIP_CALL( SCIPprimalexAddSolFree(conshdlrdata->primal, SCIPblkmem(scip), scip->set, scip->transprob,
               &solex, &stored) );

         *result = SCIP_CUTOFF;

         mpq_clear(tmp);
      }
      else
      {
         SCIP_Bool solvelpagain;

         SCIPdebugMessage("   exact LP solution is not integral\n");

         /* if the inexact LP was claimed to be infeasible, i.e., a pseudo solution is enforced here, but the exact LP
          * has been proved to be feasible, try to solve the inexact LP again with starting basis equal to optimal
          * basis of the exact LP; in order to avoid an infinite loop, this is only done once at a node
          */
         solvelpagain = FALSE;
         if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE
            && conshdlrdata->lastenfopsnode != SCIPgetCurrentNode(scip) )
         {
            SCIP_CALL( loadLPEXState(scip, conshdlrdata, &solvelpagain) );

            if( solvelpagain )
            {
               *result = SCIP_SOLVELP;
#ifdef SCIP_DEBUG
               SCIPdebugMessage("   ---> LP was infeasible but LPEX not: solve LP once again with basis of LPEX\n");
#endif
            }
         }

         /* branch on the first fractional variable in the exact LP solution */
         if( !solvelpagain )
         {
            SCIP_NODE* node;
            SCIP_Real solval;
            SCIP_Real downub;
            SCIP_Real uplb;

            /** @todo exiptodo: branching extension
             *  - implement more sophisticated rule for branching on the exact LP solution
             */

            assert(branchvar >= 0 && branchvar < nvars);

            /* create left child: add x_i <= floor(x_i^*) */
            solval = mpqGetRealRelax(scip, primsol[branchvar], GMP_RNDD);
            downub = floor(solval);
            SCIP_CALL( SCIPcreateChild(scip, &node, SCIPcalcNodeselPriority(scip, vars[branchvar],
                     SCIP_BRANCHDIR_DOWNWARDS, downub), SCIPcalcChildEstimate(scip, vars[branchvar], downub)) );
            SCIP_CALL( SCIPchgVarUbNode(scip, node, vars[branchvar], downub) );

            /* create right child: add x_i >= ceil(x_i^*) */
            solval = mpqGetRealRelax(scip, primsol[branchvar], GMP_RNDU);
            uplb = ceil(solval);
            SCIP_CALL( SCIPcreateChild(scip, &node, SCIPcalcNodeselPriority(scip, vars[branchvar],
                     SCIP_BRANCHDIR_UPWARDS, uplb), SCIPcalcChildEstimate(scip, vars[branchvar], uplb)) );
            SCIP_CALL( SCIPchgVarLbNode(scip, node, vars[branchvar], uplb) );

            assert(uplb == downub + 1.0);

#ifndef NDEBUG
            {
               mpq_t resultuplb;
               mpq_t resultdownub;

               mpq_init(resultuplb);
               mpq_init(resultdownub);

               mpq_set_d(resultuplb, uplb);
               mpq_set_d(resultdownub, downub);

               assert(mpq_cmp(resultdownub, primsol[branchvar]) < 0);
               assert(mpq_cmp(resultuplb, primsol[branchvar]) > 0);

               mpq_clear(resultuplb);
               mpq_clear(resultdownub);
            }
#endif
            *result = SCIP_BRANCHED;

#ifdef SCIP_DEBUG
            SCIPdebugMessage("   ---> branched on exact LP solution\n");
#endif
         }
      }
   }

   /* free temporary arrays */
   for( v = 0; v < ncols; ++v )
   {
      mpq_clear(primsol[v]);
   }
   mpq_clear(lpobjval);
   SCIPfreeBufferArray(scip, &primsol);

   if( integral && !inrange )
   {
      SCIPerrorMessage("storing optimal solutions of subproblems that is out of FP-range is not supported yet\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** evaluates the result of the exact LP */
static
SCIP_RETCODE evaluateLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_RESULT*          result              /**< pointer to store the result of the lp enforcement call */
   )
{
   /** @todo exiptodo:
    *  - in a similar function for LP in lp.c the case SCIPlpiIsPrimalUnbounded() is not explicitely handled. why? if case
    *    should be added, include it here as well.
    */
   /* evaluate solution status */
   if( SCIPlpiexIsOptimal(conshdlrdata->lpiex) )
   {
      SCIPdebugMessage("   exact LP solved to optimality\n");

#ifndef NDEBUG
      {
         SCIP_Bool primalfeasible;
         SCIP_Bool dualfeasible;

         SCIP_CALL( SCIPlpiexGetSolFeasibility(conshdlrdata->lpiex, &primalfeasible, &dualfeasible) );
         assert(primalfeasible);
         assert(dualfeasible);
      }
#endif
      /* check whether exact LP solution is feasible for the MIP; if it is feasible the solution is stored
       * and the current node is cut off otherwise a branching is created
       */
      SCIP_CALL( checkIntegrality(scip, cons, conshdlrdata, consdata, result) );
      assert(*result == SCIP_CUTOFF || *result == SCIP_BRANCHED || *result == SCIP_SOLVELP );
   }
   else if( SCIPlpiexIsObjlimExc(conshdlrdata->lpiex) )
   {
#ifdef USEOBJLIM
      SCIPdebugMessage("   exact LP exceeds upper objective limit\n");
      *result = SCIP_CUTOFF;
#else
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
      return SCIP_ERROR;
#endif
   }
   else if( SCIPlpiexIsPrimalInfeasible(conshdlrdata->lpiex) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n");

      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiexExistsPrimalRay(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiexIsIterlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiexIsTimelimExc(conshdlrdata->lpiex) )
   {
#ifdef USETIMELIM
      SCIPdebugMessage("   exact LP exceeds time limit\n");
      *result = SCIP_INFEASIBLE;
#else
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
      return SCIP_ERROR;
#endif
   }
   else
   {
      SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT") error or unknown return status in current exact LP (internal status: %d)\n",
         SCIPgetNNodes(scip), SCIPlpiexGetInternalStatus(conshdlrdata->lpiex));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}


/*
 * local methods for project-and-shift dual bounding method
 */

/** subroutine of constructPSdata(); chooses which columns of the matrix are designated as set S, used for projections */
static
SCIP_RETCODE psChooseS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   int i;
   int nconss;
   int nvars;
   int nextendedconss;
   int nrows;

   /* solution information for exact root LP */
   mpq_t* rootactivity;
   mpq_t* rootprimal;
   SCIP_Bool lperror;

   SCIP_COL** cols;
   SCIP_ROW** rows;

   nconss = consdata->nconss;
   nvars = consdata->nvars;
   nextendedconss = conshdlrdata->nextendedconss;

   rootactivity = NULL;
   rootprimal = NULL;


   /* build includedcons vector based on psdualcolselection, this determines the matrix D */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->includedcons, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      conshdlrdata->includedcons[i] = 0;
   if( conshdlrdata->psdualcolselection == 'n' || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      /* determine which dual variables to included in the problem
       * (ones with finite dual objective coef. in [lhs',-rhs',lb',-ub'])
       */
      for( i = 0; i < nconss; i++ )
      {
         if( !isNegInfinity(conshdlrdata, consdata->lhs[i]) )
            conshdlrdata->includedcons[i] = 1;
         if( !isPosInfinity(conshdlrdata, consdata->rhs[i]) )
            conshdlrdata->includedcons[nconss + i] = 1;
      }
      for( i = 0; i < nvars; i++ )
      {
         if( !isNegInfinity(conshdlrdata, consdata->lb[i]) )
            conshdlrdata->includedcons[2*nconss + i] = 1;
         if( !isPosInfinity(conshdlrdata, consdata->ub[i]) )
            conshdlrdata->includedcons[2*nconss + nvars + i] = 1;
      }
   }
   else if( conshdlrdata->psdualcolselection == 'a' )
   {
      /* determone which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      /* constructs exact LP of current node */
      constructCurrentLPEX(scip, conshdlrdata, consdata);

      /* load LP state from inexact LP into exact LP solver */
      SCIP_CALL( loadLPState(scip, conshdlrdata) );

      /* solve exact LP */
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, FALSE, &lperror) );
      if( lperror )
      {
         SCIPerrorMessage("Error solving root node LP in project-and-shift method.\n");
         return SCIP_ERROR;
      }

      /* allocate space for root lp return data */
      SCIP_CALL( SCIPallocBufferArray(scip, &rootprimal, nvars) );
      for( i = 0; i < nvars ; i++ )
         mpq_init(rootprimal[i]);

      SCIP_CALL( SCIPallocBufferArray(scip, &rootactivity, nconss) );
      for( i = 0; i < nconss ; i++ )
         mpq_init(rootactivity[i]);

      /* get the primal solution and activity */
      SCIP_CALL( SCIPlpiexGetSol(conshdlrdata->lpiex, NULL, rootprimal, NULL, rootactivity, NULL) );

      /* determine which dual variables to include in the problem
       * (primal constraints active at optimal solution found at root node)
       */
      for( i = 0; i < nconss; i++ )
      {
         if( mpq_equal(rootactivity[i], consdata->lhs[i]) )
            conshdlrdata->includedcons[i] = 1;
         if( mpq_equal(rootactivity[i], consdata->rhs[i]) )
            conshdlrdata->includedcons[nconss + i] = 1;
      }
      for( i = 0; i < nvars; i++ )
      {
         if( mpq_equal(rootprimal[i], consdata->lb[i]) )
            conshdlrdata->includedcons[2*nconss + i] = 1;
         if( mpq_equal(rootprimal[i], consdata->ub[i]) )
            conshdlrdata->includedcons[2*nconss + nvars + i] = 1;
      }

      /* free locally used memory */
      for( i = 0; i < nconss ; i++ )
         mpq_clear(rootactivity[i]);
      SCIPfreeBufferArray(scip, &rootactivity );
      for( i = 0; i < nvars ; i++ )
         mpq_clear(rootprimal[i]);
      SCIPfreeBufferArray(scip, &rootprimal);
   }
   else if( conshdlrdata->psdualcolselection == 'A' )
   {
      /* determine which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      /* if the active columns have already been determiend, assign them */
      if( conshdlrdata->rootactivitycon )
      {
         for( i = 0; i < nextendedconss; i++ )
         {
            if( conshdlrdata->rootactivity[i] )
               conshdlrdata->includedcons[i] = 1;
            else
               conshdlrdata->includedcons[i] = 0;
         }
      }
      else
      {
         SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
         assert(nrows == nconss);
         for( i = 0; i < nconss; i++ )
         {
            if( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, rows[i]), SCIProwGetLhs(rows[i])) )
               conshdlrdata->includedcons[i] = 1;
            if( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, rows[i]), SCIProwGetRhs(rows[i])) )
               conshdlrdata->includedcons[nconss + i] = 1;
         }
         cols = SCIPgetLPCols(scip);
         for( i = 0; i < nvars; i++ )
         {
            if( SCIPisFeasEQ(scip, SCIPcolGetPrimsol(cols[i]), SCIPcolGetLb(cols[i])) )
               conshdlrdata->includedcons[2*nconss + i] = 1;
            if( SCIPisFeasEQ(scip, SCIPcolGetPrimsol(cols[i]), SCIPcolGetUb(cols[i])) )
               conshdlrdata->includedcons[2*nconss + nvars + i] = 1;
         }
      }
   }
   else
   {
      SCIPerrorMessage("Invald value for parameter psdualcolselection\n");
   }
   return SCIP_OKAY;
}

/** subroutine of constructPSdata(); computes the LU factorization used by the project-and-shift method */
static
SCIP_RETCODE psFactorizeD(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{

   int i;
   int j;
   int rval;
   int pos;
   int nconss;
   int nvars;
   int nextendedconss;
   int nnonz;

   /* sparse representation of the matrix used for the LU factorization */
   int* projbeg;
   int* projlen;
   int* projind;
   mpq_t* projval;

   nconss = consdata->nconss;
   nvars = consdata->nvars;
   nextendedconss = conshdlrdata->nextendedconss;
   nnonz = consdata->nnonz;

   projbeg = NULL;
   projlen = NULL;
   projind = NULL;
   projval = NULL;

   /* allocate memory for the projection factorization */
   SCIP_CALL( SCIPallocBufferArray(scip, &projbeg, nextendedconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &projlen, nextendedconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &projind, 2*nnonz + 2*nvars) );
   for( i = 0; i < 2*nnonz + 2*nvars; i++)
      projind[i]=0;
   SCIP_CALL( SCIPallocBufferArray(scip, &projval, 2*nnonz + 2*nvars) );
   for( i = 0; i < 2*nnonz + 2*nvars; i++)
      mpq_init(projval[i]);

   /* allocate memory for the basis mapping */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->psbasis, nextendedconss) );

   /* use includedcons to construct psbasis, a description/mapping for D it has length npsbasis and psbasis[i] tells what
    * column (out of the original nextendecons) the ith column in D is
    */
   pos = 0;
   for( i = 0; i < nextendedconss; i++ )
   {
      if( conshdlrdata->includedcons[i] )
      {
         conshdlrdata->psbasis[pos] = i;
         pos++;
      }
   }
   conshdlrdata->npsbasis = pos;

   /* build the sparse representation of D that will be passed to the RECTLU code for factorization */
   pos = 0;
   for( i = 0; i < nextendedconss; i++ )
   {
      /* A part (lhs constraints) */
      if(i < nconss)
      {
         projlen[i] = consdata->len[i];
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = consdata->ind[ consdata->beg[i] + j];
            mpq_set( projval[ projbeg[i] + j], consdata->val[ consdata->beg[i] + j ] );
         }
         pos += consdata->len[i];
      }
      /* -A part (rhs constraints) */
      else if(i < 2 * nconss)
      {
         projlen[i] = consdata->len[i - nconss];
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = consdata->ind[ consdata->beg[i - nconss] + j ];
            mpq_neg( projval[ projbeg[i] + j], consdata->val[ consdata->beg[i - nconss] + j ] );
         }
         pos += consdata->len[i - nconss];
      }
      /* I part (lb constraints) */
      else if (i < 2*nconss + nvars)
      {
         projbeg[i] = pos;
         projlen[i] = 1;
         projind[pos] = i - 2*nconss;
         mpq_set_si(projval[pos], 1, 1);
         pos ++;
      }
      /* -I part (ub constraints) */
      else
      {
         projbeg[i] = pos;
         projlen[i] = 1;
         projind[pos] = i - (2*nconss + nvars);
         mpq_set_si(projval[pos], -1, 1);
         pos ++;
      }
   }

#ifdef PS_OUT
   printf("factoring matrix: nvars=%d, npsbasis=%d\n", nvars, conshdlrdata->npsbasis);
   for( i = 0; i < nextendedconss; i++ )
      printf("   j=%d:\t projbeg=<%d>,\t projlen=<%d>\n", i, projbeg[i], projlen[i]);

   for( i = 0; i < 2*nnonz + 2*nvars; i++ )
   {
      printf("   i=%d:\t projind=<%d>,\t projval=<", i, projind[i]);
      mpq_out_str(stdout, 10, projval[i]);
      printf(">\n");
   }
#endif

   /* factorize projection matrix D
    * - psbasis stores a mapping to tell us what D is, i.e. the dual columns corresponding to dual valuse that have a
    *   strictly positive value in the relative interior point
    * - D is equal to a subset of [A',-A',I,-I] and is given to the factor code in sparse column representation
    */
   rval = RECTLUbuildFactorization(&conshdlrdata->rectfactor, nvars, conshdlrdata->npsbasis,
      conshdlrdata->psbasis, projval, projind, projbeg, projlen);

   /* if rval != 0 then RECTLUbuildFactorization has failed. In this case the project-and-shift method will not work and
    * we will return failure
    */
   if( rval )
   {
      conshdlrdata->psdatafail = TRUE;
      SCIPdebugMessage("factorization of matrix for project-and-shift method failed.\n");
   }

#ifdef PS_OUT
   printf("   matrix factorization complete: %s\n", rval ? "failed" : "correct termination");
#endif

   for( i = 0; i < 2*nnonz + 2*nvars; i++ )
      mpq_clear(projval[i]);
   SCIPfreeBufferArray(scip, &projval);
   SCIPfreeBufferArray(scip, &projind);
   SCIPfreeBufferArray(scip, &projlen);
   SCIPfreeBufferArray(scip, &projbeg);

   return SCIP_OKAY;
}

/** prints error related to the current lpiex status, if there is one */
static
SCIP_RETCODE printlpiexerr(
   SCIP_LPIEX*           lpiex              /**< lpiex interface */
   )
{
   if( SCIPlpiexIsOptimal(lpiex) )
   {
      return SCIP_OKAY;
   }
   else if( SCIPlpiexIsObjlimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
   }
   else if( SCIPlpiexIsPrimalInfeasible(lpiex) )
   {
      SCIPerrorMessage(" Exact LP infeas.\n");
   }
   else if( SCIPlpiexExistsPrimalRay(lpiex) )
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
   }
   else if( SCIPlpiexIsIterlimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
   }
   else if( SCIPlpiexIsTimelimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
   }
   else
   {
      SCIPerrorMessage("lpiex not solved, or other error\n");
   }
   return SCIP_OKAY;
}

/** using the given rational data, an approximate LP is constructed using the LPI interface, solved and then its LP state
 *  is assigned to lpiex for the purpose of warm starting the next lpiex solve data for exact problem
 */
static
SCIP_RETCODE warmStartlpiex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPIEX*           lpiex,              /**< lpiex interface */
   int                   nvars,              /**< nvars for exact problem */
   int                   nconss,             /**< nconss for exact problem */
   int                   nnonz,              /**< nnonz for exact problem */
   mpq_t*                obj,                /**< obj for exact problem */
   mpq_t*                lb,                 /**< lb for exact problem */
   mpq_t*                ub,                 /**< ub for exact problem */
   mpq_t*                lhs,                /**< lhs for exact problem */
   mpq_t*                rhs,                /**< rhs for exact problem */
   char**                colnames,           /**< colnames for exact problem */
   int*                  ind,                /**< ind for exact problem */
   int*                  beg,                /**< beg for exact problem */
   mpq_t*                val                 /**< val for exact problem */
   )
{
   SCIP_LPI* lpi;
   SCIP_LPISTATE* lpistate;
   SCIP_Real* obj_real;
   SCIP_Real* lb_real;
   SCIP_Real* ub_real;
   SCIP_Real* lhs_real;
   SCIP_Real* rhs_real;
   SCIP_Real* val_real;
   int i;

   /* warm start the exact LP by solving the approximate LP first */

   /* allocate and copy aux problem using SCIP_Real arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj_real, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb_real, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub_real, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs_real, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs_real, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val_real, nnonz) );
   for( i = 0; i < nvars; i++)
   {
      obj_real[i] = mpq_get_d(obj[i]);
      lb_real[i] = mpq_get_d(lb[i]);
      ub_real[i] = mpq_get_d(ub[i]);
   }
   for( i = 0; i < nconss; i++)
   {
      lhs_real[i] = mpq_get_d(lhs[i]);
      rhs_real[i] = mpq_get_d(rhs[i]);
   }
   for( i = 0; i < nnonz; i++)
      val_real[i] = mpq_get_d(val[i]);

   /* build and solve approximate aux. problem */
   SCIP_CALL( SCIPlpiCreate(&lpi, SCIPgetMessagehdlr(scip), "problem" , SCIP_OBJSEN_MAXIMIZE) );
   SCIP_CALL( SCIPlpiAddCols(lpi, nvars, obj_real, lb_real, ub_real, colnames, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddRows(lpi, nconss, lhs_real, rhs_real, NULL, nnonz, beg, ind, val_real) );

   SCIP_CALL( SCIPlpiSolveDual(lpi) );

   /* load lp state into exact LP */
   SCIP_CALL( SCIPlpiGetState(lpi, SCIPblkmem(scip), &lpistate) );
   SCIP_CALL( SCIPlpiexSetState(lpiex, SCIPblkmem(scip), lpistate) );

   /* free memory used for approx LP */
   SCIP_CALL( SCIPlpiFreeState(lpi, SCIPblkmem(scip), &lpistate) );
   SCIPfreeBufferArray(scip, &val_real);
   SCIPfreeBufferArray(scip, &rhs_real);
   SCIPfreeBufferArray(scip, &lhs_real);
   SCIPfreeBufferArray(scip, &ub_real);
   SCIPfreeBufferArray(scip, &lb_real);
   SCIPfreeBufferArray(scip, &obj_real);
   if( lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }
   assert(lpi == NULL);
   return SCIP_OKAY;
}

/** subroutine of constructPSdata(); computes S-interior point or ray which is used to do the shifting step */
static
SCIP_RETCODE psComputeSintPointRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_Bool             findintpoint        /**< if true(false) compute S-interior point(ray) */
   )
{
   int i;
   int j;
   int pos;
   int nconss;
   int nvars;
   int nextendedconss; /* number of extended constraints, # of cols in [A',-A',I,-I] */
   int nnonz;
   int indx;
   mpq_t mpqtemp;
   mpq_t alpha;
   mpq_t beta;
   int nobjnz;

   /* lpiex and data used for the aux. problem */
   SCIP_LPIEX* pslpiex;
   SCIP_LPISTATE* lpistate;
   int psnvars;
   mpq_t* psobj;
   mpq_t* pslb;
   mpq_t* psub;
   int psnconss;
   mpq_t* pslhs;
   mpq_t* psrhs;
   int psnnonz;
   int* psbeg;
   int* pslen;
   int* psind;
   mpq_t* psval;
   mpq_t* dualsol;
   mpq_t* primalsol;
   char ** colnames;
   mpq_t objval;

   /* mapping between variables used in the aux. problem and the original problem */
   int ndvarmap;
   int* dvarmap;
   int* dvarincidence;

   nconss = consdata->nconss;
   nvars = consdata->nvars;
   nnonz = consdata->nnonz;
   nextendedconss = conshdlrdata->nextendedconss;

   psnvars = 0;
   psnconss = 0;
   psnnonz = 0;
   pslpiex = NULL;
   psobj = NULL;
   pslb = NULL;
   psub = NULL;
   pslhs = NULL;
   psrhs = NULL;
   psbeg = NULL;
   pslen = NULL;
   psind = NULL;
   psval = NULL;
   dvarmap = NULL;
   dualsol = NULL;
   primalsol = NULL;
   mpq_init(mpqtemp);
   mpq_init(alpha);
   mpq_init(beta);

   /* set up dvarmap - mapping between variables and original problem
    * - use the rows that are used for aux. problem
    * - dvarmap[i] is the index in the original problem of the i^th constraint in the reduced size problem
    *   (reduced from nextendedconss to ndvarmap)
    * - dvarincidence gives the incidence vector of variables used in aux problem
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &dvarmap, nextendedconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dvarincidence, nextendedconss) );
   if( conshdlrdata->psreduceauxlp )
   {
      /* if the aux. lp is reduced, this means that dvarmap is based on the conshdlrdata->includedcons */
      for( i = 0; i < nextendedconss; i++ )
         dvarincidence[i] = conshdlrdata->includedcons[i];
   }
   else
   {
      /* if the aux. lp is not reduced then expand the selection for dvarmap to include all dual vars with finite cost */
      for( i = 0; i < nextendedconss; i++ )
         dvarincidence[i] = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( !isNegInfinity(conshdlrdata, consdata->lhs[i]) )
            dvarincidence[i] = 1;
         if( !isPosInfinity(conshdlrdata, consdata->rhs[i]) )
            dvarincidence[nconss + i] = 1;
      }
      for( i = 0; i < nvars; i++ )
      {
         if( !isNegInfinity(conshdlrdata, consdata->lb[i]) )
            dvarincidence[2*nconss + i] = 1;
         if( !isPosInfinity(conshdlrdata, consdata->ub[i]) )
            dvarincidence[2*nconss + nvars + i] = 1;
      }
   }
   pos = 0;
   for( i = 0; i < nextendedconss; i++ )
   {
      if(dvarincidence[i])
      {
         dvarmap[pos] = i;
         pos++;
      }
   }
   ndvarmap = pos;

   /* if we are finding an interior ray, always use the optimized selection */
   if( conshdlrdata->psintpointselection == 'o' || !findintpoint )
   {
      /* in this case we will find an optimized interior point for which we will try to push it interior and
       * optimize over its objective value.  To do this we will solve the following problem
       * max \alpha * [lhs,-rhs,lb,ub] * y + \beta d
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                  M >= d >= 0
       * M is a bound on how interior we will let the point be, S is the set of dual columns chosen earlier
       * which could have nonzero values for the S-interior point.  The parameter psreduceauxlp=TRUE then we
       * exclude all dual variables y_i that are not in S to not be present in this problem.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */

      psnvars =  ndvarmap + 1;
      psnconss = nvars + conshdlrdata->npsbasis;
      psnnonz = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
            psnnonz += consdata->len[i];
         if( dvarincidence[nconss + i] )
            psnnonz += consdata->len[i];
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
            psnnonz++;
         if( dvarincidence[2*nconss + nvars + i] )
            psnnonz++;
      }
      psnnonz += 2*conshdlrdata->npsbasis;

      /* allocate memory for aux problem */
      SCIP_CALL( SCIPallocBufferArray(scip, &psobj, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psobj[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslb, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(pslb[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psub, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psub[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(pslhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psrhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(psrhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psbeg, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pslen, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psind, psnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psval, psnnonz) );
      for( i = 0; i < psnnonz; i++ )
         mpq_init(psval[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &primalsol, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(primalsol[i]);
      mpq_init(objval);
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, psnvars) );
      for( i = 0; i < psnvars; i++ )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &colnames[i],SCIP_MAXSTRLEN ) );
         (void) SCIPsnprintf( colnames[i] , SCIP_MAXSTRLEN, "var%d",i);
      }

      /* the representation of the problem will be:
       * max:  [\alpha*OBJ, \beta]*[y,d]'
       * s.t.: [c] <= [ A~ |  0]   [y] <= [  c   ]
       *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
       * bounds:     0 <= y <= \infty
       *             0 <= d <= M
       * y is a vector of length (ndvarmap) and d is a single variable
       * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
       * and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
       *
       * beta is set equal to the param psobjweight and alpha is set equal to
       * alpha := (1-beta)/||OBJ||
       */

      /* set up the objective */
      pos = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
         {
            mpq_set(psobj[pos], consdata->lhs[i]);
            pos++;
         }
      }
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[nconss + i] )
         {
            mpq_neg(psobj[pos], consdata->rhs[i]);
            pos++;
         }
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
         {
            mpq_set(psobj[pos], consdata->lb[i]);
            pos++;
         }
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + nvars + i])
         {
            mpq_neg(psobj[pos], consdata->ub[i]);
            pos++;
         }
      }
      assert(pos == ndvarmap);

      /* set alpha and beta. */
      mpq_set_d(alpha, conshdlrdata->psobjweight);
      mpq_set_ui(beta, 1, 1);

      if( mpq_sgn(alpha) > 0 )
      {
         mpq_sub(beta, beta, alpha);

         /*  beta = (1-alpha)*|OBJ|   Where OBJ = optimal objective value of root LP, if |OBJ|<1 use 1 instead */
         if( fabs(SCIPgetLPObjval(scip)) > 1 )
         {
            mpq_set_d(mpqtemp, fabs(SCIPgetLPObjval(scip)));
            mpq_mul(beta, beta, mpqtemp);
         }
         /* divide through by alpha and round beta to be a power of 2 */
         mpq_div(beta, beta, alpha);
         mpq_set_ui(alpha, 1, 1);
         mpq_set_d(beta, pow(2, (int) (log(mpq_get_d(beta))/log(2))));
      }

      /* set objective to normalized value */
      for( i = 0; i < ndvarmap; i ++)
         mpq_mul(psobj[i], psobj[i], alpha);
      mpq_set(psobj[ndvarmap], beta);

      /* set variable bounds */
      for( i = 0; i < ndvarmap; i++ )
      {
         mpq_set(psub[i], conshdlrdata->posinfinity);
         mpq_set_ui(pslb[i], 0, 1);
      }
      mpq_set_ui(psub[ndvarmap], PSBIGM, 1);
      mpq_set_ui(pslb[ndvarmap], 0 ,1);

      /* set up constraint bounds */
      for( i = 0; i < nvars; i++ )
      {
         mpq_set(pslhs[i], consdata->obj[i]);
         mpq_set(psrhs[i], consdata->obj[i]);
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         mpq_set_si(pslhs[nvars + i], 0, 1);
         mpq_set(psrhs[nvars + i], conshdlrdata->posinfinity);
      }

      /* set up constraint matrix: this involves transposing the constraint matrix */

      /* count the length of each constraint */
      for( i = 0; i < psnconss; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j <consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pslen[indx]++;
         }
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         pslen[nvars + i] = 2;
      }
      /* set up the beg array */
      pos = 0;
      for( i = 0; i < psnconss; i++ )
      {
         psbeg[i] = pos;
         pos += pslen[i];
      }
      assert(pos == psnnonz);

      /* reset the length array and build it up as entries are added one by one by scanning through matrix. */
      for( i = 0; i < nvars; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j < consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pos = psbeg[consdata->ind[j]] + pslen[consdata->ind[j]];
               psind[pos] = i;
               if(dvarmap[i]<nconss)
                  mpq_set(psval[pos], consdata->val[j]);
               else
                  mpq_neg(psval[pos], consdata->val[j]);
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pos = psbeg[indx] + pslen[indx];
            psind[pos] = i;
            if( dvarmap[i] < 2*nconss + nvars)
               mpq_set_ui(psval[pos], 1, 1);
            else
               mpq_set_si(psval[pos], -1, 1);
            pslen[indx]++;
         }
      }
      /* set up the last npsbasis rows */
      pos = nvars;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( conshdlrdata->includedcons[indx] )
         {
            psind[psbeg[pos]] = i;
            mpq_set_ui(psval[psbeg[pos]], 1, 1);
            psind[psbeg[pos] + 1] = psnvars - 1;
            mpq_set_si(psval[psbeg[pos] + 1], -1, 1);
            pos++;
         }
      }
      assert( pos == psnconss);

      if( !findintpoint )
      {
         /* in this case we want to find an interior ray instead of an interior point
          * the problem will be modified to the following problem:
          * max:  [OBJ, 0]*[y,d]'
          * s.t.: [0] <= [ A~ |  0]   [y] <= [  0   ]
          *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
          * bounds:     0 <= y <= \infty
          *             1 <= d <= \infty
          * y is a vector of length (ndvarmap) and d is a single variable
          * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
          * and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
          *
          * the parts that change are the objective function, the RHS/LHS of the first constraint set
          * and the lower bound for d
          */

         mpq_set_ui(psobj[ndvarmap], 0, 1);

         /* update the rhs/lhs */
         for( i = 0; i < nvars; i++ )
         {
            mpq_set_ui(pslhs[i], 0, 1);
            mpq_set_ui(psrhs[i], 0, 1);
         }

         /* update bounds on d */
         mpq_set(psub[ndvarmap], conshdlrdata->posinfinity);
         mpq_set_ui(pslb[ndvarmap], 1 ,1);
      }

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, SCIP_OBJSEN_MAXIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psnvars, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnconss, (const mpq_t*) pslhs, (const mpq_t*) psrhs,
            NULL, psnnonz, psbeg, pslen, psind, psval) );

      if( PSWARMSTARTAUXPROB )
      {
         /* warm start the exact LP by solving the approximate LP first */
         SCIP_CALL( warmStartlpiex(scip, pslpiex, psnvars, psnconss, psnnonz, psobj, pslb,
               psub, pslhs, psrhs, colnames, psind, psbeg, psval ));
      }

      /* solve the LP */
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");
         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, primalsol, NULL, NULL, NULL) );

         mpq_set(conshdlrdata->commonslack, primalsol[psnvars - 1]);
         if( mpq_sgn(conshdlrdata->commonslack) == 0 )
         {
            /* if commonslack == 0, point/ray is not interior */
            SCIPdebugMessage("   --> project-and-shift failed to find interior point/ray\n");
         }
         else
         {
            /* assign interior point solution to constraint handler data */
            for( i = 0; i < ndvarmap; i++ )
            {
               if( findintpoint )
                  mpq_set( conshdlrdata->interiorpt[dvarmap[i]], primalsol[i]);
               else
                  mpq_set( conshdlrdata->interiorray[dvarmap[i]], primalsol[i]);

            }
            if( findintpoint )
               conshdlrdata->pshaspoint = TRUE;
            else
               conshdlrdata->pshasray = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr( pslpiex ) );
      }
      if(primalsol != NULL)
      {
         for( i = 0; i < psnvars; i++ )
            mpq_clear(primalsol[i]);
      }
      SCIPfreeBufferArray(scip, &primalsol);
   }
   /* use 'a'rbitrary interior point */
   else if( conshdlrdata->psintpointselection == 'a' )
   {
      SCIPdebugMessage("building aux. problem with arbitrary interior point\n");

      /* the aux problem here that we want to solve can be written in the following way.
       * First let A# be the submatrix of [A',-A',I,-I] defined by dvarmap.  Then we want to solve:
       *
       * max   \sum \delta_i
       * s.t.:  A# * y - c*\lambda = 0
       *             y_i >= \delta_i for each i in S
       *               y_i >= 0
       *           1 >= \delta_i >=0
       *                \lambda >= 1
       *
       * solving this problem determines an interior point to the dual problem (which is y/\lambda)
       * it maximizes the number of components which are interior using the \delta_i's.
       *
       * However, instead of solving it in this form, we construct and solve the dual of this problem
       * which can be written like this:
       *
       *   min      [ 0 | 0 |-1 | 1 ] * [x,y,z,w]'
       *   s.t 0 <= [A#'|-I | 0 | 0 ]              <= 0
       *       1 <= [ 0 | I | 0 | I ] * [x,y,z,w]' <= \infty
       *       0 <= [-c'| 0 |-1 | 0 ]              <= \infty
       *             x free, y,z,w >= 0
       *
       * This problem is solved and the dual multipliers for the first set of rows give us the values of y and the
       * next block of rows tell us which components were nonzero (\delta_i) and the final row tells us what the
       * scale factor \lambda of the c in the original problem was.
       */

      /* figure out the dimensions of the aux problem */
      psnvars = nvars + 2 * ndvarmap + 1;
      psnconss = 2 * ndvarmap + 1;
      nobjnz = 0;

      /* count the number of nonzeros of the aux problem: psnnonz*/
      for( i = 0; i < nvars; i++ )
      {
         if( mpq_sgn(consdata->obj[i]) )
         {
            nobjnz++;
         }
      }
      psnnonz = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
            psnnonz += consdata->len[i];
         if( dvarincidence[nconss + i] )
            psnnonz += consdata->len[i];
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
            psnnonz++;
         if( dvarincidence[2*nconss + nvars + i] )
            psnnonz++;
      }
      psnnonz += nobjnz + 1 + 3 * ndvarmap;

      /* allocate memory for aux problem */
      SCIP_CALL( SCIPallocBufferArray(scip, &psobj, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psobj[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslb, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(pslb[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psub, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psub[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(pslhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psrhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(psrhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psbeg, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pslen, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psind, psnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psval, psnnonz) );
      for( i = 0; i < psnnonz; i++ )
         mpq_init(psval[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &dualsol, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(dualsol[i]);
      mpq_init(objval);
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, psnvars) );
      for( i = 0; i < psnvars; i++ )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &colnames[i],SCIP_MAXSTRLEN ) );
         (void) SCIPsnprintf( colnames[i] , SCIP_MAXSTRLEN, "var%d",i);
      }

      /* set objective */
      for( i = 0; i < nvars + ndvarmap; i++ )
         mpq_set_si( psobj[i ], 0, 1);
      mpq_set_si( psobj[nvars + ndvarmap],-1,1);
      for( i = nvars + ndvarmap + 1; i < psnvars; i++ )
         mpq_set_si( psobj[i], 1, 1);

      /* set variable bounds */
      for( i = 0; i < psnvars; i++ )
         mpq_set(psub[i], conshdlrdata->posinfinity);

      for( i = 0; i < psnvars; i++ )
      {
         if( i < nvars)
            mpq_set( pslb[i], conshdlrdata->neginfinity );
         else
            mpq_set_si( pslb[i], 0, 1);
      }

      /* set up constraint bounds */
      for( i = 0; i < psnconss; i++ )
      {
         if( i < ndvarmap )
         {
            mpq_set_si(pslhs[i], 0, 1);
            mpq_set_si(psrhs[i], 0, 1);
         }
         else if( i == psnconss - 1 )
         {
            mpq_set_si(pslhs[i], 0, 1);
            mpq_set(psrhs[i], conshdlrdata->posinfinity);
         }
         else
         {
            mpq_set_si(pslhs[i], 1, 1);
            mpq_set(psrhs[i], conshdlrdata->posinfinity);
         }
      }

      /* set up constraint matrix */

      /* set up the first ndvarmap rows*/
      pos = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         /* current row comes from lhs/rhs constraints of original problem */
         if( indx < 2*nconss )
         {
            if(indx >= nconss)
               indx -= nconss;
            pslen[i] = consdata->len[indx] + 1;
            psbeg[i] = pos;

            /* set A,-A part of row */
            for(j = 0; j < pslen[i] - 1; j++)
            {
               psind[ psbeg[i] + j ] = consdata->ind[ consdata->beg[indx] + j];
               /* if it is an LHS constraint */
               if( dvarmap[i] < nconss )
                  mpq_set( psval[ psbeg[i] + j], consdata->val[ consdata->beg[indx] + j ] );
               /* otherwise it is the RHS of the constraint */
               else
                  mpq_neg( psval[ psbeg[i] + j], consdata->val[ consdata->beg[indx] + j ] );
            }
            /* set I part of row */
            psind[ psbeg[i] + pslen[i] - 1] = nvars + i ;
            mpq_set_si( psval[psbeg[i] + pslen[i] - 1 ], -1, 1);

            /* update pos */
            pos += consdata->len[indx] + 1;
         }
         /* current row comes from lower bound constraints of original problem */
         else if ( indx < 2*nconss + nvars )
         {
            indx -= 2 * nconss;
            psbeg[i] = pos;
            pslen[i] = 2;
            psind[pos] = indx;
            psind[pos + 1] = nvars + i;
            mpq_set_si(psval[pos], 1, 1);
            mpq_set_si(psval[pos+1], -1, 1);
            pos += 2;
         }
         /* current row comes from upper bound constraints of original problem */
         else
         {
            indx -= (2*nconss + nvars);
            psbeg[i] = pos;
            pslen[i] = 2;
            psind[pos] = indx;
            psind[pos + 1] = nvars + i;
            mpq_set_si(psval[pos], -1, 1);
            mpq_set_si(psval[pos+1], -1, 1);
            pos += 2;
         }
      }

      /* set up the next ndvarmap rows */
      for( i = 0; i < ndvarmap; i++ )
      {
         psbeg[ndvarmap + i] = pos;
         pslen[ndvarmap + i] = 2;
         psind[pos] = nvars + i;
         psind[pos + 1] = nvars + ndvarmap + 1 + i;
         mpq_set_si(psval[pos], 1,1);
         mpq_set_si(psval[pos+1],1,1);
         pos += 2;
      }

      /* last row */
      psbeg[ psnconss - 1 ] = pos;
      pslen[ psnconss - 1 ] = nobjnz + 1; /* this = objective length + 1 */
      j = psbeg[ 2 * nextendedconss ];
      for( i = 0; i < nvars; i++ )
      {
         if( mpq_sgn(consdata->obj[i]) )
         {
            mpq_neg(psval[pos], consdata->obj[i]);
            psind[pos] = i;
            pos++;
         }
      }
      mpq_set_si(psval[pos], -1,1);
      psind[pos] = nvars + ndvarmap;
      pos++;
      assert(pos == psnnonz);

      SCIPdebugMessage("building LPIEX for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psnvars, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnconss, (const mpq_t*) pslhs, (const mpq_t*) psrhs,
            NULL, psnnonz, psbeg, pslen, psind, psval) );

      if( PSWARMSTARTAUXPROB )
      {
         /* warm start the exact LP by solving the approximate LP first */
         SCIP_CALL( warmStartlpiex(scip, pslpiex, psnvars, psnconss, psnnonz, psobj, pslb,
               psub, pslhs, psrhs, colnames, psind, psbeg, psval ));
      }
      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, NULL, dualsol, NULL, NULL) );
         if( mpq_sgn(dualsol[psnconss - 1]) )
            mpq_inv(conshdlrdata->commonslack,dualsol[psnconss - 1]);
         else
         {
            mpq_set_si(conshdlrdata->commonslack, 0, 1);
         }
         if( mpq_sgn(conshdlrdata->commonslack) == 0 )
         {
            SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
         }

         /* interior point is y/lambda */
         conshdlrdata->pshaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( conshdlrdata->includedcons[dvarmap[i]] && mpq_sgn(dualsol[i]) == 0)
            {
               SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
               conshdlrdata->pshaspoint = FALSE;
               i = ndvarmap;
            }
            else
               mpq_div( conshdlrdata->interiorpt[dvarmap[i]], dualsol[i], dualsol[psnconss - 1]);
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr( pslpiex ) );
      }

      /* free the dual solution that was allocated for this problem */
      if( dualsol != NULL )
      {
         for( i = 0; i < psnconss; i++ )
            mpq_clear(dualsol[i]);
         SCIPfreeBufferArray(scip, &dualsol);
      }
   }
   /* use 'A'rbitrary interior point in transposed form*/
   else if( conshdlrdata->psintpointselection == 'A' )
   {
      SCIPdebugMessage("building new version of arbitrary interior point aux. problem\n");

      /* the aux problem here that we want to solve can be written in the following way.
       * First let A# be the submatrix of [A',-A',I,-I] defined by dvarmap.  Then we want to solve:
       *
       * max   \sum \delta_i
       * s.t.:  A# * y - c*\lambda = 0
       *             y_i >= \delta_i for each i in S
       *               y_i >= 0
       *           1 >= \delta_i >=0
       *                \lambda >= 1
       *
       *  the representation of the problem will be:
       * min:         [  0 | 0 | -1 ] * [y,z,w]'
       * s.t.: [0] <= [ A~ | -c|  0 ]   [y] <= [  0   ]
       *       [0] <= [ I* | 0 | -I*] * [z] <= [\infty] <-- only for dual vars from includecons
       *                                [w]
       * bounds:     0 <= y <= \infty
       *             1 <= z <= \infty
       *             0 <= w <= 1
       * y is a vector of length (ndvarmap) and d is a single variable
       * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
       */
      psnvars =  ndvarmap + 1 + conshdlrdata->npsbasis;
      psnconss = nvars + conshdlrdata->npsbasis;
      psnnonz = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
            psnnonz += consdata->len[i];
         if( dvarincidence[nconss + i] )
            psnnonz += consdata->len[i];
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
            psnnonz++;
         if( dvarincidence[2*nconss + nvars + i] )
            psnnonz++;
      }
      psnnonz += 2*conshdlrdata->npsbasis + nvars;

      /* allocate memory for aux problem */
      SCIP_CALL( SCIPallocBufferArray(scip, &psobj, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psobj[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslb, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(pslb[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psub, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psub[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(pslhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psrhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(psrhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psbeg, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pslen, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psind, psnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psval, psnnonz) );
      for( i = 0; i < psnnonz; i++ )
         mpq_init(psval[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &primalsol, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(primalsol[i]);
      mpq_init(objval);
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, psnvars) );
      for( i = 0; i < psnvars; i++ )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &colnames[i],SCIP_MAXSTRLEN ) );
         (void) SCIPsnprintf( colnames[i] , SCIP_MAXSTRLEN, "var%d",i);
      }

      /* set up the objective*/
      for( i = 0; i < ndvarmap + 1; i++ )
      {
         mpq_set_ui(psobj[i], 0, 1);
      }
      for( i = ndvarmap + 1; i < psnvars; i++ )
      {
         mpq_set_si(psobj[i], -1, 1);
      }

      /* set variable bounds */
      for( i = 0; i < ndvarmap; i++ )
      {
         mpq_set(psub[i], conshdlrdata->posinfinity);
         mpq_set_ui(pslb[i], 0, 1);
      }
      mpq_set(psub[ndvarmap], conshdlrdata->posinfinity);
      mpq_set_ui(pslb[ndvarmap], 1 ,1);
      for( i = ndvarmap + 1; i < psnvars; i++ )
      {
         mpq_set_ui(psub[i], 1, 1);
         mpq_set_ui(pslb[i], 0, 1);
      }

      /* set up constraint bounds */
      for( i = 0; i < nvars; i++ )
      {
         mpq_set_ui(pslhs[i], 0, 1);
         mpq_set_ui(psrhs[i], 0, 1);
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         mpq_set_ui(pslhs[nvars + i], 0, 1);
         mpq_set(psrhs[nvars + i], conshdlrdata->posinfinity);
      }

      /* set up constraint matrix: this involves transposing the constraint matrix */
      SCIPdebugMessage("setting up constraint matrix\n");

      /* count the length of each constraint */
      for( i = 0; i < psnconss; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j <consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pslen[indx]++;
         }
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         pslen[nvars + i] = 2;
      }

      /* add another element to the first nvar rows for the c vector */
      for( i = 0; i < nvars; i++ )
      {
         pslen[i]++;
      }

      /* set up the beg array */
      pos = 0;
      for( i = 0; i < psnconss; i++ )
      {
         psbeg[i] = pos;
         pos += pslen[i];
      }
      assert(pos == psnnonz);

      /* reset the length array and build it up as entries are added one by one by scanning through matrix. */
      for( i = 0; i < nvars; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j < consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pos = psbeg[consdata->ind[j]] + pslen[consdata->ind[j]];
               psind[pos] = i;
               if(dvarmap[i]<nconss)
                  mpq_set(psval[pos], consdata->val[j]);
               else
                  mpq_neg(psval[pos], consdata->val[j]);
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pos = psbeg[indx] + pslen[indx];
            psind[pos] = i;
            if( dvarmap[i] < 2*nconss + nvars)
               mpq_set_ui(psval[pos], 1, 1);
            else
               mpq_set_si(psval[pos], -1, 1);
            pslen[indx]++;
         }
      }
      for( i = 0; i < nvars; i++ )
      {
         mpq_neg(psval[psbeg[i] + pslen[i]], consdata->obj[i]);
         psind[psbeg[i] + pslen[i]] = ndvarmap;
         pslen[i]++;
      }

      /* set up the last npsbasis rows */
      pos = nvars;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( conshdlrdata->includedcons[indx] )
         {
            psind[psbeg[pos]] = i;
            mpq_set_ui(psval[psbeg[pos]], 1, 1);
            psind[psbeg[pos] + 1] = psnvars - psnconss + pos;
            mpq_set_si(psval[psbeg[pos] + 1], -1, 1);
            pos++;
         }
      }
      assert( pos == psnconss);

      SCIPdebugMessage("building LPIEX for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psnvars, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnconss, (const mpq_t*) pslhs, (const mpq_t*) psrhs,
            NULL, psnnonz, psbeg, pslen, psind, psval) );

      if( PSWARMSTARTAUXPROB )
      {
         /* warm start the exact LP by solving the approximate LP first */
         SCIP_CALL( warmStartlpiex(scip, pslpiex, psnvars, psnconss, psnnonz, psobj, pslb,
               psub, pslhs, psrhs, colnames, psind, psbeg, psval ));
      }

      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, primalsol, NULL, NULL, NULL) );
         if( mpq_sgn(primalsol[ndvarmap]) )
            mpq_inv(conshdlrdata->commonslack, primalsol[ndvarmap]);
         else
         {
            mpq_set_si(conshdlrdata->commonslack, 0, 1);
         }
         if( mpq_sgn(conshdlrdata->commonslack) == 0 )
         {
            SCIPdebugMessage("   --> interior point not found\n");
         }

         /* interior point is y/lambda */
         conshdlrdata->pshaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( conshdlrdata->includedcons[dvarmap[i]] && mpq_sgn(primalsol[i]) == 0)
            {
               conshdlrdata->pshaspoint = FALSE;
               SCIPdebugMessage("   --> interior point not found\n");
               i = ndvarmap;
            }
            else
               mpq_div( conshdlrdata->interiorpt[dvarmap[i]], primalsol[i], primalsol[ndvarmap]);
         }

#ifdef PS_OUT
         printf("constraints all satisfied by slack=");
         mpq_out_str(stdout, 10, conshdlrdata->commonslack);
         printf("\n");

         printf("objective value of aux problem=");
         mpq_out_str(stdout, 10, objval);
         printf("\n");

         printf("relative interior solution:\n");
         for( i = 0; i <  conshdlrdata->nextendedconss; i++ )
         {
            printf("   i=%d: ", i);
            mpq_out_str(stdout, 10, conshdlrdata->interiorpt[i]);
            printf("\n");
         }
#endif
      }
      else
      {
         printlpiexerr( pslpiex );
         conshdlrdata->psdatafail = TRUE;
      }
      /* free the dual solution that was allocated for this problem */
      if(primalsol != NULL)
      {
         for( i = 0; i < psnvars; i++ )
            mpq_clear(primalsol[i]);
         SCIPfreeBufferArray(scip, &primalsol);
      }
   }
   else if( conshdlrdata->psintpointselection == 't' )
   {
      /* in this case we will find an optimized interior point for which we will try to push it interior and
       * optimize over its objective value.  To do this we will solve the following two problems
       * max                                   d
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                  M >= d >= 0
       * max          [lhs,-rhs,lb,ub] * y
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                       d >= d* <-- where d* is the optimal solution from the first problem
       * M is a bound on how interior we will let the point be, S is the set of dual columns chosen earlier
       * which could have nonzero values for the S-interior point.  The parameter psreduceauxlp=TRUE then we
       * exclude all dual variables y_i that are not in S to not be present in this problem.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */
      psnvars =  ndvarmap + 1;
      psnconss = nvars + conshdlrdata->npsbasis;
      psnnonz = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
            psnnonz += consdata->len[i];
         if( dvarincidence[nconss + i] )
            psnnonz += consdata->len[i];
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
            psnnonz++;
         if( dvarincidence[2*nconss + nvars + i] )
            psnnonz++;
      }
      psnnonz += 2*conshdlrdata->npsbasis;

      /* allocate memory for aux problem */
      SCIP_CALL( SCIPallocBufferArray(scip, &psobj, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psobj[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslb, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(pslb[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psub, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(psub[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &pslhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(pslhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psrhs, psnconss) );
      for( i = 0; i < psnconss; i++ )
         mpq_init(psrhs[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &psbeg, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pslen, psnconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psind, psnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &psval, psnnonz) );
      for( i = 0; i < psnnonz; i++ )
         mpq_init(psval[i]);
      SCIP_CALL( SCIPallocBufferArray(scip, &primalsol, psnvars) );
      for( i = 0; i < psnvars; i++ )
         mpq_init(primalsol[i]);
      mpq_init(objval);
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, psnvars) );
      for( i = 0; i < psnvars; i++ )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &colnames[i],SCIP_MAXSTRLEN ) );
         (void) SCIPsnprintf( colnames[i] , SCIP_MAXSTRLEN, "var%d",i);
      }

      /* the representation of the problem will be:
       * max:              [0,1]*[y|d]'
       * s.t.: [c] <= [ A~ |  0]   [y] <= [  c   ]
       *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
       * bounds:     0 <= y <= \infty
       *             0 <= d <= M
       * y is a vector of length (ndvarmap) and d is a single variable and A~ is the submatrix of [A',-A',I,-I] using
       * columns in dvarmap and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
       */


      /* solve the aux problem in two stages, first maximize the interiorness of the point, and then in the second phase
       * problem, move the interiorness to the constraint bounds and optimize over the objective.
       */
      for( i = 0; i < ndvarmap; i ++)
         mpq_set_ui(psobj[i], 0, 1);
      mpq_set_ui(psobj[ndvarmap], 1, 1);

      /* set variable bounds */
      for( i = 0; i < ndvarmap; i++ )
      {
         mpq_set(psub[i], conshdlrdata->posinfinity);
         mpq_set_ui(pslb[i], 0, 1);
      }
      mpq_set_ui(psub[ndvarmap], PSBIGM, 1);
      mpq_set_ui(pslb[ndvarmap], 0 ,1);

      /* set up constraint bounds */
      for( i = 0; i < nvars; i++ )
      {
         mpq_set(pslhs[i], consdata->obj[i]);
         mpq_set(psrhs[i], consdata->obj[i]);
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         mpq_set_si(pslhs[nvars + i], 0, 1);
         mpq_set(psrhs[nvars + i], conshdlrdata->posinfinity);
      }

      /* set up constraint matrix: this involves transposing the constraint matrix */

      /* count the length of each constraint */
      for( i = 0; i < psnconss; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j <consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pslen[indx]++;
         }
      }
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         pslen[nvars + i] = 2;
      }

      /* set up the beg array */
      pos = 0;
      for( i = 0; i < psnconss; i++ )
      {
         psbeg[i] = pos;
         pos += pslen[i];
      }
      assert(pos == psnnonz);

      /* reset the length array and build it up as entries are added one by one by scanning through matrix */
      for( i = 0; i < nvars; i++ )
         pslen[i] = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( indx < 2*nconss )
         {
            if( indx >= nconss )
               indx -= nconss;
            for(j = consdata->beg[indx]; j < consdata->beg[indx] + consdata->len[indx]; j++)
            {
               pos = psbeg[consdata->ind[j]] + pslen[consdata->ind[j]];
               psind[pos] = i;
               if(dvarmap[i]<nconss)
                  mpq_set(psval[pos], consdata->val[j]);
               else
                  mpq_neg(psval[pos], consdata->val[j]);
               pslen[consdata->ind[j]]++;
            }
         }
         else
         {
            if ( indx < 2*nconss + nvars )
               indx -= 2 * nconss;
            else
               indx -= (2*nconss + nvars);
            pos = psbeg[indx] + pslen[indx];
            psind[pos] = i;
            if( dvarmap[i] < 2*nconss + nvars)
               mpq_set_ui(psval[pos], 1, 1);
            else
               mpq_set_si(psval[pos], -1, 1);
            pslen[indx]++;
         }
      }

      /* set up the last npsbasis rows */
      pos = nvars;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         if( conshdlrdata->includedcons[indx] )
         {
            psind[psbeg[pos]] = i;
            mpq_set_ui(psval[psbeg[pos]], 1, 1);
            psind[psbeg[pos] + 1] = psnvars - 1;
            mpq_set_si(psval[psbeg[pos] + 1], -1, 1);
            pos++;
         }
      }
      assert( pos == psnconss);

      /* build aux LP using the exact LP interface */
      if( pslpiex != NULL )
      {
         SCIP_CALL( SCIPlpiexFree(&pslpiex) );
      }
      pslpiex = NULL;
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, SCIP_OBJSEN_MAXIMIZE) );

#ifdef PS_OUT
      /* activate extra output from exact LP solver */
      SCIPlpiexSetIntpar(pslpiex, SCIP_LPPAR_LPINFO, TRUE);
#endif

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psnvars, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnconss, (const mpq_t*) pslhs, (const mpq_t*) psrhs,
            NULL, psnnonz, psbeg, pslen, psind, psval) );
      if( PSWARMSTARTAUXPROB )
      {
         /* warm start the exact LP by solving the approximate LP first */
         SCIP_CALL( warmStartlpiex(scip, pslpiex, psnvars, psnconss, psnnonz, psobj, pslb,
               psub, pslhs, psrhs, colnames, psind, psbeg, psval ));
      }

      /* solve the LP */
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      /* get state and solution of lpiex that was just solved */
      SCIP_CALL( SCIPlpiexGetState(pslpiex, SCIPblkmem(scip), &lpistate) );
      SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, NULL, NULL, NULL, NULL) );

      /* now reset the objective value to be the original objective */
      pos = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[i] )
         {
            mpq_set(psobj[pos], consdata->lhs[i]);
            pos++;
         }
      }
      for( i = 0; i < nconss; i++ )
      {
         if( dvarincidence[nconss + i] )
         {
            mpq_neg(psobj[pos], consdata->rhs[i]);
            pos++;
         }
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + i] )
         {
            mpq_set(psobj[pos], consdata->lb[i]);
            pos++;
         }
      }
      for( i = 0; i < nvars; i++ )
      {
         if( dvarincidence[2*nconss + nvars + i])
         {
            mpq_neg(psobj[pos], consdata->ub[i]);
            pos++;
         }
      }
      assert(pos == ndvarmap);
      mpq_set_ui(psobj[ndvarmap], 0, 1);

      /* set the lower bound on the interiorness based on the objective value */
      mpq_set(pslb[ndvarmap], objval);

      /* reuse the psind array to pass indices to updated the bounds and objective */
      for( i = 0; i < psnvars; i++ )
         psind[i] = i;
      SCIP_CALL( SCIPlpiexChgBounds(pslpiex, psnvars, psind, pslb, NULL) );
      SCIP_CALL( SCIPlpiexChgObj(pslpiex, psnvars, psind, psobj) );

      /* reload state and solve new LP */
      SCIP_CALL( SCIPlpiexSetState(pslpiex, SCIPblkmem(scip), lpistate) );

      /* reoptimizing using primal simplex seems MUCH faster here, warm start basis is primal feasible */
      SCIP_CALL( SCIPlpiexSolvePrimal(pslpiex) );
      SCIP_CALL( SCIPlpiexFreeState(pslpiex, SCIPblkmem(scip), &lpistate) );

        /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, primalsol, NULL, NULL, NULL) );

         /* assign interior point solution to constraint handler data */
         mpq_set(conshdlrdata->commonslack, primalsol[psnvars - 1]);
         if( mpq_sgn(conshdlrdata->commonslack) == 0 )
         {
            SCIPdebugMessage("   --> interior point not found \n");
         }
         else
         {
            for( i = 0; i < ndvarmap; i++ )
            {
               mpq_set( conshdlrdata->interiorpt[dvarmap[i]], primalsol[i]);
            }
            conshdlrdata->pshaspoint = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr(pslpiex) );
      }
     if( primalsol != NULL )
      {
         for( i = 0; i < psnvars; i++ )
            mpq_clear(primalsol[i]);
      }
      SCIPfreeBufferArray(scip, &primalsol);
   }
   else
   {
      SCIPerrorMessage("invalid parameter setting <%c> for selection method to compute interior point\n",
         conshdlrdata->psintpointselection);
      return SCIP_PARAMETERWRONGVAL;
   }

   mpq_clear(beta);
   mpq_clear(alpha);
   mpq_clear(mpqtemp);
   for( i = 0; i < psnvars; i++ )
      mpq_clear(psobj[i]);
   for( i = 0; i < psnvars; i++ )
      mpq_clear(pslb[i]);
   for( i = 0; i < psnvars; i++ )
      mpq_clear(psub[i]);
   for( i = 0; i < psnconss; i++ )
      mpq_clear(pslhs[i]);
   for( i = 0; i < psnconss; i++ )
      mpq_clear(psrhs[i]);
   for( i = 0; i < psnnonz; i++ )
      mpq_clear(psval[i]);
   mpq_clear(objval);
   if( pslpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&pslpiex) );
   }
   assert(pslpiex == NULL);
   for( i = psnvars - 1; i >= 0; i--)
      SCIPfreeBufferArray(scip, &colnames[i] );
   SCIPfreeBufferArray(scip, &colnames);
   SCIPfreeBufferArray(scip, &psval);
   SCIPfreeBufferArray(scip, &psind);
   SCIPfreeBufferArray(scip, &pslen);
   SCIPfreeBufferArray(scip, &psbeg);
   SCIPfreeBufferArray(scip, &psrhs);
   SCIPfreeBufferArray(scip, &pslhs);
   SCIPfreeBufferArray(scip, &psub);
   SCIPfreeBufferArray(scip, &pslb);
   SCIPfreeBufferArray(scip, &psobj);
   SCIPfreeBufferArray(scip, &dvarincidence);
   SCIPfreeBufferArray(scip, &dvarmap);
   return SCIP_OKAY;
}


/** constructs datas used to compute dual bounds by the project-and-shift method */
static
SCIP_RETCODE constructPSData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   int i;
   assert(consdata->nconss > 0);

   /* consider the primal problem as
    * min c'x
    * lhs <= Ax <= rhs
    *  lb <=  x <= ub
    *
    * and the dual of the form
    *     [  A',  -A',  I,  -I] y =  c
    *                           y >= 0
    *
    * A subset S of the dual columns are chosen to give a submatrix D of [A',-A',I,-I], which is then LU factorized using
    * rectlu code then an S-interior point is found (a dual solution that is strictly positive for each column in S).
    * this data is then reused throughout the tree where the LU factorization can be used to correct feasibility of
    * the equality constraints of the dual, and a convex combination with the S-interior point can correct any
    * infeasibility coming from negative variables.
    */

   /* if the ps data was already constructed, exit */
   if( conshdlrdata->psdatacon )
      return SCIP_OKAY;
   /* now mark that this function has been called */
   conshdlrdata->psdatacon = TRUE;

   SCIPdebugMessage("calling constructPSdata()\n");

   /* process the bound changes */
   processBoundchgs(scip, conshdlrdata, consdata);

   conshdlrdata->nextendedconss = 2*consdata->nconss + 2*consdata->nvars;

   /* call function to select the set S */
   SCIP_CALL( psChooseS(scip, conshdlrdata, consdata) );

   /* compute LU factorization of D == A|_S */
   SCIP_CALL( psFactorizeD(scip, conshdlrdata, consdata) );

   /* if no fail in LU factorization, compute S-interior point and/or ray */
   if( !conshdlrdata->psdatafail )
   {
      if( SCIPpsInfeasRay(scip) || !conshdlrdata->psuseintpoint )
      {
         /* try to compute the S-interior ray if we want to use it for bounding or infeasibility */
         SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->interiorray, conshdlrdata->nextendedconss) );
         for( i = 0; i < conshdlrdata->nextendedconss; i++ )
            mpq_init(conshdlrdata->interiorray[i]);
         SCIP_CALL( psComputeSintPointRay(scip, conshdlrdata, consdata, FALSE) );
      }
      if( conshdlrdata->psuseintpoint || !conshdlrdata->pshasray)
      {
         /* now, compute S-interior point if we need it OR if the ray construction failed */
         SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->interiorpt, conshdlrdata->nextendedconss) );
         for( i = 0; i < conshdlrdata->nextendedconss; i++ )
            mpq_init(conshdlrdata->interiorpt[i]);
         SCIP_CALL( psComputeSintPointRay(scip, conshdlrdata, consdata, TRUE) );
      }
   }

   /* if construction of both point and ray has failed, mark psdatafail as true. */
   if( !conshdlrdata->pshaspoint && !conshdlrdata->pshasray )
   {
      conshdlrdata->psdatafail = TRUE;
   }

   SCIPdebugMessage("exiting constructPSdata()\n");

   return SCIP_OKAY;
}

/** updates root node LP activities of the constraints */
static
SCIP_RETCODE copyRootActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int nrows;
   int ncols;
   int nconss;
   int nvars;
   int nextendedconss;
   int i;

   assert(SCIPgetDepth(scip) >= 0);

   /* check if the LP has been solved to optimality */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only update LP activities at the root node; keep updating after adding cutting planes */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   nconss = consdata->nconss;
   nvars = consdata->nvars;
   nextendedconss = 2*nconss + 2*nvars;

   /* allocate memory */
   if( !(conshdlrdata->rootactivitycon) )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->rootactivity, nextendedconss) );
      conshdlrdata->rootactivitycon = TRUE;
   }

   /* reset LP activities */
   for( i = 0; i < nextendedconss; i++ )
      conshdlrdata->rootactivity[i] = 0;

   /* update constraint activities */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows == nconss);
   for( i = 0; i < nconss; i++ )
   {
      if( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, rows[i]), SCIProwGetLhs(rows[i])) )
         conshdlrdata->rootactivity[i] = 1;
      if( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, rows[i]), SCIProwGetRhs(rows[i])) )
         conshdlrdata->rootactivity[nconss + i] = 1;
   }

   /* update bound activities */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert(ncols == nvars);
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPisFeasEQ(scip, SCIPcolGetPrimsol(cols[i]), SCIPcolGetLb(cols[i])) )
         conshdlrdata->rootactivity[2*nconss + i] = 1;
      if( SCIPisFeasEQ(scip, SCIPcolGetPrimsol(cols[i]), SCIPcolGetUb(cols[i])) )
         conshdlrdata->rootactivity[2*nconss + nvars + i] = 1;
   }

   return SCIP_OKAY;
}

#if 0 /* uncomment if method is needed */
/* checks the condition required for the project-and-shift method to work given the current settings of the parameters.
 * this will call constructPSData() if the psdata structures are not constructed yet.
 */
static
SCIP_RETCODE checkPSConditions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_Bool*            success             /**< do the conditions hold? */
   )
{
   /* if constructPSData was not already constructed, call constructPSData */
   if( !conshdlrdata->psdatacon )
   {
      SCIP_CALL( constructPSData(scip, conshdlrdata, consdata) );
   }

   /* check if constructPSData was successful */
   if( conshdlrdata->psdatafail )
      *success = FALSE;
   else
      *success = TRUE;

   return SCIP_OKAY;
}
#endif

/** computes safe dual bound by project-and-shift method */
static
SCIP_RETCODE getPSdualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   mpq_t*                boundval,           /**< value of dual bound */
   SCIP_Bool*            success             /**< pointer to store whether dual bound computation was successful */
   )
{
   SCIP_COL** cols;
   mpq_t* approxdualsol;
   mpq_t* costvect; /**< dual cost vector for expanded problem*/
   mpq_t* violation;
   mpq_t* correction;
   SCIP_Bool useinteriorpoint;
   mpq_t mpqtemp;
   mpq_t mpqtemp2;
   mpq_t lambda1;
   mpq_t lambda2;
   mpq_t maxv;
   mpq_t dualbound;
   int i;
   int j;
   int rval;
   int nextendedconss;
   int nconss;
   int nvars;
   int currentrow;
   int isfeas;

   /* project-and-shift method:
    * 1. projection step (to ensure that equalities are satisfied):
    *   - compute error in equalities: r=c-Ay^
    *   - backsolve system of equations to find correction of error: z with Dz=r
    *   - add corretion to approximate dual solution: bold(y)=y^+[z 0]
    * 2. shifing step (to ensure that inequalities are satisfied):
    *   - take convex combination of projected approximate point bold(y) with interior point y*
    * 3. compute dual objective value of feasible dual solution and set bound
    */

   /* constructpsdata should always be called first */
   assert(conshdlrdata->psdatacon);

   /* if data has not been constructed, or it failed, then exit */
   if( conshdlrdata->psdatafail )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   SCIPdebugMessage("calling getPSdualbound()\n");

   /* decide if we should use ray or point to compute bound */
   if( conshdlrdata->psuseintpoint && conshdlrdata->pshaspoint )
   {
      /* if we are supposed to use the interior point, and it exists use it */
      useinteriorpoint = TRUE;
   }
   else
   {
      /* in this case, since psdatafail != TRUE, pshasray should be true -- use it */
      assert( conshdlrdata->pshasray );
      useinteriorpoint = FALSE;
   }

   mpq_init(mpqtemp);
   mpq_init(mpqtemp2);
   mpq_init(lambda1);
   mpq_init(lambda2);
   mpq_init(maxv);
   mpq_init(dualbound);

   nextendedconss = conshdlrdata->nextendedconss;
   nconss = consdata->nconss;
   nvars = consdata->nvars;

   /* process bound changes */
   processBoundchgs(scip, conshdlrdata, consdata);

   /* allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIP_CALL( SCIPallocBufferArray(scip, &approxdualsol, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(approxdualsol[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &costvect, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(costvect[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &violation, nvars) );
   for( i = 0; i < nvars; i++ )
      mpq_init(violation[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &correction, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(correction[i]);

   /* recover the objective coefs and approximate solution value of dual solution;
    * dual vars of lhs constraints (including -inf) and rhs constraints (including +inf),
    * dual vars of lb constraint (including -inf) and ub constraints (including +inf)
    */
   for( i = 0; i < nconss; i++ )
   {
      /* in case we want to prove infeasibility it might be that we were not able to compute a dual solution
       * with bound exceeding objective value; in this case dual solution is set to SCIP_INVALID
       */
      if( SCIProwGetDualsol(consdata->rows[i]) == SCIP_INVALID )
      {
         SCIPdebugMessage("  no valid unbounded approx dual sol given\n");
         *success = FALSE;
         goto TERMINATE;
      }
      mpq_set_d(mpqtemp, SCIProwGetDualsol(consdata->rows[i]));

      /* lhs constraint */
      if( mpq_sgn(mpqtemp) > 0 )
         mpq_set(approxdualsol[i], mpqtemp);
      /* rhs constraint */
      else
         mpq_neg(approxdualsol[i+nconss], mpqtemp);

      mpq_set(costvect[i], consdata->lhs[i]);
      mpq_neg(costvect[i+nconss], consdata->rhs[i]);
   }
   cols = SCIPgetLPCols(scip);
   for( i = 0; i < nvars; i++ )
   {
      mpq_set_d(mpqtemp, SCIPgetColRedcost(scip, cols[i]));
      /* lb constraint */
      if( mpq_sgn(mpqtemp) > 0 )
         mpq_set(approxdualsol[i + 2*nconss], mpqtemp);
      /* ub constraint */
      else
         mpq_neg(approxdualsol[i + 2*nconss + nvars], mpqtemp);

      mpq_set(costvect[i + 2*nconss], consdata->lbloc[i]);
      mpq_neg(costvect[i + 2*nconss + nvars], consdata->ubloc[i]);
   }

#ifdef PS_OUT
   printf("approximate dual solution:\n");

   mpq_set_ui(dualbound, 0, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualsol[i]);
      printf("\n");

      mpq_mul(mpqtemp, approxdualsol[i], costvect[i]);
      mpq_add(dualbound, dualbound,mpqtemp);
   }

   printf("   objective value=%.20f (", mpq_get_d(dualbound));
   mpq_out_str(stdout, 10, dualbound);
   printf(")\n");
#endif

   /* first, ensure nonnegativity of dual solution and fix artificial dual variables to zero */
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualsol[i]) < 0)
         mpq_set_si(approxdualsol[i], 0, 1);
      else if( isNegInfinity(conshdlrdata, costvect[i]) )
      {
         assert(!conshdlrdata->includedcons[i]);
         mpq_set_si(approxdualsol[i], 0, 1);
      }
   }

   /* calculate violation of equality constraints r=c-A^ty */
   for( i = 0; i < nvars; i++ )
   {
      mpq_set(violation[i], consdata->obj[i]);
   }

   /* A^ty for y corresponding to primal constraints */
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++)
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualsol[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualsol[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }

   /* A^ty for y corresponding to bound constraints */
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualsol[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualsol[i+2*nconss+nvars]);
   }


   /* project solution */

#ifdef PS_OUT
   printf("violation of solution:\n");
   for( i = 0; i < nvars; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, violation[i]);
      printf("\n");
   }
#endif

   /* if there is no violation of the constraints, then skip the projection */
   isfeas = 1;
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) )
      {
         isfeas = 0;
         i = nvars;
      }
   }

   /* isfeas is equal to one only if approximate dual solution is already feasible for the dual */
   if( !isfeas )
   {
      /* compute [z] with Dz=r (D depends on conshdlrdata->psdualcolselection) */
      rval = RECTLUsolveSystem( conshdlrdata->rectfactor, nvars, nextendedconss, violation, correction);

#ifdef PS_OUT
      printf("correction of solution:\n");
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         printf("   i=%d: ", i);
         mpq_out_str(stdout, 10, correction[i]);
         printf(", position=%d\n", conshdlrdata->psbasis[i]);
      }
#endif

      /* projection step: compute bold(y)=y^+[z 0];
       * correct only components corresponding to D (npsbasis=#of columns in D
       */
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         mpq_add(approxdualsol[conshdlrdata->psbasis[i]], approxdualsol[conshdlrdata->psbasis[i]], correction[i]);
      }
   }

#ifdef PS_OUT
   printf("updated dual solution:\n");
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualsol[i]);
      printf("\n");
   }
#endif

   if( useinteriorpoint )
   {
      /* shifting step (scale solution with interior point to be dual feasible):
       * y' = lambda1 bold(y) + lambda2 y*, where
       *   lambda1 = ( slack of int point)/ (slack of int point + max violation) = d/m+d
       *   lambda2 = 1 - lambda1
       */

      /* compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
      mpq_set_ui(lambda1, 1, 1);
      for( i = 0; i < nextendedconss; i++ )
      {
         if( mpq_sgn(approxdualsol[i]) < 0 )
         {
            mpq_set(mpqtemp2, conshdlrdata->interiorpt[i]);
            mpq_sub(mpqtemp, conshdlrdata->interiorpt[i], approxdualsol[i]);
            mpq_div(mpqtemp2, mpqtemp2, mpqtemp);
            if( mpq_cmp(lambda1, mpqtemp2) > 0 )
               mpq_set(lambda1, mpqtemp2);
         }
      }
      mpq_set_si(lambda2, 1, 1);
      mpq_sub(lambda2, lambda2, lambda1);
   }
   else
   {
      /* in this case we are using an interior ray that can be added freely to the solution */
      mpq_set_si(lambda1, 1, 1);

      /* compute lambda values: compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
      mpq_set_ui(lambda1, 1, 1);
      for( i = 0; i < nextendedconss; i++ )
      {
         if( mpq_sgn(approxdualsol[i]) < 0 && conshdlrdata->includedcons[i] )
         {
            mpq_div(mpqtemp, approxdualsol[i], conshdlrdata->interiorray[i]);
            mpq_neg(mpqtemp, mpqtemp);
            if( mpq_cmp(lambda2, mpqtemp) < 0 ) /* changed to < */
               mpq_set(lambda2, mpqtemp);
         }
      }
   }

   /* perform shift */
   if( mpq_sgn(lambda2) != 0 )
   {
      for( i = 0; i < nextendedconss; i++ )
      {
         mpq_mul(approxdualsol[i], approxdualsol[i], lambda1);
      }
      for( i =0; i < nextendedconss; i++ )
      {
         if( useinteriorpoint )
            mpq_mul(mpqtemp, conshdlrdata->interiorpt[i], lambda2);
         else
            mpq_mul(mpqtemp, conshdlrdata->interiorray[i], lambda2);
         mpq_add(approxdualsol[i], approxdualsol[i], mpqtemp);
      }
   }

   /* postprocess dual solution to reduce values when both sides of constraint used;
    * if y(lhs) and y(rhs) are both nonzero shift them such that one becomes zero
    * this will tighten the solution and improve the objective value, there is no way this can hurt
    */
   if( PSPOSTPROCESSDUALSOL )
   {
      /* y(lhs) and y(rhs) corresponding to primal constraints */
      for( i = 0; i < nconss; i++ )
      {
         /* find the min value of y(lhs) and y(rhs) */
         if( mpq_cmp(approxdualsol[i], approxdualsol[i+nconss]) > 0 )
            mpq_set(mpqtemp, approxdualsol[i+nconss]);
         else
            mpq_set(mpqtemp, approxdualsol[i]);

         /* shift if both are nonzero */
         if( mpq_sgn(mpqtemp) > 0 )
         {
            mpq_sub(approxdualsol[i], approxdualsol[i], mpqtemp);
            mpq_sub(approxdualsol[i+nconss], approxdualsol[i+nconss], mpqtemp);
         }
      }
      /* y(lhs) and y(rhs) corresponding to bound constraints */
      for( i = 0; i < nvars; i++ )
      {
         /* find the min value of y(lhs) and y(rhs) */
         if( mpq_cmp(approxdualsol[i + 2*nconss], approxdualsol[i + 2*nconss + nvars]) > 0 )
            mpq_set(mpqtemp, approxdualsol[i + 2*nconss + nvars]);
         else
            mpq_set(mpqtemp, approxdualsol[i + 2*nconss]);

         /* shift if both are nonzero */
         if( mpq_sgn(mpqtemp) > 0 )
         {
            mpq_sub(approxdualsol[i+2*nconss], approxdualsol[i+2*nconss], mpqtemp);
            mpq_sub(approxdualsol[i+2*nconss + nvars], approxdualsol[i+2*nconss + nvars], mpqtemp);
         }
      }
   }

#ifdef PS_OUT
   printf("projected and shifted dual solution (should be an exact dual feasible solution)\n");
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualsol[i]);
      printf("\n");
   }
#endif


#ifndef NDEBUG
   SCIPdebugMessage("debug test: verifying feasibility of dual solution:\n");

   /* calculate violation of equality constraints: subtract Ax to get violation b-Ax, subtract A(approxdualsol) */
   rval = 0;
   for( i = 0; i < nvars; i++ )
   {
      mpq_set(violation[i], consdata->obj[i]);
   }
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++ )
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualsol[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualsol[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualsol[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualsol[i+2*nconss+nvars]);
   }
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) )
      {
         SCIPdebugMessage("   dual solution incorrect, violates equalties\n");
         i = nvars;
         rval = 1;
      }
   }
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualsol[i]) < 0 )
      {
         SCIPdebugMessage("   dual solution incorrect, negative components\n");
         i = nextendedconss;
         rval = 1;
      }
   }
   if( !rval )
      SCIPdebugMessage("   dual solution verified\n");
   assert(!rval);
#endif

   /* compute dual bound for constructed exact dual solution */
   mpq_set_ui(dualbound, 0, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      mpq_mul(mpqtemp, approxdualsol[i], costvect[i]);
      mpq_add(dualbound, dualbound, mpqtemp);
   }
   mpq_set(*boundval, dualbound);

#ifdef PS_OUT
   printf("   common slack=%.20f (", mpq_get_d(conshdlrdata->commonslack));
   mpq_out_str(stdout, 10, conshdlrdata->commonslack);
   printf(")\n");

   printf("   max violation=%.20f (", mpq_get_d(maxv));
   mpq_out_str(stdout, 10, maxv);
   printf(")\n");

   printf("   lambda (use of interior point)=%.20f (", mpq_get_d(lambda2));
   mpq_out_str(stdout, 10, lambda2);
   printf(")\n");

   printf("   dual objective value=%.20f (", mpq_get_d(dualbound));
   mpq_out_str(stdout, 10, dualbound);
   printf(")\n");
#endif

 TERMINATE:
   /* free memory */
   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(correction[i]);
   SCIPfreeBufferArray(scip, &correction);

   for( i = 0; i < nvars; i++ )
      mpq_clear(violation[i]);
   SCIPfreeBufferArray(scip, &violation);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(costvect[i]);
   SCIPfreeBufferArray(scip, &costvect);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(approxdualsol[i]);
   SCIPfreeBufferArray(scip, &approxdualsol);

   mpq_clear(dualbound);
   mpq_clear(maxv);
   mpq_clear(lambda2);
   mpq_clear(lambda1);
   mpq_clear(mpqtemp2);
   mpq_clear(mpqtemp);

   return SCIP_OKAY;
}

/** proves node infeasibility by project and shift method */
static
SCIP_RETCODE PScorrectdualray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_Bool*            success             /**< was corrected dual ray cost improving? */
   )
{
   int i;
   int j;
   int rval;
   int nextendedconss;
   int nconss;
   int nvars;
   int currentrow;
   int isfeas;
   mpq_t* approxdualray;
   mpq_t* costvect; /**< dual cost vector for expanded problem*/
   mpq_t* violation;
   mpq_t* correction;
   mpq_t mpqtemp;
   mpq_t mpqtemp2;
   mpq_t lambda1;
   mpq_t lambda2;
   mpq_t maxv;
   mpq_t dualbound;
   SCIP_COL** cols;


   /* project-and-shift method with interior ray correction (approximate dual ray is corrected to be exactly feasible)
    * 1. projection step (to ensure that equalities are satisfied):
    *   - compute error in equalities: r=c-Ay^
    *   - backsolve system of equations to find correction of error: z with Dz=r
    *   - add corretion to approximate dual ray: bold(y)=y^+[z 0]
    *  2. shifing step (to ensure that dual ray is feasible):
    *   - take convex combination of projected approximate ray bold(y) with S-interior ray y*
    *  3. compute objective value of ray
    */

   /* constructpsdata should always be called first */
   assert(conshdlrdata->psdatacon);

   /* if ray is not constructed, exit */
   if( !conshdlrdata->pshasray )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   SCIPdebugMessage("calling PScorrectdualray()\n");

   mpq_init(mpqtemp);
   mpq_init(mpqtemp2);
   mpq_init(lambda1);
   mpq_init(lambda2);
   mpq_init(maxv);
   mpq_init(dualbound);

   nextendedconss = conshdlrdata->nextendedconss;
   nconss = consdata->nconss;
   nvars = consdata->nvars;

   /* process bound changes */
   processBoundchgs(scip, conshdlrdata, consdata);

   /* allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIP_CALL( SCIPallocBufferArray(scip, &approxdualray, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(approxdualray[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &costvect, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(costvect[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &violation, nvars) );
   for( i = 0; i < nvars; i++ )
      mpq_init(violation[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &correction, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(correction[i]);

   /* recover the objective coefs and approximate solution value of dual solution;
    * dual vars of lhs constraints (including -inf) and rhs constraints (including +inf),
    * dual vars of lb constraint (including -inf) and ub constraints (including +inf)
    */
   for( i = 0; i < nconss; i++ )
   {
      /* in case we want to prove infeasibility it might be that we were not able to compute a dual solution
       * with bound exceeding objective value; in this case dual solution is set to SCIP_INVALID
       */
      if( SCIProwGetDualsol(consdata->rows[i]) == SCIP_INVALID )
      {
         SCIPdebugMessage("  no valid unbounded approx dual sol given\n");
         *success = FALSE;
         goto TERMINATE;
      }
      mpq_set_d(mpqtemp, SCIProwGetDualfarkas(consdata->rows[i]));

      /* lhs constraint */
      if( mpq_sgn(mpqtemp) > 0 )
         mpq_set(approxdualray[i], mpqtemp);
      /* rhs constraint */
      else
         mpq_neg(approxdualray[i+nconss], mpqtemp);

      mpq_set(costvect[i], consdata->lhs[i]);
      mpq_neg(costvect[i+nconss], consdata->rhs[i]);
   }
   cols = SCIPgetLPCols(scip);
   for( i = 0; i < nvars; i++ )
   {
      mpq_set_d(mpqtemp, -SCIPgetColFarkasCoef(scip, cols[i]));
      /* lb constraint */
      if( mpq_sgn(mpqtemp) > 0 )
         mpq_set(approxdualray[i + 2*nconss], mpqtemp);
      /* ub constraint */
      else
         mpq_neg(approxdualray[i + 2*nconss + nvars], mpqtemp);

      mpq_set(costvect[i + 2*nconss], consdata->lbloc[i]);
      mpq_neg(costvect[i + 2*nconss + nvars], consdata->ubloc[i]);
   }

#ifdef PS_OUT
   printf("approximate dual ray:\n");

   mpq_set_ui(dualbound, 0, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualray[i]);
      printf("\n");

      mpq_mul(mpqtemp, approxdualray[i], costvect[i]);
      mpq_add(dualbound, dualbound,mpqtemp);
   }
   printf("   objective value=%.20f (", mpq_get_d(dualbound));
   mpq_out_str(stdout, 10, dualbound);
   printf(")\n");
#endif

   /* first, ensure nonnegativity of dual ray and fix artificial dual variables to zero */
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualray[i]) < 0)
         mpq_set_si(approxdualray[i], 0, 1);
      else if( isNegInfinity(conshdlrdata, costvect[i]) )
      {
         assert(!conshdlrdata->includedcons[i]);
         mpq_set_si(approxdualray[i], 0, 1);
      }
   }

   /* calculate violation of equality constraints r=0-A^ty */
   for( i = 0; i < nvars; i++ )
   {
      mpq_set_si(violation[i], 0, 1);
   }

   /* A^ty for y corresponding to primal constraints */
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++)
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualray[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualray[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }

   /* A^ty for y corresponding to bound constraints */
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualray[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualray[i+2*nconss+nvars]);
   }

   /* project solution */
#ifdef PS_OUT
   printf("violation of ray: \n");
   for( i = 0; i < nvars; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, violation[i]);
      printf("\n");
   }
#endif

   /* if there is no violation of the constraints, then skip the projection */
   isfeas = 1;
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) )
      {
         isfeas = 0;
         i = nvars;
      }
   }

   /* isfeas is equal to one only if approximate dual solution is already feasible for the dual */
   if( !isfeas )
   {
      /* compute [z] with Dz=r (D depends on conshdlrdata->psdualcolselection) */
      rval = RECTLUsolveSystem( conshdlrdata->rectfactor, nvars, nextendedconss, violation, correction);

#ifdef PS_OUT
      printf("correction of ray: \n");
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         printf("   i=%d: ", i);
         mpq_out_str(stdout, 10, correction[i]);
         printf(", position=%d\n", conshdlrdata->psbasis[i]);
      }
#endif

      /* project step: compute bold(y)=y^+[z 0]; correct only components corresponding to D (npsbasis=#of columns in D */
      for( i = 0; i < conshdlrdata->npsbasis; i++ )
      {
         mpq_add(approxdualray[conshdlrdata->psbasis[i]], approxdualray[conshdlrdata->psbasis[i]], correction[i]);
      }
   }

#ifdef PS_OUT
   printf("updated dual ray:\n");
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualray[i]);
      printf("\n");
   }
#endif

   /* in this case we are using an interior ray that can be added freely to the solution */
   mpq_set_si(lambda1, 1, 1);

   /* compute lambda values: compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
   mpq_set_ui(lambda1, 1, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualray[i]) < 0 && conshdlrdata->includedcons[i] )
      {
         mpq_div(mpqtemp, approxdualray[i], conshdlrdata->interiorray[i]);
         mpq_neg(mpqtemp, mpqtemp);
         if( mpq_cmp(lambda2, mpqtemp) < 0 ) /* changed to < */
            mpq_set(lambda2, mpqtemp);
      }
   }

   /* perform shift */
   if( mpq_sgn(lambda2) != 0 )
   {
      for( i = 0; i < nextendedconss; i++ )
      {
         mpq_mul(approxdualray[i], approxdualray[i], lambda1);
      }
      for( i =0; i < nextendedconss; i++ )
      {
         mpq_mul(mpqtemp, conshdlrdata->interiorray[i], lambda2);
         mpq_add(approxdualray[i], approxdualray[i], mpqtemp);
      }
   }

   /* postprocess dual solution to reduce values when both sides of constraint used;
    * if y(lhs) and y(rhs) are both nonzero shift them such that one becomes zero
    * this will tighten the solution and improve the objective value, there is no way this can hurt
    */
   if( PSPOSTPROCESSDUALSOL )
   {
      /* y(lhs) and y(rhs) corresponding to primal constraints */
      for( i = 0; i < nconss; i++ )
      {
         /* find the min value of y(lhs) and y(rhs) */
         if( mpq_cmp(approxdualray[i], approxdualray[i+nconss]) > 0 )
            mpq_set(mpqtemp, approxdualray[i+nconss]);
         else
            mpq_set(mpqtemp, approxdualray[i]);

         /* shift if both are nonzero */
         if( mpq_sgn(mpqtemp) > 0 )
         {
            mpq_sub(approxdualray[i], approxdualray[i], mpqtemp);
            mpq_sub(approxdualray[i+nconss], approxdualray[i+nconss], mpqtemp);
         }
      }
      /* y(lhs) and y(rhs) corresponding to bound constraints */
      for( i = 0; i < nvars; i++ )
      {
         /* find the min value of y(lhs) and y(rhs) */
         if( mpq_cmp(approxdualray[i + 2*nconss], approxdualray[i + 2*nconss + nvars]) > 0 )
            mpq_set(mpqtemp, approxdualray[i + 2*nconss + nvars]);
         else
            mpq_set(mpqtemp, approxdualray[i + 2*nconss]);

         /* shift if both are nonzero */
         if( mpq_sgn(mpqtemp) > 0 )
         {
            mpq_sub(approxdualray[i+2*nconss], approxdualray[i+2*nconss], mpqtemp);
            mpq_sub(approxdualray[i+2*nconss + nvars], approxdualray[i+2*nconss + nvars], mpqtemp);
         }
      }
   }

#ifdef PS_OUT
   printf("projected and shifted dual ray (should be an exact feasible dual ray)\n");
   for( i = 0; i < nextendedconss; i++ )
   {
      printf("   i=%d: ", i);
      mpq_out_str(stdout, 10, approxdualray[i]);
      printf("\n");
   }
#endif


#ifndef NDEBUG
   SCIPdebugMessage("debug test: verifying feasibility of dual ray\n");

   /* calculate violation of equality constraints: subtract Ax to get violation 0-Ax, subtract A(approxdualray) */
   rval = 0;
   for( i = 0; i < nvars; i++ )
   {
      mpq_set_si(violation[i], 0, 1);
   }
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++ )
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualray[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualray[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualray[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualray[i+2*nconss+nvars]);
   }
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) )
      {
         SCIPdebugMessage("   dual ray incorrect, violates equalties\n");
         i = nvars;
         rval = 1;
      }
   }
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualray[i])<0 )
      {
         SCIPdebugMessage("   dual ray incorrect, negative components\n");
         i = nextendedconss;
         rval = 1;
      }
   }
   if( !rval )
   {
      SCIPdebugMessage("   dual ray verified\n");
   }
   assert(!rval);
#endif

   /* compute dual bound for constructed exact dual solution */
   mpq_set_ui(dualbound, 0, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      mpq_mul(mpqtemp, approxdualray[i], costvect[i]);
      mpq_add(dualbound, dualbound, mpqtemp);
   }

   /* if the objective value of the corrected ray is positive we can prune node, otherwise not */
   if( mpq_sgn(dualbound) > 0 )
      *success = TRUE;
   else
      *success = FALSE;

#ifdef PS_OUT
   printf("   common slack=%.20f (", mpq_get_d(conshdlrdata->commonslack));
   mpq_out_str(stdout, 10, conshdlrdata->commonslack);
   printf(")\n");

   printf("   max violation=%.20f (", mpq_get_d(maxv));
   mpq_out_str(stdout, 10, maxv);
   printf(")\n");

   printf("   lambda (use of interior point)=%.20f (", mpq_get_d(lambda2));
   mpq_out_str(stdout, 10, lambda2);
   printf(")\n");

   printf("   dual objective value=%.20f (", mpq_get_d(dualbound));
   mpq_out_str(stdout, 10, dualbound);
   printf(")\n");
#endif

 TERMINATE:
   /* free memory */
   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(correction[i]);
   SCIPfreeBufferArray(scip, &correction);

   for( i = 0; i < nvars; i++ )
      mpq_clear(violation[i]);
   SCIPfreeBufferArray(scip, &violation);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(costvect[i]);
   SCIPfreeBufferArray(scip, &costvect);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(approxdualray[i]);
   SCIPfreeBufferArray(scip, &approxdualray);

   mpq_clear(dualbound);
   mpq_clear(maxv);
   mpq_clear(lambda2);
   mpq_clear(lambda1);
   mpq_clear(mpqtemp2);
   mpq_clear(mpqtemp);

   return SCIP_OKAY;
}


/*
 * local methods for primal bound shift dual bounding methods based on interval data and rational data
 */

/** calculates y*b + min{(c - y*A)*x | lb <= x <= ub} for given vectors y and c;
 *  the vector b is defined with b[i] = lhs[i] if y[i] >= 0, b[i] = rhs[i] if y[i] < 0
 *  Calculating this value in interval arithmetics gives a proved lower LP bound for the following reason (assuming,
 *  we have only left hand sides):
 *           min{cx       |  b <=  Ax, lb <= x <= ub}
 *   >=      min{cx       | yb <= yAx, lb <= x <= ub}   (restriction in minimum is relaxed)
 *   == yb + min{cx - yb  | yb <= yAx, lb <= x <= ub}   (added yb - yb == 0)
 *   >= yb + min{cx - yAx | yb <= yAx, lb <= x <= ub}   (because yAx >= yb inside minimum)
 *   >= yb + min{cx - yAx |            lb <= x <= ub}   (restriction in minimum is relaxed)
 *   This is a modified version of the function provedBound() from lp.c, we no longer assume the problem is a relaxation.
 *   interval arithmetic is used for every operation.  y is taken exactly as it is given but A,c,b are all used as intervals.
 */
static
SCIP_RETCODE provedBoundInterval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_Bool             usefarkas,          /**< use y = dual farkas and c = 0 instead of y = dual solution and c = obj? */
   SCIP_Real*            boundval            /**< value of dual bound */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL* rhsinter;
   SCIP_INTERVAL* xinter;
   SCIP_INTERVAL* ainter;
   SCIP_INTERVAL* atyinter;
   SCIP_INTERVAL* cinter;
   SCIP_INTERVAL intertemp;
   SCIP_INTERVAL ytb;
   SCIP_INTERVAL minprod;
   SCIP_ROW* row;
   SCIP_Real* y;
   SCIP_Real* ycol;
   int i;
   int j;

   /* process bound changes */
   processBoundchgs(scip, conshdlrdata, consdata);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &y, consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhsinter, consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ycol, consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ainter, consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &atyinter, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cinter, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xinter, consdata->nvars) );

   /* calculate y^Tb */
   SCIPintervalSet(&ytb, 0.0);
   /* create y, rhs and constant vector in interval arithmetic */
   for( j = 0; j < consdata->nconss; ++j )
   {
      row = consdata->rows[j];
      assert(row != NULL);

      /* create y vector in interval arithmetic, setting near zeros to zero */
      y[j] = (usefarkas ? SCIProwGetDualfarkas(row) : SCIProwGetDualsol(row));

      if( SCIPisInfinity(scip, y[j]) )
         y[j] = SCIPinfinity(scip);

      if( SCIPisInfinity(scip, -y[j]) )
         y[j] = -1.0 * SCIPinfinity(scip);

      /** @todo exiptodo: dual bounding improvement
       *  - should we also set nonzero values of y to zero if corresponding lhs/rhs is not finite (to improve dual bound)?
       *  - do such situations come up?
       */
      /* create rhs and constant vectors in interval arithmetic */
      if( SCIPisFeasPositive(scip, y[j]) )
      {
         SCIPintervalSetBounds(&rhsinter[j],
            mpqGetRealRelax(scip, consdata->lhs[j], GMP_RNDD), mpqGetRealRelax(scip, consdata->lhs[j], GMP_RNDU));
      }
      else if( SCIPisFeasNegative(scip, y[j]) )
      {
         SCIPintervalSetBounds(&rhsinter[j],
            mpqGetRealRelax(scip, consdata->rhs[j], GMP_RNDD), mpqGetRealRelax(scip, consdata->rhs[j], GMP_RNDU));
      }
      else
      {
         y[j] = 0.0;
         SCIPintervalSet(&rhsinter[j], 0.0);
      }
   }
   SCIPintervalScalprodScalars(SCIPinfinity(scip), &ytb, consdata->nconss, rhsinter, y);

#ifndef NDEBUG
   for( j = 0; j < consdata->nconss; ++j )
   {
      row = consdata->rows[j];
      assert(row != NULL);

      if( !SCIPisFeasPositive(scip, y[j]) && !SCIPisFeasNegative(scip, y[j]) )
      {
         assert(rhsinter[j].inf == 0.0);
         assert(rhsinter[j].sup == 0.0);
      }
   }
#endif

   /* calculate min{(c^T - y^TA)x} */

   for( j = 0; j < consdata->nvars; j++)
   {
      SCIPintervalSet(&atyinter[j], 0.0);
   }

#ifdef STORE_INTERA
   /* the updated code will store the interval representation of the val matrix instead of recomputing it every time here.
    * instead of using "intertemp" we will use an array conshdlrdata->intervalval
    */
   if( !conshdlrdata->intervalvalcon )
   {

      SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->intervalval, consdata->nnonz) );
      for( j = 0; j < consdata->nconss; ++j )
      {
         for( i = consdata->beg[j]; i < consdata->beg[j] + consdata->len[j]; ++i )
         {
            SCIPintervalSetBounds(&conshdlrdata->intervalval[i],
               mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD), mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU));
         }
      }
      conshdlrdata->intervalvalcon = TRUE;
   }

   /* compute infimums of -A^Ty */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   for( j = 0; j < consdata->nconss; ++j )
   {
      for( i = consdata->beg[j]; i < consdata->beg[j] + consdata->len[j]; ++i )
      {
         SCIPintervalMulScalarInf(SCIPinfinity(scip), &intertemp, conshdlrdata->intervalval[i], -1.0 * y[j]);
         SCIPintervalAddInf(SCIPinfinity(scip), &atyinter[consdata->ind[i]], atyinter[consdata->ind[i]], intertemp);
      }

   }
   /* compute supremums of -A^Ty */
   SCIPintervalSetRoundingModeUpwards();
   for( j = 0; j < consdata->nconss; ++j )
   {
      for( i = consdata->beg[j]; i < consdata->beg[j] + consdata->len[j]; ++i )
      {
         SCIPintervalMulScalarSup(SCIPinfinity(scip), &intertemp, conshdlrdata->intervalval[i], -1.0 * y[j]);
         SCIPintervalAddSup(SCIPinfinity(scip), &atyinter[consdata->ind[i]], atyinter[consdata->ind[i]], intertemp);
      }
   }
   SCIPintervalSetRoundingMode(roundmode);

#else
   /* compute infimums of -A^Ty */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   for( j = 0; j < consdata->nconss; ++j )
   {
      for( i = consdata->beg[j]; i < consdata->beg[j] + consdata->len[j]; ++i )
      {
         SCIPintervalSetBounds(&intertemp,
            mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD), mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU));
         SCIPintervalMulScalarInf(SCIPinfinity(scip), &intertemp, intertemp, -1.0 * y[j]);
         SCIPintervalAddInf(SCIPinfinity(scip), &atyinter[consdata->ind[i]], atyinter[consdata->ind[i]], intertemp);
      }
   }
   /* compute supremums of -A^Ty */
   SCIPintervalSetRoundingModeUpwards();
   for( j = 0; j < consdata->nconss; ++j )
   {
      for( i = consdata->beg[j]; i < consdata->beg[j] + consdata->len[j]; ++i )
      {
         SCIPintervalSetBounds(&intertemp,
            mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD), mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU));
         SCIPintervalMulScalarSup(SCIPinfinity(scip), &intertemp, intertemp, -1.0 * y[j]);
         SCIPintervalAddSup(SCIPinfinity(scip), &atyinter[consdata->ind[i]], atyinter[consdata->ind[i]], intertemp);
      }
   }
   SCIPintervalSetRoundingMode(roundmode);
#endif

   /* create c vector and x vector in interval arithmetic and compute min{(c^T - y^TA)x} */
   for( j = 0; j < consdata->nvars; ++j )
   {
      if( usefarkas )
      {
         SCIPintervalSet(&cinter[j], 0.0);
      }
      else
      {
         SCIPintervalSetBounds(&cinter[j],
            mpqGetRealRelax(scip, consdata->obj[j], GMP_RNDD), mpqGetRealRelax(scip, consdata->obj[j], GMP_RNDU));
      }

      SCIPintervalSetBounds(&xinter[j],
            mpqGetRealRelax(scip, consdata->lbloc[j], GMP_RNDD), mpqGetRealRelax(scip, consdata->ubloc[j], GMP_RNDU));
   }
   SCIPintervalAddVectors(SCIPinfinity(scip), atyinter, consdata->nvars, atyinter, cinter);
   SCIPintervalScalprod(SCIPinfinity(scip), &minprod, consdata->nvars, atyinter, xinter);

   /* add y^Tb */
   SCIPintervalAdd(SCIPinfinity(scip), &minprod, minprod, ytb);

   /* free buffer for storing y in interval arithmetic */
   SCIPfreeBufferArray(scip, &xinter);
   SCIPfreeBufferArray(scip, &cinter);
   SCIPfreeBufferArray(scip, &atyinter);
   SCIPfreeBufferArray(scip, &ainter);
   SCIPfreeBufferArray(scip, &ycol);
   SCIPfreeBufferArray(scip, &rhsinter);
   SCIPfreeBufferArray(scip, &y);

   *boundval = SCIPintervalGetInf(minprod);

   return SCIP_OKAY;

}

/** compute safe dual bound by Neumaier and Shcherbina Bound using exact rational arithmetic.
 * Considering the primal of the form:
 * min c'x
 * lhs <= Ax <= rhs
 *  lb <=  x <= ub
 *
 * and the dual of the form
 *
 * max [lhs',-rhs',lb',-ub'] y
 *     [  A',  -A',  I,  -I] y =  c
 *                           y >= 0
 *
 * we will take an approximate dual solution y~ and use exact arithmetic to compute its
 * error r = c - [  A',  -A',  I,  -I] y~ and then use exact arithmetic to compute a valid
 * bound by increasing the bound variables as needed to correct this error.
 */
static
SCIP_RETCODE provedBoundRational(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_Bool             usefarkas,          /**< use y = dual farkas and c = 0 instead of y = dual solution and c = obj? */
   mpq_t*                boundval            /**< value of dual bound */
   )
{
   int i;
   int j;
#ifndef NDEBUG
   int rval;
#endif
   int nextendedconss;
   int nconss;
   int nvars;
   int currentrow;
   mpq_t* approxdualsol;
   mpq_t* costvect;
   mpq_t* violation;
   mpq_t mpqtemp;
   mpq_t dualbound;

   mpq_init(mpqtemp);
   mpq_init(dualbound);

   nconss = consdata->nconss;
   nvars = consdata->nvars;
   nextendedconss = 2*nconss + 2*nvars;

   /* process bound changes */
   processBoundchgs(scip, conshdlrdata, consdata);

   /*allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIP_CALL( SCIPallocBufferArray(scip, &approxdualsol, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(approxdualsol[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &costvect, nextendedconss) );
   for( i = 0; i < nextendedconss; i++ )
      mpq_init(costvect[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &violation, nvars) );
   for( i = 0; i < nvars; i++ )
      mpq_init(violation[i]);

   /* recover the objective coefs and approximate solution value of dual solution;
    * dual vars of lhs constraints and rhs constraints,
    * dual vars of lb constraint and ub constraints
    */
   for( i = 0; i < nconss; i++ )
   {
      if( usefarkas )
         mpq_set_d(mpqtemp, SCIProwGetDualfarkas(consdata->rows[i]));
      else
         mpq_set_d(mpqtemp, SCIProwGetDualsol(consdata->rows[i]));

      /* lhs constraint */
      if( mpq_sgn(mpqtemp) > 0 )
         mpq_set(approxdualsol[i], mpqtemp);
      /* rhs constraint */
      else
         mpq_neg(approxdualsol[i+nconss], mpqtemp);

      mpq_set(costvect[i], consdata->lhs[i]);
      mpq_neg(costvect[i+nconss], consdata->rhs[i]);
   }

   for( i = 0; i < nvars; i++ )
   {
      mpq_set(costvect[i + 2*nconss], consdata->lbloc[i]);
      mpq_neg(costvect[i + 2*nconss + nvars], consdata->ubloc[i]);
   }

   /* make sure we set components equal to zero if they have infinite rhs/lhs cost */
   for( i = 0; i < 2*nconss; i ++)
   {
      if( isNegInfinity(conshdlrdata, costvect[i]) )
      {
         mpq_set_si(approxdualsol[i], 0, 1);
      }
   }

   /* first, ensure nonnegativity of dual solution */
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualsol[i]) < 0)
         mpq_set_si(approxdualsol[i], 0, 1);
   }

   /* calculate violation of equality constraints r=c-A^ty */
   for( i = 0; i < nvars; i++ )
   {
      if( usefarkas )
         mpq_set_ui(violation[i], 0, 1);
      else
         mpq_set(violation[i], consdata->obj[i]);
   }

   /* A^ty for y corresponding to primal constraints */
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++)
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualsol[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualsol[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }
   /* A^ty for y corresponding to bound constraints */
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualsol[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualsol[i+2*nconss+nvars]);
   }


   /* correct the solution to be dual feasible by increasing the primal bound dual variables */
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) > 0 )
      {
         mpq_add(approxdualsol[i+2*nconss], approxdualsol[i+2*nconss], violation[i]);
      }
      else if( mpq_sgn(violation[i]) < 0 )
      {
         mpq_sub(approxdualsol[i+2*nconss+nvars], approxdualsol[i+2*nconss+nvars], violation[i]);
      }
   }

#ifndef NDEBUG
   SCIPdebugMessage("debug test: verifying feasibility of dual solution\n");
   /* calculate violation of equality constraints: subtract Ax to get violation b-Ax, subtract A(approxdualsol) */
   rval = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( usefarkas )
         mpq_set_ui(violation[i],0,1);
      else
         mpq_set(violation[i], consdata->obj[i]);
   }
   for( i = 0; i < nconss; i++ )
   {
      for( j = consdata->beg[i]; j < consdata->beg[i] + consdata->len[i]; j++ )
      {
         currentrow = consdata->ind[j];
         mpq_mul(mpqtemp, approxdualsol[i], consdata->val[j]);
         mpq_sub(violation[currentrow], violation[currentrow], mpqtemp);
         mpq_mul(mpqtemp, approxdualsol[i+nconss], consdata->val[j]);
         mpq_add(violation[currentrow], violation[currentrow], mpqtemp);
      }
   }
   for( i = 0; i < nvars; i++ )
   {
      mpq_sub(violation[i], violation[i], approxdualsol[i+2*nconss]);
      mpq_add(violation[i], violation[i], approxdualsol[i+2*nconss+nvars]);
   }
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(violation[i]) )
      {
         SCIPdebugMessage("   dual solution incorrect, violates equalties\n");
         i = nvars;
         rval = 1;
      }
   }
   for( i = 0; i < nextendedconss; i++ )
   {
      if( mpq_sgn(approxdualsol[i])<0 )
      {
         SCIPdebugMessage("   dual solution incorrect, negative components\n");
         i = nextendedconss;
         rval = 1;
      }
   }
   if( !rval )
      SCIPdebugMessage("   dual solution verified\n");
   assert(!rval);
#endif

   /* compute dual bound for constructed exact dual solution */
   mpq_set_ui(dualbound, 0, 1);
   for( i = 0; i < nextendedconss; i++ )
   {
      mpq_mul(mpqtemp, approxdualsol[i], costvect[i]);
      mpq_add(dualbound, dualbound, mpqtemp);
   }
   mpq_set(*boundval, dualbound);

   /* free memory */
   for( i = 0; i < nvars; i++ )
      mpq_clear(violation[i]);
   SCIPfreeBufferArray(scip, &violation);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(costvect[i]);
   SCIPfreeBufferArray(scip, &costvect);

   for( i = 0; i < nextendedconss; i++ )
      mpq_clear(approxdualsol[i]);
   SCIPfreeBufferArray(scip, &approxdualsol);

   mpq_clear(dualbound);
   mpq_clear(mpqtemp);

   return SCIP_OKAY;
}


/*
 * local methods for presolving
 */

/** removes rounding locks of exactlp constraint for all variables in the given linear constraint */
static
void unlockRoundingSingleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   int                   c                   /**< index of constraint in exactlp data structures */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int j;

   /* get problem variables */
   vars = SCIPgetVars(scip);
   assert(SCIPgetNVars(scip) == consdata->nvars);

   haslhs = !isNegInfinity(conshdlrdata, consdata->lhs[c]);
   hasrhs = !isPosInfinity(conshdlrdata, consdata->rhs[c]);

   /* check all variables of the current constraint */
   for( j = consdata->beg[c]; j < consdata->beg[c] + consdata->len[c]; ++j )
   {
      assert(consdata->ind[j] >= 0 && consdata->ind[j] < consdata->nvars);
      assert(consdata->lockup[consdata->ind[j]] >= 0 && consdata->lockdown[consdata->ind[j]] >= 0);

      assert(SCIPvarGetProbindex(vars[consdata->ind[j]]) == consdata->ind[j]);


      /* val > 0 */
      if( mpq_sgn(consdata->val[j]) > 0 )
      {
         if( haslhs )
            consdata->lockdown[consdata->ind[j]]--;

         if( hasrhs )
            consdata->lockup[consdata->ind[j]]--;

         SCIPunlockVarCons(scip, vars[consdata->ind[j]], cons,
            haslhs && consdata->lockdown[consdata->ind[j]] == 0, hasrhs && consdata->lockup[consdata->ind[j]] == 0);
      }

      /* val < 0 */
      if( mpq_sgn(consdata->val[j]) < 0 )
      {
         if( haslhs )
            consdata->lockup[consdata->ind[j]]--;

         if( hasrhs )
            consdata->lockdown[consdata->ind[j]]--;

         SCIPunlockVarCons(scip, vars[consdata->ind[j]], cons,
            hasrhs && consdata->lockdown[consdata->ind[j]] == 0, haslhs && consdata->lockup[consdata->ind[j]] == 0);
      }
      assert(consdata->lockup[consdata->ind[j]] >= 0 && consdata->lockdown[consdata->ind[j]] >= 0);
   }
}

/** removes linear constraint from exactlp constraint data */
static
void delSingleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   int                   c                   /**< index of constraint in exactlp data structures */
   )
{
   int i;

   assert(consdata->len[c] > 0);
   assert(consdata->rows == NULL); /** @todo exiptodo: otherwise the LP/LPEX rows have to be deleted as well */
   assert(!conshdlrdata->lpexconstructed);

   unlockRoundingSingleCons(scip, cons, conshdlrdata, consdata, c);

   for( i = c; i < consdata->nconss - 1; ++i )
   {
      consdata->beg[i] = consdata->beg[i+1];
      consdata->len[i] = consdata->len[i+1];
      mpq_set(consdata->lhs[i], consdata->lhs[i+1]);
      mpq_set(consdata->rhs[i], consdata->rhs[i+1]);
   }
   assert(i == consdata->nconss - 1);

   /* move dummy entry of beg array */
   consdata->beg[i] = consdata->beg[i+1];

   /* update number of linear constraint in exactlp constraint */
   consdata->nconss--;
}

/** gets activity bounds for linear constraint of exactlp constraint */
static
void singleConsGetActivityBounds(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   int                   c,                  /**< index of constraint in exactlp data structures */
   mpq_t                 minactivity,        /**< pointer to store the minimal activity */
   mpq_t                 maxactivity         /**< pointer to store the maximal activity */
   )
{
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   /** @todo exiptodo: presolving extension
    *  - currently we only work on pure bound constraints
    *  - implement more general presolving step which considers constraints with more than one variable.
    */
   if( consdata->len[c] > 1 )
   {
      mpq_set(minactivity, *negInfinity(conshdlrdata));
      mpq_set(maxactivity, *posInfinity(conshdlrdata));
   }
   else
   {
      int varind;
      int v;

      v = consdata->beg[c];
      varind = consdata->ind[v];

      /* val > 0 */
      if( mpq_sgn(consdata->val[v]) > 0 )
      {
         /* maxactivity = val*ub */
         if( isPosInfinity(conshdlrdata, consdata->ubloc[varind]) )
            mpq_set(maxactivity, *posInfinity(conshdlrdata));
         else
            mpq_mul(maxactivity, consdata->val[v], consdata->ubloc[varind]);

         /* minactivity = val*lb */
         if( isNegInfinity(conshdlrdata, consdata->lbloc[varind]) )
            mpq_set(minactivity, *negInfinity(conshdlrdata));
         else
            mpq_mul(minactivity, consdata->val[v], consdata->lbloc[varind]);
      }
      /* val < 0 */
      else
      {
         assert(mpq_sgn(consdata->val[v]) < 0);

         /* maxactivity = val*lb */
         if( isNegInfinity(conshdlrdata, consdata->lbloc[varind]) )
            mpq_set(maxactivity, *posInfinity(conshdlrdata));
         else
            mpq_mul(maxactivity, consdata->val[v], consdata->lbloc[varind]);

         /* minactivity = val*ub */
         if( isPosInfinity(conshdlrdata, consdata->ubloc[varind]) )
            mpq_set(minactivity, *negInfinity(conshdlrdata));
         else
            mpq_mul(minactivity, consdata->val[v], consdata->ubloc[varind]);
      }
   }
}

/** returns adjusted lower bound value, which is rounded for integral variable types */
static
void adjustedLb(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   const mpq_t           lb,                 /**< lower bound to adjust */
   mpq_t                 newlb               /**< pointer to store adjusted lower bound */
   )
{
   if( isNegInfinity(conshdlrdata, lb) )
      mpq_set(newlb, *negInfinity(conshdlrdata));
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      mpqCeil(newlb, lb);
   else
      mpq_set(newlb, lb);
}

/** returns adjusted upper bound value, which is rounded for integral variable types */
static
void adjustedUb(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   const mpq_t           ub,                 /**< upper bound to adjust */
   mpq_t                 newub               /**< pointer to store adjusted upper bound */
   )
{
   if( isPosInfinity(conshdlrdata, ub) )
      mpq_set(newub, *posInfinity(conshdlrdata));
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      mpqFloor(newub, ub);
   else
      mpq_set(newub, ub);
}

/** adjust lower bound to integral value, if variable is integral */
static
void varAdjustLb(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_VAR*             var,                /**< problem variable */
   const mpq_t           bound,              /**< lower bound to adjust */
   mpq_t                 newbound            /**< pointer to store adjusted lower bound */
   )
{
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
#endif
   assert(var != NULL);

#ifdef SCIP_DEBUG
   gmp_snprintf(s, SCIP_MAXSTRLEN, "adjust lower bound %Qd of <%s>\n", bound, SCIPvarGetName(var));
   SCIPdebugMessage(s);
#endif

   adjustedLb(conshdlrdata, SCIPvarGetType(var), bound, newbound);
}

/** adjust upper bound to integral value, if variable is integral */
static
void varAdjustUb(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_VAR*             var,                /**< problem variable */
   const mpq_t           bound,              /**< upper bound to adjust */
   mpq_t                 newbound            /**< pointer to store adjusted upper bound */
   )
{
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
#endif
   assert(var != NULL);

#ifdef SCIP_DEBUG
   gmp_snprintf(s, SCIP_MAXSTRLEN, "adjust upper bound %Qd of <%s>\n", bound, SCIPvarGetName(var));
   SCIPdebugMessage(s);
#endif

   adjustedUb(conshdlrdata, SCIPvarGetType(var), bound, newbound);
}

/** changes local exact lower bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;
   int probidx;

   /** @todo exiptodo: presolving extension
    *  - changing local bounds in presolving is currently only supported at the root node. Because for switching between
    *    different subtrees, keeping track of the bound changes would be required which is not supported yet.
    *    Implement this.
    */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustLb(conshdlrdata, var, newbound, adjustedbound);

   /* change the local bound and update pseudo objective value */
   probidx = consdata->ind[v];
   SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
         consdata->obj[probidx], adjustedbound, consdata->ubloc[probidx]) );
   mpq_set(consdata->lbloc[probidx], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[probidx], consdata->lb[probidx]) >= 0);
   assert(mpq_cmp(consdata->ubloc[probidx], consdata->ub[probidx]) <= 0);

#ifndef NDEBUG
   /* check whether pseudo objective value was updated correctly */
   checkPseudoobjval(scip, conshdlrdata, consdata);
#endif

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}

/** changes local exact upper bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;
   int probidx;

   /** @todo exiptodo: presolving extension
    *  - changing local bounds in presolving is currently only supported at the root node. Because for switching between
    *    different subtrees, keeping track of the bound changes would be required which is not supported yet.
    *    Implement this.
    */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustUb(conshdlrdata, var, newbound, adjustedbound);

   /* change the local bound and update pseudo objective value */
   probidx = consdata->ind[v];
   SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
         consdata->obj[probidx], consdata->lbloc[probidx], adjustedbound) );
   mpq_set(consdata->ubloc[probidx], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[probidx], consdata->lb[probidx]) >= 0);
   assert(mpq_cmp(consdata->ubloc[probidx], consdata->ub[probidx]) <= 0);

#ifndef NDEBUG
   /* check whether pseudo objective value was updated correctly */
   checkPseudoobjval(scip, conshdlrdata, consdata);
#endif

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}

/** changes global exact lower bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustLb(conshdlrdata, var, newbound, adjustedbound);

   /* change the local bound */
   if( mpq_cmp(adjustedbound, consdata->lbloc[consdata->ind[v]]) > 0 )
   {
      SCIP_CALL( chgVarLbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }

   /* change the global bound */
   mpq_set(consdata->lb[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}

/** changes global exact upper bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustUb(conshdlrdata, var, newbound, adjustedbound);

   /* change the local bound */
   if( mpq_cmp(adjustedbound, consdata->ubloc[consdata->ind[v]]) < 0 )
   {
      SCIP_CALL( chgVarUbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }

   /* change the global bound */
   mpq_set(consdata->ub[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}


/** changes exact lower bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustLb(conshdlrdata, var, newbound, adjustedbound);

   /* change the bound */
   if( SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( chgVarLbGlobal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }
   else
   {
      SCIP_CALL( chgVarLbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}

/** changes exact upper bound of variable; if possible, adjusts bound to integral value */
static
SCIP_RETCODE chgVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustUb(conshdlrdata, var, newbound, adjustedbound);

   /* change the bound */
   if( SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( chgVarUbGlobal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }
   else
   {
      SCIP_CALL( chgVarUbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound) );
   }

   mpq_clear(adjustedbound);

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  than the current bound; if possible, adjusts bound to integral value
 */
static
SCIP_RETCODE tightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound,           /**< new value for bound */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   mpq_t adjustedbound;

   assert(infeasible != NULL);

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustLb(conshdlrdata, var, newbound, adjustedbound);

   if( mpq_cmp(adjustedbound, consdata->ubloc[consdata->ind[v]]) > 0 )
   {
      *infeasible = TRUE;
      goto TERMINATE;
   }

   if( mpq_cmp(adjustedbound, consdata->lbloc[consdata->ind[v]]) <= 0 )
   {
      goto TERMINATE;
   }

#ifdef SCIP_DEBUG
   {
      char s[SCIP_MAXSTRLEN];
      gmp_snprintf(s, SCIP_MAXSTRLEN, "tighten lb of var<%s>: <%Qd> --> <%Qd>\n", SCIPvarGetName(var),
         consdata->lbloc[consdata->ind[v]], adjustedbound);
      SCIPdebugMessage(s);
   }
#endif

   SCIP_CALL( chgVarLb(scip, conshdlrdata, consdata, var, v, adjustedbound) );

   if( tightened != NULL )
      *tightened = TRUE;

 TERMINATE:
   mpq_clear(adjustedbound);
   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  than the current bound; if possible, adjusts bound to integral value
 */
static
SCIP_RETCODE tightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound,           /**< new value for bound */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   mpq_t adjustedbound;

   assert(infeasible != NULL);

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustUb(conshdlrdata, var, newbound, adjustedbound);

   if( mpq_cmp(adjustedbound, consdata->lbloc[consdata->ind[v]]) < 0 )
   {
      *infeasible = TRUE;
      goto TERMINATE;
   }

   if( mpq_cmp(adjustedbound, consdata->ubloc[consdata->ind[v]]) >= 0 )
   {
      goto TERMINATE;
   }

#ifdef SCIP_DEBUG
   {
      char s[SCIP_MAXSTRLEN];
      gmp_snprintf(s, SCIP_MAXSTRLEN, "tighten ub of var<%s>: <%Qd> --> <%Qd>\n", SCIPvarGetName(var),
         consdata->ubloc[consdata->ind[v]], adjustedbound);
      SCIPdebugMessage(s);
   }
#endif

   SCIP_CALL( chgVarUb(scip, conshdlrdata, consdata, var, v, adjustedbound) );

   if( tightened != NULL )
      *tightened = TRUE;

 TERMINATE:
   mpq_clear(adjustedbound);
   return SCIP_OKAY;
}

/** tightens bounds of a single variable due to activity bounds */
static
SCIP_RETCODE tightenVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   c,                  /**< index of constraint in exactlp data structures */
   int                   v,                  /**< index of variable in ind and val array */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_VAR* var;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   mpq_t minresactivity;
   mpq_t maxresactivity;
   mpq_t newlb;
   mpq_t newub;
   SCIP_Real newlbrelax;
   SCIP_Real newubrelax;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   assert(consdata->len[c] > 0);
   assert(consdata->beg[c] <= v && v < consdata->beg[c] + consdata->len[c]);
   assert(consdata->ind[v] >= 0 && consdata->ind[v] < consdata->nvars);
   assert(SCIPvarGetProbindex(vars[consdata->ind[v]]) == consdata->ind[v]);

   var = vars[consdata->ind[v]];
   assert(var != NULL);
   assert(consdata->ind[v] == SCIPvarGetProbindex(var));

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   assert(mpq_sgn(consdata->val[v]) != 0);
   assert(!isPosInfinity(conshdlrdata, consdata->lhs[c]));
   assert(!isNegInfinity(conshdlrdata, consdata->rhs[c]));

   /** @todo exiptodo: presolving extension
    *  - changing local bounds in presolving is currently only supported at the root node. Because for switching between
    *    different subtrees, keeping track of the bound changes would be required which is not supported yet.
    *    Implement this.
    */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /** @todo exiptodo: presolving extension
    *  - currently we only work on pure bound constraints
    *  - implement more general presolving step which considers constraints with more than one variable.
    */
   if( consdata->len[c] > 1 )
      return SCIP_OKAY;

   /* initialize gmp data */
   mpq_init(minresactivity);
   mpq_init(maxresactivity);
   mpq_init(newlb);
   mpq_init(newub);

   /* get activity bounds for linear constraint after setting variable to zero */
   mpq_set_d(minresactivity, 0.0);
   mpq_set_d(maxresactivity, 0.0);

   /* val > 0 */
   if( mpq_sgn(consdata->val[v]) > 0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !isNegInfinity(conshdlrdata, minresactivity) && !isPosInfinity(conshdlrdata, consdata->rhs[c]) )
      {
         /* newub = (rhs - minresactivity)/val */
         mpq_set(newub, consdata->rhs[c]);
         mpq_sub(newub, newub, minresactivity);
         mpq_div(newub, newub, consdata->val[v]);

         if( mpq_cmp(newub, consdata->ubloc[consdata->ind[v]]) < 0 )
         {
            /* tighten exact bound of variable */
            SCIP_CALL( tightenVarUb(scip, conshdlrdata, consdata, var, v, newub, &infeasible, &tightened) );

            if( infeasible )
            {
               *cutoff = TRUE;
               goto TERMINATE;
            }

            if( tightened )
            {
#ifdef SCIP_DEBUG
               char s[SCIP_MAXSTRLEN];
               gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d>: tighten <%s>, new bds=[%Qd,%Qd]\n",
                  c, SCIPvarGetName(var), consdata->lbloc[consdata->ind[v]], consdata->ubloc[consdata->ind[v]]);
               SCIPdebugMessage(s);
#endif
               newubrelax = mpqGetRealRelax(scip, consdata->ubloc[consdata->ind[v]], GMP_RNDU);

               /* change relaxed bound of variable */
               if( SCIPisLT(scip, newubrelax, SCIPvarGetUbLocal(var)) )
               {
                  SCIP_CALL( SCIPchgVarUb(scip, var, newubrelax) );
                  assert(SCIPvarGetUbLocal(var) >= newubrelax);

                  (*nchgbds)++;
               }
            }
         }
      }

      if( !isPosInfinity(conshdlrdata, maxresactivity) && !isNegInfinity(conshdlrdata, consdata->lhs[c]) )
      {
         /* newlb = (lhs - maxresactivity)/val */
         mpq_set(newlb, consdata->lhs[c]);
         mpq_sub(newlb, newlb, maxresactivity);
         mpq_div(newlb, newlb, consdata->val[v]);

         if( mpq_cmp(newlb, consdata->lbloc[consdata->ind[v]]) > 0 )
         {
            /* tighten exact bound of variable */
            SCIP_CALL( tightenVarLb(scip, conshdlrdata, consdata, var, v, newlb, &infeasible, &tightened) );

            if( infeasible )
            {
               *cutoff = TRUE;
               goto TERMINATE;
            }

            if( tightened )
            {
#ifdef SCIP_DEBUG
               char s[SCIP_MAXSTRLEN];

               gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d>: tighten <%s>, new bds=[%Qd,%Qd]\n",
                  c, SCIPvarGetName(var), consdata->lbloc[consdata->ind[v]], consdata->ubloc[consdata->ind[v]]);
               SCIPdebugMessage(s);
#endif
               newlbrelax = mpqGetRealRelax(scip, consdata->lbloc[consdata->ind[v]], GMP_RNDD);

               /* change relaxed bound of variable */
               if( SCIPisGT(scip, newlbrelax, SCIPvarGetLbLocal(var)) )
               {
                  SCIP_CALL( SCIPchgVarLb(scip, var, newlbrelax) );
                  assert(SCIPvarGetLbLocal(var) <= newlbrelax);

                  (*nchgbds)++;
               }
            }
         }
      }
   }
   /* val < 0 */
   else
   {
      assert(mpq_sgn(consdata->val[v]) == -1);

      /* check, if we can tighten the variable's bounds */
      if( !isNegInfinity(conshdlrdata, minresactivity) && !isPosInfinity(conshdlrdata, consdata->rhs[c]) )
      {
         /* newlb = (rhs - minresactivity)/val */
         mpq_set(newlb, consdata->rhs[c]);
         mpq_sub(newlb, newlb, minresactivity);
         mpq_div(newlb, newlb, consdata->val[v]);

         if( mpq_cmp(newlb, consdata->lbloc[consdata->ind[v]]) > 0 )
         {
            /* tighten exact bound of variable */
            SCIP_CALL( tightenVarLb(scip, conshdlrdata, consdata, var, v, newlb, &infeasible, &tightened) );

            if( infeasible )
            {
               *cutoff = TRUE;
               goto TERMINATE;
            }

            if( tightened )
            {
#ifdef SCIP_DEBUG
               char s[SCIP_MAXSTRLEN];

               gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d>: tighten <%s>, new bds=[%Qd,%Qd]\n",
                  c, SCIPvarGetName(var), consdata->lbloc[consdata->ind[v]], consdata->ubloc[consdata->ind[v]]);
               SCIPdebugMessage(s);
#endif
               newlbrelax = mpqGetRealRelax(scip, consdata->lbloc[consdata->ind[v]], GMP_RNDD);

               /* change relaxed bound of variable */
               if( SCIPisGT(scip, newlbrelax, SCIPvarGetLbLocal(var)) )
               {
                  SCIP_CALL( SCIPchgVarLb(scip, var, newlbrelax) );
                  assert(SCIPvarGetLbLocal(var) <= newlbrelax);

                  (*nchgbds)++;
               }
            }
         }
      }

      if( !isPosInfinity(conshdlrdata, maxresactivity) && !isNegInfinity(conshdlrdata, consdata->lhs[c]) )
      {
         /* newub = (lhs - maxresactivity)/val */
         mpq_set(newub, consdata->lhs[c]);
         mpq_sub(newub, newub, maxresactivity);
         mpq_div(newub, newub, consdata->val[v]);

         if( mpq_cmp(newub, consdata->ubloc[consdata->ind[v]]) < 0 )
         {
            /* tighten exact bound of variable */
            SCIP_CALL( tightenVarUb(scip, conshdlrdata, consdata, var, v, newub, &infeasible, &tightened) );

            if( infeasible )
            {
               *cutoff = TRUE;
               goto TERMINATE;
            }

            if( tightened )
            {
#ifdef SCIP_DEBUG
               char s[SCIP_MAXSTRLEN];

               gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d>: tighten <%s>, new bds=[%Qd,%Qd]\n",
                  c, SCIPvarGetName(var), consdata->lbloc[consdata->ind[v]], consdata->ubloc[consdata->ind[v]]);
               SCIPdebugMessage(s);
#endif
               newubrelax = mpqGetRealRelax(scip, consdata->ubloc[consdata->ind[v]], GMP_RNDU);

               /* change relaxed bound of variable */
               if( SCIPisLT(scip, newubrelax, SCIPvarGetUbLocal(var)) )
               {
                  SCIP_CALL( SCIPchgVarUb(scip, var, newubrelax) );
                  assert(SCIPvarGetUbLocal(var) >= newubrelax);

                  (*nchgbds)++;
               }
            }
         }
      }
   }

 TERMINATE:
   /* free gmp data */
   mpq_clear(newub);
   mpq_clear(newlb);
   mpq_clear(maxresactivity);
   mpq_clear(minresactivity);

   return SCIP_OKAY;
}

/** tightens bounds of variables in linear constraint of exactlp constraint due to activity bounds */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   c,                  /**< index of constraint in exactlp data structures */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   int v;

   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   for( v = consdata->beg[c]; v < consdata->beg[c] + consdata->len[c]; ++v )
   {
      SCIP_CALL( tightenVarBounds(scip, conshdlrdata, consdata, vars, c, v, cutoff, nchgbds) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyExactlp NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExactlp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPprimalexCreate(&conshdlrdata->primal) );

   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->provedfeaslptime) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->provedinfeaslptime) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->exactfeaslptime) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->exactinfeaslptime) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->exactunsollptime) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExactlp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPprimalexFree(&conshdlrdata->primal, SCIPblkmem(scip)) );

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->provedfeaslptime) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->provedinfeaslptime) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->exactfeaslptime) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->exactinfeaslptime) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->exactunsollptime) );

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreExactlp)
{
   SCIP_CONSDATA* consdata;

   assert(nconss >= 0 && nconss <= 1);

   if( nconss == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   /* update objscale of transforemd problem;
    * note that this value is not reliable and should never be used to compute dual and primal bounds
    */
   SCIPsetTransObjscale(scip, mpqGetRealApprox(scip, consdata->objscale));

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreExactlp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(nconss >= 0 && nconss <= 1);

   if( nconss == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   /* check, wheter objective value is always integral by inspecting the problem */
   checkObjIntegral(scip, consdata);

   /* if possible, scale objective function such that it becomes integral with gcd 1 */
   SCIP_CALL( scaleObj(scip, conshdlrdata, consdata) );

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolExactlp NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolExactlp)
{  /*lint --e{715}*/
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* release and free the rows */
      SCIP_CALL( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMessage("ConsDelete method of exactlp constraints\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* free exactlp constraint */
   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   SCIP_Bool success;
   mpq_t* newobj;
   mpq_t intscalar;
   int i;

   SCIPdebugMessage("Trans method of exactlp constraints\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &newobj, sourcedata->nvars) );

   /* initialize objective coefficients for transformed constraint */
   for( i = 0; i < sourcedata->nvars; ++i )
   {
      mpq_init(newobj[i]);
      mpq_set(newobj[i], sourcedata->obj[i]);
   }
   mpq_init(intscalar);
   mpq_set_d(intscalar, 1.0);

   success = TRUE;

   /* in case of maximization, use negative objective coefficients in the transformed constraint */
   if( sourcedata->objsense == SCIP_OBJSENSE_MAXIMIZE )
   {
      mpq_t objsense;

      mpq_init(objsense);
      mpq_set_d(objsense, -1.0);

      /* update exact objective values */
      for( i = 0; i < sourcedata->nvars; ++i )
         mpq_mul(newobj[i], newobj[i], objsense);

      mpq_clear(objsense);
   }

   /* in case of working with an FP relaxation, all objective coefficients need to be FP representable;
    * if this is not the case yet, scale all objective coefficients to integral values
    */
   if( sourcedata->objneedscaling )
   {
      SCIP_VAR** origvars;
      SCIP_VAR** vars;

      assert(SCIPuseFPRelaxation(scip));

      /** @todo exiptodo: check which value is suitable for maxscale such that there is no overflow */
      SCIP_CALL( SCIPmpqCalcIntegralScalar((const mpq_t*) newobj, sourcedata->nvars, SCIPinfinity(scip), intscalar,
            &success) );

      if( success )
      {
         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, sourcedata->nvars) );

         /* get transformed problem variables */
         origvars = SCIPgetOrigVars(scip);
         assert(SCIPgetNOrigVars(scip) == sourcedata->nvars);
         SCIP_CALL( SCIPgetTransformedVars(scip, sourcedata->nvars, origvars, vars) );

         /* scale exact objective values s.t. all become FP representable and store them in the transformed problem */
         for( i = 0; i < sourcedata->nvars && success; ++i )
         {
#ifdef SCIP_DEBUG
	    char s[SCIP_MAXSTRLEN];
#endif
	    assert(SCIPvarIsOriginal(origvars[i]));
	    assert(SCIPvarIsTransformed(vars[i]));

#ifdef SCIP_DEBUG
	    gmp_snprintf(s, SCIP_MAXSTRLEN, " -> i=%d: old val=%Qd --> \n", i, newobj[i]);
	    SCIPdebugMessage(s);
#endif

	    mpq_mul(newobj[i], newobj[i], intscalar);

#ifdef SCIP_DEBUG
	    gmp_snprintf(s, SCIP_MAXSTRLEN, "                           scaled val=%Qd\n", newobj[i]);
	    SCIPdebugMessage(s);
#endif

	    SCIP_CALL( SCIPchgVarObj(scip, vars[i], mpqGetRealApprox(scip, newobj[i])) );

	    /* large integral values might not be FP representable */
	    if( !mpqIsReal(scip, newobj[i]) )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMessage("                           --> does not fit into a double\n");
#endif
               success = FALSE;
            }
         }

#ifndef NDEBUG
         /* check for correct status of original variables */
         checkOrigVars(scip);
#endif

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &vars);
      }
   }

   if( success )
   {
      /* create exactlp constraint data for target constraint */
      SCIP_CALL( consdataCreate(scip, conshdlrdata, &targetdata, conshdlrdata->eventhdlr, sourcedata->objsense,
            sourcedata->nvars, sourcedata->ninfbounds, sourcedata->ninfintbounds, sourcedata->nlargebounds, newobj,
            sourcedata->lb, sourcedata->ub, sourcedata->nconss, sourcedata->nsplitconss, sourcedata->conssize,
            sourcedata->lhs, sourcedata->rhs, sourcedata->nnonz, sourcedata->nintegral, sourcedata->beg, sourcedata->len,
            sourcedata->ind, sourcedata->val, sourcedata->minabsval, sourcedata->maxabsval, FALSE) );

      /* create target constraint */
      SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
            SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
            SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

      /* update objscale of transformed constraint */
      mpq_div(targetdata->objscale, targetdata->objscale, intscalar);

      /* check, wheter objective value is always integral by inspecting the problem */
      checkObjIntegral(scip, targetdata);
   }

   /* free temporary memory */
   mpq_clear(intscalar);

   for( i = 0; i < sourcedata->nvars; ++i )
      mpq_clear(newobj[i]);
   SCIPfreeBufferArray(scip, &newobj);

   if( !success )
   {
      SCIPerrorMessage("FP-relaxation not supported (could not scale objective funtion to FP-numbers) --> switch to FP-approximation via \"set misc usefprelax FALSE\"\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpExactlp)
{  /*lint --e{715}*/
   int i;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   for( i = 0; i < nconss; ++i )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i]) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool safedualboundcomputed;
   SCIP_Real safedualbound;
   char dualboundmethod;
   int ncolsex;
   int nrowsex;

#ifdef SCIP_DEBUG
   SCIP_Real oldlb;
   oldlb = SCIPgetLocalLowerbound(scip);
#endif

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("separating exactlp constraint <%s> on LP solution (LP solstat=%d)\n", SCIPconsGetName(conss[0]),
      SCIPgetLPSolstat(scip));

   /* update root node LP activities of the constraints for project and shift method */
   if( SCIPdualBoundMethod(scip) == 'p' || SCIPdualBoundMethod(scip) == 'a' )
   {
      SCIPdebugMessage("   recording the root LP activity for later use by project-and-shift\n");
      copyRootActivity( scip, conshdlrdata, consdata );
   }

   /** @todo exiptodo: separation extension
    *  - as soon as we actually add cutting planes here, we might want to compute a proved bound after the separation loop
    *    in solve.c. currently we disabled it inorder to avoid a second call of proved bound, the first one is done before
    *    separation in order to know whether Neumaier Shcherbina worked for the selection of the dual bounding method
    */

   /* select dual bounding method to apply */
   dualboundmethod = SCIPselectDualBoundMethod(scip, FALSE);

   /* in case the FP problem is a relaxation of the original problem and we use Neumaier and Shcherbinas
    * dual bounding method, we have already calculated a proved lower bound via postprocessing the LP solution
    * of the FP problem
    */
   if( SCIPuseFPRelaxation(scip) && (dualboundmethod == 'n' || dualboundmethod == 'a') )
      return SCIP_OKAY;

   /* dual bound will be calculated in enfops method, as we can not branch here */
   if( dualboundmethod == 'e' )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* update lower bound of current node wrt the pseudo objective value */
   if( !SCIPuseFPRelaxation(scip) )
   {
      SCIP_CALL( SCIPupdateLocalLowerbound(scip,
            mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) );

#ifdef SCIP_DEBUG
      if( oldlb < SCIPgetLocalLowerbound(scip) )
      {
         SCIPdebugMessage("   by pseudosol: local lower bound improved: %.20f --> %.20f\n", oldlb,
            SCIPgetLocalLowerbound(scip));
      }
      else
      {
         SCIPdebugMessage("   by pseudosol: local lower bound did not improve: %.20f -/-> %.20f\n", oldlb,
            SCIPgetLocalLowerbound(scip));
      }
      oldlb = SCIPgetLocalLowerbound(scip);
#endif
   }

   /* check for infeasible node by bounding */
   if( SCIPgetLocalLowerbound(scip) >= SCIPgetCutoffbound(scip) )
   {
      SCIPdebugMessage("   --> node will be cut off by bounding (lower=%g, upper=%g)\n", SCIPgetLocalLowerbound(scip),
         SCIPgetCutoffbound(scip));

      return SCIP_OKAY;
   }

   /* start timing */
   SCIPstartClock(scip, conshdlrdata->provedfeaslptime);
   safedualboundcomputed = FALSE;
   safedualbound = -SCIPinfinity(scip);

   switch( dualboundmethod )
   {
   case 'v':
      {
         SCIP_COL** cols;
         SCIP_ROW** rows;
         SCIP_Real* dualsol;
         SCIP_Real* primalsol;
         int ncols;
         int nrows;
         int j;

         /* get inexact LP data */
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
         SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

         /* constructs exact LP of current node */
         SCIP_CALL( constructCurrentLPEX(scip, conshdlrdata, consdata) );
         SCIP_CALL( SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncolsex) );
         SCIP_CALL( SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrowsex) );

         /* test whether LP state (i.e. basis information) of inexact LP is dual feasbile in rational arithmetic */
         if( ncolsex == ncols && nrowsex == nrows )
         {
            SCIP_LPISTATE* lpistate;
            SCIP_Bool dualfeasible;
            SCIP_Bool useprestep;
            mpq_t dualobjval;

            SCIP_CALL( SCIPallocBufferArray(scip, &primalsol, ncols + nrows) );
            SCIP_CALL( SCIPallocBufferArray(scip, &dualsol, nrows) );

            mpq_init(dualobjval);

            SCIP_CALL( SCIPgetLPState(scip, &lpistate) );

#ifdef TESTDBVERIFY
            {
               SCIP_Bool lperror;

               printf("\n\n................................. Test basis verification method ..........................\n");
               printf("1. solve LP exactly:\n");
               SCIP_CALL( loadLPState(scip, conshdlrdata) );
               SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, FALSE, &lperror) );
            }
#endif

            /* get approximate primal and dual solution of inexact LP solver */
            for( j = 0; j < ncols; ++j )
            {
               primalsol[j] = SCIPcolGetPrimsol(cols[j]);
               if( SCIPisFeasZero(scip, primalsol[j]) )
                  primalsol[j] = 0.0;
            }
            for( j = 0; j < nrows; ++j )
            {
               primalsol[ncols+j] = SCIPgetRowFeasibility(scip, rows[j]);
               if( SCIPisFeasZero(scip, primalsol[ncols+j]) )
                  primalsol[ncols+j] = 0.0;

               dualsol[j] = SCIProwGetDualsol(rows[j]);
               if( SCIPisFeasZero(scip, dualsol[j]) )
                  dualsol[j] = 0.0;
            }

#ifdef TESTDBVERIFY
            printf("2. apply basis verification:\n");
#endif
            /* verify dual feasibility of basis of inexact LP solver */
            useprestep = TRUE;
            SCIP_CALL( SCIPlpiexStateDualFeasible(conshdlrdata->lpiex, SCIPblkmem(scip), lpistate, useprestep,
                  primalsol, dualsol, &dualfeasible, &dualobjval) );

#ifdef TESTDBVERIFY
            printf("  --> LP basis %s dual feasible (%.50f)\n", dualfeasible ? "IS" : "IS NOT",
               mpqGetRealRelax(scip, dualobjval, GMP_RNDD));
            printf("3. by hand: Compare result of basis verify <-> status of exact LP solver with double prescion.\n");
            printf("...........................................................................................\n\n");
#endif

            /* free LP state */
            SCIP_CALL( SCIPfreeLPState(scip, &lpistate) );

            SCIPdebugMessage("   db method <v>: LP basis %s dual feasible (%.50f)\n", dualfeasible ? "is" : "is not",
               mpqGetRealRelax(scip, dualobjval, GMP_RNDD));

            if( dualfeasible )
            {
               SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) );

               safedualboundcomputed = TRUE;
               safedualbound = mpqGetRealRelax(scip, dualobjval, GMP_RNDD);
            }
            else
            {
               /* update number of fails */
               conshdlrdata->nfailprovedfeaslp++;
            }

            mpq_clear(dualobjval);

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &dualsol);
            SCIPfreeBufferArray(scip, &primalsol);
         }
         break;
      }

   case 'p':
      {
         SCIP_Bool success;
         mpq_t dualobjval;

         if( conshdlrdata->psdatafail )
         {
            /* update number of fails */
            conshdlrdata->nfailprovedfeaslp++;
            break;
         }

         SCIP_CALL( constructPSData(scip, conshdlrdata, consdata) );
         if( conshdlrdata->psdatafail )
         {
            /* update number of fails */
            conshdlrdata->nfailprovedfeaslp++;
            break;
         }
         mpq_init(dualobjval);

         SCIP_CALL( getPSdualbound(scip, conshdlrdata, consdata, &dualobjval, &success) );

         if( success )
         {
            SCIP_Real dualbound;

            dualbound = mpqGetRealRelax(scip, dualobjval, GMP_RNDD);

            /** @todo exiptodo: if this assert is not correct we have to change the way we update nfailprovedfeaslp */
            assert(!SCIPisInfinity(scip, -dualbound));

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualbound) );

            safedualboundcomputed = TRUE;
            safedualbound = dualbound;
         }
         else
         {
            /* update number of fails */
            conshdlrdata->nfailprovedfeaslp++;
         }
         mpq_clear(dualobjval);

         break;
      }
   case 'i':
      {
         SCIP_Real dualobjval;

         SCIP_CALL( provedBoundInterval(scip, conshdlrdata, consdata, FALSE, &dualobjval) );

         SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualobjval) );
         if( !SCIPisInfinity(scip, -dualobjval) )
         {
            safedualboundcomputed = TRUE;
            safedualbound = dualobjval;
         }
         else
         {
            /* update number of fails */
            conshdlrdata->nfailprovedfeaslp++;
         }
         break;
      }
   case 'x':
      {
         mpq_t dualobjval;
         SCIP_Real dualbound;

         mpq_init(dualobjval);
         SCIP_CALL( provedBoundRational(scip, conshdlrdata, consdata, FALSE, &dualobjval) );
         dualbound = mpqGetRealRelax(scip, dualobjval, GMP_RNDD);

         SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualbound) );
         if( !SCIPisInfinity(scip, -dualbound) )
         {
            safedualboundcomputed = TRUE;
            safedualbound = dualbound;
         }
         else
         {
            /* update number of fails */
            conshdlrdata->nfailprovedfeaslp++;
         }
         mpq_clear(dualobjval);

         break;
      }

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for dual bounding method\n", dualboundmethod);
      return SCIP_PARAMETERWRONGVAL;
   }

   /* stop timing and update number of calls */
   SCIPstopClock(scip, conshdlrdata->provedfeaslptime);
   conshdlrdata->nprovedfeaslp++;

   /* compute dual bound quality if dual bounding method did not fail */
   if( safedualboundcomputed )
   {
#ifdef SCIP_DEBUG
      if( oldlb < SCIPgetLocalLowerbound(scip) )
      {
         SCIPdebugMessage("   by dbmethod <%c>: local lower bound improved: %.20f -----> %.20f\n", dualboundmethod, oldlb,
            SCIPgetLocalLowerbound(scip));
      }
      else
      {
         SCIPdebugMessage("   by dbmethod <%c>: local lower bound did not improved: %.20f --/--> %.20f\n", dualboundmethod,
            oldlb, SCIPgetLocalLowerbound(scip));
      }
#endif

      SCIP_CALL( SCIPcomputeDualboundQuality(scip, conss[0], safedualbound) );;
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolExactlp NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool lperror;
   char dualboundmethod;
   char algo;


   assert(SCIPhasCurrentNodeLP(scip));
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIPdebugMessage("enforcing exactlp constraint <%s> on LP solution (LP solstat=%d)\n", SCIPconsGetName(conss[0]),
      SCIPgetLPSolstat(scip));

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_INFEASIBLE;

   /* select dual bounding method to apply */
   dualboundmethod = SCIPselectDualBoundMethod(scip, FALSE);
   assert(SCIPgetNLPBranchCands(scip) == 0 || dualboundmethod == 'e');

   /* start timing for dual bounding method 'e' (is performed here, as we cannot branch in sepalp) */
   if( dualboundmethod == 'e' )
      SCIPstartClock(scip, conshdlrdata->provedfeaslptime);
   else
      SCIPstartClock(scip, conshdlrdata->exactfeaslptime);

   /* constructs exact LP of current node */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   /* updates lower bound of current node wrt the pseudo objective value */
   if( !SCIPuseFPRelaxation(scip) )
   {
      SCIP_CALL( SCIPupdateLocalLowerbound(scip,
            mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) );
   }

   /* load LP state from inexact LP into exact LP solver */
   SCIP_CALL( loadLPState(scip, conshdlrdata) );

   /** @todo exiptodo:
    *  - which simplex algorithm sould be used here? selection by warmstartinfo is not useful; maybe use the same
    *    algorithm as in the inexact LP solver.
    */
   /* solve exact LP */
   algo = 'd';
   switch( algo )
   {
   case 'd':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, TRUE, &lperror) );
      break;

   case 'p':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_PRIMALSIMPLEX, TRUE, &lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for exact LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }

   /* stop timing and update number of calls for dual bounding method 'e'
    * (is performed here, as we cannot branch in sepalp)
    */
   if( dualboundmethod == 'e' )
   {
      SCIPstopClock(scip, conshdlrdata->provedfeaslptime);
      conshdlrdata->nprovedfeaslp++;
      if( !lperror )
      {
         mpq_t dualobjval;

         mpq_init(dualobjval);
         SCIP_CALL( SCIPlpiexGetObjval(conshdlrdata->lpiex, &dualobjval) );
         SCIP_CALL( SCIPcomputeDualboundQuality(scip, conss[0], mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) );
         mpq_clear(dualobjval);
      }
   }
   else
   {
      SCIPstopClock(scip, conshdlrdata->exactfeaslptime);
      conshdlrdata->nexactfeaslp++;
   }

   if( lperror )
   {
      SCIPerrorMessage("exact LP solver returns error: case not handled yet\n");

      if( dualboundmethod != 'e' || SCIPgetNLPBranchCands(scip) == 0 )
         conshdlrdata->nwrongexactfeaslp++;
      else
         conshdlrdata->nfailprovedfeaslp++;

      return SCIP_ERROR;
   }

   /* evaluate result of exact LP solver */
   SCIP_CALL( evaluateLPEX(scip, conss[0], conshdlrdata, consdata, result) );

   /* update number of wrong integral LP claims */
   if( *result == SCIP_BRANCHED && ( dualboundmethod != 'e' || SCIPgetNLPBranchCands(scip) == 0 ) )
      conshdlrdata->nwrongexactfeaslp++;

   /* update number of wrong integral LP claims, where even the LP was infeasible */
   if( *result == SCIP_CUTOFF && SCIPlpiexIsPrimalInfeasible(conshdlrdata->lpiex)
      && ( dualboundmethod != 'e' || SCIPgetNLPBranchCands(scip) == 0 ) )
      conshdlrdata->nwronginfeasexactfeaslp++;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool lperror;
   char dualboundmethod;
   char algo;

   SCIPdebugMessage("enforcing exactlp constraint <%s> on pseudo solution (LP solstat=%d)\n", SCIPconsGetName(conss[0]),
      SCIPgetLPSolstat(scip));

   assert(SCIPisExactSolve(scip));
   assert(!objinfeasible);
   assert(!SCIPhasCurrentNodeLP(scip));
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_INFEASIBLE;

   /** @todo exiptodo:
    *  - this is just a workaround to avoid stopping at the hard time limit. it is only needed because currently the exact
    *    LP solver QSopt_ex does not correctly handle time and iteration limits.
    */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_TIMELIMIT )
      return SCIP_OKAY;

   /* select dual bounding method to apply */
   dualboundmethod = SCIPselectDualBoundMethod(scip, TRUE);

   /* try to prove infeasibility */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE &&
      (!SCIPuseFPRelaxation(scip) || dualboundmethod != 'n') )
   {
      SCIPdebugMessage("   proving infeasibility by dbmethod <%c>\n", dualboundmethod);

      switch( dualboundmethod )
      {
      case 'p':
         {
            if( conshdlrdata->psinfeasversion )
            {
               mpq_t dualobjval;
               SCIP_Bool success;

               /* start timing */
               SCIPstartClock(scip, conshdlrdata->provedinfeaslptime);

               /* nothing can be done if requirements for project-and-shift method are not satisfied */
               if( conshdlrdata->psdatafail )
               {
                  /* update number of fails */
                  conshdlrdata->nfailprovedinfeaslp++;
                  break;
               }
               /* construct data needed for project-and-shift method if not done yet */
               SCIP_CALL( constructPSData(scip, conshdlrdata, consdata) );
               if( conshdlrdata->psdatafail )
               {
                  SCIPdebugMessage("      constructPSData failed\n");

                  /* update number of fails */
                  conshdlrdata->nfailprovedinfeaslp++;
                  break;
               }

               /* try to prove infeasibility by reparing solution with cutoff bound exceeding objective value */
               /* if we are trying to cut off node via primal bound, compute valid bound */
               /* call this if we are not using the ray for infeasibility OR the ray failed to construct */
               if( !SCIPpsInfeasRay(scip) || !conshdlrdata->pshasray  )
               {
                  mpq_init(dualobjval);
                  SCIP_CALL( getPSdualbound(scip, conshdlrdata, consdata, &dualobjval, &success) );
                  if( success )
                  {
#ifdef SCIP_DEBUG
                     SCIP_Real oldlb;
                     oldlb = SCIPgetLocalLowerbound(scip);
#endif

                     SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) );

#ifdef SCIP_DEBUG
                     if( oldlb < SCIPgetLocalLowerbound(scip) )
                     {
                        SCIPdebugMessage("      local lower bound improved: %.20f --> %.20f\n", oldlb,
                           SCIPgetLocalLowerbound(scip));
                     }
                     else
                     {
                        SCIPdebugMessage("      local lower bound did not improve: %.20f -/-> %.20f\n", oldlb,
                           SCIPgetLocalLowerbound(scip));
                     }
#endif

                     /* check if node can be cut off by bounding, i.e., whether construction was successful */
                     if( SCIPgetLocalLowerbound(scip) >= SCIPgetCutoffbound(scip) )
                     {

                        SCIPdebugMessage("      node can be cut off by bounding (lower=%g, upper=%g)\n",
                           SCIPgetLocalLowerbound(scip), SCIPgetCutoffbound(scip));
                        *result = SCIP_CUTOFF;
                     }
                     else
                     {
                        /* update number of fails */
                        conshdlrdata->nfailprovedinfeaslp++;
                     }
                  }
                  else
                  {
                     SCIPdebugMessage("      project-and-shift was not able to proof infeasibility\n");

                     /* update number of fails */
                     conshdlrdata->nfailprovedinfeaslp++;
                  }
                  mpq_clear(dualobjval);
               }
               else
               {
                  /* else, try to prove infeasibility of node by correcting a dual ray to be exactly feasible */

                  SCIP_CALL( PScorrectdualray(scip, conshdlrdata, consdata, &success) );
                  if( success )
                  {
                     SCIPdebugMessage("node is cut off by unbounded dual LP.  Exactly feasible cost-improving dual ray found.\n");
                     *result = SCIP_CUTOFF;
                  }
                  else
                  {
                     /* update number of fails */
                     SCIPdebugMessage("failed to correct dual ray \n");
                     conshdlrdata->nfailprovedinfeaslp++;
                  }
               }

               /* stop timing and update number of calls */
               SCIPstopClock(scip, conshdlrdata->provedinfeaslptime);
               conshdlrdata->nprovedinfeaslp++;
            }
            break;
         }
      case 'x':
         {
            mpq_t dualobjval;

            /* start timing */
            SCIPstartClock(scip, conshdlrdata->provedinfeaslptime);

            mpq_init(dualobjval);
            SCIP_CALL( provedBoundRational(scip, conshdlrdata, consdata, TRUE, &dualobjval) );

            /* check if node can be cut off because safe infeasibility proof was found */
            if( mpqGetRealRelax(scip, dualobjval, GMP_RNDD) > 0.0 )
            {
               SCIPdebugMessage("   --> node is cut off by bounding (farkasproof<%.20f> <= 0)\n",
                  mpqGetRealRelax(scip, dualobjval, GMP_RNDD));
               *result = SCIP_CUTOFF;
            }
            else
            {
               SCIPdebugMessage("   --> node NOT cut off by provedboundrational (farkasproof<%.20f> > 0)\n",
                  mpqGetRealRelax(scip, dualobjval, GMP_RNDD));

               /* update number of fails */
               conshdlrdata->nfailprovedinfeaslp++;
            }
            mpq_clear(dualobjval);

            /* stop timing and update number of calls */
            SCIPstopClock(scip, conshdlrdata->provedinfeaslptime);
            conshdlrdata->nprovedinfeaslp++;

            break;
         }
      case 'i':
         {
            SCIP_Real dualobjval;
            /* start timing */
            SCIPstartClock(scip, conshdlrdata->provedinfeaslptime);

            SCIP_CALL( provedBoundInterval(scip, conshdlrdata, consdata, TRUE, &dualobjval) );

            /* check if node can be cut off because safe infeasibility proof was found */
            if( dualobjval > 0.0 )
            {
               SCIPdebugMessage("   --> node is cut off by bounding (ray cost =%g)\n", dualobjval);
               *result = SCIP_CUTOFF;
            }
            else
            {
               SCIPdebugMessage("   --> node not cut off by provedbound interval (ray cost =%g)\n",dualobjval);

               /* update number of fails */
               conshdlrdata->nfailprovedinfeaslp++;
            }
            /* stop timing and update number of calls */
            SCIPstopClock(scip, conshdlrdata->provedinfeaslptime);
            conshdlrdata->nprovedinfeaslp++;
            break;
         }
      case 'n':
      case 'v':
      case 'e':
         break;
      default:
         SCIPerrorMessage("invalid parameter setting <%c> for dual bounding method\n", dualboundmethod);
         return SCIP_PARAMETERWRONGVAL;
      }
   }

   /* solve LP exactly in order to show that node can be pruned or to find a solution to branch on */
   if( *result != SCIP_CUTOFF )
   {
      SCIPdebugMessage("   as a last resort: proving infeasibility by dbmethod <e>\n");

      /* start timing */
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         if( dualboundmethod == 'e' )
            SCIPstartClock(scip, conshdlrdata->provedinfeaslptime);
         else
            SCIPstartClock(scip, conshdlrdata->exactinfeaslptime);
      }
      else
         SCIPstartClock(scip, conshdlrdata->exactunsollptime);

      /* constructs exact LP of current node */
      constructCurrentLPEX(scip, conshdlrdata, consdata);

      /* updates lower bound of current node wrt the pseudo objective value */
      if( !SCIPuseFPRelaxation(scip) )
      {
         SCIP_CALL( SCIPupdateLocalLowerbound(scip,
               mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) );
      }

      /* if the inexact LP was solved at the node, load LP state from inexact LP into exact LP solver */
      if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_NOTSOLVED )
      {
         SCIP_CALL( loadLPState(scip, conshdlrdata) );
      }

      /** @todo exiptodo:
       *  - which simplex algorithm sould be used here? selection by warmstartinfo is not useful; maybe use the same
       *    algorithm as in the inexact LP solver.
       */
      /* solve exact LP */
      algo = 'd';
      switch( algo )
      {
      case 'd':
         SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, TRUE, &lperror) );
         break;

      case 'p':
         SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_PRIMALSIMPLEX, TRUE, &lperror) );
         break;

      default:
         SCIPerrorMessage("invalid parameter setting <%c> for exact LP algorithm\n", algo);
         return SCIP_PARAMETERWRONGVAL;
      }

      /* stop timing and update number of calls */
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         if( dualboundmethod == 'e' )
         {
            SCIPstopClock(scip, conshdlrdata->provedinfeaslptime);
            conshdlrdata->nprovedinfeaslp++;
            if( lperror )
               conshdlrdata->nfailprovedinfeaslp++;
         }
         else
         {
            SCIPstopClock(scip, conshdlrdata->exactinfeaslptime);
            conshdlrdata->nexactinfeaslp++;
         }
      }
      else
      {
         SCIPstopClock(scip, conshdlrdata->exactunsollptime);
         conshdlrdata->nexactunsollp++;
         if( !lperror )
         {
            mpq_t dualobjval;

            mpq_init(dualobjval);
            SCIP_CALL( SCIPlpiexGetObjval(conshdlrdata->lpiex, &dualobjval) );
            SCIP_CALL( SCIPcomputeDualboundQuality(scip, conss[0], mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) );
            mpq_clear(dualobjval);
         }
      }

      if( lperror )
      {
         SCIPerrorMessage("exact LP solver returns error: case not handled yet\n");
         return SCIP_ERROR;
      }

      /* evaluate result of exact LP solver */
      SCIP_CALL( evaluateLPEX(scip, conss[0], conshdlrdata, consdata, result) );

      /* update number of wrong infeasible LP claims */
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE && ( *result == SCIP_BRANCHED || *result ==  SCIP_SOLVELP) )
         conshdlrdata->nwrongexactinfeaslp++;
   }

   SCIPdebugMessage("   --> enforcing pseudo solution returned result <%d>\n", *result);

   /* remember that current node is the one at which a pseudo solution was enforced last */
   conshdlrdata->lastenfopsnode = SCIPgetCurrentNode(scip);

   /** @todo exiptodo: extension
    *  - is it possible and useful to work on the pseudosolution in additon (on inexact pseudosolution from scip
    *    or exact pseudosolution via exactlp conshdlr)?
    */

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExactlp)
{  /*lint --e{715}*/
   SCIP_Bool feasible;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);
   assert(sol != NULL);
   assert(checklprows); /* in exact mode, constraint handler cannot draw conclusions from the LP's (inexact) feasibility */

   /*debugMessage("Check method of exactlp constraint\n");*/

   SCIP_CALL( checkCons(scip, SCIPconshdlrGetData(conshdlr), conss[0], sol, NULL, printreason, &feasible) );

   if( feasible )
      *result = SCIP_FEASIBLE;
   else
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropExactlp NULL
#endif


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolExactlp)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool cutoff;
   SCIP_Bool delay;
   mpq_t minactivity;
   mpq_t maxactivity;
   int oldnchgbds;
   int oldndelconss;
   int c;
   int i;
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
#endif

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Presol method of exactlp constraints\n");*/

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);

   /* remember old preprocessing counters */
   cutoff = FALSE;
   delay = FALSE;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;

   /* initialize gmp data */
   mpq_init(minactivity);
   mpq_init(maxactivity);

   /* process single exactlp constraints */
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(SCIPconsIsActive(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(SCIPgetNVars(scip) == consdata->nvars);

      /* we can only presolve exactlp constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;

      SCIPdebugMessage("presolving exactlp constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

      /* process single linear constraints of exactlp constraint */
      for( i = 0; i < consdata->nconss && !SCIPisStopped(scip); ++i )
      {
         SCIPdebugMessage("presolving linear constraint <%d> of exactlp constraint\n", i);

         /* check bounds */
         if( mpq_cmp(consdata->lhs[i], consdata->rhs[i]) > 0 )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d> of exactlp constraint is infeasible: sides=[%Qd,%Qd]\n",
               consdata->lhs[i], consdata->rhs[i]);
            SCIPdebugMessage(s);
#endif
            cutoff = TRUE;
            break;
         }

         /* tighten variable's bounds */
         SCIP_CALL( tightenBounds(scip, conshdlrdata, consdata, vars, i, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* check constraint for infeasibility and redundancy */
         singleConsGetActivityBounds(conshdlrdata, consdata, i, minactivity, maxactivity);

         if( mpq_cmp(minactivity, consdata->rhs[i]) > 0 || mpq_cmp(maxactivity, consdata->lhs[i]) < 0 )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "lincons<%d> of exactlp cons is infeas: activity=[%Qd,%Qd], sides=[%Qd,%Qd]\n",
               i, minactivity, maxactivity, consdata->lhs[i], consdata->rhs[i]);
            SCIPdebugMessage(s);
#endif
            cutoff = TRUE;
            break;
         }
         else if( mpq_cmp(minactivity, consdata->lhs[i]) >= 0 && mpq_cmp(maxactivity, consdata->rhs[i]) <= 0 )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "lincons<%d> of exactlp cons is redundant: activity=[%Qd,%Qd], sides=[%Qd,%Qd]\n",
               i, minactivity, maxactivity, consdata->lhs[i], consdata->rhs[i]);
            SCIPdebugMessage(s);
#endif
            delSingleCons(scip, cons, conshdlrdata, consdata, i);
            i--;
            (*ndelconss)++;
         }
      }

      if( cutoff )
         break;
   }

   /* free gmp data */
   mpq_clear(maxactivity);
   mpq_clear(minactivity);

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nchgbds > oldnchgbds || *ndelconss > oldndelconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExactlp NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExactlp)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int c;
   int j;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(consdata->lockup != NULL);
   assert(consdata->lockdown != NULL);

   /* get problem variables; in transforming stage we can only access original variables,
    * however, SCIPaddVarLocks() will use the transformed variable if it exists
    */
   vars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);

   /* for every variable, check whether rounding up or down could lead to the violation of some constraint */
   for( c = 0; c < consdata->nconss; ++c )
   {
      haslhs = !isNegInfinity(conshdlrdata, consdata->lhs[c]);
      hasrhs = !isPosInfinity(conshdlrdata, consdata->rhs[c]);

      /* check all variables of the current constraint */
      for( j = consdata->beg[c]; j < consdata->beg[c] + consdata->len[c]; ++j )
      {
         assert(consdata->ind[j] >= 0 && consdata->ind[j] < consdata->nvars);
         assert(consdata->lockup[consdata->ind[j]] >= 0 && consdata->lockdown[consdata->ind[j]] >= 0);

         /* val > 0 */
         if( mpq_sgn(consdata->val[j]) > 0 )
         {
            if( haslhs )
               consdata->lockdown[consdata->ind[j]]++;

            if( hasrhs )
               consdata->lockup[consdata->ind[j]]++;
         }

         /* val < 0 */
         if( mpq_sgn(consdata->val[j]) < 0 )
         {
            if( haslhs )
               consdata->lockup[consdata->ind[j]]++;

            if( hasrhs )
               consdata->lockdown[consdata->ind[j]]++;
         }
      }
   }

   /* set rounding locks for all variables */
   for( j = 0; j < consdata->nvars; ++j )
   {
      int probindex;

      probindex = SCIPvarGetProbindex(vars[j]);

      if( consdata->lockdown[probindex] > 0 && consdata->lockup[probindex] > 0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      else if( consdata->lockdown[probindex] > 0 )
      {
         assert(consdata->lockup[probindex] == 0);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos, nlocksneg) );
      }
      else if( consdata->lockup[probindex] > 0 )
      {
         assert(consdata->lockdown[probindex] == 0);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveExactlp NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveExactlp NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableExactlp NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableExactlp NULL
#endif


/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsExactlp NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExactlp)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdataPrint(scip, SCIPconshdlrGetData(conshdlr), SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
#define consCopyExactlp NULL

/** constraint parsing method of constraint handler */
#define consParseExactlp NULL

/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsExactlp NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsExactlp NULL
#endif

/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a exactlp constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdExactlp)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to exactlp constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to exactlp constraint\n", SCIPconsGetName(cons));

      /* create the exactlp constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsExactlp(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif




/*
 * Event handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real newbound;
   int varind;

   assert(eventdata != NULL);
   assert(eventdata->consdata != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);
   assert(0 <= consdata->nbndchglb && consdata->nbndchglb <= consdata->nvars);
   assert(0 <= consdata->nbndchgub && consdata->nbndchgub <= consdata->nvars);

   varind = eventdata->varind;
   assert(0 <= varind && varind < consdata->nvars);
   assert(SCIPvarGetProbindex(SCIPeventGetVar(event)) == varind);

   newbound = SCIPeventGetNewbound(event);

   conshdlrdata = SCIPconshdlrGetData(SCIPfindConshdlr(scip, CONSHDLR_NAME));

   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:

      /** @todo exiptodo: presolving extension
       *  - check whether the following argumentation is correct (currently we only have branching decisions anyway)
       */
      /* catch only bound changes that are not originated by the exactlp constraint handler (e.g., exact presolving),
       * i.e., those that come from branching decisions. the bound changes found by the exactlp constraint handler are
       * already stored here and the given bound is only a relaxation and probably weaker.
       * bound changes on integral variables can be catched in any case as they are identical to the ones already stored.
       */
      if( SCIPvarGetType(SCIPeventGetVar(event)) == SCIP_VARTYPE_CONTINUOUS )
         break;

#ifdef BOUNDCHG_OUT
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "event<%x>: store unprocessed bndchg: <%s>[%Qd,%Qd] >= %g (old sum: nchglb=%d,nchgub=%d)",
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lbloc[varind],
            consdata->ubloc[varind], newbound, consdata->nbndchglb, consdata->nbndchgub);
         SCIPdebugMessage(s);
      }
#endif

      /* check whether variable already has unprocessed lower bound change events */
      if( consdata->bndchglbpos[varind] == -1 )
      {
#ifdef BOUNDCHG_OUT
         printf("-->     first lbchg of var");
#endif
         /* store position of variable in unprocessed lower bound change arrays */
         consdata->bndchglbpos[varind] = consdata->nbndchglb;

         /* store new lower bound and index of variable */
         if( SCIPisInfinity(scip, -newbound) )
            mpq_set(consdata->bndchglb[consdata->nbndchglb], *negInfinity(conshdlrdata));
         else if( SCIPisInfinity(scip, newbound) )
            mpq_set(consdata->bndchglb[consdata->nbndchglb], *posInfinity(conshdlrdata));
         else
            mpq_set_d(consdata->bndchglb[consdata->nbndchglb], newbound);
         consdata->bndchglbind[consdata->nbndchglb] = varind;

         /* update number of unprocessed lower bound change events */
         consdata->nbndchglb++;
      }
      else
      {
#ifdef BOUNDCHG_OUT
         printf("--> NOT first lbchg of var");
#endif
         assert(0 <= consdata->bndchglbpos[varind] && consdata->bndchglbpos[varind] < consdata->nbndchglb);
         if( SCIPisInfinity(scip, -newbound) )
            mpq_set(consdata->bndchglb[consdata->bndchglbpos[varind]], *negInfinity(conshdlrdata));
         else if( SCIPisInfinity(scip, newbound) )
            mpq_set(consdata->bndchglb[consdata->bndchglbpos[varind]], *posInfinity(conshdlrdata));
         else
            mpq_set_d(consdata->bndchglb[consdata->bndchglbpos[varind]], newbound);
         assert(consdata->bndchglbind[consdata->bndchglbpos[varind]] == varind);
      }
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:

      /** @todo exiptodo: presolving extension
       *  - check whether the following argumentation is correct (currently we only have branching decisions anyway)
       */
      /* catch only bound changes that are not originated by the exactlp constraint handler (e.g., exact presolving),
       * i.e., those that come from branching decisions. the bound changes found by the exactlp constraint handler are
       * already stored here and the given bound is only a relaxation and probaly weaker.
       * bound changes on integral variables can be catched in any case as they are identical to the ones already stored.
       */
      if( SCIPvarGetType(SCIPeventGetVar(event)) == SCIP_VARTYPE_CONTINUOUS )
         break;

#ifdef BOUNDCHG_OUT
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "event<%x>: store unprocessed bndchg: <%s>[%Qd,%Qd] <= %g (old sum: nchglb=%d,nchgub=%d)",
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lbloc[varind],
            consdata->ubloc[varind], newbound, consdata->nbndchglb, consdata->nbndchgub);
         SCIPdebugMessage(s);
      }
#endif

      /* check whether variable already has unprocessed upper bound change events */
      if( consdata->bndchgubpos[varind] == -1 )
      {
#ifdef BOUNDCHG_OUT
         printf("-->     first ubchg of var");
#endif
         /* store position of variable in unprocessed upper bound change arrays */
         consdata->bndchgubpos[varind] = consdata->nbndchgub;

         /* store new upper bound and index of variable */
         if( SCIPisInfinity(scip, -newbound) )
            mpq_set(consdata->bndchgub[consdata->nbndchgub], *negInfinity(conshdlrdata));
         else if( SCIPisInfinity(scip, newbound) )
            mpq_set(consdata->bndchgub[consdata->nbndchgub], *posInfinity(conshdlrdata));
         else
            mpq_set_d(consdata->bndchgub[consdata->nbndchgub], newbound);
         consdata->bndchgubind[consdata->nbndchgub] = varind;

         /* update number of unprocessed upper bound change events */
         consdata->nbndchgub++;
      }
      else
      {
#ifdef BOUNDCHG_OUT
         printf("--> NOT first ubchg of var");
#endif
         assert(0 <= consdata->bndchgubpos[varind] && consdata->bndchgubpos[varind] < consdata->nbndchgub);

         if( SCIPisInfinity(scip, -newbound) )
            mpq_set(consdata->bndchgub[consdata->bndchgubpos[varind]], *negInfinity(conshdlrdata));
         else if( SCIPisInfinity(scip, newbound) )
            mpq_set(consdata->bndchgub[consdata->bndchgubpos[varind]], *posInfinity(conshdlrdata));
         else
            mpq_set_d(consdata->bndchgub[consdata->bndchgubpos[varind]], newbound);
         assert(consdata->bndchgubind[consdata->bndchgubpos[varind]] == varind);
      }
      break;
   default:
      SCIPerrorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= consdata->nbndchglb && consdata->nbndchglb <= consdata->nvars);
   assert(0 <= consdata->nbndchgub && consdata->nbndchgub <= consdata->nvars);

#ifdef BOUNDCHG_OUT
   printf(" (new sum: nchglb=%d, nchgub=%d)\n", consdata->nbndchglb, consdata->nbndchgub);
#endif

#ifndef NDEBUG
   {
      SCIP_VAR** vars;
      SCIP_VAR* var;


      /* get problem variables */
      var = SCIPeventGetVar(event);
      vars = SCIPgetVars(scip);
      assert(var != NULL);
      assert(vars != NULL);
      assert(SCIPgetNVars(scip) == consdata->nvars);
      assert(vars[varind] == var);
      assert(SCIPvarGetProbindex(var) == varind);
   }
#endif

   return SCIP_OKAY;
}



/*
 * constraint specific interface methods
 */

/** creates the handler for exactlp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExactlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecExactlp,
         eventhdlrdata) );

   /* create exactlp constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS, CONSHDLR_PROP_TIMING,
         conshdlrCopyExactlp, consFreeExactlp, consInitExactlp, consExitExactlp,
         consInitpreExactlp, consExitpreExactlp, consInitsolExactlp, consExitsolExactlp,
         consDeleteExactlp, consTransExactlp, consInitlpExactlp,
         consSepalpExactlp, consSepasolExactlp, consEnfolpExactlp, consEnfopsExactlp, consCheckExactlp,
         consPropExactlp, consPresolExactlp, consRespropExactlp, consLockExactlp,
         consActiveExactlp, consDeactiveExactlp,
         consEnableExactlp, consDisableExactlp, consDelvarsExactlp,
         consPrintExactlp, consCopyExactlp, consParseExactlp,
         consGetVarsExactlp, consGetNVarsExactlp, conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdExactlp, LINCONSUPGD_PRIORITY) );
#endif

   /* add exactlp constraint handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/exactlp/psobjweight",
         "weight of the original objective function in lp to compute interior point",
         &conshdlrdata->psobjweight, TRUE, DEFAULT_PSOBJWEIGHT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/exactlp/psreduceauxlp",
         "should the number of constraints in lp to compute interior point be reduced?",
         &conshdlrdata->psreduceauxlp, TRUE, DEFAULT_PSREDUCEAUXLP, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "constraints/exactlp/psdualcolselection",
         "strategy to select which dual columns to use for lp to compute interior point ('n'o selection (all cols), 'a'ctive rows from exact primal LP, 'A'ctive rows from inexact primal LP)",
         &conshdlrdata->psdualcolselection, TRUE, DEFAULT_PSDUALCOLSELECTION, "naAbB", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "constraints/exactlp/psintpointselection",
         "method to select interior point ('a'rbitrary interior point, 'o'ptimized interior point, 'A'rbitrary interior point solved in dual form, 't'wo stage optimized interior point)",
         &conshdlrdata->psintpointselection, TRUE, DEFAULT_PSINTPOINTSELECTION, "aoAt", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/exactlp/psuseintpoint",
         "should correction shift use an interior pt? (otherwise use interior ray of recession cone)",
         &conshdlrdata->psuseintpoint, TRUE, DEFAULT_PSUSEINTPOINT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/exactlp/psinfeasversion",
         "should infeasibility version of project-and-shift be used?",
         &conshdlrdata->psinfeasversion, TRUE, DEFAULT_PSINFEASVERSION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/exactlp/interleavedbfreq",
         "frequency at which dual bounding strategy is interleaved (-1: never, 0: if prommising, x: xth node)",
         &conshdlrdata->interleavedbfreq, TRUE, DEFAULT_INTERLEAVEDBFREQ, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/exactlp/computedbdiff",
         "should the quality of the safe dual bounds be computed (involves solving LPs exactly)?",
         &conshdlrdata->computedbdiff, TRUE, DEFAULT_COMPUTEDBDIFF, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/exactlp/smalldbdiff",
         "maximum value of small relativ difference of safe dual bound",
         &conshdlrdata->smalldbdiff, TRUE, DEFAULT_SMALLDBDIFF, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/exactlp/mediumdbdiff",
         "maximum value of medium relativ difference of safe dual bound",
         &conshdlrdata->mediumdbdiff, TRUE, DEFAULT_MEDIUMDBDIFF, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a exactlp constraint */
SCIP_RETCODE SCIPcreateConsExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   int                   ninfbounds,         /**< number of variables with infinite bound in safe dual bounding method */
   int                   ninfintbounds,      /**< number of integer variables with infinite bound in safe db method */
   int                   nlargebounds,       /**< number of variables with large bound in safe dual bounding method */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,             /**< number of constraints */
   int                   nsplitconss,        /**< number of constraints we would have to be split for a FP-relaxation */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int                   nintegral,          /**< number of integral nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each constraint in ind and val array */
   int*                  len,                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries (and some zeros) */
   mpq_t                 minabsval,          /**< minimum absolute nonzero constraint matrix, lhs, or rhs entry */
   mpq_t                 maxabsval,          /**< maximum absolute nonzero constraint matrix, lhs, or rhs entry */
   SCIP_Bool             objneedscaling,     /**< do objective function values need to be scaled because some are not FP representable? */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(nnonz > 0 || mpq_sgn(minabsval) > 0);
   assert(nnonz > 0 || mpq_sgn(maxabsval) > 0);

   SCIPdebugMessage("create exactlp constraint\n");

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, conshdlrdata, &consdata, conshdlrdata->eventhdlr, objsense, nvars, ninfbounds,
         ninfintbounds, nlargebounds, obj, lb, ub, nconss, nsplitconss, nconss, lhs, rhs, nnonz, nintegral, beg, len, ind,
         val, minabsval, maxabsval, objneedscaling) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   SCIPdebugMessage("create exactlp constraint (done)\n");

   return SCIP_OKAY;
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
SCIP_Bool SCIPisPosInfinityExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   const mpq_t           val                 /**< value to be compared against infinity */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return isPosInfinity(conshdlrdata, val);
}

/** checks if value is treated as negative infinite in exactlp constraint handler */
SCIP_Bool SCIPisNegInfinityExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   const mpq_t           val                 /**< value to be compared against infinity */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return isNegInfinity(conshdlrdata, val);
}

/** returns a safe external value for the given exact internal objective value, i.e., a lower or upper approximation
 *  of the exact external value depending on the objective sense of the original problem and whether the given internal
 *  value is a lower or upper bound in the transformed problem.
 */
SCIP_Real SCIPgetExternSafeObjval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   SCIP_Real             objval,             /**< safe internal objective value */
   SCIP_Bool             lowerbound          /**< TRUE if objval is lower bound in transformed prob; FALSE if upper bound */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real objextern;
   SCIP_Bool down;
   mpq_t mpqobjintern;
   mpq_t mpqobjextern;

   assert(cons != NULL);

   /* get exactlp constraint handler */
   conshdlr = SCIPconsGetHdlr(cons);

   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   mpq_init(mpqobjintern);
   mpq_init(mpqobjextern);

   if( SCIPisInfinity(scip, objval) )
      mpq_set(mpqobjintern, *posInfinity(conshdlrdata));
   else if( SCIPisInfinity(scip, -objval) )
      mpq_set(mpqobjintern, *negInfinity(conshdlrdata));
   else
      mpq_set_d(mpqobjintern, objval);

   getExternObjvalex(cons, mpqobjintern, mpqobjextern);

   down = (((SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE) && lowerbound)
      || ((SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE) && !lowerbound));

   if( down )
      objextern = mpqGetRealRelax(scip, mpqobjextern, GMP_RNDD);
   else
      objextern = mpqGetRealRelax(scip, mpqobjextern, GMP_RNDU);

   mpq_clear(mpqobjintern);
   mpq_clear(mpqobjextern);

   return objextern;
}

/** gets number of feasible exact primal solutions stored in the exact solution storage */
int SCIPgetNSolexs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->primal != NULL )
      return conshdlrdata->primal->nsols;
   else
      return 0;
}

/** gets best feasible exact primal solution found so far, or NULL if no solution has been found */
SCIP_SOLEX* SCIPgetBestSolex(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return NULL;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->primal != NULL && conshdlrdata->primal->nsols > 0 )
   {
      assert(conshdlrdata->primal->sols != NULL);
      assert(conshdlrdata->primal->sols[0] != NULL);
      return conshdlrdata->primal->sols[0];
   }

   return NULL;
}

/** returns objective value of exact primal CIP solution w.r.t. original problem */
void SCIPgetSolexOrigObj(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_CONS*           cons,               /**< exactlp constraint data */
   SCIP_SOLEX*          sol,                /**< exact primal solution */
   mpq_t                obj                 /**< pointer to store objective value */
   )
{
   mpq_t objintern;

   assert(sol != NULL);

   mpq_init(objintern);
   SCIPsolexGetObj(sol, objintern);

   getExternObjvalex(cons, objintern, obj);

   mpq_clear(objintern);
}

/** returns transformed objective value of exact primal CIP solution */
void SCIPgetSolexTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol,                /**< exact primal solution */
   mpq_t                 obj                 /**< pointer to store objective value */
   )
{
   assert(sol != NULL);

   SCIPsolexGetObj(sol, obj);
}

/** returns objective value of best exact primal CIP solution found so far w.r.t. original problem */
void SCIPgetBestSolexObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   mpq_t                 obj                 /**< pointer to store objective value */
   )
{
   SCIP_SOLEX* sol;

   /* get best exact solution found so far */
   sol = SCIPgetBestSolex(scip);
   if( sol != NULL )
   {
      SCIPgetSolexOrigObj(scip, cons, sol, obj);
   }
   else
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSHDLR* conshdlr;

      /* find the exactlp constraint handler */
      conshdlr = SCIPconsGetHdlr(cons);

      if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
      {
         SCIPerrorMessage("constraint is not of type exactlp\n");
         SCIPABORT();
      }

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);

      assert(conshdlrdata != NULL);
      getExternObjvalex(cons, *posInfinity(conshdlrdata), obj);
   }
}

/** returns transformed objective value of best exact primal CIP solution found so far */
void SCIPgetBestSolexTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   mpq_t                 obj                 /**< pointer to store objective value */
   )
{
   SCIP_SOLEX* sol;

   /* get best exact solution found so far */
   sol = SCIPgetBestSolex(scip);
   if( sol != NULL )
   {
      SCIPgetSolexTransObj(scip, sol, obj);
   }
   else
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSHDLR* conshdlr;

      /* find the exactlp constraint handler */
      conshdlr = SCIPconsGetHdlr(cons);

      if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
      {
         SCIPerrorMessage("constraint is not of type exactlp\n");
         SCIPABORT();
      }

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);

      assert(conshdlrdata != NULL);
      mpq_set(obj, *posInfinity(conshdlrdata));
   }
}

/** outputs non-zero variables of exact solution in original problem space to file stream */
SCIP_RETCODE SCIPprintSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   SCIP_SOLEX*           sol,                /**< exact primal solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   char s[SCIP_MAXSTRLEN];
   mpq_t obj;
   int n;

   mpq_init(obj);

   SCIPgetSolexOrigObj(scip, cons, sol, obj);
   n = gmp_snprintf(s, SCIP_MAXSTRLEN, "objective value:                 %20Qd\n", obj);
   if( n >= SCIP_MAXSTRLEN )
   {
      char* bigs;

      SCIP_CALL( SCIPallocMemorySize(scip, &bigs, n+1) );
      gmp_snprintf(bigs, n+1, "objective value:                 %20Qd\n", obj);
      SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), bigs);
      SCIPfreeMemory(scip, &bigs);
   }
   else
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, s);

   /** @todo exiptodo: using scip->origprob and scip->transprob is not SCIP conform
    *  - think about a more general way to get information from cons_exactlp.c.
    *  - maybe, implement all methods that operate/use exact problem data (like solex methods) in cons_exactlp.c
    */
   SCIP_CALL( SCIPsolexPrint(scip, sol, scip->origprob, scip->transprob, file, printzeros) );

   mpq_clear(obj);

   return SCIP_OKAY;
}

/** outputs non-zero variables of exact solution in original problem space in transformed problem space to file stream */
SCIP_RETCODE SCIPprintTransSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol,                /**< exact primal solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   char s[SCIP_MAXSTRLEN];
   mpq_t obj;

   mpq_init(obj);

   if( SCIPsolexGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIPerrorMessage("cannot print original space solution as transformed solution\n");
      return SCIP_INVALIDCALL;
   }

   SCIPgetSolexTransObj(scip, sol, obj);
   gmp_snprintf(s, SCIP_MAXSTRLEN, "objective value:                 %20Qd\n", obj);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, s);

   /** @todo exiptodo: using scip->origprob and scip->transprob is not SCIP conform
    *  - think about a more general way to get information from cons_exactlp.c.
    *  - maybe, implement all methods that operate/use exact problem data (like solex methods) in cons_exactlp.c
    */
   SCIP_CALL( SCIPsolexPrint(scip, sol, scip->transprob, scip->transprob, file, printzeros) );

   mpq_clear(obj);

   return SCIP_OKAY;
}

/** outputs best feasible exact primal solution found so far to file stream */
SCIP_RETCODE SCIPprintBestSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< exactlp constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOLEX* sol;

   sol = SCIPgetBestSolex(scip);

   if( sol == NULL )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "no solution available\n");
   else
   {
      SCIP_CALL( SCIPprintSolex(scip, cons, sol, file, printzeros) );
   }
   return SCIP_OKAY;
}

/** outputs best feasible exact primal solution found so far in transformed problem space to file stream */
SCIP_RETCODE SCIPprintBestTransSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOLEX* sol;

   sol = SCIPgetBestSolex(scip);

   if( sol == NULL )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "no exact solution available\n");
   else
   {
      if( SCIPsolexGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
      {
         SCIPdialogMessage(scip, NULL, "best exact solution exists only in original problem space\n");
      }
      else
      {
         SCIP_CALL( SCIPprintTransSolex(scip, sol, file, printzeros) );
      }
   }

   return SCIP_OKAY;
}

/** outputs value of variable in best feasible exact primal solution found so far to file stream */
SCIP_RETCODE SCIPprintBestSolexVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable for which solution value should be printed */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_SOLEX* sol;
   char s[SCIP_MAXSTRLEN];
   mpq_t solval;

   sol = SCIPgetBestSolex(scip);

   if( sol == NULL )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "no exact solution available\n");
   else
   {
      mpq_init(solval);

      SCIPsolexGetVal(sol, var, solval);

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%-32s", SCIPvarGetName(var));
      gmp_snprintf(s, SCIP_MAXSTRLEN, " %20Qd\n", solval);
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, s);

      mpq_clear(solval);
   }

   return SCIP_OKAY;
}


/** checks best exact primal solution for feasibility without adding it to the solution store;
 *  called for original exactlp constraints the method is used to double check the best exact solution in order to
 *  validate the presolving process
 */
SCIP_RETCODE SCIPcheckBestSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOLEX* sol;

   assert(scip != NULL);
   assert(feasible != NULL);

   *feasible = TRUE;

   /* find the exactlp constraint handler */
   conshdlr = SCIPconsGetHdlr(cons);

   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get best exact primal solution */
   sol = SCIPgetBestSolex(scip);
   if( sol == NULL )
   {
      SCIPdialogMessage(scip, NULL, "no feasible exact solution available\n");
   }
   else
   {
      SCIP_CALL( checkCons(scip, conshdlrdata, cons, NULL, sol, printreason, feasible) );
   }

   return SCIP_OKAY;
}

/** gets exact objective function value of variable */
void SCIPvarGetObjExactlp(
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   mpq_t                 obj                 /**< pointer to store objective value */
   )
{
   SCIP_CONSDATA* consdata;

   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_set(obj, consdata->obj[SCIPvarGetProbindex(var)]);
}

/** gets exact global lower bound of variable */
void SCIPvarGetLbGlobalExactlp(
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   mpq_t                 lb                  /**< pointer to store global lower bound */
   )
{
   SCIP_CONSDATA* consdata;

   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_set(lb, consdata->lb[SCIPvarGetProbindex(var)]);
}

/** gets exact global upper bound of variable */
void SCIPvarGetUbGlobalExactlp(
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   mpq_t                 ub                  /**< pointer to store global upper bound */
   )
{
   SCIP_CONSDATA* consdata;

   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_set(ub, consdata->ub[SCIPvarGetProbindex(var)]);
}

/** gets worst global bound of variable with respect to the objective function */
void SCIPvarGetWorstBoundGlobalExactlp(
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   mpq_t                 bound               /**< pointer to store worst bound */
   )
{
   SCIP_CONSDATA* consdata;

   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( mpq_sgn(consdata->obj[SCIPvarGetProbindex(var)]) >= 0.0 )
      mpq_set(bound, consdata->ub[SCIPvarGetProbindex(var)]);
   else
      mpq_set(bound, consdata->lb[SCIPvarGetProbindex(var)]);
}

/** returns safe dual bounding method to be applied; if user want's the solver to decide the most promising one is
 *  selected, otherwise the one the user wanted is returned
 *  note, for the automatic mode, this method asumes that Neumaier and Shcherbina was already tested,
 *  i.e., that we know whether it suceeded.
 */
char SCIPselectDualBoundMethod(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             infeaslp            /**< will dual bound method be applied to safely verify infeasible LP? */
   )
{
   SCIP_CONS** conss;
   char dualboundmethod;
   SCIP_Bool skip;
   SCIP_Bool projectshift;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get exactlp constraints */
   conss = SCIPgetConss(scip);
   assert(conss != NULL);
   assert(SCIPgetNConss(scip) == 1);

   /* get exactlp constraint data */
   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   skip = FALSE;
   projectshift = TRUE;

   if( infeaslp )
   {
      /* user did choose a dual bounding method */
      if( SCIPdualBoundMethod(scip) != 'a' )
         dualboundmethod = SCIPdualBoundMethod(scip);
      else
      {
         /** @todo exiptodo: using scip->lp->hasprovedbound is not SCIP conform.
          *  - implement method in lp.c or scip.c which returns lp->hasprovedbound
          *  - there is at least one other position where hasprovedbound accessed directly
          */
         /* enfops which calls this methods for infeaslp is never called if Neumair Shcherbina proved infeasibility */
         assert(!scip->lp->hasprovedbound);

         /* decide whether to apply project-and-shift or exact LP */

         /* we already know whether project-and-shift works */
         if( conshdlrdata->psdatacon )
         {
            if( conshdlrdata->psdatafail )
               projectshift = FALSE;
            else
               projectshift = TRUE;
         }
         /* we don't know whether it will work and thus predict it
          * criterium: small number of nonzeros in the constraint matrix
          */
         else
         {
            assert(SCIPgetNNonzExactlp(conss[0]) >= 0);
            if( SCIPgetNNonzExactlp(conss[0]) <= SPARSEMAXNNONZ )
               projectshift = TRUE;
            else
               projectshift = FALSE;
         }

         if( projectshift )
            dualboundmethod = 'p';
         else
            dualboundmethod = 'e';
      }
   }
   else
   {
      /* decide whether we want to interleave with exact LP call/basis verification
       * - given freq
       * or
       * - Neumair Shcherbina bound only nearly able to cutoff node
       */
      if( (conshdlrdata->interleavedbfreq > 0 && !SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) && SCIPgetDepth(scip) > 0
            && SCIPgetDepth(scip) % (conshdlrdata->interleavedbfreq) == 0)
         || (conshdlrdata->interleavedbfreq == 0 && SCIPisGE(scip, SCIPgetLocalLowerbound(scip), SCIPgetCutoffbound(scip))
            && SCIPgetLocalLowerbound(scip) < SCIPgetCutoffbound(scip)) )
      {
         dualboundmethod = 'e';
      }
      else
      {
         /* user did choose a dual bounding method */
         if( SCIPdualBoundMethod(scip) != 'a' )
            dualboundmethod = SCIPdualBoundMethod(scip);
         /* do not recompute the dual bound if Neumair Shcherbina did succeed */
         else if( scip->lp->hasprovedbound )
         {
            dualboundmethod = 'n';
         }
         else
         {
            /* decide whether it is worse to compute a dual bound:
             * - never skip
             */
            skip = FALSE;
            if( skip )
               dualboundmethod = 'n';
            else
            {
               /* decide whether to apply project-and-shift or basis verification:
                * - we already know whether project-and-shift works
                * or
                * - small number of nonzeros in the constraint matrix
                */
               if( conshdlrdata->psdatacon )
               {
                  if( conshdlrdata->psdatafail )
                     projectshift = FALSE;
                  else
                     projectshift = TRUE;
               }
               else
               {
                  assert(SCIPgetNNonzExactlp(conss[0]) >= 0);
                  if( SCIPgetNNonzExactlp(conss[0]) <= SPARSEMAXNNONZ )
                     projectshift = TRUE;
                  else
                     projectshift = FALSE;
               }

               if( projectshift )
                  dualboundmethod = 'p';
               else
                  dualboundmethod = 'v';
            }
         }
      }
   }
   SCIPdebugMessage("   selected dual bounding method for <%s> LP: %c\n", infeaslp ? "infeasible" : "feasible",
      dualboundmethod);
#ifdef DBAUTO_OUT
   printf("selected dual bounding method for <%s> LP: %c (skip safedb:%d, many largebounds:%d [<%d>/<%d>=%.2f], many nonzeros:%d [<%d>])\n",
      infeaslp ? "infeasible" : "feasible", dualboundmethod, skip,
      (SCIP_Real)(SCIPgetNInfiniteBounds(conss[0]) + SCIPgetNLargeBounds(conss[0]))/(SCIP_Real)SCIPgetNVars(scip) > FEWLBOUNDSRATIO,
      SCIPgetNInfiniteBounds(conss[0]) + SCIPgetNLargeBounds(conss[0]), SCIPgetNVars(scip),
      (SCIP_Real)(SCIPgetNInfiniteBounds(conss[0]) + SCIPgetNLargeBounds(conss[0]))/(SCIP_Real)SCIPgetNVars(scip),
      SCIPgetNNonzExactlp(conss[0]) > SPARSEMAXNNONZ, SCIPgetNNonzExactlp(conss[0]));
#endif

   return dualboundmethod;
}

/*
 * problem and solving process statistics methods
 */

/** gets total number of linear constraints handled by exactlp constraint */
int SCIPgetNConssExactlp(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nconss;
}

/** gets total number of constraints handled by exactlp constraint that were split */
int SCIPgetNSplitconssExactlp(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nsplitconss;
}

/** gets number of nonzero entries in linear constraints handled by exactlp constraint */
int SCIPgetNNonzExactlp(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nnonz;
}

/** gets number of integral nonzero entries in linear constraints handled by exactlp constraint */
int SCIPgetNIntegralExactlp(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nintegral;
}

/** gets ration of maximum and minimum absolute nonzero constraint matrix, lhs, or rhs entry */
SCIP_Real SCIPgetCoefRatioExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;
   mpq_t ratio;
   mpq_t intratio;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_init(ratio);
   mpq_init(intratio);

   assert(mpq_sgn(consdata->minabsval) > 0);
   assert(mpq_sgn(consdata->maxabsval) > 0);
   mpq_div(ratio, consdata->maxabsval, consdata->minabsval);

   mpqCeil(intratio, ratio);
   assert(mpqIsIntegral(intratio));

   return mpqGetRealApprox(scip, intratio);
}

/** gets total number variables with infinite bound needed in safe dual bounding method */
int SCIPgetNInfiniteBounds(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->ninfbounds;
}

/** gets total number integer variables with infinite bound needed in safe dual bounding method */
int SCIPgetNInfiniteIntegerBounds(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->ninfintbounds;
}

/** gets total number variables with infinite bound needed in safe dual bounding method */
int SCIPgetNLargeBounds(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
      return 0;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nlargebounds;
}


/** gets time needed for safe dual bound computation for feasible LPs */
SCIP_Real SCIPgetProvedfeaslpTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIPgetClockTime(scip, conshdlrdata->provedfeaslptime);
}

/** gets time needed for safe verification for infeasible LPs */
SCIP_Real SCIPgetProvedinfeaslpTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIPgetClockTime(scip, conshdlrdata->provedinfeaslptime);
}

/** gets number of times, safe dual bound computation for feasible LPs was called */
SCIP_Longint SCIPgetNProvedfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nprovedfeaslp;
}

/** number of times, safe verification for infeasible LPs was called */
SCIP_Longint SCIPgetNProvedinfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nprovedinfeaslp;
}

/** gets number of times, safe dual bound computation for feasible LPs failed */
SCIP_Longint SCIPgetNFailProvedfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nfailprovedfeaslp;
}

/** number of times, safe verification for infeasible LPs failed */
SCIP_Longint SCIPgetNFailProvedinfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nfailprovedinfeaslp;
}

/** number of times, safe verification for infeasible LPs was aborted because no approximate certificate was given */
SCIP_Longint SCIPgetNAbortProvedinfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nabortprovedinfeaslp;
}

/** gets time needed for exact integrality verification for feasible LPs */
SCIP_Real SCIPgetExactfeaslpTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIPgetClockTime(scip, conshdlrdata->exactfeaslptime);
}

/** gets time needed for exact verification for infeasible LPs */
SCIP_Real SCIPgetExactinfeaslpTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIPgetClockTime(scip, conshdlrdata->exactinfeaslptime);
}

/** gets time needed for exact LP solve for unsolved LP */
SCIP_Real SCIPgetExactunsollpTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0.0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIPgetClockTime(scip, conshdlrdata->exactunsollptime);
}

/** gets number of times, exact integrality verification for feasible LPs was called */
SCIP_Longint SCIPgetNExactfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nexactfeaslp;
}

/** gets number of times, exact verification for infeasible LPs was called */
SCIP_Longint SCIPgetNExactinfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nexactinfeaslp;
}

/** gets number of times, exact LP solver is called for unsolved LP */
SCIP_Longint SCIPgetNExactunsollp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nexactunsollp;
}

/** gets number of times, claimed integral LP solution was fractional */
SCIP_Longint SCIPgetNWrongExactfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nwrongexactfeaslp;
}

/** gets number of times, claimed integral LP solution was infeasible and even the LP was infeasible */
SCIP_Longint SCIPgetNWronginfeasExactfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nwronginfeasexactfeaslp;
}

/** gets number of times, claimed infeasible LP was feasible */
SCIP_Longint SCIPgetNWrongExactinfeaslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return 0;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nwrongexactinfeaslp;
}

/** computes relative difference of given safe dual bound and exact dual bound at current node and updates statistics */
SCIP_RETCODE SCIPcomputeDualboundQuality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             safedualbound       /**< safe dual bound to compute quality for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool lperror;

   assert(!SCIPisInfinity(scip, -safedualbound));

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   /* only compute relative dual bound difference if wanted */
   if( !conshdlrdata->computedbdiff )
      return SCIP_OKAY;

   SCIPdebugMessage("  computing dual bound quality for exactlp constraint <%s>:\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* constructs exact LP of current node */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   /* load LP state from inexact LP into exact LP solver */
   SCIP_CALL( loadLPState(scip, conshdlrdata) );

   /* solve exact LP */
   SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, FALSE, &lperror) );

   if( !lperror )
   {
      mpq_t dualobjval;
      mpq_t absdualobjval;
      mpq_t dbdiff;
      mpq_t dualbound;
      mpq_t smalldbdiff;
      mpq_t mediumdbdiff;

      mpq_init(dualobjval);
      mpq_init(absdualobjval);
      mpq_init(dbdiff);
      mpq_init(dualbound);
      mpq_init(smalldbdiff);
      mpq_init(mediumdbdiff);

      SCIP_CALL( SCIPlpiexGetObjval(conshdlrdata->lpiex, &dualobjval) );
      mpq_set_d(dualobjval, mpqGetRealRelax(scip, dualobjval, GMP_RNDD));

      mpq_abs(absdualobjval, dualobjval);
      mpq_set_d(dualbound, safedualbound);
      mpq_set_d(smalldbdiff, conshdlrdata->smalldbdiff);
      mpq_set_d(mediumdbdiff, conshdlrdata->mediumdbdiff);

      mpq_sub(dbdiff, dualobjval, dualbound);

#ifdef SCIP_DEBUG
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "     exactdb<%Qd> (%.14f) - safedb<%Qd> (%.14f) = <%Qd> (%.14f)\n",
            dualobjval, mpqGetRealApprox(scip, dualobjval), dualbound, mpqGetRealApprox(scip, dualbound),
            dbdiff, mpqGetRealApprox(scip, dbdiff));
         SCIPdebugMessage(s);
      }
#endif

      /* compute relative difference (exact db - safe db)/max(1.0, |exact db|) */
      if( mpq_cmp_ui(absdualobjval, 1, 1) < 0 )
         mpq_set_ui(absdualobjval, 1, 1);
      mpq_div(dbdiff, dbdiff, absdualobjval);

      assert(mpq_sgn(dbdiff) >= 0);


      if( mpq_sgn(dbdiff) == 0 )
      {
         assert(mpq_equal(dualobjval, dualbound) != 0);
         conshdlrdata->nzerodbdiff++;
      }
      else if( mpq_cmp(dbdiff, smalldbdiff) <= 0 )
         conshdlrdata->nsmalldbdiff++;
      else if( mpq_cmp(dbdiff, mediumdbdiff) <= 0 )
         conshdlrdata->nmediumdbdiff++;
      else
         conshdlrdata->nlargedbdiff++;

      mpq_clear(dualobjval);
      mpq_clear(absdualobjval);
      mpq_clear(dbdiff);
      mpq_clear(dualbound);
      mpq_clear(smalldbdiff);
      mpq_clear(mediumdbdiff);
   }

   return SCIP_OKAY;
}

/** gets numbers of safe dual bounds with zero, small, medium, and large relative difference */
void SCIPgetNDualboundDiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         nzerodbdiff,        /**< pointer to store number of dual bounds with zero relativ difference */
   SCIP_Longint*         nsmalldbdiff,       /**< pointer to store number of dual bounds with small relativ difference */
   SCIP_Longint*         nmediumdbdiff,      /**< pointer to store number of dual bounds with medium relativ difference */
   SCIP_Longint*         nlargedbdiff        /**< pointer to store number of dual bounds with large relativ difference */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert(nzerodbdiff != NULL);
   assert(nsmalldbdiff != NULL);
   assert(nmediumdbdiff != NULL);
   assert(nlargedbdiff != NULL);

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *nzerodbdiff = conshdlrdata->nzerodbdiff;
   *nsmalldbdiff = conshdlrdata->nsmalldbdiff;
   *nmediumdbdiff = conshdlrdata->nmediumdbdiff;
   *nlargedbdiff = conshdlrdata->nlargedbdiff;

   assert(*nzerodbdiff >= 0);
   assert(*nsmalldbdiff >= 0);
   assert(*nmediumdbdiff >= 0);
   assert(*nlargedbdiff >= 0);
}

/** returns whether the quality of the safe dual bounds was computed */
SCIP_Bool SCIPhasDualboundDiff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      SCIPABORT();
      return FALSE;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->computedbdiff;
}

#endif
