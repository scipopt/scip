/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_exactlp.c,v 1.1.2.2 2009/07/24 12:52:51 bzfwolte Exp $"
//#define SCIP_DEBUG /*??????????????*/
//#define LP_OUT /* only for debugging ???????????????? */
//#define BOUNDCHG_OUT /* only for debugging ?????????? */
//#define BASIS_OUT /* only for debugging ???????????????? */
//#define EXLPSOLVER_OUT /* only for debugging ???????????????? */
//#define DETAILED_DEBUG /* only for debugging ???????????????? */

/**@file   cons_exactlp.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for exactlp constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <gmp.h> 

#include "scip/cons_exactlp.h"
#include "scip/struct_lp.h" /* only for debugging ??????????*/
#include "scip/struct_scip.h" /* only for debugging ??????????*/
#include "EGlib.h" 
#include "QSopt_ex.h" 


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

#define EVENTHDLR_NAME        "exactlp"
#define EVENTHDLR_DESC        "bound change event handler for exactlp constraints"



/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   /* ???? todo: this should be moved to an appropriat place and has to be assigned in a consistent way (scipex/def.h oder defex.h) */
   mpq_t                 posinfinity;        /**< value considered to be infinity */
   mpq_t                 neginfinity;        /**< value considered to be infinity */
   SCIP_LPIEX*           lpiex;              /**< Exact LP solver interface */
   SCIP_Bool             lpexconstructed;    /**< was the exact LP of some prior node already constructed (constraints)? */
   SCIP_NODE*            lastenfopsnode;     /**< last node at which enfops was called */ 
};

/** constraint data for exactlp constraints */
struct SCIP_ConsData
{
   SCIP_OBJSEN           objsense;           /**< objective sense */
   int                   nvars;              /**< number of variables */
   mpq_t*                obj;                /**< objective function values of variables */
   mpq_t*                lb;                 /**< lower bounds of variables */
   mpq_t*                ub;                 /**< upper bounds of variables */
   int*                  bndchglbpos;        /**< positions of variables in unprocessed lower bound change arrays, or -1 */
   int*                  bndchgubpos;        /**< positions of variables in unprocessed upper bound change arrays, or -1 */
   int                   nconss;             /**< number of constraints */
   mpq_t*                lhs;                /**< left hand sides of constraints */
   mpq_t*                rhs;                /**< right hand sides of constraints */
   int                   nnonz;              /**< number of nonzero elements in the constraint matrix */
   int*                  beg;                /**< start index of each constraint in ind- and val-array */
   int*                  ind;                /**< indices of variables (probindex) corresponding to nonzero matrix entries */
   mpq_t*                val;                /**< values of nonzero constraint matrix entries */
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

/* just for debugging ????????????? */
#ifdef DETAILED_DEBUG /* just for debugging ?????????? */
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

static
SCIP_RETCODE printCurrentLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_ROW** rows;
   int nrows;
   int i;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   
   SCIPdebugMessage("current inexact LP rows:\n");

   for( i = 0; i < nrows; ++i )
   {
      SCIP_CALL( SCIPprintRow(scip, rows[i], NULL) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE checkLoadState(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata       /**< exactlp constraint handler data */
   )
{
   int i;
   int* cstat;
   int* rstat;
   int* cstatex;
   int* rstatex;
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

   SCIP_CALL( SCIPlpiGetBase(scip->lp->lpi, cstat, rstat) );
   SCIP_CALL( SCIPlpiexGetBase(conshdlrdata->lpiex, cstatex, rstatex) );

#ifdef BASIS_OUT /* ???????????????? */
   //printCurrentLP(scip);  /* ????????????? */
   for( i = 0; i < ncols; ++i)
   {
      if( cstat[i] == SCIP_BASESTAT_LOWER)
         printf(" c=%d: %d (at lower)\n", i, cstat[i]);
      else if( cstat[i] == SCIP_BASESTAT_BASIC)
         printf(" c=%d: %d (basic)\n", i, cstat[i]);
      else if( cstat[i] ==  SCIP_BASESTAT_UPPER)
         printf(" c=%d: %d (at upper)\n", i, cstat[i]);
      else if( cstat[i] ==  SCIP_BASESTAT_ZERO)
         printf(" c=%d: %d (at zero)\n", i, cstat[i]);
      else 
         printf(" c=%d: %d (error)\n", i, cstat[i]);
   }
   for( i = 0; i < nrows; ++i)
   {
      if( rstat[i] == SCIP_BASESTAT_LOWER)
         printf(" r=%d: %d (at lower)\n", i, rstat[i]);
      else if( rstat[i] == SCIP_BASESTAT_BASIC)
         printf(" r=%d: %d (basic)\n", i, rstat[i]);
      else if( rstat[i] ==  SCIP_BASESTAT_UPPER)
         printf(" r=%d: %d (at upper)\n", i, rstat[i]);
      else 
         printf(" r=%d: %d (error)\n", i, rstat[i]);
   }
#endif
 

   assert(ncols == ncolsex);
   for( i = 0; i < ncolsex; ++i)
   {
      assert(cstatex[i] == cstat[i]);
#ifdef BASIS_OUT /* ???????????????? */
      if( cstatex[i] == SCIP_BASESTAT_LOWER)
         printf("ex c=%d: %d (at lower)\n", i, cstatex[i]);
      else if( cstatex[i] == SCIP_BASESTAT_BASIC)
         printf("ex c=%d: %d (basic)\n", i, cstatex[i]);
      else if( cstatex[i] ==  SCIP_BASESTAT_UPPER)
         printf("ex c=%d: %d (at upper)\n", i, cstatex[i]);
      else if( cstatex[i] ==  SCIP_BASESTAT_ZERO)
         printf("ex c=%d: %d (at zero)\n", i, cstatex[i]);
      else 
         printf("ex c=%d: %d (error)\n", i, cstatex[i]);
#endif
   }

   assert(nrows == nrowsex);
   for( i = 0; i < nrowsex; ++i)
   {
      assert(rstatex[i] == rstat[i]);

#ifdef BASIS_OUT /* ???????????????? */
      if( rstatex[i] == SCIP_BASESTAT_LOWER)
         printf("ex r=%d: %d (at lower)\n", i, rstatex[i]);
      else if( rstatex[i] == SCIP_BASESTAT_BASIC)
         printf("ex r=%d: %d (basic)\n", i, rstatex[i]);
      else if( rstatex[i] ==  SCIP_BASESTAT_UPPER)
         printf("ex r=%d: %d (at upper)\n", i, rstatex[i]);
      else 
         printf("ex r=%d: %d (error)\n", i, rstatex[i]);
#endif
   }

#if 0 /*?????????????????*/
   {
      SCIP_LPI* lpi;
      SCIPgetLPI(scip, &lpi); /*????????????????*/
      SCIP_CALL( SCIPlpiWriteState(lpi, "testset/state_lp.bas") ); /*????????????????*/
      SCIP_CALL( SCIPlpiexWriteState(conshdlrdata->lpiex, "testset/state_lpex.bas") ); /* ????????????*/
      SCIPABORT(); /*?????????????*/
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
   mpq_t                 val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsNegInfinity(conshdlrdata->lpiex, val);
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
SCIP_Bool isPosInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   mpq_t                 val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsPosInfinity(conshdlrdata->lpiex, val);
}

/** converts given rational number into an FP number; uses given rounding mode during conversion 
 * (should be used to construct an FP relaxation of a constraint) 
 */
SCIP_Real mpqGetRealRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val,                /**< given rational number */
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
#ifdef READER_OUT  /*???????????????*/
      if( roundmode == GMP_RNDU )
      {
         gmp_printf("given <%Qd> <=! conv <%Qd | ", val, result); 
         printf(" %g>\n", mpq_get_d(result)); 
      }
      if( roundmode == GMP_RNDD )
      {
         gmp_printf("given <%Qd> >=! conv <%Qd | ", val, result); 
         printf(" %g>\n", mpq_get_d(result)); 
      }
#endif
      mpq_clear(result);
   }
#endif

   /* todo: check whether this way to treat infinity is ok (in particular if we want to construct a relaxation) ???????? */
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
   mpq_t                 val                 /**< given rational number */
   )
{
   SCIP_Real valapprox;

   valapprox = (SCIP_Real) mpq_get_d(val);
 
   /* todo: check whether this way to treat infinity is ok (in particular if we want to construct a relaxation) ????????? */
   if( SCIPisInfinity(scip, valapprox) )
      valapprox = SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -valapprox) )
      valapprox = -SCIPinfinity(scip);

   return valapprox;
}

/** creates constaint handler data for exactlp constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   QSexactStart(); /* todo: find a good way/place to call this method ?????????????? */

   /* open exact LP Solver interface */
   SCIP_CALL( SCIPlpiexCreate(&(*conshdlrdata)->lpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

   mpq_init((*conshdlrdata)->posinfinity);
   mpq_init((*conshdlrdata)->neginfinity);

   mpq_set((*conshdlrdata)->posinfinity, *SCIPlpiexPosInfinity((*conshdlrdata)->lpiex));
   mpq_set((*conshdlrdata)->neginfinity, *SCIPlpiexNegInfinity((*conshdlrdata)->lpiex));

   (*conshdlrdata)->lpexconstructed = FALSE;
   (*conshdlrdata)->lastenfopsnode = NULL;

   /* get event handler for updating bounds of variables in the exact LP */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for exactlp constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   if( (*conshdlrdata)->lpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*conshdlrdata)->lpiex) );
   }
   assert((*conshdlrdata)->lpiex == NULL);

   mpq_clear((*conshdlrdata)->posinfinity);
   mpq_clear((*conshdlrdata)->neginfinity);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

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
   
   /* get transforemd problem variables */
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

#ifdef DETAILED_DEBUG /*????????????????*/
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
   SCIP_VAR** vars; 
   int i;
   
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(consdata->nvars == 0 || consdata->eventdatas != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   assert(SCIPgetNVars(scip) == consdata->nvars);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPdropVarEvent(scip, vars[i], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, consdata->eventdatas[i], 
            consdata->eventdatas[i]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdatas[i]) );
   }

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

/** creates exactlp constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_OBJSEN           objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,              /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each variable in ind- and val-array */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val                 /**< values of nonzero constraint matrix entries */
   )
{
   int j;

   assert(consdata != NULL);
   assert(nvars > 0 || (nconss == 0 && nconss == 0) );
   assert(nconss > 0 || (nvars >= 0 && nnonz == 0) );
   assert(nnonz > 0 || (nconss >= 0 && nvars >= 0) );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   
   /* store variable specific information */ 
   if( nvars > 0 )
   {
      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->obj, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lb, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ub, nvars) );

      for( j = 0; j < nvars; ++j )
      {
         mpq_init((*consdata)->obj[j]);
         mpq_init((*consdata)->lb[j]);
         mpq_init((*consdata)->ub[j]);

         mpq_set((*consdata)->obj[j], obj[j]);
         mpq_set((*consdata)->lb[j], lb[j]);
         mpq_set((*consdata)->ub[j], ub[j]);
      }
   }
   else
   {
      (*consdata)->obj = NULL;
      (*consdata)->lb = NULL;
      (*consdata)->ub = NULL;
   }
   
   /* store constraint specific information */ 
   if( nconss > 0 )
   {
      /* allocate and copy integer array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->beg, beg, nconss+1) );

      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lhs, nconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->rhs, nconss) );

      for( j = 0; j < nconss; ++j )
      {
         mpq_init((*consdata)->lhs[j]);
         mpq_init((*consdata)->rhs[j]);

         mpq_set((*consdata)->lhs[j], lhs[j]);
         mpq_set((*consdata)->rhs[j], rhs[j]);
      }
   }
   else
   {
      (*consdata)->lhs = NULL;
      (*consdata)->rhs = NULL;
      (*consdata)->beg = NULL;
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
   (*consdata)->nvars = nvars;
   (*consdata)->nconss = nconss;
   (*consdata)->nnonz = nnonz;
   (*consdata)->rows = NULL;
   (*consdata)->eventdatas = NULL;
   (*consdata)->nbndchglb = 0;
   (*consdata)->nbndchgub = 0;

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
         mpq_init((*consdata)->bndchglb[j]);
         mpq_init((*consdata)->bndchgub[j]);
      }

      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgubind, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglbind, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgub, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglb, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchgubpos, (*consdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->bndchglbpos, (*consdata)->nvars);
   }

   /* free matrix specific information */ 
   for( j = 0; j < (*consdata)->nnonz; ++j )
   {
      mpq_clear((*consdata)->val[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ind, (*consdata)->nnonz);

   /* free constraint specific information */ 
   for( j = 0; j < (*consdata)->nconss; ++j )
   {
      mpq_clear((*consdata)->rhs[j]);
      mpq_clear((*consdata)->lhs[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->rhs, (*consdata)->nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lhs, (*consdata)->nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->beg, (*consdata)->nconss+1);

   /* free variable specific information */ 
   for( j = 0; j < (*consdata)->nvars; ++j )
   {
      mpq_clear((*consdata)->obj[j]);
      mpq_clear((*consdata)->lb[j]);
      mpq_clear((*consdata)->ub[j]);
   }
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
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   i                   /**< index of constraint in exactlp datas tructures */    
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
   if( consdata->beg[i] == consdata->beg[i+1] )
      SCIPinfoMessage(scip, file, "0 ");
   else
   {
      for( v = consdata->beg[i]; v < consdata->beg[i+1]; v++ )
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

/** prints exactlp constraint in CIP format to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   char s[SCIP_MAXSTRLEN];
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

         /* name */
         SCIPmessageFPrintInfo(file, "    <%s>:", SCIPvarGetName(vars[i]));
         
         /* objective value */
         gmp_snprintf(s, SCIP_MAXSTRLEN, " obj=%Qd", consdata->obj[i]);
         SCIPinfoMessage(scip, file, s);

         /* bounds (global bounds for transformed variables, original bounds for original variables) */
         SCIPmessageFPrintInfo(file, ", bounds=");

         if( isPosInfinity(conshdlrdata, consdata->lb[i]) )
            SCIPmessageFPrintInfo(file, "[+inf,");
         else if( isNegInfinity(conshdlrdata, consdata->lb[i]) )
            SCIPmessageFPrintInfo(file, "[-inf,");
         else
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "[%Qd,", consdata->lb[i]);
            SCIPinfoMessage(scip, file, s);
         }

         if( isPosInfinity(conshdlrdata, consdata->ub[i]) )
            SCIPmessageFPrintInfo(file, "+inf]");
         else if( isNegInfinity(conshdlrdata, consdata->ub[i]) )
            SCIPmessageFPrintInfo(file, "-inf]");
         else
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "%Qd]", consdata->ub[i]);
            SCIPinfoMessage(scip, file, s);
         }

         SCIPmessageFPrintInfo(file, "\n");
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
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   char s[SCIP_MAXSTRLEN];
   mpq_t activity;
   mpq_t solval;
   mpq_t prod;
   mpq_t violation;
   int c;
   int i;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking exactlp constraint <%s> for feasibility of solution %p\n",
      SCIPconsGetName(cons), (void*)sol);

   *violated = FALSE;

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
   
   /* check each linear constraint for exact feasibility */
   for( c = 0; c < consdata->nconss && !(*violated); ++c )
   {
      mpq_set_d(activity, 0.0);
  
      /* calculate activity of current constraint */
      for( i = consdata->beg[c]; i < consdata->beg[c+1]; ++i )
      {
         assert(consdata->ind[i] >= 0 && consdata->ind[i] < SCIPgetNVars(scip));

         /** @todo: ????????????
          *  - This only works if presolving is disabled (solval may already be an approximation since 
          *    solution values of aggregated variables are calculated in floating point arithmetic in SCIPgetSolVal()) 
          */ 
         mpq_set_d(solval, SCIPgetSolVal(scip, sol, vars[consdata->ind[i]]));
         
         mpq_mul(prod, consdata->val[i], solval);
         mpq_add(activity, activity, prod);
      }
      
      /* current constraint is violated */
      if( mpq_cmp(activity, consdata->rhs[c]) > 0 || mpq_cmp(activity, consdata->lhs[c]) < 0 )
      { 
         *violated = TRUE;

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

/** creates LP rows corresponding to exactlp constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR** rowvars; 
   SCIP_Real* rowvals; 
   SCIP_Real rowlhs;
   SCIP_Real rowrhs;
   int nrowvars;
   
   char rowname[SCIP_MAXSTRLEN];
   int nrows;
   int c;
   int i;

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   
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
   
#ifdef LP_OUT /* only for debugging ???????????????? */
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); 
#endif

   if( SCIPuseFPRelaxation(scip) ) 
   {
     /* for each row of the exactlp constraint, create a row with FP data that defines a relaxation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         int v;

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealRelax(scip, consdata->lhs[c], GMP_RNDD);
         rowrhs = mpqGetRealRelax(scip, consdata->rhs[c], GMP_RNDU);

         if( !SCIPisInfinity(scip, -rowlhs) && !SCIPisInfinity(scip, rowrhs) )
         {
            /* todo: create two rows and modify consdata structures in order to be able to transform a basis 
             * to the original problem later ??????? 
             */
            SCIPerrorMessage("consinitlp: for ranged rows, creating a FP relaxation is not supported yet\n");
            SCIPABORT(); /*lint --e{527}*/
            goto TERMINATE;
         }            
        
         v = 0;
         nrowvars = consdata->beg[c+1] - consdata->beg[c];
         assert(nrowvars >= 0 && nrowvars <= SCIPgetNVars(scip));

         /* calculate coefficients of all variables in the row */
         for( i = consdata->beg[c]; i < consdata->beg[c+1]; ++i )
         {
            int probidx;

            probidx = consdata->ind[i];
            
            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbLocal(vars[probidx]));
            assert(mpq_cmp(consdata->lb[probidx], consdata->ub[probidx]) <= 0);

            /* x_j >= 0 holds always ==> 
             *   cons with rhs: underestimate(a_j) * x_j <= a_j * x_j <= rhs
             *   cons with lhs:  overestimate(a_j) * x_j >= a_j * x_j >= lhs 
             */
            if( mpq_sgn(consdata->lb[probidx]) >= 0 )
            {
               if( !SCIPisInfinity(scip, rowrhs) )
               {
                  assert(SCIPisInfinity(scip, -rowlhs));
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
                  assert(SCIPisInfinity(scip, -rowlhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -rowlhs) && SCIPisInfinity(scip, rowrhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD);
               }
            }
            /* x_j <= 0 and x_j >= 0 may hold ==> 
             *   split x_j into negative and positive part 
             */
            else
            {
               /* todo: split variable into positive and negative part and modify consdata structures in order to be able 
                * to transform a basis to the original problem later ??????? 
                */
               SCIPerrorMessage("consinitlp: for variable that are neither nonnegative nor nonpositive, creating a FP relaxation is not supported yet\n");
               SCIPABORT(); /*lint --e{527}*/
               goto TERMINATE;
            }

            rowvars[v] = vars[probidx];
            v++;
         }
         
         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_relax_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons), 
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], nrowvars, rowvars, rowvals) );

#ifdef LP_OUT /* only for debugging ???????????????? */
         SCIPdebug(SCIPprintRow(scip, consdata->rows[c], NULL)); 
#endif
      }
   }
   else
   {
      /* for each row of the exactlp constraint, create a row with FP data that defines an approximation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         int v;

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealApprox(scip, consdata->lhs[c]);
         rowrhs = mpqGetRealApprox(scip, consdata->rhs[c]);

         v = 0;
         nrowvars = consdata->beg[c+1] - consdata->beg[c];
         assert(nrowvars >= 0 && nrowvars <= SCIPgetNVars(scip));

         /* add all variables to the row */
         for( i = consdata->beg[c]; i < consdata->beg[c+1]; ++i )
         {
            int probidx;
            probidx = consdata->ind[i];

            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbLocal(vars[probidx]));
            
            rowvals[v] = mpqGetRealApprox(scip, consdata->val[i]);
            rowvars[v] = vars[probidx];
            v++;
         }

         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_approx_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons), 
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], nrowvars, rowvars, rowvals) );
         
#ifdef LP_OUT /* only for debugging ???????????????? */
         SCIPdebug(SCIPprintRow(scip, consdata->rows[c], NULL)); 
#endif
      }
   }

 TERMINATE:
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
         /* todo: check whether it is ok/better not to force the cut to enter the LP ???????????? */
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], TRUE) ); 
      }
   }

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
      /* stores LP state (like basis information) of inexact LP solver into LP state object */
      /* todo: check whether this gives the basis corresponding to the node from which the unbounded ray was 
       *       found in case of an infeasible LP ???????????? 
       */
      SCIP_CALL( SCIPgetLPState(scip, &lpistate) ); 
      
      /* loads LP state (like basis information) into exact LP solver */
      SCIP_CALL( SCIPlpiexSetState(conshdlrdata->lpiex, SCIPblkmem(scip), lpistate) );
      
#ifdef DETAILED_DEBUG
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif
   }
#ifdef BASIS_OUT
   else
   {
      SCIPdebugMessage("loading LP basis into LPEX failed because dimensions don't match (LP: %d x %d, LPEX: %d x %d)\n",
         SCIPgetNLPRows(scip), SCIPgetNLPCols(scip), nrowsex, ncolsex);
   }
#endif
   return SCIP_OKAY;
}

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
      
#ifdef DETAILED_DEBUG
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif

      assert(*success);
   }
   else
   {
      *success = FALSE;

#ifdef BASIS_OUT
      SCIPdebugMessage("loading LP basis into LPEX failed because dimensions don't match (LP: %d x %d, LPEX: %d x %d)\n",
         SCIPgetNLPRows(scip), SCIPgetNLPCols(scip), nrowsex, ncolsex);
#endif
   }

   return SCIP_OKAY;
}

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

   SCIPdebugMessage("constructing initial exact LP\n");

   /* exact LP is constructed when it is needed for the first time */
   if( !conshdlrdata->lpexconstructed )
   {
      SCIP_VAR** vars;
      char** colnames; /* todo: is this implemented in a correct way ????????????*/
       
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
         assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[i]) == SCIPvarGetObj(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD) <= SCIPvarGetLbGlobal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU) >= SCIPvarGetUbGlobal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD) <= SCIPvarGetLbLocal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU) >= SCIPvarGetUbLocal(vars[i]));
      }

      /* add all columns to the exact LP */
      /* todo: check whether I implement SCIPlpiexAddCols() correctly (handling of case: beg=ind=val=NULL) ????????????? */
      SCIP_CALL( SCIPlpiexAddCols(conshdlrdata->lpiex, consdata->nvars, consdata->obj, consdata->lb, consdata->ub, colnames,
            0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(conshdlrdata->lpiex, consdata->nconss, consdata->lhs, consdata->rhs, NULL, 
            consdata->nnonz, consdata->beg, consdata->ind, consdata->val) );

      //SCIP_CALL( SCIPlpiexWriteLP(conshdlrdata->lpiex, "testset/debug2.lp") ); /* ????????????*/
      //SCIPABORT(); /*lint --e{527}*/ /* ????????????*/

      conshdlrdata->lpexconstructed = TRUE;

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &colnames);
   }

   /* apply and delete unprocessed lower bound changes */
   assert(consdata->eventdatas != NULL);
   if( consdata->nbndchglb > 0 )
   {
      SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchglb, consdata->bndchglbind, 
            consdata->bndchglb, NULL) );

      /* delete unprocessed lower bound change information */
      for( i = 0; i < consdata->nvars; ++i )
         consdata->bndchglbpos[i] = -1;
      consdata->nbndchglb = 0;
   }

   /* apply and delete unprocessed upper bound changes */
   if( consdata->nbndchgub > 0 )
   {
      SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchgub, consdata->bndchgubind, 
            NULL, consdata->bndchgub) );

      /* delete unprocessed upper bound change information */
      for( i = 0; i < consdata->nvars; ++i )
         consdata->bndchgubpos[i] = -1;
      consdata->nbndchgub = 0;
   }

   return SCIP_OKAY;
}

/** checks whether primal solution of exact LP solver satisfies all integrality restrictions */
static
SCIP_RETCODE checkIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
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
   SCIP_Bool fpvalue;
   SCIP_Bool integral;
   int branchvar;
#ifdef DETAILED_DEBUG /*????????? */
   SCIP_Real oldlb;
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

   /* get all problem variables and integer region in vars array */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbin, &nint, NULL, NULL) );
   
   integral = TRUE;
   branchvar = -1;
   /* check whether primal solution satisfies all integrality restrictions */
   for( v = 0; v < nbin + nint && integral; ++v )
   {
      assert(SCIPvarGetProbindex(vars[v]) == v);
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER );

      if( mpz_get_si(mpq_denref(primsol[v])) != 1 ) 
      {
         integral = FALSE;
         branchvar = v;
         SCIPdebugMessage("   exact LP value of intvar %s ", SCIPvarGetName(vars[v]));
#ifndef NDEBUG
         gmp_printf("<%Qd> is not integral\n", primsol[v]); 
#endif
      }
   }

   /* try to improve current local lower bound by using the result of the exact LP solver */
#ifdef DETAILED_DEBUG /*????????? */
   oldlb = SCIPgetLocalLowerbound(scip);
#endif
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, lpobjval, GMP_RNDD)) ); /* todo: check whether it is ok to use this function instead of SCIPupdateLocalDualbound() ?????????? */ 
#ifdef DETAILED_DEBUG /*????????? */
   if( oldlb < SCIPgetLocalLowerbound(scip) )
   {
      SCIPdebugMessage("lower bound improved: %.50f --> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
   }
   else
   {
      SCIPdebugMessage("lower bound did not improve: %.50f -/-> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
   }
#endif

   /* store exact LP solution as feasible MIP solution and cut off current node */
   if( integral )
   {
      SCIP_SOL* sol;
      mpq_t tmp;

      assert(branchvar == -1);

      mpq_init(tmp);

      /* create primal solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
      
      fpvalue = TRUE;

      /* check whether it is likely that the lpobjval is treated as infinite in scip, i.e., is out of range in scip */
      if( SCIPisInfinity(scip, REALABS(mpqGetRealApprox(scip, lpobjval))) )
         fpvalue = FALSE;
      
      for( v = 0; v < nvars && fpvalue; ++v )
      {
         SCIP_Real scipsolval;

         scipsolval = mpqGetRealApprox(scip, primsol[v]);

         /* check whether we can stored the primal solution without FP errors */
         mpq_set_d(tmp, scipsolval);
         if( !mpq_equal(primsol[v], tmp) || SCIPisInfinity(scip, REALABS(mpqGetRealApprox(scip, primsol[v]))) )
            fpvalue = FALSE;

         /* set solution value */
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[v], scipsolval) );
      }
      mpq_set_d(tmp, SCIPgetSolTransObj(scip, sol));
      assert(!fpvalue || mpq_cmp(lpobjval, tmp) <= 0);

      /* @todo (later): ????????? 
       * - maybe improve objval
       * - store all primal solutions in the exactlp conshdlr with the exact objval and 
       *   hand over an FP approximation of the exact sol with safe objval to scip (step 3c)
       */
      if( !fpvalue )   
         SCIPfreeSol(scip, &sol);
      else   
      {
         SCIP_Bool stored;
         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         assert(stored); /* todo: check whether this assert is correct ????????? */
      }

      *result = SCIP_CUTOFF;

      mpq_clear(tmp);
   }
   else
   {
      SCIP_Bool solvelpagain; 

      /* if the inexact LP was claimed to be infeasible, i.e., a pseudo solution is enforced here, but the exact LP 
       * has been proved to be feasible, try to solve the inexact LP again with starting basis equal to optimal 
       * basis of the exact LP; in order to avoid an infinite loop, this is only done once at a node 
       */
      solvelpagain = FALSE;
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE && conshdlrdata->lastenfopsnode != SCIPgetCurrentNode(scip) )
      {
         SCIP_CALL( loadLPEXState(scip, conshdlrdata, &solvelpagain) );

         if( solvelpagain )
         {
            solvelpagain = TRUE;
            *result = SCIP_SOLVELP;
         }
      }

      /* branch on the first fractional variable in the exact LP solution */
      if( !solvelpagain )
      {
         SCIP_NODE* node;
         SCIP_Real downub;
         SCIP_Real uplb;

         /* @todo (later): ???????????
          *  - implement more sophisticated branching rules 
          */

         assert(branchvar >= 0 && branchvar < nvars);
      
         /* create left child: add x_i <= floor(x_i^*) */
         downub = floor(mpqGetRealRelax(scip, primsol[branchvar], GMP_RNDD));
         SCIP_CALL( SCIPcreateChild(scip, &node, SCIPcalcNodeselPriority(scip, vars[branchvar], downub),
               SCIPcalcChildEstimate(scip, vars[branchvar], downub)) );
         SCIP_CALL( SCIPchgVarUbNode(scip, node, vars[branchvar], downub) );
      
         /* create right child: add x_i >= ceil(x_i^*) */
         uplb = ceil(mpqGetRealRelax(scip, primsol[branchvar], GMP_RNDU)); 
         SCIP_CALL( SCIPcreateChild(scip, &node, SCIPcalcNodeselPriority(scip, vars[branchvar], uplb),
               SCIPcalcChildEstimate(scip, vars[branchvar], uplb)) );
         SCIP_CALL( SCIPchgVarLbNode(scip, node, vars[branchvar], uplb) );
 
         assert(uplb = downub + 1.0);

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
      }
   }

   /* free temporary arrays */
   for( v = 0; v < ncols; ++v )
   {
      mpq_clear(primsol[v]);
   }
   mpq_clear(lpobjval);
   SCIPfreeBufferArray(scip, &primsol);

   if( !fpvalue )   
   {
      SCIPerrorMessage("storing optimal solutions of subproblems that are not FP representable is not supported yet\n");
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** solves the exact LP with the given algorithm and evaluates return status */
static
SCIP_RETCODE solveLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
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

#ifdef EXLPSOLVER_OUT /* only for debugging ???????????????? */
   SCIP_CALL( SCIPlpiexSetIntpar(conshdlrdata->lpiex, SCIP_LPPAR_LPINFO, TRUE) ); 
#endif

   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      SCIPdebugMessage("solving current primal exact LP (%d cols, %d rows):\n", ncols, nrows);

      /* call primal simplex */
      retcode = SCIPlpiexSolvePrimal(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("   (node %"SCIP_LONGINT_FORMAT") primal simplex solving error in current exact LP\n", 
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }
 
      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("   solved primal exact LP in %d iterations\n", iterations);
      break;

   case SCIP_LPALGO_DUALSIMPLEX:
      SCIPdebugMessage("solving current dual exact LP (%d cols, %d rows)\n", ncols, nrows);
      
      /* call dual simplex */
      retcode = SCIPlpiexSolveDual(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("   (node %"SCIP_LONGINT_FORMAT") dual simplex solving error in current exact LP\n", 
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }

      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("   solved dual exact LP in %d iterations\n", iterations);
      break;

   default:
      SCIPerrorMessage("invalid exact LP algorithm\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** evaluates the result of the exact LP */
static
SCIP_RETCODE evaluateLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_RESULT*          result              /**< pointer to store the result of the lp enforcement call */
   )
{
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
      checkIntegrality(scip, conshdlrdata, result);
      assert(*result == SCIP_CUTOFF || SCIP_BRANCHED);

   }
   else if( SCIPlpiexIsObjlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsPrimalInfeasible(conshdlrdata->lpiex) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n"); 
      
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiexExistsPrimalRay(conshdlrdata->lpiex) ) /* todo: check why in lp.c SCIPlpiIsPrimalUnbounded() is not used ????*/
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsIterlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsTimelimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
      SCIPABORT();
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
 * Callback methods of constraint handler
 */

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
#define consInitExactlp NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitExactlp NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreExactlp NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreExactlp NULL
#endif


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
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolExactlp NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

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

   /*debugMessage("Trans method of exactlp constraints\n");*/

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

   /* create exactlp constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr, sourcedata->objsense, sourcedata->nvars, 
         sourcedata->obj, sourcedata->lb, sourcedata->ub, sourcedata->nconss, sourcedata->lhs, sourcedata->rhs, 
         sourcedata->nnonz, sourcedata->beg, sourcedata->ind, sourcedata->val) );
 
   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

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
/* todo: implement this (because of sepafreq = 1, there will be an error ??????? */
static
SCIP_DECL_CONSSEPALP(consSepalpExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN; /* todo: find out which result should be returend ?????????*/

   /* in case the FP problem is a relaxation of the original problem, we have already calculated a proved lower bound
    * via postprocessing the LP solution of the FP problem 
    */
   if( SCIPuseFPRelaxation(scip) ) 
      return SCIP_OKAY;

   /* todo: CONTINUE here ?????????
    * - fill lpiex
    * - calculate safe lower bound and hand it to SCIP
    */
   //   SCIP_CALL( constructCurrentLPEX(scip, conshdlrdata, consdata) );
   
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
   char algo;


   assert(SCIPhasCurrentNodeLP(scip));
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIPdebugMessage("enforcing exactlp constraint <%s> on LP solution\n", SCIPconsGetName(conss[0]));

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_INFEASIBLE;

   /* constructs exact LP of current node */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   /* load LP state from inexact LP into exact LP solver */
   SCIP_CALL( loadLPState(scip, conshdlrdata) );

   /* solve exact LP */
   algo = 'd'; /* todo: what should be do here (selection by warmstartinfo is not useful, algo used in inexact LP solver better?) ?????????? */
   switch( algo )
   {
   case 'd':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, &lperror) );
      break;

   case 'p':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_PRIMALSIMPLEX, &lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for exact LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }

   if( lperror )
   {
      SCIPerrorMessage("exact LP solver returns error: case not handled yet\n");
      SCIPABORT();
   }

   /* evaluate result of exact LP solver */
   SCIP_CALL( evaluateLPEX(scip, conshdlrdata, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool lperror; 
   char algo;

   SCIPdebugMessage("enforcing exactlp constraint <%s> on pseudo solution (LP solstat=%d)\n", SCIPconsGetName(conss[0]), 
      SCIPgetLPSolstat(scip));

   assert(SCIPisExactSolve(scip));
   assert(!objinfeasible);

   assert(!SCIPhasCurrentNodeLP(scip));
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE 
      || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED
      || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_ERROR);

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_INFEASIBLE;

   /* constructs exact LP of current node */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   /* load LP state from inexact LP into exact LP solver */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      SCIP_CALL( loadLPState(scip, conshdlrdata) );
   }

   /* solve exact LP */
   algo = 'd'; /* todo: what should be do here (selection by warmstartinfo is not useful, algo used in inexact LP solver better?) ?????????? */
   switch( algo )
   {
   case 'd':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, &lperror) );
      break;

   case 'p':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_PRIMALSIMPLEX, &lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for exact LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }

   if( lperror )
   {
      SCIPerrorMessage("exact LP solver returns error: case not handled yet\n");
      SCIPABORT();
   }

   /* evaluate result of exact LP solver */
   SCIP_CALL( evaluateLPEX(scip, conshdlrdata, result) );

   SCIPdebugMessage(" -> enforcing pseudo solution returned result <%d>\n", *result);

   /* todo (later): ????????????????? 
    *  - decide whether it is possible and useful to support working on a pseudo solution 
    *    (inexact from scip or exact from exactlp conshdlr) in addition 
    */
   
   /* remember that current node is the one at which a pseudo solution was enforced last */
   conshdlrdata->lastenfopsnode = SCIPgetCurrentNode(scip);

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExactlp)
{  /*lint --e{715}*/
   SCIP_Bool violated;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);
   assert(sol != NULL);
   assert(checklprows); /* in exact mode, constraint handler cannot draw conclusions from the LP's (inexact) feasibility */
   
   /*debugMessage("Check method of exactlp constraint\n");*/

   SCIP_CALL( checkCons(scip, SCIPconshdlrGetData(conshdlr), conss[0], sol, printreason, &violated) );

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;
   
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
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolExactlp NULL
#endif


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
   SCIP_Bool* roundinguplocked;
   SCIP_Bool* roundingdownlocked;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int j;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get problem variables; in transforming stage we can only access original variables, 
    * however, SCIPaddVarLocks() will use the transformed variable if it exists 
    */
   vars = SCIPgetOrigVars(scip); 
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);

   /* allocate and initialize temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &roundingdownlocked, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roundinguplocked, consdata->nvars) );
   BMSclearMemoryArray(roundingdownlocked, consdata->nvars);
   BMSclearMemoryArray(roundinguplocked, consdata->nvars);

   /* for every variable, check whether rounding up or down could lead to the violation of some constraint */
   j = 0; 
   for( c = 0; c < consdata->nconss; ++c )
   {
      haslhs = !isNegInfinity(conshdlrdata, consdata->lhs[c]);
      hasrhs = !isPosInfinity(conshdlrdata, consdata->rhs[c]);

      /* check all variables of the current constraint */
      for( ; j < consdata->beg[c+1]; ++j )
      {
         assert(consdata->ind[j] >= 0 && consdata->ind[j] < consdata->nvars);

         /* val > 0 */
         if( mpq_sgn(consdata->val[j]) > 0 )
         {
            if( haslhs )
               roundingdownlocked[consdata->ind[j]] = TRUE;
            
            if( hasrhs )
               roundinguplocked[consdata->ind[j]] = TRUE;
         }

         /* val < 0 */
         if( mpq_sgn(consdata->val[j]) < 0 )
         {
            if( haslhs )
               roundinguplocked[consdata->ind[j]] = TRUE;
            
            if( hasrhs )
               roundingdownlocked[consdata->ind[j]] = TRUE;
         }
      } 
   }
   assert(j == consdata->nnonz);

   /* set rounding locks for all variables */
   for( j = 0; j < consdata->nvars; ++j )
   {
      int probindex;

      probindex = SCIPvarGetProbindex(vars[j]);

      if( roundingdownlocked[probindex] && roundinguplocked[probindex] )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      else if( roundingdownlocked[probindex] )
      { 
         assert(!roundinguplocked[probindex]);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos, nlocksneg) );
      } 
      else if( roundinguplocked[probindex] )
      { 
         assert(!roundingdownlocked[probindex]);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlocksneg, nlockspos) );
      } 
   }
   
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &roundinguplocked);
   SCIPfreeBufferArray(scip, &roundingdownlocked);

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
      
      /* create the bin Exactlp constraint (an automatically upgraded constraint is always unmodifiable) */
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

/** execution methode of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecExactlp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real newbound;
   int varind;

   assert(eventdata != NULL);
   assert(eventdata->consdata != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);
   assert(0 <= consdata->nbndchglb && consdata->nbndchglb < consdata->nvars);
   assert(0 <= consdata->nbndchgub && consdata->nbndchgub < consdata->nvars);

   varind = eventdata->varind;
   assert(0 <= varind && varind < consdata->nvars);
   assert(SCIPvarGetProbindex(SCIPeventGetVar(event)) == varind);

   newbound = SCIPeventGetNewbound(event);

      
   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
#ifdef BOUNDCHG_OUT /* ?????????? */
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "store unprocessed bound change (%x): <%s> [%Qd,%Qd] >= %g\n", 
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lb[varind], consdata->ub[varind], 
            newbound);
         SCIPdebugMessage(s);
      }
#endif
      /* check whether variable already has unprocessed lower bound change events */ 
      if( consdata->bndchglbpos[varind] == -1 )
      {
         /* store position of variable in unprocessed lower bound change arrays */
         consdata->bndchglbpos[varind] = consdata->nbndchglb; 

         /* store new lower bound and index of variable */
         mpq_set_d(consdata->bndchglb[consdata->nbndchglb], newbound);
         consdata->bndchglbind[consdata->nbndchglb] = varind;

         /* update number of unprocessed lower bound change events */
         consdata->nbndchglb++;
      }
      else
      {
         mpq_set_d(consdata->bndchglb[consdata->bndchglbpos[varind]], newbound);
         assert(consdata->bndchglbind[consdata->bndchglbpos[varind]] == varind);
      }
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
#ifdef BOUNDCHG_OUT /* ?????????? */
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "store unprocessed bound change (%x): <%s> [%Qd,%Qd] <= %g\n", 
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lb[varind], consdata->ub[varind], 
            newbound);
         SCIPdebugMessage(s);
      }
#endif

      /* check whether variable already has unprocessed upper bound change events */ 
      if( consdata->bndchgubpos[varind] == -1 )
      {
         /* store position of variable in unprocessed upper bound change arrays */
         consdata->bndchgubpos[varind] = consdata->nbndchgub;

         /* store new upper bound and index of variable */
         mpq_set_d(consdata->bndchgub[consdata->nbndchgub], newbound);
         consdata->bndchgubind[consdata->nbndchgub] = varind;

         /* update number of unprocessed upper bound change events */
         consdata->nbndchgub++;
      }
      else
      {
         mpq_set_d(consdata->bndchgub[consdata->bndchgubpos[varind]], newbound);
         assert(consdata->bndchgubind[consdata->bndchgubpos[varind]] == varind);
      }
      break;
   default:
      SCIPerrorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

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
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecExactlp,
         eventhdlrdata) );

   /* create exactlp constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeExactlp, consInitExactlp, consExitExactlp, 
         consInitpreExactlp, consExitpreExactlp, consInitsolExactlp, consExitsolExactlp,
         consDeleteExactlp, consTransExactlp, consInitlpExactlp,
         consSepalpExactlp, consSepasolExactlp, consEnfolpExactlp, consEnfopsExactlp, consCheckExactlp, 
         consPropExactlp, consPresolExactlp, consRespropExactlp, consLockExactlp,
         consActiveExactlp, consDeactiveExactlp, 
         consEnableExactlp, consDisableExactlp,
         consPrintExactlp,
         conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdExactlp, LINCONSUPGD_PRIORITY) );
#endif

   return SCIP_OKAY;
}

/** creates and captures a exactlp constraint */
SCIP_RETCODE SCIPcreateConsExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_OBJSEN           objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,              /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each variable in ind- and val-array */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries */
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
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, objsense, nvars, obj, lb, ub, nconss, lhs, rhs, 
         nnonz, beg, ind, val) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}
