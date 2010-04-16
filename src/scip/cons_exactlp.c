/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_exactlp.c,v 1.1.2.13 2010/04/16 15:41:14 bzfwolte Exp $"
//#define SCIP_DEBUG /*??????????????*/
//#define LP_OUT /* only for debugging ???????????????? */
//#define BOUNDCHG_OUT /* only for debugging ?????????? */
//#define PSEUDOOBJ_OUT /* only for debugging ?????????? */
//#define BASIS_OUT /* only for debugging ???????????????? */
//#define DBVERIFY_TEST /* only for testing the methode ?????????????????? */
//#define EXLPSOLVER_OUT /* only for debugging ???????????????? */
//#define DETAILED_DEBUG /* only for debugging ???????????????? */
//#define PS_OUT /* only for debugging ???????????????? */
//#define PS_OUT_2 /* only for debugging ???????????????? */
//#define USE_QSFACTOR 		/* use qsoptex for projection factorization */
//#define READER_OUT /* only for debugging ???????????????? */
//#define PRESOL_OUT /* only for debugging ?????????? */
//#define BOUNDING_OUT /* only for debugging ?????????? */

/**@file   cons_exactlp.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for exactlp constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h> 

#include "scip/cons_exactlp.h"
#include "scip/struct_lp.h" /* only for debugging ??????????*/
#include "scip/struct_scip.h" 
#include "scip/solex.h" 
#include "scip/primalex.h" 
#include "EGlib.h" 
#include "QSopt_ex.h" 
#include "scip/misc.h" 
#include "rectlu/rectlu.h"

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

#define USE_RECTFACTOR		        /**< use rectfactor for projection factorization */


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
   mpq_t                 pseudoobjval;       /**< pseudo solution value at node where all unprocessed bound changes were
                                              *   applied last with all variables set to their best bounds, ignoring 
                                              *   variables with infinite best bound */
   int                   pseudoobjvalinf;    /**< number of variables with infinite best bound in pseudo solution at node 
                                              *   where all unprocessed bound changes were applied last */
   SCIP_LPIEX*           psfactor;           /**< stores factorized matrix for project and scale */ 
   qsnum_factor_work*    rectfactor;         /**< stores factorized matrix for project and scale */
   mpq_t*                interiorpt;         /**< stores relative interior point for root node dual problem */
   int*                  impliedeq;          /**< 1 if constraints dual variable is always zero, 0 otherwise */
   int                   psdim;              /**< dimension of interior point */
   mpq_t                 commonslack;        /**< slack by which interior point satisfies all dual inequalities */
   SCIP_PRIMALEX*        primal;             /**< exact primal data and solution storage */
   SCIP_Bool             psdatacon;          /**< was project and scale data structure constructed? */  
   SCIP_Bool             psfactorfail;       /**< did factorization of projection matrix fail? */
};

/** constraint data for exactlp constraints */
struct SCIP_ConsData
{
   SCIP_OBJSEN           objsense;           /**< objective sense */
   mpq_t                 objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   mpq_t                 objscale;           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   int                   nvars;              /**< number of variables */
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
   int                   conssize;           /**< size of constraint specific array in exactlp data structure */
   mpq_t*                lhs;                /**< left hand sides of constraints */
   mpq_t*                rhs;                /**< right hand sides of constraints */
   int                   nnonz;              /**< number of nonzero lements in the constraint matrix */
   int*                  beg;                /**< start index of each constraint in ind and val array */
   int*                  len;                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind;                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val;                /**< values of nonzero constraint matrix entries (and some zeros) */
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
   //   printCurrentLP(scip);  /* ????????????? */
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

#ifndef NDEBUG
/** checks whether pseudo objective value has been updated correctly */
static
void checkPseudoobjval(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,      /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata           /**< constraint data */
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
 
   for( i = 0; i < consdata->nvars; ++i )
   {
      /* add new pseudo objective value */
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
   const mpq_t           val                 /**< given rational number */
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

#ifdef DETAILED_DEBUG /*????????????????*/
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

   QSexactStart(); /* todo: find a good way/place to call this method ?????????????? */

   /* open exact LP Solver interface */
   SCIP_CALL( SCIPlpiexCreate(&(*conshdlrdata)->lpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

   /* open exact LP Solver interface for factorization */
#ifdef USE_QSFACTOR
   SCIP_CALL( SCIPlpiexCreate(&(*conshdlrdata)->psfactor, NULL, SCIP_OBJSEN_MINIMIZE) );
#endif
#ifdef USE_RECTFACTOR
   (*conshdlrdata)->rectfactor = (qsnum_factor_work*) NULL;
#endif
   (*conshdlrdata)->psdim = 0;

   mpq_init((*conshdlrdata)->commonslack);
   
   mpq_init((*conshdlrdata)->posinfinity);
   mpq_init((*conshdlrdata)->neginfinity);

   mpq_set((*conshdlrdata)->posinfinity, *SCIPlpiexPosInfinity((*conshdlrdata)->lpiex));
   mpq_set((*conshdlrdata)->neginfinity, *SCIPlpiexNegInfinity((*conshdlrdata)->lpiex));

   (*conshdlrdata)->lpexconstructed = FALSE;
   (*conshdlrdata)->psdatacon = FALSE;
   (*conshdlrdata)->lastenfopsnode = NULL;

   mpq_init((*conshdlrdata)->pseudoobjval);
   mpq_set_d((*conshdlrdata)->pseudoobjval, 0.0);
   (*conshdlrdata)->pseudoobjvalinf = 0;

   (*conshdlrdata)->interiorpt = NULL;
   (*conshdlrdata)->impliedeq = NULL;

   /* get event handler for updating bounds of variables in the exact LP */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for exactlp constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   (*conshdlrdata)->primal = NULL;
   (*conshdlrdata)->psfactorfail = FALSE;

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

   if( (*conshdlrdata)->psdim > 0 )
   {
      for( i = 0; i < (*conshdlrdata)->psdim; i++ )
         mpq_clear((*conshdlrdata)->interiorpt[i]);
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->interiorpt);
      SCIPfreeMemoryArray(scip, &(*conshdlrdata)->impliedeq);
   }

#ifdef USE_QSFACTOR
   if( (*conshdlrdata)->psfactor != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*conshdlrdata)->psfactor) );
   }
   assert((*conshdlrdata)->psfactor == NULL);
#endif

#ifdef USE_RECTFACTOR
   if( (*conshdlrdata)->rectfactor != NULL)
      RECTLUfreeFactorization((*conshdlrdata)->rectfactor);
#endif

   if( (*conshdlrdata)->lpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*conshdlrdata)->lpiex) );
   }
   assert((*conshdlrdata)->lpiex == NULL);

   mpq_clear((*conshdlrdata)->pseudoobjval);
   mpq_clear((*conshdlrdata)->commonslack);
   mpq_clear((*conshdlrdata)->posinfinity);
   mpq_clear((*conshdlrdata)->neginfinity);

   SCIPfreeMemory(scip, conshdlrdata);

   QSexactClear();/* todo: find a good way/place to call this method ?????????????? */

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
   int                   nconss,             /**< number of constraints in exactlp data structure */
   int                   conssize,           /**< size of constraint specific array in exactlp data structure */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each constraint in ind and val array */
   int*                  len,                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries (and some zeros) */
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
      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->obj, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lb, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ub, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lbloc, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ubloc, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lockdown, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lockup, nvars) );

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
      }
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
   (*consdata)->nconss = nconss;
   (*consdata)->conssize = conssize;
   (*consdata)->nnonz = nnonz;
   (*consdata)->rows = NULL;
   (*consdata)->eventdatas = NULL;
   (*consdata)->nbndchglb = 0;
   (*consdata)->nbndchgub = 0;
   (*consdata)->objneedscaling = objneedscaling;

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
   SCIPmessageFPrintInfo(file, "<%s>:", SCIPvarGetName(var));
   
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

         SCIPmessageFPrintInfo(file, "    ");
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
         /** @todo: ????????????
          *  - This only works if presolving is disabled (solval may already be an approximation since 
          *    solution values of aggregated variables are calculated in floating point arithmetic in SCIPgetSolVal()) 
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
            /** @todo: ????????????
             *  - This only works if presolving is disabled (solval may already be an approximation since 
             *    solution values of aggregated variables are calculated in floating point arithmetic in SCIPgetSolVal()) 
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

/** updates current pseudo objective values for a change in a variable's  objective value or bounds */
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

#ifdef PSEUDOOBJ_OUT /* only for debugging ?????????? */
   gmp_printf("obj: %Qd -> %Qd, lb:  %Qd -> %Qd, ub:  %Qd -> %Qd: pseudoobjval: %Qd (%d) --> ", 
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

#ifdef PSEUDOOBJ_OUT /* only for debugging ?????????? */
   gmp_printf("%Qd (%d)\n", conshdlrdata->pseudoobjval, conshdlrdata->pseudoobjvalinf); 
#endif

   assert(conshdlrdata->pseudoobjvalinf >= 0);

   mpq_clear(deltaval);
   mpq_clear(prod);

   return SCIP_OKAY;
}

/** applies all unprocessed bounds changes of variables, i.e., updates exact local bound information in constraint
 *  data, exact LP, and in project and scale data structure
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
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("apply lower bound changes to LPEX (nbndchglb=%d, nbndchgub=%d):\n", consdata->nbndchglb, 
            consdata->nbndchgub);
#endif

         SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchglb, consdata->bndchglbind, 
               consdata->bndchglb, NULL) );
      }

      /* update lower bounds in project and scale data structure */
      if( conshdlrdata->psdatacon )
      {
         /* todo Dan: update your data structures here and call processBoundchgs() in constructPSData() in order to update bounds; 
          * NOTE: Now you have to use lbloc/ubloc when you initialize your data ????????? */ 
      }

      /* update pseudo objective value and local lower bounds */
      for( i = 0; i < consdata->nbndchglb; ++i )
      {
         int probidx;
         
         probidx = consdata->bndchglbind[i];

#ifdef PSEUDOOBJ_OUT /* only for debugging ?????????? */
         printf("update lb of var with probidx %d: ", probidx); 
#endif
         if( conshdlrdata->lpexconstructed )
         {
            /* todo: currently (since it was easier to implement) the pseudoobjval is only computed/updated, 
             * if an exact LP exists; this should be more general ?????????????
             */
            SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
                  consdata->obj[probidx], consdata->bndchglb[i], consdata->ubloc[probidx]) );
         }         
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
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("apply upper bound changes to LPEX (nbndchglb=%d, nbndchgub=%d):\n", consdata->nbndchglb, 
            consdata->nbndchgub);
#endif
         SCIP_CALL( SCIPlpiexChgBounds(conshdlrdata->lpiex, consdata->nbndchgub, consdata->bndchgubind, 
               NULL, consdata->bndchgub) );
      }

      /* update upper bounds in project and scale data structure */
      if( conshdlrdata->psdatacon )
      {
         /* todo Dan: update your data structures here and call processBoundchgs() in constructPSData() in order to update bounds; 
          * NOTE: Now you have to use lbloc/ubloc when you initialize your data ????????? */ 
      }

      /* update pseudo objective value and local lower bounds */
      for( i = 0; i <  consdata->nbndchgub; ++i )
      {
         int probidx;
         
         probidx = consdata->bndchgubind[i];

#ifdef PSEUDOOBJ_OUT /* only for debugging ?????????? */
         printf("update ub of var with probidx %d: ", probidx); 
#endif
         if( conshdlrdata->lpexconstructed )
         {
            /* todo: currently (since it was easier to implement) the pseudoobjval is only computed/updated, 
             * if an exact LP exists; this should be more general ?????????????
             */
            SCIP_CALL( updateVar(conshdlrdata, consdata->obj[probidx], consdata->lbloc[probidx], consdata->ubloc[probidx],
                  consdata->obj[probidx], consdata->lbloc[probidx], consdata->bndchgub[i]) );
         }         
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
   assert(conshdlrdata->lpexconstructed && consdata->nbndchgub == 0 && consdata->nbndchglb == 0);
   
   if( conshdlrdata->pseudoobjvalinf > 0 || SCIPignorePseudosol(scip) )
      return negInfinity(conshdlrdata);
   else
      return (const mpq_t*) (&conshdlrdata->pseudoobjval);
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
   SCIP_CONSHDLRDATA* conshdlrdata;
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
               /* todo: split variable into positive and negative part and modify consdata structures in order to be able 
                * to transform a basis to the original problem later ??????? 
                */
               SCIPerrorMessage("consinitlp: for variables that are neither nonnegative nor nonpositive, creating a FP relaxation is not supported yet\n");

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
         
#ifdef LP_OUT /* only for debugging ???????????????? */
         SCIPdebug(SCIPprintRow(scip, consdata->rows[c], NULL)); 
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
         /* todo: check whether it is ok/better not to force the cut to enter the LP ???????????? */
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], TRUE) ); 
      }
   }

   return SCIP_OKAY;
}

#ifdef BASIS_OUT /* ???????????????? */
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
      
#ifdef DETAILED_DEBUG
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif

      /* frees LP state */
      SCIP_CALL( SCIPlpiexFreeState(conshdlrdata->lpiex, SCIPblkmem(scip), &lpistate) );

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
      char** colnames; /* todo: is this implemented in a correct way ????????????*/
      mpq_t tmpzero;
       
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
      /* todo: check whether I implement SCIPlpiexAddCols() correctly (handling of case: beg=ind=val=NULL) ????????????? */
      SCIP_CALL( SCIPlpiexAddCols(conshdlrdata->lpiex, consdata->nvars, consdata->obj, consdata->lbloc, consdata->ubloc, 
            colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(conshdlrdata->lpiex, consdata->nconss, (const mpq_t*) consdata->lhs, 
            (const mpq_t*) consdata->rhs, NULL, consdata->nnonz, consdata->beg, consdata->len, consdata->ind, 
            consdata->val) );

      //SCIP_CALL( SCIPlpiexWriteLP(conshdlrdata->lpiex, "testset/debug2.lp") ); /* ????????????*/
      //SCIPABORT(); /*lint --e{527}*/ /* ????????????*/

      conshdlrdata->lpexconstructed = TRUE;

      /* calculate pseudo objective value */
      mpq_init(tmpzero);
      mpq_set_d(tmpzero, 0.0);
      for( i = 0; i < consdata->nvars; ++i )
      {
         SCIP_CALL( updateVar(conshdlrdata, tmpzero, tmpzero, tmpzero, consdata->obj[i], consdata->lb[i], 
               consdata->ub[i]) );
      }
      mpq_clear(tmpzero);

#ifndef NDEBUG
      /* check whether pseudo objective value was updated correctly */
      checkPseudoobjval(scip, conshdlrdata, consdata);
#endif
     
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
      /* stores LP state (like basis information) of inexact LP solver into LP state object */
      /* todo: check whether this gives the basis corresponding to the node from which the unbounded ray was 
       *       found in case of an infeasible LP ???????????? 
       */
      SCIP_CALL( SCIPgetLPState(scip, &lpistate) ); 
      
      /* loads LP state (like basis information) into exact LP solver */
      //      printCurrentLP(scip);  /* ????????????? */
      SCIP_CALL( SCIPlpiexSetState(conshdlrdata->lpiex, SCIPblkmem(scip), lpistate) );
      
#ifdef DETAILED_DEBUG
      SCIP_CALL( checkLoadState(scip, conshdlrdata) );
#endif

      /* free LP state */
      SCIP_CALL( SCIPfreeLPState(scip, &lpistate) );
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

/** checks whether primal solution of exact LP solver satisfies all integrality restrictions */
static
SCIP_RETCODE checkIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIP_Bool fpvalue;
   SCIP_Bool inrange;
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
   
   /* check whether primal solution satisfies all integrality restrictions */
   integral = TRUE;
   branchvar = -1;
   for( v = 0; v < nbin + nint && integral; ++v )
   {
      assert(SCIPvarGetProbindex(vars[v]) == v);
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER );

      if( mpz_get_si(mpq_denref(primsol[v])) != 1 ) 
      {
         integral = FALSE;
         branchvar = v;
#ifdef DETAILED_DEBUG  /*???????????????*/
         SCIPdebugMessage("   exact LP value of intvar %s ", SCIPvarGetName(vars[v]));
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
      SCIP_SOLEX* solex;
      mpq_t tmp;
      SCIP_Bool stored;

      assert(branchvar == -1);

      mpq_init(tmp);

#ifdef DETAILED_DEBUG  /*???????????????*/
      SCIPdebugMessage("   current exact LP solution is integral\n");
#endif

      /* create exact and approximate primal solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
      SCIP_CALL( SCIPsolexCreate(&solex, SCIPblkmem(scip), NULL) );
      
      fpvalue = TRUE;
      inrange = TRUE;

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
         SCIP_CALL( SCIPsetSolVal(scip, sol, vars[v], scipsolval) );
         /* todo: using scip->set is ugly! think about a general way to get information from the constraint handler 
          *       implement all the methods in cons_exactlp.c that operate and use exact problem data, like exact solutions
          *       ??????????? 
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
         /* @todo (later): ????????? 
          * 1. maybe improve objval (using SCIPsetSolTransObj() is only a workaround, and this method is not very safe, 
          *    i.e., should be deleted from the source code later when point 2. was implemented)
          * 2. store all primal solutions in the exactlp conshdlr with the exact objval and 
          *    hand over an FP approximation of the exact sol with safe objval to scip (step 3c)
          */
         /* store safe objective value for approximate primal solution */ 
         SCIP_CALL( SCIPsetSolTransObj(scip, sol, mpqGetRealRelax(scip, lpobjval, GMP_RNDU)) );

         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
         
#if 0 /* since exact solutions can be stored now, this is not needed anymore */
         if( stored && !fpvalue )
         {
            SCIPwarningMessage("Note: Primal solution found is NOT FP representable (primal bound stored is safe, but primal solution stored is only an FP approximation)!\n");
         }
#endif
      }

      /* add exact primal solution */ 
      /* todo: using scip->set is ugly! think about a general way get information from the constraint handler 
       *       implement all the methods in cons_exactlp.c that operate and use exact problem data, like exact solutions
       *       ??????????? 
       */
      SCIP_CALL( SCIPprimalexAddSolFree(conshdlrdata->primal, SCIPblkmem(scip), scip->set, scip->transprob, 
            &solex, &stored) );

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
#ifdef DETAILED_DEBUG  /*???????????????*/
            SCIPdebugMessage("   ---> LP was infeasible but LPEX not: solve LP once again with basis of LPEX\n");
#endif
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
#ifdef DETAILED_DEBUG  /*???????????????*/
         SCIPdebugMessage("   ---> LP was not solved: branch on result of LPEX\n");
#endif
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
      SCIPerrorMessage("storing optimal solutions of subproblems that is out of FP range is not supported yet\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** evaluates the result of the exact LP */
static
SCIP_RETCODE evaluateLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
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
      SCIP_CALL( checkIntegrality(scip, conshdlrdata, consdata, result) );
      assert(*result == SCIP_CUTOFF || SCIP_BRANCHED);

   }
   else if( SCIPlpiexIsObjlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiexIsPrimalInfeasible(conshdlrdata->lpiex) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n"); 
      
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiexExistsPrimalRay(conshdlrdata->lpiex) ) /* todo: check why in lp.c SCIPlpiIsPrimalUnbounded() is not used ????*/
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
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
      return SCIP_ERROR;
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
 * local methods for generating safe dual bounds
 */

/* todo: tell Dan that length array is present now because val array can have zero entries (dummies). 
 *       he has to adapt his code ???????? 
 */  

/** constructs the exact LP for the project and scale dual bound method  */
static
SCIP_RETCODE constructPSData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   int                   i;
   int                   j;
   int                   rval;
   int                   pos;
   mpq_t                 mpqtemp;
   int                   nobjnz;             /**< number of nonzeros in cost vector */ 
   SCIP_LPIEX*           pslpiex;            /**< modified exact LP to find interior point of dual */ 
   int                   psnvars;            /**< number of variables in modified problem */
   mpq_t*                psobj;              /**< objective function values of variables in modified problem */
   mpq_t*                pslb;               /**< lower bounds of variables in modified problem */
   mpq_t*                psub;               /**< upper bounds of variables in modified problem */
   int                   psnconss;           /**< number of constraints in modified problem */
   mpq_t*                pslhs;              /**< left hand sides of constraints in modified problem */
   mpq_t*                psrhs;              /**< right hand sides of constraints in modified problem */
   int                   psnnonz;            /**< number of nonzero elements in the constraint matrix in modified prob */
   int*                  psbeg;              /**< start index of each constraint in ind- and val-array in modified prob */
   int*                  pslen;              /**< length array for modified problem */
   int*                  psind;              /**< indices of variables (probindex) corresponding to nonzero matrix entries */
   mpq_t*                psval;              /**< values of nonzero constraint matrix entries in modified problem */   

   mpq_t*                dualsol;            /**< stores the dual solution vector of the modified problem */  
   char **               colnames;           /**< meaningless variable names */
   mpq_t                 objval;             /**< objective value of modified problem */

   int*                  projbasis;          /**< dual variables modified to do projections */   
   int                   projbasislen;       /**< number of columns to use for projection */
   int*                  projbeg;            /**< projection submatrix begin */
   int*                  projlen;            /**< projection column length */
   int*                  projind;            /**< dual row indices of projection submatrix */
   mpq_t*                projval;            /**< values of projection submatrix */     

   assert(!conshdlrdata->psfactorfail);
   assert(consdata->nconss > 0);

   if( conshdlrdata->psdatacon )
      return SCIP_OKAY;

   mpq_init(mpqtemp);

   /* store dimension of interior point */
   conshdlrdata->psdim = consdata->nconss;

   /* allocate memory for the interior point solution */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->impliedeq, conshdlrdata->psdim) ); 
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->interiorpt, conshdlrdata->psdim) );

   for( i = 0; i < conshdlrdata->psdim ; i++)
      mpq_init(conshdlrdata->interiorpt[i]);

   /* figure out the dimensions of the updated problem */
   psnvars = consdata->nvars + 2 * consdata->nconss + 1;
   psnconss = 2 * consdata->nconss + 1;
   nobjnz = 0;
   for( i = 0; i < consdata->nvars; i++)
   {
      if(mpq_sgn(consdata->obj[i]))
      {
         nobjnz++;
      }
   }
   psnnonz = consdata->nnonz + 3 * consdata->nconss + 1 + nobjnz;

   /* allocate memory for new problem */
   SCIP_CALL( SCIPallocBufferArray(scip, &psobj, psnvars) );
   for( i = 0; i < psnvars; i++)
      mpq_init(psobj[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &pslb, psnvars) );
   for( i = 0; i < psnvars; i++)
      mpq_init(pslb[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &psub, psnvars) );
   for( i = 0; i < psnvars; i++)
      mpq_init(psub[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &pslhs, psnconss) );
   for( i = 0; i < psnconss; i++)
      mpq_init(pslhs[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &psrhs, psnconss) );
   for( i = 0; i < psnconss; i++)
      mpq_init(psrhs[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &psbeg, psnconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pslen, psnconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &psind, psnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &psval, psnnonz) );
   for( i = 0; i < psnnonz; i++)
      mpq_init(psval[i]);
   SCIP_CALL( SCIPallocBufferArray(scip, &dualsol, psnconss) );
   for( i = 0; i < psnconss; i++)
      mpq_init(dualsol[i]);

   mpq_init(objval);

   SCIP_CALL( SCIPallocBufferArray(scip, &colnames, psnvars) );

   /*    printf( "psnvars = %d\n",psnvars); */
   
   for( i = 0; i < psnvars; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames[i],SCIP_MAXSTRLEN ) );
      (void) SCIPsnprintf( colnames[i] , SCIP_MAXSTRLEN, "var%d",i); 
   }


   /* allocate memory for the projection factorization */
   SCIP_CALL( SCIPallocBufferArray(scip, &projbeg,  consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &projlen,  consdata->nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &projind,  consdata->nnonz) ); /* this may be an overestimate, but that is fine */
   for( i = 0; i < consdata->nnonz; i++ )
      projind[i]=0; 
   SCIP_CALL( SCIPallocBufferArray(scip, &projval,  consdata->nnonz) );
   for( i = 0; i < consdata->nnonz; i++ )
      mpq_init(projval[i]);   
   projbasislen = consdata->nconss;
   SCIP_CALL( SCIPallocBufferArray(scip, &projbasis,  projbasislen) );


   /* construct modified LP problem from root node LP */

   /* set objective */
   for( i = 0; i < consdata->nvars + consdata->nconss; i++)
      mpq_set_si(psobj[i ], 0, 1);
   mpq_set_si( psobj[consdata->nvars + consdata->nconss],-1, 1);
   for( i = consdata->nvars + consdata->nconss + 1; i < psnvars; i++)
      mpq_set_si(psobj[i], 1, 1);   
    

   /* we assume all primal variables are unbounded, this is checked here */
    for( i = 0; i < consdata->nvars; i++)
    {
 #ifdef PS_OUT_2 /*??????????*/
       gmp_printf("i=%d> lb=%Qd,ub=%Qd\n", i, consdata->lb[i], consdata->ub[i]);
 #endif
      assert(isNegInfinity(conshdlrdata,consdata->lb[i]));
      assert(isPosInfinity(conshdlrdata,consdata->ub[i]));
    }

   /* set variable bounds */
   for( i = 0; i < psnvars; i++)
      mpq_set(psub[i],conshdlrdata->posinfinity);
   for( i = 0; i < psnvars; i++)
   {
      if( i < consdata->nvars)
         mpq_set( pslb[i], conshdlrdata->neginfinity );
      else
         mpq_set_si( pslb[i], 0, 1);      
   }
   
   /* we assume original constraints are all inequalities (min cx: Ax >= (or <=)  b), this is checked here */
   for( i = 0; i < consdata->nconss; i++)
   {
#ifdef PS_OUT_2 /*??????????*/
      gmp_printf("i=%d> lhs=%Qd, rhs=%Qd\n", i, consdata->lhs[i], consdata->rhs[i]);
#endif
      assert(isNegInfinity(conshdlrdata,consdata->lhs[i]) || isPosInfinity(conshdlrdata,consdata->rhs[i]));
   }
   
   /* set up constraint bounds */
   for( i = 0; i < psnconss; i++)
   {
      if(i >= psnconss - 1 )
      {
         mpq_set_si(pslhs[i],0,1);
         mpq_set_si(psrhs[i],0,1);
      }
      else if( i < consdata->nconss)
      {
         mpq_set_si(pslhs[i],0,1);
         mpq_set(psrhs[i],conshdlrdata->posinfinity);
      }
      else
      {
         mpq_set_si(pslhs[i],1,1);
         mpq_set(psrhs[i],conshdlrdata->posinfinity);
      }    
   }
   
   /* set up constraint matrix */
   for( i = 0; i < consdata->nconss; i++)
   {  
      psbeg[i] = consdata->beg[i] + i;
      pslen[i] = consdata->len[i] + 1;
   }
   for( i = 0; i < psnconss - consdata->nconss; i++)
   {  
      psbeg[i +  consdata->nconss] = consdata->nnonz + consdata->nconss + 2 * i;
      pslen[i +  consdata->nconss] = 2;
   }
   pslen[ psnconss - 1 ] = 0; 	/* this = objective length + 1, and is correctly set later */

   /* first m rows */
   for( i = 0; i < consdata->nconss; i++)
   {  
      if(isNegInfinity(conshdlrdata,consdata->lhs[i])) /* ? */
      {
         for( j = psbeg[i]; j < psbeg[i] + pslen[i] - 1; j++)
         {
            mpq_neg( psval[j], consdata->val[ j - i ] );
            psind[j] = consdata->ind[ j - i ];
         }
      } 
      else
      {
         for( j = psbeg[i]; j < psbeg[i] + pslen[i] - 1; j++)
         {
            mpq_set( psval[j], consdata->val[ j - i ] );
            psind[j] = consdata->ind[ j - i ];
         }
      }
      mpq_set_si( psval[ psbeg[i] + pslen[i] - 1 ], -1,1);
      psind[ psbeg[i] + pslen[i] - 1 ] = consdata->nvars + i;
   }
   /* next m rows */
   for( i = 0; i < consdata->nconss; i++)
   {  
      j = psbeg[i + consdata->nconss];
      mpq_set_si(psval[ j ],1,1);
      psind[ j ] = consdata->nvars + i ;
      mpq_set_si(psval[ j + 1 ],1,1);
      psind[ j + 1 ] = consdata->nvars + consdata->nconss + 1 + i;     
   }   
   /* last row */

   j = psbeg[ 2 * consdata->nconss ];
   /*    printf( "j = %d\n",j); */
   for(i = 0; i < consdata->nvars; i++)
   { 
      if(mpq_sgn(consdata->obj[i]))
      {
         mpq_neg(psval[j], consdata->obj[i]);
         psind[j] = i;
         j++;
      }
   }
   mpq_set_si(psval[j],- 1,1);
   psind[j] = consdata->nconss + consdata->nvars;
   pslen[psnconss - 1] = j - psbeg[ 2 * consdata->nconss ] + 1; 
   assert( j == psnnonz - 1);
   

   /* build modified LP using the exact LP interface */
   SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

   /* add all columns to the exact LP */
   SCIP_CALL( SCIPlpiexAddCols(pslpiex, psnvars, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

   /* add all constraints to the exact LP */
   SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnconss, (const mpq_t*) pslhs, (const mpq_t*) psrhs, NULL, psnnonz, psbeg, pslen, psind, psval) );

   /* write LP to file */
   /*   SCIP_CALL( SCIPlpiexWriteLP(pslpiex, "prob/psdebug.lp") );*/  /* ????????????*/

   /* solve the LP */
   SCIP_CALL( SCIPlpiexSolveDual(pslpiex));

   if( SCIPlpiexIsOptimal(pslpiex) )
   {
      SCIPdebugMessage("   exact LP solved to optimality\n"); 
      /* get optimal dual solution and print it */
      SCIP_CALL( SCIPlpiexGetSol(pslpiex, &objval, NULL, dualsol, NULL, NULL) );
#ifdef PS_OUT /*???????????????????*/
      printf("Dual solution: \n");
      for( i = 0; i < psnconss; i++)
      {
         mpq_out_str(stdout,10,dualsol[i]);
         printf(" \n");
      }   
      printf("Objective value: ");
      mpq_out_str(stdout,10,objval);
      printf(" \n");   
#endif

      /* assign interior point solution to constraint handler data */
      if(mpq_sgn(dualsol[psnconss - 1]))
         mpq_inv(conshdlrdata->commonslack,dualsol[psnconss - 1]);
      else
         mpq_set_si(conshdlrdata->commonslack, 0, 1);

      for( i = 0; i < conshdlrdata->psdim; i++)                                                                                                               
      {                                                                                                                                                       
         mpq_div( conshdlrdata->interiorpt[i], dualsol[i],dualsol[psnconss - 1]);                                                                             
      }     

#ifdef PS_OUT /*???????????????????*/
      printf("Constraints all satisfied by slack of:  ");
      mpq_out_str(stdout,10,conshdlrdata->commonslack);
      printf(" \n"); 
      printf("Relative interior solution: \n");
      for( i = 0; i <  conshdlrdata->psdim; i++)
      {
         mpq_out_str(stdout,10,conshdlrdata->interiorpt[i]);
         printf(" \n");
      } 
#endif


      /* assign implied equality characteristic vector to constraint handler data */
      for( i = 0; i < conshdlrdata->psdim; i++)
      {
         if(mpq_sgn(dualsol[i +conshdlrdata->psdim ]))
            conshdlrdata->impliedeq[i] = 0;
         else
            conshdlrdata->impliedeq[i] = 1;
      }  


#ifdef PS_OUT /*???????????????????*/
      printf("Implied equality: \n");
      for( i = 0; i <  conshdlrdata->psdim; i++)
      {
         printf("%d \n",conshdlrdata->impliedeq[i]);
      } 
#endif
   }
   else if( SCIPlpiexIsObjlimExc(pslpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
   }
   else if( SCIPlpiexIsPrimalInfeasible(pslpiex) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n"); 
   }
   else if( SCIPlpiexExistsPrimalRay(pslpiex) ) 
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
   }
   else if( SCIPlpiexIsIterlimExc(pslpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
   }
   else if( SCIPlpiexIsTimelimExc(pslpiex) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
   }
   else
   {
      SCIPerrorMessage("Other Error\n");
   }


   /* factorize projection matrix */

#ifdef USE_QSFACTOR
   /* select basis variables */                                                                                        
   /* for now, just use first n rows */                                                                                
   for(i = 0; i < projbasislen; i++)                                                                                
      projbasis[i] = i;                                                                                                
                                                                                                                       
   /* construct matrix to pass to solver interface */                                                                  
   projbeg[0] = 0;                                                                                                     
   projlen[0] = consdata->len[projbasis[0]];                                                                           
   for(i = 1; i < consdata->nvars; i++)                                                                                
   {                                                                                                                   
      j =  projbasis[i];                                                                                               
      projlen[i] = consdata->len[j];                                                                                   
      projbeg[i] = projbeg[i-1] + projlen[i-1];                                                                        
   }                                                                                                                   
                                                                                                                       
   for( i = 0; i < consdata->nvars; i++)                                                                               
   {                                                                                                                   
      if(isNegInfinity(conshdlrdata,consdata->lhs[i]))                                                                  
      {                                                                                                                
         for( j = 0 ; j < projlen[ i]; j++)                                                                            
         {                                                                                                             
            mpq_neg( projval[j + projbeg[i]], consdata->val[ consdata->beg[ projbasis[i]] + j  ] );                    
            projind[j + projbeg[i]] = consdata->ind[  consdata->beg[ projbasis[i]] + j  ];                             
         }                                                                                                             
      }                                                                                                                
      else                                                                                                             
      {                                                                                                                
         for( j = 0; j < projlen[i]; j++)                                                                              
         {                                                                                                             
            mpq_set( projval[j + projbeg[i]], consdata->val[ consdata->beg[ projbasis[i]] + j  ] );                    
            projind[j + projbeg[i]] = consdata->ind[  consdata->beg[ projbasis[i]] + j  ];                             
         }                                                                                                             
      }                                                                                                                
   }  
   /* factorize matrix */                                                                   
   SCIPdebugMessage("Factorizing matrix for projection \n"); 
   SCIP_CALL(   SCIPlpiexCreateFactor(conshdlrdata->psfactor,            /**< LP interface structure */
         consdata->nvars,            /**< dimension of matrix */
         projbeg,                    /**< column indices of matrix */
         projlen,                    /**< column lengths of matrix */
         projind,                    /**< row index of entries */
         projval                     /**< coef values of matrix */
         ));
#endif

#ifdef USE_RECTFACTOR

   /* this just uses all columns */
   for(i = 0; i < projbasislen; i++)
      projbasis[i] = i;

   for(i = 0; i < consdata->nconss; i++)
   {
      projlen[i] = consdata->len[i];
      projbeg[i] = consdata->beg[i];
   }
   for( i = 0; i < consdata->nconss; i++)
   { 
      if(isNegInfinity(conshdlrdata,consdata->lhs[i]))
      {
         for( j = 0 ; j < projlen[i]; j++)                                                                                                                
         {
            pos = projbeg[i] + j;
            mpq_neg( projval[pos], consdata->val[pos] );
            projind[pos] = consdata->ind[pos];  
         }
      }      
      else
      {
         for( j = 0 ; j < projlen[i]; j++)                                 
         {            
            pos = projbeg[i] + j;                                                                              
            mpq_set( projval[pos], consdata->val[pos] );        
            projind[pos] = consdata->ind[pos];    
         }
      } 
   }

     
   rval = RECTLUbuildFactorization(
      & conshdlrdata->rectfactor,   /**< pointer to store factor work*/
      consdata->nvars,            /**< number of rows in matrix */
      consdata->nconss,            /**< number of columns in matrix*/
      projbasis,                  /**< basis columns of matrix */
      projval,                    /**< values of matrix entries */
      projind,                    /**< row index of matrix entries */
      projbeg,                    /**< start of columns in sparse matrix */
      projlen                     /**< length of column in sparse matrix */
      );

   if(rval)
      conshdlrdata->psfactorfail= TRUE;
#endif

#ifdef PS_OUT
   if(rval)
      printf("Factorization failed!! \n");
#endif


#ifdef PS_OUT /*???????????????????*/
   printf("Printing matrix used for factorization:\n");
   printf("projbeg: \n");    
   for( i = 0; i < consdata->nconss; i++)
      printf("%d \n",projbeg[i]);                                                                                      
   printf("projlen: \n");                                                                                              
   for( i = 0; i < consdata->nconss; i++)                                                                               
      printf("%d \n",projlen[i]);                                                                                      
   printf("projind: \n");                                                                                              
   for( i = 0; i < consdata->nnonz; i++)                                                                               
      printf("%d \n",projind[i]);                                                                                      
   printf("projval: \n");                                                                                              
   for( i = 0; i < consdata->nnonz; i++)                                                                               
   {                                                                                                                   
      mpq_out_str(stdout,10,projval[i]);                                                                               
      printf(" \n");                                                                                                   
   }                                                                                                                   
#endif  

   conshdlrdata->psdatacon = TRUE; 

   /* free memory */
   SCIPfreeBufferArray(scip, &projbasis);
   for( i = 0; i < consdata->nnonz; i++)
      mpq_clear(projval[i]); 
   SCIPfreeBufferArray(scip, &projval);
   SCIPfreeBufferArray(scip, &projind);
   SCIPfreeBufferArray(scip, &projlen);
   SCIPfreeBufferArray(scip, &projbeg);

   mpq_clear(objval);
   for( i = 0; i < psnvars; i++)
      mpq_clear(psobj[i]);
   for( i = 0; i < psnvars; i++)
      mpq_clear(pslb[i]);
   for( i = 0; i < psnvars; i++)
      mpq_clear(psub[i]);
   for( i = 0; i < psnconss; i++)
      mpq_clear(pslhs[i]);
   for( i = 0; i < psnconss; i++)
      mpq_clear(psrhs[i]);
   for( i = 0; i < psnnonz; i++)
      mpq_clear(psval[i]);
   for( i = 0; i < psnconss; i++)
      mpq_clear(dualsol[i]);
   for( i = psnvars - 1; i >= 0; i--)
      SCIPfreeBufferArray(scip, &colnames[i] );   
   SCIPfreeBufferArray(scip, &colnames);
   SCIPfreeBufferArray(scip, &dualsol);
   SCIPfreeBufferArray(scip, &psval);
   SCIPfreeBufferArray(scip, &psind);
   SCIPfreeBufferArray(scip, &pslen); 
   SCIPfreeBufferArray(scip, &psbeg);
   SCIPfreeBufferArray(scip, &psrhs);
   SCIPfreeBufferArray(scip, &pslhs);
   SCIPfreeBufferArray(scip, &psub);
   SCIPfreeBufferArray(scip, &pslb);
   SCIPfreeBufferArray(scip, &psobj);

   if( pslpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&pslpiex) );
   }
   assert(pslpiex == NULL);

   mpq_clear(mpqtemp);

   return SCIP_OKAY;
}

/** constructs the exact LP for the project and scale dual bound method  */
static
SCIP_RETCODE getPSdualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   mpq_t*                boundval            /**< value of dual bound */
   )
{
   int                   i;
   int                   j;
   int                   rval;               
   int                   currentrow;
   mpq_t*                approxdualsol;      /**< stores the approximate dual solution vector */  
   int                   dim;                /**< column dimension of current dual problem */
   int                   psdim;              /**< column dimension of the root node dual problem */
   mpq_t*                violation;          /**< amount by which the approx dual solution violates constraints */  
   mpq_t*                correction;         /**< stores the correction to the approximate dual solution */
   int                   ndualcons;          /**< number of constraints in dual problem */
   mpq_t                 mpqtemp;
   mpq_t                 lambda1;            /**< multiplier for convex comb. of int. pt and proj sol */
   mpq_t                 lambda2;            /**< multiplier for convex comb. of int. pt and proj sol */ 
   mpq_t                 maxv;               /**< max violation of projected point */   
   mpq_t                 dualbound;          /**< lp bound generated */

   assert(conshdlrdata->psdatacon);
   assert(!conshdlrdata->psfactorfail);

   mpq_init(mpqtemp);
   mpq_init(lambda1);
   mpq_init(lambda2);
   mpq_init(maxv);
   mpq_init(dualbound);

   dim = consdata->nconss;
   psdim = conshdlrdata->psdim;
   ndualcons = consdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &approxdualsol, dim) );
   for( i = 0; i < dim; i++)
      mpq_init(approxdualsol[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &violation, ndualcons) );
   for( i = 0; i < ndualcons; i++)
      mpq_init(violation[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &correction, psdim) );
   for( i = 0; i < psdim; i++)
      mpq_init(correction[i]);

   /* recover the approximate dual solution */
   for( i = 0; i < dim; i++)
   {
      if(isNegInfinity(conshdlrdata,consdata->lhs[i]))
         mpq_set_d(approxdualsol[i],-SCIProwGetDualsol(consdata->rows[i]));
      else
         mpq_set_d(approxdualsol[i],SCIProwGetDualsol(consdata->rows[i]));      
   }
#ifdef PS_OUT 
   printf("Approximate dual solution:\n");
   for( i = 0; i <  dim; i++)
   {
      mpq_out_str(stdout,10,approxdualsol[i]);
      printf(" \n");
   }
   mpq_set_ui(dualbound,0,1);                                                             
   for( i = 0; i < dim; i++)                                                             
   {                                                                                     
      if(isNegInfinity(conshdlrdata,consdata->lhs[i]))                                   
      {                                                                                  
         mpq_mul(mpqtemp,consdata->rhs[i],approxdualsol[i]);                             
         mpq_sub(dualbound,dualbound,mpqtemp);                                                      
      }                                                                                  
      else                                                                               
      {                                                                                  
         mpq_mul(mpqtemp,consdata->lhs[i],approxdualsol[i]);                             
         mpq_add(dualbound,dualbound,mpqtemp);                                                      
      }                                                                                  
   }                                                                                     
   printf("objective value of approx solution is: \n");         
   mpq_out_str(stdout,10,dualbound);                                                      
   printf(" \n");                                                                        
#endif  



   /* first, ensure implied equalities are satisfied and set all negative dual values to zero */ 
   for( i = 0; i < dim; i++) 
   { 
      if( mpq_sgn(approxdualsol[i])<0)
         mpq_set_si(approxdualsol[i],0,1);
   } 
   
   /* calculate violation of equality constraints */
   for( i = 0; i < ndualcons; i++)
   {
      mpq_set(violation[i],consdata->obj[i]);
      /* SUBTRACT Ax to get violation b-Ax */
      /* subtract A(approxdualsol) */
   }
   for( i = 0; i < consdata->nconss; i++)
   {
      if(isPosInfinity(conshdlrdata,consdata->rhs[i]))
      {
         for( j = consdata->beg[i] ; j < consdata->beg[i] + consdata->len[i] ; j++)
         {
            currentrow = consdata->ind[j];
            mpq_mul(mpqtemp,approxdualsol[i],consdata->val[j]);
            mpq_sub(violation[currentrow],violation[currentrow],mpqtemp);
         }
      } 
      else
      {
         for( j = consdata->beg[i] ; j < consdata->beg[i] + consdata->len[i] ; j++)
         {
            currentrow = consdata->ind[j];
            mpq_mul(mpqtemp,approxdualsol[i],consdata->val[j]);
            mpq_add(violation[currentrow],violation[currentrow],mpqtemp);
         }
      }
   }


   /* project solution */
#ifdef PS_OUT 
   printf("violation: \n");
   for( i = 0; i < ndualcons; i++)
   {
      mpq_out_str(stdout,10,violation[i]);
      printf(" \n");
   }     
#endif


#ifdef USE_QSFACTOR
   SCIP_CALL( SCIPlpiexFactorSolve(conshdlrdata->psfactor,ndualcons,correction,violation));
#endif


#ifdef USE_RECTFACTOR
   rval = RECTLUsolveSystem( conshdlrdata->rectfactor,ndualcons ,psdim ,violation ,correction);
#endif


#ifdef PS_OUT 
   printf("correction: \n");
   for( i = 0; i < psdim; i++)
   {
      mpq_out_str(stdout,10,correction[i]);
      printf(" \n");
   }   
#endif

   /* Update approximate dual solution with calculated correction */
   /* THIS NEEDS TO BE FIXED, it is only correct if basis is first n columns!!????? */
   for( i = 0; i < psdim; i++)
   {
      mpq_add(approxdualsol[i],approxdualsol[i],correction[i]);
   }   

#ifdef PS_OUT 
   printf("updated dual solution: \n");
   for( i = 0; i < dim; i++)
   {
      mpq_out_str(stdout,10,approxdualsol[i]);
      printf(" \n");
   } 
#endif

   /* calculate max violation of inequality constraints */

   mpq_set_ui(maxv,0,1);

   for( i = 0; i < dim; i++)
   {
      if(mpq_cmp(maxv,approxdualsol[i]) > 0)
         mpq_set(maxv,approxdualsol[i]);
   }    
 
#ifdef PS_OUT 
   printf("maximum violation is: \n");  
   mpq_out_str(stdout,10,maxv);
   printf(" \n");
#endif

   /* scale solution with interior point to be dual feasible */
   /* compute scale factor */
   /* lambda1 = ( slack of int point)/ (slack of int point + max violation) */
   /* lambda2 = 1 - lambda1 */
   mpq_set(lambda1, conshdlrdata->commonslack);
   mpq_sub(mpqtemp, conshdlrdata->commonslack,maxv);
   mpq_div(lambda1,lambda1,mpqtemp);
   mpq_set_si(lambda2,1,1);
   mpq_sub(lambda2,lambda2,lambda1);


   /* perform shift */
   /* dual feasible solution = lambda1 * projected sol + lambda2* interior sol */
   if( mpq_sgn(lambda2))
   {
      for( i = 0; i < dim; i++)
      {
         mpq_mul(approxdualsol[i],approxdualsol[i],lambda1);
      }  
      /* interior point might have smaller dimension than approx sol at node */
      for( i =0; i < psdim; i++)
      {
         mpq_mul(mpqtemp,conshdlrdata->interiorpt[i],lambda2);
         mpq_add(approxdualsol[i],approxdualsol[i],mpqtemp);
      }
   }

#ifdef PS_OUT 
   printf("projected and shifted dual solution: \n");
   printf("(Should be an exact feasible solution) \n");
   for( i = 0; i < dim; i++)            
   {                                                                                     
      mpq_out_str(stdout,10,approxdualsol[i]);           
      printf(" \n");             
   }                                             
#endif 


   /* determine dual objective value of solution */

   mpq_set_ui(dualbound,0,1);
   for( i = 0; i < dim; i++)
   {
      if(isNegInfinity(conshdlrdata,consdata->lhs[i])) 
      {
         mpq_mul(mpqtemp,consdata->rhs[i],approxdualsol[i]);
         mpq_sub(dualbound,dualbound,mpqtemp);
      }
      else 
      {                              
         mpq_mul(mpqtemp,consdata->lhs[i],approxdualsol[i]);
         mpq_add(dualbound,dualbound,mpqtemp);
      }
   }

#ifdef PS_OUT 
   printf("objective value of dual feasible solution is: \n");  
   mpq_out_str(stdout,10,dualbound);
   printf(" \n");
#endif

   /* assign bound value */
   mpq_set(*boundval,dualbound);

   /* free memory */
   for( i = 0; i < psdim; i++)
      mpq_clear(correction[i]);
   SCIPfreeBufferArray(scip, &correction);

   for( i = 0; i < ndualcons; i++)
      mpq_clear(violation[i]);
   SCIPfreeBufferArray(scip, &violation);

   for( i = 0; i < dim; i++)
      mpq_clear(approxdualsol[i]);
   SCIPfreeBufferArray(scip, &approxdualsol);

   mpq_clear(dualbound);
   mpq_clear(maxv);
   mpq_clear(lambda2);
   mpq_clear(lambda1);
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
   assert(consdata->rows == NULL); /* @todo: otherwise the LP/LPEX rows have to be deleted as well */ 
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

   /* @todo: implement more general case with more than one variable */
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
void chgVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;
   
   /* @todo: changing local bounds is currently only supported at the root 
    * (because keeping track of bound changes which is necessary for switching between different subtrees 
    * is not supported yet). implement this. ??????????? */
   if( SCIPgetDepth(scip) > 0 )
      return;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustLb(conshdlrdata, var, newbound, adjustedbound);
   
   /* change the local bound */
   mpq_set(consdata->lbloc[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);
}

/** changes local exact upper bound of variable; if possible, adjusts bound to integral value */
static
void chgVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< exactlp constraint data */
   SCIP_VAR*             var,                /**< problem variable */
   int                   v,                  /**< index of variable in ind and val array to change the bound for */
   const mpq_t           newbound            /**< new value for bound */
   )
{
   mpq_t adjustedbound;
   
   /* @todo: changing local bounds is currently only supported at the root 
    * (because keeping track of bound changes which is necessary for switching between different subtrees 
    * is not supported yet). implement this. ??????????? */
   if( SCIPgetDepth(scip) > 0 )
      return;

   /* adjust bound for integral variables */
   mpq_init(adjustedbound);
   varAdjustUb(conshdlrdata, var, newbound, adjustedbound);
   
   /* change the local bound */
   mpq_set(consdata->ubloc[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);
}

/** changes global exact lower bound of variable; if possible, adjusts bound to integral value */
static
void chgVarLbGlobal(
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
      chgVarLbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound);
   
   /* change the global bound */
   mpq_set(consdata->lb[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);
}

/** changes global exact upper bound of variable; if possible, adjusts bound to integral value */
static
void chgVarUbGlobal(
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
      chgVarUbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound);
   
   /* change the global bound */
   mpq_set(consdata->ub[consdata->ind[v]], adjustedbound);
   assert(mpq_cmp(consdata->lbloc[consdata->ind[v]], consdata->lb[consdata->ind[v]]) >= 0);
   assert(mpq_cmp(consdata->ubloc[consdata->ind[v]], consdata->ub[consdata->ind[v]]) <= 0);

   mpq_clear(adjustedbound);
}


/** changes exact lower bound of variable; if possible, adjusts bound to integral value */
static
void chgVarLb(
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
      chgVarLbGlobal(scip, conshdlrdata, consdata, var, v, adjustedbound);
   else
      chgVarLbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound);

   mpq_clear(adjustedbound);
}

/** changes exact upper bound of variable; if possible, adjusts bound to integral value */
static
void chgVarUb(
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
      chgVarUbGlobal(scip, conshdlrdata, consdata, var, v, adjustedbound);
   else
      chgVarUbLocal(scip, conshdlrdata, consdata, var, v, adjustedbound);

   mpq_clear(adjustedbound);
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  than the current bound; if possible, adjusts bound to integral value
 */
static
void tightenVarLb(
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

   chgVarLb(scip, conshdlrdata, consdata, var, v, adjustedbound);
   
   if( tightened != NULL )
      *tightened = TRUE;

 TERMINATE:
   mpq_clear(adjustedbound);
   return;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  than the current bound; if possible, adjusts bound to integral value
 */
static
void tightenVarUb(
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

   chgVarUb(scip, conshdlrdata, consdata, var, v, adjustedbound);
   
   if( tightened != NULL )
      *tightened = TRUE;

 TERMINATE:
   mpq_clear(adjustedbound);
   return;
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
   
   /* @todo: changing local bounds is currently only supported at the root 
    * (because keeping track of bound changes which is necessary for switching between different subtrees 
    * is not supported yet). implement this. ??????????? */
   if( SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY; 

   /* @todo: implement general case with more than one variable ?????? */
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
            tightenVarUb(scip, conshdlrdata, consdata, var, v, newub, &infeasible, &tightened);

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
            tightenVarLb(scip, conshdlrdata, consdata, var, v, newlb, &infeasible, &tightened);

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

      /* todo: CONTINUE here ???????*/
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
            tightenVarLb(scip, conshdlrdata, consdata, var, v, newlb, &infeasible, &tightened);

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
            tightenVarUb(scip, conshdlrdata, consdata, var, v, newub, &infeasible, &tightened);

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
   if( sourcedata->objsense == SCIP_OBJSEN_MAXIMIZE )
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

      /* todo: check what value is suitable for maxscale such that there is no overflow ?????????? */
      SCIP_CALL( SCIPmpqCalcIntegralScalar((const mpq_t*) newobj, sourcedata->nvars, SCIPinfinity(scip), intscalar, &success) );
      
      if( success )
      {
         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, sourcedata->nvars) );
      
         /* get transformed problem variables */
         origvars = SCIPgetOrigVars(scip);
         assert(SCIPgetNOrigVars(scip) == sourcedata->nvars);
         SCIP_CALL( SCIPgetTransformedVars(scip, sourcedata->nvars, origvars, vars) );
         
         /* scale exact objective values s.t. all become FP representable and store these values in the transformed problem */
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
               success = FALSE;
         }
      
#ifdef DETAILED_DEBUG /*????????????????*/
         checkOrigVars(scip);
#endif
         
         /* free temporary memory */
         SCIPfreeBufferArray(scip, &vars);
      }
   }

   if( success )
   {
      /* create exactlp constraint data for target constraint */
      SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr, sourcedata->objsense, sourcedata->nvars, 
            newobj, sourcedata->lb, sourcedata->ub, sourcedata->nconss, sourcedata->conssize, sourcedata->lhs, sourcedata->rhs,
            sourcedata->nnonz, sourcedata->beg, sourcedata->len, sourcedata->ind, sourcedata->val, FALSE) );
     
      /* create target constraint */
      SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
            SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
            SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

      /* update objscale of transformed constraint */
      mpq_div(targetdata->objscale, targetdata->objscale, intscalar); 
   }

   /* free temporary memory */
   mpq_clear(intscalar);

   for( i = 0; i < sourcedata->nvars; ++i )
      mpq_clear(newobj[i]);
   SCIPfreeBufferArray(scip, &newobj);
   
   if( !success )
   {
      SCIPerrorMessage("given obj coefficient of var could not be scaled to FP representable number; which is required for working with an FP relaxation\n");
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
   int ncolsex;
   int nrowsex;
#ifdef DETAILED_DEBUG /*????????? */
   SCIP_Real oldlb;
#endif

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   SCIPdebugMessage("separating exactlp constraint <%s> on LP solution (LP solstat=%d)\n", SCIPconsGetName(conss[0]),
      SCIPgetLPSolstat(scip));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN; 

   /* in case the FP problem is a relaxation of the original problem and we use Neumaier and Shcherbinas 
    * dual bounding method, we have already calculated a proved lower bound via postprocessing the LP solution 
    * of the FP problem 
    */
   if( SCIPuseFPRelaxation(scip) && SCIPdualBoundMethod(scip) == 'n' )
      return SCIP_OKAY;

   switch( SCIPdualBoundMethod(scip) )
   {
   case 'e':
      /* dual bound will be calculated in enfops methode, as we can not branch here */
      break;

   case 'v':
      /* constructs exact LP of current node */
      SCIP_CALL( constructCurrentLPEX(scip, conshdlrdata, consdata) );

      /* update lower bound of current node wrt the pseudo objective value */
      if( !SCIPuseFPRelaxation(scip) )
      {
#ifdef DETAILED_DEBUG /*????????? */
         oldlb = SCIPgetLocalLowerbound(scip);
#endif

         SCIP_CALL( SCIPupdateLocalLowerbound(scip, 
               mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) ); /* todo: check whether it is ok to use this function instead of SCIPupdateLocalDualbound() ?????????? */ 

#ifdef DETAILED_DEBUG /*????????? */
         if( oldlb < SCIPgetLocalLowerbound(scip) )
            SCIPdebugMessage("by pseudosol: lower bound improved: %.50f --> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
         else
            SCIPdebugMessage("by pseudosol: lower bound did not improve: %.50f -/-> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
#endif
      }

      /* test whether LP state (i.e. basis information) of inexact LP is dual feasbile in rational arithmetic */
      SCIP_CALL( SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncolsex) );
      SCIP_CALL( SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrowsex) );
      if( ncolsex == SCIPgetNLPCols(scip) && nrowsex == SCIPgetNLPRows(scip) )
      {
         SCIP_LPISTATE* lpistate;
         SCIP_Bool dualfeasible;
         mpq_t dualobjval;
         
         mpq_init(dualobjval);

         SCIP_CALL( SCIPgetLPState(scip, &lpistate) ); 

#ifdef DBVERIFY_TEST /* only for testing the db method ?????????????????? */
         {
            SCIP_Bool lperror;

            printf("\n<<<<<<<<<<<<<<<<<< Load LP basis into LPEX solver and solve LPEX exactly - START >>>>>>>>>>>>>>>>\n");
            SCIP_CALL( loadLPState(scip, conshdlrdata) );
            SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, &lperror) );
            printf("<<<<<<<<<<<<<<<<<< Load LP basis into LPEX solver and solve LPEX exactly - END   >>>>>>>>>>>>>>>>\n\n");
         }
#endif

#ifdef EXLPSOLVER_OUT /* only for debugging ????????? */ 
         printf("<<<<<<<<<<<<<<<<<< Load LP basis into LPEX solver and only check basis   - START >>>>>>>>>>>>>>>>\n");
#endif
         SCIP_CALL( SCIPlpiexStateDualFeasible(conshdlrdata->lpiex, SCIPblkmem(scip), lpistate, &dualfeasible, 
               &dualobjval) );
         
         /* free LP state */
         SCIP_CALL( SCIPfreeLPState(scip, &lpistate) );
         
#ifdef EXLPSOLVER_OUT /* only for debugging ????????? */ 
         printf("<<<<<<<<<<<<<<<<<< Load LP basis into LPEX solver and only check basis   - END   >>>>>>>>>>>>>>>>\n\n");
#endif

         SCIPdebugMessage("DB method <v>: LP basis %s dual feasible\n", dualfeasible ? "is" : "is not");
      
         if( dualfeasible )
         {
#ifdef DETAILED_DEBUG /*????????? */
            oldlb = SCIPgetLocalLowerbound(scip);
#endif

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) ); /* todo: check whether it is ok to use this function instead of SCIPupdateLocalDualbound() ?????????? */ 

#ifdef DETAILED_DEBUG /*????????? */
            if( oldlb < SCIPgetLocalLowerbound(scip) )
               SCIPdebugMessage("by db method (verify): lower bound improved: %.50f --> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
#endif
         }

         mpq_clear(dualobjval);
      }
      break;

   case 'r':
      SCIPerrorMessage("Dual bounding method <%c> has not been implemented yet\n", SCIPdualBoundMethod(scip));
      return SCIP_ERROR;
      //break;

   case 'p':
      {
         mpq_t dualobjval;                                     

         if( conshdlrdata->psfactorfail )
            break;
         SCIPdebugMessage("Computing bound by project and scale\n");
         SCIP_CALL( constructPSData(scip, conshdlrdata, consdata) );

         if( conshdlrdata->psfactorfail )
            break;
         mpq_init(dualobjval);
         SCIP_CALL( getPSdualbound(scip, conshdlrdata, consdata, & dualobjval) );

#ifdef DETAILED_DEBUG /*????????? */          
         oldlb = SCIPgetLocalLowerbound(scip);                                                                     
#endif 

         SCIP_CALL( SCIPupdateLocalLowerbound(scip, mpqGetRealRelax(scip, dualobjval, GMP_RNDD)) ); 

#ifdef DETAILED_DEBUG /*????????? */  
         if( oldlb < SCIPgetLocalLowerbound(scip) )               
            SCIPdebugMessage("by db method (project): lower bound improved: %.50f --> %.50f\n", oldlb, SCIPgetLocalLowerbound(scip));
#endif                                                                                                                
 
      
         mpq_clear(dualobjval);
         /*    SCIPerrorMessage("Dual bounding method <%c> has not been completely implemented yet\n", SCIPdualBoundMethod(scip));*/
         /*return SCIP_ERROR;*/
         break;
      }
   default:
      SCIPerrorMessage("invalid parameter setting <%c> for dual bounding method\n", SCIPdualBoundMethod(scip));
      return SCIP_PARAMETERWRONGVAL;
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

   /* updates lower bound of current node wrt the pseudo objective value */
   if( !SCIPuseFPRelaxation(scip) )
   {
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, 
            mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) ); /* todo: check whether it is ok to use this function instead of SCIPupdateLocalDualbound() ?????????? */ 
   }
   
   /* load LP state from inexact LP into exact LP solver */
   SCIP_CALL( loadLPState(scip, conshdlrdata) );

   /* solve exact LP */
   algo = 'd'; /* todo: what should be done here (selection by warmstartinfo is not useful, algo used in inexact LP solver better?) ?????????? */
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
      return SCIP_ERROR;
   }

   /* evaluate result of exact LP solver */
   SCIP_CALL( evaluateLPEX(scip, conshdlrdata, consdata, result) );

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
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss == 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_INFEASIBLE;

   /* todo: (later) this is ugly and only used to avoid stopping with the hard time limit. 
    *       Should be improved later (set a time and iteration limit on solving LPEX anyway) ?????????????? */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_TIMELIMIT )
      return SCIP_OKAY;

   /* constructs exact LP of current node */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   /* updates lower bound of current node wrt the pseudo objective value */
   if( !SCIPuseFPRelaxation(scip) )
   {
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, 
            mpqGetRealRelax(scip, *getPseudoObjval(scip, conshdlrdata, consdata), GMP_RNDD)) ); /* todo: check whether it is ok to use this function instead of SCIPupdateLocalDualbound() ?????????? */ 
   }

   /* if the inexact LP was solved at the node, load LP state from inexact LP into exact LP solver */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_NOTSOLVED )
   {
      SCIP_CALL( loadLPState(scip, conshdlrdata) );
   }

   /* solve exact LP */
   algo = 'd'; /* todo: what should be done here (selection by warmstartinfo is not useful, algo used in inexact LP solver better?) ?????????? */
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
      return SCIP_ERROR;
   }

   /* evaluate result of exact LP solver */
   SCIP_CALL( evaluateLPEX(scip, conshdlrdata, consdata, result) );

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
         
         /* todo: CONTINUE HERE ???????? 
          *  - search for redundant constraints 
          *    (in particulat delete constraints with only one variable if it is redundant) 
          */
         /* check constraint for infeasibility and redundancy */
         singleConsGetActivityBounds(conshdlrdata, consdata, i, minactivity, maxactivity);
         
         if( mpq_cmp(minactivity, consdata->rhs[i]) > 0 || mpq_cmp(maxactivity, consdata->lhs[i]) < 0 )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d> of exactlp constraint is infeasible: activitybounds=[%Qd,%Qd], sides=[%Qd,%Qd]\n",
               i, minactivity, maxactivity, consdata->lhs[i], consdata->rhs[i]);
            SCIPdebugMessage(s);
#endif
            cutoff = TRUE;
            break;
         }
         else if( mpq_cmp(minactivity, consdata->lhs[i]) >= 0 && mpq_cmp(maxactivity, consdata->rhs[i]) <= 0 )
         {
#ifdef SCIP_DEBUG
            gmp_snprintf(s, SCIP_MAXSTRLEN, "linear constraint <%d> of exactlp constraint is redundant: activitybounds=[%Qd,%Qd], sides=[%Qd,%Qd]\n",
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

      /* catch only bound changes that are not originated by the exactlp constraint handler (e.g., exact presolving),
       * i.e., those that come from branching decisions. the bound changes found by the exactlp constraint handler are 
       * already stored here and the given bound is only a relaxation and probably weaker. 
       * bound changes on integral variables can be catched in any case as they are identical to the ones already stored. 
       */
      /* todo: check whether this argumentation is correct ???????? */
      if( SCIPvarGetType(SCIPeventGetVar(event)) == SCIP_VARTYPE_CONTINUOUS )
         break;

#ifdef BOUNDCHG_OUT /* ?????????? */
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "store unprocessed bound change (%x): <%s> [%Qd,%Qd] >= %g (old: nbndchglb=%d,nbndchgub=%d)\n", 
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lb[varind], consdata->ub[varind], 
            newbound, consdata->nbndchglb, consdata->nbndchgub);
         SCIPdebugMessage(s);
      }
#endif

      /* check whether variable already has unprocessed lower bound change events */ 
      if( consdata->bndchglbpos[varind] == -1 )
      {
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("   first unprocessed boundchange of variable!\n");
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
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("   NOT the first unprocessed boundchange of variable!\n");
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

      /* catch only bound changes that are not originated by the exactlp constraint handler (e.g., exact presolving),
       * i.e., those that come from branching decisions. the bound changes found by the exactlp constraint handler are 
       * already stored here and the given bound is only a relaxation and probaly weaker. 
       * bound changes on integral variables can be catched in any case as they are identical to the ones already stored. 
       */
      /* todo: check whether this argumentation is correct ???????? */
      if( SCIPvarGetType(SCIPeventGetVar(event)) == SCIP_VARTYPE_CONTINUOUS )
         break;

#ifdef BOUNDCHG_OUT /* ?????????? */
      {
         char s[SCIP_MAXSTRLEN];
         gmp_snprintf(s, SCIP_MAXSTRLEN, "store unprocessed bound change (%x): <%s> [%Qd,%Qd] <= %g (old: nbndchglb=%d,nbndchgub=%d)\n", 
            SCIPeventGetType(event), SCIPvarGetName(SCIPeventGetVar(event)), consdata->lb[varind], consdata->ub[varind], 
            newbound, consdata->nbndchglb, consdata->nbndchgub);
         SCIPdebugMessage(s);
      }
#endif

      /* check whether variable already has unprocessed upper bound change events */ 
      if( consdata->bndchgubpos[varind] == -1 )
      {
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("   first unprocessed boundchange of variable!\n");
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
#ifdef BOUNDCHG_OUT /* ?????????? */
         SCIPdebugMessage("   NOT the first unprocessed boundchange of variable!\n");
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

#ifdef BOUNDCHG_OUT /* ?????????? */
   SCIPdebugMessage("   --> new: nbndchglb=%d, nbndchgub=%d\n", consdata->nbndchglb, consdata->nbndchgub);
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
         consPrintExactlp, consCopyExactlp, consParseExactlp,
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
   int                   nconss,             /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each constraint in ind and val array */
   int*                  len,                /**< number of nonzeros in val array corresponding to constraint */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries (and some zeros) */
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
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, objsense, nvars, obj, lb, ub, nconss, nconss, lhs, rhs, 
         nnonz, beg, len, ind, val, objneedscaling) );

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
   
#ifdef BOUNDING_OUT /* ???????????????? */
   printf("SCIPgetExternSafeObjval() for %s=%g\n", lowerbound ? "lowerbound" : "upperbound", objval);
#endif
   if( SCIPisInfinity(scip, objval) )
      mpq_set(mpqobjintern, *posInfinity(conshdlrdata));
   else if( SCIPisInfinity(scip, -objval) )
      mpq_set(mpqobjintern, *negInfinity(conshdlrdata));
   else
      mpq_set_d(mpqobjintern, objval);
   
   getExternObjvalex(cons, mpqobjintern, mpqobjextern); 

   down = (((SCIPgetObjsense(scip) == SCIP_OBJSEN_MINIMIZE) && lowerbound) 
      || ((SCIPgetObjsense(scip) == SCIP_OBJSEN_MAXIMIZE) && !lowerbound));

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

   mpq_init(obj);
   
   SCIPgetSolexOrigObj(scip, cons, sol, obj);
   gmp_snprintf(s, SCIP_MAXSTRLEN, "objective value:                 %20Qd\n", obj);
   SCIPmessageFPrintInfo(file, s);

   /* todo: usinf scip->orig/transprob is ugly! think about a general way get information from the constraint handler 
    *       implement all the methods in cons_exactlp.c that operate and use exact problem data, like exact solutions
    *       ??????????? 
    */
   SCIP_CALL( SCIPsolexPrint(sol, scip->origprob, scip->transprob, file, printzeros) );

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
   SCIPmessageFPrintInfo(file, s);

   /* todo: usinf scip->orig/transprob is ugly! think about a general way get information from the constraint handler 
    *       implement all the methods in cons_exactlp.c that operate and use exact problem data, like exact solutions
    *       ??????????? 
    */
   SCIP_CALL( SCIPsolexPrint(sol, scip->transprob, scip->transprob, file, printzeros) );

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
      SCIPmessageFPrintInfo(file, "no solution available\n");
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
      SCIPmessageFPrintInfo(file, "no exact solution available\n");
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
      SCIPmessageFPrintInfo(file, "no exact solution available\n");
   else
   {
      mpq_init(solval);
      
      SCIPsolexGetVal(sol, var, solval);

      SCIPmessageFPrintInfo(file, "%-32s", SCIPvarGetName(var));
      gmp_snprintf(s, SCIP_MAXSTRLEN, " %20Qd\n", solval);
      SCIPmessageFPrintInfo(file, s);
      
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

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlp\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   mpq_set(ub, consdata->ub[SCIPvarGetProbindex(var)]);
}

