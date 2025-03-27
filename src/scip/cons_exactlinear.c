/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_exactlinear.c
 * @brief  Constraint handler for exact linear constraints in their most general form, \f$lhs <= a^T x <= rhs\f$.
 * @author Leon Eifler
 * @author Sander Borst
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/clock.h"
#include "scip/def.h"
#include "scip/struct_stat.h"
#include "scip/type_retcode.h"
#include "blockmemshell/memory.h"
#include "scip/certificate.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_exactlinear.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"
#include "scip/debug.h"
#include "scip/intervalarith.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_event.h"
#include "scip/pub_lp.h"
#include "scip/pub_lpexact.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_event.h"
#include "scip/scip_exact.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/var.h"
#include "scip/sepastoreexact.h"
#include <ctype.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#define CONSHDLR_NAME          "exactlinear"
#define CONSHDLR_DESC          "exact linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY   +100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "exactlinear"
#define EVENTHDLR_DESC         "bound change event handler for exact linear constraints"

#define DEFAULT_TIGHTENBOUNDSFREQ       1 /**< multiplier on propagation frequency, how often the bounds are tightened */
#define DEFAULT_MAXROUNDS               5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT          -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS            50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT       200 /**< maximal number of cuts separated per separation round in root node */
#define DEFAULT_SORTVARS             TRUE /**< should variables be sorted after presolve w.r.t their coefficient absolute for faster
                                           *  propagation? */
#define DEFAULT_SEPARATEALL         FALSE /**< should all constraints be subject to cardinality cut generation instead of only
                                           *   the ones with non-zero dual value? */
#define DEFAULT_LIMITDENOM          FALSE /**< should denominator sizes for continuous variables be controlled?*/
#define DEFAULT_BOUNDMAXDENOM        256L /**< maximal denominator for rational bounds on continuous variables after propagation */


/** constraint data for linear constraints */
struct SCIP_ConsData
{
   SCIP_RATIONAL*        lhs;                /**< left hand side of row (for ranged rows) */
   SCIP_RATIONAL*        rhs;                /**< right hand side of row */
   SCIP_Real             lhsreal;            /**< real relaxation of lhs */
   SCIP_Real             rhsreal;            /**< real relaxation of rhs */
   SCIP_RATIONAL*        violation;          /**< used to store violation */
   SCIP_RATIONAL*        activity;           /**< used to store activity */
   SCIP_Real             maxabsval;          /**< maximum absolute value of all coefficients */
   SCIP_Real             minabsval;          /**< minimal absolute value of all coefficients */
   SCIP_Real             minactivity;        /**< minimal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             maxactivity;        /**< maximal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastminactivity;    /**< last minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastmaxactivity;    /**< last maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             glbminactivity;     /**< minimal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             glbmaxactivity;     /**< maximal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastglbminactivity; /**< last global minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastglbmaxactivity; /**< last global maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             maxactdelta;        /**< maximal activity contribution of a single variable, or SCIP_INVALID if invalid */
   SCIP_VAR*             maxactdeltavar;     /**< variable with maximal activity contribution, or NULL if invalid */
   SCIP_RATIONAL*        maxabsvalexact;     /**< exact maximum absolute value of all coefficients */
   SCIP_RATIONAL*        minabsvalexact;     /**< exact minimal absolute value of all coefficients */
   SCIP_ROW*             rowlhs;             /**< LP row, if constraint is already stored in LP row format; represents fp-relaxation of lhs-part of rowexact;
                                                  only this row will be added to the exact LP, rowrhs is used for safe aggregation of rows */
   SCIP_ROW*             rowrhs;             /**< LP row, if constraint is already stored in LP row format; represents fp-relaxation of rhs-part of rowexact */
   SCIP_ROWEXACT*        rowexact;           /**< Exact rational lp row */
   SCIP_VAR**            vars;               /**< variables of constraint entries */
   SCIP_RATIONAL**       vals;               /**< coefficients of constraint entries */
   SCIP_INTERVAL*        valsreal;           /**< values of val rounded up/down to closest fp-representable numbers */
   SCIP_EVENTDATA**      eventdata;          /**< event data for bound change events of the variables */
   SCIP_EVENTDATA**      roweventdata;       /**< event data for bound change events of the variables in the row */
   int                   minactivityneginf;  /**< number of coefficients contributing with neg. infinite value to minactivity */
   int                   minactivityposinf;  /**< number of coefficients contributing with pos. infinite value to minactivity */
   int                   maxactivityneginf;  /**< number of coefficients contributing with neg. infinite value to maxactivity */
   int                   maxactivityposinf;  /**< number of coefficients contributing with pos. infinite value to maxactivity */
   int                   minactivityneghuge; /**< number of coefficients contributing with huge neg. value to minactivity */
   int                   minactivityposhuge; /**< number of coefficients contributing with huge pos. value to minactivity */
   int                   maxactivityneghuge; /**< number of coefficients contributing with huge neg. value to maxactivity */
   int                   maxactivityposhuge; /**< number of coefficients contributing with huge pos. value to maxactivity */
   int                   glbminactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbminactivity */
   int                   glbminactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbminactivity */
   int                   glbmaxactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbmaxactivity */
   int                   glbmaxactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbmaxactivity */
   int                   glbminactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbminactivity */
   int                   glbminactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbminactivity */
   int                   glbmaxactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbmaxactivity */
   int                   glbmaxactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbmaxactivity */
   int                   varssize;           /**< size of the vars- and vals-arrays */
   int                   nvars;              /**< number of nonzeros in constraint */
   int                   nbinvars;           /**< the number of binary variables in the constraint, only valid after
                                              *   sorting in stage >= SCIP_STAGE_INITSOLVE
                                              */
   unsigned int          boundstightened:2;  /**< is constraint already propagated with bound tightening? */
   unsigned int          rangedrowpropagated:2; /**< did we perform ranged row propagation on this constraint?
                                                 *   (0: no, 1: yes, 2: with potentially adding artificial constraint */
   unsigned int          validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int          validminabsval:1;   /**< is the minimum absolute value valid? */
   unsigned int          validactivities:1;  /**< are the activity bounds (local and global) valid? */
   unsigned int          validminact:1;      /**< is the local minactivity valid? */
   unsigned int          validmaxact:1;      /**< is the local maxactivity valid? */
   unsigned int          validglbminact:1;   /**< is the global minactivity valid? */
   unsigned int          validglbmaxact:1;   /**< is the global maxactivity valid? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int          coefsorted :1;      /**< are the constraint's variables sorted? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
   unsigned int          implsadded:1;       /**< were the implications of the constraint already extracted? */
   unsigned int          indexsorted:1;      /**< are binary variables sorted w.r.t. the absolute value of their coefficient? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
   unsigned int          hascontvar:1;       /**< does the constraint contain at least one continuous variable? */
   unsigned int          hasnonbinvar:1;     /**< does the constraint contain at least one non-binary variable? */
   unsigned int          hasnonbinvalid:1;   /**< is the information stored in hasnonbinvar and hascontvar valid? */
   unsigned int          onerowrelax:1;      /**< is one floating-point row enough for the fp-relaxation? if so only rowlhs is used */
   unsigned int          hasfprelax:1;       /**< is the constraint possible to be represented as a fp relaxation (only false if var without bound is present) */
};

/** event data for bound change event */
struct SCIP_EventData
{
   SCIP_CONS*            cons;               /**< linear constraint to process the bound change for */
   int                   varpos;             /**< position of variable in vars array */
   bool                  rowvar;             /**< is the event a row event? */
   int                   filterpos;          /**< position of event in variable's event filter */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_EXLINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   SCIP_RATIONAL*        maxaggrnormscale;   /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                              *   (0.0: disable constraint aggregation) */
   SCIP_RATIONAL*        maxcardbounddist;   /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for separating knapsack cardinality cuts */
   SCIP_RATIONAL*        mingainpernmincomp; /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
   SCIP_RATIONAL*        maxeasyactivitydelta;/**< maximum activity delta to run easy propagation on linear constraint
                                               *   (faster, but numerically less stable) */
   int                   linconsupgradessize;/**< size of linconsupgrade array */
   int                   nlinconsupgrades;   /**< number of linear constraint upgrade methods */
   int                   tightenboundsfreq;  /**< multiplier on propagation frequency, how often the bounds are tightened */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
   int                   naddconss;          /**< number of added constraints */
   SCIP_Longint          ncheckserrorbound;  /**< number of times running error analyis activity computation was called */
   SCIP_Longint          nsuccesserrorbound; /**< number of times running error analyis activity computation could determine feasibility */
   SCIP_Longint          nabotserrorbound;    /**< number of times running error analysis activity computation not appliccable (e.g. row->len != fprow->len) */
   SCIP_Longint          nconsprop;         /**< number of times a constraint was propagated */
   SCIP_Longint          nconspropnoninit;  /**< number of times a non-initial (conflict) constraint was propagated */
   SCIP_Longint          propnonzeros;       /**< number of nonzeros in propagated rows */
   SCIP_Longint          propnonzerosnoninit;/**< number of nonzeros in propagated rows in non-initial (conflict) propagations */
   SCIP_Bool             separateall;        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   SCIP_Bool             sortvars;           /**< should binary variables be sorted for faster propagation? */
   SCIP_Bool             propcont;           /**< should bounds on continuous variables be tightened by propagation?*/
   SCIP_Bool             limitdenom;         /**< should denominator sizes for continuous variables be controlled?*/
   SCIP_Longint          boundmaxdenom;      /**< maximal denominator for rational bounds on continuous variables after propagation */
};


/*
 * Propagation rules
 */

/*lint --e{749} */
enum Proprule
{
   PROPRULE_1_RHS        = 1,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the right hand side of the inequality */
   PROPRULE_1_LHS        = 2,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the left hand side of the inequality */
   PROPRULE_1_RANGEDROW  = 3,                /**< fixed variables and gcd of all left variables tighten bounds of a
                                              *   single variable in this reanged row */
   PROPRULE_INVALID      = 0                 /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    proprule:8;         /**< propagation rule that was applied */
         unsigned int    pos:24;             /**< variable position, the propagation rule was applied at */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};

typedef struct InferInfo INFERINFO;


/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}


/** constructs an inference information out of a propagation rule and a position number */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   INFERINFO inferinfo;

   assert(pos >= 0);
   /* in the inferinfo struct only 24 bits for 'pos' are reserved */
   assert(pos < (1<<24));

   inferinfo.val.asbits.proprule = (unsigned int) proprule; /*lint !e641*/
   inferinfo.val.asbits.pos = (unsigned int) pos; /*lint !e732*/

   return inferinfo;
}

/** constructs an inference information out of a propagation rule and a position number, returns info as int */
static
int getInferInt(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   return inferInfoToInt(getInferInfo(proprule, pos));
}

/** ensures, that vars and vals arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int k;
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vals, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->valsreal, consdata->varssize, newsize) );
      for( k = consdata->varssize; k < newsize; ++k )
         SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &consdata->vals[k]) );

      if( consdata->eventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize, newsize) );
      }
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}


/*
 * local methods for managing linear constraint update methods
 */


/** creates constraint handler data for linear constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );
   (*conshdlrdata)->linconsupgrades = NULL;
   (*conshdlrdata)->linconsupgradessize = 0;
   (*conshdlrdata)->nlinconsupgrades = 0;
   (*conshdlrdata)->naddconss = 0;
   (*conshdlrdata)->ncheckserrorbound = 0;
   (*conshdlrdata)->nabotserrorbound = 0;
   (*conshdlrdata)->nsuccesserrorbound = 0;
   (*conshdlrdata)->nconsprop = 0;
   (*conshdlrdata)->nconspropnoninit = 0;
   (*conshdlrdata)->propnonzeros = 0;
   (*conshdlrdata)->propnonzerosnoninit = 0;
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxaggrnormscale) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxcardbounddist) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxeasyactivitydelta) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*conshdlrdata)->mingainpernmincomp) );

   /* set event handler for updating linear constraint activity bounds */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*conshdlrdata)->linconsupgrades, (*conshdlrdata)->linconsupgradessize);

   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxaggrnormscale);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxcardbounddist);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*conshdlrdata)->maxeasyactivitydelta);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*conshdlrdata)->mingainpernmincomp);

   SCIPfreeBlockMemory(scip, conshdlrdata);
}

/*
 * local methods
 */

/** installs rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPrationalIsZero(val));

   if( SCIPrationalIsPositive(val) )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
         !SCIPrationalIsNegInfinity(consdata->lhs), !SCIPrationalIsInfinity(consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
         !SCIPrationalIsInfinity(consdata->rhs), !SCIPrationalIsNegInfinity(consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPrationalIsZero(val));

   if( SCIPrationalIsPositive(val) )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPrationalIsNegInfinity(consdata->lhs),
         !SCIPrationalIsInfinity(consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPrationalIsInfinity(consdata->rhs),
         !SCIPrationalIsNegInfinity(consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
/**! [SnippetDebugAssertions] */
static
SCIP_RETCODE consCatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);
   assert(consdata->vars[pos] != NULL);
   assert(SCIPvarIsTransformed(consdata->vars[pos]));
   assert(consdata->eventdata != NULL);
   assert(consdata->eventdata[pos] == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &(consdata->eventdata[pos])) ); /*lint !e866*/
   consdata->eventdata[pos]->cons = cons;
   consdata->eventdata[pos]->varpos = pos;
   consdata->eventdata[pos]->rowvar = false;

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED,
         eventhdlr, consdata->eventdata[pos], &consdata->eventdata[pos]->filterpos) );

   consdata->removedfixings = consdata->removedfixings && SCIPvarIsActive(consdata->vars[pos]);

   return SCIP_OKAY;
}

/**! [SnippetDebugAssertions] */

/** deletes event data for variable at given position, and drops events */
static
SCIP_RETCODE consDropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);
   assert(consdata->eventdata != NULL);
   assert(consdata->eventdata[pos] != NULL);
   assert(consdata->eventdata[pos]->cons == cons);
   assert(consdata->eventdata[pos]->varpos == pos);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED,
         eventhdlr, consdata->eventdata[pos], consdata->eventdata[pos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->eventdata[pos]); /*lint !e866*/

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consCatchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->eventdata == NULL);

   /* allocate eventdata array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize) );
   assert(consdata->eventdata != NULL);
   BMSclearMemoryArray(consdata->eventdata, consdata->nvars);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consCatchEvent(scip, cons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consDropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->eventdata != NULL);

   /* drop event of every single variable */
   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      SCIP_CALL( consDropEvent(scip, cons, eventhdlr, i) );
   }

   /* free eventdata array */
   SCIPfreeBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize);
   assert(consdata->eventdata == NULL);

   return SCIP_OKAY;
}

/** creates a linear constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of constraint entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of row */
   SCIP_RATIONAL*        rhs                 /**< right hand side of row */
   )
{
   int v;
   SCIP_RATIONAL* constant;
   SCIP_Real lhsrel;
   SCIP_Real rhsrel;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPrationalIsGT(lhs, rhs) )
   {
      SCIPwarningMessage(scip, "left hand side of linear constraint greater than right hand side\n");
      SCIPwarningMessage(scip, " -> lhs=%g, rhs=%g\n", SCIPrationalGetReal(lhs), SCIPrationalGetReal(rhs));
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->varssize = 0;
   (*consdata)->nvars = nvars;
   (*consdata)->hascontvar = FALSE;
   (*consdata)->hasnonbinvar = FALSE;
   (*consdata)->hasnonbinvalid = TRUE;
   (*consdata)->vars = NULL;
   (*consdata)->vals = NULL;
   (*consdata)->valsreal = NULL;

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &constant) );
   if( nvars > 0 )
   {
      int k;

      SCIP_VAR** varsbuffer;
      SCIP_RATIONAL** valsbuffer;
      SCIP_INTERVAL* valsrealbuffer;

      /* copy variables into temporary buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &varsbuffer, nvars) );
      SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &valsbuffer, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &valsrealbuffer, nvars) );
      k = 0;

      /* loop over variables and sort out fixed ones */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_VAR* var;

         var = vars[v];

         assert(var != NULL);
         if( !SCIPrationalIsZero(vals[v]) )
         {
            /* treat fixed variable as a constant if problem compression is enabled */
            if( SCIPisConsCompressionEnabled(scip) && SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetUbGlobalExact(var)) )
            {
               SCIPrationalAddProd(constant, vals[v], SCIPvarGetLbGlobalExact(var));
            }
            else
            {
               varsbuffer[k] = var;
               SCIPrationalSetRational(valsbuffer[k], vals[v]);
               SCIPintervalSetRational(&(valsrealbuffer[k]), vals[v]);
               k++;

               /* update hascontvar and hasnonbinvar flags */
               if( !(*consdata)->hascontvar )
               {
                  SCIP_VARTYPE vartype = SCIPvarGetType(var);

                  if( vartype != SCIP_VARTYPE_BINARY )
                  {
                     (*consdata)->hasnonbinvar = TRUE;

                     if( vartype == SCIP_VARTYPE_CONTINUOUS )
                        (*consdata)->hascontvar = TRUE;
                  }
               }
            }
         }
      }
      (*consdata)->nvars = k;

      if( k > 0 )
      {
         /* copy the possibly reduced buffer arrays into block */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, varsbuffer, k) );
         SCIP_CALL( SCIPrationalCopyBlockArray(SCIPblkmem(scip), &(*consdata)->vals, valsbuffer, k) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->valsreal, valsrealbuffer, k) );
         (*consdata)->varssize = k;
      }

      SCIPrationalFreeBufferArray(SCIPbuffer(scip), &valsbuffer, nvars);
      SCIPfreeBufferArray(scip, &varsbuffer);
      SCIPfreeBufferArray(scip, &valsrealbuffer);
   }

   (*consdata)->eventdata = NULL;

   lhsrel = SCIPrationalRoundReal(lhs, SCIP_R_ROUND_DOWNWARDS);
   rhsrel = SCIPrationalRoundReal(rhs, SCIP_R_ROUND_UPWARDS);

   /* due to compressed copying, we may have fixed variables contributing to the left and right hand side */
   if( !SCIPrationalIsZero(constant) )
   {
      if( !SCIPrationalIsAbsInfinity(lhs) )
         SCIPrationalDiff(lhs, lhs, constant);

      if( !SCIPrationalIsAbsInfinity(rhs) )
         SCIPrationalDiff(rhs, rhs, constant);
   }

   (*consdata)->rowlhs = NULL;
   (*consdata)->rowrhs = NULL;
   (*consdata)->rowexact = NULL;
   SCIP_CALL( SCIPrationalCopyBlock(SCIPblkmem(scip), &(*consdata)->lhs, lhs) );
   SCIP_CALL( SCIPrationalCopyBlock(SCIPblkmem(scip), &(*consdata)->rhs, rhs) );
   (*consdata)->lhsreal = lhsrel;
   (*consdata)->rhsreal = rhsrel;
   SCIP_CALL( SCIPrationalCreateString(SCIPblkmem(scip), &(*consdata)->maxabsvalexact, "inf") );
   SCIP_CALL( SCIPrationalCreateString(SCIPblkmem(scip), &(*consdata)->minabsvalexact, "inf") );
   (*consdata)->maxabsval = SCIP_INVALID;
   (*consdata)->minabsval = SCIP_INVALID;
   (*consdata)->minactivity = SCIP_INVALID;
   (*consdata)->maxactivity = SCIP_INVALID;
   (*consdata)->lastminactivity = SCIP_INVALID;
   (*consdata)->lastmaxactivity = SCIP_INVALID;
   (*consdata)->maxactdelta = SCIP_INVALID;
   (*consdata)->maxactdeltavar = NULL;
   (*consdata)->minactivityneginf = -1;
   (*consdata)->minactivityposinf = -1;
   (*consdata)->maxactivityneginf = -1;
   (*consdata)->maxactivityposinf = -1;
   (*consdata)->minactivityneghuge = -1;
   (*consdata)->minactivityposhuge = -1;
   (*consdata)->maxactivityneghuge = -1;
   (*consdata)->maxactivityposhuge = -1;
   (*consdata)->glbminactivity = SCIP_INVALID;
   (*consdata)->glbmaxactivity = SCIP_INVALID;
   (*consdata)->lastglbminactivity = SCIP_INVALID;
   (*consdata)->lastglbmaxactivity = SCIP_INVALID;
   (*consdata)->glbminactivityneginf = -1;
   (*consdata)->glbminactivityposinf = -1;
   (*consdata)->glbmaxactivityneginf = -1;
   (*consdata)->glbmaxactivityposinf = -1;
   (*consdata)->glbminactivityneghuge = -1;
   (*consdata)->glbminactivityposhuge = -1;
   (*consdata)->glbmaxactivityneghuge = -1;
   (*consdata)->glbmaxactivityposhuge = -1;
   (*consdata)->validmaxabsval = FALSE;
   (*consdata)->validminabsval = FALSE;
   (*consdata)->validactivities = FALSE;
   (*consdata)->validminact = FALSE;
   (*consdata)->validmaxact = FALSE;
   (*consdata)->validglbminact = FALSE;
   (*consdata)->validglbmaxact = FALSE;
   (*consdata)->boundstightened = 0;
   (*consdata)->presolved = FALSE;
   (*consdata)->removedfixings = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->normalized = FALSE;
   (*consdata)->upgradetried = FALSE;
   (*consdata)->upgraded = FALSE;
   (*consdata)->indexsorted = (nvars <= 1);
   (*consdata)->merged = (nvars <= 1);
   (*consdata)->cliquesadded = FALSE;
   (*consdata)->implsadded = FALSE;
   (*consdata)->coefsorted = FALSE;
   (*consdata)->nbinvars = -1;
   (*consdata)->varsdeleted = FALSE;
   (*consdata)->rangedrowpropagated = 0;
   (*consdata)->onerowrelax = FALSE;
   (*consdata)->hasfprelax = FALSE;

   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*consdata)->activity) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(*consdata)->violation) );

   if( SCIPisTransformed(scip) )
   {
      /* get transformed variables */
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   /* capture variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPrationalIsZero((*consdata)->vals[v]));
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
   }

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &constant);

   return SCIP_OKAY;
}

/** frees a linear constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->varssize >= 0);

   /* release the row */
   if( (*consdata)->rowlhs != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowlhs) );
   }
   if( (*consdata)->rowrhs != NULL && !(*consdata)->onerowrelax )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowrhs) );
   }

   /* release variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPrationalIsZero((*consdata)->vals[v]));
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }

   SCIPrationalFreeBlockArray(SCIPblkmem(scip), &(*consdata)->vals, (*consdata)->varssize);

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vals, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->valsreal, (*consdata)->varssize);

   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->lhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->rhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->maxabsvalexact);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->minabsvalexact);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->violation);
   SCIPrationalFreeBlock(SCIPblkmem(scip), &(*consdata)->activity);

   SCIPfreeBlockMemory(scip, consdata);
   return SCIP_OKAY;
}
/** prints linear constraint in CIP format to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPrationalIsNegInfinity(consdata->lhs)
      && !SCIPrationalIsInfinity(consdata->rhs)
      && !SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
   {
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
      SCIPinfoMessage(scip, file, " <= ");
   }

   /* print coefficients and variables */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0");
   else
   {
      /* post linear sum of the linear constraint */
      SCIP_CALL( SCIPwriteVarsLinearsumExact(scip, file, consdata->vars, consdata->vals, consdata->nvars, TRUE) );
   }

   /* print right hand side */
   if( SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
   }
   else if( !SCIPrationalIsInfinity(consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->rhs);
   }
   else if( !SCIPrationalIsNegInfinity(consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
   }
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}

/** prints linear constraint and contained solution values of variables to file stream */
static
SCIP_RETCODE consPrintConsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< solution to print */
   SCIP_Bool             useexactsol,        /**< should the exact sol be used */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  [%s] <%s>: ", SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), SCIPconsGetName(cons));

   /* print left hand side for ranged rows */
   /* print left hand side for ranged rows */
   if( !SCIPrationalIsNegInfinity(consdata->lhs)
      && !SCIPrationalIsInfinity(consdata->rhs)
      && !SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
   {
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
      SCIPinfoMessage(scip, file, " <= ");
   }

   /* print coefficients and variables */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0");
   else
   {
      int v;

      /* post linear sum of the linear constraint */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( consdata->vals != NULL )
         {
            if( SCIPrationalIsEQReal(consdata->vals[v], 1.0) )
            {
               if( v > 0 )
                  SCIPinfoMessage(scip, file, " +");
            }
            else if( SCIPrationalIsEQReal(consdata->vals[v], -1.0) )
               SCIPinfoMessage(scip, file, " -");
            else
            {
               SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->vals[v]);
            }
         }
         else if( consdata->nvars > 0 )
            SCIPinfoMessage(scip, file, " +");

         /* print variable name */
         SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[v], TRUE) );

         if( sol != NULL )
         {
            if( !useexactsol )
               SCIPinfoMessage(scip, file, " (%+.9g)", SCIPgetSolVal(scip, sol, consdata->vars[v]));
            else
            {
               SCIP_RATIONAL* tmp;
               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmp) );
               SCIPgetSolValExact(scip, sol, consdata->vars[v], tmp);
               SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, tmp);
               SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmp);
            }
         }
      }
   }

   /* print right hand side */
   if( SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
   }
   else if( !SCIPrationalIsInfinity(consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->rhs);
   }
   else if( !SCIPrationalIsNegInfinity(consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, consdata->lhs);
   }
   else
      SCIPinfoMessage(scip, file, " [free]");

   SCIPinfoMessage(scip, file, ";\n");

   return SCIP_OKAY;
}

/** invalidates activity bounds, such that they are recalculated in next get */
static
void consdataInvalidateActivities(
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   consdata->validactivities = FALSE;
   consdata->validminact = FALSE;
   consdata->validmaxact = FALSE;
   consdata->validglbminact = FALSE;
   consdata->validglbmaxact = FALSE;
   consdata->validmaxabsval = FALSE;
   consdata->validminabsval = FALSE;
   consdata->hasnonbinvalid = FALSE;
   consdata->minactivity = SCIP_INVALID;
   consdata->maxactivity = SCIP_INVALID;
   consdata->lastminactivity = SCIP_INVALID;
   consdata->lastmaxactivity = SCIP_INVALID;
   consdata->maxabsval = SCIP_INVALID;
   consdata->minabsval = SCIP_INVALID;
   consdata->maxactdelta = SCIP_INVALID;
   SCIPrationalSetInfinity(consdata->maxabsvalexact);
   SCIPrationalSetInfinity(consdata->minabsvalexact);
   consdata->maxactdeltavar = NULL;
   consdata->minactivityneginf = -1;
   consdata->minactivityposinf = -1;
   consdata->maxactivityneginf = -1;
   consdata->maxactivityposinf = -1;
   consdata->minactivityneghuge = -1;
   consdata->minactivityposhuge = -1;
   consdata->maxactivityneghuge = -1;
   consdata->maxactivityposhuge = -1;
   consdata->glbminactivity = SCIP_INVALID;
   consdata->glbmaxactivity = SCIP_INVALID;
   consdata->lastglbminactivity = SCIP_INVALID;
   consdata->lastglbmaxactivity = SCIP_INVALID;
   consdata->glbminactivityneginf = -1;
   consdata->glbminactivityposinf = -1;
   consdata->glbmaxactivityneginf = -1;
   consdata->glbmaxactivityposinf = -1;
   consdata->glbminactivityneghuge = -1;
   consdata->glbminactivityposhuge = -1;
   consdata->glbmaxactivityneghuge = -1;
   consdata->glbmaxactivityposhuge = -1;
}

/** computes the pseudo activity of a constraint */
static
void consdataComputePseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_RATIONAL*        pseudoactivity
   )
{
   int i;
   int pseudoactivityposinf;
   int pseudoactivityneginf;
   SCIP_RATIONAL* bound;
   SCIP_RATIONAL* val;

   SCIPrationalSetInt(pseudoactivity, 0L, 1L);

   pseudoactivityposinf = 0;
   pseudoactivityneginf = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      val = consdata->vals[i];
      bound = (SCIPvarGetBestBoundType(consdata->vars[i]) == SCIP_BOUNDTYPE_LOWER) ? SCIPvarGetLbLocalExact(consdata->vars[i]) : SCIPvarGetUbLocalExact(consdata->vars[i]);

      if( SCIPrationalIsInfinity(bound) )
      {
         if( SCIPrationalIsPositive(val) )
            pseudoactivityposinf++;
         else
            pseudoactivityneginf++;
      }
      else
      {
         if( SCIPrationalIsNegInfinity(bound) )
         {
            if( SCIPrationalIsPositive(val) )
               pseudoactivityneginf++;
            else
               pseudoactivityposinf++;
         }
         else
         {
            SCIPrationalAddProd(pseudoactivity, val, bound);
         }
      }
   }

   if( pseudoactivityneginf > 0 && pseudoactivityposinf > 0 )
      /** @todo introduce a rational equivalent of SCIP_INVALID (maybe an additional flag in SCIP_RATIONAL) */
      return;
   else if( pseudoactivityneginf > 0 )
      SCIPrationalSetNegInfinity(pseudoactivity);
   else if( pseudoactivityposinf > 0 )
      SCIPrationalSetInfinity(pseudoactivity);
}

/** recompute the minactivity of a constraint */
static
void consdataRecomputeMinactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;
   SCIP_ROUNDMODE prevmode;

   consdata->minactivity = 0;

   prevmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->valsreal[i].inf > 0.0 ) ? SCIPvarGetLbLocal(consdata->vars[i]) : SCIPvarGetUbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->valsreal[i].inf * bound) && !SCIPisHugeValue(scip, -consdata->valsreal[i].inf * bound) )
         consdata->minactivity += consdata->valsreal[i].inf * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastminactivity = consdata->minactivity;
   SCIPintervalSetRoundingMode(prevmode);
}

/** recompute the maxactivity of a constraint */
static
void consdataRecomputeMaxactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;
   SCIP_ROUNDMODE prevmode;

   consdata->maxactivity = 0.0;
   prevmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->valsreal[i].sup > 0.0 ) ? SCIPvarGetUbLocal(consdata->vars[i]) : SCIPvarGetLbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->valsreal[i].sup * bound) && !SCIPisHugeValue(scip, -consdata->valsreal[i].sup * bound) )
         consdata->maxactivity += consdata->valsreal[i].sup * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastmaxactivity = consdata->maxactivity;

   SCIPintervalSetRoundingMode(prevmode);
}

/** recompute the global minactivity of a constraint */
static
void consdataRecomputeGlbMinactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;
   SCIP_ROUNDMODE prevmode;

   prevmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   consdata->glbminactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      assert(consdata->vars[i] == consdata->vars[i]);
      bound = (consdata->valsreal[i].inf > 0.0 ) ? SCIPvarGetLbGlobal(consdata->vars[i]) : SCIPvarGetUbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->valsreal[i].inf * bound) && !SCIPisHugeValue(scip, -consdata->valsreal[i].inf * bound) )
         consdata->glbminactivity += consdata->valsreal[i].inf * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbminactivity = consdata->glbminactivity;

   SCIPintervalSetRoundingMode(prevmode);
}

/** recompute the global maxactivity of a constraint */
static
void consdataRecomputeGlbMaxactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;
   SCIP_ROUNDMODE prevmode;

   prevmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   consdata->glbmaxactivity = 0.0;
   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->valsreal[i].sup > 0.0 ) ? SCIPvarGetUbGlobal(consdata->vars[i]) : SCIPvarGetLbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->valsreal[i].sup * bound) && !SCIPisHugeValue(scip, -consdata->valsreal[i].sup * bound) )
         consdata->glbmaxactivity += consdata->valsreal[i].sup * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
   SCIPintervalSetRoundingMode(prevmode);
}

/** calculates minimum absolute value of coefficients */
static
void consdataCalcMinAbsvalEx(
   SCIP*                 scip,
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   assert(consdata != NULL);
   assert(!consdata->validminabsval);

   consdata->validminabsval = TRUE;

   if( consdata->nvars > 0 )
   {
      SCIPrationalAbs(consdata->minabsvalexact, consdata->vals[0]);
      assert(!SCIPrationalIsZero(consdata->vals[0]));
   }
   else
      SCIPrationalSetReal(consdata->minabsvalexact, 0.0);

   for( i = 1; i < consdata->nvars; ++i )
   {
      assert(!SCIPrationalIsZero(consdata->vals[i]));

      if( SCIPrationalIsAbsGT(consdata->minabsvalexact, consdata->vals[i]) )
         SCIPrationalAbs(consdata->minabsvalexact, consdata->vals[i]);
   }
}

/** checks the type of all variables of the constraint and sets hasnonbinvar and hascontvar flags accordingly */
static
void consdataCheckNonbinvar(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int v;

   assert(!consdata->hasnonbinvalid);
   consdata->hasnonbinvar = FALSE;
   consdata->hascontvar = FALSE;

   for( v = consdata->nvars - 1; v >= 0; --v )
   {
      SCIP_VARTYPE vartype = SCIPvarGetType(consdata->vars[v]);

      if( vartype != SCIP_VARTYPE_BINARY )
      {
         consdata->hasnonbinvar = TRUE;

         if( vartype == SCIP_VARTYPE_CONTINUOUS )
         {
            consdata->hascontvar = TRUE;
            break;
         }
      }
   }
   assert(consdata->hascontvar || v < 0);

   consdata->hasnonbinvalid = TRUE;
}

#ifdef CHECKMAXACTDELTA
/* checks that the stored maximal activity delta (if not invalid) is correct */
static
void checkMaxActivityDelta(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   if( consdata->maxactdelta != SCIP_INVALID )
   {
      SCIP_Ratoinal* SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &maxactdelta );
      SCIP_Ratoinal* SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &domain );
      SCIP_Ratoinal* SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &delta );
      SCIP_RATIONAL* lb;
      SCIP_RATIONAL* ub;
      int v;

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         lb = SCIPvarGetLbLocalExact(consdata->vars[v]);
         ub = SCIPvarGetUbLocalExact(consdata->vars[v]);

         if( SCIPrationalIsNegInfinity(lb) || SCIPrationalIsInfinity(ub) )
         {
            SCIPrationalSetInfinity(maxactdelta);
            break;
         }

         SCIPrationalDiff(domain, ub, lb);
         SCIPrationalAbs(delta, consdata->vals[v]);
         SCIPrationalMult(delta, delta, domain);

         if( SCIPrationalisGT(delta,maxactdelta) )
         {
            SCIPrationalSetRational(maxactdelta, delta);
         }
      }
      assert(SCIPrationalIsEQ(maxactdelta, consdata->maxactdelta));

      SCIPrationalFreeBuffer(SCIPbuffer(scip), delta);
      SCIPrationalFreeBuffer(SCIPbuffer(scip), domain);
      SCIPrationalFreeBuffer(SCIPbuffer(scip), maxactdelta);   }
}
#else
#define checkMaxActivityDelta(scip, consdata) /**/
#endif

/** recompute maximal activity contribution for a single variable */
static
void consdataRecomputeMaxActivityDelta(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real delta;
   int v;
   consdata->maxactdelta = 0.0;

   if( !consdata->hasnonbinvalid )
      consdataCheckNonbinvar(consdata);

   /* easy case, the problem consists only of binary variables */
   if( !consdata->hasnonbinvar )
   {
      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         if( SCIPvarGetLbLocal(consdata->vars[v]) < 0.5 && SCIPvarGetUbLocal(consdata->vars[v]) > 0.5 )
         {
            delta = SCIPintervalAbsMax(consdata->valsreal[v]);

            if( delta > consdata->maxactdelta )
            {
               consdata->maxactdelta = delta;
               consdata->maxactdeltavar = consdata->vars[v];
            }
         }
      }
      return;
   }

   for( v = consdata->nvars - 1; v >= 0; --v )
   {
      SCIP_Real domain;
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(consdata->vars[v]);
      ub = SCIPvarGetUbLocal(consdata->vars[v]);

      if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
      {
         consdata->maxactdelta = SCIPinfinity(scip);
         consdata->maxactdeltavar = consdata->vars[v];
         break;
      }

      domain = ub - lb;
      delta = SCIPintervalAbsMax(consdata->valsreal[v]) * domain;

      if( delta > consdata->maxactdelta )
      {
         consdata->maxactdelta = delta;
         consdata->maxactdeltavar = consdata->vars[v];
      }
   }
}

/** updates activities for a change in a bound */
static
void consdataUpdateActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed; can be NULL for global bound changes */
   SCIP_Real             oldbound,           /**< old bound of variable */
   SCIP_Real             newbound,           /**< new bound of variable */
   SCIP_INTERVAL         valrange,           /**< coefficient of constraint entry */
   SCIP_BOUNDTYPE        boundtype,          /**< type of the bound change */
   SCIP_Bool             global              /**< is it a global or a local bound change? */
   )
{
   SCIP_Real* activity;
   SCIP_Real* lastactivity;
   int* activityposinf;
   int* activityneginf;
   int* activityposhuge;
   int* activityneghuge;
   SCIP_Real oldcontribution;
   SCIP_Real newcontribution;
   SCIP_Real delta;
   SCIP_Bool validact;
   SCIP_Bool finitenewbound;
   SCIP_Bool hugevalnewcont;
   SCIP_Real val;
   SCIP_ROUNDMODE prevmode;

   prevmode = SCIPintervalGetRoundingMode();

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(global || (var != NULL));
   assert(consdata->validactivities);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->lastminactivity < SCIP_INVALID);
   assert(consdata->lastmaxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);
   assert(consdata->minactivityneghuge >= 0);
   assert(consdata->minactivityposhuge >= 0);
   assert(consdata->maxactivityneghuge >= 0);
   assert(consdata->maxactivityposhuge >= 0);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->lastglbminactivity < SCIP_INVALID);
   assert(consdata->lastglbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);
   assert(consdata->glbminactivityneghuge >= 0);
   assert(consdata->glbminactivityposhuge >= 0);
   assert(consdata->glbmaxactivityneghuge >= 0);
   assert(consdata->glbmaxactivityposhuge >= 0);

   delta = 0.0;

   /* we are updating global activities */
   if( global )
   {
      /* depending on the boundtype and the coefficient, we choose the activity to be updated:
       * lower bound + pos. coef: update minactivity
       * lower bound + neg. coef: update maxactivity, positive and negative infinity counters have to be switched
       * upper bound + pos. coef: update maxactivity
       * upper bound + neg. coef: update minactivity, positive and negative infinity counters have to be switched
       */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         if( valrange.sup > 0.0 )
         {
            activity = &(consdata->glbminactivity);
            lastactivity = &(consdata->lastglbminactivity);
            activityposinf = &(consdata->glbminactivityposinf);
            activityneginf = &(consdata->glbminactivityneginf);
            activityposhuge = &(consdata->glbminactivityposhuge);
            activityneghuge = &(consdata->glbminactivityneghuge);
            validact = consdata->validglbminact;
            SCIPintervalSetRoundingModeDownwards();
            val = valrange.inf;
         }
         else
         {
            activity = &(consdata->glbmaxactivity);
            lastactivity = &(consdata->lastglbmaxactivity);
            activityposinf = &(consdata->glbmaxactivityneginf);
            activityneginf = &(consdata->glbmaxactivityposinf);
            activityposhuge = &(consdata->glbmaxactivityposhuge);
            activityneghuge = &(consdata->glbmaxactivityneghuge);
            validact = consdata->validglbmaxact;
            SCIPintervalSetRoundingModeUpwards();
            val = valrange.sup;
         }
      }
      else
      {
         if( valrange.sup > 0.0 )
         {
            activity = &(consdata->glbmaxactivity);
            lastactivity = &(consdata->lastglbmaxactivity);
            activityposinf = &(consdata->glbmaxactivityposinf);
            activityneginf = &(consdata->glbmaxactivityneginf);
            activityposhuge = &(consdata->glbmaxactivityposhuge);
            activityneghuge = &(consdata->glbmaxactivityneghuge);
            validact = consdata->validglbmaxact;
            SCIPintervalSetRoundingModeUpwards();
            val = valrange.sup;
         }
         else
         {
            activity = &(consdata->glbminactivity);
            lastactivity = &(consdata->lastglbminactivity);
            activityposinf = &(consdata->glbminactivityneginf);
            activityneginf = &(consdata->glbminactivityposinf);
            activityposhuge = &(consdata->glbminactivityposhuge);
            activityneghuge = &(consdata->glbminactivityneghuge);
            validact = consdata->validglbminact;
            SCIPintervalSetRoundingModeDownwards();
            val = valrange.inf;
         }
      }
   }
   /* we are updating local activities */
   else
   {
      /* depending on the boundtype and the coefficient, we choose the activity to be updated:
       * lower bound + pos. coef: update minactivity
       * lower bound + neg. coef: update maxactivity, positive and negative infinity counters have to be switched
       * upper bound + pos. coef: update maxactivity
       * upper bound + neg. coef: update minactivity, positive and negative infinity counters have to be switched
       */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         if( valrange.sup > 0.0 )
         {
            activity = &(consdata->minactivity);
            lastactivity = &(consdata->lastminactivity);
            activityposinf = &(consdata->minactivityposinf);
            activityneginf = &(consdata->minactivityneginf);
            activityposhuge = &(consdata->minactivityposhuge);
            activityneghuge = &(consdata->minactivityneghuge);
            validact = consdata->validminact;
            SCIPintervalSetRoundingModeDownwards();
            val = valrange.inf;
         }
         else
         {
            activity = &(consdata->maxactivity);
            lastactivity = &(consdata->lastmaxactivity);
            activityposinf = &(consdata->maxactivityneginf);
            activityneginf = &(consdata->maxactivityposinf);
            activityposhuge = &(consdata->maxactivityposhuge);
            activityneghuge = &(consdata->maxactivityneghuge);
            validact = consdata->validmaxact;
            SCIPintervalSetRoundingModeUpwards();
            val = valrange.sup;
         }
      }
      else
      {
         if( valrange.sup > 0.0 )
         {
            activity = &(consdata->maxactivity);
            lastactivity = &(consdata->lastmaxactivity);
            activityposinf = &(consdata->maxactivityposinf);
            activityneginf = &(consdata->maxactivityneginf);
            activityposhuge = &(consdata->maxactivityposhuge);
            activityneghuge = &(consdata->maxactivityneghuge);
            validact = consdata->validmaxact;
            SCIPintervalSetRoundingModeUpwards();
            val = valrange.sup;
         }
         else
         {
            activity = &(consdata->minactivity);
            lastactivity = &(consdata->lastminactivity);
            activityposinf = &(consdata->minactivityneginf);
            activityneginf = &(consdata->minactivityposinf);
            activityposhuge = &(consdata->minactivityposhuge);
            activityneghuge = &(consdata->minactivityneghuge);
            validact = consdata->validminact;
            SCIPintervalSetRoundingModeDownwards();
            val = valrange.inf;
         }
      }
   }

   oldcontribution = SCIPintervalNegateReal(val) * oldbound;
   newcontribution = val * newbound;
   hugevalnewcont = SCIPisHugeValue(scip, REALABS(newcontribution));
   finitenewbound = !SCIPisInfinity(scip, REALABS(newbound));

   if( SCIPisInfinity(scip, REALABS(oldbound)) )
   {
      /* old bound was +infinity */
      if( oldbound > 0.0 )
      {
         assert((*activityposinf) >= 1);

         /* we only have to do something if the new bound is not again +infinity */
         if( finitenewbound || newbound < 0.0 )
         {
            /* decrease the counter for positive infinite contributions */
            (*activityposinf)--;

            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            if( !finitenewbound && newbound < 0.0 )
               (*activityneginf)++;
            else if( hugevalnewcont )
            {
               /* if the contribution of this variable is too large, increase the counter for huge values */
               if( newcontribution > 0.0 )
                  (*activityposhuge)++;
               else
                  (*activityneghuge)++;
            }
            /* "normal case": just add the contribution to the activity */
            else
               delta = newcontribution;
         }
      }
      /* old bound was -infinity */
      else
      {
         assert(oldbound < 0.0);
         assert((*activityneginf) >= 1);

         /* we only have to do something ig the new bound is not again -infinity */
         if( finitenewbound || newbound > 0.0 )
         {
            /* decrease the counter for negative infinite contributions */
            (*activityneginf)--;

            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( !finitenewbound && newbound > 0.0 )
               (*activityposinf)++;
            else if( hugevalnewcont )
            {
               /* if the contribution of this variable is too large, increase the counter for huge values */
               if( newcontribution > 0.0 )
                  (*activityposhuge)++;
               else
                  (*activityneghuge)++;
            }
            /* "normal case": just add the contribution to the activity */
            else
               delta = newcontribution;
         }
      }
   }
   else if( SCIPisHugeValue(scip, REALABS(oldcontribution)) )
   {
      /* old contribution was too large and positive */
      if( -oldcontribution > 0.0 )
      {
         assert((*activityposhuge) >= 1);

         /* decrease the counter for huge positive contributions; it might be increased again later,
          * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
          */
         (*activityposhuge)--;

         if( !finitenewbound )
         {
            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( newbound > 0.0 )
               (*activityposinf)++;
            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            else
               (*activityneginf)++;
         }
         else if( hugevalnewcont )
         {
            /* if the contribution of this variable is too large and positive, increase the corresponding counter */
            if( newcontribution > 0.0 )
               (*activityposhuge)++;
            /* if the contribution of this variable is too large and negative, increase the corresponding counter */
            else
               (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = newcontribution;
      }
      /* old contribution was too large and negative */
      else
      {
         assert(-oldcontribution < 0.0);
         assert((*activityneghuge) >= 1);

         /* decrease the counter for huge negative contributions; it might be increased again later,
          * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
          */
         (*activityneghuge)--;

         if( !finitenewbound )
         {
            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( newbound > 0.0 )
               (*activityposinf)++;
            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            else
               (*activityneginf)++;
         }
         else if( hugevalnewcont )
         {
            /* if the contribution of this variable is too large and positive, increase the corresponding counter */
            if( newcontribution > 0.0 )
               (*activityposhuge)++;
            /* if the contribution of this variable is too large and negative, increase the corresponding counter */
            else
               (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = newcontribution;
      }
   }
   /* old bound was finite and not too large */
   else
   {
      if( !finitenewbound )
      {
         /* if the new bound is +infinity, the old contribution has to be subtracted
          * and the counter for positive infinite contributions has to be increased
          */
         if( newbound > 0.0 )
         {
            (*activityposinf)++;
            delta = oldcontribution;
         }
         /* if the new bound is -infinity, the old contribution has to be subtracted
          * and the counter for negative infinite contributions has to be increased
          */
         else
         {
            assert(newbound < 0.0 );

            (*activityneginf)++;
            delta = oldcontribution;
         }
      }
      /* if the contribution of this variable is too large, increase the counter for huge values */
      else if( hugevalnewcont )
      {
         if( newcontribution > 0.0 )
         {
            (*activityposhuge)++;
            delta = oldcontribution;
         }
         else
         {
            (*activityneghuge)++;
            delta = oldcontribution;
         }
      }
      /* "normal case": just update the activity */
      else
         delta = newcontribution + oldcontribution;
   }

   /* update the activity, if the current value is valid and there was a change in the finite part */
   if( validact && (delta != 0.0) )
   {
      /* if the absolute value of the activity is increased, this is regarded as reliable,
       * otherwise, we check whether we can still trust the updated value
       */
      (*activity) = (*activity) + delta;
      assert(!SCIPisInfinity(scip, -(*activity)) && !SCIPisInfinity(scip, *activity));

      if( REALABS((*lastactivity)) < REALABS(*activity) )
      {
         (*lastactivity) = (*activity);
      }
   }

   SCIPintervalSetRoundingMode(prevmode);
}
/** updates minimum and maximum activity for a change in lower bound */
static
void consdataUpdateActivitiesLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, FALSE);

      assert(!SCIPisInfinity(scip, -consdata->minactivity) && !SCIPisInfinity(scip, consdata->minactivity));
      assert(!SCIPisInfinity(scip, -consdata->maxactivity) && !SCIPisInfinity(scip, consdata->maxactivity));
   }
}

/** updates minimum and maximum activity for a change in upper bound */
static
void consdataUpdateActivitiesUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, FALSE);

      assert(!SCIPisInfinity(scip, -consdata->minactivity) && !SCIPisInfinity(scip, consdata->minactivity));
      assert(!SCIPisInfinity(scip, -consdata->maxactivity) && !SCIPisInfinity(scip, consdata->maxactivity));
   }
}

/** updates minimum and maximum global activity for a change in the global lower bound */
static
void consdataUpdateActivitiesGlbLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, TRUE);

      assert(!SCIPisInfinity(scip, -consdata->glbminactivity) && !SCIPisInfinity(scip, consdata->glbminactivity));
      assert(!SCIPisInfinity(scip, -consdata->glbmaxactivity) && !SCIPisInfinity(scip, consdata->glbmaxactivity));
   }
}

/** updates minimum and maximum global activity for a change in global upper bound */
static
void consdataUpdateActivitiesGlbUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, TRUE);

      assert(!SCIPisInfinity(scip, -consdata->glbminactivity) && !SCIPisInfinity(scip, consdata->glbminactivity));
      assert(!SCIPisInfinity(scip, -consdata->glbmaxactivity) && !SCIPisInfinity(scip, consdata->glbmaxactivity));
   }
}

/** updates minimum and maximum activity and maximum absolute value for coefficient addition */
static
void consdataUpdateAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        valExact,           /**< coefficient of constraint entry */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      assert(consdata->maxabsval < SCIP_INVALID);

      absval = MAX(REALABS(val.inf), REALABS(val.sup)); /*lint !e777 !e666*/
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }

   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      assert(consdata->minabsval < SCIP_INVALID);

      absval = MAX(REALABS(val.inf), REALABS(val.sup)); /*lint !e777  !e666*/
      consdata->minabsval = MIN(consdata->minabsval, absval);
   }

      /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      if( SCIPrationalIsAbsEQ(valExact, consdata->maxabsvalexact) )
      {
         consdata->validmaxabsval = FALSE;
         SCIPrationalSetInfinity(consdata->maxabsvalexact);
      }
   }

   /* invalidate minimum absolute value, if this coefficient was the minimum */
   if( consdata->validminabsval )
   {
      if( SCIPrationalIsAbsEQ(valExact, consdata->minabsvalexact) )
      {
         consdata->validminabsval = FALSE;
         SCIPrationalSetInfinity(consdata->minabsvalexact);
      }
   }

   /* update minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, 0.0, SCIPvarGetLbLocal(var), val);
      consdataUpdateActivitiesUb(scip, consdata, var, 0.0, SCIPvarGetUbLocal(var), val);
      consdataUpdateActivitiesGlbLb(scip, consdata, 0.0, SCIPvarGetLbGlobal(var), val);
      consdataUpdateActivitiesGlbUb(scip, consdata, 0.0, SCIPvarGetUbGlobal(var), val);
   }
}

/** updates minimum and maximum activity for coefficient deletion, invalidates maximum absolute value if necessary */
static
void consdataUpdateDelCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        valExact,           /**< exact coefficient of constraint entry */
   SCIP_INTERVAL         val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = SCIPintervalAbsMax(val);

      if( SCIPisEQ(scip, absval, consdata->maxabsval) )
      {
         consdata->validmaxabsval = FALSE;
         consdata->maxabsval = SCIP_INVALID;
      }
   }

   /* invalidate minimum absolute value, if this coefficient was the minimum */
   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      absval = SCIPintervalAbsMax(val);

      if( SCIPisEQ(scip, absval, consdata->minabsval) )
      {
         consdata->validminabsval = FALSE;
         consdata->minabsval = SCIP_INVALID;
      }
   }

   /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      if( SCIPrationalIsAbsEQ(valExact, consdata->maxabsvalexact) )
      {
         consdata->validmaxabsval = FALSE;
         SCIPrationalSetInfinity(consdata->maxabsvalexact);
      }
   }

   /* invalidate minimum absolute value, if this coefficient was the minimum */
   if( consdata->validminabsval )
   {
      if( SCIPrationalIsAbsEQ(valExact, consdata->minabsvalexact) )
      {
         consdata->validminabsval = FALSE;
         SCIPrationalSetInfinity(consdata->minabsvalexact);
      }
   }

   /* update minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, SCIPvarGetLbLocal(var), 0.0, val);
      consdataUpdateActivitiesUb(scip, consdata, var, SCIPvarGetUbLocal(var), 0.0, val);
      consdataUpdateActivitiesGlbLb(scip, consdata, SCIPvarGetLbGlobal(var), 0.0, val);
      consdataUpdateActivitiesGlbUb(scip, consdata, SCIPvarGetUbGlobal(var), 0.0, val);
   }
}

/** returns the minimum absolute value of all coefficients in the constraint */
static
SCIP_RATIONAL* consdataGetMinAbsvalEx(
   SCIP*                 scip,
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validminabsval )
      consdataCalcMinAbsvalEx(scip, consdata);
   assert(consdata->validminabsval);

   return consdata->minabsvalexact;
}


/** updates minimum and maximum activity for coefficient change, invalidates maximum absolute value if necessary */
static
void consdataUpdateChgCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_INTERVAL         oldval,             /**< old coefficient of constraint entry */
   SCIP_RATIONAL*        oldvalExact,        /**< old exact coefficient of constraint entry */
   SCIP_INTERVAL         newval,             /**< new coefficient of constraint entry */
   SCIP_RATIONAL*        newvalExact         /**< new coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = SCIPintervalAbsMax(newval);

      if( absval >= consdata->maxabsval )
      {
         consdata->maxabsval = absval;
      }
      else
      {
         absval = SCIPintervalAbsMax(oldval);

         /* invalidate maximum absolute value */
         if( SCIPisEQ(scip, absval, consdata->maxabsval) )
         {
            consdata->validmaxabsval = FALSE;
            consdata->maxabsval = SCIP_INVALID;
         }
      }
   }

   /* update minimum absolute value */
   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      absval = SCIPintervalAbsMax(newval);

      if( absval <= consdata->minabsval )
      {
         consdata->minabsval = absval;
      }
      else
      {
         absval = SCIPintervalAbsMax(oldval);

         /* invalidate minimum absolute value */
         if( SCIPisEQ(scip, absval, consdata->minabsval) )
         {
            consdata->validminabsval = FALSE;
            consdata->minabsval = SCIP_INVALID;
         }
      }
   }
   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      if( SCIPrationalIsAbsGT(newvalExact, consdata->maxabsvalexact) )
      {
         SCIPrationalAbs(consdata->maxabsvalexact, newvalExact);
      }
      else
      {
         /* invalidate maximum absolute value */
         if( SCIPrationalIsAbsEQ(oldvalExact, consdata->maxabsvalexact) )
         {
            consdata->validmaxabsval = FALSE;
            SCIPrationalSetInfinity(consdata->maxabsvalexact);
         }
      }
   }
      /* update minimum absolute value */
   if( consdata->validminabsval )
   {
      if( SCIPrationalIsAbsGT(consdata->minabsvalexact, newvalExact) )
      {
         SCIPrationalAbs(consdata->minabsvalexact, newvalExact);
      }
      else
      {
         /* invalidate minimum absolute value */
         if( SCIPrationalIsAbsEQ(oldvalExact, consdata->minabsvalexact) )
         {
            consdata->validminabsval = FALSE;
            SCIPrationalSetInfinity(consdata->minabsvalexact);
         }
      }
   }

   /* update maximum activity delta */
   if( !SCIPisInfinity(scip, consdata->maxactdelta ) )
   {
      SCIP_Real domain;
      SCIP_Real delta;

      assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
      assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));

      domain = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
      delta = SCIPintervalAbsMax(newval) * domain;

      if( delta > consdata->maxactdelta )
      {
         consdata->maxactdelta = delta;
         consdata->maxactdeltavar = var;
      }
      else
      {
         /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
         if( consdata->maxactdeltavar == var )
            consdata->maxactdelta = SCIP_INVALID;
      }
   }

   /* @todo as in cons_linear, do something more clever here, e.g. if oldval * newval >= 0, do the update directly */
   consdataUpdateDelCoef(scip, consdata, var, oldvalExact, oldval);
   consdataUpdateAddCoef(scip, consdata, var, newvalExact, newval);
}

/** ensures that every nonzero is a least minval so that we don't get problem with SCIPs 0 in floating point representation */
static
void consdataScaleMinValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             minval              /**< minmimal value for coefficients in constraint */
   )
{
   int i;
   SCIP_RATIONAL* scalingfactor;
   SCIP_RATIONAL* minabsval;

   assert(scip != NULL);
   assert(consdata != NULL);

   minabsval = consdataGetMinAbsvalEx(scip, consdata);

   assert(!SCIPrationalIsZero(minabsval) || consdata->nvars == 0);

   (void) SCIPrationalCreateBuffer(SCIPbuffer(scip), &scalingfactor);

   if( SCIPrationalIsLTReal(minabsval, minval) )
   {
      SCIPrationalSetReal(scalingfactor, minval);
      SCIPrationalDiv(scalingfactor, scalingfactor, minabsval);

      for( i = 0; i < consdata->nvars; i++ )
      {
         SCIPrationalMult(consdata->vals[i], consdata->vals[i], scalingfactor);
         SCIPintervalSetRational(&(consdata->valsreal[i]), consdata->vals[i]);
      }

      SCIPrationalMult(consdata->rhs, consdata->rhs, scalingfactor);
      consdata->rhsreal = SCIPrationalRoundReal(consdata->rhs, SCIP_R_ROUND_UPWARDS);

      SCIPrationalMult(consdata->lhs, consdata->lhs, scalingfactor);
      consdata->lhsreal = SCIPrationalRoundReal(consdata->lhs, SCIP_R_ROUND_DOWNWARDS);
   }

   consdataInvalidateActivities(consdata);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &scalingfactor);
}

/** calculates minimum and maximum local and global activity for constraint from scratch;
 *  additionally recalculates maximum absolute value of coefficients
 */
static
void consdataCalcActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(!consdata->validactivities);
   assert(consdata->minactivity >= SCIP_INVALID || consdata->validminact);
   assert(consdata->maxactivity >= SCIP_INVALID || consdata->validmaxact);
   assert(consdata->glbminactivity >= SCIP_INVALID || consdata->validglbminact);
   assert(consdata->glbmaxactivity >= SCIP_INVALID || consdata->validglbmaxact);

   consdata->validmaxabsval = TRUE;
   consdata->validminabsval = TRUE;
   consdata->validactivities = TRUE;
   consdata->validminact = TRUE;
   consdata->validmaxact = TRUE;
   consdata->validglbminact = TRUE;
   consdata->validglbmaxact = TRUE;
   consdata->maxabsval = 0.0;
   consdata->minabsval = (consdata->nvars == 0 ? 0.0 : SCIPintervalAbsMax(consdata->valsreal[0]));
   consdata->minactivity = 0.0;
   consdata->maxactivity = 0.0;
   consdata->lastminactivity = 0.0;
   consdata->lastmaxactivity = 0.0;
   consdata->minactivityneginf = 0;
   consdata->minactivityposinf = 0;
   consdata->maxactivityneginf = 0;
   consdata->maxactivityposinf = 0;
   consdata->minactivityneghuge = 0;
   consdata->minactivityposhuge = 0;
   consdata->maxactivityneghuge = 0;
   consdata->maxactivityposhuge = 0;
   consdata->glbminactivity = 0.0;
   consdata->glbmaxactivity = 0.0;
   consdata->lastglbminactivity = 0.0;
   consdata->lastglbmaxactivity = 0.0;
   consdata->glbminactivityneginf = 0;
   consdata->glbminactivityposinf = 0;
   consdata->glbmaxactivityneginf = 0;
   consdata->glbmaxactivityposinf = 0;
   consdata->glbminactivityneghuge = 0;
   consdata->glbminactivityposhuge = 0;
   consdata->glbmaxactivityneghuge = 0;
   consdata->glbmaxactivityposhuge = 0;

   for( i = 0; i < consdata->nvars; ++i )
      consdataUpdateAddCoef(scip, consdata, consdata->vars[i], consdata->vals[i], consdata->valsreal[i]);
   consdata->lastminactivity = consdata->minactivity;
   consdata->lastmaxactivity = consdata->maxactivity;
   consdata->lastglbminactivity = consdata->glbminactivity;
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
}

/** computes the activity of a row for a given solution plus a bound on the floating-point error using running error analysis */
static
SCIP_Bool consdataComputeSolActivityWithErrorbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity,           /**< buffer to return floating-point activity */
   SCIP_Real*            errorbound          /**< buffer to return bound on absolute floating-point error */
   )
{
   SCIP_Real solval;
   SCIP_Real sum;
   SCIP_Real mu;
   SCIP_Real inf;
   SCIP_Bool success;
   int v;

   assert(activity != NULL);
   assert(errorbound != NULL);

   inf = SCIPinfinity(scip);
   *activity = SCIP_UNKNOWN;
   *errorbound = inf;

   sum = 0.0;
   mu = 0.0;
   /* normally we want to use the row since all fixed/aggregated variables do not appear there */
   if( consdata->rowlhs == NULL )
   {
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_LOOSE )
            solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);
         else
            return FALSE;

         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            return FALSE;

         sum += consdata->valsreal[v].inf * solval;
         mu += REALABS(sum);
         /* the factor 3 + eps is needed to account for rounding errors in valsreal[v]/solval */
         mu += (3.0 + SCIP_REAL_UNITROUNDOFF) * REALABS(consdata->valsreal[v].inf * solval);
      }
   }
   else
   {
      success = SCIPgetRowSolActivityWithErrorboundExact(scip, consdata->rowexact, sol, &sum, &mu);

      if( !success )
         return FALSE;
   }

   sum = MAX(sum, -inf);
   sum = MIN(sum, +inf);
   *activity = sum;

   if( SCIPisInfinity(scip, sum) || SCIPisInfinity(scip, -sum) )
      *errorbound = inf;
   else
      *errorbound = mu * 1.1 * SCIP_REAL_UNITROUNDOFF;

   return TRUE;
}

/** gets minimal activity for constraint and given values of counters for infinite and huge contributions
 *  and (if needed) delta to subtract from stored finite part of activity in case of a residual activity
 */
static
void getMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   posinf,             /**< number of coefficients contributing pos. infinite value */
   int                   neginf,             /**< number of coefficients contributing neg. infinite value */
   int                   poshuge,            /**< number of coefficients contributing huge pos. value */
   int                   neghuge,            /**< number of coefficients contributing huge neg. value */
   SCIP_Real             delta,              /**< value to subtract from stored minactivity
                                              *   (contribution of the variable set to zero when getting residual activity) */
   SCIP_Bool             global,             /**< should the global or local minimal activity be returned? */
   SCIP_Bool             goodrelax,          /**< should a good relaxation be computed or are relaxed acticities ignored, anyway? */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Bool*            isrelax,            /**< pointer to store whether the activity is a relaxation,
                                              *   i.e. is <= the exact minactivity (in case of huge contributing values) */
   SCIP_Bool*            issettoinfinity     /**< pointer to store whether minactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(posinf >= 0);
   assert(neginf >= 0);
   assert(poshuge >= 0);
   assert(neghuge >= 0);
   assert(minactivity != NULL);
   assert(isrelax != NULL);
   assert(issettoinfinity != NULL);

   /* if we have pos. infinite contributions, the minactivity is +infty */
   if( posinf > 0 )
   {
      *minactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have neg. (and no pos.) infinite contributions, the minactivity is -infty */
   else if( neginf > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have neg. huge contributions, we only know that -infty is a relaxation of the minactivity */
   else if( neghuge > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   /* we do not need a good relaxation and we have positive huge contributions, so we just return -infty as activity */
   else if( !goodrelax && poshuge > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   else
   {
      SCIP_Real tmpactivity;

      /* recompute minactivity if it is not valid */
      if( global )
      {
         if( !consdata->validglbminact )
            consdataRecomputeGlbMinactivity(scip, consdata);
         assert(consdata->validglbminact);

         tmpactivity = consdata->glbminactivity;
      }
      else
      {
         if( !consdata->validminact )
            consdataRecomputeMinactivity(scip, consdata);
         assert(consdata->validminact);

         tmpactivity = consdata->minactivity;
      }

      /* we have no infinite and no neg. huge contributions, but pos. huge contributions;
       * a feasible relaxation of the minactivity is the number of positive huge contributions
       * times the minimum value counting as "huge" plus finite (and non-huge) part of minactivity - delta
       */
      if( poshuge > 0 )
      {
         *minactivity = 1.0 * poshuge * SCIPgetHugeValue(scip) + (tmpactivity - delta);
         *issettoinfinity = FALSE;
         *isrelax = TRUE;
      }
      /* all counters are zero, so the minactivity is just stored and we subtract the delta */
      else
      {
         *minactivity = tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = FALSE;
      }
   }
}

/** gets maximal activity for constraint and given values of counters for infinite and huge contributions
 *  and (if needed) delta to subtract from stored finite part of activity in case of a residual activity
 */
static
void getMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   posinf,             /**< number of coefficients contributing pos. infinite value */
   int                   neginf,             /**< number of coefficients contributing neg. infinite value */
   int                   poshuge,            /**< number of coefficients contributing huge pos. value */
   int                   neghuge,            /**< number of coefficients contributing huge neg. value */
   SCIP_Real             delta,              /**< value to subtract from stored maxactivity
                                              *   (contribution of the variable set to zero when getting residual activity) */
   SCIP_Bool             global,             /**< should the global or local maximal activity be returned? */
   SCIP_Bool             goodrelax,          /**< should a good relaxation be computed or are relaxed acticities ignored, anyway? */
   SCIP_Real*            maxactivity,        /**< pointer to store the maximal activity */
   SCIP_Bool*            isrelax,            /**< pointer to store whether the activity is a relaxation,
                                              *   i.e. is >= the exact maxactivity (in case of huge contributing values) */
   SCIP_Bool*            issettoinfinity     /**< pointer to store whether maxactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(posinf >= 0);
   assert(neginf >= 0);
   assert(poshuge >= 0);
   assert(neghuge >= 0);
   assert(maxactivity != NULL);
   assert(isrelax != NULL);
   assert(issettoinfinity != NULL);

   /* if we have neg. infinite contributions, the maxactivity is -infty */
   if( neginf > 0 )
   {
      *maxactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have pos. (and no neg.) infinite contributions, the maxactivity is +infty */
   else if( posinf > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have pos. huge contributions, we only know that +infty is a relaxation of the maxactivity */
   else if( poshuge > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   /* we do not need a good relaxation and we have positve huge contributions, so we just return +infty as activity */
   else if( !goodrelax && neghuge > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   else
   {
      SCIP_Real tmpactivity;

      /* recompute maxactivity if it is not valid */
      if( global )
      {
         if( !consdata->validglbmaxact )
            consdataRecomputeGlbMaxactivity(scip, consdata);
         assert(consdata->validglbmaxact);

         tmpactivity = consdata->glbmaxactivity;
      }
      else
      {
         if( !consdata->validmaxact )
            consdataRecomputeMaxactivity(scip, consdata);
         assert(consdata->validmaxact);

         tmpactivity = consdata->maxactivity;
      }

      /* we have no infinite, and no pos. huge contributions, but neg. huge contributions;
       * a feasible relaxation of the maxactivity is minus the number of negative huge contributions
       * times the minimum value counting as "huge" plus the finite (and non-huge) part of maxactivity minus delta
       */
      if( neghuge > 0 )
      {
         *maxactivity = -1.0 * neghuge * SCIPgetHugeValue(scip) + tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = TRUE;
      }
      /* all counters are zero, so the maxactivity is just stored and we subtract the delta */
      else
      {
         *maxactivity = tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = FALSE;
      }
   }
}

/** gets activity bounds for constraint */
static
void consdataGetActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed activities ignored, anyway? */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Real*            maxactivity,        /**< pointer to store the maximal activity */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned minactivity is just a relaxation,
                                              *   i.e. <= the exact minactivity (in case of huge contributions),
                                              *   or equal to the exact minimal activity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned maxactivity is just a relaxation,
                                              *   i.e. >= the exact maxactivity (in case of huge contributions),
                                              *   or equal to the exact maximal activity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxactivity was set to infinity or calculated */

   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validminact);
      assert(consdata->validmaxact);
   }
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);

   getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
      consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
      minactivity, minisrelax, isminsettoinfinity);

   getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
      consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
      maxactivity, maxisrelax, ismaxsettoinfinity);
}

/** gets activity bounds for constraint after setting variable to zero */
static
void consdataGetActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_INTERVAL         val,                /**< coefficient value of variable in linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   SCIP_Real*            maxresactivity,     /**< pointer to store the maximal residual activity */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned residual minactivity is just a
                                              *   relaxation, i.e. <= the exact residual minactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned residual maxactivity is just a
                                              *   relaxation, i.e. <= the exact residual maxactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minresactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxresactivity was set to infinity or calculated */
   )
{
   SCIP_Real minactbound;
   SCIP_Real maxactbound;
   SCIP_Real absval;
   SCIP_ROUNDMODE prevmode;
   prevmode = SCIPintervalGetRoundingMode();

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(minisrelax != NULL);
   assert(maxisrelax != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validminact);
      assert(consdata->validmaxact);
   }
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);
   assert(consdata->minactivityneghuge >= 0);
   assert(consdata->minactivityposhuge >= 0);
   assert(consdata->maxactivityneghuge >= 0);
   assert(consdata->maxactivityposhuge >= 0);

   if( val.sup > 0.0 )
   {
      minactbound = SCIPvarGetLbLocal(var);
      maxactbound = SCIPvarGetUbLocal(var);
      absval = val.sup;
   }
   else
   {
      minactbound = -SCIPvarGetUbLocal(var);
      maxactbound = -SCIPvarGetLbLocal(var);
      absval = -val.inf;
   }

   /* get/compute minactivity by calling getMinActivity() with updated counters for infinite and huge values
    * and contribution of variable set to zero that has to be subtracted from finite part of activity
    */
   if( SCIPisInfinity(scip, minactbound) )
   {
      assert(consdata->minactivityposinf >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf - 1, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisInfinity(scip, -minactbound) )
   {
      assert(consdata->minactivityneginf >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf - 1,
         consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, minactbound * absval) )
   {
      assert(consdata->minactivityposhuge >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge - 1, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, -minactbound * absval) )
   {
      assert(consdata->minactivityneghuge >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge - 1, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else
   {
      SCIP_Real delta;
      delta = absval * minactbound;
      SCIPintervalSetRoundingModeUpwards();
      SCIPintervalSetRoundingModeDownwards();
      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge, delta, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }

   /* get/compute maxactivity by calling getMaxActivity() with updated counters for infinite and huge values
    * and contribution of variable set to zero that has to be subtracted from finite part of activity
    */
   if( SCIPisInfinity(scip, -maxactbound) )
   {
      assert(consdata->maxactivityneginf >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf - 1,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisInfinity(scip, maxactbound) )
   {
      assert(consdata->maxactivityposinf >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf - 1, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, absval * maxactbound) )
   {
      assert(consdata->maxactivityposhuge >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge - 1, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, -absval * maxactbound) )
   {
      assert(consdata->maxactivityneghuge >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge - 1, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else
   {
      SCIP_Real delta;
      delta = absval * maxactbound;
      SCIPintervalSetRoundingModeDownwards();
      SCIPintervalSetRoundingModeUpwards();
      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, delta, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   SCIPintervalSetRoundingMode(prevmode);
}

/** calculates the activity of the linear constraint for given solution */
static
void consdataGetActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol,                /**< solution to get activity for, NULL to current solution */
   SCIP_Bool             useexact,           /**< should the exact solution be used */
   SCIP_RATIONAL*        activity            /**< pointer to store the activity */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( (sol == NULL) && !SCIPhasCurrentNodeLP(scip) )
      consdataComputePseudoActivity(scip, consdata, activity);
   else
   {
      SCIP_RATIONAL* solval;
      int nposinf;
      int nneginf;
      SCIP_Bool negsign;
      int v;

      (void) SCIPrationalCreateBuffer(SCIPbuffer(scip), &solval);

      SCIPrationalSetInt(activity, 0L, 1L);
      nposinf = 0;
      nneginf = 0;

      for( v = 0; v < consdata->nvars; ++v )
      {
         if( useexact )
            SCIPgetSolValExact(scip, sol, consdata->vars[v], solval);
         else
            SCIPrationalSetReal(solval, SCIPgetSolVal(scip, sol, consdata->vars[v]));

         if( SCIPrationalIsNegative(consdata->vals[v]) )
            negsign = TRUE;
         else
            negsign = FALSE;

         if( (SCIPrationalIsInfinity(solval) && !negsign) || (SCIPrationalIsNegInfinity(solval) && negsign) )
            ++nposinf;
         else if( (SCIPrationalIsInfinity(solval) && negsign) || (SCIPrationalIsNegInfinity(solval) && !negsign) )
            ++nneginf;
         else
         {
            SCIPrationalAddProd(activity, solval, consdata->vals[v]);
         }
      }
      assert(nneginf >= 0 && nposinf >= 0);

      SCIPdebugMsg(scip, "activity of linear constraint: %.15g, %d positive infinity values, %d negative infinity values \n", SCIPrationalGetReal(activity), nposinf, nneginf);

      /* check for amount of infinity values and correct the activity */
      if( nposinf > 0 && nneginf > 0 )
         /** @todo introduce a rational equivalent of SCIP_INVALID (maybe an additional flag in SCIP_RATIONAL) */
         return;
      else if( nposinf > 0 )
         SCIPrationalSetInfinity(activity);
      else if( nneginf > 0 )
         SCIPrationalSetNegInfinity(activity);

      SCIPrationalDebugMessage("corrected activity of linear constraint: %q\n", activity);
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &solval);
   }
}

/** calculates the feasibility of the linear constraint for given solution */
static
void consdataGetFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol,                /**< solution to get feasibility for, NULL to current solution */
   SCIP_RATIONAL*        ret                 /**< pointer to store the result */
   )
{
   SCIP_RATIONAL* activity;
   SCIP_RATIONAL* op1;
   SCIP_RATIONAL* op2;

   assert(scip != NULL);
   assert(consdata != NULL);

   (void) SCIPrationalCreateBuffer(SCIPbuffer(scip), &activity);
   (void) SCIPrationalCreateBuffer(SCIPbuffer(scip), &op1);
   (void) SCIPrationalCreateBuffer(SCIPbuffer(scip), &op2);

   consdataGetActivity(scip, consdata, sol, FALSE, activity);
   SCIPrationalDiff(op1, consdata->rhs, activity);
   SCIPrationalDiff(op2, activity, consdata->lhs);

   SCIPrationalMin(ret, op1, op2);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &activity);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &op1);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &op2);
}

/** creates an LP row in a linear constraint data */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   );

/** prints the certificate for a given original exact linear constraint */
SCIP_RETCODE SCIPconsPrintCertificateOrigExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_CONSDATA* consdata;
   int* varsindex;
   int i;

   /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   /* print constraint into certificate output */
   if( SCIPisCertificateActive(scip) )
   {
      certificate = SCIPgetCertificate(scip);
      consdata = SCIPconsGetData(cons);

      SCIP_CALL( SCIPallocBufferArray(scip, &varsindex, consdata->nvars) );
      for( i = 0; i < consdata->nvars; ++i )
         varsindex[i] = SCIPvarGetCertificateIndex(consdata->vars[i]);

      /* print constraint */
      if( SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
      {
         assert(!SCIPrationalIsAbsInfinity(consdata->lhs));
         SCIP_CALL( SCIPcertificatePrintCons(certificate, TRUE, NULL, 'E', consdata->lhs, consdata->nvars, varsindex, consdata->vals) );
      }
      else
      {
         if( !SCIPrationalIsNegInfinity(consdata->lhs) )
         {
            SCIP_CALL( SCIPcertificatePrintCons(certificate, TRUE, NULL, 'G', consdata->lhs, consdata->nvars, varsindex, consdata->vals) );
         }
         if( !SCIPrationalIsInfinity(consdata->rhs) )
         {
            SCIP_CALL( SCIPcertificatePrintCons(certificate, TRUE, NULL, 'L', consdata->rhs, consdata->nvars, varsindex, consdata->vals) );
         }
      }

      SCIPfreeBufferArray(scip, &varsindex);
   }

   return SCIP_OKAY;
}

/** index comparison method of linear constraints: compares two indices of the variable set in the linear constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   return SCIPvarCompare(consdata->vars[ind1], consdata->vars[ind2]);
}

/** index comparison method of linear constraints: compares two indices of the variable set in the linear constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVarProp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   var1 = consdata->vars[ind1];
   var2 = consdata->vars[ind2];

   /* exactly one variable is binary */
   if( SCIPvarIsBinary(var1) != SCIPvarIsBinary(var2) )
   {
      return (SCIPvarIsBinary(var1) ? -1 : +1);
   }
   /* both variables are binary */
   else if( SCIPvarIsBinary(var1) )
   {
      if( SCIPrationalIsAbsEQ(consdata->vals[ind1], consdata->vals[ind2]) ) {
         return (SCIPvarGetProbindex(var1) - SCIPvarGetProbindex(var2));
      }
      if( SCIPrationalIsAbsGT(consdata->vals[ind1], consdata->vals[ind2]) )
         return -1;
      else
         return +1;
   }
   else
   {
      SCIP_VARTYPE vartype1 = SCIPvarGetType(var1);
      SCIP_VARTYPE vartype2 = SCIPvarGetType(var2);

      if( vartype1 < vartype2 )
      {
         return -1;
      }
      else if( vartype1 > vartype2 )
      {
         return +1;
      }
      else
      {
         /* both variables are continuous */
         if( !SCIPvarIsIntegral(var1) )
         {
            assert(!SCIPvarIsIntegral(var2));
            return (SCIPvarGetProbindex(var1) - SCIPvarGetProbindex(var2));
         }
         else
         {
            SCIP_RATIONAL* abscont1;
            SCIP_RATIONAL* abscont2;

            (void) SCIPrationalCreate(&abscont1);
	    (void) SCIPrationalCreate(&abscont2);

            SCIPrationalDiff(abscont1, SCIPvarGetUbGlobalExact(var1), SCIPvarGetLbGlobalExact(var1));
            SCIPrationalMult(abscont1, consdata->vals[ind1], abscont1);

            SCIPrationalDiff(abscont2, SCIPvarGetUbGlobalExact(var2), SCIPvarGetLbGlobalExact(var2));
            SCIPrationalMult(abscont2, consdata->vals[ind2], abscont2);

            if( SCIPrationalIsAbsEQ(abscont1, abscont2) ) {
               SCIPrationalFree(&abscont1);
               SCIPrationalFree(&abscont2);
               return (SCIPvarGetProbindex(var1) - SCIPvarGetProbindex(var2));
            }
            if( SCIPrationalIsAbsGT(abscont2, abscont1) )
            {
               SCIPrationalFree(&abscont1);
               SCIPrationalFree(&abscont2);
               return 1;
            }
            else
            {
               SCIPrationalFree(&abscont1);
               SCIPrationalFree(&abscont2);
               return -1;
            }
         }
      }
   }
}

/** permutes the constraint's variables according to a given permutation. */
static
void permSortConsdata(
   SCIP_CONSDATA*        consdata,           /**< the constraint data */
   int*                  perm,               /**< the target permutation */
   int                   nvars               /**< the number of variables */
   )
{  /*lint --e{715}*/
   SCIP_VAR* varv;
   SCIP_EVENTDATA* eventdatav;
   SCIP_INTERVAL valrealv;
   SCIP_RATIONAL* valv;
   int v;
   int i;
   int nexti;

   assert(perm != NULL);
   assert(consdata != NULL);

   /* permute the variables in the linear constraint according to the target permutation */
   eventdatav = NULL;
   for( v = 0; v < nvars; ++v )
   {
      if( perm[v] != v )
      {
         varv = consdata->vars[v];
         valv = consdata->vals[v];
         valrealv = consdata->valsreal[v];
         if( consdata->eventdata != NULL )
            eventdatav = consdata->eventdata[v];
         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < nvars);
            assert(perm[i] != i);
            consdata->vars[i] = consdata->vars[perm[i]];
            consdata->vals[i] = consdata->vals[perm[i]];
            consdata->valsreal[i] = consdata->valsreal[perm[i]];
            if( consdata->eventdata != NULL )
            {
               consdata->eventdata[i] = consdata->eventdata[perm[i]];
               consdata->eventdata[i]->varpos = i;
            }
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->vars[i] = varv;
         consdata->vals[i] = valv;
         consdata->valsreal[i] = valrealv;
         if( consdata->eventdata != NULL )
         {
            consdata->eventdata[i] = eventdatav;
            consdata->eventdata[i]->varpos = i;
         }
         perm[i] = i;
      }
   }
#ifdef SCIP_DEBUG
   /* check sorting */
   for( v = 0; v < nvars; ++v )
   {
      assert(perm[v] == v);
      assert(consdata->eventdata == NULL || consdata->eventdata[v]->varpos == v);
   }
#endif
}

/** sorts linear constraint's variables depending on the stage of the solving process:
 * - during PRESOLVING
 *       sorts variables by binaries, integers, implicit integers, and continuous variables,
 *       and the variables of the same type by non-decreasing variable index
 *
 * - during SOLVING
 *       sorts variables of the remaining problem by binaries, integers, implicit integers, and continuous variables,
 *       and binary and integer variables by their global max activity delta (within each group),
 *       ties within a group are broken by problem index of the variable.
 *
 *       This fastens the propagation time of the constraint handler.
 */
static
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   /* check if there are variables for sorting */
   if( consdata->nvars <= 1 )
   {
      consdata->indexsorted = TRUE;
      consdata->coefsorted = TRUE;
      consdata->nbinvars = (consdata->nvars == 1 ? (int)SCIPvarIsBinary(consdata->vars[0]) : 0);
   }
   else if( (!consdata->indexsorted && SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE)
      || (!consdata->coefsorted && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE) )
   {
      int* perm;
      int v;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call sorting method  */
      if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
         SCIPsort(perm, consdataCompVar, (void*)consdata, consdata->nvars);
      else
         SCIPsort(perm, consdataCompVarProp, (void*)consdata, consdata->nvars);

      permSortConsdata(consdata, perm, consdata->nvars);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);

      if( SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE )
      {
         consdata->indexsorted = FALSE;
         consdata->coefsorted = TRUE;

         /* count binary variables in the sorted vars array */
         consdata->nbinvars = 0;
         for( v = 0; v < consdata->nvars; ++v )
         {
            if( SCIPvarIsBinary(consdata->vars[v]) )
               ++consdata->nbinvars;
            else
               break;
         }
      }
      else
      {
         consdata->indexsorted = TRUE;
         consdata->coefsorted = FALSE;
      }
   }

   return SCIP_OKAY;
}


/*
 * local linear constraint handler methods
 */

/** sets left hand side of linear constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_RATIONAL*        lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool locked;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPrationalIsInfinity(lhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPrationalIsInfinity(consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPrationalIsEQ(consdata->lhs, lhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPrationalIsEQ(lhs, consdata->rhs) )
   {
      SCIPrationalSetRational(consdata->rhs, lhs);
      assert(consdata->rowlhs == NULL);
   }

   locked = FALSE;
   for( i = 0; i < NLOCKTYPES && !locked; i++ )
      locked = SCIPconsIsLockedType(cons, (SCIP_LOCKTYPE) i);

   /* if necessary, update the rounding locks of variables */
   if( locked )
   {
      if( SCIPrationalIsNegInfinity(consdata->lhs) && !SCIPrationalIsNegInfinity(lhs) )
      {
         SCIP_VAR** vars;
         SCIP_RATIONAL** vals;
         int v;

         /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPrationalIsZero(vals[v]));

            if( SCIPrationalIsPositive(vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
      else if( !SCIPrationalIsNegInfinity(consdata->lhs) && SCIPrationalIsNegInfinity(lhs) )
      {
         SCIP_VAR** vars;
         SCIP_RATIONAL** vals;
         int v;

         /* the left hand side switched from a non-infinite value to -infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPrationalIsZero(vals[v]));

            if( SCIPrationalIsPositive(vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
   }

   /* check whether the left hand side is increased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( !SCIPrationalIsNegInfinity(lhs) && SCIPrationalIsGT(lhs, consdata->lhs) )
   {
      consdata->boundstightened = 0;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
      consdata->implsadded = FALSE;

      /* mark the constraint for propagation */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      }
   }

   /* set new left hand side and update constraint data */
   SCIPrationalSetRational(consdata->lhs, lhs);
   consdata->lhsreal = SCIPrationalRoundReal(lhs, SCIP_R_ROUND_DOWNWARDS);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->rangedrowpropagated = 0;

   /* update the lhs of the LP row */
   if( consdata->rowexact != NULL )
   {
      SCIP_CALL( SCIPchgRowExactLhs(scip, consdata->rowexact, lhs) );
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_RATIONAL*        rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool locked;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPrationalIsNegInfinity(rhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPrationalIsNegInfinity(consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPrationalIsEQ(consdata->rhs, rhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPrationalIsEQ(rhs, consdata->lhs) )
   {
      SCIPrationalSetRational(consdata->rhs, rhs);
      assert(consdata->rowlhs == NULL);
   }

   locked = FALSE;
   for( i = 0; i < NLOCKTYPES && !locked; i++ )
      locked = SCIPconsIsLockedType(cons, (SCIP_LOCKTYPE) i);

   /* if necessary, update the rounding locks of variables */
   if( locked )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPrationalIsInfinity(consdata->rhs) && !SCIPrationalIsInfinity(rhs) )
      {
         SCIP_VAR** vars;
         SCIP_RATIONAL** vals;
         int v;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPrationalIsZero(vals[v]));

            if( SCIPrationalIsPositive(vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
      else if( !SCIPrationalIsInfinity(consdata->rhs) && SCIPrationalIsInfinity(rhs) )
      {
         SCIP_VAR** vars;
         SCIP_RATIONAL** vals;
         int v;

         /* the right hand side switched from a non-infinite value to infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPrationalIsZero(vals[v]));

            if( SCIPrationalIsPositive(vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
   }

   /* check whether the right hand side is decreased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( !SCIPrationalIsInfinity(rhs) && SCIPrationalIsLT(rhs, consdata->rhs) )
   {
      consdata->boundstightened = 0;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
      consdata->implsadded = FALSE;

      /* mark the constraint for propagation */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      }
   }

   /* set new right hand side and update constraint data */
   SCIPrationalSetRational(consdata->rhs, rhs);
   consdata->rhsreal = SCIPrationalRoundReal(rhs, SCIP_R_ROUND_UPWARDS);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->rangedrowpropagated = 0;

   /* update the rhs of the LP row */
   if( consdata->rowexact != NULL )
   {
      SCIP_CALL( SCIPchgRowExactRhs(scip, consdata->rowexact, rhs) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPrationalIsZero(val) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   SCIPrationalSetRational(consdata->vals[consdata->nvars], val);
   SCIPintervalSetRational(&(consdata->valsreal[consdata->nvars]), val);
   consdata->nvars++;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, the variable needs an additional event data */
   if( transformed )
   {
      if( consdata->eventdata != NULL )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONSHDLRDATA* conshdlrdata;

         /* check for event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* initialize eventdata array */
         consdata->eventdata[consdata->nvars-1] = NULL;

         /* catch bound change events of variable */
         SCIP_CALL( consCatchEvent(scip, cons, conshdlrdata->eventhdlr, consdata->nvars-1) );
      }

      /* update minimum and maximum activities */
      consdataUpdateAddCoef(scip, consdata, var, consdata->vals[consdata->nvars - 1], consdata->valsreal[consdata->nvars - 1]);
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockRounding(scip, cons, var, val) );

   /* mark the constraint for propagation */
   if( transformed )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->removedfixings = consdata->removedfixings && SCIPvarIsActive(var);

   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   if( consdata->nvars == 1 )
   {
     consdata->indexsorted = TRUE;
     consdata->coefsorted = TRUE;
     consdata->merged = TRUE;
   }
   else
   {
      consdata->merged = FALSE;

      if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
      {
         consdata->indexsorted = consdata->indexsorted && (consdataCompVar((void*)consdata, consdata->nvars-2, consdata->nvars-1) <= 0);
         consdata->coefsorted = FALSE;
      }
      else
      {
         consdata->indexsorted = FALSE;
         consdata->coefsorted = consdata->coefsorted && (consdataCompVarProp((void*)consdata, consdata->nvars-2, consdata->nvars-1) <= 0);
      }
   }

   /* update hascontvar and hasnonbinvar flags */
   if( consdata->hasnonbinvalid && !consdata->hascontvar )
   {
      SCIP_VARTYPE vartype = SCIPvarGetType(var);

      if( vartype != SCIP_VARTYPE_BINARY )
      {
         consdata->hasnonbinvar = TRUE;

         if( vartype == SCIP_VARTYPE_CONTINUOUS )
            consdata->hascontvar = TRUE;
      }
   }

   /* add the new coefficient to the LP row */
   if( consdata->rowexact != NULL )
   {
     SCIP_CALL( SCIPaddVarsToRowExact(scip, consdata->rowexact, 1, &var, &val) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from linear constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_RATIONAL* val;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, var, val) );

   /* if we are in transformed problem, delete the event data of the variable */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* check for event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      if( consdata->eventdata != NULL )
      {
         SCIP_CALL( consDropEvent(scip, cons, conshdlrdata->eventhdlr, pos) );
         assert(consdata->eventdata[pos] == NULL);
      }
   }

   /* move the last variable to the free slot */
   if( pos != consdata->nvars - 1 )
   {
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      SCIPrationalSetRational(consdata->vals[pos], consdata->vals[consdata->nvars - 1]);
      consdata->valsreal[pos] = consdata->valsreal[consdata->nvars -1];

      if( consdata->eventdata != NULL )
      {
         consdata->eventdata[pos] = consdata->eventdata[consdata->nvars-1];
         assert(consdata->eventdata[pos] != NULL);
         consdata->eventdata[pos]->varpos = pos;
      }

      consdata->indexsorted = consdata->indexsorted && (pos + 2 >= consdata->nvars);
      consdata->coefsorted = consdata->coefsorted && (pos + 2 >= consdata->nvars);
   }
   consdata->nvars--;

   /* mark the constraint for propagation */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   /* check if hasnonbinvar flag might be incorrect now */
   if( consdata->hasnonbinvar && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
   {
      consdata->hasnonbinvalid = FALSE;
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** changes coefficient value at given position of linear constraint data */
static
SCIP_RETCODE chgCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of coefficient to delete */
   SCIP_RATIONAL*        newval              /**< new value of coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_RATIONAL* val;
   SCIP_Bool locked;
   SCIP_INTERVAL newvalfp;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPrationalIsZero(newval));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   locked = FALSE;
   for( i = 0; i < NLOCKTYPES && !locked; i++ )
      locked = SCIPconsIsLockedType(cons, (SCIP_LOCKTYPE) i);

   /* if necessary, update the rounding locks of the variable */
   if( locked && ((SCIPrationalIsNegative(newval) && SCIPrationalIsPositive(val)) || (SCIPrationalIsNegative(val) && SCIPrationalIsPositive(newval))) )
   {
      assert(SCIPconsIsTransformed(cons));

      /* remove rounding locks for variable with old coefficient */
      SCIP_CALL( unlockRounding(scip, cons, var, val) );

      /* install rounding locks for variable with new coefficient */
      SCIP_CALL( lockRounding(scip, cons, var, newval) );
   }
   SCIPintervalSetRational(&newvalfp, newval);
   /* update minimum and maximum activities */
   if( SCIPconsIsTransformed(cons) )
      consdataUpdateChgCoef(scip, consdata, var, consdata->valsreal[pos], val, newvalfp, newval);

      /* change the value */
   SCIPrationalSetRational(consdata->vals[pos], newval);
   consdata->valsreal[pos] = newvalfp;
   if( consdata->coefsorted )
   {
      if( pos > 0 )
         consdata->coefsorted = (consdataCompVarProp((void*)consdata, pos - 1, pos) <= 0);
      if( consdata->coefsorted && pos < consdata->nvars - 1 )
         consdata->coefsorted = (consdataCompVarProp((void*)consdata, pos, pos + 1) <= 0);
   }
   /* mark the constraint for propagation */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   return SCIP_OKAY;
}

/* perform deletion of variables in all constraints of the constraint handler */
static
SCIP_RETCODE performVarDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int v;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* iterate over all constraints */
   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      /* constraint is marked, that some of its variables were deleted */
      if( consdata->varsdeleted )
      {
         /* iterate over all variables of the constraint and delete them from the constraint */
         for( v = consdata->nvars - 1; v >= 0; --v )
         {
            if( SCIPvarIsDeleted(consdata->vars[v]) )
            {
               SCIP_CALL( delCoefPos(scip, conss[i], v) );
            }
         }
         consdata->varsdeleted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_RATIONAL* valsum;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &valsum) );

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata) );

   /* go backwards through the constraint looking for multiple occurrences of the same variable;
    * backward direction is necessary, since delCoefPos() modifies the given position and
    * the subsequent ones
    */
   v = consdata->nvars-1;
   while( v >= 1 )
   {
      var = consdata->vars[v];
      if( consdata->vars[v-1] == var )
      {
         SCIPrationalSetRational(valsum, consdata->vals[v]);
         do
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );
            --v;
            SCIPrationalAdd(valsum, valsum, consdata->vals[v]);
         }
         while( v >= 1 && consdata->vars[v-1] == var );

         /* modify the last existing occurrence of the variable */
         assert(consdata->vars[v] == var);
         if( SCIPrationalIsZero(valsum) )
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );

            /* if the variable defining the maximal activity delta was removed from the constraint, the maximal activity
             * delta needs to be recalculated on the next real propagation
             */
            if( consdata->maxactdeltavar == var )
            {
               consdata->maxactdelta = SCIP_INVALID;
               consdata->maxactdeltavar = NULL;
            }
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, v, valsum) );
         }
      }
      --v;
   }

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &valsum);
   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables by their non-fixed counterparts */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility is detected; or NULL if this
                                              *   information is not needed; in this case, we apply all fixings
                                              *   instead of stopping after the first infeasible one */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_RATIONAL* val;
   SCIP_RATIONAL** aggrscalars;
   SCIP_RATIONAL* fixedval;
   SCIP_RATIONAL* aggrconst;
   SCIP_Real negconst;
   int v;
   int naggrvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->eventdata == NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consCatchAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      assert(consdata->eventdata != NULL);
   }

   if( !consdata->removedfixings )
   {
      SCIP_RATIONAL* lhssubtrahend;
      SCIP_RATIONAL* rhssubtrahend;
      SCIP_RATIONAL* tmpval;

      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &lhssubtrahend) );
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &rhssubtrahend) );
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmpval) );

      SCIPdebugMsg(scip, "applying fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      v = 0;
      while( v < consdata->nvars )
      {
         var = consdata->vars[v];
         val = consdata->vals[v];
         assert(SCIPvarIsTransformed(var));

         switch( SCIPvarGetStatus(var) )
         {
         case SCIP_VARSTATUS_ORIGINAL:
            SCIPerrorMessage("original variable in transformed linear constraint\n");
            return SCIP_INVALIDDATA;

         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
            ++v;
            break;

         case SCIP_VARSTATUS_FIXED:
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
            fixedval = SCIPvarGetLbGlobalExact(var);
            if( !SCIPrationalIsNegInfinity(consdata->lhs) )
            {
               if( SCIPrationalIsAbsInfinity(fixedval) )
               {
                  if( SCIPrationalGetSign(val) == SCIPrationalGetSign(fixedval) )
                  {
                     SCIPrationalSetNegInfinity(tmpval);
                     SCIP_CALL( chgLhs(scip, cons, tmpval) );
                  }
                  else
                  {
                     if( infeasible != NULL )
                     {
                        /* if lhs gets infinity it means that the problem is infeasible */
                        *infeasible = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        SCIPrationalSetInfinity(tmpval);
                        SCIP_CALL( chgLhs(scip, cons, tmpval) );
                     }
                  }
               }
               else
                  SCIPrationalAddProd(lhssubtrahend, val, fixedval);
            }
            if( !SCIPrationalIsInfinity(consdata->rhs) )
            {
               if( SCIPrationalIsAbsInfinity(fixedval) )
               {
                  if( SCIPrationalGetSign(val) == SCIPrationalGetSign(fixedval) )
                  {
                     if( infeasible != NULL )
                     {
                        /* if rhs gets -infinity it means that the problem is infeasible */
                        *infeasible = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        SCIPrationalSetNegInfinity(tmpval);
                        SCIP_CALL( chgRhs(scip, cons, tmpval) );
                     }
                  }
                  else
                  {
                     SCIPrationalSetInfinity(tmpval);
                     SCIP_CALL( chgRhs(scip, cons, tmpval) );
                  }
               }
               else
                  SCIPrationalAddProd(rhssubtrahend, val, fixedval);
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_AGGREGATED:
         {
            SCIP_VAR* activevar = SCIPvarGetAggrVar(var);
            SCIP_RATIONAL* activescalar;
            SCIP_RATIONAL* activeconstant;

            SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &activescalar) );
            SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &activeconstant) );

            SCIPrationalMult(activescalar, val, SCIPvarGetAggrScalarExact(var));
            SCIPrationalMult(activeconstant, val, SCIPvarGetAggrConstantExact(var));

            assert(activevar != NULL);
            SCIP_CALL( SCIPgetProbvarSumExact(scip, &activevar, activescalar, activeconstant) );
            assert(activevar != NULL);

            if( !SCIPrationalIsZero(activescalar) )
            {
               SCIP_CALL( addCoef(scip, cons, activevar, activescalar) );
            }

            if( !SCIPrationalIsZero(activeconstant) )
            {
               if( !SCIPrationalIsNegInfinity(consdata->lhs) )
                  SCIPrationalAdd(lhssubtrahend, lhssubtrahend, activeconstant);
               if( !SCIPrationalIsInfinity(consdata->rhs) )
                  SCIPrationalAdd(rhssubtrahend, rhssubtrahend, activeconstant);
            }

            SCIP_CALL( delCoefPos(scip, cons, v) );

            SCIPrationalFreeBuffer(SCIPbuffer(scip), &activescalar);
            SCIPrationalFreeBuffer(SCIPbuffer(scip), &activeconstant);
            break;
         }
         case SCIP_VARSTATUS_MULTAGGR:
            SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );
            naggrvars = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalarsExact(var);
            for( i = 0; i < naggrvars; ++i )
            {
               SCIPrationalMult(tmpval, val, aggrscalars[i]);
               SCIP_CALL( addCoef(scip, cons, aggrvars[i], tmpval) );
            }
            aggrconst = SCIPvarGetMultaggrConstantExact(var);

            if( !SCIPrationalIsNegInfinity(consdata->lhs) )
            {
               SCIPrationalMult(tmpval, val, aggrconst);
               SCIPrationalAdd(lhssubtrahend, lhssubtrahend, tmpval);
            }
            if( !SCIPrationalIsInfinity(consdata->rhs) )
            {
               SCIPrationalMult(tmpval, val, aggrconst);
               SCIPrationalAdd(rhssubtrahend, rhssubtrahend, tmpval);
            }

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_NEGATED:
            SCIPrationalNegate(tmpval, val);
            SCIP_CALL( addCoef(scip, cons, SCIPvarGetNegationVar(var), tmpval) );
            negconst = SCIPvarGetNegationConstant(var);

            if( !SCIPrationalIsNegInfinity(consdata->lhs) )
            {
               SCIPrationalMultReal(tmpval, val, negconst);
               SCIPrationalAdd(lhssubtrahend, lhssubtrahend, tmpval);
            }
            if( !SCIPrationalIsInfinity(consdata->rhs) )
            {
               SCIPrationalMultReal(tmpval, val, negconst);
               SCIPrationalAdd(rhssubtrahend, rhssubtrahend, tmpval);
            }

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         default:
            SCIPerrorMessage("unknown variable status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA;  /*lint !e527*/
         }
      }

      if( !SCIPrationalIsAbsInfinity(consdata->lhs) )
      {
         SCIPrationalDiff(tmpval, consdata->lhs, lhssubtrahend);
         SCIP_CALL( chgLhs(scip, cons, tmpval) );
      }
      if( !SCIPrationalIsAbsInfinity(consdata->rhs) )
      {
         SCIPrationalDiff(tmpval, consdata->rhs, rhssubtrahend);
         SCIP_CALL( chgRhs(scip, cons, tmpval) );
      }

      consdata->removedfixings = TRUE;

      SCIPdebugMsg(scip, "after fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
       * to clean up the constraint
       */
      SCIP_CALL( mergeMultiples(scip, cons) );

      SCIPdebugMsg(scip, "after merging:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmpval);
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &rhssubtrahend);
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &lhssubtrahend);
   }
   assert(consdata->removedfixings);

#ifndef NDEBUG
   /* check, if all fixings are applied */
   for( v = 0; v < consdata->nvars; ++v )
      assert(SCIPvarIsActive(consdata->vars[v]));
#endif

   return SCIP_OKAY;
}

/** prints activity conflict to  certificate file */
static
SCIP_RETCODE certificatePrintActivityConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool rhs
   )
{
   SCIP_Real side;
   SCIP_Real activity;
   SCIP_RATIONAL* diff;
   int nvals;
   SCIP_RATIONAL** vals;

   if( !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;
   SCIP_CALL(SCIPrationalCreateBuffer(SCIPbuffer(scip), &diff));

   if( rhs )
   {
      consdataRecomputeMinactivity(scip, consdata);
      side = consdata->rhsreal;
      activity = consdata->minactivity;
      assert( activity > side );
   }
   else
   {
      consdataRecomputeMaxactivity(scip, consdata);
      side = consdata->lhsreal;
      activity = consdata->maxactivity;
      assert( activity < side );
   }

   if( consdata->rowexact != NULL )
   {
      nvals = SCIProwExactGetNNonz(consdata->rowexact);
      vals = SCIProwExactGetVals(consdata->rowexact);
   }
   else
   {
      nvals = consdata->nvars;
      vals = consdata->vals;
   }
   SCIPrationalSetReal(diff, activity);
   SCIPrationalDiffReal(diff, diff, side);

   SCIP_CALL( SCIPcertificatePrintActivityConflict(scip, cons, consdata->rowexact, consdata->lhs, consdata->rhs,
      nvals, vals, consdata->vars, diff, rhs) );

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &diff);

   return SCIP_OKAY;
}

/** tightens bounds of a single variable due to activity bounds */
static
SCIP_RETCODE tightenVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of the variable in the vars array */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool             force               /**< should a possible bound change be forced even if below bound strengthening tolerance */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_INTERVAL valrange;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool minisrelax;
   SCIP_Bool maxisrelax;
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_ROUNDMODE prevmode;
   SCIP_RATIONAL* tmpbound;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   prevmode = SCIPintervalGetRoundingMode();

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      goto RETURN_SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   var = consdata->vars[pos];

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
   {
      return SCIP_OKAY;
   }
   else
   {
      SCIP_VAR* tmpVar;
      SCIP_Real tmpBound;
      SCIP_BOUNDTYPE tmpBoundtype;
      tmpVar = var;
      SCIP_CALL( SCIPvarGetProbvarBound(&tmpVar, &tmpBound, &tmpBoundtype) );
      if( SCIPvarGetStatus(tmpVar) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(tmpVar) == SCIP_VARSTATUS_FIXED ) {
         goto RETURN_SCIP_OKAY;
      }
   }

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && !conshdlrdata->propcont )
      return SCIP_OKAY;

   valrange = consdata->valsreal[pos];
   lhs = consdata->lhsreal;
   rhs = consdata->rhsreal;
   consdataGetActivityResiduals(scip, consdata, var, valrange, FALSE, &minresactivity, &maxresactivity,
      &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);
   assert(var != NULL);
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPisLE(scip, lb, ub));

   if( valrange.sup > 0.0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) && !minisrelax )
      {
         SCIP_Real newub;
         SCIP_INTERVAL ubinterval;

         SCIPintervalSetRoundingModeUpwards();

         /* newub = (rhs + SCIPintervalNegateReal(minresactivity))/valrange.inf; */
         SCIPintervalSet(&ubinterval, rhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &ubinterval, ubinterval, minresactivity);
         SCIPintervalDiv(SCIPinfinity(scip), &ubinterval, ubinterval, valrange);
         newub = ubinterval.sup;

         if( !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            /* activity is never unreliable in exact solving */

            /* tighten upper bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newub=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, valrange.inf, minresactivity, maxresactivity, lhs, rhs, newub);

            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               SCIP_Longint boundmaxdenom;

               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmpbound) );
               SCIPrationalSetReal(tmpbound, newub);

               if( conshdlrdata->limitdenom )
               {
                  boundmaxdenom = conshdlrdata->boundmaxdenom;
                  SCIPrationalComputeApproximation(tmpbound, tmpbound, boundmaxdenom, 1);
               }

               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBoundEx(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_UPPER, tmpbound, false, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars,  consdata->nvars) );

               SCIP_CALL( SCIPinferVarUbConsExact(scip, var, tmpbound, cons, getInferInt(PROPRULE_1_RHS, pos),
                     &infeasible, &tightened) );
               SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmpbound);
            }
            else
            {
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBound(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_UPPER, newub, false, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               newub = SCIPadjustedVarUbExactFloat(scip, var, newub);
               SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
            }

            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               /* analyze conflict */
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, TRUE) );
               *cutoff = TRUE;
               goto RETURN_SCIP_OKAY;
            }
            if( tightened )
            {
               ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               (*nchgbds)++;

               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newlb;
         SCIP_INTERVAL lbinterval;

         SCIPintervalSetRoundingModeDownwards();
         /* newlb = (lhs + SCIPintervalNegateReal(maxresactivity))/valrange.sup; */
         SCIPintervalSet(&lbinterval, lhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &lbinterval, lbinterval, maxresactivity);
         SCIPintervalDiv(SCIPinfinity(scip), &lbinterval, lbinterval, valrange);
         newlb = lbinterval.inf;

         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            /* tighten lower bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, valrange.inf, minresactivity, maxresactivity, lhs, rhs, newlb);

            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               SCIP_Longint boundmaxdenom;

               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmpbound) );
               SCIPrationalSetReal(tmpbound, newlb);

               if( conshdlrdata->limitdenom )
               {
                  boundmaxdenom = conshdlrdata->boundmaxdenom;
                  SCIPrationalComputeApproximation(tmpbound, tmpbound, boundmaxdenom, -1);
               }
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBoundEx(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_LOWER, tmpbound, true, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               SCIP_CALL( SCIPinferVarLbConsExact(scip, var, tmpbound, cons, getInferInt(PROPRULE_1_LHS, pos),
                     &infeasible, &tightened) );
               SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmpbound);
            }
            else
            {
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBound(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_LOWER, newlb, true, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               newlb = SCIPadjustedVarLbExactFloat(scip, var, newlb);
               SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                     &infeasible, &tightened) );
            }

            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, FALSE) );

               *cutoff = TRUE;
               goto RETURN_SCIP_OKAY;
            }
            if( tightened )
            {
               (*nchgbds)++;
               SCIPdebug(lb = SCIPvarGetLbLocal(var)); /* get bound again: it may be additionally modified due to integrality */
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's bounds */
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) && !minisrelax )
      {
         SCIP_Real newlb;
         SCIP_INTERVAL lbinterval;

         SCIPintervalSetRoundingModeDownwards();

         SCIPintervalSet(&lbinterval, rhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &lbinterval, lbinterval, minresactivity);
         SCIPintervalDiv(SCIPinfinity(scip), &lbinterval, lbinterval, valrange);
         newlb = lbinterval.inf;

         assert(newlb <= lbinterval.inf);

         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            /* tighten lower bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, valrange.sup, minresactivity, maxresactivity, lhs, rhs, newlb);

            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               SCIP_Longint boundmaxdenom;

               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmpbound) );
               SCIPrationalSetReal(tmpbound, newlb);

               if( conshdlrdata->limitdenom )
               {
                  boundmaxdenom = conshdlrdata->boundmaxdenom;
                  SCIPrationalComputeApproximation(tmpbound, tmpbound, boundmaxdenom, -1);
               }
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBoundEx(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_LOWER, tmpbound, false, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               SCIP_CALL( SCIPinferVarLbConsExact(scip, var, tmpbound, cons, getInferInt(PROPRULE_1_RHS, pos),
                     &infeasible, &tightened) );
               SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmpbound);
            }
            else
            {
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBound(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_LOWER, newlb, false, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               newlb = SCIPadjustedVarLbExactFloat(scip, var, newlb);
               SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
            }

            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, TRUE) );

               /**@todo analyze conflict detected in exactlinear constraint handler */
               *cutoff = TRUE;
               goto RETURN_SCIP_OKAY;
            }
            if( tightened )
            {
               lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               (*nchgbds)++;
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newub;
         SCIP_INTERVAL ubinterval;

         SCIPintervalSetRoundingModeUpwards();

         /* newub = (maxresactivity + SCIPintervalNegateReal(lhs))/SCIPintervalNegateReal(valrange.inf); */
         SCIPintervalSet(&ubinterval, lhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &ubinterval, ubinterval, maxresactivity);
         SCIPintervalDiv(SCIPinfinity(scip), &ubinterval, ubinterval, valrange);
         newub = ubinterval.sup;

         if( !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            /* tighten upper bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g], newub=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, valrange.sup, minresactivity, maxresactivity, lhs, rhs, newub);

            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               SCIP_Longint boundmaxdenom;

               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmpbound) );
               SCIPrationalSetReal(tmpbound, newub);

               if( conshdlrdata->limitdenom )
               {
                  boundmaxdenom = conshdlrdata->boundmaxdenom;
                  SCIPrationalComputeApproximation(tmpbound, tmpbound, boundmaxdenom, 1);
               }
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBoundEx(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_UPPER, tmpbound, true, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               SCIP_CALL( SCIPinferVarUbConsExact(scip, var, tmpbound, cons, getInferInt(PROPRULE_1_LHS, pos),
                     &infeasible, &tightened) );
               SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmpbound);
            }
            else
            {
               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( SCIPcertificatePrintActivityVarBound(scip, SCIPgetCertificate(scip), NULL,
                     SCIP_BOUNDTYPE_UPPER, newub, true, cons, var, consdata->rowexact, consdata->vals, consdata->lhs, consdata->rhs, consdata->vars, consdata->nvars) );

               newub = SCIPadjustedVarUbExactFloat(scip, var, newub);
               SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                     &infeasible, &tightened) );
            }

            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               if( SCIPcertificateShouldTrackBounds(scip) )
                  SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, FALSE) );

               *cutoff = TRUE;
               goto RETURN_SCIP_OKAY;
            }
            if( tightened )
            {
               (*nchgbds)++;
               SCIPdebug(ub = SCIPvarGetUbLocal(var)); /* get bound again: it may be additionally modified due to integrality */
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }
   RETURN_SCIP_OKAY:
   SCIPintervalSetRoundingMode(prevmode);
   return SCIP_OKAY;
}

#define MAXTIGHTENROUNDS 10

/** tightens bounds of variables in constraint due to activity bounds */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool             sortvars,           /**< should variables be used in sorted order? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   unsigned int tightenmode;
   int nvars;
   int nrounds;
   int lastchange;
   int oldnchgbds;
#ifndef SCIP_DEBUG
   int oldnchgbdstotal;
#endif
   int v;
   SCIP_Bool force;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* currently, we do not need to call applyFixings() as in cons_linear.c */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   force = (nvars == 1) && !SCIPconsIsModifiable(cons);

   /* we are at the root node or during presolving */
   if( SCIPgetDepth(scip) < 1 )
      tightenmode = 2;
   else
      tightenmode = 1;

   /* stop if we already tightened the constraint and the tightening is not forced */
   if( !force && (consdata->boundstightened >= tightenmode) ) /*lint !e574*/
      return SCIP_OKAY;

   /* ensure that the variables are properly sorted */
   if( sortvars && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->coefsorted )
   {
      SCIP_CALL( consdataSort(scip, consdata) );
      assert(consdata->coefsorted);
   }

   /* update maximal activity delta if necessary */
   if( consdata->maxactdelta == SCIP_INVALID ) /*lint !e777*/
      consdataRecomputeMaxActivityDelta(scip, consdata);

   assert(consdata->maxactdelta != SCIP_INVALID); /*lint !e777*/
   assert(!SCIPisFeasNegative(scip, consdata->maxactdelta));
   checkMaxActivityDelta(scip, consdata);

   /* this may happen if all variables are fixed */
   if( SCIPisFeasZero(scip, consdata->maxactdelta) )
      return SCIP_OKAY;

   if( !SCIPisInfinity(scip, consdata->maxactdelta) )
   {
      SCIP_Real slack;
      SCIP_Real surplus;
      SCIP_Real minactivity;
      SCIP_Real maxactivity;
      SCIP_Bool minisrelax;
      SCIP_Bool maxisrelax;
      SCIP_Bool isminsettoinfinity;
      SCIP_Bool ismaxsettoinfinity;

      /* use maximal activity delta to skip propagation (cannot deduce anything) */
      consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minisrelax, &maxisrelax,
         &isminsettoinfinity, &ismaxsettoinfinity);

      assert(!SCIPisInfinity(scip, minactivity));
      assert(!SCIPisInfinity(scip, -maxactivity));

      slack = (SCIPisInfinity(scip, consdata->rhsreal) || isminsettoinfinity) ? SCIPinfinity(scip) : (consdata->rhsreal - minactivity);
      surplus = (SCIPisInfinity(scip, -consdata->lhsreal) || ismaxsettoinfinity) ? SCIPinfinity(scip) : (maxactivity - consdata->lhsreal);

      /* check if the constraint will propagate */
      if( consdata->maxactdelta <= MIN(slack, surplus) )
         return SCIP_OKAY;
   }

   /* as long as the bounds might be tightened again, try to tighten them; abort after a maximal number of rounds */
   lastchange = -1;
   oldnchgbds = 0;

#ifndef SCIP_DEBUG
   oldnchgbdstotal = *nchgbds;
#endif

   for( nrounds = 0; (force || consdata->boundstightened < tightenmode) && nrounds < MAXTIGHTENROUNDS; ++nrounds ) /*lint !e574*/
   {
      /* ensure that the variables are properly sorted
       *
       * note: it might happen that integer variables become binary during bound tightening at the root node
       */
      if( sortvars && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->coefsorted )
      {
         SCIP_CALL( consdataSort(scip, consdata) );
         assert(consdata->coefsorted);
      }

      /* mark the constraint to have the variables' bounds tightened */
      consdata->boundstightened = (unsigned int)tightenmode;
      /* try to tighten the bounds of each variable in the constraint. During solving process, the binary variable
       * sorting enables skipping variables
       */
      v = 0;
      while( v < nvars && v != lastchange && !(*cutoff) )
      {
         oldnchgbds = *nchgbds;

         SCIP_CALL( tightenVarBounds(scip, cons, v, cutoff, nchgbds, force) );


         /* if there was no progress, skip the rest of the binary variables */
         if( *cutoff )
         {
            break;
         }
         else if( *nchgbds > oldnchgbds )
         {
            lastchange = v;
            ++v;
         }
         else if( consdata->coefsorted && v < consdata->nbinvars - 1
            && !SCIPisFeasEQ(scip, SCIPvarGetUbLocal(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v])) )
            v = consdata->nbinvars;
         else
            ++v;
      }

#ifndef SCIP_DEBUG
      SCIPdebugMessage("linear constraint <%s> found %d bound changes in round %d\n", SCIPconsGetName(cons),
         *nchgbds - oldnchgbdstotal, nrounds);
      oldnchgbdstotal += oldnchgbds;
#endif
   }

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or current solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool             useexactsol,        /**< should the sol or solex be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_RATIONAL* activity;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMsg(scip, "checking linear constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug(consPrintConsSol(scip, cons, sol, useexactsol, NULL));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;
   activity = consdata->activity;

   /* only check exact constraint if fp cons is feasible enough */
   if( (consdata->rowexact == NULL || checklprows) && !SCIPrationalIsEQ(consdata->lhs, consdata->rhs) )
   {
      SCIP_Real activityfp;
      SCIP_Real mu;

      success = consdataComputeSolActivityWithErrorbound(scip, consdata, sol, &activityfp, &mu);

      conshdlrdata->ncheckserrorbound++;

      if( !success )
         conshdlrdata->nabotserrorbound++;

      if( success )
      {
         if( activityfp - mu > consdata->rhsreal || activityfp + mu < consdata->lhsreal )
         {
            SCIPdebugMsg(scip, "discarding solution due to fp check: activityfp=%g, lhsreal=%g, rhsreal=%g, mu=%g\n",
               activityfp, consdata->lhsreal, consdata->rhsreal, mu);
            *violated = TRUE;
            conshdlrdata->nsuccesserrorbound++;
            return SCIP_OKAY;
         }
         else if( activityfp + mu < consdata->rhsreal && activityfp - mu >= consdata->lhsreal )
         {
            SCIPdebugMsg(scip, "skipping exact check due to fp check: activityfp=%g, lhsreal=%g, rhsreal=%g, mu=%g\n",
               activityfp, consdata->lhsreal, consdata->rhsreal, mu);
            *violated = FALSE;
            conshdlrdata->nsuccesserrorbound++;
            return SCIP_OKAY;
         }
         else
         {
            SCIPdebugMsg(scip, "no decision due to fp check: activityfp=%g, lhsreal=%g, rhsreal=%g, mu=%g\n",
               activityfp, consdata->lhsreal, consdata->rhsreal, mu);
         }
      }
   }

   if( consdata->rowexact != NULL )
   {
      if( !checklprows && SCIProwExactIsInLP(consdata->rowexact) && SCIPlpExactIsSolved(scip) )
         return SCIP_OKAY;
      else if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
         consdataComputePseudoActivity(scip, consdata, activity);
      else
      {
         SCIP_CALL( SCIPgetRowSolActivityExact(scip, consdata->rowexact, sol, useexactsol, activity) );
      }
   }
   else
      consdataGetActivity(scip, consdata, sol, useexactsol, activity);

   SCIPrationalDebugMessage("consdata activity=%q (lhs=%q, rhs=%q, row=%p, checklprows=%u, rowinlp=%u, sol=%p, hascurrentnodelp=%u)\n",
      activity, consdata->lhs, consdata->rhs, (void*)consdata->rowexact, checklprows,
      consdata->rowexact == NULL ? 0 : SCIProwExactIsInLP(consdata->rowexact), (void*)sol,
      consdata->rowexact == NULL ? FALSE : SCIPhasCurrentNodeLP(scip));

   /* the activity of pseudo solutions may be invalid if it comprises positive and negative infinity contributions; we
    * return infeasible for safety
    */
   if( ((!SCIPrationalIsNegInfinity(consdata->lhs) && SCIPrationalIsLT(activity, consdata->lhs)) ||
      (!SCIPrationalIsInfinity(consdata->rhs) && SCIPrationalIsGT(activity, consdata->rhs))) )
   {
      *violated = TRUE;

      /* only reset constraint age if we are in enforcement */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }
   else
   {
      /* only increase constraint age if we are in enforcement */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** creates an LP row in a linear constraint data */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool onerowrelax;
   SCIP_Bool hasfprelax;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   assert(consdata != NULL);
   assert(consdata->rowexact == NULL);

   /** create empty fp-rows */
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(consdata->rowrhs), cons, SCIPconsGetName(cons), -SCIPinfinity(scip), SCIPinfinity(scip),
      SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(consdata->rowlhs), cons, SCIPconsGetName(cons), -SCIPinfinity(scip), SCIPinfinity(scip),
      SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   /** create exact row */
   SCIP_CALL( SCIPcreateEmptyRowConsExact(scip, &consdata->rowexact, consdata->rowlhs, consdata->rowrhs,
      consdata->lhs, consdata->rhs, consdata->hasfprelax) );

   SCIP_CALL( SCIPcaptureRowExact(scip, consdata->rowexact) );

   SCIP_CALL( SCIPaddVarsToRowExact(scip, consdata->rowexact, consdata->nvars, consdata->vars, consdata->vals) );

   onerowrelax = TRUE;
   hasfprelax = TRUE;

   SCIP_CALL( SCIPgenerateFpRowsFromRowExact(scip, consdata->rowexact, consdata->rowlhs,
      consdata->rowrhs, &onerowrelax, &hasfprelax) );

   consdata->onerowrelax = onerowrelax;
   consdata->hasfprelax = hasfprelax;
   consdataInvalidateActivities(consdata);
   if( !(consdata->hasfprelax) || consdata->onerowrelax )
      consdata->rowrhs = NULL;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->roweventdata, SCIProwExactGetNNonz(consdata->rowexact)) );
   BMSclearMemoryArray(consdata->roweventdata, SCIProwExactGetNNonz(consdata->rowexact));

   return SCIP_OKAY;
}

/** adds linear constraint as cut to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rowexact == NULL )
   {
      /* convert consdata object into LP row and exact lp row */
      SCIP_CALL( createRows(scip, cons) );
   }
   assert(consdata->rowlhs != NULL);
   assert(consdata->rowexact != NULL);

   if( consdata->nvars == 0 )
   {
      SCIPdebugMsg(scip, "Empty linear constraint enters LP: <%s>\n", SCIPconsGetName(cons));
   }

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->rowlhs) )
   {
      SCIPdebugMsg(scip, "adding relaxation of linear constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, consdata->rowlhs, NULL)) );
      SCIPdebug( SCIP_CALL( SCIPprintRowExact(scip, consdata->rowexact, NULL)) );

      /* if presolving is turned off, the row might be trivial */
      if( !SCIPrationalIsNegInfinity(consdata->lhs) || !SCIPrationalIsInfinity(consdata->rhs) )
      {
         SCIP_CALL( SCIPaddRow(scip, consdata->rowlhs, FALSE, cutoff) );
         SCIP_CALL( SCIPaddRowExact(scip, consdata->rowexact) );
      }
#ifndef NDEBUG
      else
      {
         int pr;
         int cr;
         SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrounds", &pr) );
         SCIP_CALL( SCIPgetIntParam(scip, "constraints/linear/maxprerounds", &cr) );
         assert( pr == 0 || cr == 0 );
      }
#endif
   }

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by given solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             separatecards,      /**< should knapsack cardinality cuts be generated? */
   SCIP_Bool             separateall,        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   int*                  ncuts,              /**< pointer to add up the number of found cuts */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   )
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int oldncuts;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   assert(ncuts != NULL);

   oldncuts = *ncuts;
   *cutoff = FALSE;

   SCIP_CALL( checkCons(scip, cons, conshdlrdata, sol, FALSE, (sol != NULL), &violated) );

   if( violated )
   {
      /* insert LP row as cut */
      SCIP_CALL( addRelaxation(scip, cons, cutoff) );
      (*ncuts)++;
   }

   if( *ncuts > oldncuts )
   {
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** propagation method for linear constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool             tightenbounds,      /**< should the variable's bounds be tightened? */
   SCIP_Bool             sortvars,           /**< should variable sorting for faster propagation be used? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Bool minactisrelax;
   SCIP_Bool maxactisrelax;
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /*SCIPdebugMsg(scip, "propagating linear constraint <%s>\n", SCIPconsGetName(cons));*/

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->eventdata == NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consCatchAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      assert(consdata->eventdata != NULL);
   }

   *cutoff = FALSE;

   /* we can only infer activity bounds of the linear constraint, if it is not modifiable */
   if( !SCIPconsIsModifiable(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      if( !SCIPconsIsInitial(cons) )
      {
         conshdlrdata->nconspropnoninit++;
         conshdlrdata->propnonzerosnoninit += consdata->nvars;
      }
      else
      {
         conshdlrdata->nconsprop++;
         conshdlrdata->propnonzeros += consdata->nvars;
      }

      /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
      if( !SCIPinRepropagation(scip) )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      /* tighten the variable's bounds */
      if( tightenbounds )
      {
         int oldnchgbds;

         oldnchgbds = *nchgbds;

         SCIP_CALL( tightenBounds(scip, cons, sortvars, cutoff, nchgbds) );

         if( *nchgbds > oldnchgbds )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
      }

      /* check constraint for infeasibility and redundancy */
      if( !(*cutoff) )
      {
         consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax,
            &isminsettoinfinity, &ismaxsettoinfinity);

         if( minactivity > consdata->rhsreal )
         {
            SCIPrationalDebugMessage("linear constraint <%s> is infeasible (rhs): activitybounds=[%.15g,%.15g], sides=[%q,%q]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, TRUE) );

            /**@todo analyze conflict detected in exactlinear constraint handler */
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( maxactivity < consdata->lhsreal )
         {
            SCIPrationalDebugMessage("linear constraint <%s> is infeasible (lhs): activitybounds=[%.15g,%.15g], sides=[%q,%q]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhsreal, consdata->rhsreal);
            SCIP_CALL( certificatePrintActivityConflict(scip, cons, consdata, FALSE) );

            /**@todo analyze conflict detected in exactlinear constraint handler */
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( minactivity >= consdata->lhsreal && maxactivity <= consdata->rhsreal )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhsreal, consdata->rhsreal);

            /* remove the constraint locally unless it has become empty, in which case it is removed globally */
            if( consdata->nvars > 0 )
               SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            else
               SCIP_CALL( SCIPdelCons(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Presolving methods
 */

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   SCIP_Bool checkexact;
   SCIP_Bool cutoff = FALSE;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( sol == NULL )
      checkexact =  SCIPlpExactIsSolved(scip);
   else
      checkexact = SCIPisExactSol(scip, sol);

   SCIPdebugMsg(scip, "Enforcement method of linear constraints for %s solution\n", sol == NULL ? "LP" : "relaxation");
   SCIPdebug( SCIPprintSol(scip, sol, NULL, FALSE));

   /* check for violated constraints
    * LP is processed at current node -> we can add violated linear constraints to the SCIP_LP
    */
   *result = SCIP_FEASIBLE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], conshdlrdata, sol, checkexact, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], &cutoff) );
         if( cutoff )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
      }
   }

   /* check all obsolete linear constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], conshdlrdata, sol, checkexact, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], &cutoff) );
         if( cutoff )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
      }
   }

   SCIPdebugMsg(scip, "-> constraints checked, %s\n", *result == SCIP_FEASIBLE ? "all constraints feasible" : "infeasibility detected");

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrExactLinear(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExactLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   /* check for event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   assert(nconss == 0 || conss != NULL);

   /* catch events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* catch all events */
      SCIP_CALL( consCatchAllEvents(scip, conss[c], conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExactLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   /* check for event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* drop events for the constraints */
   for( c = nconss - 1; c >= 0; --c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->eventdata != NULL )
      {
         /* drop all events */
         SCIP_CALL( consDropAllEvents(scip, conss[c], conshdlrdata->eventhdlr) );
         assert(consdata->eventdata == NULL);
      }
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreExactLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   /* delete all linear constraints that were upgraded to a more specific constraint type;
    * make sure, only active variables remain in the remaining constraints
    */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
      {
         SCIPerrorMessage("exact linear constraint upgrade not implemented yet\n");
         return SCIP_ERROR;
      }
      else
      {
         /* since we are not allowed to detect infeasibility in the exitpre stage, we dont give an infeasible pointer */
         SCIP_CALL( applyFixings(scip, conss[c], NULL) );
      }
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolExactLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->rowlhs != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowlhs) );
         SCIPfreeBlockMemoryArray(scip, &consdata->roweventdata, SCIProwExactGetNNonz(consdata->rowexact));

         SCIP_CALL( SCIPreleaseRowExact(scip, &consdata->rowexact) );
         if( consdata->rowrhs != NULL )
         {
            assert(!consdata->onerowrelax);
            SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowrhs) );
         }
      }
   }

   /**@todo when enabling restarts, extend SCIPconvertCutsToConss() in order to convert exact cuts to exactlinear
    *       constraints and call here
    */

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(SCIPisExact(scip));
   assert(cons != NULL);

   if( SCIPconsIsDeleted(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSDATA* consdata;

      assert(conshdlr != NULL);
      assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

      /* get constraint data */
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* check for event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* free event data */
      if( consdata->eventdata != NULL )
      {
         /* drop bound change events of variables */
         SCIP_CALL( consDropAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      }
      assert(consdata->eventdata == NULL);
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(SCIPisExact(scip));
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   if( (*consdata)->eventdata != NULL )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* drop all events */
      SCIP_CALL( consDropAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      assert((*consdata)->eventdata == NULL);
   }
   /* free linear constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(scip != NULL);
   assert(SCIPisExact(scip));
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->rowlhs == NULL && sourcedata->rowexact == NULL);  /* in original problem, there cannot be LP rows */

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars, sourcedata->vals, sourcedata->lhs, sourcedata->rhs) );

   if( sourcedata->nvars > 0 )
      consdataScaleMinValue(scip, targetdata, 2 * SCIPepsilon(scip));

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpExactLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   *infeasible = FALSE;

   for( c = 0; c < nconss && !(*infeasible); ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      /* add both the relaxation to the fp-lp as well as the correct constraint to the exact lp */
      SCIP_CALL( addRelaxation(scip, conss[c], infeasible) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real cutoffbound;
   SCIP_Real maxbound;
   SCIP_Bool separatecards;
   SCIP_Bool cutoff;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   /* check if we want to produce knapsack cardinality cuts at this node */
   loclowerbound = SCIPgetLocalLowerbound(scip);
   glblowerbound = SCIPgetLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);
   maxbound = glblowerbound + SCIPrationalGetReal(conshdlrdata->maxcardbounddist) * (cutoffbound - glblowerbound);

   separatecards = SCIPisLE(scip, loclowerbound, maxbound);
   separatecards = separatecards && (SCIPgetNLPBranchCands(scip) > 0);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts && !cutoff; ++c )
   {
      SCIPdebugMsg(scip, "separating exact linear constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIP_CALL( separateCons(scip, conss[c], conshdlrdata, NULL, separatecards, conshdlrdata->separateall, &ncuts, &cutoff) );
   }

   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts && !cutoff; ++c )
   {
      SCIPdebugMsg(scip, "separating exact linear constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIP_CALL( separateCons(scip, conss[c], conshdlrdata, sol, TRUE, conshdlrdata->separateall, &ncuts, &cutoff) );
   }

   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMsg(scip, "Enfops method of linear constraints\n");

   if( !SCIPisExact(scip) )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the solution is infeasible anyway due to objective value, skip the enforcement */
   if( objinfeasible )
   {
      SCIPdebugMsg(scip, "-> pseudo solution is objective infeasible, return.\n");

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], conshdlrdata, NULL, FALSE, TRUE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   SCIPdebugMsg(scip, "-> constraints checked, %s\n", *result == SCIP_FEASIBLE ? "all constraints feasible" : "infeasibility detected");

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool checkexact;
   int c;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if the fp-solution has a stand-in exact solution we check that instead */
   checkexact = SCIPisExactSol(scip, sol);

   /* check all linear constraints for feasibility */
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      SCIP_Bool violated = FALSE;
      SCIP_CALL( checkCons(scip, conss[c], conshdlrdata, sol, checkexact, checklprows, &violated) );

      if( violated )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIP_CONSDATA* consdata;
            SCIP_RATIONAL* activity;

            SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &activity) );

            consdata = SCIPconsGetData(conss[c]);
            assert( consdata != NULL);

            consdataGetActivity(scip, consdata, sol, checkexact, activity);

            SCIP_CALL( consPrintConsSol(scip, conss[c], sol, checkexact, NULL ) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( SCIPrationalIsAbsInfinity(activity) ) /*lint !e777*/
               SCIPinfoMessage(scip, NULL, "activity invalid due to positive and negative infinity contributions\n");
            else if( SCIPrationalIsLT(activity, consdata->lhs) )
            {
               SCIPrationalDiff(activity, consdata->lhs, activity);
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by ");
               SCIPrationalPrint(activity);
               SCIPinfoMessage(scip, NULL, "\n");
            }
            else if( SCIPrationalIsGT(activity, consdata->rhs) )
            {
               SCIPrationalDiff(activity, activity, consdata->rhs);
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by ");
               SCIPrationalPrint(activity);
               SCIPinfoMessage(scip, NULL, "\n");
            }

            SCIPrationalFreeBuffer(SCIPbuffer(scip), &activity);
         }
      }
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool tightenbounds;
   SCIP_Bool cutoff;

   int nchgbds;
   int i;

   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check, if we want to tighten variable's bounds (in probing, we always want to tighten the bounds) */
   if( SCIPinProbing(scip) )
      tightenbounds = TRUE;
   else
   {
      int depth;
      int propfreq;
      int tightenboundsfreq;

      depth = SCIPgetDepth(scip);
      propfreq = SCIPconshdlrGetPropFreq(conshdlr);
      tightenboundsfreq = propfreq * conshdlrdata->tightenboundsfreq;
      tightenbounds = (conshdlrdata->tightenboundsfreq >= 0)
         && ((tightenboundsfreq == 0 && depth == 0) || (tightenboundsfreq >= 1 && (depth % tightenboundsfreq == 0)));
   }

   cutoff = FALSE;
   nchgbds = 0;

   /* process constraints marked for propagation */
   for( i = 0; i < nmarkedconss && !cutoff; i++ )
   {
      SCIP_CALL( SCIPunmarkConsPropagate(scip, conss[i]) );
      SCIP_CALL( propagateCons(scip, conss[i], tightenbounds,
            conshdlrdata->sortvars, &cutoff, &nchgbds) );
   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int i;

   assert(scip != NULL);
   assert(SCIPisExact(scip));
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslhs = !SCIPrationalIsNegInfinity(consdata->lhs);
   hasrhs = !SCIPrationalIsInfinity(consdata->rhs);

   /* update rounding locks of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPrationalIsPositive(consdata->vals[i]) )
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[i], locktype, nlockspos, nlocksneg) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[i], locktype, nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[i], locktype, nlocksneg, nlockspos) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[i], locktype, nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsExactLinear)
{
   assert(scip != NULL);
   assert(SCIPisExact(scip) || nconss == 0);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   if( nconss > 0 )
   {
      SCIP_CALL( performVarDeletions(scip, conshdlr, conss, nconss) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExactLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyExactLinear)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_INTERVAL* sourcecoefs;
   const char* consname;
   int nvars;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsExactLinear(sourcescip, sourcecons);
   sourcecoefs = SCIPgetValsRealExactLinear(sourcescip, sourcecons);
   nvars = SCIPgetNVarsExactLinear(sourcescip, sourcecons);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIP_CALL( SCIPcopyConsExactLinear(scip, cons, sourcescip, consname, nvars, sourcevars, sourcecoefs,
         SCIPrationalGetReal(SCIPgetLhsExactLinear(sourcescip, sourcecons)), SCIPrationalGetReal(SCIPgetRhsExactLinear(sourcescip, sourcecons)), varmap, consmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL || *valid == FALSE);

   return SCIP_OKAY;
}

/* find operators '<=', '==', '>=', [free] in input string and return those places. There should only be one operator,
 * except for ranged rows for which exactly two operators '<=' must be present
 */
static
SCIP_RETCODE findOperators(
   const char*           str,                /**< null terminated input string */
   char**                firstoperator,      /**< pointer to store the string starting at the first operator */
   char**                secondoperator,     /**< pointer to store the string starting at the second operator */
   SCIP_Bool*            success             /**< pointer to store if the line contains a valid operator order */
   )
{
   char* curr;

   assert(str != NULL);
   assert(firstoperator != NULL);
   assert(secondoperator != NULL);

   *firstoperator = NULL;
   *secondoperator = NULL;

   curr = (char*)str;
   *success = TRUE;

   /* loop over the input string to find all operators */
   while( *curr && *success )
   {
      SCIP_Bool found = FALSE;
      int increment = 1;

      /* try if we found a possible operator */
      switch( *curr )
      {
      case '<':
      case '=':
      case '>':

         /* check if the two characters curr[0,1] form an operator together */
         if( curr[1] == '=' )
         {
            found = TRUE;

            /* update increment to continue after this operator */
            increment = 2;
         }
         break;
      case '[':
         if( strncmp(curr, "[free]", 6) == 0 )
         {
            found = TRUE;

            /* update increment to continue after this operator */
            increment = 6;
         }
         break;
      default:
         break;
      }

      /* assign the found operator to the first or second pointer and check for violations of the linear constraint grammar */
      if( found )
      {
         if( *firstoperator == NULL )
         {
            *firstoperator = curr;
         }
         else
         {
            if( *secondoperator != NULL )
            {
               SCIPerrorMessage("Found more than two operators in line %s\n", str);
               *success = FALSE;
            }
            else if( strncmp(*firstoperator, "<=", 2) != 0 )
            {
               SCIPerrorMessage("Two operators in line that is not a ranged row: %s", str);
               *success = FALSE;
            }
            else if( strncmp(curr, "<=", 2) != 0 )
            {
               SCIPerrorMessage("Bad second operator, expected ranged row specification: %s", str);
               *success = FALSE;
            }

            *secondoperator = curr;
         }
      }

      curr += increment;
   }

   /* check if we did find at least one operator */
   if( *success )
   {
      if( *firstoperator == NULL )
      {
         SCIPerrorMessage("Could not find any operator in line %s\n", str);
         *success = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseExactLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_RATIONAL** coefs;
   int        nvars;
   int        coefssize;
   int        requsize;
   SCIP_RATIONAL*  lhs;
   SCIP_RATIONAL*  rhs;
   char*      endptr;
   char*      firstop;
   char*      secondop;
   SCIP_Bool  operatorsuccess;
   char*      lhsstrptr;
   char*      rhsstrptr;
   char*      varstrptr;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* set left and right hand side to their default values */
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &lhs) );
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &rhs) );

   SCIPrationalSetNegInfinity(lhs);
   SCIPrationalSetInfinity(rhs);

   (*success) = FALSE;

   /* return of string empty */
   if( !*str )
      return SCIP_OKAY;

   /* ignore whitespace */
   while( isspace((unsigned char)*str) )
      ++str;

   /* find operators in the line first, all other remaining parsing depends on occurence of the operators '<=', '>=', '==',
    * and the special word [free]
    */
   SCIP_CALL( findOperators(str, &firstop, &secondop, &operatorsuccess) );

   /* if the grammar is not valid for parsing a linear constraint, return */
   if( ! operatorsuccess )
      return SCIP_OKAY;

   varstrptr = (char *)str;
   lhsstrptr = rhsstrptr = NULL;
   assert(firstop != NULL);

   /* assign the strings for parsing the left hand side, right hand side, and the linear variable sum */
   switch( *firstop )
   {
      case '<':
         assert(firstop[1] == '=');
         /* we have ranged row lhs <= a_1 x_1 + ... + a_n x_n <= rhs */
         if( secondop != NULL )
         {
            assert(secondop[0] == '<' && secondop[1] == '=');
            lhsstrptr = (char *)str;
            varstrptr = firstop + 2;
            rhsstrptr = secondop + 2;
         }
         else
         {
            /* we have an inequality with infinite left hand side a_1 x_1 + ... + a_n x_n <= rhs */
            lhsstrptr = NULL;
            varstrptr = (char *)str;
            rhsstrptr = firstop + 2;
         }
         break;
      case '>':
         assert(firstop[1] == '=');
         assert(secondop == NULL);
         /* we have a_1 x_1 + ... + a_n x_n >= lhs */
         lhsstrptr = firstop + 2;
         break;
      case '=':
         assert(firstop[1] == '=');
         assert(secondop == NULL);
         /* we have a_1 x_1 + ... + a_n x_n == lhs (rhs) */
         rhsstrptr = firstop + 2;
         lhsstrptr = firstop + 2;
         break;
      case '[':
         assert(strncmp(firstop, "[free]", 6) == 0);
         assert(secondop == NULL);
         /* nothing to assign in case of a free a_1 x_1 + ... + a_n x_n [free] */
         break;
      default:
         /* it should not be possible that a different character appears in that position */
         SCIPerrorMessage("Parsing has wrong operator character '%c', should be one of <=>[", *firstop);
         return SCIP_READERROR;
   }

   /* parse left hand side, if necessary */
   if( lhsstrptr != NULL )
   {
      if( ! SCIPparseRational(scip, lhsstrptr, lhs, &endptr) )
      {
         SCIPerrorMessage("error parsing left hand side number from <%s>\n", lhsstrptr);
         return SCIP_OKAY;
      }

      /* in case of an equation, assign the left also to the right hand side */
      if( rhsstrptr == lhsstrptr )
         SCIPrationalSetRational(rhs, lhs);
   }

   /* parse right hand side, if different from left hand side */
   if( rhsstrptr != NULL && rhsstrptr != lhsstrptr )
   {
      if( ! SCIPparseRational(scip, rhsstrptr, rhs, &endptr) )
      {
         SCIPerrorMessage("error parsing right hand side number from <%s>\n", lhsstrptr);
         return SCIP_OKAY;
      }
   }

   /* initialize buffers for storing the variables and coefficients */
   coefssize = 100;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  coefssize) );
   SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &coefs, coefssize) );

   assert(varstrptr != NULL);

   /* parse linear sum to get variables and coefficients */
   SCIP_CALL( SCIPparseVarsLinearsumExact(scip, varstrptr, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );

   if( *success && requsize > coefssize )
   {
      /* realloc buffers and try again */
      SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  requsize) );
      SCIP_CALL( SCIPrationalReallocBufferArray(SCIPbuffer(scip), &coefs, coefssize, requsize) );

      coefssize = requsize;

      SCIP_CALL( SCIPparseVarsLinearsumExact(scip, varstrptr, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );
      assert(!*success || requsize <= coefssize); /* if successful, then should have had enough space now */
   }

   if( !*success )
   {
      SCIPerrorMessage("no luck in parsing linear sum '%s'\n", varstrptr);
   }
   else
   {
      SCIP_CALL( SCIPcreateConsExactLinear(scip, cons, name, nvars, vars, coefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPrationalFreeBufferArray(SCIPbuffer(scip), &coefs, coefssize);
   SCIPfreeBufferArray(scip, &vars);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &rhs);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &lhs);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/**! [Callback for the number of variables]*/
/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExactLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}
/**! [Callback for the number of variables]*/

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecExactLinear)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_EVENTTYPE eventtype;
   SCIP_Bool updateActivities;
   assert(scip != NULL);
   assert(SCIPisExact(scip));
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   cons = eventdata->cons;
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   if( consdata == NULL )
      return SCIP_OKAY;
   /* we can skip events dropped for deleted constraints */
   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);
   updateActivities = ((consdata->rowexact != NULL) == eventdata->rowvar) && consdata->validactivities;
   assert(!consdata->validactivities || (consdata->validminact && consdata->validmaxact && consdata->validglbminact && consdata->validglbmaxact));

   if( ((eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0) )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_INTERVAL valrange;
      int varpos;
      varpos = eventdata->varpos;

      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      valrange = consdata->valsreal[varpos];

      /* we only need to update the activities if the constraint is active,
       * otherwise we mark them to be invalid
       */
      if( SCIPconsIsActive(cons) )
      {
         /* update the activity values */
         if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
            consdataUpdateActivitiesLb(scip, consdata, var, oldbound, newbound, valrange);
         else
         {
            assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
            consdataUpdateActivitiesUb(scip, consdata, var, oldbound, newbound, valrange);
         }
      }
      else
         consdataInvalidateActivities(consdata);

      consdata->presolved = FALSE;
      consdata->rangedrowpropagated = 0;

      /* bound change can turn the constraint infeasible or redundant only if it was a tightening */
      if( (eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0 )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

         /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
         if( consdata->maxactdeltavar == var )
         {
            consdata->maxactdelta = SCIP_INVALID;
            consdata->maxactdeltavar = NULL;
         }

         /* check whether bound tightening might now be successful */
         if( consdata->boundstightened > 0)
         {
            switch( eventtype )
            {
            case SCIP_EVENTTYPE_LBTIGHTENED:
               if( (valrange.sup > 0.0 ? !SCIPisInfinity(scip, consdata->rhsreal) : !SCIPisInfinity(scip, -consdata->lhsreal)) )
                  consdata->boundstightened = 0;
               break;
            case SCIP_EVENTTYPE_UBTIGHTENED:
               if( (valrange.sup > 0.0 ? !SCIPisInfinity(scip, -consdata->lhsreal) : !SCIPisInfinity(scip, consdata->rhsreal)) )
                  consdata->boundstightened = 0;
               break;
            default:
               SCIPerrorMessage("invalid event type %" SCIP_EVENTTYPE_FORMAT "\n", eventtype);
               return SCIP_INVALIDDATA;
            }
         }
      }
      /* update maximal activity delta if a bound was relaxed */
      else if( !SCIPisInfinity(scip, consdata->maxactdelta) )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real domain;
         SCIP_Real delta;

         assert((eventtype & SCIP_EVENTTYPE_BOUNDRELAXED) != 0);

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         domain = ub - lb;
         delta = SCIPintervalAbsMax(valrange) * domain;

         if( delta > consdata->maxactdelta )
         {
            consdata->maxactdelta = delta;
            consdata->maxactdeltavar = var;
         }
      }
   }
   else if( (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0 )
   {
      /* we want to remove the fixed variable */
      consdata->presolved = FALSE;
      consdata->removedfixings = FALSE;
      consdata->rangedrowpropagated = 0;

      /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
      if( consdata->maxactdeltavar == var )
      {
         consdata->maxactdelta = SCIP_INVALID;
         consdata->maxactdeltavar = NULL;
      }
   }
   else if( (eventtype & SCIP_EVENTTYPE_VARUNLOCKED) != 0 )
   {
      /* there is only one lock left: we may multi-aggregate the variable as slack of an equation */
      assert(SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) <= 1);
      assert(SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) <= 1);
      consdata->presolved = FALSE;
   }
   else if( (eventtype & SCIP_EVENTTYPE_GBDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_INTERVAL valrange;
      int varpos;

      varpos = eventdata->varpos;

      if( updateActivities )
      {
         oldbound = SCIPeventGetOldbound(event);
         newbound = SCIPeventGetNewbound(event);
         assert(var != NULL);
         assert(consdata->vars[varpos] == var);
         valrange = consdata->valsreal[varpos];

         consdata->rangedrowpropagated = 0;

         /* update the activity values */
         if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
            consdataUpdateActivitiesGlbLb(scip, consdata, oldbound, newbound, valrange);
         else
         {
            assert((eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0);
            consdataUpdateActivitiesGlbUb(scip, consdata, oldbound, newbound, valrange);
         }
      }

      /* if the variable is binary but not fixed it had to become binary due to this global change */
      if( SCIPvarIsBinary(var) && SCIPisGT(scip, SCIPvarGetUbGlobal(var), SCIPvarGetLbGlobal(var)) )
      {
         if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
            consdata->indexsorted = FALSE;
         else
            consdata->coefsorted = FALSE;
      }
   }
   else if( ((eventtype & SCIP_EVENTTYPE_TYPECHANGED) != 0) )
   {
      assert(SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED);

      /* for presolving it only matters if a variable type changed from continuous to some kind of integer */
      consdata->presolved = (consdata->presolved && SCIPeventGetOldtype(event) < SCIP_VARTYPE_CONTINUOUS);

      /* the ordering is preserved if the type changes from something different to binary to binary but SCIPvarIsBinary() is true */
      consdata->indexsorted = (consdata->indexsorted && SCIPeventGetNewtype(event) == SCIP_VARTYPE_BINARY && SCIPvarIsBinary(var));
   }
   else if( (eventtype & SCIP_EVENTTYPE_VARDELETED) )
   {
      consdata->varsdeleted = TRUE;
   }
   return SCIP_OKAY;
}


/*
 * Callback methods of conflict handler
 */

/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExactLinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecExactLinear, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExactLinear, consEnfopsExactLinear, consCheckExactLinear, consLockExactLinear,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* mark constraint handler as exact */
   SCIPconshdlrMarkExact(conshdlr);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyExactLinear, consCopyExactLinear) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveExactLinear) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteExactLinear) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsExactLinear) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitExactLinear) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreExactLinear) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolExactLinear) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeExactLinear) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsExactLinear) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsExactLinear) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitExactLinear) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpExactLinear) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseExactLinear) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintExactLinear) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropExactLinear, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpExactLinear, consSepasolExactLinear, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransExactLinear) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxExactLinear) );

   /* add constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/tightenboundsfreq",
         "multiplier on propagation frequency, how often the bounds are tightened (-1: never, 0: only at root)",
         &conshdlrdata->tightenboundsfreq, TRUE, DEFAULT_TIGHTENBOUNDSFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/separateall",
         "should all constraints be subject to cardinality cut generation instead of only the ones with non-zero dual value?",
         &conshdlrdata->separateall, FALSE, DEFAULT_SEPARATEALL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/sortvars", "apply binaries sorting in decr. order of coeff abs value?",
         &conshdlrdata->sortvars, TRUE, DEFAULT_SORTVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/propcont",
         "should bounds on continuous variables be tightened by propagation?",
         &conshdlrdata->propcont, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/limitdenom",
         "should denominators of rational bounds on continuous variables be controlled?",
         &conshdlrdata->limitdenom, TRUE, DEFAULT_LIMITDENOM, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/" CONSHDLR_NAME "/boundmaxdenom",
         "maximal denominator for rational bounds on continuous variables after propagation",
         &conshdlrdata->boundmaxdenom, TRUE, DEFAULT_BOUNDMAXDENOM, 1L, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of constraint entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of constraint */
   SCIP_RATIONAL*        rhs,                /**< right hand side of constraint */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE && nvars > 0 )
   {
      SCIP_VAR** consvars;
      SCIP_RATIONAL** consvals;
      SCIP_RATIONAL* constant;
      int nconsvars;
      int requiredsize;

      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &constant) );

      nconsvars = nvars;
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvars, vars, nconsvars) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvals, vals, nconsvars) );

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSumExact(scip, consvars, consvals, &nconsvars, nconsvars, constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSumExact(scip, consvars, consvals, &nconsvars, requiredsize, constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      /* adjust sides and check that we do not subtract infinity values */
      if( SCIPrationalIsAbsInfinity(constant) )
      {
         if( SCIPrationalIsNegative(constant) )
         {
            if( SCIPrationalIsInfinity(lhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            if( SCIPrationalIsInfinity(rhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }

            SCIPrationalSetNegInfinity(lhs);
            SCIPrationalSetNegInfinity(rhs);
         }
         else
         {
            if( SCIPrationalIsNegInfinity(lhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            if( SCIPrationalIsNegInfinity(rhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }

            SCIPrationalSetInfinity(lhs);
            SCIPrationalSetInfinity(lhs);
         }
      }

      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, nconsvars, consvars, consvals, lhs, rhs) );
      assert(consdata != NULL);

      SCIPrationalFreeBuffer(SCIPbuffer(scip), &constant);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, vals, lhs, rhs) );
      assert(consdata != NULL);
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a linear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of constraint entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of constraint */
   SCIP_RATIONAL*        rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsExactLinear(scip, cons, name, nvars, vars, vals, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates a linear constraint from an exact linear constraint by rounding values to floating-point and captures it */
SCIP_RETCODE SCIPcopyConsExactLinear(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in source variable array */
   SCIP_VAR**            sourcevars,         /**< source variables of the linear constraints */
   SCIP_INTERVAL*        sourcecoefs,        /**< coefficient array of the linear constraint, or NULL if all coefficients are one */
   SCIP_Real             lhs,                /**< left hand side of the linear constraint */
   SCIP_Real             rhs,                /**< right hand side of the linear constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* coefs;

   SCIP_Real constant;
   int requiredsize;
   int v;
   SCIP_Bool success;

   /**@todo This method is currently only used for subSCIPs in floating-point heuristics, but should be extended to be
    *       able to perform an exact copy in the future.  This would allow application of the cons_components presolver,
    *       for example.  In this case, whether an exact or an fp copy is created, could probably be decided by checking
    *       SCIPisExact() for the target SCIP.
    */
   assert(!SCIPisExact(scip));
   (*valid) = FALSE;

   if( SCIPisGT(scip, lhs, rhs) )
   {
      return SCIP_OKAY;
   }

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, 0, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      return SCIP_OKAY;
   }

   /* duplicate variable array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, sourcevars, nvars) );

   /* duplicate coefficient array */
   if( sourcecoefs != NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
      for( int i = 0; i < nvars; i++ )
      {
         coefs[i] = SCIPintervalGetSup(sourcecoefs[i]);
         assert(!SCIPisInfinity(scip, coefs[i]) && !SCIPisInfinity(scip, -coefs[i]));
      }
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
      for( v = 0; v < nvars; ++v )
         coefs[v] = 1.0;
   }

   constant = 0.0;

   /* transform source variable to active variables of the source SCIP since only these can be mapped to variables of
    * the target SCIP
    */
   if( !SCIPvarIsOriginal(vars[0]) )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, nvars, &constant, &requiredsize) );

      if( requiredsize > nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, requiredsize, &constant, &requiredsize) );
         assert(requiredsize <= nvars);
      }
   }
   else
   {
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsOriginal(vars[v]));
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &coefs[v], &constant) );
         assert(vars[v] != NULL);
      }
   }

   success = TRUE;
   /* map variables of the source constraint to variables of the target SCIP */
   for( v = 0; v < nvars && success; ++v )
   {
      SCIP_VAR* var;
      var = vars[v];

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &vars[v], varmap, consmap, global, &success) );
      assert(!(success) || vars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( success )
   {
      if( !SCIPisInfinity(scip, -lhs) )
         lhs -= constant;

      if( !SCIPisInfinity(scip, rhs) )
         rhs -= constant;

      SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, coefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** adds coefficient to linear constraint (if it is not zero) */
SCIP_RETCODE SCIPaddCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      return SCIP_INVALIDDATA;
   }

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_RATIONAL** consvals;
      SCIP_RATIONAL* constant;
      SCIP_RATIONAL* rhs;
      SCIP_RATIONAL* lhs;
      int nconsvars;
      int requiredsize;
      int v;

      SCIPerrorMessage("adding coefficients after presolving not supported yet in exact solving mode \n");
      SCIPABORT();

      nconsvars = 1;
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
      consvars[0] = var;
      SCIP_CALL( SCIPrationalCopyBlock(SCIPblkmem(scip), &consvals[0], val) );
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &constant) );

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSumExact(scip, consvars, consvals, &nconsvars, nconsvars, constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSumExact(scip, consvars, consvals, &nconsvars, requiredsize, constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      SCIP_CALL( SCIPrationalCopyBlock(SCIPblkmem(scip), &lhs, consdata->lhs) );
      SCIP_CALL( SCIPrationalCopyBlock(SCIPblkmem(scip), &rhs, consdata->rhs) );

      /* adjust sides and check that we do not subtract infinity values */
      if( SCIPrationalIsAbsInfinity(constant) )
      {
         SCIPfreeBufferArray(scip, &consvals);
         SCIPfreeBufferArray(scip, &consvars);

         SCIPerrorMessage("adding variable <%s> to constraint <%s> leads to infinite constant and cannot be handled safely\n",
            SCIPvarGetName(var), SCIPconsGetName(cons));

         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
      /* constant is not infinite */
      else
      {
         if( !SCIPrationalIsAbsInfinity(lhs) )
            SCIPrationalDiff(lhs, lhs, constant);
         if( !SCIPrationalIsAbsInfinity(rhs) )
            SCIPrationalDiff(rhs, rhs, constant);
      }

      /* add all active variables to constraint */
      for( v = nconsvars - 1; v >= 0; --v )
      {
         SCIP_CALL( addCoef(scip, cons, consvars[v], consvals[v]) );
      }

      /* update left and right hand sides */
      SCIP_CALL( chgLhs(scip, cons, lhs));
      SCIP_CALL( chgRhs(scip, cons, rhs));

      SCIPrationalFreeBuffer(SCIPbuffer(scip), &constant);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      SCIP_CALL( addCoef(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** changes coefficient of variable in linear constraint; deletes the variable if coefficient is zero; adds variable if
 *  not yet contained in the constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_RETCODE SCIPchgCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< new coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool found;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) || !SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints and variables\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   found = FALSE;
   i = 0;
   while( i < consdata->nvars )
   {
      if( vars[i] == var )
      {
         if( found || SCIPrationalIsZero(val) )
         {
            SCIP_CALL( delCoefPos(scip, cons, i) );

            /* decrease i by one since otherwise we would skip the coefficient which has been switched to position i */
            i--;
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, i, val) );
         }
         found = TRUE;
      }
      i++;
   }

   if( !found && !SCIPrationalIsZero(val) )
   {
      SCIP_CALL( SCIPaddCoefExactLinear(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** deletes variable from linear constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_RETCODE SCIPdelCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   SCIP_RATIONAL* temp;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &temp) );

   SCIP_CALL( SCIPchgCoefExactLinear(scip, cons, var, temp) );

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &temp);

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
SCIP_RATIONAL* SCIPgetLhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of linear constraint */
SCIP_RATIONAL* SCIPgetRhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** changes left hand side of linear constraint */
SCIP_RETCODE SCIPchgLhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        lhs                 /**< new left hand side */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
SCIP_RETCODE SCIPchgRhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

/** gets the number of variables in the linear constraint */
int SCIPgetNVarsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the linear constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_INTERVAL* SCIPgetValsRealExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->valsreal;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_RATIONAL** SCIPgetValsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vals;
}

/** gets the activity of the linear constraint in the given solution
 *
 *  @note if the activity comprises positive and negative infinity contributions, the result is currently undefined
 */
SCIP_RETCODE SCIPgetActivityExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_RATIONAL*        ret
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      SCIPrationalSetInfinity(ret);  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rowexact != NULL )
   {
      SCIP_CALL( SCIPgetRowSolActivityExact(scip, consdata->rowexact, sol, FALSE, ret) );
   }
   else
      consdataGetActivity(scip, consdata, sol, TRUE, ret);

   return SCIP_OKAY;
}

/** gets the feasibility of the linear constraint in the given solution */
SCIP_RETCODE SCIPgetFeasibilityExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_RATIONAL*        ret                 /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rowexact != NULL )
      SCIP_CALL( SCIPgetRowSolFeasibilityExact(scip, consdata->rowexact, sol, ret) );
   else
      consdataGetFeasibility(scip, consdata, sol, ret);

   return SCIP_OKAY;
}

/** gets the dual solution of the linear constraint in the current LP
 *
 *  @note this method currently returns the value from the floating-point LP
 */
void SCIPgetFpDualsolExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        ret                 /**< result pointer */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(!SCIPconsIsOriginal(cons)); /* original constraints would always return 0 */

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rowlhs != NULL )
      SCIPrationalSetReal(ret, SCIProwGetDualsol(consdata->rowlhs));
   else
      SCIPrationalSetReal(ret, 0.0);
}

/** gets the dual Farkas value of the linear constraint in the current infeasible LP
 *
 *  @note this method currently returns an approximate value from the floating-point LP
 */
void SCIPgetFpDualfarkasExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        ret                 /**< result pointer */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(!SCIPconsIsOriginal(cons)); /* original constraints would always return 0 */

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rowlhs != NULL )
      SCIPrationalSetReal(ret, SCIProwGetDualfarkas(consdata->rowlhs));
   else
      SCIPrationalSetReal(ret, 0.0);
}

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rowlhs;
}

/** returns the exact linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROWEXACT* SCIPgetRowExactExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not of type exactlinear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rowexact;
}
