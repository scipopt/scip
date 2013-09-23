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

/**@file   cons_quadratic.c
 * @brief  constraint handler for quadratic constraints \f$\textrm{lhs} \leq \sum_{i,j=1}^n a_{i,j} x_ix_j + \sum_{i=1}^n b_i x_i \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 * 
 * @todo SCIP might fix linear variables on +/- infinity; remove them in presolve and take care later
 * @todo round constraint sides to integers if all coefficients and variables are (impl.) integer
 * @todo constraints in one variable should be replaced by linear or bounddisjunction constraint
 * @todo check if some quadratic terms appear in several constraints and try to simplify (e.g., nous1)
 * @todo skip separation in enfolp if for current LP (check LP id) was already separated
 * @todo watch unbounded variables to enable/disable propagation
 * @todo sort order in bilinvar1/bilinvar2 such that the var which is involved in more terms is in bilinvar1, and use this info propagate and AddLinearReform
 * @todo catch/drop events in consEnable/consDisable, do initsol/exitsol stuff also when a constraint is enabled/disabled during solve
 * @todo underestimate for multivariate concave quadratic terms as in cons_nonlinear
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h> /* for strcmp */ 
#include <ctype.h>  /* for isspace */
#include <math.h>

#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_and.h"
#include "scip/cons_varbound.h"
#include "scip/intervalarith.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpi_ipopt.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "quadratic"
#define CONSHDLR_DESC          "quadratic constraints of the form lhs <= b' x + x' A x <= rhs"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -50 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING SCIP_PROPTIMING_BEFORELP

#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#define NONLINCONSUPGD_PRIORITY   40000 /**< priority of upgrading nonlinear constraints */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/* Activating this define enables reformulation of bilinear terms x*y with implications from x to y into linear terms.
 * However, implications are not enforced by SCIP. Thus, if, e.g., the used implication was derived from this constraint and we then reformulate the constraint,
 * then the implication may not be enforced in a solution.
 * This issue need to be fixed before this feature can be enabled.
 */
/* #define CHECKIMPLINBILINEAR */

/* enable new propagation for bivariate quadratic terms */
#define PROPBILINNEW

/*
 * Data structures
 */

/** eventdata for variable bound change events in quadratic constraints */
struct SCIP_QuadVarEventData
{
   SCIP_CONSDATA*        consdata;           /**< the constraint data */
   int                   varidx;             /**< the index of the variable which bound change is caught, positive for linear variables, negative for quadratic variables */
   int                   filterpos;          /**< position of eventdata in SCIP's event filter */
};

/** Data of a quadratic constraint. */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   int                   nlinvars;           /**< number of linear variables */
   int                   linvarssize;        /**< length of linear variable arrays */
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */
   SCIP_QUADVAREVENTDATA** lineventdata;       /**< eventdata for bound change of linear variable */

   int                   nquadvars;          /**< number of variables in quadratic terms */
   int                   quadvarssize;       /**< length of quadratic variable terms arrays */
   SCIP_QUADVARTERM*     quadvarterms;       /**< array with quadratic variable terms */

   int                   nbilinterms;        /**< number of bilinear terms */
   int                   bilintermssize;     /**< length of bilinear term arrays */
   SCIP_BILINTERM*       bilinterms;         /**< bilinear terms array */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   unsigned int          linvarssorted:1;    /**< are the linear variables already sorted? */
   unsigned int          linvarsmerged:1;    /**< are equal linear variables already merged? */
   unsigned int          quadvarssorted:1;   /**< are the quadratic variables already sorted? */
   unsigned int          quadvarsmerged:1;   /**< are equal quadratic variables already merged? */
   unsigned int          bilinsorted:1;      /**< are the bilinear terms already sorted? */
   unsigned int          bilinmerged:1;      /**< are equal bilinear terms already merged? */

   unsigned int          isconvex:1;         /**< is quadratic function is convex ? */
   unsigned int          isconcave:1;        /**< is quadratic function is concave ? */
   unsigned int          iscurvchecked:1;    /**< is quadratic function checked on convexity or concavity ? */
   unsigned int          isremovedfixings:1; /**< did we removed fixed/aggr/multiaggr variables ? */
   unsigned int          ispropagated:1;     /**< was the constraint propagated with respect to the current bounds ? */
   unsigned int          ispresolved:1;      /**< did we checked for possibilities of upgrading or implicit integer variables ? */
   unsigned int          initialmerge:1;     /**< did we perform an initial merge and clean in presolving yet ? */
#ifdef CHECKIMPLINBILINEAR
   unsigned int          isimpladded:1;      /**< has there been an implication added for a binary variable in a bilinear term? */
#endif

   SCIP_Real             minlinactivity;     /**< sum of minimal activities of all linear terms with finite minimal activity */
   SCIP_Real             maxlinactivity;     /**< sum of maximal activities of all linear terms with finite maximal activity */
   int                   minlinactivityinf;  /**< number of linear terms with infinite minimal activity */
   int                   maxlinactivityinf;  /**< number of linear terms with infinity maximal activity */
   SCIP_INTERVAL         quadactivitybounds; /**< bounds on the activity of the quadratic term, if up to date, otherwise empty interval */
   SCIP_Real             activity;           /**< activity of quadratic function w.r.t. current solution */
   SCIP_Real             lhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   int                   linvar_maydecrease; /**< index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
   int                   linvar_mayincrease; /**< index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */

   SCIP_VAR**            sepaquadvars;       /**< variables corresponding to quadvarterms to use in separation, only available in solving stage */
   int*                  sepabilinvar2pos;   /**< position of second variable in bilinear terms to use in separation, only available in solving stage */
   SCIP_Real             lincoefsmin;        /**< minimal absolute value of coefficients in linear part, only available in solving stage */
   SCIP_Real             lincoefsmax;        /**< maximal absolute value of coefficients in linear part, only available in solving stage */

   SCIP_Real*            factorleft;         /**< coefficients of left factor if constraint function is factorable */
   SCIP_Real*            factorright;        /**< coefficients of right factor if constraint function is factorable */
};

/** quadratic constraint update method */
struct SCIP_QuadConsUpgrade
{
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd));  /**< method to call for upgrading quadratic constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};
typedef struct SCIP_QuadConsUpgrade SCIP_QUADCONSUPGRADE; /**< quadratic constraint update method */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   replacebinaryprodlength; /**< length of linear term which when multiplied with a binary variable is replaced by an auxiliary variable and an equivalent linear formulation */
   int                   empathy4and;        /**< how much empathy we have for using the AND constraint handler: 0 avoid always; 1 use sometimes; 2 use as often as possible */
   SCIP_Bool             binreforminitial;   /**< whether to make constraints added due to replacing products with binary variables initial */
   SCIP_Real             mincutefficacysepa; /**< minimal efficacy of a cut in order to add it to relaxation during separation */
   SCIP_Real             mincutefficacyenfofac; /**< minimal target efficacy of a cut in order to add it to relaxation during enforcement as factor of feasibility tolerance (may be ignored) */
   SCIP_Bool             doscaling;          /**< should constraints be scaled in the feasibility check ? */
   SCIP_Real             cutmaxrange;        /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linearizeheursol;   /**< whether linearizations of convex quadratic constraints should be added to cutpool when some heuristics finds a new solution */
   SCIP_Bool             checkcurvature;     /**< whether functions should be checked for convexity/concavity */
   SCIP_Bool             checkfactorable;    /**< whether functions should be checked to be factorable */
   SCIP_Bool             linfeasshift;       /**< whether to make solutions in check feasible if possible */
   SCIP_Bool             disaggregate;       /**< whether to disaggregate quadratic constraints */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a single constraint within one round of SCIP propagation during solve */
   int                   maxproproundspresolve; /**< limit on number of propagation rounds for a single constraint within one presolving round */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */
   SCIP_Bool             enfocutsremovable;  /**< are cuts added during enforcement removable from the LP in the same node? */

   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subnlp heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */
   SCIP_EVENTHDLR*       eventhdlr;          /**< our handler for variable bound change events */
   int                   newsoleventfilterpos; /**< filter position of new solution event handler, if caught */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */
   SCIP_QUADCONSUPGRADE** quadconsupgrades;  /**< quadratic constraint upgrade methods for specializing quadratic constraints */
   int                   quadconsupgradessize; /**< size of quadconsupgrade array */
   int                   nquadconsupgrades;  /**< number of quadratic constraint upgrade methods */
};


/*
 * local methods for managing quadratic constraint update methods
 */


/** checks whether a quadratic constraint upgrade method has already be registered */
static
SCIP_Bool conshdlrdataHasUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(quadconsupgd != NULL);
   assert(conshdlrname != NULL);

   for( i = conshdlrdata->nquadconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->quadconsupgrades[i]->quadconsupgd == quadconsupgd )
      {
         SCIPwarningMessage(scip, "Try to add already known upgrade message for constraint handler <%s>.\n", conshdlrname);
         return TRUE;
      }
   }

   return FALSE;
}

/*
 * Local methods
 */

/** translate from one value of infinity to another
 * 
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/* catches variable bound change events on a linear variable in a quadratic constraint */
static
SCIP_RETCODE catchLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSDATA*  consdata;
   SCIP_QUADVAREVENTDATA* eventdata;
   SCIP_EVENTTYPE  eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventdata->consdata = consdata;
   eventdata->varidx = linvarpos;

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvars[linvarpos], eventtype, eventhdlr, (SCIP_EVENTDATA*)eventdata, &eventdata->filterpos) );

   consdata->lineventdata[linvarpos] = eventdata;

   return SCIP_OKAY;
}

/* drops variable bound change events on a linear variable in a quadratic constraint */
static
SCIP_RETCODE dropLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE  eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);
   assert(consdata->lineventdata[linvarpos] != NULL);
   assert(consdata->lineventdata[linvarpos]->consdata == consdata);
   assert(consdata->lineventdata[linvarpos]->varidx == linvarpos);
   assert(consdata->lineventdata[linvarpos]->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvars[linvarpos], eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata->lineventdata[linvarpos], consdata->lineventdata[linvarpos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->lineventdata[linvarpos]);  /*lint !e866 */

   return SCIP_OKAY;
}

/* catches variable bound change events on a quadratic variable in a quadratic constraint */
static
SCIP_RETCODE catchQuadVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   quadvarpos          /**< position of variable in quadratic variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUADVAREVENTDATA* eventdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(quadvarpos >= 0);
   assert(quadvarpos < consdata->nquadvars);
   assert(consdata->quadvarterms[quadvarpos].eventdata == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;
#ifdef CHECKIMPLINBILINEAR
   eventtype |= SCIP_EVENTTYPE_IMPLADDED;
#endif
   eventdata->consdata = consdata;
   eventdata->varidx   = -quadvarpos-1;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->quadvarterms[quadvarpos].var, eventtype, eventhdlr, (SCIP_EVENTDATA*)eventdata, &eventdata->filterpos) );

   consdata->quadvarterms[quadvarpos].eventdata = eventdata;

   return SCIP_OKAY;
}

/* catches variable bound change events on a quadratic variable in a quadratic constraint */
static
SCIP_RETCODE dropQuadVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   quadvarpos          /**< position of variable in quadratic variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(quadvarpos >= 0);
   assert(quadvarpos < consdata->nquadvars);
   assert(consdata->quadvarterms[quadvarpos].eventdata != NULL);
   assert(consdata->quadvarterms[quadvarpos].eventdata->consdata == consdata);
   assert(consdata->quadvarterms[quadvarpos].eventdata->varidx == -quadvarpos-1);
   assert(consdata->quadvarterms[quadvarpos].eventdata->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;
#ifdef CHECKIMPLINBILINEAR
   eventtype |= SCIP_EVENTTYPE_IMPLADDED;
#endif

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->quadvarterms[quadvarpos].var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata->quadvarterms[quadvarpos].eventdata, consdata->quadvarterms[quadvarpos].eventdata->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->quadvarterms[quadvarpos].eventdata);

   return SCIP_OKAY;
}

/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lineventdata == NULL);

   /* we will update isremovedfixings, so reset it to TRUE first */
   consdata->isremovedfixings = TRUE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize) );
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( catchLinearVarEvents(scip, eventhdlr, cons, i) );

      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(consdata->linvars[i]);
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      assert(consdata->quadvarterms[i].eventdata == NULL);

      SCIP_CALL( catchQuadVarEvents(scip, eventhdlr, cons, i) );

      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(consdata->quadvarterms[i].var);
   }

   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->lineventdata != NULL )
   {
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         if( consdata->lineventdata[i] != NULL )
         {
            SCIP_CALL( dropLinearVarEvents(scip, eventhdlr, cons, i) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize);
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].eventdata != NULL )
      {
         SCIP_CALL( dropQuadVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var,                /**< variable to lock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** unlocks a linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var,                /**< variable to unlock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** locks a quadratic variable in a constraint */
static
SCIP_RETCODE lockQuadraticVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var                 /**< variable to lock */
   )
{
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** unlocks a quadratic variable in a constraint */
static
SCIP_RETCODE unlockQuadraticVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var                 /**< variable to unlock */
   )
{
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** computes the minimal and maximal activity for the linear part in a constraint data
 *  only sums up terms that contribute finite values
 *  gives the number of terms that contribute infinite values
 *  only computes those activities where the corresponding side of the constraint is finite
 */
static
void consdataUpdateLinearActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty       /**< infinity value used in interval operations */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   SCIP_Real bnd;

   assert(scip != NULL);
   assert(consdata != NULL);

   /* if variable bounds are not strictly consistent, then the activity update methods may yield inconsistent activities
    * in this case, we also recompute the activities
    */
   if( consdata->minlinactivity != SCIP_INVALID && consdata->maxlinactivity != SCIP_INVALID &&  /*lint !e777 */
      (consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity) )
   {
      /* activities should be up-to-date */
      assert(consdata->minlinactivityinf >= 0);
      assert(consdata->maxlinactivityinf >= 0);
      return;
   }

   consdata->minlinactivityinf = 0;
   consdata->maxlinactivityinf = 0;

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   consdata->minlinactivity = SCIPisInfinity(scip,  consdata->rhs) ? -intervalinfty : 0.0;
   consdata->maxlinactivity = SCIPisInfinity(scip, -consdata->lhs) ?  intervalinfty : 0.0;

   if( consdata->nlinvars == 0 )
      return;

   /* if the activities computed here should be still up-to-date after bound changes,
    * variable events need to be caught */
   assert(consdata->lineventdata != NULL);

   prevroundmode = SCIPintervalGetRoundingMode();

   if( !SCIPisInfinity(scip,  consdata->rhs) )
   {
      /* compute minimal activity only if there is a finite right hand side */
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = MIN(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         else
         {
            bnd = MAX(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         consdata->minlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* compute maximal activity only if there is a finite left hand side */
      SCIPintervalSetRoundingModeUpwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = MAX(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         else
         {
            bnd = MIN(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         consdata->maxlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   SCIPintervalSetRoundingMode(prevroundmode);

   assert(consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity);
}

/** update the linear activities after a change in the lower bound of a variable */
static
void consdataUpdateLinearActivityLbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with lower bounds at infinity */
   assert(!SCIPisInfinity(scip, oldbnd));
   assert(!SCIPisInfinity(scip, newbnd));

   /* @todo since we check the linear activity for consistency later anyway, we may skip changing the rounding mode here */

   /* assume lhs <= a*x + y <= rhs, then the following bound changes can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      assert(!SCIPisInfinity(scip, -consdata->minlinactivity));

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      assert(!SCIPisInfinity(scip, consdata->maxlinactivity));

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** update the linear activities after a change in the upper bound of a variable */
static
void consdataUpdateLinearActivityUbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with upper bounds at -infinity */
   assert(!SCIPisInfinity(scip, -oldbnd));
   assert(!SCIPisInfinity(scip, -newbnd));

   /* @todo since we check the linear activity for consistency later anyway, we may skip changing the rounding mode here */

   /* assume lhs <= a*x + y <= rhs, then the following bound changes can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      assert(!SCIPisInfinity(scip, consdata->maxlinactivity));

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      assert(!SCIPisInfinity(scip, -consdata->minlinactivity));

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   int varidx;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   consdata = ((SCIP_QUADVAREVENTDATA*)eventdata)->consdata;
   assert(consdata != NULL);

   varidx = ((SCIP_QUADVAREVENTDATA*)eventdata)->varidx;
   assert(varidx <  0 ||  varidx   < consdata->nlinvars);
   assert(varidx >= 0 || -varidx-1 < consdata->nquadvars);

   eventtype = SCIPeventGetType(event);

   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      if( varidx < 0 )
      {
         /* mark activity bounds for quad term as not up to date anymore */
         SCIPintervalSetEmpty(&consdata->quadactivitybounds);
      }
      else
      {
         /* update activity bounds for linear terms */
         if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
            consdataUpdateLinearActivityLbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
         else
            consdataUpdateLinearActivityUbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
      }

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         consdata->ispropagated = FALSE;
   }

   if( eventtype & SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->isremovedfixings = FALSE;
   }

#ifdef CHECKIMPLINBILINEAR
   if( eventtype & SCIP_EVENTTYPE_IMPLADDED )
   {
      assert(varidx < 0); /* we catch impladded events only for quadratic variables */
      /* if variable is binary (quite likely if an implication has been added) and occurs in a bilinear term, then mark that we should check implications */
      if( SCIPvarIsBinary(SCIPeventGetVar(event)) && consdata->quadvarterms[-varidx-1].nadjbilin > 0 )
         consdata->isimpladded = TRUE;
   }
#endif

   return SCIP_OKAY;
}

/** ensures, that linear vars and coefs arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureLinearVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nlinvars <= consdata->linvarssize);

   if( num > consdata->linvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->linvars,  consdata->linvarssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lincoefs, consdata->linvarssize, newsize) );
      if( consdata->lineventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize, newsize) );
      }
      consdata->linvarssize = newsize;
   }
   assert(num <= consdata->linvarssize);

   return SCIP_OKAY;
}

/** ensures, that quadratic variable terms array can store at least num entries */
static
SCIP_RETCODE consdataEnsureQuadVarTermsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nquadvars <= consdata->quadvarssize);

   if( num > consdata->quadvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->quadvarterms, consdata->quadvarssize, newsize) );
      consdata->quadvarssize = newsize;
   }
   assert(num <= consdata->quadvarssize);

   return SCIP_OKAY;
}

/** ensures, that adjacency array can store at least num entries */
static
SCIP_RETCODE consdataEnsureAdjBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADVARTERM*     quadvarterm,        /**< quadratic variable term */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(quadvarterm != NULL);
   assert(quadvarterm->nadjbilin <= quadvarterm->adjbilinsize);

   if( num > quadvarterm->adjbilinsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &quadvarterm->adjbilin, quadvarterm->adjbilinsize, newsize) );
      quadvarterm->adjbilinsize = newsize;
   }
   assert(num <= quadvarterm->adjbilinsize);

   return SCIP_OKAY;
}

/** ensures, that bilinear term arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nbilinterms <= consdata->bilintermssize);

   if( num > consdata->bilintermssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->bilinterms, consdata->bilintermssize, newsize) );
      consdata->bilintermssize = newsize;
   }
   assert(num <= consdata->bilintermssize);

   return SCIP_OKAY;
}

/** creates empty constraint data structure */
static
SCIP_RETCODE consdataCreateEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< a buffer to store pointer to new constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->lhs = -SCIPinfinity(scip);
   (*consdata)->rhs =  SCIPinfinity(scip);

   (*consdata)->linvarssorted  = TRUE;
   (*consdata)->linvarsmerged  = TRUE;
   (*consdata)->quadvarssorted = TRUE;
   (*consdata)->quadvarsmerged = TRUE;
   (*consdata)->bilinsorted    = TRUE;
   (*consdata)->bilinmerged    = TRUE;

   (*consdata)->isremovedfixings = TRUE;
   (*consdata)->ispropagated     = TRUE;
   (*consdata)->initialmerge     = FALSE;

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   return SCIP_OKAY;
}

/** creates constraint data structure */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< a buffer to store pointer to new constraint data */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            lincoefs,           /**< array of coefficients of linear variables */
   int                   nquadvars,          /**< number of quadratic variables */
   SCIP_QUADVARTERM*     quadvarterms,       /**< array of quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_BILINTERM*       bilinterms,         /**< array of bilinear terms */
   SCIP_Bool             capturevars         /**< whether we should capture variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);

   assert(nlinvars == 0 || linvars  != NULL);
   assert(nlinvars == 0 || lincoefs != NULL);
   assert(nquadvars == 0 || quadvarterms != NULL);
   assert(nbilinterms == 0 || bilinterms != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;

   if( nlinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linvars,  linvars,  nlinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->lincoefs, lincoefs, nlinvars) );
      (*consdata)->nlinvars = nlinvars;
      (*consdata)->linvarssize = nlinvars;

      if( capturevars )
         for( i = 0; i < nlinvars; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, linvars[i]) );
         }
   }
   else
   {
      (*consdata)->linvarssorted = TRUE;
      (*consdata)->linvarsmerged = TRUE;
      (*consdata)->minlinactivity = 0.0;
      (*consdata)->maxlinactivity = 0.0;
      (*consdata)->minlinactivityinf = 0;
      (*consdata)->maxlinactivityinf = 0;
   }

   if( nquadvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->quadvarterms, quadvarterms, nquadvars) );

      for( i = 0; i < nquadvars; ++i )
      {
         (*consdata)->quadvarterms[i].eventdata = NULL;
         if( quadvarterms[i].nadjbilin )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->quadvarterms[i].adjbilin, quadvarterms[i].adjbilin, quadvarterms[i].nadjbilin) );
            (*consdata)->quadvarterms[i].adjbilinsize = quadvarterms[i].nadjbilin;
         }
         else
         {
            assert((*consdata)->quadvarterms[i].nadjbilin == 0);
            (*consdata)->quadvarterms[i].adjbilin = NULL;
            (*consdata)->quadvarterms[i].adjbilinsize = 0;
         }
         if( capturevars )
         {
            SCIP_CALL( SCIPcaptureVar(scip, quadvarterms[i].var) );
         }
      }

      (*consdata)->nquadvars = nquadvars;
      (*consdata)->quadvarssize = nquadvars;
      SCIPintervalSetEmpty(&(*consdata)->quadactivitybounds);
   }
   else
   {
      (*consdata)->quadvarssorted = TRUE;
      (*consdata)->quadvarsmerged = TRUE;
      SCIPintervalSet(&(*consdata)->quadactivitybounds, 0.0);
   }

   if( nbilinterms > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->bilinterms, bilinterms, nbilinterms) );
      (*consdata)->nbilinterms = nbilinterms;
      (*consdata)->bilintermssize = nbilinterms;
   }
   else
   {
      (*consdata)->bilinsorted = TRUE;
      (*consdata)->bilinmerged = TRUE;
   }

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->activity = SCIP_INVALID;
   (*consdata)->lhsviol  = SCIPisInfinity(scip, -lhs) ? 0.0 : SCIP_INVALID;
   (*consdata)->rhsviol  = SCIPisInfinity(scip,  rhs) ? 0.0 : SCIP_INVALID;

   return SCIP_OKAY;
}

/** frees constraint data structure */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data to free */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free sepa arrays, may exists if constraint is deleted in solving stage */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->sepaquadvars,     (*consdata)->nquadvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->sepabilinvar2pos, (*consdata)->nbilinterms);

   /* release linear variables and free linear part */
   if( (*consdata)->linvarssize > 0 )
   {
      for( i = 0; i < (*consdata)->nlinvars; ++i )
      {
         assert((*consdata)->lineventdata == NULL || (*consdata)->lineventdata[i] == NULL); /* variable events should have been dropped earlier */
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvars[i]) );
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->linvars,  (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->lincoefs, (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lineventdata, (*consdata)->linvarssize);
   }
   assert((*consdata)->linvars == NULL);
   assert((*consdata)->lincoefs == NULL);
   assert((*consdata)->lineventdata == NULL);

   /* release quadratic variables and free quadratic variable term part */
   for( i = 0; i < (*consdata)->nquadvars; ++i )
   {
      assert((*consdata)->quadvarterms[i].eventdata == NULL); /* variable events should have been dropped earlier */
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->quadvarterms[i].adjbilin, (*consdata)->quadvarterms[i].adjbilinsize);
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->quadvarterms[i].var) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->quadvarterms, (*consdata)->quadvarssize);

   /* free bilinear terms */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bilinterms, (*consdata)->bilintermssize);

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}

/** sorts linear part of constraint data */
static
void consdataSortLinearVars(
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->linvarssorted )
      return;

   if( consdata->nlinvars <= 1 )
   {
      consdata->linvarssorted = TRUE;
      return;
   }

   if( consdata->lineventdata == NULL )
   {
      SCIPsortPtrReal((void**)consdata->linvars, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);
   }
   else
   {
      int i;

      SCIPsortPtrPtrReal((void**)consdata->linvars, (void**)consdata->lineventdata, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);

      /* update variable indices in event data */
      for( i = 0; i < consdata->nlinvars; ++i )
         if( consdata->lineventdata[i] != NULL )
            consdata->lineventdata[i]->varidx = i;
   }

   consdata->linvarssorted = TRUE;
}

#if 0 /* no-one needs this routine currently */
/** returns the position of variable in the linear coefficients array of a constraint, or -1 if not found */
static
int consdataFindLinearVar(
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int pos;

   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->nlinvars == 0 )
      return -1;

   consdataSortLinearVars(consdata);

   if( !SCIPsortedvecFindPtr((void**)consdata->linvars, SCIPvarComp, (void*)var, consdata->nlinvars, &pos) )
      pos = -1;

   return pos;
}
#endif

/** index comparison method for quadratic variable terms: compares two indices of the quadratic variable set in the quadratic constraint */
static
SCIP_DECL_SORTINDCOMP(quadVarTermComp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nquadvars);
   assert(0 <= ind2 && ind2 < consdata->nquadvars);

   return SCIPvarCompare(consdata->quadvarterms[ind1].var, consdata->quadvarterms[ind2].var);
}

/** sorting of quadratic variable terms */
static
SCIP_RETCODE consdataSortQuadVarTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   int* perm;
   int  i;
   int  nexti;
   int  v;
   SCIP_QUADVARTERM quadterm;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->quadvarssorted )
      return SCIP_OKAY;

   if( consdata->nquadvars == 0 )
   {
      consdata->quadvarssorted = TRUE;
      return SCIP_OKAY;
   }

   /* get temporary memory to store the sorted permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nquadvars) );

   /* call bubble sort */
   SCIPsort(perm, quadVarTermComp, (void*)consdata, consdata->nquadvars);

   /* permute the quadratic variable terms according to the resulting permutation */
   for( v = 0; v < consdata->nquadvars; ++v )
   {
      if( perm[v] != v )
      {
         quadterm = consdata->quadvarterms[v];

         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < consdata->nquadvars);
            assert(perm[i] != i);
            consdata->quadvarterms[i] = consdata->quadvarterms[perm[i]];
            if( consdata->quadvarterms[i].eventdata != NULL )
            {
               consdata->quadvarterms[i].eventdata->varidx = -i-1;
            }
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->quadvarterms[i] = quadterm;
         if( consdata->quadvarterms[i].eventdata != NULL )
         {
            consdata->quadvarterms[i].eventdata->varidx = -i-1;
         }
         perm[i] = i;
      }
   }
   consdata->quadvarssorted = TRUE;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/** returns the position of variable in the quadratic variable terms array of a constraint, or -1 if not found */
static
SCIP_RETCODE consdataFindQuadVarTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   SCIP_VAR*             var,                /**< variable to search for */
   int*                  pos                 /**< buffer where to store position of var in quadvarterms array, or -1 if not found */
   )
{
   int left;
   int right;
   int cmpres;

   assert(consdata != NULL);
   assert(var != NULL);
   assert(pos != NULL);

   if( consdata->nquadvars == 0 )
   {
      *pos = -1;
      return SCIP_OKAY;
   }

   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   left = 0;
   right = consdata->nquadvars - 1;
   while( left <= right )
   {
      int middle;

      middle = (left+right)/2;
      assert(0 <= middle && middle < consdata->nquadvars);

      cmpres = SCIPvarCompare(var, consdata->quadvarterms[middle].var);

      if( cmpres < 0 )
         right = middle - 1;
      else if( cmpres > 0 )
         left  = middle + 1;
      else
      {
         *pos = middle;
         return SCIP_OKAY;
      }
   }
   assert(left == right+1);

   *pos = -1;

   return SCIP_OKAY;
}

/** index comparison method for bilinear terms: compares two index pairs of the bilinear term set in the quadratic constraint */
static
SCIP_DECL_SORTINDCOMP(bilinTermComp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;
   int var1cmp;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nbilinterms);
   assert(0 <= ind2 && ind2 < consdata->nbilinterms);

   var1cmp = SCIPvarCompare(consdata->bilinterms[ind1].var1, consdata->bilinterms[ind2].var1);
   if( var1cmp != 0 )
      return var1cmp;

   return SCIPvarCompare(consdata->bilinterms[ind1].var2, consdata->bilinterms[ind2].var2);
}

/** sorting of bilinear terms */
static
SCIP_RETCODE consdataSortBilinTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   int* perm;
   int* invperm;
   int  i;
   int  nexti;
   int  v;
   SCIP_BILINTERM bilinterm;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->bilinsorted )
      return SCIP_OKAY;

   if( consdata->nbilinterms == 0 )
   {
      consdata->bilinsorted = TRUE;
      return SCIP_OKAY;
   }

   /* get temporary memory to store the sorted permutation and the inverse permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm,    consdata->nbilinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &invperm, consdata->nbilinterms) );

   /* call bubble sort */
   SCIPsort(perm, bilinTermComp, (void*)consdata, consdata->nbilinterms);

   /* compute inverted permutation */
   for( v = 0; v < consdata->nbilinterms; ++v )
   {
      assert(0 <= perm[v] && perm[v] < consdata->nbilinterms);
      invperm[perm[v]] = v;
   }

   /* permute the bilinear terms according to the resulting permutation */
   for( v = 0; v < consdata->nbilinterms; ++v )
   {
      if( perm[v] != v )
      {
         bilinterm = consdata->bilinterms[v];

         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < consdata->nbilinterms);
            assert(perm[i] != i);
            consdata->bilinterms[i] = consdata->bilinterms[perm[i]];
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->bilinterms[i] = bilinterm;
         perm[i] = i;
      }
   }

   /* update the adjacency information in the quadratic variable terms */
   for( v = 0; v < consdata->nquadvars; ++v )
      for( i = 0; i < consdata->quadvarterms[v].nadjbilin; ++i )
         consdata->quadvarterms[v].adjbilin[i] = invperm[consdata->quadvarterms[v].adjbilin[i]];

   consdata->bilinsorted = TRUE;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &invperm);

   return SCIP_OKAY;
}

/** moves a linear variable from one position to another */
static
void consdataMoveLinearVar(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nlinvars);
   assert(newpos >= 0);
   assert(newpos < consdata->linvarssize);

   if( newpos == oldpos )
      return;

   consdata->linvars [newpos] = consdata->linvars [oldpos];
   consdata->lincoefs[newpos] = consdata->lincoefs[oldpos];

   if( consdata->lineventdata != NULL )
   {
      assert(newpos >= consdata->nlinvars || consdata->lineventdata[newpos] == NULL);

      consdata->lineventdata[newpos] = consdata->lineventdata[oldpos];
      consdata->lineventdata[newpos]->varidx = newpos;

      consdata->lineventdata[oldpos] = NULL;
   }

   consdata->linvarssorted = FALSE;
}   

/** moves a quadratic variable from one position to another */
static
void consdataMoveQuadVarTerm(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nquadvars);
   assert(newpos >= 0);
   assert(newpos < consdata->quadvarssize);

   if( newpos == oldpos )
      return;

   assert(newpos >= consdata->nquadvars || consdata->quadvarterms[newpos].eventdata == NULL);

   consdata->quadvarterms[newpos] = consdata->quadvarterms[oldpos];

   if( consdata->quadvarterms[newpos].eventdata != NULL )
   {
      consdata->quadvarterms[newpos].eventdata->varidx = -newpos-1;
      consdata->quadvarterms[oldpos].eventdata = NULL;
   }

   consdata->quadvarssorted = FALSE;
}   

/** adds linear coefficient in quadratic constraint */
static
SCIP_RETCODE addLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             coef                /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, coef) )
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

   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars+1) );
   consdata->linvars [consdata->nlinvars] = var;
   consdata->lincoefs[consdata->nlinvars] = coef;

   ++consdata->nlinvars;

   /* catch variable events */
   if( consdata->lineventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      consdata->lineventdata[consdata->nlinvars-1] = NULL;

      /* catch bound change events of variable */
      SCIP_CALL( catchLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, consdata->nlinvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockLinearVariable(scip, cons, var, coef) );

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;
   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var);
   if( consdata->nlinvars == 1 )
      consdata->linvarssorted = TRUE;
   else
      consdata->linvarssorted = consdata->linvarssorted && (SCIPvarCompare(consdata->linvars[consdata->nlinvars-2], consdata->linvars[consdata->nlinvars-1]) == -1);
   /* always set too FALSE since the new linear variable should be checked if already existing as quad var term */ 
   consdata->linvarsmerged = FALSE;

   return SCIP_OKAY;
}

/** deletes linear coefficient at given position from quadratic constraint data */
static
SCIP_RETCODE delLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nlinvars);

   var  = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );

   /* if we catch variable events, drop the events on the variable */
   if( consdata->lineventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvars[pos]) );

   /* move the last variable to the free slot */
   consdataMoveLinearVar(consdata, consdata->nlinvars-1, pos);

   --consdata->nlinvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;

   return SCIP_OKAY;
}

/** changes linear coefficient value at given position of quadratic constraint */
static
SCIP_RETCODE chgLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos,                /**< position of linear coefficient to change */
   SCIP_Real             newcoef             /**< new value of linear coefficient */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newcoef));

   conshdlrdata = NULL;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos);
   assert(pos < consdata->nlinvars);
   assert(!SCIPisZero(scip, newcoef));

   var = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if necessary, remove the rounding locks and event catching of the variable */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         assert(SCIPconsIsTransformed(cons));

         /* remove rounding locks for variable with old coefficient */
         SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );
      }

      if( consdata->lineventdata[pos] != NULL )
      {
         /* get event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* drop bound change events of variable */
         SCIP_CALL( dropLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
      }
   }

   /* change the coefficient */
   consdata->lincoefs[pos] = newcoef;

   /* if necessary, install the rounding locks and event catching of the variable again */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         /* install rounding locks for variable with new coefficient */
         SCIP_CALL( lockLinearVariable(scip, cons, var, newcoef) );
      }

      if( conshdlrdata != NULL )
      {
         /* catch bound change events of variable */
         SCIP_CALL( catchLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
      }
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   return SCIP_OKAY;
}

/** adds quadratic variable term to quadratic constraint */
static
SCIP_RETCODE addQuadVarTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef,            /**< square coefficient of variable */
   SCIP_Bool             catchevents         /**< whether we should catch variable events */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   SCIP_QUADVARTERM* quadvarterm;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

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

   SCIP_CALL( consdataEnsureQuadVarTermsSize(scip, consdata, consdata->nquadvars+1) );

   quadvarterm = &consdata->quadvarterms[consdata->nquadvars];
   quadvarterm->var       = var;
   quadvarterm->lincoef   = lincoef;
   quadvarterm->sqrcoef   = sqrcoef;
   quadvarterm->adjbilinsize = 0;
   quadvarterm->nadjbilin = 0;
   quadvarterm->adjbilin  = NULL;
   quadvarterm->eventdata = NULL;

   ++consdata->nquadvars;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* catch variable events, if we do so */
   if( catchevents )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variable */
      SCIP_CALL( catchQuadVarEvents(scip, conshdlrdata->eventhdlr, cons, consdata->nquadvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockQuadraticVariable(scip, cons, var) );

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;
   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var);
   if( consdata->nquadvars == 1 )
      consdata->quadvarssorted = TRUE;
   else
      consdata->quadvarssorted = consdata->quadvarssorted && (SCIPvarCompare(consdata->quadvarterms[consdata->nquadvars-2].var, consdata->quadvarterms[consdata->nquadvars-1].var) == -1);
   /* also set to FALSE if nquadvars == 1, since the new variable should be checked for linearity and other stuff in mergeAndClean ... */ 
   consdata->quadvarsmerged = FALSE;

   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** deletes quadratic variable term at given position from quadratic constraint data */
static
SCIP_RETCODE delQuadVarTermPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos                 /**< position of term to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nquadvars);

   var = consdata->quadvarterms[pos].var;
   assert(var != NULL);
   assert(consdata->quadvarterms[pos].nadjbilin == 0);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockQuadraticVariable(scip, cons, var) );

   /* if we catch variable events, drop the events on the variable */
   if( consdata->quadvarterms[pos].eventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropQuadVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->quadvarterms[pos].var) );

   /* free adjacency array */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->quadvarterms[pos].adjbilin, consdata->quadvarterms[pos].adjbilinsize);

   /* move the last variable term to the free slot */
   consdataMoveQuadVarTerm(consdata, consdata->nquadvars-1, pos);

   --consdata->nquadvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** replace variable in quadratic variable term at given position of quadratic constraint data 
 * allows to replace x by coef*y+offset, thereby maintaining linear and square coefficients and bilinear terms */
static
SCIP_RETCODE replaceQuadVarTermPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos,                /**< position of term to replace */
   SCIP_VAR*             var,                /**< new variable */
   SCIP_Real             coef,               /**< linear coefficient of new variable */
   SCIP_Real             offset              /**< offset of new variable */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real constant;

   int i;
   SCIP_VAR* var2;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(pos >= 0);
   assert(pos <  consdata->nquadvars);

   quadvarterm = &consdata->quadvarterms[pos];

   /* remove rounding locks for old variable */
   SCIP_CALL( unlockQuadraticVariable(scip, cons, quadvarterm->var) );

   /* if we catch variable events, drop the events on the old variable */
   if( quadvarterm->eventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      eventhdlr = conshdlrdata->eventhdlr;

      /* drop bound change events of variable */
      SCIP_CALL( dropQuadVarEvents(scip, eventhdlr, cons, pos) );
   }
   else
   {
      eventhdlr = NULL;
   }

   /* compute constant and put into lhs/rhs */
   constant = quadvarterm->lincoef * offset + quadvarterm->sqrcoef * offset * offset;
   if( constant != 0.0 )
   {
      /* maintain constant part */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->lhs -= constant;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         consdata->rhs -= constant;
   }

   /* update linear and square coefficient */
   quadvarterm->lincoef *= coef;
   quadvarterm->lincoef += 2.0 * quadvarterm->sqrcoef * coef * offset;
   quadvarterm->sqrcoef *= coef * coef;

   /* update bilinear terms */
   for( i = 0; i < quadvarterm->nadjbilin; ++i )
   {
      bilinterm = &consdata->bilinterms[quadvarterm->adjbilin[i]];

      if( bilinterm->var1 == quadvarterm->var )
      {
         bilinterm->var1 = var;
         var2 = bilinterm->var2;
      }
      else
      {
         assert(bilinterm->var2 == quadvarterm->var);
         bilinterm->var2 = var;
         var2 = bilinterm->var1;
      }

      if( var == var2 )
      {
         /* looks like we actually have a square term here */
         quadvarterm->lincoef += bilinterm->coef * offset;
         quadvarterm->sqrcoef += bilinterm->coef * coef;
         /* deleting bilinear terms is expensive, since it requires updating adjacency information
          * thus, for now we just set the coefficient to 0.0 and clear in later when the bilinear terms are merged */
         bilinterm->coef = 0.0;
         continue;
      }

      /* swap var1 and var2 if they are in wrong order */
      if( SCIPvarCompare(bilinterm->var1, bilinterm->var2) > 0 )
      {
         SCIP_VAR* tmp;
         tmp = bilinterm->var1;
         bilinterm->var1 = bilinterm->var2;
         bilinterm->var2 = tmp;
      }
      assert(SCIPvarCompare(bilinterm->var1, bilinterm->var2) == -1);

      if( offset != 0.0 )
      {
         /* need to find var2 and add offset*bilinterm->coef to linear coefficient */
         int var2pos;

         var2pos = 0;
         while( consdata->quadvarterms[var2pos].var != var2 )
         {
            ++var2pos;
            assert(var2pos < consdata->nquadvars);
         }

         consdata->quadvarterms[var2pos].lincoef += bilinterm->coef * offset;
      }

      bilinterm->coef *= coef;
   }

   /* release old variable */
   SCIP_CALL( SCIPreleaseVar(scip, &quadvarterm->var) );

   /* set new variable */
   quadvarterm->var = var;

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, quadvarterm->var) );

   /* catch variable events, if we do so */
   if( eventhdlr != NULL )
   {
      /* catch bound change events of variable */
      SCIP_CALL( catchQuadVarEvents(scip, eventhdlr, cons, pos) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockQuadraticVariable(scip, cons, var) );

   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var);
   consdata->quadvarssorted = (consdata->nquadvars == 1);
   consdata->quadvarsmerged = FALSE;
   consdata->bilinmerged &= (quadvarterm->nadjbilin == 0);  /*lint !e514*/

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** adds a bilinear term to quadratic constraint */
static
SCIP_RETCODE addBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   var1pos,            /**< position of first variable in quadratic variables array */
   int                   var2pos,            /**< position of second variable in quadratic variables array */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;

   assert(scip != NULL);
   assert(cons != NULL);

   if( var1pos == var2pos )
   {
      SCIPerrorMessage("tried to add bilinear term where both variables are the same\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(var1pos >= 0);
   assert(var1pos < consdata->nquadvars);
   assert(var2pos >= 0);
   assert(var2pos < consdata->nquadvars);

   SCIP_CALL( consdataEnsureBilinSize(scip, consdata, consdata->nbilinterms + 1) );

   bilinterm = &consdata->bilinterms[consdata->nbilinterms];
   if( SCIPvarCompare(consdata->quadvarterms[var1pos].var, consdata->quadvarterms[var2pos].var) < 0 )
   {
      bilinterm->var1 = consdata->quadvarterms[var1pos].var;
      bilinterm->var2 = consdata->quadvarterms[var2pos].var;
   }
   else
   {
      bilinterm->var1 = consdata->quadvarterms[var2pos].var;
      bilinterm->var2 = consdata->quadvarterms[var1pos].var;
   }
   bilinterm->coef = coef;

   if( bilinterm->var1 == bilinterm->var2 )
   {
      SCIPerrorMessage("tried to add bilinear term where both variables are the same, but appear at different positions in quadvarterms array\n");
      return SCIP_INVALIDDATA;
   }
   assert(SCIPvarCompare(bilinterm->var1, bilinterm->var2) == -1);

   SCIP_CALL( consdataEnsureAdjBilinSize(scip, &consdata->quadvarterms[var1pos], consdata->quadvarterms[var1pos].nadjbilin + 1) );
   SCIP_CALL( consdataEnsureAdjBilinSize(scip, &consdata->quadvarterms[var2pos], consdata->quadvarterms[var2pos].nadjbilin + 1) );

   consdata->quadvarterms[var1pos].adjbilin[consdata->quadvarterms[var1pos].nadjbilin] = consdata->nbilinterms;
   consdata->quadvarterms[var2pos].adjbilin[consdata->quadvarterms[var2pos].nadjbilin] = consdata->nbilinterms;
   ++consdata->quadvarterms[var1pos].nadjbilin;
   ++consdata->quadvarterms[var2pos].nadjbilin;

   ++consdata->nbilinterms;

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;
   if( consdata->nbilinterms == 1 )
   {
      consdata->bilinsorted = TRUE;
      consdata->bilinmerged = TRUE;
   }
   else
   {
      consdata->bilinsorted = consdata->bilinsorted
         && (bilinTermComp(consdata, consdata->nbilinterms-2, consdata->nbilinterms-1) >= 0);
      consdata->bilinmerged = FALSE;
   }

   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** removes a set of bilinear terms and updates adjacency information in quad var terms
 * Note: this function sorts the given array termposs */
static
SCIP_RETCODE removeBilinearTermsPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   nterms,             /**< number of terms to delete */
   int*                  termposs            /**< indices of terms to delete */
   )
{
   SCIP_CONSDATA* consdata;
   int* newpos;
   int i;
   int j;
   int offset;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nterms == 0 || termposs != NULL);

   if( nterms == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPsortInt(termposs, nterms);

   SCIP_CALL( SCIPallocBufferArray(scip, &newpos, consdata->nbilinterms) );

   i = 0;
   offset = 0;
   for( j = 0; j < consdata->nbilinterms; ++j )
   {
      /* if j'th term is deleted, increase offset and continue */
      if( i < nterms && j == termposs[i] )
      {
         ++offset;
         ++i;
         newpos[j] = -1;
         continue;
      }

      /* otherwise, move it forward and remember new position */
      if( offset > 0 )
         consdata->bilinterms[j-offset] = consdata->bilinterms[j];
      newpos[j] = j - offset;
   }
   assert(offset == nterms);

   /* update adjacency and activity information in quad var terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      offset = 0;
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         assert(consdata->quadvarterms[i].adjbilin[j] < consdata->nbilinterms);
         if( newpos[consdata->quadvarterms[i].adjbilin[j]] == -1 )
         {
            /* corresponding bilinear term was deleted, thus increase offset */
            ++offset;
         }
         else
         {
            /* update index of j'th bilinear term and store at position j-offset */
            consdata->quadvarterms[i].adjbilin[j-offset] = newpos[consdata->quadvarterms[i].adjbilin[j]];
         }
      }
      consdata->quadvarterms[i].nadjbilin -= offset;
      /* some bilinear term was removed, so invalidate activity bounds */
   }

   consdata->nbilinterms -= nterms;

   SCIPfreeBufferArray(scip, &newpos);

   /* some quad vars may be linear now */
   consdata->quadvarsmerged = FALSE;

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   return SCIP_OKAY;
}

/* merges quad var terms that correspond to the same variable and does additional cleanup
 * if a quadratic variable terms is actually linear, makes a linear term out of it
 * also replaces squares of binary variables by the binary variables, i.e., adds sqrcoef to lincoef
 */
static
SCIP_RETCODE mergeAndCleanQuadVarTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->quadvarsmerged )
      return SCIP_OKAY;

   if( consdata->nquadvars == 0 )
   {
      consdata->quadvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nquadvars )
   {
      /* make sure quad var terms are sorted (do this in every round, since we may move variables around) */
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

      quadvarterm = &consdata->quadvarterms[i];

      for( j = i+1; j < consdata->nquadvars && consdata->quadvarterms[j].var == quadvarterm->var; ++j )
      {
         /* add quad var term j to current term i */
         quadvarterm->lincoef += consdata->quadvarterms[j].lincoef;
         quadvarterm->sqrcoef += consdata->quadvarterms[j].sqrcoef;
         if( consdata->quadvarterms[j].nadjbilin > 0 )
         {
            /* move adjacency information from j to i */
            SCIP_CALL( consdataEnsureAdjBilinSize(scip, quadvarterm, quadvarterm->nadjbilin + consdata->quadvarterms[j].nadjbilin) );
            BMScopyMemoryArray(&quadvarterm->adjbilin[quadvarterm->nadjbilin], consdata->quadvarterms[j].adjbilin, consdata->quadvarterms[j].nadjbilin);  /*lint !e866*/
            quadvarterm->nadjbilin += consdata->quadvarterms[j].nadjbilin;
            consdata->quadvarterms[j].nadjbilin = 0;
         }
         consdata->quadvarterms[j].lincoef = 0.0;
         consdata->quadvarterms[j].sqrcoef = 0.0;
         /* mark that activity information in quadvarterm is not up to date anymore */
      }

      /* remove quad var terms i+1..j-1 backwards */
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delQuadVarTermPos(scip, cons, j) );
      }

      /* for binary variables, x^2 = x
       * however, we may destroy convexity of a quadratic term that involves also bilinear terms
       * thus, we do this step only if the variable does not appear in any bilinear term */
      if( quadvarterm->sqrcoef != 0.0 && SCIPvarIsBinary(quadvarterm->var) && quadvarterm->nadjbilin == 0 )
      {
         SCIPdebugMessage("replace square of binary variable by itself: <%s>^2 --> <%s>\n", SCIPvarGetName(quadvarterm->var), SCIPvarGetName(quadvarterm->var));
         quadvarterm->lincoef += quadvarterm->sqrcoef;
         quadvarterm->sqrcoef = 0.0;

         /* invalidate nonlinear row */
         if( consdata->nlrow != NULL )
         {
            SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
         }
      }

      /* if its 0.0 or linear, get rid of it */
      if( SCIPisZero(scip, quadvarterm->sqrcoef) && quadvarterm->nadjbilin == 0 )
      {
         if( !SCIPisZero(scip, quadvarterm->lincoef) )
         {
            /* seem to be a linear term now, thus add as linear term */
            SCIP_CALL( addLinearCoef(scip, cons, quadvarterm->var, quadvarterm->lincoef) );
         }
         /* remove term at pos i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
      }
      else
      {
         ++i;
      }
   }

   consdata->quadvarsmerged = TRUE;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   return SCIP_OKAY;
}

/* merges entries with same linear variable into one entry and cleans up entries with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanLinearVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newcoef;
   int i;
   int j;
   int qvarpos;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->linvarsmerged )
      return SCIP_OKAY;

   if( consdata->nlinvars == 0 )
   {
      consdata->linvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nlinvars )
   {
      /* make sure linear variables are sorted (do this in every round, since we may move variables around) */
      consdataSortLinearVars(consdata);

      /* sum up coefficients that correspond to variable i */
      newcoef = consdata->lincoefs[i];
      for( j = i+1; j < consdata->nlinvars && consdata->linvars[i] == consdata->linvars[j]; ++j )
         newcoef += consdata->lincoefs[j];
      /* delete the additional variables in backward order */ 
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, j) );
      }

      /* check if there is already a quadratic variable term with this variable */
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->linvars[i], &qvarpos) );
      if( qvarpos >= 0)
      {
         /* add newcoef to linear coefficient of quadratic variable and mark linear variable as to delete */
         assert(qvarpos < consdata->nquadvars);
         assert(consdata->quadvarterms[qvarpos].var == consdata->linvars[i]);
         consdata->quadvarterms[qvarpos].lincoef += newcoef;
         newcoef = 0.0;
         SCIPintervalSetEmpty(&consdata->quadactivitybounds);
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, newcoef) )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );
      }
      else
      {
         SCIP_CALL( chgLinearCoefPos(scip, cons, i, newcoef) );
         ++i;
      }
   }

   consdata->linvarsmerged = TRUE;

   return SCIP_OKAY;
}

/* merges bilinear terms with same variables into a single term, removes bilinear terms with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanBilinearTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   int i;
   int j;
   int* todelete;
   int ntodelete;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->bilinmerged )
      return SCIP_OKAY;

   if( consdata->nbilinterms == 0 )
   {
      consdata->bilinmerged = TRUE;
      return SCIP_OKAY;
   }

   /* alloc memory for array of terms that need to be deleted finally */
   ntodelete = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &todelete, consdata->nbilinterms) );

   /* make sure bilinear terms are sorted */
   SCIP_CALL( consdataSortBilinTerms(scip, consdata) );

   i = 0;
   while( i < consdata->nbilinterms )
   {
      bilinterm = &consdata->bilinterms[i];

      /* sum up coefficients that correspond to same variables as term i */
      for( j = i+1; j < consdata->nbilinterms && bilinterm->var1 == consdata->bilinterms[j].var1 && bilinterm->var2 == consdata->bilinterms[j].var2; ++j )
      {
         bilinterm->coef += consdata->bilinterms[j].coef;
         todelete[ntodelete++] = j;
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, bilinterm->coef) )
      {
         todelete[ntodelete++] = i;
      }

      /* continue with term after the current series */
      i = j;
   }

   /* delete bilinear terms */
   SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );

   SCIPfreeBufferArray(scip, &todelete);

   consdata->bilinmerged = TRUE;

   return SCIP_OKAY;
}

/** removes fixes (or aggregated) variables from a quadratic constraint */
static
SCIP_RETCODE removeFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real coef;
   SCIP_Real offset;
   SCIP_VAR* var;
   SCIP_VAR* var2;
   int var2pos;
   int i;
   int j;
   int k;

   SCIP_Bool have_change;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   have_change = FALSE;
   i = 0;
   while( i < consdata->nlinvars )
   {
      var = consdata->linvars[i];

      if( SCIPvarIsActive(var) )
      {
         ++i;
         continue;
      }

      have_change = TRUE;

      coef = consdata->lincoefs[i];
      offset = 0.0;

      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &offset) );

      SCIPdebugMessage("  linear term %g*<%s> is replaced by %g * <%s> + %g\n", consdata->lincoefs[i], SCIPvarGetName(consdata->linvars[i]), coef, SCIPvarGetName(var), offset);

      /* delete previous variable (this will move another variable to position i) */
      SCIP_CALL( delLinearCoefPos(scip, cons, i) );

      /* put constant part into bounds */
      if( offset != 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs -= offset;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhs -= offset;
      }

      /* nothing left to do if variable had been fixed */
      if( coef == 0.0 )
         continue;

      /* if GetProbvar gave a linear variable, just add it
       * if it's a multilinear variable, add it's disaggregated variables */
      if( SCIPvarIsActive(var) )
      {
         SCIP_CALL( addLinearCoef(scip, cons, var, coef) );
      }
      else
      {
         int        naggrs;
         SCIP_VAR** aggrvars;
         SCIP_Real* aggrscalars;
         SCIP_Real  aggrconstant;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         naggrs = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrscalars = SCIPvarGetMultaggrScalars(var);
         aggrconstant = SCIPvarGetMultaggrConstant(var);

         SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + naggrs) );

         for( j = 0; j < naggrs; ++j )
         {
            SCIP_CALL( addLinearCoef(scip, cons, aggrvars[j], coef * aggrscalars[j]) );
         }

         if( aggrconstant != 0.0 )
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= coef * aggrconstant;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= coef * aggrconstant;
         }
      }
   }

   i = 0;
   while( i < consdata->nquadvars )
   {
      var = consdata->quadvarterms[i].var;

      if( SCIPvarIsActive(var) )
      {
         ++i;
         continue;
      }

      have_change = TRUE;

      coef   = 1.0;
      offset = 0.0;
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &offset) );

      SCIPdebugMessage("  quadratic variable <%s> with status %d is replaced by %g * <%s> + %g\n", SCIPvarGetName(consdata->quadvarterms[i].var), SCIPvarGetStatus(consdata->quadvarterms[i].var), coef, SCIPvarGetName(var), offset);

      /* handle fixed variable */
      if( coef == 0.0 )
      {
         /* if not fixed to 0.0, add to linear coefs of vars in bilinear terms, and deal with linear and square term as constant */ 
         if( offset != 0.0 )
         {
            for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
            {
               bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[j]];

               var2 = bilinterm->var1 == consdata->quadvarterms[i].var ? bilinterm->var2 : bilinterm->var1;
               assert(var2 != consdata->quadvarterms[i].var);

               var2pos = 0;
               while( consdata->quadvarterms[var2pos].var != var2 )
               {
                  ++var2pos;
                  assert(var2pos < consdata->nquadvars);
               }
               consdata->quadvarterms[var2pos].lincoef += bilinterm->coef * offset;
            }

            offset = consdata->quadvarterms[i].lincoef * offset + consdata->quadvarterms[i].sqrcoef * offset * offset;
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= offset;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= offset;
         }

         /* remove bilinear terms */
         SCIP_CALL( removeBilinearTermsPos(scip, cons, consdata->quadvarterms[i].nadjbilin, consdata->quadvarterms[i].adjbilin) );

         /* delete quad. var term i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );

         continue;
      }

      assert(var != NULL);

      /* if GetProbvar gave an active variable, replace the quad var term so that it uses the new variable */
      if( SCIPvarIsActive(var) )
      {
         /* replace x by coef*y+offset */
         SCIP_CALL( replaceQuadVarTermPos(scip, cons, i, var, coef, offset) );

         continue;
      }
      else
      {
         /* if GetProbVar gave a multi-aggregated variable, add new quad var terms and new bilinear terms
          * x is replaced by coef * (sum_i a_ix_i + b) + offset
          * lcoef * x + scoef * x^2 + bcoef * x * y ->
          *   (b*coef + offset) * (lcoef + (b*coef + offset) * scoef)
          * + sum_i a_i*coef * (lcoef + 2 (b*coef + offset) * scoef) x_i
          * + sum_i (a_i*coef)^2 * scoef * x_i^2
          * + 2 sum_{i,j, i<j} (a_i a_j coef^2 scoef) x_i x_j
          * + bcoef * (b*coef + offset + coef * sum_i a_ix_i) y 
          */
         int        naggrs;
         SCIP_VAR** aggrvars;     /* x_i */
         SCIP_Real* aggrscalars;  /* a_i */
         SCIP_Real  aggrconstant; /* b */
         int nquadtermsold;

         SCIP_Real lcoef;
         SCIP_Real scoef;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         naggrs = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrscalars = SCIPvarGetMultaggrScalars(var);
         aggrconstant = SCIPvarGetMultaggrConstant(var);

         lcoef = consdata->quadvarterms[i].lincoef;
         scoef = consdata->quadvarterms[i].sqrcoef;

         nquadtermsold = consdata->nquadvars;

         SCIP_CALL( consdataEnsureQuadVarTermsSize(scip, consdata, consdata->nquadvars + naggrs) );

         /* take care of constant part */
         if( aggrconstant != 0.0 || offset != 0.0 )
         {
            SCIP_Real constant;
            constant = (aggrconstant * coef + offset) * (lcoef + (aggrconstant * coef + offset) * scoef);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= constant;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= constant;
         }

         /* add x_i's with linear and square coefficients */
         for( j = 0; j < naggrs; ++j )
         {
            SCIP_CALL( addQuadVarTerm(scip, cons, aggrvars[j],
                  coef * aggrscalars[j] * (lcoef + 2.0 * scoef * (coef * aggrconstant + offset)),
                  coef * coef * aggrscalars[j] * aggrscalars[j] * scoef,
                  TRUE) );
         }

         /* ensure space for bilinear terms */
         SCIP_CALL( consdataEnsureBilinSize(scip, consdata, consdata->nquadvars + (scoef != 0.0 ? (naggrs * (naggrs-1))/2 : 0) + consdata->quadvarterms[j].nadjbilin * naggrs) );

         /* add x_j*x_k's */
         if( scoef != 0.0 )
         {
            for( j = 0; j < naggrs; ++j )
               for( k = 0; k < j; ++k )
               {
                  assert(aggrvars[j] != aggrvars[k]);
                  SCIP_CALL( addBilinearTerm(scip, cons, nquadtermsold + j, nquadtermsold + k, 
                        2.0 * aggrscalars[j] * aggrscalars[k] * coef * coef * scoef) );
               }
         }

         /* add x_i*y's */
         for( k = 0; k < consdata->quadvarterms[i].nadjbilin; ++k )
         {
            bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[k]];
            var2 = (bilinterm->var1 == consdata->quadvarterms[i].var) ? bilinterm->var2 : bilinterm->var1;
            assert(var2 != consdata->quadvarterms[i].var);

            /* this is not efficient, but we cannot sort the quadratic terms here, since we currently iterate over them */
            var2pos = 0;
            while( consdata->quadvarterms[var2pos].var != var2 )
            {
               ++var2pos;
               assert(var2pos < consdata->nquadvars);
            }

            for( j = 0; j < naggrs; ++j )
            {
               if( aggrvars[j] == var2 )
               { /* x_i == y, so we have a square term here */
                  consdata->quadvarterms[var2pos].sqrcoef += bilinterm->coef * coef * aggrscalars[j];
               }
               else
               { /* x_i != y, so we need to add a bilinear term here */
                  SCIP_CALL( addBilinearTerm(scip, cons, nquadtermsold + j, var2pos,
                        bilinterm->coef * coef * aggrscalars[j]) );
               }
            }

            consdata->quadvarterms[var2pos].lincoef += bilinterm->coef * (aggrconstant * coef + offset);
         }

         /* remove bilinear terms */
         SCIP_CALL( removeBilinearTermsPos(scip, cons, consdata->quadvarterms[i].nadjbilin, consdata->quadvarterms[i].adjbilin) );

         /* delete quad. var term i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
      }
   }

   consdata->isremovedfixings = TRUE;

   SCIPdebugMessage("removed fixations from <%s>\n  -> ", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

#ifndef NDEBUG
   for( i = 0; i < consdata->nlinvars; ++i )
      assert(SCIPvarIsActive(consdata->linvars[i]));

   for( i = 0; i < consdata->nquadvars; ++i )
      assert(SCIPvarIsActive(consdata->quadvarterms[i].var));
#endif

   if( !have_change )
      return SCIP_OKAY;

   /* some quadratic variable may have been replaced by an already existing linear variable
    * in this case, we want the linear variable to be removed, which happens in mergeAndCleanLinearVars
    */ 
   consdata->linvarsmerged = FALSE;

   SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanQuadVarTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanLinearVars(scip, cons) );

#ifndef NDEBUG
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      assert(consdata->bilinterms[i].var1 != consdata->bilinterms[i].var2);
      assert(consdata->bilinterms[i].coef != 0.0);
      assert(SCIPvarCompare(consdata->bilinterms[i].var1, consdata->bilinterms[i].var2) < 0);
   }
#endif

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int        nquadvars;     /* number of variables in quadratic terms */
   SCIP_VAR** quadvars;      /* variables in quadratic terms */
   int        nquadelems;    /* number of quadratic elements (square and bilinear terms) */
   SCIP_QUADELEM* quadelems; /* quadratic elements (square and bilinear terms) */
   int        nquadlinterms; /* number of linear terms using variables that are in quadratic terms */
   SCIP_VAR** quadlinvars;   /* variables of linear terms using variables that are in quadratic terms */
   SCIP_Real* quadlincoefs;  /* coefficients of linear terms using variables that are in quadratic terms */
   int i;
   int idx1;
   int idx2;
   int lincnt;
   int elcnt;
   SCIP_VAR* lastvar;
   int lastvaridx;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   nquadvars = consdata->nquadvars;
   nquadelems = consdata->nbilinterms;
   nquadlinterms = 0;
   for( i = 0; i < nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].sqrcoef != 0.0 )
         ++nquadelems;
      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) )
         ++nquadlinterms;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars,  nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadlinvars,  nquadlinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadlincoefs, nquadlinterms) );

   lincnt = 0;
   elcnt = 0;
   for( i = 0; i < nquadvars; ++i )
   {
      quadvars[i] = consdata->quadvarterms[i].var;

      if( consdata->quadvarterms[i].sqrcoef != 0.0 )
      {
         assert(elcnt < nquadelems);
         quadelems[elcnt].idx1 = i;
         quadelems[elcnt].idx2 = i;
         quadelems[elcnt].coef = consdata->quadvarterms[i].sqrcoef;
         ++elcnt;
      }

      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) )
      {
         assert(lincnt < nquadlinterms);
         quadlinvars [lincnt] = consdata->quadvarterms[i].var;
         quadlincoefs[lincnt] = consdata->quadvarterms[i].lincoef;
         ++lincnt;
      }
   }
   assert(lincnt == nquadlinterms);

   /* bilinear terms are sorted first by first variable, then by second variable
    * thus, it makes sense to remember the index of the previous first variable for the case a series of bilinear terms with the same first var appears */
   lastvar = NULL;
   lastvaridx = -1;
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      if( lastvar == consdata->bilinterms[i].var1 )
      {
         assert(lastvaridx >= 0);
         assert(consdata->quadvarterms[lastvaridx].var == consdata->bilinterms[i].var1);
      }
      else
      {
         lastvar = consdata->bilinterms[i].var1;
         SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, lastvar, &lastvaridx) );
      }
      idx1 = lastvaridx;

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[i].var2, &idx2) );

      assert(elcnt < nquadelems);
      quadelems[elcnt].idx1 = MIN(idx1, idx2);
      quadelems[elcnt].idx2 = MAX(idx1, idx2);
      quadelems[elcnt].coef = consdata->bilinterms[i].coef;
      ++elcnt;
   }
   assert(elcnt == nquadelems);

   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
         consdata->nlinvars, consdata->linvars, consdata->lincoefs,
         nquadvars, quadvars, nquadelems, quadelems,
         NULL, consdata->lhs, consdata->rhs) );

   SCIP_CALL( SCIPaddLinearCoefsToNlRow(scip, consdata->nlrow, nquadlinterms, quadlinvars, quadlincoefs) );

   SCIPfreeBufferArray(scip, &quadvars);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &quadlinvars);
   SCIPfreeBufferArray(scip, &quadlincoefs);

   return SCIP_OKAY;
}

static
SCIP_RETCODE presolveSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RESULT*          result,             /**< to store result of solve: cutoff, success, or do-not-find */
   SCIP_Bool*            redundant,          /**< to store whether constraint is redundant now (should be deleted) */
   int*                  naggrvars           /**< counter on number of variable aggregations */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(redundant != NULL);

   *result = SCIP_DIDNOTFIND;
   *redundant = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if constraint is an equality with two variables, at least one of them binary,
    * and linear after fixing the binary, then we can aggregate the variables */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) && consdata->nlinvars == 0 && consdata->nquadvars == 2 &&
      ((SCIPvarIsBinary(consdata->quadvarterms[0].var) && consdata->quadvarterms[1].sqrcoef == 0.0) ||
       (SCIPvarIsBinary(consdata->quadvarterms[1].var) && consdata->quadvarterms[0].sqrcoef == 0.0)) )
   {
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;
      SCIP_Real a;
      SCIP_Real b;
      SCIP_Real c;
      SCIP_VAR* x;
      SCIP_VAR* y;
      int binvaridx;

      /* constraint is a*(x+x^2) + b*y + c*x*y = rhs, with x binary variable
       * x = 0 -> b*y == rhs
       * x = 1 -> (b+c)*y == rhs - a
       *
       * if b != 0 and b+c != 0, then y = (rhs-a)/(b+c) * x + rhs/b * (1-x) = ((rhs-a)/(b+c) - rhs/b) * x + rhs/b
       */

      binvaridx = (SCIPvarIsBinary(consdata->quadvarterms[0].var) && consdata->quadvarterms[1].sqrcoef == 0.0) ? 0 : 1;

      x = consdata->quadvarterms[binvaridx].var;
      a = consdata->quadvarterms[binvaridx].sqrcoef + consdata->quadvarterms[binvaridx].lincoef;

      y = consdata->quadvarterms[1-binvaridx].var;
      b = consdata->quadvarterms[1-binvaridx].lincoef;

      assert(consdata->nbilinterms <= 1);  /* should actually be 1, since constraint is otherwise linear */
      c = (consdata->nbilinterms == 1) ? consdata->bilinterms[0].coef : 0.0;

      if( !SCIPisZero(scip, b) && !SCIPisZero(scip, b+c) )
      {
         SCIPdebugMessage("<%s> = 0 -> %g*<%s> = %g  and  <%s> = 1 -> %g*<%s> = %g\n", SCIPvarGetName(x), b, SCIPvarGetName(y), consdata->rhs, SCIPvarGetName(x), b+c, SCIPvarGetName(y), consdata->rhs - a);
         SCIPdebugMessage("=> attempt aggregation <%s> = %g*<%s> + %g\n", SCIPvarGetName(y), (consdata->rhs-a)/(b+c) - consdata->rhs/b, SCIPvarGetName(x), consdata->rhs/b);

         SCIP_CALL( SCIPaggregateVars(scip, x, y, (consdata->rhs-a)/(b+c) - consdata->rhs/b, -1.0, -consdata->rhs/b, &infeasible, redundant, &aggregated) );
         if( infeasible )
            *result = SCIP_CUTOFF;
         else if( *redundant || aggregated )
         {
            /* aggregated (or were already aggregated), so constraint is now redundant */
            *result = SCIP_SUCCESS;
            *redundant = TRUE;

            if( aggregated )
               ++*naggrvars;
         }
      }

      /* @todo if b is 0 or b+c is 0, or lhs != rhs, then could replace by varbound constraint */
   }

   return SCIP_OKAY;
}


/** reformulates products of binary variables as AND constraint
 *  For a product x*y, with x and y binary variables, the product is replaced by a new auxiliary variable z and the constraint z = {x and y} is added.
 */
static
SCIP_RETCODE presolveTryAddAND(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  naddconss           /**< buffer where to add the number of AND constraints added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   char               name[SCIP_MAXSTRLEN];
   SCIP_VAR*          vars[2];
   SCIP_VAR*          auxvar;
   SCIP_CONS*         andcons;
   int                i;
   int                ntodelete;
   int*               todelete;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if user does not like AND very much, then return */
   if( conshdlrdata->empathy4and < 2 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nbilinterms == 0 )
      return SCIP_OKAY;

   /* get array to store indices of bilinear terms that shall be deleted */
   SCIP_CALL( SCIPallocBufferArray(scip, &todelete, consdata->nbilinterms) );
   ntodelete = 0;

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      vars[0] = consdata->bilinterms[i].var1;
      if( !SCIPvarIsBinary(vars[0]) )
         continue;

      vars[1] = consdata->bilinterms[i].var2;
      if( !SCIPvarIsBinary(vars[1]) )
         continue;

      /* create auxiliary variable */ 
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s*%s_%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]), SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, 
            TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );
#ifdef SCIP_DEBUG_SOLUTION
      {
         SCIP_Real var0val;
         SCIP_Real var1val;
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[0], &var0val) );
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[1], &var1val) );
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, var0val * var1val) );
      }
#endif

      /* create and constraint auxvar = x and y */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%sAND%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));
      SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, auxvar, 2, vars,
            SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
            SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, andcons) );
      SCIPdebugMessage("added AND constraint: ");
      SCIPdebugPrintCons(scip, andcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &andcons) );
      ++*naddconss;

      /* add bilincoef * auxvar to linear terms */
      SCIP_CALL( addLinearCoef(scip, cons, auxvar, consdata->bilinterms[i].coef) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

      /* remember that we have to delete this bilinear term */
      assert(ntodelete < consdata->nbilinterms);
      todelete[ntodelete++] = i;
   }

   /* remove bilinear terms that have been replaced */
   SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );
   SCIPfreeBufferArray(scip, &todelete);

   return SCIP_OKAY;
}

/** gets bounds of variable y if x takes a certain value
 * checks whether x = xval has implications on y
 */
static
void getImpliedBounds(
   SCIP_VAR*             x,                  /**< variable which implications to check */
   SCIP_Bool             xval,               /**< value of x to check for (TRUE for 1, FALSE for 0) */
   SCIP_VAR*             y,                  /**< variable to check if bounds can be reduced */
   SCIP_INTERVAL*        resultant           /**< buffer to store bounds on y */
   )
{
   int nimpls;
   int nbinimpls;
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;
   SCIP_Real* implbounds;
   int pos;

   assert(x != NULL);
   assert(y != NULL);
   assert(resultant != NULL);

   SCIPintervalSetBounds(resultant, MIN(SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y)), MAX(SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y)));  /*lint !e666 */

   if( !SCIPvarIsBinary(x) || !SCIPvarIsActive(x) )
      return;

   /* analyze implications for x = xval */
   nimpls = SCIPvarGetNImpls(x, xval);
   if( nimpls == 0 )
      return;

   nbinimpls  = SCIPvarGetNBinImpls (x, xval);
   implvars   = SCIPvarGetImplVars  (x, xval);
   impltypes  = SCIPvarGetImplTypes (x, xval);
   implbounds = SCIPvarGetImplBounds(x, xval);

   assert(implvars != NULL);
   assert(impltypes != NULL);
   assert(implbounds != NULL);

   if( SCIPvarGetType(y) == SCIP_VARTYPE_BINARY )
   {
      if( !SCIPsortedvecFindPtr((void**)implvars, SCIPvarComp, (void*)y, nbinimpls, &pos) )
         return;
   }
   else if( nbinimpls < nimpls )
   {
      if( !SCIPsortedvecFindPtr((void**)&implvars[nbinimpls], SCIPvarComp, (void*)y, nimpls - nbinimpls, &pos) )
         return;
   }
   else
      return;

   /* if there are several implications on y, go to the first one */
   while( pos > 0 && implvars[pos-1] == y )
      --pos;

   /* update implied lower and upper bounds on y
    * but make sure that resultant will not be empty, due to tolerances
    */
   while( pos < nimpls && implvars[pos] == y )
   {
      if( impltypes[pos] == SCIP_BOUNDTYPE_LOWER )
         resultant->inf = MAX(resultant->inf, MIN(resultant->sup, implbounds[pos]));
      else
         resultant->sup = MIN(resultant->sup, MAX(resultant->inf, implbounds[pos]));
      ++pos;
   }

   assert(!SCIPintervalIsEmpty(*resultant));
}

/** Reformulates products of binary times bounded continuous variables as system of linear inequalities (plus auxiliary variable).
 * 
 *  For a product x*y, with y a binary variable and x a continous variable with finite bounds,
 *  an auxiliary variable z and the inequalities \f$ x^L * y \leq z \leq x^U * y \f$ and \f$ x - (1-y)*x^U \leq z \leq x - (1-y)*x^L \f$ are added.
 * 
 *  If x is a linear term consisting of more than one variable, it is split up in groups of linear terms of length at most maxnrvar.
 *  For each product of linear term of length at most maxnrvar with y, an auxiliary z and linear inequalities are added.
 * 
 *  If y is a binary variable, the AND constraint \f$ z = x \wedge y \f$ may be added instead of linear constraints.
 */
static
SCIP_RETCODE presolveTryAddLinearReform(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  naddconss           /**< buffer where to add the number of auxiliary constraints added */
   )
{  /*lint --e{666} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR**         xvars;
   SCIP_Real*         xcoef;
   SCIP_INTERVAL      xbndszero;
   SCIP_INTERVAL      xbndsone;
   SCIP_INTERVAL      act0;
   SCIP_INTERVAL      act1;
   int                nxvars;
   SCIP_VAR*          y;
   SCIP_VAR*          bvar;
   char               name[SCIP_MAXSTRLEN];
   int                nbilinterms;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         auxcons;
   int                i;
   int                j;
   int                k;
   int                bilinidx;
   SCIP_Real          bilincoef;
   SCIP_Real          mincoef;
   SCIP_Real          maxcoef;
   int*               todelete;
   int                ntodelete;
   int                maxnrvar;
   SCIP_Bool          integral;
   SCIP_Longint       gcd;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxnrvar = conshdlrdata->replacebinaryprodlength;
   if( maxnrvar == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   xvars = NULL;
   xcoef = NULL;
   todelete = NULL;
   gcd = 0;

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      y = consdata->quadvarterms[i].var;
      if( !SCIPvarIsBinary(y) )
         continue;

      nbilinterms = consdata->quadvarterms[i].nadjbilin;
      if( nbilinterms == 0 )
         continue;

      SCIP_CALL( SCIPreallocBufferArray(scip, &xvars, MIN(maxnrvar, nbilinterms)+2) ); /* add 2 for later use when creating linear constraints */
      SCIP_CALL( SCIPreallocBufferArray(scip, &xcoef, MIN(maxnrvar, nbilinterms)+2) );

      /* alloc array to store indices of bilinear terms that shall be deleted */
      SCIP_CALL( SCIPreallocBufferArray(scip, &todelete, nbilinterms) );
      ntodelete = 0;

      /* setup a list of bounded variables x_i with coefficients a_i that are multiplied with binary y: y*(sum_i a_i*x_i)
       * and compute range of sum_i a_i*x_i for the cases y = 0 and y = 1
       * we may need several rounds of maxnrvar < nbilinterms
       */
      j = 0;
      do
      {
         nxvars = 0;
         SCIPintervalSet(&xbndszero, 0.0);
         SCIPintervalSet(&xbndsone, 0.0);

         mincoef = SCIPinfinity(scip);
         maxcoef = 0.0;
         integral = TRUE;

         /* collect at most maxnrvar variables for x term */
         for( ; j < nbilinterms && nxvars < maxnrvar; ++j )
         {
            bilinidx = consdata->quadvarterms[i].adjbilin[j];
            assert(bilinidx >= 0);
            assert(bilinidx < consdata->nbilinterms);

            bvar = consdata->bilinterms[bilinidx].var1;
            if( bvar == y )
               bvar = consdata->bilinterms[bilinidx].var2;
            assert(bvar != y);

            /* skip products with unbounded variables */
            if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bvar)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bvar)) )
               continue;

            bilincoef = consdata->bilinterms[bilinidx].coef;
            assert(bilincoef != 0.0);

            /* get activity of bilincoef * x if y = 0 */
            getImpliedBounds(y, FALSE, bvar, &act0);
            SCIPintervalMulScalar(SCIPinfinity(scip), &act0, act0, bilincoef);

            /* get activity of bilincoef * x if y = 1 */
            getImpliedBounds(y,  TRUE, bvar, &act1);
            SCIPintervalMulScalar(SCIPinfinity(scip), &act1, act1, bilincoef);

            /* skip products that give rise to very large coefficients (big big-M's)
             * we just reuse cutmaxrange as threshold, which is similarly motivated
             */
            if( SCIPintervalGetInf(act0) <= -conshdlrdata->cutmaxrange || SCIPintervalGetSup(act0) >= conshdlrdata->cutmaxrange )
            {
               SCIPdebugMessage("skip reform of %g<%s><%s> due to huge activity [%g,%g] for <%s> = 0.0\n",
                  bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar), SCIPintervalGetInf(act0), SCIPintervalGetSup(act0), SCIPvarGetName(y));
               continue;
            }
            if( SCIPintervalGetInf(act1) <= -conshdlrdata->cutmaxrange || SCIPintervalGetSup(act1) >= conshdlrdata->cutmaxrange )
            {
               SCIPdebugMessage("skip reform of %g<%s><%s> due to huge activity [%g,%g] for <%s> = 1.0\n",
                  bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar), SCIPintervalGetInf(act1), SCIPintervalGetSup(act1), SCIPvarGetName(y));
               continue;
            }

            /* add bvar to x term */  
            xvars[nxvars] = bvar;
            xcoef[nxvars] = bilincoef;
            ++nxvars;

            /* update bounds on x term */
            SCIPintervalAdd(SCIPinfinity(scip), &xbndszero, xbndszero, act0);
            SCIPintervalAdd(SCIPinfinity(scip), &xbndsone,  xbndsone,  act1);

            if( REALABS(bilincoef) < mincoef )
               mincoef = ABS(bilincoef);
            if( REALABS(bilincoef) > maxcoef )
               maxcoef = ABS(bilincoef);

            /* update whether all coefficients will be integral and if so, compute their gcd */
            integral &= (SCIPvarGetType(bvar) < SCIP_VARTYPE_CONTINUOUS) && SCIPisIntegral(scip, bilincoef);  /*lint !e514 */
            if( integral )
            {
               if( nxvars == 1 )
                  gcd = (SCIP_Longint)SCIPround(scip, REALABS(bilincoef));
               else
                  gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)SCIPround(scip, REALABS(bilincoef)));
            }

            /* remember that we have to remove this bilinear term later */
            assert(ntodelete < nbilinterms);
            todelete[ntodelete++] = bilinidx;
         }

         if( nxvars == 0 ) /* all (remaining) x_j seem to be unbounded */
            break;

         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbndszero)));
         assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(xbndszero)));
         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbndsone)));
         assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(xbndsone)));

#ifdef SCIP_DEBUG
         if( SCIPintervalGetInf(xbndszero) != SCIPintervalGetInf(xbndsone) || /*lint !e777*/
            +SCIPintervalGetSup(xbndszero) != SCIPintervalGetSup(xbndsone) ) /*lint !e777*/
         {
            SCIPdebugMessage("got different bounds for y = 0: [%g, %g] and y = 1: [%g, %g]\n", xbndszero.inf, xbndszero.sup, xbndsone.inf, xbndsone.sup);
         }
#endif

         if( nxvars == 1 && conshdlrdata->empathy4and >= 1 && SCIPvarIsBinary(xvars[0]) )
         {
            /* product of two binary variables, replace by auxvar and AND constraint */
            /* add auxiliary variable z */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s*%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT,
                  TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

#ifdef SCIP_DEBUG_SOLUTION
            {
               SCIP_Real var0val;
               SCIP_Real var1val;
               SCIP_CALL( SCIPdebugGetSolVal(scip, xvars[0], &var0val) );
               SCIP_CALL( SCIPdebugGetSolVal(scip, y, &var1val) );
               SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, var0val * var1val) );
            }
#endif

            /* add constraint z = x and y */
            xvars[1] = y;
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%sAND%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsAnd(scip, &auxcons, name, auxvar, 2, xvars,
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added AND constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add linear term coef*auxvar */
            SCIP_CALL( addLinearCoef(scip, cons, auxvar, xcoef[0]) );

            /* forget about auxvar */
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         else
         {
            /* product of binary variable with more than one binary or with continuous variables or with binary and user did not like AND -> replace by auxvar and linear constraints */
            SCIP_Real scale;

            /* scale auxiliary constraint by some nice value,
             * if all coefficients are integral, take a value that preserves integrality (-> gcd), so we can make the auxiliary variable impl. integer
             */
            if( integral )
            {
               scale = (SCIP_Real)gcd;
               assert(scale >= 1.0);
            }
            else if( nxvars == 1 )
            {
               /* scaling by the only coefficient gives auxiliary variable = x * y, which thus will be implicit integral provided y is not continuous */
               assert(mincoef == maxcoef);  /*lint !e777 */
               scale = mincoef;
               integral = SCIPvarGetType(xvars[0]) < SCIP_VARTYPE_CONTINUOUS;
            }
            else
            {
               scale = 1.0;
               if( maxcoef < 0.5 )
                  scale = maxcoef;
               if( mincoef > 2.0 )
                  scale = mincoef;
               if( scale != 1.0 )
                  scale = SCIPselectSimpleValue(scale / 2.0, 1.5 * scale, MAXDNOM);
            }
            assert(scale > 0.0);
            assert(!SCIPisInfinity(scip, scale));

            /* if x-term is always negative for y = 1, negate scale so we get a positive auxiliary variable; maybe this is better sometimes? */
            if( !SCIPisPositive(scip, SCIPintervalGetSup(xbndsone)) )
               scale = -scale;

            SCIPdebugMessage("binary reformulation using scale %g, nxvars = %d, integral = %u\n", scale, nxvars, integral);
            if( scale != 1.0 )
            {
               SCIPintervalDivScalar(SCIPinfinity(scip), &xbndszero, xbndszero, scale);
               SCIPintervalDivScalar(SCIPinfinity(scip), &xbndsone,  xbndsone, scale);
               for( k = 0; k < nxvars; ++k )
                  xcoef[k] /= scale;
            }

            /* add auxiliary variable z */
            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s*%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s*%s*more_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, MIN(0., SCIPintervalGetInf(xbndsone)), MAX(0., SCIPintervalGetSup(xbndsone)),
                  0.0, integral ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

            /* compute value of auxvar in debug solution */
#ifdef SCIP_DEBUG_SOLUTION
            {
               SCIP_Real debugval;
               SCIP_Real varval;

               SCIP_CALL( SCIPdebugGetSolVal(scip, y, &varval) );
               if( SCIPisZero(scip, varval) )
               {
                  SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, 0.0) );
               }
               else
               {
                  assert(SCIPisEQ(scip, varval, 1.0));

                  debugval = 0.0;
                  for( k = 0; k < nxvars; ++k )
                  {
                     SCIP_CALL( SCIPdebugGetSolVal(scip, xvars[k], &varval) );
                     debugval += xcoef[k] * varval;
                  }
                  SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
               }
            }
#endif

            /* add auxiliary constraints
             * it seems to be advantageous to make the varbound constraints initial and the linear constraints not initial
             * maybe because it is more likely that a binary variable takes value 0 instead of 1, and thus the varbound constraints
             * are more often active, compared to the linear constraints added below
             * also, the varbound constraints are more sparse than the linear cons
             */
            if( SCIPisNegative(scip, SCIPintervalGetInf(xbndsone)) )
            {
               /* add 0 <= z - xbndsone.inf * y constraint (as varbound constraint) */
               if( nxvars == 1 )
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_1", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               else
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_1", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetInf(xbndsone), 0.0, SCIPinfinity(scip),
                     SCIPconsIsInitial(cons) /*&& conshdlrdata->binreforminitial*/,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
               SCIPdebugPrintCons(scip, auxcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*naddconss;
            }
            if( SCIPisPositive(scip, SCIPintervalGetSup(xbndsone)) )
            {
               /* add z - xbndsone.sup * y <= 0 constraint (as varbound constraint) */
               if( nxvars == 1 )
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_2", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               else
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_2", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetSup(xbndsone), -SCIPinfinity(scip), 0.0,
                     SCIPconsIsInitial(cons) /*&& conshdlrdata->binreforminitial*/,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
               SCIPdebugPrintCons(scip, auxcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*naddconss;
            }

            /* add xbndszero.inf <= sum_i a_i*x_i + xbndszero.inf * y - z constraint */
            xvars[nxvars]   = y;
            xvars[nxvars+1] = auxvar;
            xcoef[nxvars]   = SCIPintervalGetInf(xbndszero);
            xcoef[nxvars+1] = -1;

            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_3", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_3", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, SCIPintervalGetInf(xbndszero), SCIPinfinity(scip),
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add sum_i a_i*x_i + xbndszero.sup * y - z <= xbndszero.sup constraint */
            xcoef[nxvars] = SCIPintervalGetSup(xbndszero);

            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_4", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_4", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, -SCIPinfinity(scip), SCIPintervalGetSup(xbndszero),
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add linear term scale*auxvar to this constraint */
            SCIP_CALL( addLinearCoef(scip, cons, auxvar, scale) );

            /* forget about auxvar */
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
      }
      while( j < nbilinterms );

      /* remove bilinear terms that have been replaced */
      SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );
   }

   SCIPfreeBufferArrayNull(scip, &xvars);
   SCIPfreeBufferArrayNull(scip, &xcoef);
   SCIPfreeBufferArrayNull(scip, &todelete);

   return SCIP_OKAY;
}

/** tries to automatically convert a quadratic constraint (or a part of it) into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss           /**< buffer to increase with number of additional constraints created during upgrade */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real lincoef;
   SCIP_Real quadcoef;
   SCIP_Real lb;
   SCIP_Real ub;
   int nbinlin;
   int nbinquad;
   int nintlin;
   int nintquad;
   int nimpllin;
   int nimplquad;
   int ncontlin;
   int ncontquad;
   SCIP_Bool integral;
   int i;
   int j;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsModifiable(cons));
   assert(upgraded   != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss  != NULL);

   *upgraded = FALSE;

   nupgdconss_ = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if there are no upgrade methods, we can also stop */
   if( conshdlrdata->nquadconsupgrades == 0 )
      return SCIP_OKAY;

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* calculate some statistics on quadratic constraint */
   nbinlin   = 0;
   nbinquad  = 0;
   nintlin   = 0;
   nintquad  = 0;
   nimpllin  = 0;
   nimplquad = 0;
   ncontlin  = 0;
   ncontquad = 0;
   integral  = TRUE;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      var = consdata->linvars[i];
      lincoef = consdata->lincoefs[i];
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, lincoef));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nbinlin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nintlin++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nimpllin++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisRelEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb);
         ncontlin++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvarterms[i].var;
      lincoef  = consdata->quadvarterms[i].lincoef;
      quadcoef = consdata->quadvarterms[i].sqrcoef;
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nbinquad++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nintquad++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nimplquad++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisRelEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb + quadcoef * lb * lb);
         ncontquad++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }

   if( integral )
   {
      for( i = 0; i < consdata->nbilinterms && integral; ++i )
      {
         if( SCIPvarGetType(consdata->bilinterms[i].var1) < SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(consdata->bilinterms[i].var2) < SCIP_VARTYPE_CONTINUOUS )
            integral = integral && SCIPisIntegral(scip, consdata->bilinterms[i].coef);
         else
            integral = FALSE;
      }
   }

   /* call the upgrading methods */

   SCIPdebugMessage("upgrading quadratic constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nquadconsupgrades);
   SCIPdebugMessage(" binlin=%d binquad=%d intlin=%d intquad=%d impllin=%d implquad=%d contlin=%d contquad=%d integral=%u\n",
      nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral);
   SCIPdebugPrintCons(scip, cons, NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      if( !conshdlrdata->quadconsupgrades[i]->active )
         continue;

      SCIP_CALL( conshdlrdata->quadconsupgrades[i]->quadconsupgd(scip, cons,
            nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral,
            &nupgdconss_, upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->quadconsupgrades[i]->quadconsupgd(scip, cons,
               nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral,
               &nupgdconss_, upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIPdebugMessage(" -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugPrintf("\t");
            SCIPdebugPrintCons(scip, upgdconss[j], NULL);

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMessage("delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** helper function for presolveDisaggregate */
static
SCIP_RETCODE presolveDisaggregateMarkComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   quadvaridx,         /**< index of quadratic variable to mark */
   SCIP_HASHMAP*         var2component,      /**< variables to components mapping */
   int                   componentnr         /**< the component number to mark to */ 
   )
{
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_VAR* othervar;
   int othervaridx;
   int i;

   assert(consdata != NULL);
   assert(quadvaridx >= 0);
   assert(quadvaridx < consdata->nquadvars);
   assert(var2component != NULL);
   assert(componentnr >= 0);

   quadvarterm = &consdata->quadvarterms[quadvaridx];

   if( SCIPhashmapExists(var2component, quadvarterm->var) )
   {
      /* if we saw the variable before, then it should have the same component number */
      assert((int)(size_t)SCIPhashmapGetImage(var2component, quadvarterm->var) == componentnr);
      return SCIP_OKAY;
   }

   /* assign component number to variable */
   SCIP_CALL( SCIPhashmapInsert(var2component, quadvarterm->var, (void*)(size_t)componentnr) );

   /* assign same component number to all variables this variable is multiplied with */
   for( i = 0; i < quadvarterm->nadjbilin; ++i )
   {
      othervar = consdata->bilinterms[quadvarterm->adjbilin[i]].var1 == quadvarterm->var ?
         consdata->bilinterms[quadvarterm->adjbilin[i]].var2 : consdata->bilinterms[quadvarterm->adjbilin[i]].var1;
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, othervar, &othervaridx) );
      assert(othervaridx >= 0);
      SCIP_CALL( presolveDisaggregateMarkComponent(scip, consdata, othervaridx, var2component, componentnr) );
   }

   return SCIP_OKAY;
}

/** for quadratic constraints that consists of a sum of quadratic terms, disaggregates the sum into a set of constraints by introducing auxiliary variables */
static
SCIP_RETCODE presolveDisaggregate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   int*                  naddconss           /**< pointer to counter of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* var2component;
   int ncomponents;
   int i;
   int comp;
   SCIP_CONS** auxconss;
   SCIP_VAR** auxvars;
   SCIP_Real* auxcoefs;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* make sure there are no quadratic variables without coefficients */
   SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanQuadVarTerms(scip, cons) );

   if( consdata->nquadvars <= 1 )
      return SCIP_OKAY;

   /* sort quadratic variable terms here, so we can later search in it without reordering the array */
   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   /* check how many quadratic terms with non-overlapping variables we have
    * in other words, the number of components in the sparsity graph of the quadratic term matrix */
   ncomponents = 0;
   SCIP_CALL( SCIPhashmapCreate(&var2component, SCIPblkmem(scip), SCIPcalcHashtableSize(consdata->nquadvars)) );
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* if variable was marked already, skip it */
      if( SCIPhashmapExists(var2component, (void*)consdata->quadvarterms[i].var) )
         continue;

      SCIP_CALL( presolveDisaggregateMarkComponent(scip, consdata, i, var2component, ncomponents) );
      ++ncomponents;
   }
   assert(ncomponents >= 1);

   /* if there is only one component, we cannot disaggregate
    * @todo we could still split the constraint into several while keeping the number of variables sharing several constraints as small as possible
    */
   if( ncomponents == 1 )
   {
      SCIPhashmapFree(&var2component);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &auxconss, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxvars,  ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxcoefs, ncomponents) );

   /* create auxiliary variables and empty constraints for each component */
   for( comp = 0; comp < ncomponents; ++comp )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp%d", SCIPconsGetName(cons), comp);

      SCIP_CALL( SCIPcreateVar(scip, &auxvars[comp], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPcreateConsQuadratic2(scip, &auxconss[comp], name, 0, NULL, NULL, 0, NULL, 0, NULL,
            (SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : 0.0),
            (SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : 0.0),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

      auxcoefs[comp] = SCIPinfinity(scip);
   }

   /* add quadratic variables to each component constraint
    * delete adjacency information */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      comp = (int)(size_t) SCIPhashmapGetImage(var2component, consdata->quadvarterms[i].var);
      assert(comp >= 0);
      assert(comp < ncomponents);

      /* add variable term to corresponding constraint */
      SCIP_CALL( SCIPaddQuadVarQuadratic(scip, auxconss[comp], consdata->quadvarterms[i].var, consdata->quadvarterms[i].lincoef, consdata->quadvarterms[i].sqrcoef) );

      /* reduce coefficient of aux variable */
      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) && ABS(consdata->quadvarterms[i].lincoef) < auxcoefs[comp] )
         auxcoefs[comp] = REALABS(consdata->quadvarterms[i].lincoef);
      if( !SCIPisZero(scip, consdata->quadvarterms[i].sqrcoef) && ABS(consdata->quadvarterms[i].sqrcoef) < auxcoefs[comp] )
         auxcoefs[comp] = REALABS(consdata->quadvarterms[i].sqrcoef);

      SCIPfreeBlockMemoryArray(scip, &consdata->quadvarterms[i].adjbilin, consdata->quadvarterms[i].adjbilinsize);
      consdata->quadvarterms[i].nadjbilin = 0;
      consdata->quadvarterms[i].adjbilinsize = 0;
   }

   /* add bilinear terms to each component constraint */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      comp = (int)(size_t) SCIPhashmapGetImage(var2component, consdata->bilinterms[i].var1);
      assert(comp == (int)(size_t) SCIPhashmapGetImage(var2component, consdata->bilinterms[i].var2));
      assert(!SCIPisZero(scip, consdata->bilinterms[i].coef));

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, auxconss[comp], 
            consdata->bilinterms[i].var1, consdata->bilinterms[i].var2, consdata->bilinterms[i].coef) );

      if( ABS(consdata->bilinterms[i].coef) < auxcoefs[comp] )
         auxcoefs[comp] = ABS(consdata->bilinterms[i].coef);
   }

   /* forget about bilinear terms in cons */
   SCIPfreeBlockMemoryArray(scip, &consdata->bilinterms, consdata->bilintermssize);
   consdata->nbilinterms = 0;
   consdata->bilintermssize = 0;

   /* remove quadratic variable terms from cons */
   for( i = consdata->nquadvars - 1; i >= 0; --i )
   {
      SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
   }
   assert(consdata->nquadvars == 0);

   /* add auxiliary variables to auxiliary constraints
    * add aux vars and constraints to SCIP 
    * add aux vars to this constraint
    * @todo compute debug solution values and set for auxvars
    */
   SCIPdebugMessage("add %d constraints for disaggregation of quadratic constraint <%s>\n", ncomponents, SCIPconsGetName(cons));
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + ncomponents) );
   for( comp = 0; comp < ncomponents; ++comp )
   {
      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, auxconss[comp], auxvars[comp], -auxcoefs[comp]) );

      SCIP_CALL( SCIPaddVar(scip, auxvars[comp]) );

      SCIP_CALL( SCIPaddCons(scip, auxconss[comp]) );
      SCIPdebugPrintCons(scip, auxconss[comp], NULL);

      SCIP_CALL( addLinearCoef(scip, cons, auxvars[comp], auxcoefs[comp]) );

      SCIP_CALL( SCIPreleaseCons(scip, &auxconss[comp]) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvars[comp]) );
   }
   *naddconss += ncomponents;

   SCIPdebugPrintCons(scip, cons, NULL);

   SCIPfreeBufferArray(scip, &auxconss);
   SCIPfreeBufferArray(scip, &auxvars);
   SCIPfreeBufferArray(scip, &auxcoefs);
   SCIPhashmapFree(&var2component);

   return SCIP_OKAY;
}

#ifdef CHECKIMPLINBILINEAR
/** checks if there are bilinear terms x*y with a binary variable x and an implication x = {0,1} -> y = 0
 * in this case, the bilinear term can be removed (x=0 case) or replaced by y (x=1 case)
 */
static
SCIP_RETCODE presolveApplyImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int*                  nbilinremoved       /**< buffer to store number of removed bilinear terms */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_INTERVAL implbnds;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nbilinremoved != NULL);

   *nbilinremoved = 0;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("apply implications in <%s>\n", SCIPconsGetName(cons));

   /* sort quadvarterms in case we need to search */
   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      x = consdata->quadvarterms[i].var;
      assert(x != NULL);

      if( consdata->quadvarterms[i].nadjbilin == 0 )
         continue;

      if( !SCIPvarIsBinary(x) )
         continue;

      if( !SCIPvarIsActive(x) )
         continue;

      if( SCIPvarGetNImpls(x, TRUE) == 0 && SCIPvarGetNImpls(x, FALSE) == 0 )
         continue;

      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         k = consdata->quadvarterms[i].adjbilin[j];
         assert(k >= 0);
         assert(k < consdata->nbilinterms);

         if( consdata->bilinterms[k].coef == 0.0 )
            continue;

         y = consdata->bilinterms[k].var1 == x ? consdata->bilinterms[k].var2 : consdata->bilinterms[k].var1;
         assert(x != y);

         getImpliedBounds(x, TRUE, y, &implbnds);
         if( SCIPisZero(scip, implbnds.inf) && SCIPisZero(scip, implbnds.sup) )
         {
            /* if x = 1 implies y = 0, then we can remove the bilinear term x*y, since it is always 0
             * we only set the coefficient to 0.0 here and mark the bilinterms as not merged */
            SCIPdebugMessage("remove bilinear term %g<%s><%s> from <%s> due to implication\n", consdata->bilinterms[k].coef, SCIPvarGetName(x), SCIPvarGetName(y), SCIPconsGetName(cons));
            consdata->bilinterms[k].coef = 0.0;
            consdata->bilinmerged = FALSE;
            ++*nbilinremoved;
            continue;
         }

         getImpliedBounds(x, FALSE, y, &implbnds);
         if( SCIPisZero(scip, implbnds.inf) && SCIPisZero(scip, implbnds.sup) )
         {
            /* if x = 0 implies y = 0, then we can replace the bilinear term x*y by y
             * we only move the coefficient to the linear coef of y here and mark the bilinterms as not merged */
            SCIPdebugMessage("replace bilinear term %g<%s><%s> by %g<%s> in <%s> due to implication\n", consdata->bilinterms[k].coef, SCIPvarGetName(x), SCIPvarGetName(y), consdata->bilinterms[k].coef, SCIPvarGetName(y), SCIPconsGetName(cons));
            assert(consdata->quadvarssorted);
            SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, cons, y, consdata->bilinterms[k].coef) );
            consdata->bilinterms[k].coef = 0.0;
            consdata->bilinmerged = FALSE;
            ++*nbilinremoved;
         }
      }
   }

   if( *nbilinremoved > 0 )
   {
      SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );

      /* invalidate nonlinear row */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      consdata->ispropagated  = FALSE;
      consdata->ispresolved   = FALSE;
      consdata->iscurvchecked = FALSE;
   }

   consdata->isimpladded = FALSE;

   return SCIP_OKAY;
}
#endif

/** checks a quadratic constraint for convexity and/or concavity without checking multivariate functions */
static
void checkCurvatureEasy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_Bool*            determined,         /**< pointer to store whether the curvature could be determined */
   SCIP_Bool             checkmultivariate   /**< whether curvature will be checked later on for multivariate functions */
   )
{
   SCIP_CONSDATA* consdata;
   int nquadvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(determined != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nquadvars = consdata->nquadvars;
   *determined = TRUE;

   if( consdata->iscurvchecked )
      return;

   SCIPdebugMessage("Checking curvature of constraint <%s> without multivariate functions\n", SCIPconsGetName(cons));

   if( nquadvars == 1 )
   {
      assert(consdata->nbilinterms == 0);
      consdata->isconvex      = !SCIPisNegative(scip, consdata->quadvarterms[0].sqrcoef);
      consdata->isconcave     = !SCIPisPositive(scip, consdata->quadvarterms[0].sqrcoef);
      consdata->iscurvchecked = TRUE;
   }
   else if( nquadvars == 0 )
   {
      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;
      consdata->iscurvchecked = TRUE;
   }
   else if( consdata->nbilinterms == 0 )
   {
      int v;

      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;

      for( v = nquadvars - 1; v >= 0; --v )
      {
         consdata->isconvex  = consdata->isconvex  && !SCIPisNegative(scip, consdata->quadvarterms[v].sqrcoef);
         consdata->isconcave = consdata->isconcave && !SCIPisPositive(scip, consdata->quadvarterms[v].sqrcoef);
      }

      consdata->iscurvchecked = TRUE;
   }
   else if( !checkmultivariate )
   {
      consdata->isconvex  = FALSE;
      consdata->isconcave = FALSE;
      consdata->iscurvchecked = TRUE;
   }
   else
      *determined = FALSE;
}

/** checks a quadratic constraint for convexity and/or concavity */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_Bool             checkmultivariate   /**< whether curvature should also be checked for multivariate functions */
   )
{
   SCIP_CONSDATA* consdata;
   double*        matrix;
   SCIP_HASHMAP*  var2index;
   int            i;
   int            n;
   int            nn;
   int            row;
   int            col;
   double*        alleigval;
   SCIP_Bool      determined;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   n = consdata->nquadvars;

   if( consdata->iscurvchecked )
      return SCIP_OKAY;

   /* easy checks for curvature detection */
   checkCurvatureEasy(scip, cons, &determined, checkmultivariate);

   /* if curvature was already detected stop */
   if( determined )
      return SCIP_OKAY;

   SCIPdebugMessage("Checking curvature of constraint <%s> with multivariate functions\n", SCIPconsGetName(cons));

   if( n == 2 )
   {
      /* compute eigenvalues by hand */
      assert(consdata->nbilinterms == 1);
      consdata->isconvex =
         consdata->quadvarterms[0].sqrcoef >= 0 &&
         consdata->quadvarterms[1].sqrcoef >= 0 &&
         4 * consdata->quadvarterms[0].sqrcoef * consdata->quadvarterms[1].sqrcoef >= consdata->bilinterms[0].coef * consdata->bilinterms[0].coef;
      consdata->isconcave = 
         consdata->quadvarterms[0].sqrcoef <= 0 &&
         consdata->quadvarterms[1].sqrcoef <= 0 &&
         4 * consdata->quadvarterms[0].sqrcoef * consdata->quadvarterms[1].sqrcoef >= consdata->bilinterms[0].coef * consdata->bilinterms[0].coef;
      consdata->iscurvchecked = TRUE;
      return SCIP_OKAY;
   }

   /* lower triangular of quadratic term matrix */
   nn = n * n;
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, nn) );
   BMSclearMemoryArray(matrix, nn);

   consdata->isconvex  = TRUE;
   consdata->isconcave = TRUE;

   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * n)) );
   for( i = 0; i < n; ++i )
   {
      if( consdata->quadvarterms[i].nadjbilin > 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvarterms[i].var, (void*)(size_t)i) );
         matrix[i*n + i] = consdata->quadvarterms[i].sqrcoef;
      }
      /* nonzero elements on diagonal tell a lot about convexity/concavity */
      if( SCIPisNegative(scip, consdata->quadvarterms[i].sqrcoef) )
         consdata->isconvex  = FALSE;
      if( SCIPisPositive(scip, consdata->quadvarterms[i].sqrcoef) )
         consdata->isconcave = FALSE;
   }

   if( !consdata->isconvex && !consdata->isconcave )
   {
      SCIPfreeBufferArray(scip, &matrix);
      SCIPhashmapFree(&var2index);
      consdata->iscurvchecked = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPisIpoptAvailableIpopt() )
   {
      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var1));
         assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var2));
         row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var1);
         col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var2);
         if( row < col )
            matrix[row * n + col] = consdata->bilinterms[i].coef/2;
         else
            matrix[col * n + row] = consdata->bilinterms[i].coef/2;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
      /* @todo Can we compute only min and max eigen value?
       * @todo Can we estimate the numerical error? 
       * @todo Trying a cholesky factorization may be much faster. 
       */
      if( LapackDsyev(FALSE, n, matrix, alleigval) != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Failed to compute eigenvalues of quadratic coefficient matrix of constraint %s. Assuming matrix is indefinite.\n", SCIPconsGetName(cons));
         consdata->isconvex = FALSE;
         consdata->isconcave = FALSE;
      }
      else
      {
         /* deconvexification reformulates a stricly convex quadratic function in binaries such that it becomes not-strictly convex
          * by adding the -lambda*(x^2-x) terms for lambda the smallest eigenvalue of the matrix
          * the result is still a convex form "but less so" (ref. papers by Guignard et.al.), but with hopefully tighter value for the continuous relaxation
          */
#ifdef DECONVEXIFY
         SCIP_Bool allbinary;
         printf("cons <%s>[%g,%g] spectrum = [%g,%g]\n", SCIPconsGetName(cons), consdata->lhs, consdata->rhs, alleigval[0], alleigval[n-1]);
#endif
         consdata->isconvex  &= !SCIPisNegative(scip, alleigval[0]);   /*lint !e514*/
         consdata->isconcave &= !SCIPisPositive(scip, alleigval[n-1]); /*lint !e514*/
         consdata->iscurvchecked = TRUE;
#ifdef DECONVEXIFY
         for( i = 0; i < consdata->nquadvars; ++i )
            if( !SCIPvarIsBinary(consdata->quadvarterms[i].var) )
               break;
         allbinary = i == consdata->nquadvars;

         if( !SCIPisInfinity(scip, consdata->rhs) && alleigval[0] > 0.1 && allbinary )
         {
            printf("deconvexify cons <%s> by shifting hessian by %g\n", SCIPconsGetName(cons), alleigval[0]);
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               consdata->quadvarterms[i].sqrcoef -= alleigval[0];
               consdata->quadvarterms[i].lincoef += alleigval[0];
            }
         }

         if( !SCIPisInfinity(scip, consdata->lhs) && alleigval[n-1] < -0.1 && allbinary )
         {
            printf("deconcavify cons <%s> by shifting hessian by %g\n", SCIPconsGetName(cons), alleigval[n-1]);
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               consdata->quadvarterms[i].sqrcoef -= alleigval[n-1];
               consdata->quadvarterms[i].lincoef += alleigval[n-1];
            }
         }
#endif
      }

      SCIPfreeBufferArray(scip, &alleigval);
   }
   else
   {
      consdata->isconvex = FALSE;
      consdata->isconcave = FALSE;
      consdata->iscurvchecked = TRUE; /* set to TRUE since it does not help to repeat this procedure again and again (that will not bring Ipopt in) */
   }

   SCIPhashmapFree(&var2index);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** check whether indefinite constraint function is factorable and store corresponding coefficients */
static
SCIP_RETCODE checkFactorable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_BILINTERM* bilinterm;
   SCIP_CONSDATA* consdata;
   SCIP_Real* a;
   SCIP_Real* eigvals;
   SCIP_Real sigma1;
   SCIP_Real sigma2;
   int n;
   int i;
   int idx1;
   int idx2;
   int posidx;
   int negidx;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->factorleft == NULL);
   assert(consdata->factorright == NULL);

   /* we don't need this if there are no bilinear terms */
   if( consdata->nbilinterms == 0 )
      return SCIP_OKAY;

   /* write constraint as lhs <= linear + x'^T A x' <= rhs where x' = (x,1) and
    * A = ( Q     b/2 )
    *     ( b^T/2  0  )
    * compute an eigenvalue factorization of A and check if there are one positive and one negative eigenvalue
    * if so, then let sigma1^2 and -sigma2^2 be these eigenvalues and v1 and v2 be the first two rows of the inverse eigenvector matrix
    * thus, x'^T A x' = sigma1^2 (v1^T x')^2 - sigma2^2 (v2^T x')^2
    *                 = (sigma1 (v1^T x') - sigma2 (v2^T x')) * (sigma1 (v1^T x') + sigma2 (v2^T x'))
    * we then store sigma1 v1^T - sigma2 v2^T as left factor coef, and sigma1 v1^T + sigma2 v2^T as right factor coef
    */

   /* if we already know that there are only positive or only negative eigenvalues, then don't try */
   if( consdata->iscurvchecked && (consdata->isconvex || consdata->isconcave) )
      return SCIP_OKAY;

   n = consdata->nquadvars + 1;

   /* @todo handle case n=3 explicitly */

   /* skip too large matrices */
   if( n > 50 )
      return SCIP_OKAY;

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &a, n*n) );
   BMSclearMemoryArray(a, n*n);

   /* set lower triagonal entries of A corresponding to bilinear terms */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      bilinterm = &consdata->bilinterms[i];

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var1, &idx1) );
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var2, &idx2) );
      assert(idx1 >= 0);
      assert(idx2 >= 0);
      assert(idx1 != idx2);

      a[MIN(idx1,idx2) * n + MAX(idx1,idx2)] = bilinterm->coef / 2.0;
   }

   /* set lower triangle entries of A corresponding to square and linear terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      a[i*n + i]   = consdata->quadvarterms[i].sqrcoef;
      a[i*n + n-1] = consdata->quadvarterms[i].lincoef / 2.0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, n) );
   if( LapackDsyev(TRUE, n, a, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMessage("Failed to compute eigenvalues and eigenvectors of augmented quadratic form matrix for constraint <%s>.\n", SCIPconsGetName(cons));
      goto CLEANUP;
   }

   /* check if there is exactly one positive and one negative eigenvalue */
   posidx = -1;
   negidx = -1;
   for( i = 0; i < n; ++i )
   {
      if( SCIPisPositive(scip, eigvals[i]) )
      {
         if( posidx == -1 )
            posidx = i;
         else
            break;
      }
      else if( SCIPisNegative(scip, eigvals[i]) )
      {
         if( negidx == -1 )
            negidx = i;
         else
            break;
      }
   }
   if( i < n || posidx == -1 || negidx == -1 )
   {
      SCIPdebugMessage("Augmented quadratic form of constraint <%s> is not factorable.\n", SCIPconsGetName(cons));
      goto CLEANUP;
   }
   assert(SCIPisPositive(scip, eigvals[posidx]));
   assert(SCIPisNegative(scip, eigvals[negidx]));

   /* compute factorleft and factorright */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->factorleft,  consdata->nquadvars + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->factorright, consdata->nquadvars + 1) );

   /* eigenvectors are stored in a, inverse eigenvector matrix is transposed of a
    * it seems that v1 and v2 are at &a[posidx*n] and &a[negidx*n]
    */
   sigma1 = sqrt( eigvals[posidx]);
   sigma2 = sqrt(-eigvals[negidx]);
   for( i = 0; i < n; ++i )
   {
      consdata->factorleft[i]  = sigma1 * a[posidx * n + i] - sigma2 * a[negidx * n + i];
      consdata->factorright[i] = sigma1 * a[posidx * n + i] + sigma2 * a[negidx * n + i];
      if( SCIPisZero(scip, consdata->factorleft[i]) )
         consdata->factorleft[i] = 0.0;
      if( SCIPisZero(scip, consdata->factorright[i]) )
         consdata->factorright[i] = 0.0;
   }

   SCIPdebugMessage("constraint <%s> has factorable quadratic form: (%g", SCIPconsGetName(cons), consdata->factorleft[n-1]);
   for( i = 0; i < consdata->nquadvars; ++i )
      if( consdata->factorleft[i] != 0.0 )
         SCIPdebugPrintf(" %+g<%s>", consdata->factorleft[i], SCIPvarGetName(consdata->quadvarterms[i].var));
   SCIPdebugPrintf(") * (%g", consdata->factorright[n-1]);
   for( i = 0; i < consdata->nquadvars; ++i )
      if( consdata->factorright[i] != 0.0 )
         SCIPdebugPrintf(" %+g<%s>", consdata->factorright[i], SCIPvarGetName(consdata->quadvarterms[i].var));
   SCIPdebugPrintf(")\n");

#ifndef NDEBUG
   /* check whether factorleft * factorright^T is matrix of augmented quadratic form */
   BMSclearMemoryArray(a, n*n);

   /* set lower triagonal entries of A corresponding to bilinear terms */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      bilinterm = &consdata->bilinterms[i];

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var1, &idx1) );
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var2, &idx2) );
      assert(idx1 >= 0);
      assert(idx2 >= 0);
      assert(idx1 != idx2);

      a[MIN(idx1,idx2) * n + MAX(idx1,idx2)] = bilinterm->coef / 2.0;
   }

   /* set lower triagonal entries of A corresponding to square and linear terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      a[i*n + i]   = consdata->quadvarterms[i].sqrcoef;
      a[i*n + n-1] = consdata->quadvarterms[i].lincoef / 2.0;
   }

   /* compare matrix entries */
   for( i = 0; i < n; ++i )
   {
      int j;

      /* on off-diagonal elements, only the sum of corresponding entries need to coincide */
      for( j = 0; j < i; ++j )
         assert(SCIPisRelEQ(scip, consdata->factorleft[i] * consdata->factorright[j] + consdata->factorleft[j] * consdata->factorright[i], 2*a[j*n+i]));

      /* check diagonal elements */
      assert(SCIPisRelEQ(scip, consdata->factorleft[i] * consdata->factorright[i], a[i*n+i]));
   }
#endif

 CLEANUP:
   SCIPfreeBufferArray(scip, &a);
   SCIPfreeBufferArray(scip, &eigvals);

   return SCIP_OKAY;
}

/** gets maximal absolute value in gradient of quadratic function */
static
SCIP_Real getGradientMaxElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      maxelem;
   SCIP_Real      g;
   int            i, j, k;
   SCIP_VAR*      var;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
   {
      maxelem = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
         if( REALABS(consdata->lincoefs[i]) > maxelem )
            maxelem = REALABS(consdata->lincoefs[i]);
   }
   else
   {
      maxelem = consdata->lincoefsmax;
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvarterms[i].var;
      assert(!SCIPisInfinity(scip,  SCIPgetSolVal(scip, sol, var)));
      assert(!SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, var)));
      g  =       consdata->quadvarterms[i].lincoef;
      g += 2.0 * consdata->quadvarterms[i].sqrcoef * SCIPgetSolVal(scip, sol, var);
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         k = consdata->quadvarterms[i].adjbilin[j];
         if( consdata->bilinterms[k].var1 == var )
            g += consdata->bilinterms[k].coef * SCIPgetSolVal(scip, sol, consdata->bilinterms[k].var2);
         else
            g += consdata->bilinterms[k].coef * SCIPgetSolVal(scip, sol, consdata->bilinterms[k].var1);
      }
      if( REALABS(g) > maxelem )
         maxelem = REALABS(g);
   }

   return maxelem;
}

/** computes activity and violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool             doscaling           /**< should we scale the violation by the gradient of the quadratic function ? */ 
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real varval;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->activity = 0.0;
   varval = 0.0;

   /* @todo Take better care of variables at +/- infinity: e.g., run instance waste in debug mode with a short timelimit (30s). */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( SCIPisInfinity(scip, ABS(SCIPgetSolVal(scip, sol, consdata->linvars[i]))) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      consdata->activity += consdata->lincoefs[i] * SCIPgetSolVal(scip, sol, consdata->linvars[i]);
   }

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      varval = SCIPgetSolVal(scip, sol, consdata->quadvarterms[j].var);
      if( SCIPisInfinity(scip, ABS(varval)) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      consdata->activity += (consdata->quadvarterms[j].lincoef + consdata->quadvarterms[j].sqrcoef * varval) * varval;
   }

   for( j = 0; j < consdata->nbilinterms; ++j )
      consdata->activity += consdata->bilinterms[j].coef * SCIPgetSolVal(scip, sol, consdata->bilinterms[j].var1) * SCIPgetSolVal(scip, sol, consdata->bilinterms[j].var2);

   if( consdata->activity < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhsviol = consdata->lhs - consdata->activity;
   else
      consdata->lhsviol = 0.0;

   if( consdata->activity > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhsviol = consdata->activity - consdata->rhs;
   else
      consdata->rhsviol = 0.0;

   if( doscaling && (consdata->lhsviol > 0.0 || consdata->rhsviol > 0.0) )
   {
      SCIP_Real norm;
      norm = getGradientMaxElement(scip, cons, sol);
      if( norm > 1.0 )
      {
         /* scale only if > 1.0, since LP solvers may scale also only if cut norm is > 1 */
         consdata->lhsviol /= norm;
         consdata->rhsviol /= norm;
      }
   }

   return SCIP_OKAY;
}

/** computes violation of a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool             doscaling,          /**< are we scaling when computing violation ? */
   SCIP_CONS**           maxviolcon          /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol;
   int            c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(maxviolcon != NULL);

   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, conss[c], sol, doscaling) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }

   return SCIP_OKAY;
}

/** tries to compute cut for multleft * <coefleft, x'> * multright <= rhs / (multright * <coefright, x'>) where x'=(x,1) */
static
void generateCutFactorableDo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_Real             multleft,           /**< multiplicator on lhs */
   SCIP_Real*            coefleft,           /**< coefficient for factor on lhs */
   SCIP_Real             multright,          /**< multiplicator on both sides */
   SCIP_Real*            coefright,          /**< coefficient for factor that goes to rhs */
   SCIP_Real             rightminactivity,   /**< minimal activity of <coefright, x> */
   SCIP_Real             rightmaxactivity,   /**< maximal activity of <coefright, x> */
   SCIP_Real             rhs,                /**< denominator on rhs */
   SCIP_Real*            cutcoef,            /**< array to store cut coefficients for quadratic variables */
   SCIP_Real*            cutrhs,             /**< buffer to store cut rhs */
   SCIP_Bool*            islocal,            /**< buffer to set to TRUE if local information was used */
   SCIP_Bool*            success,            /**< buffer to indicate whether a cut was successfully computed */
   char*                 name                /**< buffer to store name of cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real constant;
   int i;

   assert(cutcoef != NULL);
   assert(rightminactivity * multright > 0.0);
   assert(rightmaxactivity * multright > 0.0);
   assert(multright == 1.0 || multright == -1.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( rhs > 0.0 )
   {
      /* if rhs > 0.0, then rhs / (multright * <coefright, x'>) is convex, thus need secant:
       *  1 / multright*<coefright, x'> <= 1/minact + 1/maxact - 1/(minact * maxact) multright*<coefright, x'>
       *  where [minact, maxact] = multright * [rightminactivity, rightmaxactivity]
       *
       * assuming multright is either -1 or 1, and substituting gives
       *  multright/rightminactivity + multright/rightmaxactivity  - multright/(rightminactivity * rightmaxactivity) *<coefright, x'>
       *
       * multiplying by rhs, gives the estimate
       *  rhs / (multright * <coefright, x'>) <= rhs * multright * (1/rightminactivity + 1/rightmaxactivity - 1/(rightminactivity * rightmaxactivity) * <coefright, x'>)
       */

      /* cannot do if unbounded */
      if( SCIPisInfinity(scip, rightmaxactivity) )
      {
         *success = FALSE;
         return;
      }

      assert(SCIPisFeasLE(scip, rightminactivity, rightmaxactivity));

      constant = multleft * multright * coefleft[consdata->nquadvars];
      constant -= rhs * multright * (1.0 / rightminactivity + 1.0 / rightmaxactivity);
      constant += rhs * multright * coefright[consdata->nquadvars] / (rightminactivity * rightmaxactivity);

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         cutcoef[i] = multleft * multright * coefleft[i];
         cutcoef[i] += rhs * multright / (rightminactivity * rightmaxactivity) * coefright[i];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_factorablesecant_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
   }
   else
   {
      SCIP_Real refvalue;

      /* if rhs < 0.0, then rhs / (multright * <coefright, x'>) is convex, thus need linearization:
       *     rhs / (multright * <coefright, x'>)
       *  <= rhs / (multright * <coefright, ref'>) - rhs / (multright * <coefright, ref'>)^2 * (multright * <coefright, x'> - multright * <coefright, ref'>)
       *   = 2*rhs / (multright * <coefright, ref'>) - rhs / (multright * <coefright, ref'>)^2 * (multright * <coefright, x'>)
       *
       *  where ref' = (ref, 1)
       */

      /* compute <coefright, ref'> */
      refvalue = coefright[consdata->nquadvars];
      for( i = 0; i < consdata->nquadvars; ++i )
         refvalue += coefright[i] * ref[i];

      /* should not happen, since we checked activity of <coefright,x> before, and assume ref within bounds */
      assert(!SCIPisZero(scip, refvalue));

      constant  = multleft * multright * coefleft[consdata->nquadvars];
      constant -= 2.0 * rhs / (multright * refvalue);
      constant += rhs / (refvalue * refvalue) * multright * coefright[consdata->nquadvars];

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         cutcoef[i] = multleft * multright * coefleft[i];
         cutcoef[i] += rhs / (refvalue * refvalue) * multright * coefright[i];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_factorablelinearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
   }

   *cutrhs = -constant;

   /* @todo does not always need to be local */
   *islocal = TRUE;
   *success = TRUE;
}

/** tries to generate a cut if constraint quadratic function is factorable and there are no linear variables
 * (ax+b)(cx+d) <= rhs and cx+d >= 0 -> (ax+b) <= rhs / (cx+d), where the right hand side is concave and can be linearized
 */
static
SCIP_RETCODE generateCutFactorable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_Real*            cutcoef,            /**< array to store cut coefficients for quadratic variables */
   SCIP_Real*            cutlhs,             /**< buffer to store cut lhs */
   SCIP_Real*            cutrhs,             /**< buffer to store cut rhs */
   SCIP_Bool*            islocal,            /**< buffer to set to TRUE if local information was used */
   SCIP_Bool*            success,            /**< buffer to indicate whether a cut was successfully computed */
   char*                 name                /**< buffer to store name of cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real leftminactivity;
   SCIP_Real leftmaxactivity;
   SCIP_Real rightminactivity;
   SCIP_Real rightmaxactivity;
   SCIP_Real multleft;
   SCIP_Real multright;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(cutcoef != NULL);
   assert(cutlhs  != NULL);
   assert(cutrhs  != NULL);
   assert(islocal != NULL);
   assert(success != NULL);
   assert(name != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nlinvars == 0);
   assert(consdata->factorleft != NULL);
   assert(consdata->factorright != NULL);

   *success = FALSE;
   *cutlhs = -SCIPinfinity(scip);

   leftminactivity = consdata->factorleft[consdata->nquadvars];
   leftmaxactivity = consdata->factorleft[consdata->nquadvars];
   rightminactivity = consdata->factorright[consdata->nquadvars];
   rightmaxactivity = consdata->factorright[consdata->nquadvars];
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( !SCIPisInfinity(scip, -leftminactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, leftmaxactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }

      if( !SCIPisInfinity(scip, -rightminactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, rightmaxactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
   }

   /* write violated constraints as multleft * factorleft * factorright <= rhs */
   if( violside == SCIP_SIDETYPE_RIGHT )
   {
      rhs = consdata->rhs;
      multleft = 1.0;
   }
   else
   {
      rhs = -consdata->lhs;
      multleft = -1.0;
   }

   if( SCIPisZero(scip, rhs) )
   {
      /* @todo do something for rhs == 0.0? */
      return SCIP_OKAY;
   }

   if( !SCIPisFeasPositive(scip, leftminactivity) && !SCIPisFeasNegative(scip, leftmaxactivity) )
   {
      /* left factor has 0 within activity bounds, or is very close, at least */
      if( !SCIPisFeasPositive(scip, rightminactivity) && !SCIPisFeasNegative(scip, rightmaxactivity) )
      {
         /* right factor also has 0 within activity bounds, or is very close, at least
          * -> cannot separate
          */
         return SCIP_OKAY;
      }

      /* write violated constraint as multleft * factorleft * multright * (multright * factorright) <= rhs
       *   such that multright * factorright > 0.0
       */
      if( rightminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorleft * multright <= rhs / (factorright * multright) */
      generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorleft, multright, consdata->factorright, rightminactivity, rightmaxactivity, rhs, cutcoef, cutrhs, islocal, success, name);
   }
   else if( !SCIPisFeasPositive(scip, rightminactivity) && !SCIPisFeasNegative(scip, rightmaxactivity) )
   {
      /* left factor is bounded away from 0
       * right factor has 0 within activity bounds, or is very close, at least
       * -> so divide by left factor
       */

      /* write violated constraint as multleft * factorright * multright * (multright * factorleft) <= rhs
       *   such that multright * factorleft > 0.0
       */
      if( leftminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorright * multright <= rhs / (factorleft * multright) */
      generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorright, multright, consdata->factorleft, leftminactivity, leftmaxactivity, rhs, cutcoef, cutrhs, islocal, success, name);
   }
   else if( SCIPisInfinity(scip, -leftminactivity) || SCIPisInfinity(scip, leftmaxactivity) ||
      (!SCIPisInfinity(scip, -rightminactivity) && !SCIPisInfinity(scip, rightmaxactivity) && rightmaxactivity - rightminactivity < leftmaxactivity - leftminactivity) )
   {
      /* both factors are bounded away from 0, but the right one has a smaller activity range, so divide by that one */

      /* write violated constraint as multleft * factorleft * multright * (multright * factorright) <= rhs
       *   such that multright * factorright > 0.0
       */
      if( rightminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorleft * multright <= rhs / (factorright * multright) */
      generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorleft, multright, consdata->factorright, rightminactivity, rightmaxactivity, rhs, cutcoef, cutrhs, islocal, success, name);
   }
   else
   {
      /* both factors are bounded away from 0, but the left one has a smaller activity range, so divide by that one */

      /* write violated constraint as multleft * factorright * multright * (multright * factorleft) <= rhs
       *   such that multright * factorleft > 0.0
       */
      if( leftminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorright * multright <= rhs / (factorleft * multright) */
      generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorright, multright, consdata->factorleft, leftminactivity, leftmaxactivity, rhs, cutcoef, cutrhs, islocal, success, name);
   }

   return SCIP_OKAY;
}

/* finds intersections of a parametric line (x,y) = (x0,y0) + t [(x1,y1) - (x0,y0)] on curves x*y = wl and x*y = wu
 * returns TRUE if unsuccessful and FALSE otherwise
 */
static
SCIP_Bool generateCutLTIfindIntersection(
   SCIP*                 scip,
   SCIP_Real             x0,
   SCIP_Real             y0_,
   SCIP_Real             x1,
   SCIP_Real             y1_,
   SCIP_Real             wl,
   SCIP_Real             wu,
   SCIP_Real*            xl,
   SCIP_Real*            yl,
   SCIP_Real*            xu,
   SCIP_Real*            yu
   )
{
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real c;
   SCIP_Real tl;
   SCIP_Real tu;

   assert(wl == SCIP_INVALID || (xl != NULL && yl != NULL));  /*lint !e777 */
   assert(wu == SCIP_INVALID || (xu != NULL && yu != NULL));  /*lint !e777 */

   /* The parametric line is of the form
    *
    *  x = x0 + t (x1-x0)
    *  y = y0 + t (y1-y0)
    *
    * and for that to satisfy xy = wl and xy = wu we must have
    *
    * x0 y0 + t [x0 (y1-y0) + y0 (x1-x0)] + t^2 (x1-x0) (y1-y0) = wl
    *                                                           = wu
    *
    * or a t^2 + b t + c - wl = 0 for proper values of a,b,c.
    *    a t^2 + b t + c - wu = 0
    *
    * Because of the way this procedure will be used, one of the two
    * solutions found we must always use the minimum nonnegative one
    */

   a = (x1 - x0) * (y1_ - y0_);
   c = x0 * y0_;
   b = x0 * y1_ + y0_ * x1 - 2.0 * c;

   tl = 0.0;
   tu = 0.0;

   if( !SCIPisZero(scip, (SCIP_Real)a) )
   {
      if( wl != SCIP_INVALID )  /*lint !e777 */
      {
         SCIP_Real tl1;
         SCIP_Real tl2;
         SCIP_Real denom;

         assert(b * b - 4.0 * a * (c - wl) >= 0.0);
         denom = sqrt(b * b - 4.0 * a * (c - wl));
         tl1 = (-b - denom) / (2.0 * a);
         tl2 = (-b + denom) / (2.0 * a);
         tl = (tl1 < 0.0) ? tl2 : tl1;
      }

      if( wu != SCIP_INVALID )  /*lint !e777 */
      {
         SCIP_Real tu1;
         SCIP_Real tu2;
         SCIP_Real denom;

         assert(b * b - 4.0 * a * (c - wu) >= 0.0);
         denom = sqrt(b * b - 4.0 * a * (c - wu));
         tu1 = (-b - denom) / (2.0 * a);
         tu2 = (-b + denom) / (2.0 * a);
         tu = (tu1 < 0.0) ? tu2 : tu1;
      }
   }
   else if( !SCIPisZero(scip, (SCIP_Real)b) )
   {
      if( wl != SCIP_INVALID )  /*lint !e777 */
         tl = (wl - c) / b;
      if( wu != SCIP_INVALID )  /*lint !e777 */
         tu = (wu - c) / b;
   }
   else
   {
      /* no or infinitely many solutions */
      return TRUE;
   }

   if( wl != SCIP_INVALID )  /*lint !e777 */
   {
      *xl = (SCIP_Real)(x0  + tl * (x1  - x0 ));
      *yl = (SCIP_Real)(y0_ + tl * (y1_ - y0_));

      if( !SCIPisRelEQ(scip, *xl * *yl, wl) )
      {
         SCIPdebugMessage("probable numerical difficulties, give up\n");
         return TRUE;
      }
   }

   if( wu != SCIP_INVALID )  /*lint !e777 */
   {
      *xu = (SCIP_Real)(x0  + tu * (x1 -  x0));
      *yu = (SCIP_Real)(y0_ + tu * (y1_ - y0_));

      if( !SCIPisRelEQ(scip, *xu * *yu, wu) )
      {
         SCIPdebugMessage("probable numerical difficulties, give up\n");
         return TRUE;
      }
   }

   return FALSE;
}

/* generate coefficients for a plane through points (x1, y1_, x1*y1) and (x2, y2, x2*y2)
 * such that intersecting it with one of them (the first if whichuse is FALSE, the second otherwise)
 * gives a tangent to the curve x*y = k
 *
 * returns TRUE on error and FALSE on success
 */
static
SCIP_Bool generateCutLTIgenMulCoeff(
   SCIP*                 scip,
   SCIP_Real             x1,
   SCIP_Real             y1_,
   SCIP_Real             x2,
   SCIP_Real             y2,
   SCIP_Bool             whichuse,
   SCIP_Real*            cx,
   SCIP_Real*            cy,
   SCIP_Real*            cw
   )
{
   SCIP_Real xd;
   SCIP_Real yd;
   SCIP_Real xo;
   SCIP_Real yo;

   assert(cx != NULL);
   assert(cy != NULL);
   assert(cw != NULL);

   /* the x-y slope of this constraint must be tangent to a curve x*y = k at (xD,yD) */
   if( !whichuse )
   {
      xd = x1;
      xo = x2;
      yd = y1_;
      yo = y2;
   }
   else
   {
      xd = x2;
      xo = x1;
      yd = y2;
      yo = y1_;
   }

   *cx = yd;
   *cy = xd;

   /* lift it so that it touches the other curve */

   /* if the two points are on the same curve, then no cut */
   if( SCIPisZero(scip, xo * yo - xd * yd) )
      return TRUE;

   /* should ALWAYS be negative */
   *cw = (2.0 * xd * yd - (*cx * xo + *cy * yo)) / (xo * yo - xd * yd);

   return FALSE;
}

/** computes coefficients of a lifted-tangent inequality for x*y = w
 * The code is an adaptation of the methods in exprMul-upperHull.cpp in Couenne/stable/0.4 rev773,
 * written by P. Belotti and licensed under Eclipse Public License.
 */
static
void generateCutLTIcomputeCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xl,                 /**< lower bound on x */
   SCIP_Real             xu,                 /**< upper bound on x */
   SCIP_Real             x0,                 /**< reference point for x */
   SCIP_Real             yl,                 /**< lower bound on y */
   SCIP_Real             yu,                 /**< upper bound on y */
   SCIP_Real             y0_,                /**< reference point for y */
   SCIP_Real             wl,                 /**< lower bound on w */
   SCIP_Real             wu,                 /**< upper bound on w */
   SCIP_Real             w0,                 /**< reference point for w */
   SCIP_Real*            cx,                 /**< buffer where to store cut coefficient for x */
   SCIP_Real*            cy,                 /**< buffer where to store cut coefficient for y */
   SCIP_Real*            cw,                 /**< buffer where to store cut coefficient for w */
   SCIP_Real*            c0,                 /**< buffer where to store cut left-hand-side */
   SCIP_Bool*            success             /**< buffer where to indicate whether cut coefficients were computed */
   )
{
   SCIP_Bool flipx;
   SCIP_Bool flipy;
   SCIP_Bool flipw;
   SCIP_Real tmp;
   SCIP_Real xlow;
   SCIP_Real ylow;
   SCIP_Real xupp;
   SCIP_Real yupp;
   SCIP_Real c0x;
   SCIP_Real c0y;
   SCIP_Real c0w;

   assert(scip != NULL);
   assert(cx != NULL);
   assert(cy != NULL);
   assert(cw != NULL);
   assert(c0 != NULL);
   assert(success != NULL);

   *success = FALSE;
   *cx = 0.0;
   *cy = 0.0;
   *cw = 0.0;
   *c0 = 0.0;

   SCIPdebugMessage("entering points:\n");
   SCIPdebugMessage("x: %9g\t[%9g\t%9g]\n", x0, xl, xu);
   SCIPdebugMessage("y: %9g\t[%9g\t%9g]\n", y0_, yl, yu);
   SCIPdebugMessage("w: %9g\t[%9g\t%9g]\n", w0, wl, wu);

   /* generateCutLTI should have recognized these */
   assert(wl >= 0.0 || wu <= 0.0);
   assert(!SCIPisInfinity(scip, -wl));
   assert(!SCIPisInfinity(scip,  wu));

   assert(SCIPisFeasGE(scip, x0, xl));
   assert(SCIPisFeasLE(scip, x0, xu));
   assert(SCIPisFeasGE(scip, y0_, yl));
   assert(SCIPisFeasLE(scip, y0_, yu));

   /* preliminary bound tightening */
   if( wl >= 0.0 )
   {
      if( xl >= 0.0 || yl >= 0.0 || SCIPisLT(scip, xl * yl, wl) )
      {
         xl = MAX(xl, 0.0);
         yl = MAX(yl, 0.0);
      }
      else if( xu <= 0.0 || yu <= 0.0 || SCIPisLT(scip, xu * yu, wl) )
      {
         xu = MIN(xu, 0.0);
         yu = MIN(yu, 0.0);
      }
      else
      {
         /* both variables have mixed sign (xl < 0 && xu > 0 && yl < 0 && yu > 0) and both xl*yl and xu*yu are feasible
          * cannot generate cut for this
          */
         return;
      }
   }
   else
   {
      if( xl >= 0.0 || yu <= 0.0 || SCIPisGT(scip, xl * yu, wu) )
      {
         xl = MAX(xl, 0.0);
         yu = MIN(yu, 0.0);
      }
      else if( xu <= 0.0 || yl >= 0.0 || SCIPisGT(scip, xu * yl, wu))
      {
         xu = MIN(xu, 0.0);
         yl = MAX(yl, 0.0);
      }
      else
      {
         /* both variables have mixed sign (xl < 0 && xu > 0 && yl < 0 && yu > 0) and both xl*yu and xu*yl are feasible
          * cannot generate cut for this
          */
         return;
      }
   }

   /* if x or y is fixed now or even infeasible, then do not think about a cut */
   if( SCIPisGE(scip, xl, xu) || SCIPisGE(scip, yl, yu) )
      return;

   /* reduce to positive orthant by flipping variables */
   if( xl < 0.0 )
   {
      flipx = TRUE;
      tmp = xu;
      xu = -xl;
      xl = -tmp;
      x0 = -x0;
   }
   else
      flipx = FALSE;

   if( yl < 0.0 )
   {
      flipy = TRUE;
      tmp = yu;
      yu = -yl;
      yl = -tmp;
      y0_ = -y0_;
   }
   else
      flipy = FALSE;

   if( flipx ^ flipy )
   {
      flipw = TRUE;
      tmp = wu;
      wu = -wl;
      wl = -tmp;
      w0 = -w0;
   }
   else
      flipw = FALSE;

   /* project refpoint into box not only for numerical reasons, but also due to preliminary bound tightening above */
   x0 = MIN(xu, MAX(x0, xl));
   y0_ = MIN(yu, MAX(y0_, yl));
   w0 = MIN(wu, MAX(w0, wl));

   SCIPdebugMessage("reduced points:\n");
   SCIPdebugMessage("x: %9g\t[%9g\t%9g]\n", x0, xl, xu);
   SCIPdebugMessage("y: %9g\t[%9g\t%9g]\n", y0_, yl, yu);
   SCIPdebugMessage("w: %9g\t[%9g\t%9g]\n", w0, wl, wu);

   if( SCIPisGE(scip, xl * yl, wl) && SCIPisLE(scip, xu * yu, wu) )
   {
      SCIPdebugMessage("box for x and y inside feasible region -> nothing to separate\n");
      return;
   }
   if( SCIPisGE(scip, x0 * y0_, w0) )
   {
      SCIPdebugMessage("point to separate not below curve -> cannot separate\n");
      return;
   }

   /* find intersections of halfline from origin
    * return if no proper point could be found
    */
   if( generateCutLTIfindIntersection(scip, 0.0, 0.0, x0, y0_, wl, wu, &xlow, &ylow, &xupp, &yupp) )
      return;

   SCIPdebugMessage("intersections:\n");
   SCIPdebugMessage("lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
   SCIPdebugMessage("upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

   /* Case 1: If both are outside of bounding box, either NW or SE, then McCormick is sufficient, so return */
   if( (xlow <= xl && yupp >= yu) || (ylow <= yl && xupp >= xu) )
      return;

   /* There will be at least one cut. Define coefficients and rhs ---will have to change them back if (flipX || flipY) */
   if( xlow >= xl && xupp <= xu && ylow >= yl && yupp <= yu )
   {
      /* Case 2: both are inside. Easy lifting... */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
         return;

      c0x = *cx * xlow;
      c0y = *cy * ylow;
      c0w = *cw * wl;
   }
   else if( xlow >= xl && ylow >= yl && (xupp > xu || yupp > yu) )
   {
      /* Case 3a and 3b: through lower curve, but not upper. */
      if( yupp > yu )
      {
         /* upper intersect is North; place it within box */
         assert(!SCIPisInfinity(scip, yu));
         yupp = yu;
         xupp = wu / yu;
      }
      else
      {
         /* upper intersect is East; place it within box */
         assert(!SCIPisInfinity(scip, xu));
         xupp = xu;
         yupp = wu / xu;
      }

      /* find intersection on low curve on half line through new point and (x0,y0_) */
      if( generateCutLTIfindIntersection(scip, xupp, yupp, x0, y0_, wl, SCIP_INVALID, &xlow, &ylow, NULL, NULL) )
         return;

      /* check whether McCormick is sufficient */
      if( xlow < xl || ylow < yl )
         return;

      /* lift inequality on lower point */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
         return;

      c0x = *cx * xlow;
      c0y = *cy * ylow;
      c0w = *cw * wl;
   }
   else if( xupp <= xu && yupp <= yu && (xlow < xl || ylow < yl) )
   {
      /* Case 4a and 4b: viceversa (lift for validity) */
      if( ylow < yl )
      {
         /* upper intersect is South; place it within box */
         assert(!SCIPisZero(scip, yl));
         ylow = yl;
         xlow = wl / yl;
      }
      else
      {
         /* upper intersect is West; place it within box */
         assert(!SCIPisZero(scip, xl));
         xlow = xl;
         ylow = wl / xl;
      }

      /* find intersection on low curve on half line through new point and (x0,y0) */
      if( generateCutLTIfindIntersection(scip, xlow, ylow, x0, y0_, SCIP_INVALID, wu, NULL, NULL, &xupp, &yupp) )
         return;

      /* check whether McCormick is sufficient */
      if( xupp > xu || yupp > yu )
         return;

      /* lift inequality on UPPER point */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, TRUE, cx, cy, cw) )
         return;

      c0x = *cx * xupp;
      c0y = *cy * yupp;
      c0w = *cw * wu;
   }
   else if( (xlow < xl && xupp > xu) || (ylow < yl && yupp > yu) )
   {
      /* Case 5: both outside of bounding box, N and S or W and E. */
#if 0
      SCIP_Real xlow2;
      SCIP_Real ylow2;
      SCIP_Real xupp2;
      SCIP_Real yupp2;
#endif

      if( ylow < yl )
      {
         /* upper intersect is South; place it within box */
         assert(!SCIPisZero(scip, yl));
         assert(!SCIPisZero(scip, yu));
         ylow = yl;
         yupp = yu;
         xlow = wl / yl;
         xupp = wu / yu;
      }
      else
      {
         /* upper intersect is West; place it within box */
         assert(!SCIPisZero(scip, xl));
         assert(!SCIPisZero(scip, xu));
         xlow = xl;
         xupp = xu;
         ylow = wl / xl;
         yupp = wu / xu;
      }

      SCIPdebugMessage("New intersections:\n");
      SCIPdebugMessage("lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
      SCIPdebugMessage("upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

#if 1
      /* Nothing to find. Just separate two inequalities at the same point, just using different support */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
      {
         if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, TRUE, cx, cy, cw) )
            return;

         c0x = *cx * xupp;
         c0y = *cy * yupp;
         c0w = *cw * wu;
      }
      else
      {
         c0x = *cx * xlow;
         c0y = *cy * ylow;
         c0w = *cw * wl;
      }

#else
      /* find the intersection on the lower (upper) curve on the line through xLP and the upper (lower) point
       * this does not seem to work (cuts off solution at nous2), so it is disabled for now
       */
      if( generateCutLTIfindIntersection(scip, xlow, ylow, x0, y0_, SCIP_INVALID, wu, NULL, NULL, &xupp2, &yupp2) ||
         generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp2, yupp2, FALSE, cx, cx, cw) )
      {
         if( generateCutLTIfindIntersection(scip, xupp, yupp, x0, y0_, wl, SCIP_INVALID, &xlow2, &ylow2, NULL, NULL) ||
            generateCutLTIgenMulCoeff(scip, xlow2, ylow2, xupp, yupp, TRUE, cx, cy, cw) )
            return;

         c0x = *cx * xupp;
         c0y = *cy * yupp;
         c0w = *cw * wu;
      }
      else
      {
         c0x = *cx * xlow;
         c0y = *cy * ylow;
         c0w = *cw * wl;
      }
#endif
   }
   else
   {
      SCIPdebugMessage("points are in a weird position:\n");
      SCIPdebugMessage("lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
      SCIPdebugMessage("upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

      return;
   }

   SCIPdebugMessage("cut w.r.t. reduced points: %gx-%g %+gy-%g %+gw-%g >= 0\n",
      *cx, c0x, *cy, c0y, *cw, c0w);

   /* re-transform back into original variables */
   if( flipx )
      *cx = -*cx;
   if( flipy )
      *cy = -*cy;
   if( flipw )
      *cw = -*cw;

   *c0 = c0x + c0y + c0w;

   *success = TRUE;
}

/** tries to generate a cut if constraint quadratic function is factorable and there are linear variables
 * computes what is called a lifted tangent inequality described in
 *   Belotti, Miller, Namazifar, Lifted inequalities for bounded products of variables, SIAG/OPT Views-and-News 22:1, 2011
 */
static
SCIP_RETCODE generateCutLTI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_SOL*             sol,                /**< solution that shall be cutoff, NULL for LP solution */
   SCIP_Real**           cutcoeflin,         /**< buffer to store pointer to array with coefficients for linear variables */
   SCIP_Real*            cutcoefquad,        /**< array to store cut coefficients for quadratic variables */
   SCIP_Real*            cutlhs,             /**< buffer to store cut lhs */
   SCIP_Real*            cutrhs,             /**< buffer to store cut rhs */
   SCIP_Bool*            islocal,            /**< buffer to set to TRUE if local information was used */
   SCIP_Bool*            success,            /**< buffer to indicate whether a cut was successfully computed */
   char*                 name                /**< buffer to store name of cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real leftminactivity;
   SCIP_Real leftmaxactivity;
   SCIP_Real leftrefactivity;
   SCIP_Real rightminactivity;
   SCIP_Real rightmaxactivity;
   SCIP_Real rightrefactivity;
   SCIP_Real rhsminactivity;
   SCIP_Real rhsmaxactivity;
   SCIP_Real rhsrefactivity;
   SCIP_Real coefleft;
   SCIP_Real coefright;
   SCIP_Real coefrhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(cutcoeflin != NULL);
   assert(cutcoefquad != NULL);
   assert(cutlhs  != NULL);
   assert(cutrhs  != NULL);
   assert(islocal != NULL);
   assert(success != NULL);
   assert(name != NULL);
   /* currently only separate LP solution or solutions given as SCIP_SOL, i.e., no cutgeneration during initlp */
   assert(sol != NULL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nlinvars > 0);
   assert(consdata->factorleft != NULL);
   assert(consdata->factorright != NULL);

   *success = FALSE;
   *cutlhs = -SCIPinfinity(scip); /* for compiler */

   /* write violated constraints as factorleft * factorright '==' rhs
    * where rhs are constraint sides - activity bound of linear part
    */
   rhsminactivity = consdata->lhs;
   rhsmaxactivity = consdata->rhs;
   rhsrefactivity = (violside == SCIP_SIDETYPE_LEFT ? consdata->lhs : consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( !SCIPisInfinity(scip, -rhsminactivity) )
      {
         if( consdata->lincoefs[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i])) )
               rhsminactivity = -SCIPinfinity(scip);
            else
               rhsminactivity -= consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }
         else
         {
            assert(consdata->lincoefs[i] > 0.0);
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i])) )
               rhsminactivity = -SCIPinfinity(scip);
            else
               rhsminactivity -= consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }
      }
      if( !SCIPisInfinity(scip, rhsmaxactivity) )
      {
         if( consdata->lincoefs[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i])) )
               rhsmaxactivity = SCIPinfinity(scip);
            else
               rhsmaxactivity -= consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }
         else
         {
            assert(consdata->lincoefs[i] > 0.0);
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i])) )
               rhsmaxactivity = SCIPinfinity(scip);
            else
               rhsmaxactivity -= consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }
      }
      rhsrefactivity -= consdata->lincoefs[i] * SCIPgetSolVal(scip, sol, consdata->linvars[i]);
   }

   if( SCIPisInfinity(scip, -rhsminactivity) || SCIPisInfinity(scip, rhsmaxactivity) )
   {
      /* if right hand side is unbounded, then cannot do LTI */
      return SCIP_OKAY;
   }

   if( !SCIPisFeasPositive(scip, rhsminactivity) && !SCIPisFeasNegative(scip, rhsmaxactivity) )
   {
      /* if right hand side has 0 inside activity, then cannot do anything
       * if it has 0.0 as min or max activity, then a usual McCormick should be sufficient, too
       */
      return SCIP_OKAY;
   }

   leftminactivity = consdata->factorleft[consdata->nquadvars];
   leftmaxactivity = consdata->factorleft[consdata->nquadvars];
   leftrefactivity = consdata->factorleft[consdata->nquadvars];
   rightminactivity = consdata->factorright[consdata->nquadvars];
   rightmaxactivity = consdata->factorright[consdata->nquadvars];
   rightrefactivity = consdata->factorright[consdata->nquadvars];
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( !SCIPisInfinity(scip, -leftminactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, leftmaxactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
      leftrefactivity += consdata->factorleft[i] * ref[i];

      if( !SCIPisInfinity(scip, -rightminactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, rightmaxactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
      rightrefactivity += consdata->factorright[i] * ref[i];
   }

   /* success can only be expected for separation of violated x*y <= w, assuming x>=0, y>=0
    * @todo we should check this early? */

   /* call Couenne magic */
   generateCutLTIcomputeCoefs(scip,
      leftminactivity, leftmaxactivity, leftrefactivity,
      rightminactivity, rightmaxactivity, rightrefactivity,
      rhsminactivity, rhsmaxactivity, rhsrefactivity,
      &coefleft, &coefright, &coefrhs, cutlhs,
      success);

   if( !*success )
      return SCIP_OKAY;

   SCIPdebugMessage("LTI for x[%g,%g] * y[%g,%g] = w[%g,%g]: %gx %+gy %+gw >= %g;  feas: %g\n",
      leftminactivity, leftmaxactivity, rightminactivity, rightmaxactivity, rhsminactivity, rhsmaxactivity,
      coefleft, coefright, coefrhs, *cutlhs,
      coefleft * leftrefactivity + coefright * rightrefactivity + coefrhs * rhsrefactivity - *cutlhs
      );

   if( coefleft * leftrefactivity + coefright * rightrefactivity + coefrhs * rhsrefactivity >= *cutlhs )
   {
      SCIPdebugMessage("does not cutoff point? :-(\n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* setup cut coefs for
    *   coefleft * leftfactor + coefright * rightfactor + coefrhs * w >= cutlhs, where conslhs - lincoefs <= w <= consrhs - lincoefs
    */
   for( i = 0; i < consdata->nquadvars; ++i )
      cutcoefquad[i] = coefleft * consdata->factorleft[i] + coefright * consdata->factorright[i];
   assert(i == consdata->nquadvars);
   *cutlhs -= coefleft * consdata->factorleft[i] + coefright * consdata->factorright[i];

   SCIP_CALL( SCIPallocBufferArray(scip, cutcoeflin, consdata->nlinvars) );
   for( i = 0; i < consdata->nlinvars; ++i )
      (*cutcoeflin)[i] = -coefrhs * consdata->lincoefs[i];
   if( coefrhs > 0.0 )
   {
      /* use coefrhs * w <= coefrhs * (consrhs - lincoefs) */
      assert(!SCIPisInfinity(scip, consdata->rhs));
      *cutlhs -= coefrhs * consdata->rhs;
   }
   else
   {
      /* use coefrhs * w <= coeflhs * (conslhs - lincoefs) */
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      *cutlhs -= coefrhs * consdata->lhs;
   }

   *cutrhs = SCIPinfinity(scip);
   *islocal = TRUE;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_lti_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   *success = TRUE;

   return SCIP_OKAY;
}

/** computes coefficients of linearization of a square term in a reference point */
static
void addSquareLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             refpoint,           /**< point where to linearize */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpoint)) )
   {
      *success = FALSE;
      return;
   }

   if( !isint || SCIPisIntegral(scip, refpoint) )
   {
      SCIP_Real tmp;

      /* sqrcoef * x^2  ->  tangent in refpoint = sqrcoef * 2 * refpoint * (x - refpoint) */

      tmp = sqrcoef * refpoint;

      if( SCIPisInfinity(scip, 2.0 * REALABS(tmp)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef += 2.0 * tmp;
      tmp *= refpoint;
      *linconstant -= tmp;
   }
   else
   {
      /* sqrcoef * x^2 ->  secant between f=floor(refpoint) and f+1 = sqrcoef * (f^2 + ((f+1)^2 - f^2) * (x-f)) = sqrcoef * (-f*(f+1) + (2*f+1)*x) */
      SCIP_Real f;
      SCIP_Real coef;
      SCIP_Real constant;

      f = SCIPfloor(scip, refpoint);

      coef     =  sqrcoef * (2.0 * f + 1.0);
      constant = -sqrcoef * f * (f + 1.0);

      if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef     += coef;
      *linconstant += constant;
   }
}

/** computes coefficients of secant of a square term */
static
void addSquareSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real             refpoint,           /**< point for which to compute value of linearization */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   )
{
   SCIP_Real coef;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lb));
   assert(!SCIPisInfinity(scip, -ub));
   assert(SCIPisLE(scip, lb, ub));
   assert(SCIPisLE(scip, lb, refpoint));
   assert(SCIPisGE(scip, ub, refpoint));
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
   {
      /* unboundedness */
      *success = FALSE;
      return;
   }

   /* sqrcoef * x^2 -> sqrcoef * (lb * lb + (ub*ub - lb*lb)/(ub-lb) * (x-lb)) = sqrcoef * (lb*lb + (ub+lb)*(x-lb)) = sqrcoef * ((lb+ub)*x - lb*ub) */

   coef     =  sqrcoef * (lb + ub);
   constant = -sqrcoef * lb * ub;
   if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoef     += coef;
   *linconstant += constant;
}

/** computes coefficients of linearization of a bilinear term in a reference point */
static
void addBilinLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             refpointx,          /**< point where to linearize first  variable */
   SCIP_Real             refpointy,          /**< point where to linearize second variable */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;

   assert(scip != NULL);
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpointx)) || SCIPisInfinity(scip, REALABS(refpointy)) )
   {
      *success = FALSE;
      return;
   }

   /* bilincoef * x * y ->  bilincoef * (refpointx * refpointy + refpointy * (x - refpointx) + refpointx * (y - refpointy))
    *                    = -bilincoef * refpointx * refpointy + bilincoef * refpointy * x + bilincoef * refpointx * y
    */

   constant = -bilincoef * refpointx * refpointy;

   if( SCIPisInfinity(scip, REALABS(bilincoef * refpointx)) || SCIPisInfinity(scip, REALABS(bilincoef * refpointy)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoefx    += bilincoef * refpointy;
   *lincoefy    += bilincoef * refpointx;
   *linconstant += constant;
}

/** computes coefficients of McCormick under- or overestimation of a bilinear term */
static
void addBilinMcCormick(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;
   SCIP_Real coefx;
   SCIP_Real coefy;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( SCIPisRelEQ(scip, lbx, ubx) )
   {
      /* x is fixed, so bilinear term is at most linear */
      if( SCIPisRelEQ(scip, lby, uby) )
      {
         /* also y is fixed, so bilinear term is constant */
         coefx    = 0.0;
         coefy    = 0.0;
         constant = bilincoef * refpointx * refpointy;
      }
      else
      {
         coefx    = 0.0;
         coefy    = bilincoef * refpointx;
         constant = 0.0;
      }
   }
   else if( SCIPisRelEQ(scip, lby, uby) )
   {
      /* y is fixed, so bilinear term is linear */
      coefx    = bilincoef * refpointy;
      coefy    = 0.0;
      constant = 0.0;
   }
   else
   {
      /* both x and y are not fixed */
      if( overestimate )
         bilincoef = -bilincoef;

      if( bilincoef > 0.0 )
      {
         if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, -lby) &&
            (SCIPisInfinity(scip,  ubx) || SCIPisInfinity(scip,  uby) || (uby - refpointy) * (ubx - refpointx) >= (refpointy - lby) * (refpointx - lbx)) )
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * lby;
         }
         else if( !SCIPisInfinity(scip, ubx) && !SCIPisInfinity(scip, uby) )
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * uby;
         }
         else
         {
            *success = FALSE;
            return;
         }
      }
      else
      {
         /* bilincoef < 0.0 */
         if( !SCIPisInfinity(scip,  ubx) && !SCIPisInfinity(scip, -lby) &&
            (SCIPisInfinity(scip, -lbx) || SCIPisInfinity(scip,  uby) || (ubx - lbx) * (refpointy - lby) <= (uby - lby) * (refpointx - lbx)) )
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * lby;
         }
         else if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, uby) )
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * uby;
         }
         else
         {
            *success = FALSE;
            return;
         }
      }

      if( overestimate )
      {
         coefx    = -coefx;
         coefy    = -coefy;
         constant = -constant;
      }
   }

   if( SCIPisInfinity(scip, REALABS(coefx)) || SCIPisInfinity(scip, REALABS(coefy)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoefx    += coefx;
   *lincoefy    += coefy;
   *linconstant += constant;
}

/** computes cut coefficients by linearizing a quadratic function */
static
SCIP_RETCODE generateCutConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< side for which to generate cut */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_Real*            coef,               /**< array to store cut coefficients w.r.t. quadratic variables */
   SCIP_Real*            lhs,                /**< buffer to store left-hand-side of cut */
   SCIP_Real*            rhs,                /**< buffer to store right-hand-side of cut */
   SCIP_Bool*            islocal,            /**< buffer to set to TRUE if local bounds were used */
   SCIP_Bool*            success,            /**< buffer to indicate whether a cut was successfully computed */
   char*                 name                /**< buffer to store name for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real constant;
   SCIP_VAR* var;
   int var2pos;
   int j;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(coef != NULL);
   assert(lhs  != NULL);
   assert(rhs  != NULL);
   assert(islocal != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   constant = 0.0;
   BMSclearMemoryArray(coef, consdata->nquadvars);
   *success = TRUE;

   /* do first-order Taylor for each term */
   for( j = 0; j < consdata->nquadvars && *success; ++j )
   {
      /* initialize coefficients to linear coefficients of quadratic variables */
      coef[j] += consdata->quadvarterms[j].lincoef;

      /* add linearization of square term */
      var = consdata->quadvarterms[j].var;
      addSquareLinearization(scip, consdata->quadvarterms[j].sqrcoef, ref[j], consdata->quadvarterms[j].nadjbilin == 0 && SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef[j], &constant, success);

      /* add linearization of bilinear terms that have var as first variable */
      for( k = 0; k < consdata->quadvarterms[j].nadjbilin && *success; ++k )
      {
         bilinterm = &consdata->bilinterms[consdata->quadvarterms[j].adjbilin[k]];
         if( bilinterm->var1 != var )
            continue;
         assert(bilinterm->var2 != var);
         assert(consdata->sepabilinvar2pos != NULL);

         var2pos = consdata->sepabilinvar2pos[consdata->quadvarterms[j].adjbilin[k]];
         assert(var2pos >= 0);
         assert(var2pos < consdata->nquadvars);
         assert(consdata->quadvarterms[var2pos].var == bilinterm->var2);

         addBilinLinearization(scip, bilinterm->coef, ref[j], ref[var2pos], &coef[j], &coef[var2pos], &constant, success);
      }
   }

   if( !*success )
   {
      SCIPdebugMessage("no success in linearization of <%s> in reference point\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   if( violside == SCIP_SIDETYPE_LEFT )
   {
      *lhs = consdata->lhs - constant;
      *rhs = SCIPinfinity(scip);
   }
   else
   {
      *lhs = -SCIPinfinity(scip);
      *rhs = consdata->rhs - constant;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_side%d_linearization_%d", SCIPconsGetName(cons), violside, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** computes cut coefficients for a nonconvex quadratic function */
static
SCIP_RETCODE generateCutNonConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< side for which to generate cut */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_Real*            coef,               /**< array to store cut coefficients w.r.t. quadratic variables */
   SCIP_Real*            lhs,                /**< buffer to store left-hand-side of cut */
   SCIP_Real*            rhs,                /**< buffer to store right-hand-side of cut */
   SCIP_Bool*            islocal,            /**< buffer to set to TRUE if local bounds were used */
   SCIP_Bool*            success,            /**< buffer to indicate whether a cut was successfully computed */
   char*                 name                /**< buffer to store name for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real constant;
   SCIP_Real sqrcoef;
   SCIP_VAR* var;
   int var2pos;
   int j;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(coef != NULL);
   assert(lhs  != NULL);
   assert(rhs  != NULL);
   assert(islocal != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *islocal = TRUE;
   constant = 0.0;
   BMSclearMemoryArray(coef, consdata->nquadvars);
   *success = TRUE;

   /* underestimate (secant, McCormick) or linearize each term separately */
   for( j = 0; j < consdata->nquadvars && *success; ++j )
   {
      /* initialize coefficients to linear coefficients of quadratic variables */
      coef[j] += consdata->quadvarterms[j].lincoef;

      var = consdata->quadvarterms[j].var;

      sqrcoef = consdata->quadvarterms[j].sqrcoef;
      if( sqrcoef != 0.0 )
      {
         if( (violside == SCIP_SIDETYPE_LEFT  && sqrcoef <= 0.0) || (violside == SCIP_SIDETYPE_RIGHT && sqrcoef > 0.0) )
         {
            /* convex -> linearize */
            addSquareLinearization(scip, sqrcoef, ref[j], SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef[j], &constant, success);
         }
         else
         {
            /* not convex -> secant approximation */
            addSquareSecant(scip, sqrcoef, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), ref[j], &coef[j], &constant, success);
         }
      }

      for( k = 0; k < consdata->quadvarterms[j].nadjbilin && *success; ++k )
      {
         bilinterm = &consdata->bilinterms[consdata->quadvarterms[j].adjbilin[k]];
         if( bilinterm->var1 != var )
            continue;
         assert(bilinterm->var2 != var);
         assert(consdata->sepabilinvar2pos != NULL);

         var2pos = consdata->sepabilinvar2pos[consdata->quadvarterms[j].adjbilin[k]];
         assert(var2pos >= 0);
         assert(var2pos < consdata->nquadvars);
         assert(consdata->quadvarterms[var2pos].var == bilinterm->var2);

         addBilinMcCormick(scip, bilinterm->coef,
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), ref[j],
            SCIPvarGetLbLocal(bilinterm->var2), SCIPvarGetUbLocal(bilinterm->var2), ref[var2pos],
            violside == SCIP_SIDETYPE_LEFT, &coef[j], &coef[var2pos], &constant, success);
      }
   }

   if( !*success )
   {
      SCIPdebugMessage("no success to find estimator for nonconvex <%s>\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   if( violside == SCIP_SIDETYPE_LEFT )
   {
      *lhs = consdata->lhs - constant;
      *rhs = SCIPinfinity(scip);
   }
   else
   {
      *lhs = -SCIPinfinity(scip);
      *rhs = consdata->rhs - constant;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_side%d_estimation_%d", SCIPconsGetName(cons), violside, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex) or McCormick (if nonconvex) in a given reference point */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_SOL*             sol,                /**< point that we aim to separate, or NULL for LP solution */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            efficacy,           /**< buffer to store efficacy of row in reference solution, or NULL if not of interest */
   SCIP_Real             maxrange,           /**< maximal range allowed */
   SCIP_Bool             checkcurvmultivar,  /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   SCIP_Real             minefficacy         /**< minimal required efficacy (violation scaled by maximal absolute coefficient) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      islocal;
   char           cutname[SCIP_MAXSTRLEN];
   SCIP_Real*     lincoefs;
   SCIP_Real*     coef;
   SCIP_Real      lhs;
   SCIP_Real      rhs;
   SCIP_Bool      success;
   SCIP_Real      mincoef;
   SCIP_Real      maxcoef;
   SCIP_Real      lincoefsmax;
   SCIP_Real      lincoefsmin;
   SCIP_Real      viol;
   SCIP_VAR*      var;
   int            j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref != NULL);
   assert(row != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(violside != SCIP_SIDETYPE_LEFT  || !SCIPisInfinity(scip, -consdata->lhs));
   assert(violside != SCIP_SIDETYPE_RIGHT || !SCIPisInfinity(scip,  consdata->rhs));

   *row = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &coef, consdata->nquadvars) );
   lincoefs = consdata->lincoefs;
   lincoefsmax = consdata->lincoefsmax;
   lincoefsmin = consdata->lincoefsmin;
   islocal = SCIPconsIsLocal(cons);
   success = FALSE;
   lhs = -SCIPinfinity(scip);
   rhs = SCIPinfinity(scip);
   cutname[0] = '\0';

   /* if constraint function is factorable, then try to use factorable form to generate cut */
   if( consdata->factorleft != NULL )
   {
      if( consdata->nlinvars == 0 )
      {
         SCIP_CALL( generateCutFactorable(scip, cons, violside, ref, coef, &lhs, &rhs, &islocal, &success, cutname) );
      }
      else if( sol != NULL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         int i;

         /* generateCutLTI needs reference values also for the linear variables, which we only have if sol is given or LP has been solved */
         SCIP_CALL( generateCutLTI(scip, cons, violside, ref, sol, &lincoefs, coef, &lhs, &rhs, &islocal, &success, cutname) );

         /* in case of LTI cuts, we have to recompute the min and max of lincoefs, since they may have been modified */
         for( i = 0; i < consdata->nlinvars; ++i )
         {
            if( REALABS(lincoefs[i]) > lincoefsmax )
               lincoefsmax = REALABS(lincoefs[i]);
            if( REALABS(lincoefs[i]) < lincoefsmin )
               lincoefsmin = REALABS(lincoefs[i]);
         }
      }
   }

   /* if constraint is not factorable or failed to generate cut, try default method */
   if( !success )
   {
      SCIP_CALL( checkCurvature(scip, cons, checkcurvmultivar) );

      if( (violside == SCIP_SIDETYPE_LEFT && consdata->isconcave) || (violside == SCIP_SIDETYPE_RIGHT && consdata->isconvex) )
      {
         SCIP_CALL( generateCutConvex(scip, cons, violside, ref, coef, &lhs, &rhs, &islocal, &success, cutname) );
      }
      else
      {
         SCIP_CALL( generateCutNonConvex(scip, cons, violside, ref, coef, &lhs, &rhs, &islocal, &success, cutname) );
      }
   }

   /* cut should be one-sided, if any found */
   assert(!success || SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs));

   /* check if range of cut coefficients is ok
    * compute cut activity and violation in sol
    */
   mincoef = 0.0; /* only for lint */
   maxcoef = 0.0; /* only for compiler */
   viol = 0.0;    /* only for compiler */
   if( success )
   {
      SCIP_Real constant;
      SCIP_Real abscoef;
      int       mincoefidx;
      SCIP_Real refactivity;
      SCIP_Real refactivitylinpart;

      /* compute activity of linear part in sol, if required
       * it is required if we need to check or return cut efficacy and for some debug output below
       * round almost integral coefficients in integers, since this will happen when adding coefs to row (see comments below)
       */
      refactivitylinpart = 0.0;
#ifndef SCIP_DEBUG
      if( !SCIPisInfinity(scip, -minefficacy) || efficacy != NULL )
#endif
         for( j = 0; j < consdata->nlinvars; ++j )
            /* Loose variable have the best bound as LP solution value.
             * HOWEVER, they become column variables when they are added to a row (via SCIPaddVarsToRow below).
             * When this happens, their LP solution value changes to 0.0!
             * So when calculating the row activity, we treat loose variable as if they were already column variables.
             */
            if( SCIPvarGetStatus(consdata->linvars[j]) != SCIP_VARSTATUS_LOOSE )
               refactivitylinpart += (SCIPisIntegral(scip, lincoefs[j]) ? SCIPround(scip, lincoefs[j]) : lincoefs[j]) * SCIPgetSolVal(scip, sol, consdata->linvars[j]);

      assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

      constant = 0.0;
      do
      {
         refactivity = refactivitylinpart;
         mincoefidx = -1;
         mincoef = lincoefsmin;
         maxcoef = lincoefsmax;

         for( j = 0; j < consdata->nquadvars; ++j )
         {
            /* coefficients smaller than epsilon are rounded to 0.0 when added to row, this can be problematic if variable value is very large (bad numerics)
             * in this case, set coef to 0.0 here, but modify constant so that cut is still valid (if possible)
             * i.e., estimate coef[i]*x >= coef[i] * bound(x) or coef[i]*x <= coef[i] * bound(x), depending on whether lhs or rhs is finite
             * if required bound of x is not finite, then do nothing
             */
            if( coef[j] != 0.0 && SCIPisZero(scip, coef[j]) )
            {
               SCIP_Real xbnd;

               var = consdata->quadvarterms[j].var;
               if( !SCIPisInfinity(scip, rhs) )
                  if( islocal )
                     xbnd = coef[j] > 0.0 ? SCIPvarGetLbLocal(var)  : SCIPvarGetUbLocal(var);
                  else
                     xbnd = coef[j] > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
               else
                  if( islocal )
                     xbnd = coef[j] > 0.0 ? SCIPvarGetUbLocal(var)  : SCIPvarGetLbLocal(var);
                  else
                     xbnd = coef[j] > 0.0 ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);

               if( !SCIPisInfinity(scip, REALABS(xbnd)) )
               {
                  SCIPdebugMessage("var <%s> [%g,%g] has tiny coef %g, replace coefficient by constant %g\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), coef[j], coef[j] * xbnd);
                  constant += coef[j] * xbnd;
                  coef[j] = 0.0;
               }

               continue;
            }

            /* coefficients very close to integral values are rounded to integers when added to LP
             * we do this here already to be sure to compute an accurate activity/violation/efficacy
             * (otherwise, if variables are at huge values, cut may look efficient here, but not anymore after rounding coefs)
             */
            if( SCIPisIntegral(scip, coef[j]) )
               coef[j] = SCIPround(scip, coef[j]);

            if( coef[j] == 0.0 )
               continue;

            /* As above: When calculating the row activity, we treat loose variable as if they were already column variables. */
            if( SCIPvarGetStatus(consdata->quadvarterms[j].var) != SCIP_VARSTATUS_LOOSE )
               refactivity += coef[j] * SCIPgetSolVal(scip, sol, consdata->quadvarterms[j].var);

            abscoef = REALABS(coef[j]);
            if( abscoef < mincoef )
            {
               mincoef = abscoef;
               mincoefidx = j;
            }
            if( abscoef > maxcoef )
               maxcoef = abscoef;
         }

         if( maxcoef < mincoef )
         {
            /* if all coefficients are zero, then mincoef and maxcoef are still at their initial values
             * thus, skip cut generation if its boring
             */
            assert(maxcoef == 0.0);  /*lint !e777 */
            assert(mincoef == SCIPinfinity(scip));  /*lint !e777 */

            if( !SCIPisFeasPositive(scip, lhs) && !SCIPisFeasNegative(scip, rhs) )
            {
               SCIPdebugMessage("skip cut for constraint <%s> since all coefficients are zero and it's always satisfied\n", SCIPconsGetName(cons));
               success = FALSE;
            }
            else
            {
               /* cut will cutoff node */
            }

            break;
         }

         if( maxcoef / mincoef > maxrange  )
         {
            SCIPdebugMessage("cut coefficients for constraint <%s> have very large range: mincoef = %g maxcoef = %g\n", SCIPconsGetName(cons), mincoef, maxcoef);
            if( mincoefidx >= 0 )
            {
               var = consdata->quadvarterms[mincoefidx].var;
               /* try to eliminate coefficient with minimal absolute value by weakening cut and try again
                * since we use local bounds, we need to make the row local if they are different from their global counterpart
                */
               if( ((coef[mincoefidx] > 0.0 && !SCIPisInfinity(scip, rhs)) || (coef[mincoefidx] < 0.0 && !SCIPisInfinity(scip, -lhs))) &&
                  !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  SCIPdebugMessage("eliminate coefficient %g for <%s> [%g, %g]\n", coef[mincoefidx], SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
                  constant += coef[mincoefidx] * SCIPvarGetLbLocal(var);
                  coef[mincoefidx] = 0.0;
                  islocal |= SCIPisGT(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var));
                  continue;
               }
               else if( ((coef[mincoefidx] < 0.0 && !SCIPisInfinity(scip, rhs)) || (coef[mincoefidx] > 0.0 && !SCIPisInfinity(scip, -lhs))) &&
                  !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  SCIPdebugMessage("eliminate coefficient %g for <%s> [%g, %g]\n", coef[mincoefidx], SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
                  constant += coef[mincoefidx] * SCIPvarGetUbLocal(var);
                  coef[mincoefidx] = 0.0;
                  islocal |= SCIPisLT(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var));
                  continue;
               }
            }

            SCIPdebugMessage("skip cut\n");
            success = FALSE;
         }
         break;
      }
      while( TRUE );  /*lint !e506 */

      if( !SCIPisInfinity(scip, -lhs) )
      {
         lhs -= constant;
         viol = lhs - refactivity;
      }
      if( !SCIPisInfinity(scip,  rhs) )
      {
         rhs -= constant;
         viol = refactivity - rhs;
      }
   }

   if( success && SCIPisInfinity(scip, REALABS(lhs)) && SCIPisInfinity(scip, REALABS(rhs)) )
   {
      SCIPdebugMessage("skip cut for constraint <%s> because both sides are not finite: lhs = %g, rhs = %g\n", SCIPconsGetName(cons), lhs, rhs);
      success = FALSE;
   }

   /* check if reference point violates cut sufficiently
    * in difference to SCIPgetCutEfficacy, we scale by norm only if the norm is > 1.0
    * this avoid finding cuts efficient which are only very slightly violated
    * CPLEX does not seem to scale row coefficients up too
    * also we use infinity norm, since that seem to be the usual scaling strategy in LP solvers (equilibrium scaling)
    */
   if( success && !SCIPisInfinity(scip, -minefficacy) && viol / MAX(1.0, maxcoef) < minefficacy )
   {
      SCIPdebugMessage("skip cut for constraint <%s> because efficacy %g/%g too low (< %g)\n", SCIPconsGetName(cons), viol, MAX(1.0, maxcoef), minefficacy);
      success = FALSE;
   }

   /* generate row */
   if( success )
   {
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), cutname, lhs, rhs, islocal && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

      /* add coefficients from linear part */
      SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, lincoefs) );

      /* add coefficients from quadratic part */
      assert(consdata->sepaquadvars != NULL || consdata->nquadvars == 0);
      SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nquadvars, consdata->sepaquadvars, coef) );

      SCIPdebugMessage("found cut <%s>, lhs=%g, rhs=%g, mincoef=%g, maxcoef=%g, range=%g, nnz=%d, violation=%g, efficacy=%g\n",
         SCIProwGetName(*row), lhs, rhs,
         mincoef, maxcoef, maxcoef/mincoef,
         SCIProwGetNNonz(*row), viol, viol / MAX(1.0, maxcoef));  /*lint !e414 */

      if( efficacy != NULL )
      {
         *efficacy = viol / MAX(1.0, maxcoef);
         /* check that our computed efficacy is > feastol, iff efficacy computed by row is > feastol
          * (they should be equal, but due to numerics...)
          */
         assert(SCIPisFeasPositive(scip, *efficacy) == SCIPisFeasPositive(scip, -SCIPgetRowSolFeasibility(scip, *row, sol)/MAX(1.0,SCIPgetRowMaxCoef(scip, *row))));  /*lint !e666*/
      }
   }

   SCIPfreeBufferArray(scip, &coef);

   /* if coefficients for linear variables are different than those in constraint, then free array */
   if( lincoefs != consdata->lincoefs )
   {
      SCIPfreeBufferArray(scip, &lincoefs);
   }

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex) or McCormick (if nonconvex) in a solution
 */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution where to generate cut, or NULL if LP solution should be used */
   SCIP_SOL*             refsol,             /**< reference point where to generate cut, or NULL if sol should be used */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            efficacy,           /**< buffer to store efficacy of row in reference solution, or NULL if not of interest */
   SCIP_Real             maxrange,           /**< maximal range allowed */
   SCIP_Bool             checkcurvmultivar,  /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   SCIP_Real             minefficacy         /**< minimal required efficacy (violation scaled by maximal absolute coefficient) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*  var;
   SCIP_Real  lb;
   SCIP_Real  ub;
   SCIP_Real* ref;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( refsol == NULL )
      refsol = sol;

   /* get reference point */
   SCIP_CALL( SCIPallocBufferArray(scip, &ref, consdata->nquadvars) );
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      var = consdata->quadvarterms[j].var;
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      /* do not like variables at infinity */
      assert(!SCIPisInfinity(scip,  lb));
      assert(!SCIPisInfinity(scip, -ub));

      ref[j] = SCIPgetSolVal(scip, refsol, var);
      ref[j] = MIN(ub, MAX(lb, ref[j])); /* project value into bounds */
   }

   SCIP_CALL( generateCut(scip, cons, ref, sol, violside, row, efficacy, maxrange, checkcurvmultivar, minefficacy) );

   SCIPfreeBufferArray(scip, &ref);

   return SCIP_OKAY;
}

/** tries to find a cut that intersects with an unbounded ray of the LP
 * for convex functions, we do this by linearizing in the feasible solution of the LPI
 * for nonconvex functions, we just call generateCutSol with the unbounded solution as reference point */
static
SCIP_RETCODE generateCutUnboundedLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            rowrayprod,         /**< buffer to store product of ray with row coefficients, or NULL if not of interest */
   SCIP_Real             maxrange,           /**< maximal range allowed */
   SCIP_Bool             checkcurvmultivar   /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_VAR*  var;
   SCIP_Real* ref;
   SCIP_Real  matrixrayprod;
   SCIP_Real  linrayprod;
   SCIP_Real  quadrayprod;
   SCIP_Real  rayval;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *row = NULL;

   if( !SCIPhasPrimalRay(scip) )
   {
      SCIPdebugMessage("do not have primal ray, thus cannot resolve unboundedness\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( checkCurvature(scip, cons, checkcurvmultivar) );
   if( (!consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) ||
      (!consdata->isconcave && violside == SCIP_SIDETYPE_LEFT) )
   {
      /* if not convex, just call generateCut and hope it's getting something useful */
      SCIP_CALL( generateCutSol(scip, cons, NULL, NULL, violside, row, NULL, maxrange, FALSE, -SCIPinfinity(scip)) );

      /* compute product of cut coefficients with ray, if required */
      if( *row != NULL && rowrayprod != NULL )
      {
         *rowrayprod = 0.0;
         for( i = 0; i < SCIProwGetNNonz(*row); ++i )
         {
            assert(SCIProwGetCols(*row)[i] != NULL);
            var = SCIPcolGetVar(SCIProwGetCols(*row)[i]);
            assert(var != NULL);

            *rowrayprod += SCIProwGetVals(*row)[i] * SCIPgetPrimalRayVal(scip, var);
         }
      }

      return SCIP_OKAY;
   }

   /* we seek for a linearization of the quadratic function such that it intersects with the unbounded ray
    * that is, we need a reference point ref such that for the gradient g of xAx+bx in ref, we have
    *   <g, ray> > 0.0 if rhs is finite and <g, ray> < 0.0 if lhs is finite
    * Since g = 2*A*ref + b, we have <g, ray> = <2*A*ref + b, ray> = <ref, 2*A*ray> + <b,ray>
    * initially, for finite rhs, we set ref_i = 1.0 if (A*ray)_i > 0.0 and ref_i = -1.0 if (A*ray)_i < 0.0 (for finite lhs analog)
    * <ref, 2*A*ray> + <b,ray> is sufficiently larger 0.0, we call generateCut for this point, otherwise, we scale up ref
    */

   quadrayprod = 0.0; /* <ref, 2*A*ray> */
   linrayprod = 0.0;  /* <b, ray> */
   SCIP_CALL( SCIPallocBufferArray(scip, &ref, consdata->nquadvars) );
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvarterms[i].var;
      rayval = SCIPgetPrimalRayVal(scip, var);

      /* compute i-th entry of (2*A*ray) */
      matrixrayprod = 2.0 * consdata->quadvarterms[i].sqrcoef * rayval;
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[j]];
         matrixrayprod += bilinterm->coef * SCIPgetPrimalRayVal(scip, bilinterm->var1 == var ? bilinterm->var2 : bilinterm->var1);
      }

      if( SCIPisPositive(scip, matrixrayprod) )
         ref[i] = (violside == SCIP_SIDETYPE_RIGHT ?  1.0 : -1.0);
      else if( SCIPisNegative(scip, matrixrayprod) )
         ref[i] = (violside == SCIP_SIDETYPE_RIGHT ? -1.0 :  1.0);
      else
         ref[i] = 0.0;

      quadrayprod += matrixrayprod * ref[i];
      linrayprod += consdata->quadvarterms[i].lincoef * rayval;
   }
   assert((violside == SCIP_SIDETYPE_RIGHT && quadrayprod >= 0.0) || (violside == SCIP_SIDETYPE_LEFT && quadrayprod <= 0.0));

   if( SCIPisZero(scip, quadrayprod) )
   {
      SCIPdebugMessage("ray is zero along cons <%s>\n", SCIPconsGetName(cons));
      SCIPfreeBufferArray(scip, &ref);
      return SCIP_OKAY;
   }

   /* add linear part to linrayprod */
   for( i = 0; i < consdata->nlinvars; ++i )
      linrayprod += consdata->lincoefs[i] * SCIPgetPrimalRayVal(scip, consdata->linvars[i]);

   SCIPdebugMessage("initially have <b,ray> = %g and <ref, 2*A*ref> = %g\n", linrayprod, quadrayprod);

   /* we scale the refpoint up, such that <ref, 2*A*ray> >= -2*<b, ray> (rhs finite) or <ref, 2*A*ray> <= -2*<b, ray> (lhs finite), if <b,ray> is not zero
    * if <b,ray> is zero, then we scale refpoint up if |<ref, 2*A*ray>| < 1.0
    */
   if( (!SCIPisZero(scip, linrayprod) && violside == SCIP_SIDETYPE_RIGHT && quadrayprod < -2*linrayprod) ||
      ( !SCIPisZero(scip, linrayprod) && violside == SCIP_SIDETYPE_LEFT  && quadrayprod > -2*linrayprod) ||
      (SCIPisZero(scip, linrayprod) && REALABS(quadrayprod) < 1.0) )
   {
      SCIP_Real scale;

      if( !SCIPisZero(scip, linrayprod) )
         scale = 2*REALABS(linrayprod/quadrayprod);  /*lint !e795 */
      else
         scale = 1.0/REALABS(quadrayprod);

      SCIPdebugMessage("scale refpoint by %g\n", scale);
      for( i = 0; i < consdata->nquadvars; ++i )
         ref[i] *= scale;
      quadrayprod *= scale;
   }

   if( rowrayprod != NULL )
      *rowrayprod = quadrayprod + linrayprod;

   SCIPdebugMessage("calling generateCut, expecting ray product %g\n", quadrayprod + linrayprod);
   SCIP_CALL( generateCut(scip, cons, ref, NULL, violside, row, NULL, maxrange, FALSE, -SCIPinfinity(scip)) );

   SCIPfreeBufferArray(scip, &ref);

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 * 
 *  assumes that constraint violations have been computed 
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_Real          actminefficacy;
   SCIP_SIDETYPE      violside;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         /* we are not feasible anymore */
         if( *result == SCIP_FEASIBLE )
            *result = SCIP_DIDNOTFIND;

         violside = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT;

         /* actual minimal efficacy */
         actminefficacy = inenforcement && ((violside == SCIP_SIDETYPE_RIGHT && consdata->isconvex ) || (violside == SCIP_SIDETYPE_LEFT && consdata->isconcave)) ? SCIPfeastol(scip) : minefficacy;

         /* generate cut */
         if( sol == NULL && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            /* if the LP is unbounded, then we need a cut that cuts into the direction of a hopefully existing primal ray
             * that is, assume a ray r is given such that p + t*r is feasible for the LP for all t >= t_0 and some p
             * given a cut lhs <= <c,x> <= rhs, we check whether it imposes an upper bound on t and thus bounds the ray
             * this is given if rhs < infinity and <c,r> > 0, since we then enforce <c,p+t*r> = <c,p> + t<c,r> <= rhs, i.e., t <= (rhs - <c,p>)/<c,r>
             * similar, lhs > -infinity and <c,r> < 0 is good
             */
            SCIP_Real rayprod;
            SCIP_Real feasibility;
            SCIP_Real norm;

            rayprod = 0.0; /* for compiler */
            SCIP_CALL( generateCutUnboundedLP(scip, conss[c], violside, &row, &rayprod, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature) );

            if( row != NULL )
            {
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPisPositive(scip, rayprod) )
                  feasibility = -rayprod;
               else if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) && SCIPisNegative(scip, rayprod) )
                  feasibility =  rayprod;
               else
                  feasibility = 0.0;

               norm = SCIPgetRowMaxCoef(scip, row);
               if( norm > 1.0 )
                  efficacy = -feasibility / norm;
               else
                  efficacy = -feasibility;
            }
         }
         else
         {
            /* @todo If convex, can we easily move the refpoint closer to the feasible region to get a stronger cut?
             * E.g., use bisection on the line between LP solution and best primal (or LP interior)
             */
            SCIP_CALL( generateCutSol(scip, conss[c], sol, NULL, violside, &row, &efficacy, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature, actminefficacy) );
            /* @todo If generation failed not because of low efficacy, then probably because of numerical issues;
             * if the constraint is convex and we are desperate to get a cut, then we may try again with a better chosen reference point
             */
         }

         if( row == NULL ) /* failed to generate cut */
            continue;

         if( SCIPisGT(scip, efficacy, actminefficacy) )  /*lint !e644 */
         {
            /* cut cuts off solution */
            SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE /* forcecut */) );

            *result = SCIP_SEPARATED;

            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );

            SCIPdebugMessage("add cut with efficacy %g for constraint <%s> violated by %g\n", efficacy,
               SCIPconsGetName(conss[c]), consdata->lhsviol+consdata->rhsviol);

            /* mark row as not removable from LP for current node, if in enforcement */
            if( inenforcement && !conshdlrdata->enfocutsremovable )
               SCIPmarkRowNotRemovableLocal(scip, row);
         }
         if( bestefficacy != NULL && efficacy > *bestefficacy )
            *bestefficacy = efficacy;

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */ 
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 * if separatedlpsol is not NULL, then a cut that separates the LP solution is added to the sepastore and is forced to enter the LP
 * if separatedlpsol is not NULL, but cut does not separate the LP solution, then it is added to the cutpool only
 * if separatedlpsol is NULL, then cut is added to cutpool only
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found and added to LP, or NULL if adding to cutpool only */
   SCIP_Real             minefficacy         /**< minimal efficacy of a cut when checking for separation of LP solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool addedtolp;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613 */
         continue;

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      if( consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( generateCutSol(scip, conss[c], NULL, ref, SCIP_SIDETYPE_RIGHT, &row, NULL, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
      }
      else if( consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIP_CALL( generateCutSol(scip, conss[c], NULL, ref, SCIP_SIDETYPE_LEFT,  &row, NULL, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
      }
      else
         continue;

      if( row == NULL )
         continue;

      addedtolp = FALSE;

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         SCIP_Real feasibility;
         SCIP_Real norm;

         feasibility = SCIPgetRowLPFeasibility(scip, row);
         norm = SCIPgetRowMaxCoef(scip, row);

         if( -feasibility / MAX(1.0, norm) >= minefficacy )
         {
            *separatedlpsol = TRUE;
            addedtolp = TRUE;
            SCIP_CALL( SCIPaddCut(scip, NULL, row, TRUE) );
            SCIPdebugMessage("added linearization cut <%s> to LP, efficacy = %g\n", SCIProwGetName(row), -feasibility / MAX(1.0, norm));
         }
      }

      if( !SCIProwIsLocal(row) && !addedtolp )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
         SCIPdebugMessage("added linearization cut <%s> to cutpool\n", SCIProwGetName(row));
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_SOL*      sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMessage("caught new sol event %x from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0) );

   return SCIP_OKAY;
}

/** computes the infeasibilities of variables from the convexification gaps in the constraints and notifies the branching rule about them
 */
static
SCIP_RETCODE registerVariableInfeasibilities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA*     consdata;
   int                c;
   int                j;
   SCIP_Bool          xbinary;
   SCIP_Bool          ybinary;
   SCIP_Bool          xunbounded;
   SCIP_Bool          yunbounded;
   SCIP_VAR*          x;
   SCIP_VAR*          y;
   SCIP_Real          xlb;
   SCIP_Real          xub;
   SCIP_Real          xval;
   SCIP_Real          ylb;
   SCIP_Real          yub;
   SCIP_Real          yval;
   SCIP_Real          gap;
   SCIP_Real          coef_;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;
   yval = SCIP_INVALID;
   xval = SCIP_INVALID;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->nquadvars )
         continue;

      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || consdata->isconcave) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || consdata->isconvex ) )
         continue;
      SCIPdebugMessage("cons %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);

      /* square terms */
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         x = consdata->quadvarterms[j].var;
         if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && consdata->quadvarterms[j].sqrcoef < 0) ||
            ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && consdata->quadvarterms[j].sqrcoef > 0) )
         {
            xlb = SCIPvarGetLbLocal(x);
            xub = SCIPvarGetUbLocal(x);
            if( SCIPisRelEQ(scip, xlb, xub) )
            {
               SCIPdebugMessage("ignore fixed variable <%s>[%g, %g], diff %g\n", SCIPvarGetName(x), xlb, xub, xub-xlb);
               continue;
            }

            xval = SCIPgetSolVal(scip, NULL, x);

            /* if variable is at bounds, then no need to branch, since secant is exact there */
            if( SCIPisLE(scip, xval, xlb) || SCIPisGE(scip, xval, xub) )
               continue;

            if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
               gap = SCIPinfinity(scip);
            else
               gap = (xval-xlb)*(xub-xval)/(1+2*ABS(xval));
            assert(!SCIPisFeasNegative(scip, gap));
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, MAX(gap, 0.0), SCIP_INVALID) );
            ++*nnotify;
         }
      }

      /* bilinear terms */
      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         /* if any of the variables if fixed, then it actually behaves like a linear term, so we don't need to branch on it */
         x = consdata->bilinterms[j].var1;
         xlb = SCIPvarGetLbLocal(x);
         xub = SCIPvarGetUbLocal(x);
         if( SCIPisRelEQ(scip, xlb, xub) )
            continue;

         y = consdata->bilinterms[j].var2;
         ylb = SCIPvarGetLbLocal(y);
         yub = SCIPvarGetUbLocal(y);
         if( SCIPisRelEQ(scip, ylb, yub) )
            continue;

         xunbounded = SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub);
         yunbounded = SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub);

         /* compute gap, if both variable are bounded */
         gap = SCIPinfinity(scip);
         if( !xunbounded && !yunbounded )
         {
            xval = SCIPgetSolVal(scip, NULL, x);
            yval = SCIPgetSolVal(scip, NULL, y);

            /* if both variables are at one of its bounds, then no need to branch, since McCormick is exact there */
            if( (SCIPisLE(scip, xval, xlb) || SCIPisGE(scip, xval, xub)) &&
               ( SCIPisLE(scip, yval, ylb) || SCIPisGE(scip, yval, yub)) )
               continue;

            xval = MAX(xlb, MIN(xval, xub));
            yval = MAX(ylb, MIN(yval, yub));

            coef_ = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? -consdata->bilinterms[j].coef : consdata->bilinterms[j].coef;
            if( coef_ > 0.0 )
            {
               if( (xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb )
                  gap = (xval*yval - xlb*yval - ylb*xval + xlb*ylb) / (1+sqrt(xval*xval + yval*yval));
               else
                  gap = (xval*yval - xval*yub - yval*xub + xub*yub) / (1+sqrt(xval*xval + yval*yval));
            }
            else
            { /* coef_ < 0 */
               if( (xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub )
                  gap = -(xval*yval - xval*ylb - yval*xub + xub*ylb) / (1+sqrt(xval*xval + yval*yval));
               else
                  gap = -(xval*yval - xval*yub - yval*xlb + xlb*yub) / (1+sqrt(xval*xval + yval*yval));
            }

            assert(!SCIPisNegative(scip, gap / MAX3(REALABS(xval), REALABS(yval), 1.0)));  /*lint !e666*/
            if( gap < 0.0 )
               gap = 0.0;
         }

         /* if one of the variables is binary or integral with domain width 1, then branching on this makes the term linear, so prefer this */
         xbinary = SCIPvarIsBinary(x) || (SCIPvarIsIntegral(x) && xub - xlb < 1.5);
         ybinary = SCIPvarIsBinary(y) || (SCIPvarIsIntegral(y) && yub - ylb < 1.5);
         if( xbinary )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( ybinary )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( xbinary || ybinary )
            continue;

         /* if one of the variables is unbounded, then branch on it first */
         if( xunbounded )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( yunbounded )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( xunbounded || yunbounded )
            continue;

         /* if both variables are integral, prefer the one with the smaller domain, so variable gets fixed soon
          * does not seem to work well on tln instances, so disable for now and may look at it later again
          */
#ifdef BRANCHTOLINEARITY
         if( SCIPvarIsIntegral(x) && SCIPvarIsIntegral(y) )
         {
            if( SCIPisLT(scip, xub-xlb, yub-ylb) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
               ++*nnotify;
               continue;
            }
            if( SCIPisGT(scip, xub-xlb, yub-ylb) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
               ++*nnotify;
               continue;
            }
         }
#endif

         /* in the regular case, suggest those variables which are not at its bounds for branching
          * this is, because after branching both variables will be one the bounds, and McCormick will be exact then */
         if( !SCIPisLE(scip, xval, xlb) && !SCIPisGE(scip, xval, xub) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( !SCIPisLE(scip, yval, ylb) && !SCIPisGE(scip, yval, yub) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMessage("registered %d branching candidates\n", *nnotify);

   return SCIP_OKAY;
}

/** registers a quadratic variable from a violated constraint as branching candidate that has a large absolute value in the LP relaxation */
static
SCIP_RETCODE registerLargeLPValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int                 i;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   *brvar = NULL;
   brvarval = -1.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         /* do not propose fixed variables */
         if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
            continue;
         val = SCIPgetSolVal(scip, NULL, consdata->quadvarterms[i].var);
         if( ABS(val) > brvarval )
         {
            brvarval = ABS(val);
            *brvar = consdata->quadvarterms[i].var;
         }
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** replaces violated quadratic constraints where all quadratic variables are fixed by linear constraints */
static
SCIP_RETCODE replaceByLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            addedcons,          /**< buffer to store whether a linear constraint was added */
   SCIP_Bool*            reduceddom,         /**< whether a domain has been reduced */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONS*          cons;
   SCIP_CONSDATA*      consdata;
   SCIP_RESULT         checkresult;
   SCIP_Real           constant;
   SCIP_Real           val1;
   SCIP_Real           val2;
   int                 i;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(addedcons != NULL);
   assert(reduceddom != NULL);
   assert(infeasible != NULL);

   *addedcons = FALSE;
   *reduceddom = FALSE;
   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      constant = 0.0;

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         /* variables should be fixed if constraint is violated */
         assert(SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var)));

         val1 = (SCIPvarGetUbLocal(consdata->quadvarterms[i].var) + SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) / 2.0;
         constant += (consdata->quadvarterms[i].lincoef + consdata->quadvarterms[i].sqrcoef * val1) * val1;
      }

      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         val1 = (SCIPvarGetUbLocal(consdata->bilinterms[i].var1) + SCIPvarGetLbLocal(consdata->bilinterms[i].var1)) / 2.0;
         val2 = (SCIPvarGetUbLocal(consdata->bilinterms[i].var2) + SCIPvarGetLbLocal(consdata->bilinterms[i].var2)) / 2.0;
         constant += consdata->bilinterms[i].coef * val1 * val2;
      }

      /* check if we have a bound change */
      if ( consdata->nlinvars == 1 )
      {
         SCIP_Bool tightened;
         SCIP_Real coef;
         SCIP_Real lhs;
         SCIP_Real rhs;

         coef = *consdata->lincoefs;

         /* compute lhs/rhs */
         if ( SCIPisInfinity(scip, -consdata->lhs) )
            lhs = -SCIPinfinity(scip);
         else
            lhs = consdata->lhs - constant;

         if ( SCIPisInfinity(scip, consdata->rhs) )
            rhs = SCIPinfinity(scip);
         else
            rhs = consdata->rhs - constant;

         SCIPdebugMessage("Linear constraint with one variable: %g <= %g <%s> <= %g\n", lhs, coef, SCIPvarGetName(*consdata->linvars), rhs);

         /* possibly correct lhs/rhs */
         assert( ! SCIPisZero(scip, coef) );
         if ( coef >= 0.0 )
         {
            if ( ! SCIPisInfinity(scip, -lhs) )
               lhs /= coef;
            if ( ! SCIPisInfinity(scip, rhs) )
               rhs /= coef;
         }
         else
         {
            SCIP_Real h;
            h = rhs;
            if ( ! SCIPisInfinity(scip, -lhs) )
               rhs = lhs/coef;
            else
               rhs = SCIPinfinity(scip);

            if ( ! SCIPisInfinity(scip, h) )
               lhs = h/coef;
            else
               lhs = -SCIPinfinity(scip);
         }
         SCIPdebugMessage("Linear constraint is a bound: %g <= <%s> <= %g\n", lhs, SCIPvarGetName(*consdata->linvars), rhs);

         if ( ! SCIPisInfinity(scip, -lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, *consdata->linvars, lhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMessage("Lower bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMessage("Lower boundx changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }

         if ( ! SCIPisInfinity(scip, rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, *consdata->linvars, rhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMessage("Upper bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMessage("Upper bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }
      }
      else
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPconsGetName(conss[c]),
               consdata->nlinvars, consdata->linvars, consdata->lincoefs,
               (SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : (consdata->lhs - constant)),
               (SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : (consdata->rhs - constant)),
               SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
               SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  TRUE,
               SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
               SCIPconsIsStickingAtNode(conss[c])) );

         SCIPdebugMessage("replace quadratic constraint <%s> by linear constraint after all quadratic vars have been fixed\n", SCIPconsGetName(conss[c]) );
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPcheckCons(scip, cons, NULL, FALSE, FALSE, FALSE, &checkresult) );

         if( checkresult != SCIP_INFEASIBLE )
         {
            SCIPdebugMessage("linear constraint is feasible, thus do not add\n");
         }
         else
         {
            SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
            *addedcons = TRUE;
         }
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
   }

   return SCIP_OKAY;
}

/* tightens a lower bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new lower bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(intervalinfty > 0.0);
   assert(bnd > -intervalinfty);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   /* new bound is no improvement */
   if( SCIPisHugeValue(scip, -bnd) || SCIPisLE(scip, bnd, SCIPvarGetLbLocal(var)) )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   /* new lower bound is very low (between -intervalinfty and -SCIPinfinity()) */
   if( SCIPisInfinity(scip, -bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarLb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarLb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMessage("%s found constraint <%s> infeasible due to tightened lower bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing" : "", SCIPconsGetName(cons), bnd, SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMessage("%s tightened lower bound of variable <%s> in constraint <%s> to %g\n", SCIPinProbing(scip) ? "in probing" : "", SCIPvarGetName(var), SCIPconsGetName(cons), bnd);
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/* tightens an upper bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new upper bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(intervalinfty > 0.0);
   assert(bnd < intervalinfty);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   /* new bound is no improvement */
   if( SCIPisHugeValue(scip, bnd) || SCIPisGE(scip, bnd, SCIPvarGetUbLocal(var)) )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, -bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   /* new upper bound is very high (between SCIPinfinity() and intervalinfty) */
   if( SCIPisInfinity(scip, bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarUb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarUb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMessage("%s found constraint <%s> infeasible due to tightened upper bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing" : "", SCIPconsGetName(cons), bnd, SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMessage("%s tightened upper bound of variable <%s> in constraint <%s> to %g\n", SCIPinProbing(scip) ? "in probing" : "", SCIPvarGetName(var), SCIPconsGetName(cons), bnd);
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** solves a quadratic equation \f$ a x^2 + b x \in rhs \f$ (with b an interval) and reduces bounds on x or deduces infeasibility if possible
 */
static
SCIP_RETCODE propagateBoundsQuadVar(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which bounds with might tighten */
   SCIP_Real             a,                  /**< coefficient in square term */
   SCIP_INTERVAL         b,                  /**< coefficient in linear term */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< result of propagation */
   int*                  nchgbds             /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   /* compute solution of a*x^2 + b*x \in rhs */
   if( a == 0.0 && SCIPintervalGetInf(b) == 0.0 && SCIPintervalGetSup(b) == 0.0 )
   { /* relatively easy case: 0.0 \in rhs, thus check if infeasible or just redundant */
      if( SCIPintervalGetInf(rhs) > 0.0 || SCIPintervalGetSup(rhs) < 0.0 )
      {
         SCIPdebugMessage("found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *result = SCIP_CUTOFF;
      }
      return SCIP_OKAY;
   }
   else if( SCIPvarGetLbLocal(var) >= 0.0 )
   { 
      SCIP_INTERVAL a_;

      /* need only positive solutions */
      SCIPintervalSet(&a_, a);
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &newrange, a_, b, rhs);
   }
   else if( SCIPvarGetUbLocal(var) <= 0.0 )
   {
      /* need only negative solutions */
      SCIP_INTERVAL a_;
      SCIP_INTERVAL tmp;
      SCIPintervalSet(&a_, a);
      SCIPintervalSetBounds(&tmp, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &tmp, a_, tmp, rhs);
      if( SCIPintervalIsEmpty(tmp) )
      {
         SCIPdebugMessage("found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      SCIPintervalSetBounds(&newrange, -SCIPintervalGetSup(tmp), -SCIPintervalGetInf(tmp));
   }
   else
   {
      /* need both positive and negative solution */
      SCIP_INTERVAL a_;
      SCIPintervalSet(&a_, a);
      SCIPintervalSolveUnivariateQuadExpression(intervalinfty, &newrange, a_, b, rhs);
   }

   /* SCIPdebugMessage("%g x^2 + [%g, %g] x in [%g, %g] -> [%g, %g]\n", a, b.inf, b.sup, rhs.inf, rhs.sup, newrange.inf, newrange.sup); */

   if( SCIPisInfinity(scip, SCIPintervalGetInf(newrange)) || SCIPisInfinity(scip, -SCIPintervalGetSup(newrange)) )
   { /* domain outside [-infty, +infty] -> declare node infeasible */
      SCIPdebugMessage("found <%s> infeasible because propagated domain of quadratic variable <%s> is outside of (-infty, +infty)\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if( SCIPintervalIsEmpty(newrange) )
   {
      SCIPdebugMessage("found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(newrange)) )
   {
      SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, SCIPintervalGetInf(newrange), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip,  SCIPintervalGetSup(newrange)) )
   {
      SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, SCIPintervalGetSup(newrange), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/* the new version below computes potentially tighter bounds, but also always adds a small safety area since it is not implemented roundingsafe,
 * this may be a reason why it gives worse results on one of two instances
 * further, I have only very few instances where one can expect a difference
 */
#ifndef PROPBILINNEW
/** tries to deduce domain reductions for x in xsqrcoef x^2 + xlincoef x + ysqrcoef y^2 + ylincoef y + bilincoef x y \\in rhs
 * NOTE that domain reductions for y are not deduced 
 */
static
SCIP_RETCODE propagateBoundsBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint, where the bilinear term belongs to */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_Real             xsqrcoef,           /**< square coefficient of x */
   SCIP_Real             xlincoef,           /**< linear coefficient of x */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             ysqrcoef,           /**< square coefficient of y */
   SCIP_Real             ylincoef,           /**< linear coefficient of y */
   SCIP_Real             bilincoef,          /**< bilinear coefficient of x*y */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< pointer to store result of domain propagation */
   int*                  nchgbds             /**< counter to increment if domain reductions are found */
   )
{
   SCIP_INTERVAL myrhs;
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL lincoef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);
   assert(bilincoef != 0.0);

   if( SCIPintervalIsEntire(intervalinfty, rhs) )
      return SCIP_OKAY;

   /* try to find domain reductions for x */
   SCIPintervalSetBounds(&varbnds, MIN(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y)), MAX(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y)));  /*lint !e666 */

   /* put ysqrcoef*y^2 + ylincoef * y into rhs */
   if( SCIPintervalGetSup(rhs) >= intervalinfty )
   {
      /* if rhs is unbounded by above, it is sufficient to get an upper bound on ysqrcoef*y^2 + ylincoef * y */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real      tmp;

      SCIPintervalSet(&lincoef, ylincoef);
      tmp = SCIPintervalQuadUpperBound(intervalinfty, ysqrcoef, lincoef, varbnds);
      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();
      SCIPintervalSetBounds(&myrhs, SCIPintervalGetInf(rhs) - tmp, intervalinfty);
      SCIPintervalSetRoundingMode(roundmode);
   }
   else if( SCIPintervalGetInf(rhs) <= -intervalinfty )
   {
      /* if rhs is unbounded by below, it is sufficient to get a  lower bound on ysqrcoef*y^2 + ylincoef * y */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real      tmp;

      SCIPintervalSet(&lincoef, -ylincoef);
      tmp = -SCIPintervalQuadUpperBound(intervalinfty, -ysqrcoef, lincoef, varbnds);
      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();
      SCIPintervalSetBounds(&myrhs, -intervalinfty, SCIPintervalGetSup(rhs) - tmp);
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      /* if rhs is bounded, we need both bounds on ysqrcoef*y^2 + ylincoef * y */
      SCIP_INTERVAL tmp;

      SCIPintervalSet(&lincoef, ylincoef);
      SCIPintervalQuad(intervalinfty, &tmp, ysqrcoef, lincoef, varbnds);
      SCIPintervalSub(intervalinfty, &myrhs, rhs, tmp);
   }

   /* create equation xsqrcoef * x^2 + (xlincoef + bilincoef * [ylb, yub]) * x \in myrhs */
   SCIPintervalMulScalar(intervalinfty, &lincoef, varbnds, bilincoef);
   SCIPintervalAddScalar(intervalinfty, &lincoef, lincoef, xlincoef);

   /* propagate bounds on x */
   SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, x, xsqrcoef, lincoef, myrhs, result, nchgbds) );

   return SCIP_OKAY;
}
#else
/** tries to deduce domain reductions for x in xsqrcoef x^2 + xlincoef x + ysqrcoef y^2 + ylincoef y + bilincoef x y \\in rhs
 * NOTE that domain reductions for y are not deduced
 */
static
SCIP_RETCODE propagateBoundsBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint, where the bilinear term belongs to */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_Real             xsqrcoef,           /**< square coefficient of x */
   SCIP_Real             xlincoef,           /**< linear coefficient of x */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             ysqrcoef,           /**< square coefficient of y */
   SCIP_Real             ylincoef,           /**< linear coefficient of y */
   SCIP_Real             bilincoef,          /**< bilinear coefficient of x*y */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< pointer to store result of domain propagation */
   int*                  nchgbds             /**< counter to increment if domain reductions are found */
   )
{
   SCIP_INTERVAL xbnds;
   SCIP_INTERVAL ybnds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);
   assert(bilincoef != 0.0);

   if( SCIPintervalIsEntire(intervalinfty, rhs) )
      return SCIP_OKAY;

   SCIPintervalSetBounds(&xbnds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x))),   /*lint !e666*/
      +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x))));  /*lint !e666*/
   SCIPintervalSetBounds(&ybnds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y))),   /*lint !e666*/
      +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y))));  /*lint !e666*/

   /* try to find domain reductions for x */
   SCIPintervalSolveBivariateQuadExpressionAllScalar(intervalinfty, &xbnds, xsqrcoef, ysqrcoef, bilincoef, xlincoef, ylincoef, rhs, xbnds, ybnds);

   if( SCIPintervalIsEmpty(xbnds) )
   {
      SCIPdebugMessage("found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(x));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(xbnds)) )
   {
      SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, x, SCIPintervalGetInf(xbnds), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip,  SCIPintervalGetSup(xbnds)) )
   {
      SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, x, SCIPintervalGetSup(xbnds), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}
#endif

/** computes the minimal and maximal activity for the quadratic part in a constraint data
 *  only sums up terms that contribute finite values
 *  gives the number of terms that contribute infinite values
 *  only computes those activities where the corresponding side of the constraint is finite
 */
static
void propagateBoundsGetQuadActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_Real*            minquadactivity,    /**< minimal activity of quadratic variable terms where only terms with finite minimal activity contribute */
   SCIP_Real*            maxquadactivity,    /**< maximal activity of quadratic variable terms where only terms with finite maximal activity contribute */
   int*                  minactivityinf,     /**< number of quadratic variables that contribute -infinity to minimal activity */
   int*                  maxactivityinf,     /**< number of quadratic variables that contribute +infinity to maximal activity */
   SCIP_INTERVAL*        quadactcontr        /**< contribution of each quadratic variables to quadactivity */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   int       j;
   int       k;
   SCIP_INTERVAL tmp;
   SCIP_Real bnd;
   SCIP_INTERVAL xrng;
   SCIP_INTERVAL lincoef;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(minquadactivity != NULL);
   assert(maxquadactivity != NULL);
   assert(minactivityinf != NULL);
   assert(maxactivityinf != NULL);
   assert(quadactcontr != NULL);

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   *minquadactivity = SCIPisInfinity(scip,  consdata->rhs) ? -intervalinfty : 0.0;
   *maxquadactivity = SCIPisInfinity(scip, -consdata->lhs) ?  intervalinfty : 0.0;

   *minactivityinf = 0;
   *maxactivityinf = 0;

   if( consdata->nquadvars == 0 )
   {
      SCIPintervalSet(&consdata->quadactivitybounds, 0.0); 
      return;
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* there should be no quadratic variables fixed at -/+ infinity due to our locks */
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->quadvarterms[i].var)));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->quadvarterms[i].var)));

      SCIPintervalSetBounds(&quadactcontr[i], -intervalinfty, intervalinfty);

      SCIPintervalSetBounds(&xrng,
         -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var))),
         +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var))));

      SCIPintervalSet(&lincoef, consdata->quadvarterms[i].lincoef);
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         k = consdata->quadvarterms[i].adjbilin[j];
         if( consdata->bilinterms[k].var1 != consdata->quadvarterms[i].var )
            continue; /* handle this term later */

         SCIPintervalSetBounds(&tmp, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
            +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
         SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
      }

      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         /* compute maximal activity only if there is a finite left hand side */
         bnd = SCIPintervalQuadUpperBound(intervalinfty, consdata->quadvarterms[i].sqrcoef, lincoef, xrng);
         if( SCIPisInfinity(scip,  bnd) )
         {
            ++*maxactivityinf;
         }
         else if( SCIPisInfinity(scip, -bnd) )
         {
            /* if maximal activity is below value for -infinity, let's take -1e10 as upper bound on maximal activity
             * @todo Something better?
             */
            bnd = -sqrt(SCIPinfinity(scip));
            *maxquadactivity += bnd;
            quadactcontr[i].sup = bnd;
         }
         else
         {
            prevroundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeUpwards();
            *maxquadactivity += bnd;
            SCIPintervalSetRoundingMode(prevroundmode);
            quadactcontr[i].sup = bnd;
         }
      }

      if( !SCIPisInfinity(scip,  consdata->rhs) )
      {
         /* compute minimal activity only if there is a finite right hand side */
         SCIPintervalSetBounds(&lincoef, -SCIPintervalGetSup(lincoef), -SCIPintervalGetInf(lincoef));
         bnd = -SCIPintervalQuadUpperBound(intervalinfty, -consdata->quadvarterms[i].sqrcoef, lincoef, xrng);

         if( SCIPisInfinity(scip, -bnd) )
         {
            ++*minactivityinf;
         }
         else if( SCIPisInfinity(scip, bnd) )
         {
            /* if minimal activity is above value for infinity, let's take 1e10 as lower bound on minimal activity
             * @todo Something better?
             */
            bnd = sqrt(SCIPinfinity(scip));
            *minquadactivity += bnd;
            quadactcontr[i].inf = bnd;
         }
         else
         {
            prevroundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();
            *minquadactivity += bnd;
            SCIPintervalSetRoundingMode(prevroundmode);
            quadactcontr[i].inf = bnd;
         }
      }

   }

   SCIPintervalSetBounds(&consdata->quadactivitybounds,
      (*minactivityinf > 0 ? -intervalinfty : *minquadactivity),
      (*maxactivityinf > 0 ?  intervalinfty : *maxquadactivity));
   assert(!SCIPintervalIsEmpty(consdata->quadactivitybounds));
}

/** propagates bounds on a quadratic constraint */
static
SCIP_RETCODE propagateBoundsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   SCIP_Bool*            redundant           /**< buffer where to store whether constraint has been found to be redundant */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA*     consdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_INTERVAL      consactivity;  /* activity of linear plus quadratic part */
   SCIP_Real          intervalinfty; /* infinity used for interval computation */  
   SCIP_Real          minquadactivity; /* lower bound on finite activities of quadratic part */
   SCIP_Real          maxquadactivity; /* upper bound on finite activities of quadratic part */
   int                quadminactinf; /* number of quadratic variables that contribute -infinity to minimal activity of quadratic term */
   int                quadmaxactinf; /* number of quadratic variables that contribute +infinity to maximal activity of quadratic term */
   SCIP_INTERVAL*     quadactcontr;  /* contribution of each quadratic variable term to quadactivity */

   SCIP_VAR*          var;
   SCIP_INTERVAL      rhs;           /* right hand side of quadratic equation */
   SCIP_INTERVAL      tmp;
   SCIP_ROUNDMODE     roundmode;
   SCIP_Real          bnd;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN;
   *redundant = FALSE;

   if( consdata->ispropagated )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);

   quadactcontr = NULL;
   quadminactinf = -1;
   quadmaxactinf = -1;

   SCIPdebugMessage("start domain propagation for constraint <%s>\n", SCIPconsGetName(cons));

   consdata->ispropagated = TRUE;

   /* make sure we have activity of linear term and that they are consistent */
   consdataUpdateLinearActivity(scip, consdata, intervalinfty);
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->minlinactivityinf >= 0);
   assert(consdata->maxlinactivityinf >= 0);

   /* sort quadratic variable terms, in case we need to search for terms occuring in bilinear terms later
    * we sort here already, since we rely on a constant variable order during this method
    */
   if( consdata->nbilinterms > 0 )
   {
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );
   }

   /* compute activity of quad term part, if not up to date
    * in that case, we also collect the contribution of each quad var term for later */
   if( SCIPintervalIsEmpty(consdata->quadactivitybounds) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &quadactcontr, consdata->nquadvars) );
      propagateBoundsGetQuadActivity(scip, consdata, intervalinfty, &minquadactivity, &maxquadactivity, &quadminactinf, &quadmaxactinf, quadactcontr);
      assert(!SCIPintervalIsEmpty(consdata->quadactivitybounds));
   }

   SCIPdebugMessage("linear activity: [%g, %g]   quadratic activity: [%g, %g]\n",
      (consdata->minlinactivityinf > 0 ? -SCIPinfinity(scip) : consdata->minlinactivity),
      (consdata->maxlinactivityinf > 0 ?  SCIPinfinity(scip) : consdata->maxlinactivity),
      consdata->quadactivitybounds.inf, consdata->quadactivitybounds.sup);

   /* extend constraint bounds by epsilon to avoid some numerical difficulties */
   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -consdata->lhs+SCIPepsilon(scip)),
      +infty2infty(SCIPinfinity(scip), intervalinfty,  consdata->rhs+SCIPepsilon(scip)));

   /* check redundancy and infeasibility */
   SCIPintervalSetBounds(&consactivity, consdata->minlinactivityinf > 0 ? -intervalinfty : consdata->minlinactivity, consdata->maxlinactivityinf > 0 ? intervalinfty : consdata->maxlinactivity);
   SCIPintervalAdd(intervalinfty, &consactivity, consactivity, consdata->quadactivitybounds);
   if( SCIPintervalIsSubsetEQ(intervalinfty, consactivity, consbounds) )
   {
      SCIPdebugMessage("found constraint <%s> to be redundant: sides: [%g, %g], activity: [%g, %g]\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity));
      *redundant = TRUE;
      goto CLEANUP;
   }

   /* was SCIPintervalAreDisjoint(consbounds, consactivity), but that would allow violations up to eps only
    * we need to decide feasibility w.r.t. feastol (but still want to propagate w.r.t. eps)
    */
   if( SCIPisFeasGT(scip, consbounds.inf, consactivity.sup) || SCIPisFeasLT(scip, consbounds.sup, consactivity.inf) )
   {
      SCIPdebugMessage("found constraint <%s> to be infeasible; sides: [%g, %g], activity: [%g, %g], infeas: %g\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity),
         MAX(consdata->lhs - SCIPintervalGetSup(consactivity), SCIPintervalGetInf(consactivity) - consdata->rhs));
      *result = SCIP_CUTOFF;
      goto CLEANUP;
   }

   /* propagate linear part \in rhs = consbounds - quadactivity (use the one from consdata, since that includes infinities) */
   SCIPintervalSub(intervalinfty, &rhs, consbounds, consdata->quadactivitybounds);
   if( !SCIPintervalIsEntire(intervalinfty, rhs) )
   {
      SCIP_Real coef;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         coef = consdata->lincoefs[i];
         var  = consdata->linvars[i];

         /* skip fixed variables
          * @todo is that a good or a bad idea?
          *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
          */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            continue;

         if( coef > 0.0 )
         {
            if( SCIPintervalGetSup(rhs) < intervalinfty )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the upper bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  /* tighten upper bound on x to (rhs.sup - (minlinactivity - coef * xlb)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd += coef * SCIPvarGetLbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the minimal linear activity equal -infinity, so
                   * we tighten upper bound on x to just (rhs.sup - minlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetInf(rhs) > -intervalinfty )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the lower bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* tighten lower bound on x to (rhs.inf - (maxlinactivity - coef * xub)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd += coef * SCIPvarGetUbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (rhs.inf - maxlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is +infinity and x is not solely responsible for this */
            }
         }
         else
         {
            assert(coef < 0.0 );
            if( SCIPintervalGetInf(rhs) > -intervalinfty )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the upper bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
                  /* compute upper bound on x to (maxlinactivity - coef * xlb) - rhs.inf / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd += (-coef) * SCIPvarGetLbLocal(var);
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (maxlinactivity - rhs.inf) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetSup(rhs) < intervalinfty )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the lower bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* compute lower bound on x to (minlinactivity - coef * xub) - rhs.sup / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd += (-coef) * SCIPvarGetUbLocal(var);
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal -infinity, so
                   * we tighten lower bound on x to just (minlinactivity - rhs.sup) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }
         }
      }
      if( *result == SCIP_CUTOFF )
         goto CLEANUP;
   }

   /* propagate quadratic part \in rhs = consbounds - linactivity */
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
   consdataUpdateLinearActivity(scip, consdata, intervalinfty); /* make sure, activities of linear part did not become invalid by above bound changes, if any */
   assert(consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity);
   SCIPintervalSetBounds(&tmp,
      (consdata->minlinactivityinf > 0 ? -intervalinfty : consdata->minlinactivity),
      (consdata->maxlinactivityinf > 0 ?  intervalinfty : consdata->maxlinactivity));
   SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);
   if( !SCIPintervalIsEntire(intervalinfty, rhs) )
   {
      if( consdata->nquadvars == 1 )
      {
         /* quadratic part is just a*x^2+b*x -> a common case that we treat directly */
         SCIP_INTERVAL lincoef;    /* linear coefficient of quadratic equation */

         assert(consdata->nbilinterms == 0);

         var = consdata->quadvarterms[0].var;
         SCIPintervalSet(&lincoef, consdata->quadvarterms[0].lincoef);

         /* propagate a*x^2 + b*x \in rhs */
         SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, consdata->quadvarterms[0].sqrcoef, lincoef, rhs, result, nchgbds) );
      }
      else if( consdata->nbilinterms == 1 && consdata->nquadvars == 2 )
      {
         /* quadratic part is just ax*x^2+bx*x + ay*y^2+by*y + c*xy -> a common case that we treat directly */
         assert(consdata->bilinterms[0].var1 == consdata->quadvarterms[0].var || consdata->bilinterms[0].var1 == consdata->quadvarterms[1].var);
         assert(consdata->bilinterms[0].var2 == consdata->quadvarterms[0].var || consdata->bilinterms[0].var2 == consdata->quadvarterms[1].var);

         /* find domain reductions for x from a_x x^2 + b_x x + a_y y^2 + b_y y + c x y \in rhs */
         SCIP_CALL( propagateBoundsBilinearTerm(scip, cons, intervalinfty,
               consdata->quadvarterms[0].var, consdata->quadvarterms[0].sqrcoef, consdata->quadvarterms[0].lincoef,
               consdata->quadvarterms[1].var, consdata->quadvarterms[1].sqrcoef, consdata->quadvarterms[1].lincoef,
               consdata->bilinterms[0].coef,
               rhs, result, nchgbds) );
         if( *result != SCIP_CUTOFF )
         {
            /* find domain reductions for y from a_x x^2 + b_x x + a_y y^2 + b_y y + c x y \in rhs */
            SCIP_CALL( propagateBoundsBilinearTerm(scip, cons, intervalinfty,
                  consdata->quadvarterms[1].var, consdata->quadvarterms[1].sqrcoef, consdata->quadvarterms[1].lincoef,
                  consdata->quadvarterms[0].var, consdata->quadvarterms[0].sqrcoef, consdata->quadvarterms[0].lincoef,
                  consdata->bilinterms[0].coef,
                  rhs, result, nchgbds) );
         }
      }
      else
      {
         /* general case */

         /* compute "advanced" information on quad var term activities, if not up-to-date */
         if( quadminactinf == -1  )
         {
            assert(quadactcontr == NULL);
            SCIP_CALL( SCIPallocBufferArray(scip, &quadactcontr, consdata->nquadvars) );
            propagateBoundsGetQuadActivity(scip, consdata, intervalinfty, &minquadactivity, &maxquadactivity, &quadminactinf, &quadmaxactinf, quadactcontr);
         }
         assert(quadactcontr != NULL);
         assert(quadminactinf >= 0);
         assert(quadmaxactinf >= 0);

         /* if the quad activities are not hopelessly unbounded on useful sides, try to deduce domain reductions on quad vars */
         if( (SCIPintervalGetSup(rhs) <  intervalinfty && quadminactinf <= 1) ||
            ( SCIPintervalGetInf(rhs) > -intervalinfty && quadmaxactinf <= 1) )
         {
            SCIP_INTERVAL lincoef;
            SCIP_INTERVAL rhs2;
            int j;
            int k;

            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;

               /* skip fixed variables
                * @todo is that a good or a bad idea?
                *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
                */
               if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                  continue;

               /* compute rhs2 such that we can propagate quadvarterm(x_i) \in rhs2 */

               /* setup rhs2.sup = rhs.sup - (quadactivity.inf - quadactcontr[i].inf), if everything were finite
                * if only quadactcontr[i].inf is infinite (i.e., the other i are all finite), we just get rhs2.sup = rhs.sup
                * otherwise we get rhs2.sup = infinity */
               if( SCIPintervalGetSup(rhs) < intervalinfty )
               {
                  if( quadminactinf == 0 || (quadminactinf == 1 && SCIPintervalGetInf(quadactcontr[i]) <= -intervalinfty) )
                  {
                     /* the residual quad min activity w.r.t. quad var term i is finite */
                     assert(!SCIPisInfinity(scip, -minquadactivity));  /*lint !e644*/
                     roundmode = SCIPintervalGetRoundingMode();
                     SCIPintervalSetRoundingModeUpwards();
                     rhs2.sup = rhs.sup - minquadactivity;
                     if( quadminactinf == 0 && SCIPintervalGetInf(quadactcontr[i]) != 0.0 )
                     {
                        /* the residual quad min activity w.r.t. quad var term i is finite and nonzero, so add it to right hand side */
                        assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(quadactcontr[i])));
                        rhs2.sup += SCIPintervalGetInf(quadactcontr[i]);
                     }
                     SCIPintervalSetRoundingMode(roundmode);
                  }
                  else
                  {
                     /* there are either >= 2 quad var terms contributing -infinity, or there is one which is not i */
                     rhs2.sup = intervalinfty;
                  }
               }
               else
               {
                  rhs2.sup = intervalinfty;
               }

               /* setup rhs.inf = rhs.inf - (quadactivity.sup - quadactcontr[i].sup), see also above */
               if( SCIPintervalGetInf(rhs) > -intervalinfty )
               {
                  if( quadmaxactinf == 0 || (quadmaxactinf == 1 && SCIPintervalGetSup(quadactcontr[i]) >= intervalinfty) )
                  {
                     /* the residual quad max activity w.r.t. quad var term i is finite and nonzero, so add it to right hand side */
                     assert(!SCIPisInfinity(scip, maxquadactivity));  /*lint !e644*/
                     roundmode = SCIPintervalGetRoundingMode();
                     SCIPintervalSetRoundingModeDownwards();
                     rhs2.inf = rhs.inf - maxquadactivity;
                     /* the residual quad max activity w.r.t. quad var term i is finite */
                     if( quadmaxactinf == 0 && SCIPintervalGetSup(quadactcontr[i]) != 0.0 )
                     {
                        assert(!SCIPisInfinity(scip, SCIPintervalGetSup(quadactcontr[i])));
                        rhs2.inf += SCIPintervalGetSup(quadactcontr[i]);
                     }
                     SCIPintervalSetRoundingMode(roundmode);
                  }
                  else
                  {
                     /* there are either >= 2 quad var terms contributing infinity, or there is one which is not i */
                     rhs2.inf = -intervalinfty;
                  }
               }
               else
               {
                  rhs2.inf = -intervalinfty;
               }
               assert(!SCIPintervalIsEmpty(rhs2));

               /* if rhs2 is entire, then there is nothing we could propagate */
               if( SCIPintervalIsEntire(intervalinfty, rhs2) )
                  continue;

               /* assemble linear coefficient for quad equation a*x^2 + b*x \in rhs2 */
               SCIPintervalSet(&lincoef, consdata->quadvarterms[i].lincoef);
               for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
               {
                  k = consdata->quadvarterms[i].adjbilin[j];
#if 1
                  if( consdata->bilinterms[k].var1 == var )
                  {
                     /* bilinear term k contributes to the activity of quad var term i, so just add bounds to linear coef */
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
                     SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);

                  }
                  else
                  {
                     /* bilinear term k does not contribute to the activity of quad var term i
                      * so bounds on term k are contained in rhs2
                      * if they are finite, we try to remove them from rhs2 and update lincoef instead
                      * if the bounds on bilinear term k as added to rhs2 are old due to recent bound tightening, we may not do best possible, but still correct
                      * HOWEVER: when computing rhs2, we may not just have added the bounds for the bilinear term, but for the associated quadratic term
                      *   for this complete term, we used SCIPintervalQuad to compute the bounds
                      *   since we do not want to repeat a call to SCIPintervalQuad for that quadratic term with bilinear term k removed,
                      *   we only remove the bounds for the bilinear term k from rhs2 if the associated quadratic term consists only of this bilinear term,
                      *   i.e., the quadratic term corresponding to var1 should be only var1*var2, but have no square or linear coefs or other bilinear terms
                      *   (for efficiency reasons, we check here only if there are any other bilinear terms than var1*var2 associated with var1, even if they are not associated with the quad var term for var1)
                      */
                     SCIP_INTERVAL me;
                     SCIP_INTERVAL bilinbounds;
                     int otherpos;

                     assert(consdata->bilinterms[k].var2 == var);

                     assert(consdata->quadvarssorted);
                     SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[k].var1, &otherpos) );
                     assert(otherpos >= 0);
                     assert(consdata->quadvarterms[otherpos].var == consdata->bilinterms[k].var1);

                     if( (consdata->quadvarterms[otherpos].sqrcoef != 0.0) || consdata->quadvarterms[otherpos].lincoef != 0.0 || consdata->quadvarterms[otherpos].nadjbilin > 1 )
                        continue;

                     /* set tmp to bounds of other variable and multiply with bilin coef */
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var1), SCIPvarGetUbLocal(consdata->bilinterms[k].var1))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var1), SCIPvarGetUbLocal(consdata->bilinterms[k].var1))));
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);

                     /* set me to bounds of i'th variable */
                     SCIPintervalSetBounds(&me,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));

                     /* remove me*tmp from rhs2 */

                     roundmode = SCIPintervalGetRoundingMode();

                     if( rhs2.inf > -intervalinfty )
                     {
                        /* need upward rounding for SCIPintervalMulSup */
                        SCIPintervalSetRoundingModeUpwards();
                        SCIPintervalMulSup(intervalinfty, &bilinbounds, me, tmp);
                        /* rhs2.inf += bilinbounds.sup, but we are in upward rounding */
                        if( bilinbounds.sup < intervalinfty )
                           rhs2.inf = SCIPintervalNegateReal(SCIPintervalNegateReal(rhs2.inf) - bilinbounds.sup);
                     }

                     if( rhs2.sup <  intervalinfty )
                     {
                        /* need downward rounding for SCIPintervalMulInf */
                        SCIPintervalSetRoundingModeDownwards();
                        SCIPintervalMulInf(intervalinfty, &bilinbounds, me, tmp);
                        /* rhs2.sup += bilinbounds.inf, but we are in downward rounding */
                        if( bilinbounds.inf > -intervalinfty )
                           rhs2.sup = SCIPintervalNegateReal(SCIPintervalNegateReal(rhs2.sup) - bilinbounds.inf);
                     }

                     SCIPintervalSetRoundingMode(roundmode);

                     /* add tmp to lincoef */
                     SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
                  }
#else
                  if( consdata->bilinterms[k].var1 != var )
                     continue; /* this term does not contribute to the activity of quad var term i */

                  SCIPintervalSetBounds(&tmp,
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                     +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
                  SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
                  SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
#endif
               }

               /* deduce domain reductions for x_i */
               SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, consdata->quadvarterms[i].sqrcoef, lincoef, rhs2, result, nchgbds) );
               if( *result == SCIP_CUTOFF )
                  goto CLEANUP;
            }
         }
      }
   }

 CLEANUP:
   SCIPfreeBufferArrayNull(scip, &quadactcontr);

   return SCIP_OKAY;
}

/** calls domain propagation for a set of constraints */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds             /**< buffer where to add the the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT propresult;
   SCIP_Bool   redundant;
   int         c;
   int         roundnr;
   SCIP_Bool   success;
   int         maxproprounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);
   assert(nchgbds != NULL);

   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      maxproprounds = conshdlrdata->maxproproundspresolve;
   else
      maxproprounds = conshdlrdata->maxproprounds;

   do
   {
      success = FALSE;
      ++roundnr;

      SCIPdebugMessage("starting domain propagation round %d of %d for %d constraints\n", roundnr, maxproprounds, nconss);

      for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
      {
         assert(conss != NULL);
         if( !SCIPconsIsEnabled(conss[c]) )
            continue;

         SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );
         if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         {
            *result = propresult;
            success = TRUE;
         }
         if( redundant )
         {
            SCIPdebugMessage("deleting constraint <%s> locally\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
         }
      }

   }
   while( success && *result != SCIP_CUTOFF && roundnr < maxproprounds );

   return SCIP_OKAY;
}

/* checks for a linear variable that can be increase or decreased without harming feasibility */
static
void consdataFindUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   int poslock;
   int neglock;

   consdata->linvar_maydecrease = -1;
   consdata->linvar_mayincrease = -1;

   /* check for a linear variable that can be increase or decreased without harming feasibility */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      /* compute locks of i'th linear variable */
      assert(consdata->lincoefs[i] != 0.0);
      if( consdata->lincoefs[i] > 0.0 )
      {
         poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
      }
      else
      {
         poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - neglock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can decrease x without harming other constraints */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_maydecrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_maydecrease]) / consdata->lincoefs[consdata->linvar_maydecrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_maydecrease = i;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - poslock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can increase x without harm */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_mayincrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_mayincrease]) / consdata->lincoefs[consdata->linvar_mayincrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_mayincrease = i;
      }
   }

#ifdef SCIP_DEBUG
   if( consdata->linvar_mayincrease >= 0 )
   {
      SCIPdebugMessage("may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_mayincrease]));
   }
   if( consdata->linvar_maydecrease >= 0 )
   {
      SCIPdebugMessage("may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_maydecrease]));
   }
#endif
}

/** Given a solution where every quadratic constraint is either feasible or can be made feasible by
 * moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 * The method assumes that this is always possible and that not all constraints are feasible already.
 */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to process */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded to construct a solution that satisfies all provided constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_SOL* newsol;
   SCIP_VAR* var;
   int c;
   SCIP_Real viol;
   SCIP_Real norm;
   SCIP_Real delta;
   SCIP_Real gap;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   *success = FALSE;

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );
   SCIPdebugMessage("attempt to make solution from <%s> feasible by shifting linear variable\n",
      sol != NULL ? (SCIPsolGetHeur(sol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(sol)) : "tree") : "LP");

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* recompute violation of solution in case solution has changed
       * get absolution violation and sign */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conss[c], newsol, TRUE) );  /*lint !e613*/
         viol = consdata->lhs - consdata->activity;
      }
      else if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conss[c], newsol, TRUE) );  /*lint !e613*/
         viol = consdata->rhs - consdata->activity;
      }
      else
         continue; /* constraint is satisfied */

      assert(viol != 0.0);
      if( consdata->linvar_mayincrease >= 0 &&
         ((viol > 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) || (viol < 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0)) )
      {
         /* have variable where increasing makes the constraint less violated */
         var = consdata->linvars[consdata->linvar_mayincrease];
         /* compute how much we would like to increase var */
         delta = viol / consdata->lincoefs[consdata->linvar_mayincrease];
         assert(delta > 0.0);
         /* if var has an upper bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            gap = SCIPvarGetUbGlobal(var) - SCIPgetSolVal(scip, newsol, var);
            delta = MIN(MAX(0.0, gap), delta);
         }
         if( SCIPisPositive(scip, delta) )
         {
            /* if variable is integral, round delta up so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPceil(scip, delta);

            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMessage("increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_mayincrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvar_maydecrease >= 0 &&
         ((viol > 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) || (viol < 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0)) )
      {
         /* have variable where decreasing makes constraint less violated */
         var = consdata->linvars[consdata->linvar_maydecrease];
         /* compute how much we would like to decrease var */
         delta = viol / consdata->lincoefs[consdata->linvar_maydecrease];
         assert(delta < 0.0);
         /* if var has a lower bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            gap = SCIPgetSolVal(scip, newsol, var) - SCIPvarGetLbGlobal(var);
            delta = MAX(MIN(0.0, gap), delta);
         }
         if( SCIPisNegative(scip, delta) )
         {
            /* if variable is integral, round delta down so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPfloor(scip, delta);
            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMessage("increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_maydecrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so maybe we could not make constraint feasible due to variable bounds
       * check if we are feasible w.r.t. (relative) feasibility tolerance */
      norm = getGradientMaxElement(scip, conss[c], newsol);  /*lint !e613*/
      if( norm > 1.0 )
         viol /= norm;
      /* if still violated, we give up */
      if( SCIPisGT(scip, REALABS(viol), SCIPfeastol(scip)) )
         break;

      /* if objective value is not better than current upper bound, we give up */
      if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && !SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip)) )
         break;
   }

   /* if we have a solution that should satisfy all quadratic constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic */
   if( c == nconss )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      SCIPdebugMessage("pass solution with objective val %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->trysolheur != NULL);

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );
      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** tries to upgrade a nonlinear constraint into a quadratic constraint */
static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdQuadratic)
{
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   int i;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   /* if a quadratic expression has been simplified, then all children of the node should be variables */
   if( !SCIPexprgraphAreAllNodeChildrenVars(node) )
      return SCIP_OKAY;

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_CONST:
   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   case SCIP_EXPR_SUM:
   case SCIP_EXPR_LINEAR:
      /* these should not appear as exprgraphnodes after constraint presolving */
      return SCIP_OKAY;

   case SCIP_EXPR_DIV:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_INTPOWER:
   case SCIP_EXPR_SIGNPOWER:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_PRODUCT:
   case SCIP_EXPR_POLYNOMIAL:
      /* these do not look like an quadratic expression (assuming the expression graph simplifier did run) */
      return SCIP_OKAY;

   case SCIP_EXPR_MUL:
   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_QUADRATIC:
      /* these mean that we have something quadratic */
      break;

   case SCIP_EXPR_PARAM:
   case SCIP_EXPR_LAST:
   default:
      SCIPwarningMessage(scip, "unexpected expression operator %d in nonlinear constraint <%s>\n", SCIPexprgraphGetNodeOperator(node), SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   /* setup a quadratic constraint */

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsQuadratic(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
         0, NULL, 0, NULL,
         SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
   assert(!SCIPconsIsStickingAtNode(cons));

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));

   /* add variables from expression tree as "quadratic" variables to quadratic constraint */
   for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
   {
      assert(SCIPexprgraphGetNodeChildren(node)[i] != NULL);
      SCIP_CALL( SCIPaddQuadVarQuadratic(scip, upgdconss[0], (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[i]), 0.0, 0.0) );
   }

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_MUL:
   {
      /* expression is product of two variables, so add bilinear term to constraint */
      assert(SCIPexprgraphGetNodeNChildren(node) == 2);

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, upgdconss[0],
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[0]),
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[1]),
            1.0) );

      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      /* expression is square of a variable, so change square coefficient of quadratic variable */
      assert(SCIPexprgraphGetNodeNChildren(node) == 1);

      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, upgdconss[0],
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[0]),
            1.0) );

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      /* expression is quadratic */
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_Real* lincoefs;

      quadelems  = SCIPexprgraphGetNodeQuadraticQuadElements(node);
      nquadelems = SCIPexprgraphGetNodeQuadraticNQuadElements(node);
      lincoefs   = SCIPexprgraphGetNodeQuadraticLinearCoefs(node);

      SCIPaddConstantQuadratic(scip, upgdconss[0], SCIPexprgraphGetNodeQuadraticConstant(node));

      if( lincoefs != NULL )
         for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
            if( lincoefs[i] != 0.0 )
            {
               /* linear term */
               SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, upgdconss[0],
                     (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[i]),
                     lincoefs[i]) );
            }

      for( i = 0; i < nquadelems; ++i )
      {
         assert(quadelems[i].idx1 < SCIPexprgraphGetNodeNChildren(node));
         assert(quadelems[i].idx2 < SCIPexprgraphGetNodeNChildren(node));

         if( quadelems[i].idx1 == quadelems[i].idx2 )
         {
            /* square term */
            SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, upgdconss[0],
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx1]),
                  quadelems[i].coef) );
         }
         else
         {
            /* bilinear term */
            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, upgdconss[0],
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx1]),
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx2]),
                  quadelems[i].coef) );
         }
      }

      break;
   }

   default:
   {
      SCIPerrorMessage("you should not be here\n");
      return SCIP_ERROR;
   }
   }  /*lint !e788 */

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyQuadratic)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      assert(conshdlrdata->quadconsupgrades[i] != NULL);
      SCIPfreeMemory(scip, &conshdlrdata->quadconsupgrades[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->quadconsupgrades);

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* catch variable events */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* drop variable events */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
   }

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIP_OKAY;
}
#endif

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA*     consdata;
   int                c;
#ifndef NDEBUG
   int                i;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->isremovedfixings )
      {
         SCIP_CALL( removeFixedVariables(scip, conss[c]) );
      }

      /* make sure we do not have duplicate bilinear terms, quad var terms, or linear vars */
      SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
      SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
      SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );

      assert(consdata->isremovedfixings);
      assert(consdata->linvarsmerged);
      assert(consdata->quadvarsmerged);
      assert(consdata->bilinmerged);

#ifndef NDEBUG
      for( i = 0; i < consdata->nlinvars; ++i )
         assert(SCIPvarIsActive(consdata->linvars[i]));

      for( i = 0; i < consdata->nquadvars; ++i )
         assert(SCIPvarIsActive(consdata->quadvarterms[i].var));
#endif

      /* tell SCIP that we have something nonlinear */
      if( SCIPconsIsAdded(conss[c]) && consdata->nquadvars > 0 )
         SCIPenableNLP(scip);
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for a linear variable that can be increase or decreased without harming feasibility */
      consdataFindUnlockedLinearVar(scip, consdata);

      /* setup lincoefsmin, lincoefsmax */
      consdata->lincoefsmin = SCIPinfinity(scip);
      consdata->lincoefsmax = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         consdata->lincoefsmin = MIN(consdata->lincoefsmin, REALABS(consdata->lincoefs[i]));  /*lint !e666 */
         consdata->lincoefsmax = MAX(consdata->lincoefsmax, REALABS(consdata->lincoefs[i]));  /*lint !e666 */
      }

      /* add nlrow representation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )
      {
         if( consdata->nlrow == NULL )
         {
            SCIP_CALL( createNlRow(scip, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }

      /* setup sepaquadvars and sepabilinvar2pos */
      assert(consdata->sepaquadvars == NULL);
      assert(consdata->sepabilinvar2pos == NULL);
      if( consdata->nquadvars > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->sepaquadvars,     consdata->nquadvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->sepabilinvar2pos, consdata->nbilinterms) );

         /* make sure, quadratic variable terms are sorted */
         SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

         for( i = 0; i < consdata->nquadvars; ++i )
            consdata->sepaquadvars[i] = consdata->quadvarterms[i].var;

         for( i = 0; i < consdata->nbilinterms; ++i )
         {
            SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[i].var2, &consdata->sepabilinvar2pos[i]) );
         }
      }

      if( conshdlrdata->checkfactorable )
      {
         /* check if constraint function is factorable, i.e., can be written as product of two linear functions */
         SCIP_CALL( checkFactorable(scip, conss[c]) );
      }
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 && conshdlrdata->linearizeheursol )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      /* @todo Should we catch every new solution or only new *best* solutions */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   if( nconss != 0 && !SCIPisIpoptAvailableIpopt() && !SCIPisInRestart(scip) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Quadratic constraint handler does not have LAPACK for eigenvalue computation. Will assume that matrices (with size > 2x2) are indefinite.\n");
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfolpnode = NULL;
   conshdlrdata->nenfolprounds = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      assert(consdata->sepaquadvars     != NULL || consdata->nquadvars == 0);
      assert(consdata->sepabilinvar2pos != NULL || consdata->nquadvars == 0);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->sepaquadvars,     consdata->nquadvars);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->sepabilinvar2pos, consdata->nbilinterms);

      SCIPfreeBlockMemoryArrayNull(scip, &consdata->factorleft,  consdata->nquadvars + 1);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->factorright, consdata->nquadvars + 1);
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(SCIPconsGetData(cons) == *consdata);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   SCIP_CALL( consdataFree(scip, consdata) );

   assert(*consdata == NULL);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransQuadratic)
{  
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int            i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->lhs, sourcedata->rhs,
         sourcedata->nlinvars, sourcedata->linvars, sourcedata->lincoefs,
         sourcedata->nquadvars, sourcedata->quadvarterms,
         sourcedata->nbilinterms, sourcedata->bilinterms,
         FALSE) );

   for( i = 0; i < targetdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->linvars[i], &targetdata->linvars[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvars[i]) );
   }

   for( i = 0; i < targetdata->nquadvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->quadvarterms[i].var, &targetdata->quadvarterms[i].var) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->quadvarterms[i].var) );
   }

   for( i = 0; i < targetdata->nbilinterms; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->bilinterms[i].var1, &targetdata->bilinterms[i].var1) );
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->bilinterms[i].var2, &targetdata->bilinterms[i].var2) );

      if( SCIPvarCompare(targetdata->bilinterms[i].var1, targetdata->bilinterms[i].var2) > 0 )
      {
         SCIP_VAR* tmp;
         tmp = targetdata->bilinterms[i].var2;
         targetdata->bilinterms[i].var2 = targetdata->bilinterms[i].var1;
         targetdata->bilinterms[i].var1 = tmp;
      }
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   SCIPdebugMessage("created transformed quadratic constraint ");
   SCIPdebugPrintCons(scip, *targetcons, NULL);

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   SCIP_ROW*          row;
   SCIP_Real*         x;
   int                c;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      row = NULL;

      if( consdata->nquadvars == 0 )
      {
         /* if we are actually linear, add the constraint as row to the LP */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(conss[c]), SCIPconsGetName(conss[c]), consdata->lhs, consdata->rhs,
               SCIPconsIsLocal(conss[c]), FALSE , TRUE) );  /*lint !e613 */
         SCIP_CALL( SCIPaddVarsToRow(scip, row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
         SCIP_CALL( SCIPreleaseRow (scip, &row) );
         continue;
      }

      /* alloc memory for reference point */
      SCIP_CALL( SCIPallocBufferArray(scip, &x, consdata->nquadvars) );

      /* for convex parts, add linearizations in 5 points */
      if( (consdata->isconvex && !SCIPisInfinity(scip,  consdata->rhs)) ||
         (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)) )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real lambda;
         int k;

         for( k = 0; k < 5; ++k )
         {
            lambda = 0.1 * (k+1); /* lambda = 0.1, 0.2, 0.3, 0.4, 0.5 */
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;
               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);

               if( ub > -INITLPMAXVARVAL )
                  lb = MAX(lb, -INITLPMAXVARVAL);
               if( lb <  INITLPMAXVARVAL )
                  ub = MIN(ub,  INITLPMAXVARVAL);

               /* make bounds finite */
               if( SCIPisInfinity(scip, -lb) )
                  lb = MIN(-10.0, ub - 0.1*REALABS(ub));  /*lint !e666 */
               if( SCIPisInfinity(scip,  ub) )
                  ub = MAX( 10.0, lb + 0.1*REALABS(lb));  /*lint !e666 */

               if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
                  x[i] = lambda * ub + (1.0 - lambda) * lb;
               else
                  x[i] = lambda * lb + (1.0 - lambda) * ub;
            }

            SCIP_CALL( generateCut(scip, conss[c], x, NULL, consdata->isconvex ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT, &row, NULL, conshdlrdata->cutmaxrange, FALSE, -SCIPinfinity(scip)) );  /*lint !e613 */
            if( row != NULL )
            {
               SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
               SCIPdebugMessage("initlp adds row <%s> for lambda = %g of conss <%s>\n", SCIProwGetName(row), lambda, SCIPconsGetName(conss[c]));  /*lint !e613 */
               SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
               SCIP_CALL( SCIPreleaseRow (scip, &row) );
            }
         }
      }

      /* for concave parts, add underestimator w.r.t. at most 2 reference points */
      if( (!consdata->isconvex && !SCIPisInfinity(scip,  consdata->rhs)) ||
         (!consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)) )
      {
         SCIP_Bool unbounded;
         SCIP_Bool possquare;
         SCIP_Bool negsquare;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real lambda;
         int k;

         unbounded = FALSE; /* whether there are unbounded variables */
         possquare = FALSE; /* whether there is a positive square term */
         negsquare = FALSE; /* whether there is a negative square term */
         lambda = 0.6; /* weight of preferred bound */
         for( k = 0; k < 2; ++k )
         {
            /* set reference point to 0 projected on bounds for unbounded variables or in between lower and upper bound for bounded variables
             * in the first round, we set it closer to the best bound, in the second closer to the worst bound
             * the reason is, that for a bilinear term with bounded variables, there are always two linear underestimators
             * if the reference point is set to the middle, then rounding and luck decides which underestimator is chosen
             * we thus choose the reference point not to be the middle, so both McCormick terms are definitely chosen one time
             * of course, the possible number of cuts is something in the order of 2^nquadvars, and we choose two of them here
             */
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;
               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);

               if( SCIPisInfinity(scip, -lb) )
               {
                  if( SCIPisInfinity(scip, ub) )
                     x[i] = 0.0;
                  else
                     x[i] = MIN(0.0, ub);
                  unbounded = TRUE;
               }
               else
               {
                  if( SCIPisInfinity(scip, ub) )
                  {
                     x[i] = MAX(0.0, lb);
                     unbounded = TRUE;
                  }
                  else
                     x[i] = lambda * SCIPvarGetBestBoundLocal(var) + (1.0-lambda) * SCIPvarGetWorstBoundLocal(var);
               }

               possquare |= consdata->quadvarterms[i].sqrcoef > 0.0;  /*lint !e514 */
               negsquare |= consdata->quadvarterms[i].sqrcoef < 0.0;  /*lint !e514 */
            }

            if( !consdata->isconvex  && !SCIPisInfinity(scip,  consdata->rhs) )
            {
               SCIP_CALL( generateCut(scip, conss[c], x, NULL, SCIP_SIDETYPE_RIGHT, &row, NULL, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
               if( row != NULL )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
                  SCIPdebugMessage("initlp adds row <%s> for rhs of conss <%s>, round %d\n", SCIProwGetName(row), SCIPconsGetName(conss[c]), k);  /*lint !e613 */
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
                  SCIP_CALL( SCIPreleaseRow (scip, &row) );
               }
            }
            if( !consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( generateCut(scip, conss[c], x, NULL, SCIP_SIDETYPE_LEFT, &row, NULL, conshdlrdata->cutmaxrange, conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
               if( row != NULL )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
                  SCIPdebugMessage("initlp adds row <%s> for lhs of conss <%s>, round %d\n", SCIProwGetName(row), SCIPconsGetName(conss[c]), k);  /*lint !e613 */
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
                  SCIP_CALL( SCIPreleaseRow (scip, &row) );
               }
            }

            /* if there are unbounded variables, then there is typically only at most one possible underestimator, so don't try another round
             * similar, if there are no bilinear terms and no linearizations of square terms, then the reference point does not matter, so don't do another round */
            if( unbounded ||
               (consdata->nbilinterms == 0 && (!possquare || SCIPisInfinity(scip,  consdata->rhs))) ||
               (consdata->nbilinterms == 0 && (!negsquare || SCIPisInfinity(scip, -consdata->lhs))) )
               break;

            /* invert lambda for second round */
            lambda = 1.0 - lambda;
         }
      }

      SCIPfreeBufferArray(scip, &x);
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || (SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY && conshdlrdata->sepanlpmincont <= 1.0)) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CONSDATA*  consdata;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool       solvednlp;
      int c;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we might want to do this
          * but first check whether there is a violated constraint side which corresponds to a convex function
          */
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);  /*lint !e613 */

            consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
            assert(consdata != NULL);

            /* skip feasible constraints */
            if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
               continue;

            /* make sure curvature has been checked */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

            if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && consdata->isconvex) ||
               ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && consdata->isconcave) )
               break;
         }

         if( c < nconss )
         {
            /* try to solve NLP and update solstat */

            /* ensure linear conss are in NLP */
            if( conshdlrdata->subnlpheur != NULL )
            {
               SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, conshdlrdata->subnlpheur, TRUE, TRUE) );
            }

            /* set LP solution as starting values, if available */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );
            }

            /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
            SCIP_CALL( SCIPsolveNLP(scip) );

            solstat = SCIPgetNLPSolstat(scip);
            SCIPdebugMessage("solved NLP relax, solution status: %d\n", solstat);

            solvednlp = TRUE;
         }
      }

      conshdlrdata->sepanlp = TRUE;

      if( solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      {
         SCIPdebugMessage("NLP relaxation is globally infeasible, thus can cutoff node\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* if we have feasible NLP solution, generate linearization cuts there */
         SCIP_Bool lpsolseparated;
         SCIP_SOL* nlpsol;

         SCIP_CALL( SCIPcreateNLPSol(scip, &nlpsol, NULL) );
         assert(nlpsol != NULL);

         /* if we solved the NLP and solution is integral, then pass it to trysol heuristic */
         if( solvednlp && conshdlrdata->trysolheur != NULL )
         {
            int nfracvars;

            nfracvars = 0;
            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIPdebugMessage("pass solution with obj. value %g to trysol\n", SCIPgetSolOrigObj(scip, nlpsol));
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, conshdlrdata->mincutefficacysepa) );

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         /* if a cut that separated the LP solution was added, then return, otherwise continue with usual separation in LP solution */
         if( lpsolseparated )
         {
            SCIPdebugMessage("linearization cuts separate LP solution\n");
            *result = SCIP_SEPARATED;

            return SCIP_OKAY;
         }
      }
   }
   /* if we do not want to try solving the NLP, or have no NLP, or have no NLP solver, or solving the NLP failed,
    * or separating with NLP solution as reference point failed, then try (again) with LP solution as reference point
    */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, conshdlrdata->mincutefficacysepa, FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, conshdlrdata->mincutefficacysepa, FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Real          maxviol;
   SCIP_RESULT        propresult;
   SCIP_RESULT        separateresult;
   int                nchgbds;
   int                nnotify;
   SCIP_Real          sepaefficacy;
   SCIP_Real          minefficacy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   consdata = SCIPconsGetData(maxviolcon);
   assert(consdata != NULL);
   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   SCIPdebugMessage("enfolp with max violation %g in cons <%s>\n", maxviol, SCIPconsGetName(maxviolcon));

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
    * in this case, check if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an assert in scip.c)
    * the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may be expensive
    * we only increment nenfolprounds until 101 to avoid an overflow
    */
   if( conshdlrdata->lastenfolpnode == SCIPgetCurrentNode(scip) )
   {
      if( conshdlrdata->nenfolprounds > 100 )
      {
         if( SCIPisStopped(scip) )
         {
            SCIP_NODE* child;

            SCIP_CALL( SCIPcreateChild(scip, &child, 1.0, SCIPnodeGetEstimate(SCIPgetCurrentNode(scip))) );
            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }
      else
         ++conshdlrdata->nenfolprounds;
   }
   else
   {
      conshdlrdata->lastenfolpnode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenfolprounds = 0;
   }

   /* run domain propagation */
   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &nchgbds) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      SCIPdebugMessage("propagation succeeded (%s)\n", propresult == SCIP_CUTOFF ? "cutoff" : "reduceddom");
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>feastol)
    * however, if the maximal violation is very small, also the best cut efficacy cannot be large
    * thus, in the latter case, we are also happy if the efficacy is at least, say, 75% of the maximal violation
    * but in any case we need an efficacy that is at least feastol
    */
   minefficacy = MIN(0.75*maxviol, conshdlrdata->mincutefficacyenfofac * SCIPfeastol(scip));  /*lint !e666 */
   minefficacy = MAX(minefficacy, SCIPfeastol(scip));  /*lint !e666 */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, minefficacy, TRUE, &separateresult, &sepaefficacy) );
   if( separateresult == SCIP_SEPARATED )
   {
      SCIPdebugMessage("separation succeeded (bestefficacy = %g, minefficacy = %g)\n", sepaefficacy, minefficacy);
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a good cut
    * -> collect variables for branching
    */

   SCIPdebugMessage("separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy, minefficacy, maxviol);

   /* find branching candidates */
   SCIP_CALL( registerVariableInfeasibilities(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify == 0 && !solinfeasible && minefficacy > SCIPfeastol(scip) )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, SCIPfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
      if( separateresult == SCIP_SEPARATED )
      {
         SCIPdebugMessage("separation fallback succeeded, efficacy = %g\n", sepaefficacy);
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }

   if( nnotify == 0 && !solinfeasible )
   {
      /* fallback 2: separation probably failed because of numerical difficulties with a convex constraint;
       *  if noone declared solution infeasible yet and we had not even found a weak cut, try to resolve by branching
       */
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeLPValueVariableForBranching(scip, conss, nconss, &brvar) );
      if( brvar == NULL )
      {
         /* fallback 3: all quadratic variables seem to be fixed -> replace by linear constraint */
         SCIP_Bool addedcons;
         SCIP_Bool reduceddom;
         SCIP_Bool infeasible;

         SCIP_CALL( replaceByLinearConstraints(scip, conss, nconss, &addedcons, &reduceddom, &infeasible) );
         /* if the linear constraints are actually feasible, then adding them and returning SCIP_CONSADDED confuses SCIP
          * when it enforces the new constraints again and nothing resolves the infeasiblity that we declare here thus,
          * we only add them if considered violated, and otherwise claim the solution is feasible (but print a
          * warning) */
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else if ( addedcons )
            *result = SCIP_CONSADDED;
         else if ( reduceddom )
            *result = SCIP_REDUCEDDOM;
         else
         {
            *result = SCIP_FEASIBLE;
            SCIPwarningMessage(scip, "could not enforce feasibility by separating or branching; declaring solution with viol %g as feasible\n", maxviol);
            assert(!SCIPisInfinity(scip, maxviol));
         }
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("Could not find any usual branching variable candidate. Proposed variable <%s> with LP value %g for branching.\n",
            SCIPvarGetName(brvar), SCIPgetSolVal(scip, NULL, brvar));
         nnotify = 1;
      }
   }

   assert(*result == SCIP_INFEASIBLE && (solinfeasible || nnotify > 0));
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_VAR*          var;
   int                c;
   int                i;
   int                nchgbds;
   int                nnotify;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   SCIPdebugMessage("enfops with max violation in cons <%s>\n", SCIPconsGetName(maxviolcon));

   /* run domain propagation */
   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &nchgbds) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> collect all variables in violated constraints for branching
    */
   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         var = consdata->linvars[i];
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         var = consdata->quadvarterms[i].var;
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }
   }

   if( nnotify == 0 )
   {
      SCIPdebugMessage("All variables in violated constraints fixed (up to epsilon). Cannot find branching candidate. Forcing solution of LP.\n");
      *result = SCIP_SOLVELP;
   }

   assert(*result == SCIP_SOLVELP || (*result == SCIP_INFEASIBLE && nnotify > 0));
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropQuadratic)
{
   int         nchgbds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, result, &nchgbds) );

   return SCIP_OKAY;
}  /*lint !e715 */

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolQuadratic)
{  /*lint --e{715,788}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        solveresult;
   SCIP_Bool          redundant;
   SCIP_Bool          havechange;
   SCIP_Bool          doreformulations;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   /* if other presolvers did not find enough changes for another presolving round,
    * then try the reformulations (replacing products with binaries, disaggregation, setting default variable bounds)
    * otherwise, we wait with these
    * @todo first do all usual presolving steps, then check SCIPisPresolveFinished(scip), and if true then do reformulations (and usual steps again)
    */
   doreformulations = (nrounds > 0 || SCIPconshdlrWasPresolvingDelayed(conshdlr)) && SCIPisPresolveFinished(scip);
   SCIPdebugMessage("presolving will %swait with reformulation\n", doreformulations ? "not " : "");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIPdebugMessage("process constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIPdebugPrintCons(scip, conss[c], NULL);

      if( !consdata->initialmerge )
      {
         SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
         SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
         SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );
         consdata->initialmerge = TRUE;
      }

      havechange = FALSE;
#ifdef CHECKIMPLINBILINEAR
      if( consdata->isimpladded )
      {
         int nbilinremoved;
         SCIP_CALL( presolveApplyImplications(scip, conss[c], &nbilinremoved) );
         if( nbilinremoved > 0 )
         {
            *nchgcoefs += nbilinremoved;
            havechange = TRUE;
            *result = SCIP_SUCCESS;
         }
         assert(!consdata->isimpladded);
      }
#endif
      /* call upgrade methods if the constraint has not been presolved yet or there has been a bound tightening or possibly be a change in variable type
       * we want to do this before (multi)aggregated variables are replaced, since that may change structure, e.g., introduce bilinear terms
       */
      if( !consdata->ispresolved || !consdata->ispropagated || nnewchgvartypes > 0 )
      {
         SCIP_Bool upgraded;

         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );
         if( upgraded )
         {
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      if( !consdata->isremovedfixings )
      {
         SCIP_CALL( removeFixedVariables(scip, conss[c]) );
         assert(consdata->isremovedfixings);
         havechange = TRUE;
      }

      /* try to "solve" the constraint, e.g., reduce to a variable aggregation */
      SCIP_CALL( presolveSolve(scip, conss[c], &solveresult, &redundant, naggrvars) );
      if( solveresult == SCIP_CUTOFF )
      {
         SCIPdebugMessage("solving constraint <%s> says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( redundant )
      {
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++*ndelconss;
         *result = SCIP_SUCCESS;
         break;
      }
      if( solveresult == SCIP_SUCCESS )
      {
         *result = SCIP_SUCCESS;
         havechange = TRUE;
      }

      /* @todo divide constraint by gcd of coefficients if all are integral */

      if( doreformulations )
      {
         int naddconss_old;

         naddconss_old = *naddconss;

         SCIP_CALL( presolveTryAddAND(scip, conshdlr, conss[c], naddconss) );
         assert(*naddconss >= naddconss_old);

         if( *naddconss == naddconss_old )
         {
            /* user not so empathic about AND, or we don't have products of two binaries, so try this more general reformulation */
            SCIP_CALL( presolveTryAddLinearReform(scip, conshdlr, conss[c], naddconss) );
            assert(*naddconss >= naddconss_old);
         }

         if( conshdlrdata->disaggregate )
         {
            /* try disaggregation, if enabled */
            SCIP_CALL( presolveDisaggregate(scip, conshdlr, conss[c], naddconss) );
         }

         if( *naddconss > naddconss_old )
         {
            /* if something happened, report success and cleanup constraint */
            *result = SCIP_SUCCESS;
            havechange = TRUE;
            SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );
         }
      }

      if( consdata->nlinvars == 0 && consdata->nquadvars == 0 )
      {
         /* all variables fixed or removed, constraint function is 0.0 now */
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) ); /* well, there shouldn't be any variables left anyway */
         if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasPositive(scip, consdata->lhs)) ||
            ( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasNegative(scip, consdata->rhs)) )
         { /* left hand side positive or right hand side negative */
            SCIPdebugMessage("constraint <%s> is constant and infeasible\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            ++*ndelconss;
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* left and right hand side are consistent */
         SCIPdebugMessage("constraint <%s> is constant and feasible, deleting\n", SCIPconsGetName(conss[c]));
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++*ndelconss;
         *result = SCIP_SUCCESS;
         continue;
      }

      if( !consdata->ispropagated )
      {
         /* try domain propagation if there were bound changes or constraint has changed (in which case, processVarEvents may have set ispropagated to false) */
         SCIP_RESULT propresult;
         int roundnr;

         roundnr = 0;
         do
         {
            ++roundnr;

            SCIPdebugMessage("starting domain propagation round %d of %d\n", roundnr, conshdlrdata->maxproproundspresolve);

            SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );

            if( propresult == SCIP_CUTOFF )
            {
               SCIPdebugMessage("propagation on constraint <%s> says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            /* delete constraint if found redundant by bound tightening */
            if( redundant )
            {
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
               SCIP_CALL( SCIPdelCons(scip, conss[c]) );
               ++*ndelconss;
               *result = SCIP_SUCCESS;
               break;
            }

            if( propresult == SCIP_REDUCEDDOM )
            {
               *result = SCIP_SUCCESS;
               havechange = TRUE;
            }

         }
         while( !consdata->ispropagated && roundnr < conshdlrdata->maxproproundspresolve );

         if( redundant )
            continue;
      }

      /* check if we have a single linear continuous variable that we can make implicit integer */
      if( (nnewchgvartypes != 0 || havechange || !consdata->ispresolved)
         && (SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->lhs)) )
      {
         int       ncontvar;
         SCIP_VAR* candidate;
         SCIP_Bool fail;

         fail = FALSE;
         candidate = NULL;
         ncontvar = 0;

         for( i = 0; !fail && i < consdata->nlinvars; ++i )
         {
            if( !SCIPisIntegral(scip, consdata->lincoefs[i]) )
            {
               fail = TRUE;
            }
            else if( SCIPvarGetType(consdata->linvars[i]) == SCIP_VARTYPE_CONTINUOUS )
            {
               if( ncontvar > 0 ) /* now at 2nd continuous variable */
                  fail = TRUE;
               else if( SCIPisEQ(scip, ABS(consdata->lincoefs[i]), 1.0) )
                  candidate = consdata->linvars[i];
               ++ncontvar;
            }
         }
         for( i = 0; !fail && i < consdata->nquadvars; ++i )
            fail = SCIPvarGetType(consdata->quadvarterms[i].var) == SCIP_VARTYPE_CONTINUOUS ||
               !SCIPisIntegral(scip, consdata->quadvarterms[i].lincoef) ||
               !SCIPisIntegral(scip, consdata->quadvarterms[i].sqrcoef);
         for( i = 0; !fail && i < consdata->nbilinterms; ++i )
            fail = !SCIPisIntegral(scip, consdata->bilinterms[i].coef);

         if( !fail && candidate != NULL )
         {
            SCIP_Bool infeasible;

            SCIPdebugMessage("make variable <%s> implicit integer due to constraint <%s>\n", SCIPvarGetName(candidate), SCIPconsGetName(conss[c]));

            SCIP_CALL( SCIPchgVarType(scip, candidate, SCIP_VARTYPE_IMPLINT, &infeasible) );
            if( infeasible )
            {
               SCIPdebugMessage("infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(candidate));
               *result = SCIP_CUTOFF;

               return SCIP_OKAY;
            }

            ++(*nchgvartypes);
            *result = SCIP_SUCCESS;
            havechange = TRUE;
         }
      }

      /* call upgrade methods again if constraint has been changed */
      if( havechange )
      {
         SCIP_Bool upgraded;

         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );
         if( upgraded )
         {
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      consdata->ispresolved = TRUE;
   }

   /* if we did not try reformulations, ensure that presolving is called again even if there were only a few changes (< abortfac) */
   if( !doreformulations )
      *result = SCIP_DELAYED;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool      haslb;
   SCIP_Bool      hasub;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslb = !SCIPisInfinity(scip, -consdata->lhs);
   hasub = !SCIPisInfinity(scip, consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( consdata->lincoefs[i] > 0 )
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* @todo try to be more clever, but variable locks that depend on the bounds of other variables are not trival to maintain */
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->quadvarterms[i].var, nlockspos+nlocksneg, nlockspos+nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nlinvars == 0 && consdata->nquadvars == 0 )
   {
      SCIPinfoMessage(scip, file, "0 ");
   }
   else
   {
      SCIP_VAR*** monomialvars;
      SCIP_Real** monomialexps;
      SCIP_Real*  monomialcoefs;
      int*        monomialnvars;
      int         nmonomials;
      int         monomialssize;
      int         j;

      monomialssize = consdata->nlinvars + 2 * consdata->nquadvars + consdata->nbilinterms;
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars,  monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialexps,  monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialcoefs, monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialnvars, monomialssize) );

      nmonomials = 0;
      for( j = 0; j < consdata->nlinvars; ++j )
      {
         assert(nmonomials < monomialssize);

         SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */

         monomialvars[nmonomials][0] = consdata->linvars[j];
         monomialexps[nmonomials] = NULL;
         monomialcoefs[nmonomials] = consdata->lincoefs[j];
         monomialnvars[nmonomials] = 1;
         ++nmonomials;
      }

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         if( consdata->quadvarterms[j].lincoef != 0.0 )
         {
            assert(nmonomials < monomialssize);

            SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */

            monomialvars[nmonomials][0] = consdata->quadvarterms[j].var;
            monomialexps[nmonomials] = NULL;
            monomialcoefs[nmonomials] = consdata->quadvarterms[j].lincoef;
            monomialnvars[nmonomials] = 1;
            ++nmonomials;
         }

         if( consdata->quadvarterms[j].sqrcoef != 0.0 )
         {
            assert(nmonomials < monomialssize);

            SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */
            SCIP_CALL( SCIPallocBufferArray(scip, &monomialexps[nmonomials], 1) );  /*lint !e866 */

            monomialvars[nmonomials][0] = consdata->quadvarterms[j].var;
            monomialexps[nmonomials][0] = 2.0;
            monomialcoefs[nmonomials] = consdata->quadvarterms[j].sqrcoef;
            monomialnvars[nmonomials] = 1;
            ++nmonomials;
         }
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         assert(nmonomials < monomialssize);

         SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 2) );  /*lint !e866 */

         monomialvars[nmonomials][0] = consdata->bilinterms[j].var1;
         monomialvars[nmonomials][1] = consdata->bilinterms[j].var2;
         monomialexps[nmonomials] = NULL;
         monomialcoefs[nmonomials] = consdata->bilinterms[j].coef;
         monomialnvars[nmonomials] = 2;
         ++nmonomials;
      }

      SCIP_CALL( SCIPwriteVarsPolynomial(scip, file, monomialvars, monomialexps, monomialcoefs, monomialnvars, nmonomials, TRUE) );

      for( j = 0; j < nmonomials; ++j )
      {
         SCIPfreeBufferArray(scip, &monomialvars[j]);
         SCIPfreeBufferArrayNull(scip, &monomialexps[j]);
      }

      SCIPfreeBufferArray(scip, &monomialvars);
      SCIPfreeBufferArray(scip, &monomialexps);
      SCIPfreeBufferArray(scip, &monomialcoefs);
      SCIPfreeBufferArray(scip, &monomialnvars);
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   }
   else
   {
      /* should be ignored by parser */
      SCIPinfoMessage(scip, file, " [free]");
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;
   SCIP_Bool          maypropfeasible; /* whether we may be able to propose a feasible solution */

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   maxviol = 0.0;
   maypropfeasible = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL) &&
      SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, conshdlrdata->doscaling) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g (scaled: %.15g)\n", consdata->lhs - consdata->activity, consdata->lhsviol);
            }
            if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g (scaled: %.15g)\n", consdata->activity - consdata->rhs, consdata->rhsviol);
            }
         }
         if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible )
            return SCIP_OKAY;
         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = consdata->lhsviol + consdata->rhsviol;

         /* do not try to shift linear variables if activity is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, REALABS(consdata->activity)) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            /* update information on linear variables that may be in- or decreased */
            if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
               consdataFindUnlockedLinearVar(scip, consdata);

            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
   }

   if( *result == SCIP_INFEASIBLE && maypropfeasible )
   {
      SCIP_Bool success;

      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol, &success) );

      /* do not pass solution to NLP heuristic if we made it feasible this way */
      if( success )
         return SCIP_OKAY;
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyQuadratic)
{
   SCIP_CONSDATA*    consdata;
   SCIP_CONSDATA*    targetconsdata;
   SCIP_VAR**        linvars;
   SCIP_QUADVARTERM* quadvarterms;
   SCIP_BILINTERM*   bilinterms;
   int               i;
   int               j;
   int               k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   linvars = NULL;
   quadvarterms = NULL;
   bilinterms = NULL;

   *valid = TRUE;

   if( consdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->nlinvars) );
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->linvars[i], &linvars[i], varmap, consmap, global, valid) );
         assert(!(*valid) || linvars[i] != NULL);

         /* we do not copy, if a variable is missing */
         if( !(*valid) )
            goto TERMINATE;  
      }
   }

   if( consdata->nbilinterms != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinterms, consdata->nbilinterms) );
   }

   if( consdata->nquadvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &quadvarterms, consdata->nquadvars) );
      for( i = 0; i < consdata->nquadvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->quadvarterms[i].var, &quadvarterms[i].var, varmap, consmap, global, valid) );
         assert(!(*valid) || quadvarterms[i].var != NULL);

         /* we do not copy, if a variable is missing */
         if( !(*valid) )
            goto TERMINATE;

         quadvarterms[i].lincoef   = consdata->quadvarterms[i].lincoef;
         quadvarterms[i].sqrcoef   = consdata->quadvarterms[i].sqrcoef;
         quadvarterms[i].eventdata = NULL;
         quadvarterms[i].nadjbilin = consdata->quadvarterms[i].nadjbilin;
         quadvarterms[i].adjbilin  = consdata->quadvarterms[i].adjbilin;

         assert(consdata->nbilinterms != 0 || consdata->quadvarterms[i].nadjbilin == 0);

         for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
         {
            assert(bilinterms != NULL);

            k = consdata->quadvarterms[i].adjbilin[j];
            assert(consdata->bilinterms[k].var1 != NULL);
            assert(consdata->bilinterms[k].var2 != NULL);
            if( consdata->bilinterms[k].var1 == consdata->quadvarterms[i].var )
            {
               assert(consdata->bilinterms[k].var2 != consdata->quadvarterms[i].var);
               bilinterms[k].var1 = quadvarterms[i].var;
            }
            else
            {
               assert(consdata->bilinterms[k].var2 == consdata->quadvarterms[i].var);
               bilinterms[k].var2 = quadvarterms[i].var;
            }
            bilinterms[k].coef = consdata->bilinterms[k].coef; 
         }
      }
   }

   assert(stickingatnode == FALSE);
   SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name ? name : SCIPconsGetName(sourcecons),
         consdata->nlinvars, linvars, consdata->lincoefs,
         consdata->nquadvars, quadvarterms,
         consdata->nbilinterms, bilinterms,
         consdata->lhs, consdata->rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* copy information on curvature */
   targetconsdata = SCIPconsGetData(*cons);
   targetconsdata->isconvex      = consdata->isconvex;
   targetconsdata->isconcave     = consdata->isconcave;
   targetconsdata->iscurvchecked = consdata->iscurvchecked;

 TERMINATE:
   SCIPfreeBufferArrayNull(sourcescip, &quadvarterms);
   SCIPfreeBufferArrayNull(sourcescip, &bilinterms);
   SCIPfreeBufferArrayNull(sourcescip, &linvars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseQuadratic)
{  /*lint --e{715}*/
   SCIP_VAR*** monomialvars;
   SCIP_Real** monomialexps;
   SCIP_Real*  monomialcoefs;
   char*       endptr;
   int*        monomialnvars;
   int         nmonomials;

   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   (*success) = FALSE;

   /* return of string empty */
   if( !*str )
      return SCIP_OKAY;

   /* ignore whitespace */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, &endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_OKAY;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the first coefficient, but not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   SCIP_CALL( SCIPparseVarsPolynomial(scip, str, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, &nmonomials, &endptr, success) );

   if( *success )
   {
      /* check for right hand side */
      str = endptr;

      /* ignore whitespace */
      while( isspace((unsigned char)*str) )
         ++str;

      if( *str && str[0] == '<' && str[1] == '=' )
      {
         /* we seem to get a right-hand-side */
         str += 2;

         if( !SCIPstrToRealValue(str, &rhs, &endptr) )
         {
            SCIPerrorMessage("error parsing right-hand-side from %s\n", str);
            *success = FALSE;
         }
      }
      else if( *str && str[0] == '>' && str[1] == '=' )
      {
         /* we seem to get a left-hand-side */
         str += 2;

         /* we should not have a left-hand-side already */
         assert(SCIPisInfinity(scip, -lhs));

         if( !SCIPstrToRealValue(str, &lhs, &endptr) )
         {
            SCIPerrorMessage("error parsing left-hand-side from %s\n", str);
            *success = FALSE;
         }
      }
      else if( *str && str[0] == '=' && str[1] == '=' )
      {
         /* we seem to get a left- and right-hand-side */
         str += 2;

         /* we should not have a left-hand-side already */
         assert(SCIPisInfinity(scip, -lhs));

         if( !SCIPstrToRealValue(str, &lhs, &endptr) )
         {
            SCIPerrorMessage("error parsing left-hand-side from %s\n", str);
            *success = FALSE;
         }
         else
         {
            rhs = lhs;
         }
      }
   }

   if( *success )
   {
      int i;

      /* setup constraint */
      assert(stickingatnode == FALSE);
      SCIP_CALL( SCIPcreateConsQuadratic(scip, cons, name, 0, NULL, NULL,
            0, NULL, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

      for( i = 0; i < nmonomials; ++i )
      {
         if( monomialnvars[i] == 0 )
         {
            /* constant monomial */
            SCIPaddConstantQuadratic(scip, *cons, monomialcoefs[i]);
         }
         else if( monomialnvars[i] == 1 && monomialexps[i][0] == 1.0 )
         {
            /* linear monomial */
            SCIP_CALL( SCIPaddLinearVarQuadratic(scip, *cons, monomialvars[i][0], monomialcoefs[i]) );
         }
         else if( monomialnvars[i] == 1 && monomialexps[i][0] == 2.0 )
         {
            /* square monomial */
            SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, monomialvars[i][0], 0.0, monomialcoefs[i]) );
         }
         else if( monomialnvars[i] == 2 && monomialexps[i][0] == 1.0 && monomialexps[i][1] == 1.0 )
         {
            /* bilinear term */
            SCIP_VAR* var1;
            SCIP_VAR* var2;
            int pos;

            var1 = monomialvars[i][0];
            var2 = monomialvars[i][1];
            if( var1 == var2 )
            {
               /* actually a square term */
               SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var1, 0.0, monomialcoefs[i]) );
            }
            else
            {
               SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, var1, &pos) );
               if( pos == -1 )
               {
                  SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var1, 0.0, 0.0) );
               }

               SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, var2, &pos) );
               if( pos == -1 )
               {
                  SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var2, 0.0, 0.0) );
               }
            }

            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, *cons, var1, var2, monomialcoefs[i]) );
         }
         else
         {
            SCIPerrorMessage("polynomial in quadratic constraint does not have degree at most 2\n");
            *success = FALSE;
            SCIP_CALL( SCIPreleaseCons(scip, cons) );
            break;
         }
      }
   }

   SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nlinvars + consdata->nquadvars )
      (*success) = FALSE;
   else
   {
      int i;

      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->linvars, consdata->nlinvars);

      for( i = 0; i < consdata->nquadvars; ++i )
         vars[consdata->nlinvars+i] = consdata->quadvarterms[i].var;

      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nlinvars + consdata->nquadvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for quadratic constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create quadratic constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpQuadratic, consEnfopsQuadratic, consCheckQuadratic, consLockQuadratic,
         conshdlrdata) );
   assert(conshdlr != NULL);


   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyQuadratic, consCopyQuadratic) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolQuadratic) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeQuadratic) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsQuadratic) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpQuadratic) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseQuadratic) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolQuadratic, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintQuadratic) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropQuadratic, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpQuadratic, consSepasolQuadratic, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransQuadratic) );

   /* add quadratic constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/replacebinaryprod",
         "max. length of linear term which when multiplied with a binary variables is replaced by an auxiliary variable and a linear reformulation (0 to turn off)",
         &conshdlrdata->replacebinaryprodlength, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/empathy4and",
         "empathy level for using the AND constraint handler: 0 always avoid using AND; 1 use AND sometimes; 2 use AND as often as possible",
         &conshdlrdata->empathy4and, FALSE, 0, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/binreforminitial",
         "whether to make non-varbound linear constraints added due to replacing products with binary variables initial",
         &conshdlrdata->binreforminitial, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacysepa",
         "minimal efficacy for a cut to be added to the LP during separation; overwrites separating/efficacy",
         &conshdlrdata->mincutefficacysepa, TRUE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacyenfofac",
         "minimal target efficacy of a cut in order to add it to relaxation during enforcement as a factor of the feasibility tolerance (may be ignored)",
         &conshdlrdata->mincutefficacyenfofac, TRUE, 2.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/scaling", 
         "whether a quadratic constraint should be scaled w.r.t. the current gradient norm when checking for feasibility",
         &conshdlrdata->doscaling, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, TRUE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linearizeheursol",
         "whether linearizations of convex quadratic constraints should be added to cutpool in a solution found by some heuristic",
         &conshdlrdata->linearizeheursol, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/checkcurvature",
         "whether multivariate quadratic functions should be checked for convexity/concavity",
         &conshdlrdata->checkcurvature, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/checkfactorable",
         "whether constraint functions should be checked to be factorable",
         &conshdlrdata->checkfactorable, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting a linear variable (esp. useful if constraint was actually objective function)",
         &conshdlrdata->linfeasshift, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/disaggregate",
         "whether to disaggregate quadratic parts that decompose into a sum of non-overlapping quadratic terms",
         &conshdlrdata->disaggregate, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxproprounds",
         "limit on number of propagation rounds for a single constraint within one round of SCIP propagation during solve",
         &conshdlrdata->maxproprounds, TRUE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxproproundspresolve",
         "limit on number of propagation rounds for a single constraint within one round of SCIP presolve",
         &conshdlrdata->maxproproundspresolve, TRUE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/enfocutsremovable",
         "are cuts added during enforcement removable from the LP in the same node?",
         &conshdlrdata->enfocutsremovable, TRUE, FALSE, NULL, NULL) );

   conshdlrdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->eventhdlr),CONSHDLR_NAME"_boundchange", "signals a bound change to a quadratic constraint",
         processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
         processNewSolutionEvent, NULL) );

   /* include the quadratic constraint upgrade in the nonlinear constraint handler */
   SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdQuadratic, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );

   return SCIP_OKAY;
}

/** includes a quadratic constraint update method into the quadratic constraint handler */
SCIP_RETCODE SCIPincludeQuadconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method be active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_QUADCONSUPGRADE* quadconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(quadconsupgd != NULL);
   assert(conshdlrname != NULL );

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdataHasUpgrade(scip, conshdlrdata, quadconsupgd, conshdlrname) )
   {
      /* create a quadratic constraint upgrade data object */
      SCIP_CALL( SCIPallocMemory(scip, &quadconsupgrade) );
      quadconsupgrade->quadconsupgd = quadconsupgd;
      quadconsupgrade->priority     = priority;
      quadconsupgrade->active       = active;

      /* insert quadratic constraint upgrade method into constraint handler data */
      assert(conshdlrdata->nquadconsupgrades <= conshdlrdata->quadconsupgradessize);
      if( conshdlrdata->nquadconsupgrades+1 > conshdlrdata->quadconsupgradessize )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nquadconsupgrades+1);
         SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->quadconsupgrades, newsize) );
         conshdlrdata->quadconsupgradessize = newsize;
      }
      assert(conshdlrdata->nquadconsupgrades+1 <= conshdlrdata->quadconsupgradessize);

      for( i = conshdlrdata->nquadconsupgrades; i > 0 && conshdlrdata->quadconsupgrades[i-1]->priority < quadconsupgrade->priority; --i )
         conshdlrdata->quadconsupgrades[i] = conshdlrdata->quadconsupgrades[i-1];
      assert(0 <= i && i <= conshdlrdata->nquadconsupgrades);
      conshdlrdata->quadconsupgrades[i] = quadconsupgrade;
      conshdlrdata->nquadconsupgrades++;

      /* adds parameter to turn on and off the upgrading step */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/"CONSHDLR_NAME"/upgrade/%s", conshdlrname);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable quadratic upgrading for constraint handler <%s>", conshdlrname);
      SCIP_CALL( SCIPaddBoolParam(scip,
            paramname, paramdesc,
            &quadconsupgrade->active, FALSE, active, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 * 
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP*  quadvaridxs;
   SCIP_Real      sqrcoef;
   int i;
   int var1pos;
   int var2pos;

   int nbilinterms;

   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data and constraint */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   /* add quadratic variables and remember their indices */
   SCIP_CALL( SCIPhashmapCreate(&quadvaridxs, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nquadterms)) );
   nbilinterms = 0;
   for( i = 0; i < nquadterms; ++i )
   {
      if( SCIPisZero(scip, quadcoefs[i]) )
         continue;

      /* if it is actually a square term, remember it's coefficient */
      if( quadvars1[i] == quadvars2[i] )
         sqrcoef = quadcoefs[i];
      else
         sqrcoef = 0.0;

      /* add quadvars1[i], if not in there already */
      if( !SCIPhashmapExists(quadvaridxs, quadvars1[i]) )
      {
         SCIP_CALL( addQuadVarTerm(scip, *cons, quadvars1[i], 0.0, sqrcoef, FALSE) );
         assert(consdata->nquadvars >= 0);
         assert(consdata->quadvarterms[consdata->nquadvars-1].var == quadvars1[i]);

         SCIP_CALL( SCIPhashmapInsert(quadvaridxs, quadvars1[i], (void*)(size_t)(consdata->nquadvars-1)) );
      }
      else if( !SCIPisZero(scip, sqrcoef) )
      {
         /* if it's there already, but we got a square coefficient, add it to the previous one */ 
         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars1[i]);
         assert(consdata->quadvarterms[var1pos].var == quadvars1[i]);
         consdata->quadvarterms[var1pos].sqrcoef += sqrcoef;
      }

      if( quadvars1[i] == quadvars2[i] )
         continue;

      /* add quadvars2[i], if not in there already */
      if( !SCIPhashmapExists(quadvaridxs, quadvars2[i]) )
      {
         assert(sqrcoef == 0.0);
         SCIP_CALL( addQuadVarTerm(scip, *cons, quadvars2[i], 0.0, 0.0, FALSE) );
         assert(consdata->nquadvars >= 0);
         assert(consdata->quadvarterms[consdata->nquadvars-1].var == quadvars2[i]);

         SCIP_CALL( SCIPhashmapInsert(quadvaridxs, quadvars2[i], (void*)(size_t)(consdata->nquadvars-1)) );
      }

      ++nbilinterms;
   }

   /* add bilinear terms, if we saw any */
   if( nbilinterms > 0 )
   {
      SCIP_CALL( consdataEnsureBilinSize(scip, consdata, nbilinterms) );
      for( i = 0; i < nquadterms; ++i )
      {
         if( SCIPisZero(scip, quadcoefs[i]) )
            continue;

         /* square terms have been taken care of already */
         if( quadvars1[i] == quadvars2[i] )
            continue;

         assert(SCIPhashmapExists(quadvaridxs, quadvars1[i]));
         assert(SCIPhashmapExists(quadvaridxs, quadvars2[i]));

         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars1[i]);
         var2pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars2[i]);

         SCIP_CALL( addBilinearTerm(scip, *cons, var1pos, var2pos, quadcoefs[i]) );
      }
   }

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )
         continue;

      /* if it's a linear coefficient for a quadratic variable, add it there, otherwise add as linear variable */
      if( SCIPhashmapExists(quadvaridxs, linvars[i]) )
      {
         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, linvars[i]);
         assert(consdata->quadvarterms[var1pos].var == linvars[i]);
         consdata->quadvarterms[var1pos].lincoef += lincoefs[i];
      }
      else
      {
         SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );
      }
   }

   if( SCIPisTransformed(scip) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );
   }

   SCIPhashmapFree(&quadvaridxs);

   /* merge duplicate bilinear terms, move quad terms that are linear to linear vars */
   SCIP_CALL( mergeAndCleanBilinearTerms(scip, *cons) );
   SCIP_CALL( mergeAndCleanQuadVarTerms(scip, *cons) );
   SCIP_CALL( mergeAndCleanLinearVars(scip, *cons) );

   SCIPdebugMessage("created quadratic constraint ");
   SCIPdebugPrintCons(scip, *cons, NULL);

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint with all its
 *  flags set to their default values.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation (u) */
   )
{
   SCIP_CALL( SCIPcreateConsQuadratic(scip, cons, name, nlinvars, linvars, lincoefs,
         nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 * 
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */ 
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */ 
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs,                /**< constraint right hand side (u) */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(modifiable == FALSE); /* we do not support column generation */
   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadvarterms == 0 || quadvarterms != NULL);
   assert(nbilinterms == 0 || bilinterms != NULL);

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, lhs, rhs,
         nlinvars, linvars, lincoefs, nquadvarterms, quadvarterms, nbilinterms, bilinterms,
         TRUE) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   if( SCIPisTransformed(scip) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );
   }

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint in its most basic version, i.e.,
 *  all constraint flags are set to their default values.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs                 /**< constraint right hand side (u) */
   )
{
   SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name, nlinvars, linvars, lincoefs,
         nquadvarterms, quadvarterms, nbilinterms, bilinterms, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** Adds a constant to the constraint function, that is, subtracts a constant from both sides */
void SCIPaddConstantQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             constant            /**< constant to subtract from both sides */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lhs <= consdata->rhs);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhs -= constant;
   if( !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhs -= constant;

   if( consdata->lhs > consdata->rhs )
   {
      assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));
      consdata->lhs = consdata->rhs;
   }
}

/** Adds a linear variable with coefficient to a quadratic constraint.
 */
SCIP_RETCODE SCIPaddLinearVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   SCIP_CALL( addLinearCoef(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** Adds a quadratic variable with linear and square coefficient to a quadratic constraint.
 */
SCIP_RETCODE SCIPaddQuadVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef             /**< square coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(lincoef)));
   assert(!SCIPisInfinity(scip, REALABS(sqrcoef)));

   SCIP_CALL( addQuadVarTerm(scip, cons, var, lincoef, sqrcoef, SCIPconsIsTransformed(cons)) );

   return SCIP_OKAY;
}

/** Adds a linear coefficient for a quadratic variable.
 *
 * Variable will be added with square coefficient 0.0 if not existing yet.
 */
SCIP_RETCODE SCIPaddQuadVarLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to linear coefficient of variable */
   )
{
   SCIP_CONSDATA* consdata;
   int pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var, &pos) );
   if( pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var, coef, 0.0, SCIPconsIsTransformed(cons)) );
      return SCIP_OKAY;
   }
   assert(pos < consdata->nquadvars);
   assert(consdata->quadvarterms[pos].var == var);

   consdata->quadvarterms[pos].lincoef += coef;

   /* update flags and invalid activities */
   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = consdata->ispresolved && !SCIPisZero(scip, consdata->quadvarterms[pos].lincoef);

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** Adds a square coefficient for a quadratic variable.
 *
 * Variable will be added with linear coefficient 0.0 if not existing yet.
 */
SCIP_RETCODE SCIPaddSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to square coefficient of variable */
   )
{
   SCIP_CONSDATA* consdata;
   int pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var, &pos) );
   if( pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var, 0.0, coef, SCIPconsIsTransformed(cons)) );
      return SCIP_OKAY;
   }
   assert(pos < consdata->nquadvars);
   assert(consdata->quadvarterms[pos].var == var);

   consdata->quadvarterms[pos].sqrcoef += coef;

   /* update flags and invalid activities */
   consdata->isconvex      = FALSE;
   consdata->isconcave     = FALSE;
   consdata->iscurvchecked = FALSE;
   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = consdata->ispresolved && !SCIPisZero(scip, consdata->quadvarterms[pos].sqrcoef);

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** Adds a bilinear term to a quadratic constraint.
 *
 * Variables will be added with linear and square coefficient 0.0 if not existing yet.
 * If variables are equal, only the square coefficient of the variable is updated.
 */
SCIP_RETCODE SCIPaddBilinTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_CONSDATA* consdata;
   int            var1pos;
   int            var2pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( var1 == var2 )
   {
      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, var1, coef) );
      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var1, &var1pos) );
   if( var1pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var1, 0.0, 0.0, SCIPconsIsTransformed(cons)) );
      var1pos = consdata->nquadvars-1;
   }

   if( !consdata->quadvarssorted )
   {
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );
      /* sorting may change the position of var1 */
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var1, &var1pos) );
      assert(var1pos >= 0);
   }

   assert(consdata->quadvarssorted);
   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var2, &var2pos) );
   if( var2pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var2, 0.0, 0.0, SCIPconsIsTransformed(cons)) );
      var2pos = consdata->nquadvars-1;
   }

   assert(consdata->quadvarterms[var1pos].var == var1);
   assert(consdata->quadvarterms[var2pos].var == var2);

   SCIP_CALL( addBilinearTerm(scip, cons, var1pos, var2pos, coef) );

   return SCIP_OKAY;
}

/** Gets the quadratic constraint as a nonlinear row representation.
 */
SCIP_RETCODE SCIPgetNlRowQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** Gets the number of variables in the linear term of a quadratic constraint.
 */
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nlinvars;
}

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->linvars;
}

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_Real* SCIPgetCoefsLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lincoefs;
}

/** Gets the number of quadratic variable terms of a quadratic constraint.
 */
int SCIPgetNQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nquadvars;
}

/** Gets the quadratic variable terms of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarTermsQuadratic.
 */
SCIP_QUADVARTERM* SCIPgetQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->quadvarterms;
}

/** Ensures that quadratic variable terms are sorted. */
SCIP_RETCODE SCIPsortQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   SCIP_CALL( consdataSortQuadVarTerms(scip, SCIPconsGetData(cons)) );

   return SCIP_OKAY;
}

/** Finds the position of a quadratic variable term for a given variable.
 *
 * @note If the quadratic variable terms have not been sorted before, then a search may reorder the current order of the terms.
 */
SCIP_RETCODE SCIPfindQuadVarTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to search for */
   int*                  pos                 /**< buffer to store position of quadvarterm for var, or -1 if not found */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(var != NULL);
   assert(pos != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, SCIPconsGetData(cons), var, pos) );

   return SCIP_OKAY;
}

/** Gets the number of bilinear terms of a quadratic constraint.
 */
int SCIPgetNBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nbilinterms;
}

/** Gets the bilinear terms of a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_BILINTERM* SCIPgetBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->bilinterms;
}

/** Gets the left hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lhs;
}

/** Gets the right hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhs;
}

/** Check the quadratic function of a quadratic constraint for its semi-definiteness, if not done yet.
 */
SCIP_RETCODE SCIPcheckCurvatureQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkCurvature(scip, cons, TRUE) );

   return SCIP_OKAY;
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
 */
SCIP_Bool SCIPisConvexQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_Bool determined;

   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   checkCurvatureEasy(scip, cons, &determined, FALSE);
   assert(determined);

   return (SCIPconsGetData(cons)->isconvex);
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
 */
SCIP_Bool SCIPisConcaveQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_Bool determined;

   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   checkCurvatureEasy(scip, cons, &determined, FALSE);
   assert(determined);

   return (SCIPconsGetData(cons)->isconcave);
}

/** Computes the violation of a constraint by a solution */
SCIP_RETCODE SCIPgetViolationQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< buffer to store violation of constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violation != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, TRUE) ); /* we assume that scaling was left on */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violation = MAX(consdata->lhsviol, consdata->rhsviol);

   return SCIP_OKAY;
}

/** Indicates whether the quadratic constraint is local w.r.t. the current local bounds.
 *
 * That is, checks whether each variable with a square term is fixed and for each bilinear term at least one variable is fixed.
 */
SCIP_Bool SCIPisLinearLocalQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check all square terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].sqrcoef == 0.0 )
         continue;

      var1 = consdata->quadvarterms[i].var;
      assert(var1 != NULL);

      if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) )
         return FALSE;
   }

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      var1 = consdata->bilinterms[i].var1;
      var2 = consdata->bilinterms[i].var2;

      assert(var1 != NULL);
      assert(var2 != NULL);

      if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) &&
         ! SCIPisRelEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
         return FALSE;
   }

   return TRUE;
}

/** Adds the constraint to an NLPI problem. */
SCIP_RETCODE SCIPaddToNlpiProblemQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< NLPI problem where to add constraint */
   SCIP_HASHMAP*         scipvar2nlpivar,    /**< mapping from SCIP variables to variable indices in NLPI */
   SCIP_Bool             names               /**< whether to pass constraint names to NLPI */
   )
{
   SCIP_CONSDATA* consdata;
   int            nlininds;
   int*           lininds;
   SCIP_Real*     linvals;
   int            nquadelems;
   SCIP_QUADELEM* quadelems;
   SCIP_VAR*      othervar;
   const char*    name;
   int            j;
   int            l;
   int            lincnt;
   int            quadcnt;
   int            idx1;
   int            idx2;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(scipvar2nlpivar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* count nonzeros in quadratic part */
   nlininds = consdata->nlinvars;
   nquadelems = consdata->nbilinterms;
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      if( consdata->quadvarterms[j].sqrcoef != 0.0 )
         ++nquadelems;
      if( consdata->quadvarterms[j].lincoef != 0.0 )
         ++nlininds;
   }

   /* setup linear part */
   lininds = NULL;
   linvals = NULL;
   lincnt  = 0;
   if( nlininds > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlininds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlininds) );

      for( j = 0; j < consdata->nlinvars; ++j )
      {
         linvals[j] = consdata->lincoefs[j];
         assert(SCIPhashmapExists(scipvar2nlpivar, consdata->linvars[j]));
         lininds[j] = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->linvars[j]);
      }

      lincnt = consdata->nlinvars;
   }

   /* setup quadratic part */
   quadelems = NULL;
   if( nquadelems > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );
   }
   quadcnt = 0;

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      assert(SCIPhashmapExists(scipvar2nlpivar, consdata->quadvarterms[j].var));
      idx1 = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpivar, consdata->quadvarterms[j].var);
      if( consdata->quadvarterms[j].lincoef != 0.0 )
      {
         assert(lininds != NULL);
         assert(linvals != NULL);
         lininds[lincnt] = idx1;
         linvals[lincnt] = consdata->quadvarterms[j].lincoef;
         ++lincnt;
      }

      if( consdata->quadvarterms[j].sqrcoef != 0.0 )
      {
         assert(quadcnt < nquadelems);
         assert(quadelems != NULL);
         quadelems[quadcnt].idx1 = idx1;
         quadelems[quadcnt].idx2 = idx1;
         quadelems[quadcnt].coef = consdata->quadvarterms[j].sqrcoef;
         ++quadcnt;
      }

      for( l = 0; l < consdata->quadvarterms[j].nadjbilin; ++l )
      {
         othervar = consdata->bilinterms[consdata->quadvarterms[j].adjbilin[l]].var2;
         /* if othervar is on position 2, then we process this bilinear term later (or it was processed already) */
         if( othervar == consdata->quadvarterms[j].var )
            continue;

         assert(quadcnt < nquadelems);
         assert(quadelems != NULL);
         assert(SCIPhashmapExists(scipvar2nlpivar, othervar));
         idx2 = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpivar, othervar);
         quadelems[quadcnt].idx1 = MIN(idx1, idx2);
         quadelems[quadcnt].idx2 = MAX(idx1, idx2);
         quadelems[quadcnt].coef = consdata->bilinterms[consdata->quadvarterms[j].adjbilin[l]].coef;
         ++quadcnt;
      }
   }

   assert(quadcnt == nquadelems);
   assert(lincnt  == nlininds);

   name = names ? SCIPconsGetName(cons) : NULL;

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1,
         &consdata->lhs, &consdata->rhs,
         &nlininds, &lininds, &linvals ,
         &nquadelems, &quadelems,
         NULL, NULL, &name) );

   SCIPfreeBufferArrayNull(scip, &quadelems);
   SCIPfreeBufferArrayNull(scip, &lininds);
   SCIPfreeBufferArrayNull(scip, &linvals);

   return SCIP_OKAY;
}
