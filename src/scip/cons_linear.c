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

/**@file   cons_linear.c
 * @brief Constraint handler for linear constraints in their most general form, \f$lhs <= a^T x <= rhs\f$.
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Michael Winkler
 * @author Gerald Gamrath
 *
 *  Linear constraints are separated with a high priority, because they are easy
 *  to separate. Instead of using the global cut pool, the same effect can be
 *  implemented by adding linear constraints to the root node, such that they are
 *  separated each time, the linear constraints are separated. A constraint
 *  handler, which generates linear constraints in this way should have a lower
 *  separation priority than the linear constraint handler, and it should have a
 *  separation frequency that is a multiple of the frequency of the linear
 *  constraint handler. In this way, it can be avoided to separate the same cut
 *  twice, because if a separation run of the handler is always preceded by a
 *  separation of the linear constraints, the priorily added constraints are
 *  always satisfied.
 *
 *  Linear constraints are enforced and checked with a very low priority. Checking
 *  of (many) linear constraints is much more involved than checking the solution
 *  values for integrality. Because we are separating the linear constraints quite
 *  often, it is only necessary to enforce them for integral solutions. A constraint
 *  handler which generates pool cuts in its enforcing method should have an
 *  enforcing priority smaller than that of the linear constraint handler to avoid
 *  regenerating constraints which already exist.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc.h"

#define CONSHDLR_NAME          "linear"
#define CONSHDLR_DESC          "linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY   +100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"

#define CONFLICTHDLR_NAME      "linear"
#define CONFLICTHDLR_DESC      "conflict handler creating linear constraints"
#define CONFLICTHDLR_PRIORITY  -1000000

#define DEFAULT_TIGHTENBOUNDSFREQ       1 /**< multiplier on propagation frequency, how often the bounds are tightened */
#define DEFAULT_MAXROUNDS               5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT          -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS            50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT       200 /**< maximal number of cuts separated per separation round in root node */
#define DEFAULT_PRESOLPAIRWISE       TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_PRESOLUSEHASHING     TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define DEFAULT_NMINCOMPARISONS    200000 /**< number for minimal pairwise presolving comparisons */
#define DEFAULT_MINGAINPERNMINCOMP  1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise 
                                           *   comparison round */
#define DEFAULT_SORTVARS             TRUE /**< should variables be sorted after presolve w.r.t their coefficient absolute for faster
                                           *  propagation? */
#define DEFAULT_MAXAGGRNORMSCALE      0.0 /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                           *   (0.0: disable constraint aggregation) */
#define DEFAULT_MAXCARDBOUNDDIST      0.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                           *   to best node's dual bound for separating knapsack cardinality cuts */
#define DEFAULT_SEPARATEALL         FALSE /**< should all constraints be subject to cardinality cut generation instead of only
                                           *   the ones with non-zero dual value? */
#define DEFAULT_AGGREGATEVARIABLES   TRUE /**< should presolving search for redundant variables in equations */
#define DEFAULT_SIMPLIFYINEQUALITIES TRUE/**< should presolving try to simplify inequalities */
#define DEFAULT_DUALPRESOLVING       TRUE /**< should dual presolving steps be performed? */

#define MAXDNOM                   10000LL /**< maximal denominator for simple rational fixed values */
#define MAXSCALEDCOEF               1e+03 /**< maximal coefficient value after scaling */
#define MAXSCALEDCOEFINTEGER        1e+05 /**< maximal coefficient value after scaling if all variables are of integral
                                           *   type
                                           */

#define HASHSIZE_LINEARCONS        131101 /**< minimal size of hash table in linear constraint tables */

#define QUADCONSUPGD_PRIORITY     1000000 /**< priority of the constraint handler for upgrading of quadratic constraints */
#define NONLINCONSUPGD_PRIORITY   1000000 /**< priority of the constraint handler for upgrading of nonlinear constraints */

#ifdef WITH_PRINTORIGCONSTYPES
/** constraint type */
enum SCIP_Constype
{
   SCIP_CONSTYPE_EMPTY         =  0,         /**<  */
   SCIP_CONSTYPE_FREE          =  1,         /**<  */
   SCIP_CONSTYPE_SINGLETON     =  2,         /**<  */
   SCIP_CONSTYPE_AGGREGATION   =  3,         /**<  */
   SCIP_CONSTYPE_VARBOUND      =  4,         /**<  */
   SCIP_CONSTYPE_SETPARTITION  =  5,         /**<  */
   SCIP_CONSTYPE_SETPACKING    =  6,         /**<  */
   SCIP_CONSTYPE_SETCOVERING   =  7,         /**<  */
   SCIP_CONSTYPE_CARDINALITY   =  8,         /**<  */
   SCIP_CONSTYPE_INVKNAPSACK   =  9,         /**<  */
   SCIP_CONSTYPE_EQKNAPSACK    = 10,         /**<  */
   SCIP_CONSTYPE_BINPACKING    = 11,         /**<  */
   SCIP_CONSTYPE_KNAPSACK      = 12,         /**<  */
   SCIP_CONSTYPE_INTKNAPSACK   = 13,         /**<  */
   SCIP_CONSTYPE_MIXEDBINARY   = 14,         /**<  */
   SCIP_CONSTYPE_GENERAL       = 15          /**<  */
};
typedef enum SCIP_Constype SCIP_CONSTYPE;
#endif


/** constraint data for linear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of row (for ranged rows) */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             maxabsval;          /**< maximum absolute value of all coefficients */
   SCIP_Real             pseudoactivity;     /**< pseudo activity value in current pseudo solution */
   SCIP_Real             lastpseudoactivity; /**< last pseudo activity which was computed by complete summation
                                              *   over all contributing values */
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
   SCIP_Longint          possignature;       /**< bit signature of coefficients that may take a positive value */
   SCIP_Longint          negsignature;       /**< bit signature of coefficients that may take a negative value */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_VAR**            vars;               /**< variables of constraint entries */
   SCIP_Real*            vals;               /**< coefficients of constraint entries */
   SCIP_EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   int                   pseudoactivityneginf;/**< number of coefficients contributing with neg. infinite value to pseudoactivity */
   int                   pseudoactivityposinf;/**< number of coefficients contributing with pos. infinite value to pseudoactivity */
   int                   pseudoactivityneghuge;/**< number of coefficients contributing with huge neg. value to pseudoactivity */
   int                   pseudoactivityposhuge;/**< number of coefficients contributing with huge pos. value to pseudoactivity */
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
   unsigned int          validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int          validactivities:1;  /**< are the pseudo activity and activity bounds (local and global) valid? */
   unsigned int          validpseudoact:1;   /**< is the pseudo activity valid? */
   unsigned int          validminact:1;      /**< is the local minactivity valid? */
   unsigned int          validmaxact:1;      /**< is the local maxactivity valid? */
   unsigned int          validglbminact:1;   /**< is the global minactivity valid? */
   unsigned int          validglbmaxact:1;   /**< is the global maxactivity valid? */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   unsigned int          boundstightened:1;  /**< is constraint already propagated with bound tightening? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          validsignature:1;   /**< is the bit signature valid? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int          donotupgrade:1;     /**< should the constraint not be upgraded? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
   unsigned int          binvarssorted:1;    /**< are binary variables sorted w.r.t. the absolute of their coefficient? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
};

/** event data for bound change event */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< linear constraint data to process the bound change for */
   int                   varpos;             /**< position of variable in vars array */
   int                   filterpos;          /**< position of event in variable's event filter */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_LINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   SCIP_Real             maxaggrnormscale;   /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                              *   (0.0: disable constraint aggregation) */
   SCIP_Real             maxcardbounddist;   /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for separating knapsack cardinality cuts */
   SCIP_Real             mingainpernmincomp; /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
   int                   linconsupgradessize;/**< size of linconsupgrade array */
   int                   nlinconsupgrades;   /**< number of linear constraint upgrade methods */
   int                   tightenboundsfreq;  /**< multiplier on propagation frequency, how often the bounds are tightened */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
   int                   nmincomparisons;    /**< number for minimal pairwise presolving comparisons */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             separateall;        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   SCIP_Bool             aggregatevariables; /**< should presolving search for redundant variables in equations */
   SCIP_Bool             simplifyinequalities;/**< should presolving try to cancel down or delete coefficients in inequalities */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be performed? */
   SCIP_Bool             sortvars;           /**< should binary variables be sorted for faster propagation? */
};

/** linear constraint update method */
struct SCIP_LinConsUpgrade
{
   SCIP_DECL_LINCONSUPGD((*linconsupgd));    /**< method to call for upgrading linear constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};


/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1_RHS        = 1,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the right hand side of the inequality */
   PROPRULE_1_LHS        = 2,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the left hand side of the inequality */
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

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
int inferInfoGetProprule(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.proprule;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a propagation rule and a position number */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   INFERINFO inferinfo;
   assert( pos >= 0 );

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


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that linconsupgrades array can store at least num entries */
static
SCIP_RETCODE conshdlrdataEnsureLinconsupgradesSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nlinconsupgrades <= conshdlrdata->linconsupgradessize);

   if( num > conshdlrdata->linconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->linconsupgrades, newsize) );
      conshdlrdata->linconsupgradessize = newsize;
   }
   assert(num <= conshdlrdata->linconsupgradessize);

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vals, consdata->varssize, newsize) );
      if( consdata->eventdatas != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize, newsize) );
      }
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}


/*
 * local methods for managing linear constraint update methods
 */

/** creates a linear constraint upgrade data object */
static
SCIP_RETCODE linconsupgradeCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade,     /**< pointer to store the linear constraint upgrade */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority            /**< priority of upgrading method */
   )
{
   assert(scip != NULL);
   assert(linconsupgrade != NULL);
   assert(linconsupgd != NULL);

   SCIP_CALL( SCIPallocMemory(scip, linconsupgrade) );
   (*linconsupgrade)->linconsupgd = linconsupgd;
   (*linconsupgrade)->priority = priority;
   (*linconsupgrade)->active = TRUE;

   return SCIP_OKAY;
}

/** frees a linear constraint upgrade data object */
static
void linconsupgradeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade      /**< pointer to the linear constraint upgrade */
   )
{
   assert(scip != NULL);
   assert(linconsupgrade != NULL);
   assert(*linconsupgrade != NULL);

   SCIPfreeMemory(scip, linconsupgrade);
}

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

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );
   (*conshdlrdata)->linconsupgrades = NULL;
   (*conshdlrdata)->linconsupgradessize = 0;
   (*conshdlrdata)->nlinconsupgrades = 0;

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
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   for( i = 0; i < (*conshdlrdata)->nlinconsupgrades; ++i )
   {
      linconsupgradeFree(scip, &(*conshdlrdata)->linconsupgrades[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &(*conshdlrdata)->linconsupgrades);

   SCIPfreeMemory(scip, conshdlrdata);
}

/** creates a linear constraint upgrade data object */
static
SCIP_Bool conshdlrdataHasUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(linconsupgd != NULL);
   assert(conshdlrname != NULL);

   for( i = conshdlrdata->nlinconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->linconsupgrades[i]->linconsupgd == linconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade message %p for constraint handler %s.\n", linconsupgd, conshdlrname);
#endif
         return TRUE;
      }
   }

   return FALSE;
}

/** adds a linear constraint update method to the constraint handler's data */
static
SCIP_RETCODE conshdlrdataIncludeUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_LINCONSUPGRADE*  linconsupgrade      /**< linear constraint upgrade method */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(linconsupgrade != NULL);

   SCIP_CALL( conshdlrdataEnsureLinconsupgradesSize(scip, conshdlrdata, conshdlrdata->nlinconsupgrades+1) );
   
   for( i = conshdlrdata->nlinconsupgrades;
        i > 0 && conshdlrdata->linconsupgrades[i-1]->priority < linconsupgrade->priority; --i )
   {
      conshdlrdata->linconsupgrades[i] = conshdlrdata->linconsupgrades[i-1];
   }
   assert(0 <= i && i <= conshdlrdata->nlinconsupgrades);
   conshdlrdata->linconsupgrades[i] = linconsupgrade;
   conshdlrdata->nlinconsupgrades++;

   return SCIP_OKAY;
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
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
static
SCIP_RETCODE consdataCatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);
   assert(consdata->vars[pos] != NULL);
   assert(SCIPvarIsTransformed(consdata->vars[pos]));
   assert(consdata->eventdatas != NULL);
   assert(consdata->eventdatas[pos] == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &(consdata->eventdatas[pos])) ); /*lint !e866*/
   consdata->eventdatas[pos]->consdata = consdata;
   consdata->eventdatas[pos]->varpos = pos;

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED | SCIP_EVENTTYPE_OBJCHANGED,
         eventhdlr, consdata->eventdatas[pos], &consdata->eventdatas[pos]->filterpos) );

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
SCIP_RETCODE consdataDropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);
   assert(consdata->eventdatas != NULL);
   assert(consdata->eventdatas[pos] != NULL);
   assert(consdata->eventdatas[pos]->consdata == consdata);
   assert(consdata->eventdatas[pos]->varpos == pos);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED | SCIP_EVENTTYPE_OBJCHANGED,
         eventhdlr, consdata->eventdatas[pos], consdata->eventdatas[pos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->eventdatas[pos]); /*lint !e866*/

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consdataCatchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->eventdatas == NULL);

   /* allocate eventdatas array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize) );
   assert(consdata->eventdatas != NULL);
   BMSclearMemoryArray(consdata->eventdatas, consdata->nvars);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consdataCatchEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consdataDropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->eventdatas != NULL);

   /* drop event of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consdataDropEvent(scip, consdata, eventhdlr, i) );
   }

   /* free eventdatas array */
   SCIPfreeBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize);
   assert(consdata->eventdatas == NULL);

   return SCIP_OKAY;
}

/** returns whether we are in a stage, where the variable events should be caught */
static
SCIP_Bool needEvents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return (SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) < SCIP_STAGE_FREETRANS);
}

/** creates a linear constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs                 /**< right hand side of row */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of linear constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%g, rhs=%g\n", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->varssize = nvars;
   (*consdata)->nvars = nvars;
   if( nvars > 0 )
   {
      int k;

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vals, vals, nvars) );
      k = 0;
      for( v = 0; v < nvars; ++v )
      {
         assert((*consdata)->vars[v] != NULL);
         if( !SCIPisZero(scip, (*consdata)->vals[v]) )
         {
            (*consdata)->vars[k] = (*consdata)->vars[v];
            (*consdata)->vals[k] = (*consdata)->vals[v];
            k++;
         }
      }
      (*consdata)->nvars = k;
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->vals = NULL;
   }
   (*consdata)->eventdatas = NULL;

   (*consdata)->row = NULL;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->maxabsval = SCIP_INVALID;
   (*consdata)->pseudoactivity = SCIP_INVALID;
   (*consdata)->minactivity = SCIP_INVALID;
   (*consdata)->maxactivity = SCIP_INVALID;
   (*consdata)->lastminactivity = SCIP_INVALID;
   (*consdata)->lastmaxactivity = SCIP_INVALID;
   (*consdata)->pseudoactivityneginf = -1;
   (*consdata)->pseudoactivityposinf = -1;
   (*consdata)->pseudoactivityneghuge = -1;
   (*consdata)->pseudoactivityposhuge = -1;
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
   (*consdata)->possignature = 0;
   (*consdata)->negsignature = 0;
   (*consdata)->validmaxabsval = FALSE;
   (*consdata)->validactivities = FALSE;
   (*consdata)->validpseudoact = FALSE;
   (*consdata)->validminact = FALSE;
   (*consdata)->validmaxact = FALSE;
   (*consdata)->validglbminact = FALSE;
   (*consdata)->validglbmaxact = FALSE;
   (*consdata)->propagated = FALSE;
   (*consdata)->boundstightened = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->removedfixings = FALSE;
   (*consdata)->validsignature = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->normalized = FALSE;
   (*consdata)->upgradetried = FALSE;
   (*consdata)->upgraded = FALSE;
   (*consdata)->donotupgrade = FALSE;
   (*consdata)->sorted = (nvars <= 1);
   (*consdata)->merged = (nvars <= 1);
   (*consdata)->cliquesadded = FALSE;
   (*consdata)->binvarssorted = FALSE;
   (*consdata)->nbinvars = -1;
   (*consdata)->varsdeleted = FALSE;

   if( SCIPisTransformed(scip) )
   {
      /* get transformed variables */
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* catch bound change events of variables */
      if( needEvents(scip) )
      {
         SCIP_CALL( consdataCatchAllEvents(scip, *consdata, eventhdlr) );
         assert((*consdata)->eventdatas != NULL);
      }
   }

   /* capture variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPisZero(scip, (*consdata)->vals[v]));
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
   }

   return SCIP_OKAY;
}

/** frees a linear constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->varssize >= 0);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* free event datas */
   if( (*consdata)->eventdatas != NULL )
   {
      /* drop bound change events of variables */
      SCIP_CALL( consdataDropAllEvents(scip, *consdata, eventhdlr) );
   }
   assert((*consdata)->eventdatas == NULL);

   /* release variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPisZero(scip, (*consdata)->vals[v]));
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vals, (*consdata)->varssize);
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
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0");
   else
   {
      /* post linear sum of the linear constraint */
      SCIP_CALL( SCIPwriteVarsLinearsum(scip, file, consdata->vars, consdata->vals, consdata->nvars, TRUE) );
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}

/** invalidates pseudo activity and activity bounds, such that they are recalculated in next get */
static
void consdataInvalidateActivities(
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   consdata->validactivities = FALSE;
   consdata->validpseudoact = FALSE;
   consdata->validminact = FALSE;
   consdata->validmaxact = FALSE;
   consdata->validglbminact = FALSE;
   consdata->validglbmaxact = FALSE;
   consdata->pseudoactivity = SCIP_INVALID;
   consdata->minactivity = SCIP_INVALID;
   consdata->maxactivity = SCIP_INVALID;
   consdata->lastminactivity = SCIP_INVALID;
   consdata->lastmaxactivity = SCIP_INVALID;
   consdata->pseudoactivityneginf = -1;
   consdata->pseudoactivityposinf = -1;
   consdata->pseudoactivityneghuge = -1;
   consdata->pseudoactivityposhuge = -1;
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

/** recompute the pseudo activity of a constraint */
static
void consdataRecomputePseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;

   consdata->pseudoactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (SCIPvarGetBestBoundType(consdata->vars[i]) == SCIP_BOUNDTYPE_LOWER) ? SCIPvarGetLbLocal(consdata->vars[i]) : SCIPvarGetUbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound) )
         consdata->pseudoactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validpseudoact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastpseudoactivity = consdata->pseudoactivity;
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

   consdata->minactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetLbLocal(consdata->vars[i]) : SCIPvarGetUbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->minactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastminactivity = consdata->minactivity;
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

   consdata->maxactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetUbLocal(consdata->vars[i]) : SCIPvarGetLbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->maxactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastmaxactivity = consdata->maxactivity;
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

   consdata->glbminactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetLbGlobal(consdata->vars[i]) : SCIPvarGetUbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound) )
         consdata->glbminactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbminactivity = consdata->glbminactivity;
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

   consdata->glbmaxactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetUbGlobal(consdata->vars[i]) : SCIPvarGetLbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound) )
         consdata->glbmaxactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
}


/** updates activities for a change in a bound */
static
void consdataUpdateActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed; can be NULL for global bound changes */
   SCIP_Real             oldbound,           /**< old bound of variable */
   SCIP_Real             newbound,           /**< new bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_BOUNDTYPE        boundtype,          /**< type of the bound change */
   SCIP_Bool             global,             /**< is it a global or a local bound change? */
   SCIP_Bool             pseudo,             /**< is this the call to update the pseudo activity? */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   SCIP_Real* activity;
   SCIP_Real* lastactivity;
   int* activityposinf;
   int* activityneginf;
   int* activityposhuge;
   int* activityneghuge;
   SCIP_Real delta;
   SCIP_Bool validact;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(global || pseudo || (var != NULL));
   assert(!global || !pseudo);
   assert(consdata->validactivities);
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->lastpseudoactivity < SCIP_INVALID);
   assert(consdata->pseudoactivityneginf >= 0);
   assert(consdata->pseudoactivityposinf >= 0);
   assert(consdata->pseudoactivityneghuge >= 0);
   assert(consdata->pseudoactivityposhuge >= 0);
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

   /* adjust pseudo activity, if we are updating a local activity and the variables best bound w.r.t. the objective was changed */
   if( !global && !pseudo && (SCIPvarGetBestBoundType(var) == boundtype) )
      consdataUpdateActivities(scip, consdata, var, oldbound, newbound, val, boundtype, FALSE, TRUE, checkreliability);

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
         if( val > 0.0 )
         {
            activity = &(consdata->glbminactivity);
            lastactivity = &(consdata->lastglbminactivity);
            activityposinf = &(consdata->glbminactivityposinf);
            activityneginf = &(consdata->glbminactivityneginf);
            activityposhuge = &(consdata->glbminactivityposhuge);
            activityneghuge = &(consdata->glbminactivityneghuge);
            validact = consdata->validglbminact;
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
         }
      }
      else
      {
         if( val > 0.0 )
         {
            activity = &(consdata->glbmaxactivity);
            lastactivity = &(consdata->lastglbmaxactivity);
            activityposinf = &(consdata->glbmaxactivityposinf);
            activityneginf = &(consdata->glbmaxactivityneginf);
            activityposhuge = &(consdata->glbmaxactivityposhuge);
            activityneghuge = &(consdata->glbmaxactivityneghuge);
            validact = consdata->validglbmaxact;
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
         }
      }
   }
   /* we are updating local activities */
   else
   {
      /* adjust pseudo activity; we do not have to distinguish bound types or coefficients,
       * because we only call the pseudoactivity update, if the bound is the best bound w.r.t. the objective function
       * and the sign of the coefficient does not matter in this case
       */
      if( pseudo )
      {
         activity = &(consdata->pseudoactivity);
         lastactivity = &(consdata->lastpseudoactivity);
         activityposinf = &(consdata->pseudoactivityposinf);
         activityneginf = &(consdata->pseudoactivityneginf);
         activityposhuge = &(consdata->pseudoactivityposhuge);
         activityneghuge = &(consdata->pseudoactivityneghuge);
         validact = consdata->validpseudoact;
      }
      /* depending on the boundtype and the coefficient, we choose the activity to be updated:
       * lower bound + pos. coef: update minactivity
       * lower bound + neg. coef: update maxactivity, positive and negative infinity counters have to be switched
       * upper bound + pos. coef: update maxactivity
       * upper bound + neg. coef: update minactivity, positive and negative infinity counters have to be switched
       */
      else
      {
         if( boundtype == SCIP_BOUNDTYPE_LOWER )
         {
            if( val > 0.0 )
            {
               activity = &(consdata->minactivity);
               lastactivity = &(consdata->lastminactivity);
               activityposinf = &(consdata->minactivityposinf);
               activityneginf = &(consdata->minactivityneginf);
               activityposhuge = &(consdata->minactivityposhuge);
               activityneghuge = &(consdata->minactivityneghuge);
               validact = consdata->validminact;
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
            }
         }
         else
         {
            if( val > 0.0 )
            {
               activity = &(consdata->maxactivity);
               lastactivity = &(consdata->lastmaxactivity);
               activityposinf = &(consdata->maxactivityposinf);
               activityneginf = &(consdata->maxactivityneginf);
               activityposhuge = &(consdata->maxactivityposhuge);
               activityneghuge = &(consdata->maxactivityneghuge);
               validact = consdata->validmaxact;
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
            }
         }
      }
   }

   /* old bound was +infinity */
   if( SCIPisInfinity(scip, oldbound) )
   {
      assert((*activityposinf) >= 1);

      /* we only have to do something ig the new bound is not again +infinity */
      if( !SCIPisInfinity(scip, newbound) )
      {
         /* decrease the counter for positive infinite contributions */
         (*activityposinf)--;

         /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
         if( SCIPisInfinity(scip, -newbound) )
            (*activityneginf)++;
         /* if the contribution of this variable is too large, increase the counter for huge values */
         else if( SCIPisHugeValue(scip, val * newbound) )
         {
            (*activityposhuge)++;
         }
         else if( SCIPisHugeValue(scip, -val * newbound) )
         {
            (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = val * newbound;
      }
   }
   /* old bound was -infinity */
   else if( SCIPisInfinity(scip, -oldbound) )
   {
      assert((*activityneginf) >= 1);

      /* we only have to do something ig the new bound is not again -infinity */
      if( !SCIPisInfinity(scip, -newbound) )
      {
         /* decrease the counter for negative infinite contributions */
         (*activityneginf)--;

         /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
         if( SCIPisInfinity(scip, newbound) )
            (*activityposinf)++;
         /* if the contribution of this variable is too large, increase the counter for huge values */
         else if( SCIPisHugeValue(scip, val * newbound) )
         {
            (*activityposhuge)++;
         }
         else if( SCIPisHugeValue(scip, -val * newbound) )
         {
            (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = val * newbound;
      }
   }
   /* old contribution was too large and positive */
   else if( SCIPisHugeValue(scip, val * oldbound) )
   {
      assert((*activityposhuge) >= 1);

      /* decrease the counter for huge positive contributions; it might be increased again later,
       * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
       */
      (*activityposhuge)--;

      /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
      if( SCIPisInfinity(scip, newbound) )
         (*activityposinf)++;
      /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
      else if( SCIPisInfinity(scip, -newbound) )
         (*activityneginf)++;
      /* if the contribution of this variable is too large and positive, increase the corresponding counter */
      else if( SCIPisHugeValue(scip, val * newbound) )
      {
         (*activityposhuge)++;
      }
      /* if the contribution of this variable is too large and negative, increase the corresponding counter */
      else if( SCIPisHugeValue(scip, -val * newbound) )
      {
         (*activityneghuge)++;
      }
      /* "normal case": just add the contribution to the activity */
      else
         delta = val * newbound;
   }
   /* old contribution was too large and negative */
   else if( SCIPisHugeValue(scip, -val * oldbound) )
   {
      assert((*activityneghuge) >= 1);

      /* decrease the counter for huge negative contributions; it might be increased again later,
       * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
       */
      (*activityneghuge)--;

      /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
      if( SCIPisInfinity(scip, newbound) )
         (*activityposinf)++;
      /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
      else if( SCIPisInfinity(scip, -newbound) )
         (*activityneginf)++;
      /* if the contribution of this variable is too large and positive, increase the corresponding counter */
      else if( SCIPisHugeValue(scip, val * newbound) )
      {
         (*activityposhuge)++;
      }
      /* if the contribution of this variable is too large and negative, increase the corresponding counter */
      else if( SCIPisHugeValue(scip, -val * newbound) )
      {
         (*activityneghuge)++;
      }
      /* "normal case": just add the contribution to the activity */
      else
         delta = val * newbound;
   }
   /* old bound was finite and not too large */
   else
   {
      /* if the new bound is +infinity, the old contribution has to be subtracted
       * and the counter for positive infinite contributions has to be increased
       */
      if( SCIPisInfinity(scip, newbound) )
      {
         (*activityposinf)++;
         delta = -val * oldbound;
      }
      /* if the new bound is -infinity, the old contribution has to be subtracted
       * and the counter for negative infinite contributions has to be increased
       */
      else if( SCIPisInfinity(scip, -newbound) )
      {
         (*activityneginf)++;
         delta = -val * oldbound;
      }
      /* if the contribution of this variable is too large, increase the counter for huge values */
      else if( SCIPisHugeValue(scip, val * newbound) )
      {
         (*activityposhuge)++;
         delta = -val * oldbound;
      }
      else if( SCIPisHugeValue(scip, -val * newbound) )
      {
         (*activityneghuge)++;
         delta = -val * oldbound;
      }
      /* "normal case": just update the activity */
      else
         delta = val * (newbound - oldbound);
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
      else
      {
         if( checkreliability && SCIPisUpdateUnreliable(scip, (*activity), (*lastactivity)) )
         {
            SCIPdebugMessage("%s%s activity of linear constraint unreliable after update: %16.9g\n", (global ? "global " : ""),
               (pseudo ? "pseudo" : ((boundtype == SCIP_BOUNDTYPE_LOWER) == (val > 0.0) ? "min" : "max")), (*activity));

            /* mark the activity that was just changed and is not reliable anymore to be invalid */
            if( global )
            {
               if( (boundtype == SCIP_BOUNDTYPE_LOWER) == (val > 0.0) )
                  consdata->validglbminact = FALSE;
               else
                  consdata->validglbmaxact = FALSE;
            }
            else
            {
               if( pseudo )
                  consdata->validpseudoact = FALSE;
               else
               {
                  if( (boundtype == SCIP_BOUNDTYPE_LOWER) == (val > 0.0) )
                     consdata->validminact = FALSE;
                  else
                     consdata->validmaxact = FALSE;
               }
            }
         }
      }
   }
}

/** updates minimum and maximum activity for a change in lower bound */
static
void consdataUpdateActivitiesLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, FALSE, FALSE, checkreliability);

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
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, FALSE, FALSE, checkreliability);

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
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, TRUE, FALSE, checkreliability);

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
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, TRUE, FALSE, checkreliability);

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
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
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

      absval = REALABS(val);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }

   /* update pseudo, minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, 0.0, SCIPvarGetLbLocal(var), val, checkreliability);
      consdataUpdateActivitiesUb(scip, consdata, var, 0.0, SCIPvarGetUbLocal(var), val, checkreliability);
      consdataUpdateActivitiesGlbLb(scip, consdata, 0.0, SCIPvarGetLbGlobal(var), val, checkreliability);
      consdataUpdateActivitiesGlbUb(scip, consdata, 0.0, SCIPvarGetUbGlobal(var), val, checkreliability);
   }
}

/** updates minimum and maximum activity for coefficient deletion, invalidates maximum absolute value if necessary */
static
void consdataUpdateDelCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = REALABS(val);

      if( SCIPisEQ(scip, absval, consdata->maxabsval) )
      {
         consdata->validmaxabsval = FALSE;
         consdata->maxabsval = SCIP_INVALID;
      }
   }

   /* update pseudo, minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, SCIPvarGetLbLocal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesUb(scip, consdata, var, SCIPvarGetUbLocal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesGlbLb(scip, consdata, SCIPvarGetLbGlobal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesGlbUb(scip, consdata, SCIPvarGetUbGlobal(var), 0.0, val, checkreliability);
   }
}

/** updates minimum and maximum activity for coefficient change, invalidates maximum absolute value if necessary */
static
void consdataUpdateChgCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             oldval,             /**< old coefficient of constraint entry */
   SCIP_Real             newval,             /**< new coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* old zero coefficients should be handled by consdataUpdateAddCoef() */
   assert(!SCIPisZero(scip, oldval));

   /* new zero coefficients should be handled by consdataUpdateDelCoef() */
   assert(!SCIPisZero(scip, newval));

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = REALABS(newval);

      if( SCIPisGE(scip, absval, consdata->maxabsval) )
         consdata->maxabsval = absval;
      else
      {
         absval = REALABS(oldval);

         /* invalidate maximum absolute value */
         if( SCIPisEQ(scip, absval, consdata->maxabsval) )
         {
            consdata->validmaxabsval = FALSE;
            consdata->maxabsval = SCIP_INVALID;
         }
      }
   }

   /* @todo do something more clever here, e.g. if oldval * newval >= 0, do the update directly */
   consdataUpdateDelCoef(scip, consdata, var, oldval, checkreliability);
   consdataUpdateAddCoef(scip, consdata, var, newval, checkreliability);
}

/** calculates maximum absolute value of coefficients */
static
void consdataCalcMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real absval;
   int i;

   assert(consdata != NULL);
   assert(!consdata->validmaxabsval);
   assert(consdata->maxabsval >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->maxabsval = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      absval = consdata->vals[i];
      absval = REALABS(absval);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }
}

/** returns the maximum absolute value of all coefficients in the constraint */
static
SCIP_Real consdataGetMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validmaxabsval )
      consdataCalcMaxAbsval(consdata);
   assert(consdata->validmaxabsval);
   assert(consdata->maxabsval < SCIP_INVALID);

   return consdata->maxabsval;
}

/** calculates pseudo activity, and minimum and maximum local and global activity for constraint;
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
   assert(consdata->pseudoactivity >= SCIP_INVALID);
   assert(consdata->minactivity >= SCIP_INVALID);
   assert(consdata->maxactivity >= SCIP_INVALID);
   assert(consdata->glbminactivity >= SCIP_INVALID);
   assert(consdata->glbmaxactivity >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->validactivities = TRUE;
   consdata->validpseudoact = TRUE;
   consdata->validminact = TRUE;
   consdata->validmaxact = TRUE;
   consdata->validglbminact = TRUE;
   consdata->validglbmaxact = TRUE;
   consdata->maxabsval = 0.0;
   consdata->pseudoactivity = 0.0;
   consdata->lastpseudoactivity = 0.0;
   consdata->pseudoactivityneginf = 0;
   consdata->pseudoactivityposinf = 0;
   consdata->pseudoactivityneghuge = 0;
   consdata->pseudoactivityposhuge = 0;
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
      consdataUpdateAddCoef(scip, consdata, consdata->vars[i], consdata->vals[i], FALSE);

   consdata->lastpseudoactivity = consdata->pseudoactivity;
   consdata->lastminactivity = consdata->minactivity;
   consdata->lastmaxactivity = consdata->maxactivity;
   consdata->lastglbminactivity = consdata->glbminactivity;
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
}

/** gets pseudo activity for constraint */
static
SCIP_Real consdataGetPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validpseudoact);
   }
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);

   SCIPdebugMessage("pseudo activity of linear constraint: %.15g\n", consdata->pseudoactivity);

   if( consdata->pseudoactivityneginf > 0 )
      return -SCIPinfinity(scip);
   else if( consdata->pseudoactivityposinf > 0 )
      return SCIPinfinity(scip);

   /* recompute the pseudo activity, if it is not valid */
   if( !consdata->validpseudoact )
      consdataRecomputePseudoactivity(scip, consdata);
   assert(consdata->validpseudoact);

   return consdata->pseudoactivity;
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
   /* we do not need a good relaxation and we have positve huge contributions, so we just return -infty as activity */
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
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Real*            maxactivity,        /**< pointer to store the maximal activity */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned minactivity is just a relaxation,
                                              *   i.e. <= the exact minactivity (in case of huge contributions),
                                              *   or equal to the exact minimal activity */
   SCIP_Bool*            maxisrelax          /**< pointer to store whether the returned maxactivity is just a relaxation,
                                              *   i.e. >= the exact maxactivity (in case of huge contributions),
                                              *   or equal to the exact maximal activity */
   )
{
   SCIP_Bool issettoinfinity;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validminact);
      assert(consdata->validmaxact);
   }
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);

   getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
      consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
      minactivity, minisrelax, &issettoinfinity);

   getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
      consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
      maxactivity, maxisrelax, &issettoinfinity);
}

/** calculates activity bounds for constraint after setting variable to zero */
static
void consdataGetReliableResidualActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             cancelvar,          /**< variable to calculate activity residual for */
   SCIP_Real*            resactivity,        /**< pointer to store the residual activity */
   SCIP_Bool             isminresact,        /**< should minimal or maximal residual activity be calculated? */
   SCIP_Bool             useglobalbounds     /**< should global or local bounds be used? */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(cancelvar != NULL);
   assert(resactivity != NULL);

   *resactivity = 0.0;

   for( v = 0; v < consdata->nvars; ++v )
   {
      var = consdata->vars[v];
      assert(var != NULL);
      if( var == cancelvar )
         continue;

      val = consdata->vals[v];

      if( useglobalbounds )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
      }
      else
      {
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }

      assert(!SCIPisZero(scip, val));
      assert(SCIPisLE(scip, lb, ub));

      if( val > 0.0 )
      {
         if( isminresact )
         {
            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisHugeValue(scip, REALABS(val*lb)));
            *resactivity += val*lb;
         }
         else
         {
            assert(!SCIPisInfinity(scip, ub));
            assert(!SCIPisHugeValue(scip, REALABS(val*ub)));
            *resactivity += val*ub;
         }
      }
      else
      {
         if( isminresact)
         {
            assert(!SCIPisInfinity(scip, ub));
            assert(!SCIPisHugeValue(scip, REALABS(val*ub)));
            *resactivity += val*ub;
         }
         else
         {
            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisHugeValue(scip, REALABS(val*lb)));
            *resactivity += val*lb;
         }
      }
   }
   assert(!SCIPisInfinity(scip, *resactivity) && !SCIPisInfinity(scip, -(*resactivity)));
}

/** gets activity bounds for constraint after setting variable to zero */
static
void consdataGetActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
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
   assert(consdata->pseudoactivity < SCIP_INVALID);
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

   if( val > 0.0 )
   {
      minactbound = SCIPvarGetLbLocal(var);
      maxactbound = SCIPvarGetUbLocal(var);
      absval = val;
   }
   else
   {
      minactbound = -SCIPvarGetUbLocal(var);
      maxactbound = -SCIPvarGetLbLocal(var);
      absval = -val;
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
      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge, absval * minactbound, FALSE, goodrelax,
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
      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, absval * maxactbound, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
}

/** gets global activity bounds for constraint */
static
void consdataGetGlbActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            glbminactivity,     /**< pointer to store the minimal activity, or NULL, if not needed */
   SCIP_Real*            glbmaxactivity,     /**< pointer to store the maximal activity, or NULL, if not needed */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned minactivity is just a relaxation,
                                              *   i.e. <= the exact minactivity (in case of huge contributions),
                                              *   or equal to the exact minimal activity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned maxactivity is just a relaxation,
                                              *   i.e. >= the exact maxactivity (in case of huge contributions),
                                              *   or equal to the exact maximal activity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minresactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxresactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert((glbminactivity != NULL && minisrelax != NULL && isminsettoinfinity != NULL)
      || (glbmaxactivity != NULL && maxisrelax != NULL && ismaxsettoinfinity != NULL));

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validglbminact);
      assert(consdata->validglbmaxact);
   }
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);
   assert(consdata->glbminactivityneghuge >= 0);
   assert(consdata->glbminactivityposhuge >= 0);
   assert(consdata->glbmaxactivityneghuge >= 0);
   assert(consdata->glbmaxactivityposhuge >= 0);

   if( glbminactivity != NULL )
   {
      assert(isminsettoinfinity != NULL);
      assert(minisrelax != NULL);

      getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
         consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
         glbminactivity, minisrelax, isminsettoinfinity);
   }

   if( glbmaxactivity != NULL )
   {
      assert(ismaxsettoinfinity != NULL);
      assert(maxisrelax != NULL);

      getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
         consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
         glbmaxactivity, maxisrelax, ismaxsettoinfinity);
   }
}

/** gets global activity bounds for constraint after setting variable to zero */
static
void consdataGetGlbActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity, or NULL, if not needed */
   SCIP_Real*            maxresactivity,     /**< pointer to store the maximal residual activity, or NULL, if not needed */
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

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert((minresactivity != NULL && minisrelax != NULL && isminsettoinfinity != NULL )
      || (maxresactivity != NULL && maxisrelax != NULL && ismaxsettoinfinity != NULL));

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);

   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);

   if( val > 0.0 )
   {
      minactbound = SCIPvarGetLbGlobal(var);
      maxactbound = SCIPvarGetUbGlobal(var);
      absval = val;
   }
   else
   {
      minactbound = -SCIPvarGetUbGlobal(var);
      maxactbound = -SCIPvarGetLbGlobal(var);
      absval = -val;
   }

   if( minresactivity != NULL )
   {
      assert(isminsettoinfinity != NULL);
      assert(minisrelax != NULL);

      /* get/compute minactivity by calling getMinActivity() with updated counters for infinite and huge values
       * and contribution of variable set to zero that has to be subtracted from finite part of activity
       */
      if( SCIPisInfinity(scip, minactbound) )
      {
         assert(consdata->glbminactivityposinf >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf - 1, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisInfinity(scip, -minactbound) )
      {
         assert(consdata->glbminactivityneginf >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf - 1,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, minactbound * absval) )
      {
         assert(consdata->glbminactivityposhuge >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge - 1, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, -minactbound * absval) )
      {
         assert(consdata->glbminactivityneghuge >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge - 1, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else
      {
         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, absval * minactbound, TRUE,
            goodrelax, minresactivity, minisrelax, isminsettoinfinity);
      }
   }

   if( maxresactivity != NULL )
   {
      assert(ismaxsettoinfinity != NULL);
      assert(maxisrelax != NULL);

      /* get/compute maxactivity by calling getMaxActivity() with updated counters for infinite and huge values
       * and contribution of variable set to zero that has to be subtracted from finite part of activity
       */
      if( SCIPisInfinity(scip, -maxactbound) )
      {
         assert(consdata->glbmaxactivityneginf >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf - 1,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisInfinity(scip, maxactbound) )
      {
         assert(consdata->glbmaxactivityposinf >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf - 1, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, absval * maxactbound) )
      {
         assert(consdata->glbmaxactivityposhuge >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge - 1, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, -absval * maxactbound) )
      {
         assert(consdata->glbmaxactivityneghuge >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge - 1, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else
      {
         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, absval * maxactbound, TRUE,
            goodrelax, maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
   }
}

/** calculates the activity of the linear constraint for given solution */
static
SCIP_Real consdataGetActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get activity for, NULL to current solution */
   )
{
   SCIP_Real activity;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
   {
      /* for performance reasons, the pseudo activity is updated with each bound change, so we don't have to
       * recalculate it
       */
      activity = consdataGetPseudoActivity(scip, consdata);
   }
   else
   {
      SCIP_Real solval;
      int nposinf;
      int nneginf;
      SCIP_Bool negsign;
      int v;

      activity = 0.0;
      nposinf = 0;
      nneginf = 0;
      negsign = 0;

      for( v = 0; v < consdata->nvars; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);

         if( consdata->vals[v] < 0 )
            negsign = TRUE;
         else 
            negsign = FALSE;

         if( (SCIPisInfinity(scip, solval) && !negsign) || (SCIPisInfinity(scip, -solval) && negsign) )
            ++nposinf;
         else if( (SCIPisInfinity(scip, solval) && negsign) || (SCIPisInfinity(scip, -solval) && !negsign) )
            ++nneginf;
         else
            activity += consdata->vals[v] * solval;
      }
      assert(nneginf >= 0 && nposinf >= 0);

      SCIPdebugMessage("activity of linear constraint: %.15g, %d positive infinity values, %d negative infinity values \n", activity, nposinf, nneginf);

      /* check for amount of infinity values and correct the activity */
      if( nposinf > 0 && nneginf > 0 )
         activity = (consdata->rhs + consdata->lhs) / 2;
      else if( nposinf > 0 )
         activity = SCIPinfinity(scip);
      else if( nneginf > 0 )
         activity = -SCIPinfinity(scip);

      SCIPdebugMessage("corrected activity of linear constraint: %.15g\n", activity);
   }

   if( activity < 0 )
      activity = MAX(activity, -SCIPinfinity(scip)); /*lint !e666*/
   else
      activity = MIN(activity, SCIPinfinity(scip)); /*lint !e666*/

   return activity;
}

/** calculates the feasibility of the linear constraint for given solution */
static
SCIP_Real consdataGetFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get feasibility for, NULL to current solution */
   )
{
   SCIP_Real activity;

   assert(scip != NULL);
   assert(consdata != NULL);

   activity = consdataGetActivity(scip, consdata, sol);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}

/** returns the signature bitmask for the given variable */
static
SCIP_Longint getVarSignature(
   SCIP_VAR*             var                 /**< variable */
   )
{
   int sigidx;

   assert(var != NULL);

   sigidx = SCIPvarGetIndex(var) % (int)(8*sizeof(SCIP_Longint));
   return ((SCIP_Longint)1) << sigidx; /*lint !e703*/
}

/** updates bit signatures after adding a single coefficient */
static
void consdataUpdateSignatures(
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   pos                 /**< position of coefficient to update signatures for */
   )
{
   SCIP_Longint varsignature;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real val;

   assert(consdata != NULL);
   assert(consdata->validsignature);

   varsignature = getVarSignature(consdata->vars[pos]);
   lb = SCIPvarGetLbGlobal(consdata->vars[pos]);
   ub = SCIPvarGetUbGlobal(consdata->vars[pos]);
   val = consdata->vals[pos];
   if( (val > 0.0 && ub > 0.0) || (val < 0.0 && lb < 0.0) )
      consdata->possignature |= varsignature;
   if( (val > 0.0 && lb < 0.0) || (val < 0.0 && ub > 0.0) )
      consdata->negsignature |= varsignature;
}

/** calculates the bit signatures of the given constraint data */
static
void consdataCalcSignatures(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validsignature )
   {
      int i;

      consdata->validsignature = TRUE;
      consdata->possignature = 0;
      consdata->negsignature = 0;
      for( i = 0; i < consdata->nvars; ++i )
         consdataUpdateSignatures(consdata, i);
   }
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

/** permutes the constraint's variables according to a given permutation. */
static 
void permSortConsdata(
   SCIP_CONSDATA*        consdata,           /**< the constraint data */
   int*                  perm,               /**< the target permutation */
   int                   nvars,              /**< the number of variables */
   SCIP_Bool             isinpresolving      /**< is the scip stage before initsolve */
   )
{  /*lint --e{715}*/
   SCIP_VAR* varv;
   SCIP_EVENTDATA* eventdatav;
   SCIP_Real valv;
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
         if( consdata->eventdatas != NULL )
            eventdatav = consdata->eventdatas[v];
         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < nvars);
            assert(perm[i] != i);
            consdata->vars[i] = consdata->vars[perm[i]];
            consdata->vals[i] = consdata->vals[perm[i]];
            if( consdata->eventdatas != NULL )
            {
               consdata->eventdatas[i] = consdata->eventdatas[perm[i]];
               consdata->eventdatas[i]->varpos = i;
            }
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->vars[i] = varv;
         consdata->vals[i] = valv;
         if( consdata->eventdatas != NULL )
         {
            consdata->eventdatas[i] = eventdatav;
            consdata->eventdatas[i]->varpos = i;
         }
         perm[i] = i;
      }
   }
#ifdef SCIP_DEBUG
   /* check sorting */
   for( v = 0; v < nvars; ++v )
   {
      if( isinpresolving )
      {
         assert(v == nvars-1 || SCIPvarCompare(consdata->vars[v], consdata->vars[v+1]) <= 0);
      }
      assert(perm[v] == v);
      assert(consdata->eventdatas == NULL || consdata->eventdatas[v]->varpos == v);
   }
#endif
}

/** sorts linear constraint's variables depending on the stage of the solving process:
 *  - during PRESOLVING
 *       sorts variables by binaries, integers, implicit integers, and continuous variables,
 *       and the variables of the same type by non-decreasing variable index
 *
 * -  during SOLVING
 *       sorts binary variables of the remaining problem w.r.t the absolute of their coefficient.
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
      consdata->sorted = TRUE;
      consdata->binvarssorted = TRUE;
   }
   else if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE && !consdata->sorted )
   {
      int* perm;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call sorting method  */
      SCIPsort(perm, consdataCompVar, (void*)consdata, consdata->nvars);

      permSortConsdata(consdata, perm, consdata->nvars, TRUE);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);

      consdata->sorted = TRUE;
      consdata->binvarssorted = FALSE;
   }
   else if( SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->binvarssorted )
   {
      SCIP_EVENTDATA** eventdatas;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;
      int v;
      int lastbin;

      nvars = consdata->nvars;
      vars = consdata->vars;
      vals = consdata->vals;
      eventdatas = consdata->eventdatas;
      assert(vars != NULL || nvars == 0);
      assert(vals != NULL || nvars == 0);

      lastbin = 0;
      /* count binary variables and permute variables such that binaries appear first in the sorted vars array */
      for( v = 0; v < nvars; ++v )
      {
         assert( vars != NULL); /* for flexelint */
         assert( vals != NULL); /* for flexelint */
         if( SCIPvarIsBinary(vars[v]) )
         {
            /* swap variable at the end of the binary variables, if necessary */
            if( lastbin < v )
            {
               SCIP_VAR* tmpvar;
               SCIP_Real tmpval;
               
               tmpvar = vars[lastbin];
               tmpval = vals[lastbin];

               vars[lastbin] = vars[v];
               vals[lastbin] = vals[v];

               vars[v] = tmpvar;
               vals[v] = tmpval;

               if( eventdatas != NULL )
               {
                  SCIP_EVENTDATA* tmpeventdata;
                  
                  tmpeventdata = eventdatas[lastbin];
                  eventdatas[lastbin] = eventdatas[v];
                  eventdatas[lastbin]->varpos = lastbin;
                  eventdatas[v] = tmpeventdata;
                  eventdatas[v]->varpos = v;
               }
               assert(SCIPvarIsBinary(vars[lastbin]));
            }
#ifndef NDEBUG
            else
               assert(lastbin == v);
#endif
            ++lastbin;
         }
      }
      consdata->nbinvars = lastbin;

#ifndef NDEBUG
      /* check sorting */
      for( v = 0; v < nvars; ++v )
      {
         assert(vars != NULL); /* for flexelint */
         assert(eventdatas == NULL || eventdatas[v]->varpos == v);
         assert((v >= consdata->nbinvars && !SCIPvarIsBinary(vars[v])) || (v < consdata->nbinvars && SCIPvarIsBinary(vars[v])));
      }
#endif

      if( consdata->nbinvars > 1 )
      {
         SCIP_Real* absvals;
         int*       perm;

         assert(lastbin == consdata->nbinvars);
         assert(lastbin <= nvars);
         assert(vals != NULL);

         /* initialize absolute coefficients and the target permutation for binary variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &absvals, lastbin) );
         SCIP_CALL( SCIPallocBufferArray(scip, &perm, lastbin) );

         for( v = 0; v < lastbin; ++v )
         {
            absvals[v] = ABS(vals[v]);
            perm[v] = v;
         }

         /* execute the sorting */
         SCIPsortDownRealInt(absvals, perm, lastbin);

         permSortConsdata(consdata, perm, lastbin, FALSE);

         /* free temporary arrays */
         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &absvals);
      }
      consdata->binvarssorted = TRUE;
         
      /* presolve sorting cannot be guaranteed after binary sorting */
      consdata->sorted = (consdata->sorted && consdata->nbinvars == 0);
   }
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE || consdata->binvarssorted);
   assert(SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE || consdata->sorted);

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
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if ( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, lhs, consdata->rhs) )
   {
      consdata->rhs = lhs;
      assert(consdata->row == NULL);
   }

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from a non-infinite value to -infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
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

   if( !SCIPisInfinity(scip, -lhs) && SCIPisGT(scip, lhs, consdata->lhs) )
   {
      consdata->boundstightened = FALSE;
      consdata->propagated = FALSE;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
   }

   /* new left hand side */
   consdata->lhs = lhs;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;

   /* update the lhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, consdata->row, lhs) );
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if ( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, rhs, consdata->lhs) )
   {
      consdata->lhs = rhs;
      assert(consdata->row == NULL);
   }

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from a non-infinite value to infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
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

   if( !SCIPisInfinity(scip, rhs) && SCIPisLT(scip, rhs, consdata->rhs) )
   {
      consdata->boundstightened = FALSE;
      consdata->propagated = FALSE;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
   }

   /* set new right hand side */
   consdata->rhs = rhs;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;

   /* update the rhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, consdata->row, rhs) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, val) )
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
   consdata->vals[consdata->nvars] = val;
   consdata->nvars++;
   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, the variable needs an additional event data */
   if( transformed )
   {
      if( consdata->eventdatas != NULL )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONSHDLRDATA* conshdlrdata;

         /* get event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* initialize eventdatas array */
         consdata->eventdatas[consdata->nvars-1] = NULL;

         /* catch bound change events of variable */
         SCIP_CALL( consdataCatchEvent(scip, consdata, conshdlrdata->eventhdlr, consdata->nvars-1) );
      }

      /* update minimum and maximum activities */
      consdataUpdateAddCoef(scip, consdata, var, val, FALSE); 
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockRounding(scip, cons, var, val) );

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->removedfixings = consdata->removedfixings && SCIPvarIsActive(var);
   if( consdata->validsignature )
      consdataUpdateSignatures(consdata, consdata->nvars-1);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   if( consdata->nvars == 1 )
   {
      consdata->binvarssorted = TRUE;
      consdata->sorted = TRUE;
      consdata->merged = TRUE;
   }
   else
   {
      consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(var);
      consdata->sorted = consdata->sorted
         && (SCIPvarCompare(consdata->vars[consdata->nvars-2], consdata->vars[consdata->nvars-1]) == -1);
      consdata->merged = FALSE;
   }

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, val) );
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
   SCIP_Real val;

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

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      if( consdata->eventdatas != NULL )
      {
         SCIP_CALL( consdataDropEvent(scip, consdata, conshdlrdata->eventhdlr, pos) );
         assert(consdata->eventdatas[pos] == NULL);
      }
   }

   /* move the last variable to the free slot */
   if( pos != consdata->nvars-1 )
   {
      consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(consdata->vars[pos]);

      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->vals[pos] = consdata->vals[consdata->nvars-1];

      if( consdata->eventdatas != NULL )
      {
         consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
         assert(consdata->eventdatas[pos] != NULL);
         consdata->eventdatas[pos]->varpos = pos;
      }
      consdata->sorted = FALSE;
   }
   consdata->nvars--;

   /* if at most one variable is left, the activities should be recalculated (to correspond exactly to the bounds
    * of the remaining variable, or give exactly 0.0)
    */
   if( consdata->nvars <= 1 )
      consdataInvalidateActivities(consdata);
   else
   {
      /* if we are in transformed problem, update minimum and maximum activities */
      if( SCIPconsIsTransformed(cons) )
         consdataUpdateDelCoef(scip, consdata, var, val, TRUE);
   }

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->validsignature = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   /* delete coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, -val) );
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
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newval));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(!SCIPisZero(scip, newval));

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* if necessary, update the rounding locks of the variable */
   if( SCIPconsIsLocked(cons) && newval * val < 0.0 )
   {
      assert(SCIPconsIsTransformed(cons));

      /* remove rounding locks for variable with old coefficient */
      SCIP_CALL( unlockRounding(scip, cons, var, val) );

      /* install rounding locks for variable with new coefficient */
      SCIP_CALL( lockRounding(scip, cons, var, newval) );
   }

   /* change the value */
   consdata->vals[pos] = newval;

   consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(var);

   /* update minimum and maximum activities */
   if( SCIPconsIsTransformed(cons) )
      consdataUpdateChgCoef(scip, consdata, var, val, newval, TRUE);

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->validsignature = consdata->validsignature && (newval * val > 0.0);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   return SCIP_OKAY;
}

/** scales a linear constraint with a constant scalar */
static
SCIP_RETCODE scaleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint to scale */
   SCIP_Real             scalar              /**< value to scale constraint with */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newval;
   SCIP_Real absscalar;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);
   assert(!SCIPisEQ(scip, scalar, 1.0));

   /* scale the coefficients */
   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      newval = scalar * consdata->vals[i];

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasFloor, we add an additional 0.5 before
       * flooring down our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->vals[i], scalar) )
         newval = SCIPfeasFloor(scip, newval + 0.5);

      if( SCIPisZero(scip, newval) )
      {
         SCIPwarningMessage(scip, "coefficient %.15g of variable <%s> in linear constraint <%s> scaled to zero (scalar: %.15g)\n",
            consdata->vals[i], SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons), scalar);
         SCIP_CALL( delCoefPos(scip, cons, i) );
      }
      else
         consdata->vals[i] = newval;
   }

   /* scale the sides */
   if( scalar < 0.0 )
   {
      SCIP_Real lhs;

      lhs = consdata->lhs;
      consdata->lhs = -consdata->rhs;
      consdata->rhs = -lhs;
   }
   absscalar = REALABS(scalar);
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      newval = absscalar * consdata->lhs;

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasFloor, we add an additional 0.5 before
       * flooring down our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->lhs, absscalar) )
         consdata->lhs = SCIPfeasFloor(scip, newval + 0.5);
      else
         consdata->lhs = newval;
   }
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      newval = absscalar * consdata->rhs;

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasCeil, we subtract 0.5 before ceiling up
       * our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->rhs, absscalar) )
         consdata->rhs = SCIPfeasCeil(scip, newval - 0.5);
      else
         consdata->rhs = newval;
   }

   consdataInvalidateActivities(consdata);
   consdata->cliquesadded = FALSE;

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


/** normalizes a linear constraint with the following rules:
 *  - if all coefficients have them same absolute value, change them to (-)1.0
 *  - multiplication with +1 or -1:
 *      Apply the following rules in the given order, until the sign of the factor is determined. Later rules only apply,
 *      if the current rule doesn't determine the sign):
 *        1. the right hand side must not be negative
 *        2. the right hand side must not be infinite
 *        3. the absolute value of the right hand side must be greater than that of the left hand side
 *        4. the number of positive coefficients must not be smaller than the number of negative coefficients
 *        5. multiply with +1
 *  - rationals to integrals
 *      Try to identify a rational representation of the fractional coefficients, and multiply all coefficients
 *      by the smallest common multiple of all denominators to get integral coefficients.
 *      Forbid large denominators due to numerical stability.
 *  - division by greatest common divisor
 *      If all coefficients are integral, divide them by the greatest common divisor.
 */
static
SCIP_RETCODE normalizeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint to normalize */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Longint gcd;
   SCIP_Longint maxmult;
   SCIP_Real epsilon;
   SCIP_Real feastol;
   SCIP_Real maxabsval;
   SCIP_Bool success;
   int nvars;
   int mult;
   int nposcoeffs;
   int nnegcoeffs;
   int i;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   /* we must not change a modifiable constraint in any way */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint is already normalized */
   if( consdata->normalized )
      return SCIP_OKAY;

   /* get coefficient arrays */
   vals = consdata->vals;
   nvars = consdata->nvars;
   vars = consdata->vars;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( nvars == 0 )
   {
      consdata->normalized = TRUE;
      return SCIP_OKAY;
   }

   assert(vars != NULL);
   assert(vals != NULL);

   /* get maximal absolute coefficient */
   maxabsval = consdataGetMaxAbsval(consdata);

   /* check if all coefficients are in absolute value equal, and not 1.0 */
   if( !SCIPisEQ(scip, maxabsval, 1.0) )
   {
      SCIP_Bool abscoefsequ;

      abscoefsequ = TRUE;

      for( v = nvars - 1; v >= 0; --v )
      {
	 if( !SCIPisEQ(scip, REALABS(vals[v]), maxabsval) )
	 {
	    abscoefsequ = FALSE;
	    break;
	 }
      }

      /* all coefficients are in absolute value equal, so change them to (-)1.0 */
      if( abscoefsequ )
      {
	 SCIPdebugMessage("divide linear constraint with %g, because all coefficents are in absolute value the same\n", maxabsval);
	 SCIPdebugPrintCons(scip, cons, NULL);
	 SCIP_CALL( scaleCons(scip, cons, 1/maxabsval) );

	 if( consdata->validmaxabsval )
	 {
	    if( !SCIPisEQ(scip, consdata->maxabsval, 1.0) )
	       consdata->maxabsval = 1.0;

	    maxabsval = 1.0;
	 }
	 else
	 {
	    /* get maximal absolute coefficient */
	    maxabsval = consdataGetMaxAbsval(consdata);
	 }

	 /* get new consdata information, because scalecons() might have deleted variables */
	 vals = consdata->vals;
	 nvars = consdata->nvars;
	 vars = consdata->vars;

	 assert(nvars == 0 || vars != NULL);
	 assert(nvars == 0 || vals != NULL);
      }
   }

   /* nvars might have changed */
   if( nvars == 0 )
   {
      consdata->normalized = TRUE;
      return SCIP_OKAY;
   }

   assert(vars != NULL);
   assert(vals != NULL);

   /* calculate the maximal multiplier for common divisor calculation:
    *   |p/q - val| < epsilon  and  q < feastol/epsilon  =>  |p - q*val| < feastol
    * which means, a value of feastol/epsilon should be used as maximal multiplier;
    * additionally, we don't want to scale the constraint if this would lead to too
    * large coefficients
    */
   epsilon = SCIPepsilon(scip) * 0.9;  /* slightly decrease epsilon to be safe in rational conversion below */
   feastol = SCIPfeastol(scip);
   maxmult = (SCIP_Longint)(feastol/epsilon + feastol);
   maxmult = MIN(maxmult, (SCIP_Longint)( MAXSCALEDCOEF/MAX(maxabsval, 1.0)));

   /* if all variables are of integral type we will allow a greater multiplier */
   if( consdata->sorted )
   {
      if( SCIPvarGetType(vars[nvars - 1]) != SCIP_VARTYPE_CONTINUOUS )
      {
	 maxmult = (SCIP_Longint) (MAXSCALEDCOEFINTEGER/(MAX(maxabsval, 1.0)));
      }
   }
   else
   {
      SCIP_Bool foundcont;

      foundcont = FALSE;

      for( v = nvars - 1; v >= 0; --v )
      {
	 if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
	 {
	    foundcont = TRUE;
	    break;
	 }
      }

      if( !foundcont )
      {
	 maxmult = (SCIP_Longint) (MAXSCALEDCOEFINTEGER/(MAX(maxabsval, 1.0)));
      }
   }

   /*
    * multiplication with +1 or -1
    */
   mult = 0;

   /* 1. the right hand side must not be negative */
   if( SCIPisPositive(scip, consdata->lhs) )
      mult = +1;
   else if( SCIPisNegative(scip, consdata->rhs) )
      mult = -1;

   if( mult == 0 )
   {
      /* 2. the right hand side must not be infinite */
      if( SCIPisInfinity(scip, -consdata->lhs) )
         mult = +1;
      else if( SCIPisInfinity(scip, consdata->rhs) )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 3. the absolute value of the right hand side must be greater than that of the left hand side */
      if( SCIPisGT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
         mult = +1;
      else if( SCIPisLT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 4. the number of positive coefficients must not be smaller than the number of negative coefficients */
      nposcoeffs = 0;
      nnegcoeffs = 0;
      for( i = 0; i < nvars; ++i )
      {
         if( vals[i] > 0.0 )
            nposcoeffs++;
         else
            nnegcoeffs++;
      }
      if( nposcoeffs > nnegcoeffs )
         mult = +1;
      else if( nposcoeffs < nnegcoeffs )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 5. multiply with +1 */
      mult = +1;
   }

   assert(mult == +1 || mult == -1);
   if( mult == -1 )
   {
      /* scale the constraint with -1 */
      SCIPdebugMessage("multiply linear constraint with -1.0\n");
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( scaleCons(scip, cons, -1.0) );

      /* scalecons() can delete variables, but scaling with -1 should not do that */
      assert(nvars == consdata->nvars);
   }

   /*
    * rationals to integrals
    */
   success = TRUE;
   scm = 1;
   for( i = 0; i < nvars && success && scm <= maxmult; ++i )
   {
      if( !SCIPisIntegral(scip, vals[i]) )
      {
         /* epsilon has been slightly decreased above - to be on the safe side */
         success = SCIPrealToRational(vals[i], -epsilon, epsilon , maxmult, &nominator, &denominator);
         if( success )
            scm = SCIPcalcSmaComMul(scm, denominator);
      }
   }
   assert(scm >= 1);
   success = success && (scm <= maxmult);
   if( success && scm != 1 )
   {
      /* scale the constraint with the smallest common multiple of all denominators */
      SCIPdebugMessage("scale linear constraint with %"SCIP_LONGINT_FORMAT" to make coefficients integral\n", scm);
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( scaleCons(scip, cons, (SCIP_Real)scm) );

      if( consdata->validmaxabsval )
      {
	 consdata->maxabsval *= REALABS((SCIP_Real)scm);
      }

      /* get new consdata information, because scalecons() might have deleted variables */
      vals = consdata->vals;
      nvars = consdata->nvars;
      assert(nvars == 0 || vals != NULL);
   }

   /*
    * division by greatest common divisor
    */
   if( success && nvars >= 1 )
   {
      /* all coefficients are integral: divide them by their greatest common divisor */
      assert(SCIPisIntegral(scip, vals[0]));
      gcd = (SCIP_Longint)(REALABS(vals[0]) + feastol);
      assert(gcd >= 1);
      for( i = 1; i < nvars && gcd > 1; ++i )
      {
         assert(SCIPisIntegral(scip, vals[i]));
         gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[i]) + feastol));
      }

      if( gcd > 1 )
      {
         /* divide the constraint by the greatest common divisor of the coefficients */
         SCIPdebugMessage("divide linear constraint by greatest common divisor %"SCIP_LONGINT_FORMAT"\n", gcd);
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIP_CALL( scaleCons(scip, cons, 1.0/(SCIP_Real)gcd) );

	 if( consdata->validmaxabsval )
	 {
	    consdata->maxabsval /= REALABS((SCIP_Real)gcd);
	 }
      }
   }

   /* mark constraint to be normalized */
   consdata->normalized = TRUE;

   SCIPdebugMessage("normalized constraint:\n");
   SCIPdebugPrintCons(scip, cons, NULL);

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
   SCIP_Real valsum;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

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
         valsum = consdata->vals[v];
         do
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );
            --v;
            valsum += consdata->vals[v];
         }
         while( v >= 1 && consdata->vars[v-1] == var );

         /* modify the last existing occurrence of the variable */
         assert(consdata->vars[v] == var);
         if( SCIPisZero(scip, valsum) )
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, v, valsum) );
         }
      }
      --v;
   }

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
   SCIP_Real val;
   SCIP_Real* aggrscalars;
   SCIP_Real fixedval;
   SCIP_Real aggrconst;
   int v;
   int naggrvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !consdata->removedfixings )
   {
      SCIP_Real lhssubtrahend;
      SCIP_Real rhssubtrahend;

      lhssubtrahend = 0.0;
      rhssubtrahend = 0.0;

      SCIPdebugMessage("applying fixings:\n");
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
            fixedval = SCIPvarGetLbGlobal(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               if( SCIPisInfinity(scip, ABS(fixedval)) )
               {
                  if( val * fixedval > 0.0 )
                  {
                     SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
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
                        SCIP_CALL( chgLhs(scip, cons, SCIPinfinity(scip)) );
                     }
                  }
               }
               else
                  lhssubtrahend += val * fixedval;
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               if( SCIPisInfinity(scip, ABS(fixedval)) )
               {
                  if( val * fixedval > 0.0 )
                  {
                     if( infeasible != NULL )
                     {
                        /* if rhs gets -infinity it means that the problem is infeasible */
                        *infeasible = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        SCIP_CALL( chgRhs(scip, cons, -SCIPinfinity(scip)) );
                     }
                  }
                  else
                  {
                     SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
                  }
               }
               else
                  rhssubtrahend += val * fixedval;
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_AGGREGATED:
	 {
	    SCIP_VAR* activevar = SCIPvarGetAggrVar(var);
	    SCIP_Real activescalar = val * SCIPvarGetAggrScalar(var);
	    SCIP_Real activeconstant = val * SCIPvarGetAggrConstant(var);

	    assert(activevar != NULL);
	    SCIP_CALL( SCIPgetProbvarSum(scip, &activevar, &activescalar, &activeconstant) );
	    assert(activevar != NULL);

	    if( !SCIPisZero(scip, activescalar) )
	    {
	       SCIP_CALL( addCoef(scip, cons, activevar, activescalar) );
	    }

	    if( !SCIPisZero(scip, activeconstant) )
	    {
	       if( !SCIPisInfinity(scip, -consdata->lhs) )
		  lhssubtrahend += activeconstant;
	       if( !SCIPisInfinity(scip, consdata->rhs) )
		  rhssubtrahend += activeconstant;
	    }

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;
	 }
         case SCIP_VARSTATUS_MULTAGGR:
            SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );
            naggrvars = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            for( i = 0; i < naggrvars; ++i )
            {
               SCIP_CALL( addCoef(scip, cons, aggrvars[i], val * aggrscalars[i]) );
            }
            aggrconst = SCIPvarGetMultaggrConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_NEGATED:
            SCIP_CALL( addCoef(scip, cons, SCIPvarGetNegationVar(var), -val) );
            aggrconst = SCIPvarGetNegationConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         default:
            SCIPerrorMessage("unknown variable status\n");
            SCIPABORT();
         }
      }
      
      if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->lhs) )
      {
         /* for large numbers that are relatively equal, substraction can lead to cancellation,
          * causing wrong fixings of other variables --> better use a real zero here;
          * for small numbers, polishing the difference might lead to wrong results -->
          * better use the exact difference in this case
          */
         if( SCIPisFeasEQ(scip, lhssubtrahend, consdata->lhs) && SCIPisFeasGE(scip, REALABS(lhssubtrahend), 1.0) ) 
         {
            SCIP_CALL( chgLhs(scip, cons, 0.0) );
         }
         else
         {
            SCIP_CALL( chgLhs(scip, cons, consdata->lhs - lhssubtrahend) );
         }
      }
      if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs))
      {

         /* for large numbers that are relatively equal, substraction can lead to cancellation,
          * causing wrong fixings of other variables --> better use a real zero here;
          * for small numbers, polishing the difference might lead to wrong results -->
          * better use the exact difference in this case
          */
         if( SCIPisFeasEQ(scip, rhssubtrahend, consdata->rhs ) && SCIPisFeasGE(scip, REALABS(rhssubtrahend), 1.0) )
         {
            SCIP_CALL( chgRhs(scip, cons, 0.0) );
         }
         else
         {
            SCIP_CALL( chgRhs(scip, cons, consdata->rhs - rhssubtrahend) );
         }
      }

      consdata->removedfixings = TRUE;

      SCIPdebugMessage("after fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
       * to clean up the constraint
       */
      SCIP_CALL( mergeMultiples(scip, cons) );

      SCIPdebugMessage("after merging:\n");
      SCIPdebugPrintCons(scip, cons, NULL);
   }
   assert(consdata->removedfixings);

#ifndef NDEBUG
   /* check, if all fixings are applied */
   for( v = 0; v < consdata->nvars; ++v )
      assert(SCIPvarIsActive(consdata->vars[v]));
#endif

   return SCIP_OKAY;
}

/** for each variable in the linear constraint, except the inferred variable, adds one bound to the conflict analysis'
 *  candidate store (bound depends on sign of coefficient and whether the left or right hand side was the reason for the
 *  inference variable's bound change); the conflict analysis can be initialized with the linear constraint being the
 *  conflict detecting constraint by using NULL as inferred variable
 */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   int                   inferpos,           /**< position of the inferred variable in the vars array */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the bound change? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   vals = consdata->vals;
   nvars = consdata->nvars;
   assert(vars != NULL || nvars == 0);
   assert(vals != NULL || nvars == 0);
   assert(-1 <= inferpos && inferpos < nvars);
   assert((infervar == NULL) == (inferpos == -1));
   assert(inferpos == -1 || vars[inferpos] == infervar); /*lint !e613*/

   /* for each variable, add the bound to the conflict queue, that is responsible for the minimal or maximal
    * residual value, depending on whether the left or right hand side is responsible for the bound change:
    *  - if the right hand side is the reason, the minimal residual activity is responsible
    *  - if the left hand side is the reason, the maximal residual activity is responsible
    */

   /* if the variable is integral we only need to add reason bounds until the propagation could be applied */
   if( infervar == NULL || SCIPvarIsIntegral(infervar) )
   {
      SCIP_Real minresactivity;
      SCIP_Real maxresactivity;
      SCIP_Bool minisrelax;
      SCIP_Bool maxisrelax;
      SCIP_Bool isminsettoinfinity;
      SCIP_Bool ismaxsettoinfinity;

      minresactivity = -SCIPinfinity(scip);
      maxresactivity = SCIPinfinity(scip);

      /* calculate the minimal and maximal global activity of all other variables involved in the constraint */
      if( infervar != NULL )
      {
         assert(vals != NULL); /* for flexelint */
         if( reasonisrhs )
            consdataGetGlbActivityResiduals(scip, consdata, infervar, vals[inferpos], FALSE, &minresactivity, NULL,
               &minisrelax, NULL, &isminsettoinfinity, NULL);
         else
            consdataGetGlbActivityResiduals(scip, consdata, infervar, vals[inferpos], FALSE, NULL, &maxresactivity,
               NULL, &maxisrelax, NULL, &ismaxsettoinfinity);
      }
      else
      {
         if( reasonisrhs )
            consdataGetGlbActivityBounds(scip, consdata, FALSE, &minresactivity, NULL,
               &minisrelax, NULL, &isminsettoinfinity, NULL);
         else
            consdataGetGlbActivityBounds(scip, consdata, FALSE, NULL, &maxresactivity,
               NULL, &maxisrelax, NULL, &ismaxsettoinfinity);
      }

      /* we can only do something clever, if the residual activity is finite and not relaxed */
      if( (reasonisrhs && !isminsettoinfinity && !minisrelax) || (!reasonisrhs && !ismaxsettoinfinity && !maxisrelax) ) /*lint !e644*/
      {
         SCIP_Real rescap;
         SCIP_Bool resactisinf;

         resactisinf = FALSE;

         /* calculate the residual capacity that would be left, if the variable would be set to one more / one less
          * than its inferred bound
          */
         if( infervar != NULL )
         {
            assert(vals != NULL); /* for flexelint */

            if( reasonisrhs )
            {
               if( SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastglbminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, infervar, &minresactivity, TRUE, TRUE);
                  if( SCIPisInfinity(scip, -minresactivity) )
                     resactisinf = TRUE;
               }
               rescap = consdata->rhs - minresactivity;
            }
            else
            {
               if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastglbmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, infervar, &maxresactivity, FALSE, TRUE);
                  if( SCIPisInfinity(scip, maxresactivity) )
                     resactisinf = TRUE;
               }
               rescap = consdata->lhs - maxresactivity;
            }

            if( reasonisrhs == (vals[inferpos] > 0.0) )
               rescap -= vals[inferpos] * (SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) + 1.0);
            else
               rescap -= vals[inferpos] * (SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) - 1.0);
         }
         else
            rescap = (reasonisrhs ? consdata->rhs - minresactivity : consdata->lhs - maxresactivity);

         if( !resactisinf )
         {
            /* now add bounds as reasons until the residual capacity is exceeded */
            for( i = 0; i < nvars; ++i )
            {
               assert(vars != NULL); /* for flexelint */
               assert(vals != NULL); /* for flexelint */

               /* zero coefficients and the infered variable can be ignored */
               if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
                  continue;

               /* check if the residual capacity is exceeded */
               if( (reasonisrhs && SCIPisFeasNegative(scip, rescap))
                  || (!reasonisrhs && SCIPisFeasPositive(scip, rescap)) )
                  break;

               /* update the residual capacity due to the local bound of this variable */
               if( reasonisrhs == (vals[i] > 0.0) )
               {
                  /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound */
                  SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
                  rescap -= vals[i] * (SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE) - SCIPvarGetLbGlobal(vars[i]));
               }
               else
               {
                  /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound */
                  SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
                  rescap -= vals[i] * (SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) - SCIPvarGetUbGlobal(vars[i]));
               }
            }
            return SCIP_OKAY;
         }
      }
   }

   /* for a bound change on a continuous variable, all locally changed bounds are responsible */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars != NULL); /* for flexelint */
      assert(vals != NULL); /* for flexelint */

      /* zero coefficients and the infered variable can be ignored */
      if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
         continue;

      if( reasonisrhs == (vals[i] > 0.0) )
      {
         /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound is responsible */
         SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
      }
      else
      {
         /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound is responsible */
         SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
      }
   }

   return SCIP_OKAY;
}

/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) activity residuals of all other variables tighten bounds of single variable
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
#ifndef NDEBUG
   SCIP_Real* vals;
#endif
   int nvars;
   int inferpos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;
#ifndef NDEBUG
   vals = consdata->vals;
   assert(vars != NULL);
   assert(vals != NULL);
#endif

   /* get the position of the inferred variable in the vars array */
   inferpos = inferInfoGetPos(inferinfo);
   if( inferpos >= nvars || vars[inferpos] != infervar )
   {
      /* find inference variable in constraint */
      /**@todo use a binary search here; the variables can be sorted by variable index */
      for( inferpos = 0; inferpos < nvars && vars[inferpos] != infervar; ++inferpos )
      {}
   }
   assert(inferpos < nvars);
   assert(vars[inferpos] == infervar);
   assert(!SCIPisZero(scip, vals[inferpos]));

   switch( inferInfoGetProprule(inferinfo) )
   {
   case PROPRULE_1_RHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not exceed the right hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_UPPER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, TRUE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_1_LHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not fall below the left hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_LOWER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, FALSE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in linear constraint <%s> at position %d for %s bound of variable <%s>\n",
         inferInfoGetProprule(inferinfo), SCIPconsGetName(cons), inferInfoGetPos(inferinfo),
         boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar));
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, ";\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bounds on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conflict detecting constraint */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the conflict? */
   )
{
   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add the conflicting bound for each variable of infeasible constraint to conflict candidate queue */
   SCIP_CALL( addConflictBounds(scip, cons, NULL, NULL, -1, reasonisrhs) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

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
   SCIP_Real val;
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

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   *cutoff = FALSE;

   var = consdata->vars[pos];

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   val = consdata->vals[pos];
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
      &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);
   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPisLE(scip, lb, ub));

   if( val > 0.0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) && !minisrelax )
      {
         SCIP_Real newub;

         newub = (rhs - minresactivity)/val;

         if( !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            SCIP_Bool activityunreliable;
            activityunreliable = SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity);

            /* check minresactivities for reliability */
            if( activityunreliable )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
               newub = (rhs - minresactivity)/val;
               activityunreliable = SCIPisInfinity(scip, -minresactivity) ||
                  (!SCIPisUbBetter(scip, newub, lb, ub) && (!SCIPisFeasLT(scip, newub, ub) || !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisLT(scip, newub, ub)));
            }

            if( !activityunreliable )
            {
               /* tighten upper bound */
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newub=%.15g\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
               SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
               if( infeasible )
               {
                  SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

                  /* analyze conflict */
                  SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
                  assert(SCIPisFeasLE(scip, ub, newub));
                  (*nchgbds)++;

                  SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
               }
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newlb;

         newlb = (lhs - maxresactivity)/val;
         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            /* check maxresactivities for reliability */
            if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
               newlb = (lhs - maxresactivity)/val;

               if( SCIPisInfinity(scip, maxresactivity) || (!SCIPisLbBetter(scip, newlb, lb, ub) 
                     && (!SCIPisFeasGT(scip, newlb, lb) || !SCIPvarIsIntegral(var)) 
                     && (!force || !SCIPisGT(scip, newlb, lb))) )
                  return SCIP_OKAY;
            }

            /* tighten lower bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
            SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasGE(scip, lb, newlb));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
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

         newlb = (rhs - minresactivity)/val;
         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            SCIP_Bool activityunreliable;
            activityunreliable = SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity);
            /* check minresactivities for reliability */
            if( activityunreliable )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
               newlb = (rhs - minresactivity)/val;

               activityunreliable = SCIPisInfinity(scip, -minresactivity) 
                  || (!SCIPisLbBetter(scip, newlb, lb, ub) && (!SCIPisFeasGT(scip, newlb, lb) || !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisGT(scip, newlb, lb)));
            }

            if( !activityunreliable )
            {
               /* tighten lower bound */
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
               SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
               if( infeasible )
               {
                  SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

                  /* analyze conflict */
                  SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
                  assert(SCIPisFeasGE(scip, lb, newlb));
                  (*nchgbds)++;
                  SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
               }
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newub;

         newub = (lhs - maxresactivity)/val;
         if(  !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            /* check maxresactivities for reliability */
            if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
               newub = (lhs - maxresactivity)/val;

               if( SCIPisInfinity(scip, maxresactivity) || (!SCIPisUbBetter(scip, newub, lb, ub) 
                     && (!SCIPisFeasLT(scip, newub, ub) && !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisLT(scip, newub, ub))) )
                  return SCIP_OKAY;
            }

            /* tighten upper bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g], newub=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
            SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasLE(scip, ub, newub));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }

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
   int nvars;
   int nrounds;
   int lastchange;
   int oldnchgbds;
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

   /* if a constraint was created after presolve, then it may hold fixed variables
    * if there are even multi-aggregated variables, then we cannot do bound tightening on these
    * thus, ensure here again that variable fixings have been applied
    */
   SCIP_CALL( applyFixings(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   force = (nvars == 1) && !SCIPconsIsModifiable(cons);

   /* ensure that the variables are properly sorted */
   if( sortvars && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->binvarssorted )
   {
      SCIP_CALL( consdataSort(scip, consdata) );
      assert(consdata->binvarssorted);
   }

   /* as long as the bounds might be tightened again, try to tighten them; abort after a maximal number of rounds */
   lastchange = -1;
   for( nrounds = 0; (force || !consdata->boundstightened) && nrounds < MAXTIGHTENROUNDS; ++nrounds )
   {
      /* mark the constraint to have the variables' bounds tightened */
      consdata->boundstightened = TRUE;

      /* try to tighten the bounds of each variable in the constraint. During solving process, 
       * the binary variable sorting enables skipping variables */
      v = 0;
      while( v < nvars && v != lastchange && !(*cutoff) )
      {
         oldnchgbds = *nchgbds;

         assert(!sortvars || SCIPgetStage(scip) < SCIP_STAGE_SOLVING || consdata->binvarssorted);

         SCIP_CALL( tightenVarBounds(scip, cons, v, cutoff, nchgbds, force) );
         
         /* if there was no progress, skip the rest of the binary variables */
         if( *nchgbds > oldnchgbds )
         {
            lastchange = v;
            v++;
         }
         else if( consdata->binvarssorted && v < consdata->nbinvars - 1
            && !SCIPisFeasEQ(scip, SCIPvarGetUbLocal(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v])) )
            v = consdata->nbinvars;
         else
            ++v;
      }
   }
#ifndef NDEBUG
   if( force && SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      assert(*cutoff || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->vars[0]), SCIPvarGetUbLocal(consdata->vars[0])));
#endif

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or current solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool             checklprows,        /**< has linear constraint to be checked, if it is already in current LP? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMessage("checking linear constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   if( consdata->row != NULL )
   {
      if( !checklprows && SCIProwIsInLP(consdata->row) )
         return SCIP_OKAY;
      else if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
         activity = consdataGetPseudoActivity(scip, consdata);
      else
         activity = SCIPgetRowSolActivity(scip, consdata->row, sol);
   }
   else
      activity = consdataGetActivity(scip, consdata, sol);

   SCIPdebugMessage("  consdata activity=%.15g (lhs=%.15g, rhs=%.15g, row=%p, checklprows=%u, rowinlp=%u, sol=%p, hascurrentnodelp=%u)\n",
      activity, consdata->lhs, consdata->rhs, (void*)consdata->row, checklprows,
      consdata->row == NULL ? 0 : SCIProwIsInLP(consdata->row), (void*)sol,
      consdata->row == NULL ? FALSE : SCIPhasCurrentNodeLP(scip));
   
   if( SCIPisFeasLT(scip, activity, consdata->lhs) || SCIPisFeasGT(scip, activity, consdata->rhs) )
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
      *violated = FALSE;

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
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRow(scip, consdata->row, consdata->nvars, consdata->vars, consdata->vals) );

   return SCIP_OKAY;
}

/** adds linear constraint as cut to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      /* convert consdata object into LP row */
      SCIP_CALL( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding relaxation of linear constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, consdata->row, NULL)) );
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by given solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             separatecards,      /**< should knapsack cardinality cuts be generated? */
   SCIP_Bool             separateall,        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   int*                  ncuts,              /**< pointer to add up the number of found cuts */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool violated;
   int oldncuts;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(ncuts != NULL);
   assert(consdata != NULL);

   oldncuts = *ncuts;
   *cutoff = FALSE;

   SCIP_CALL( checkCons(scip, cons, sol, (sol != NULL), &violated) );

   if( violated )
   {
      /* insert LP row as cut */
      SCIP_CALL( addRelaxation(scip, cons, sol) );
      (*ncuts)++;
   }
   else if( !SCIPconsIsModifiable(cons) && separatecards )
   {
      /* relax linear constraint into knapsack constraint and separate lifted cardinality cuts */
      if( !separateall && sol == NULL )
      {
         /* we only want to call the knapsack cardinality cut separator for rows that have a non-zero dual solution */
         if( consdata->row != NULL && SCIProwIsInLP(consdata->row) )
         {
            SCIP_Real dualsol;

            dualsol = SCIProwGetDualsol(consdata->row);
            if( SCIPisFeasNegative(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, consdata->rhs) )
               {
                  SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                        consdata->vals, +1.0, consdata->rhs, sol, ncuts, cutoff) );
               }
            }
            else if( SCIPisFeasPositive(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) )
               {
                  SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                        consdata->vals, -1.0, -consdata->lhs, sol, ncuts, cutoff) );
               }
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                  consdata->vals, +1.0, consdata->rhs, sol, ncuts, cutoff) );
         }
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                  consdata->vals, -1.0, -consdata->lhs, sol, ncuts, cutoff) );
         }
      }
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

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /*SCIPdebugMessage("propagating linear constraint <%s>\n", SCIPconsGetName(cons));*/

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated && (!tightenbounds || consdata->boundstightened) )
      return SCIP_OKAY;

   /* mark constraint to be propagated */
   consdata->propagated = TRUE;

   /* we can only infer activity bounds of the linear constraint, if it is not modifiable */
   if( !SCIPconsIsModifiable(cons) )
   {
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
         consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible (rhs): activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible (lhs): activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisGE(scip, minactivity, consdata->lhs) && SCIPisLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Presolving methods
 */

/** extracts cliques of the constraint and adds them to SCIP */
static
SCIP_RETCODE extractCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsclique;
   SCIP_Bool rhsclique;
   int i;
   int nposcoefs;
   int nnegcoefs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we already added the cliques of the constraints */
   if( consdata->cliquesadded )
      return SCIP_OKAY;

   consdata->cliquesadded = TRUE;

   /* sort variables by variable type */
   SCIP_CALL( consdataSort(scip, consdata) );

   /* currently, we only check whether the constraint is a set packing / partitioning constraint */
   /**@todo extract more cliques from linear constraints */

   /* check if all variables are binary, if the coefficients are +1 or -1, and if the right hand side is equal
    * to 1 - number of negative coefficients, or if the left hand side is equal to number of positive coefficients - 1
    */
   nposcoefs = 0;
   nnegcoefs = 0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( !SCIPvarIsBinary(consdata->vars[i]) )
         return SCIP_OKAY;
      else if( SCIPisEQ(scip, consdata->vals[i], +1.0) )
         nposcoefs++;
      else if( SCIPisEQ(scip, consdata->vals[i], -1.0) )
         nnegcoefs++;
      else
         return SCIP_OKAY;
   }
   lhsclique = SCIPisEQ(scip, consdata->lhs, (SCIP_Real)nposcoefs - 1.0);
   rhsclique = SCIPisEQ(scip, consdata->rhs, 1.0 - (SCIP_Real)nnegcoefs);
   if( lhsclique || rhsclique )
   {
      SCIP_Bool* values;
      SCIP_Bool infeasible;
      int nbdchgs;

      SCIPdebugMessage("linear constraint <%s>: adding clique with %d vars (%d pos, %d neg)\n",
         SCIPconsGetName(cons), consdata->nvars, nposcoefs, nnegcoefs);
      SCIP_CALL( SCIPallocBufferArray(scip, &values, consdata->nvars) );
      for( i = 0; i < consdata->nvars; ++i )
         values[i] = (rhsclique == (consdata->vals[i] > 0.0));
      SCIP_CALL( SCIPaddClique(scip, consdata->vars, values, consdata->nvars, &infeasible, &nbdchgs) );
      if( infeasible )
         *cutoff = TRUE;
      *nchgbds += nbdchgs;
      SCIPfreeBufferArray(scip, &values);
   }

   return SCIP_OKAY;
}

/** tightens left and right hand side of constraint due to integrality */
static
SCIP_RETCODE tightenSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool integral;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisIntegral(scip, consdata->lhs) || !SCIPisIntegral(scip, consdata->rhs) )
   {
      integral = TRUE;
      for( i = 0; i < consdata->nvars && integral; ++i )
      {
         integral = SCIPisIntegral(scip, consdata->vals[i])
            && (SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS);
      }
      if( integral )
      {
         SCIPdebugMessage("linear constraint <%s>: make sides integral: sides=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisIntegral(scip, consdata->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons, SCIPfeasCeil(scip, consdata->lhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisIntegral(scip, consdata->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons, SCIPfeasFloor(scip, consdata->rhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         SCIPdebugMessage("linear constraint <%s>: new integral sides: sides=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
      }
   }

   return SCIP_OKAY;
}

#define MAXVALRECOMP 1e+06

/** tightens coefficients of binary, integer, and implicit integer variables due to activity bounds in presolving:
 *  given an inequality  lhs <= a*x + ai*xi <= rhs, with a non-continuous variable  li <= xi <= ui
 *  let minact := min{a*x + ai*xi}, maxact := max{a*x + ai*xi}
 *  (i) ai >= 0:
 *      if  minact + ai >= lhs  and  maxact - ai <= rhs: (**)
 *       - a deviation from the lower/upper bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := max(lhs - minact, maxact - rhs)
 *           lhs' := lhs - (ai - ai')*li
 *           rhs' := rhs - (ai - ai')*ui
 * (ii) ai < 0:
 *      if  minact - ai >= lhs  and  maxact + ai <= rhs: (***)
 *       - a deviation from the upper/lower bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := min(rhs - maxact, minact - lhs)
 *           lhs' := lhs - (ai - ai')*ui
 *           rhs' := rhs - (ai - ai')*li
 *
 *  We further try to remove redundant variable from the constraint;
 *  Variables which fulfill conditions (**) or (***) are called surely non-redundant variables.
 *  A deviation of only one from their bound makes the lhs/rhs feasible (i.e., redundant), even if all other 
 *  variables are set to their "worst" bound. If all variables which are not surely non-redundant cannot make 
 *  the lhs/rhs redundant, even if they are set to their "best" bound, they can be removed from the constraint.
 *  E.g., for binary variables and an inequality x_1 +x_2 +10y_1 +10y_2 >= 5, setting either of the y_i to one 
 *  suffices to fulfill the inequality, whereas the x_i do not contribute to feasibility and can be removed.
 */
static
SCIP_RETCODE consdataTightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real minactivity; /* minimal value w.r.t. the variable's local bounds for the constraint's
                           * activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real maxactivity; /* maximal value w.r.t. the variable's local bounds for the constraint's
                           * activity, ignoring the coefficients contributing with infinite value */
   SCIP_Bool minactisrelax; /* do huge finite values contribute to the minactivity? */
   SCIP_Bool maxactisrelax; /* do huge finite values contribute to the maxactivity? */
   SCIP_Real minleftactivity; /* minimal activity without surely non-redundant variables. */
   SCIP_Real maxleftactivity; /* maximal activity without surely non-redundant variables. */
   SCIP_Real aggrlhs; /* lhs without minimal activity of surely non-redundant variables. */
   SCIP_Real aggrrhs; /* rhs without maximal activity of surely non-redundant variables. */
   SCIP_Real lval; /* candidate for new value arising from considering the left hand side */
   SCIP_Real rval; /* candidate for new value arising from considering the left hand side */
   SCIP_Real val;
   SCIP_Real newval;
   SCIP_Real newlhs;
   SCIP_Real newrhs;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* @todo Is this still needed with automatic recomputation of activities? */
   /* if the maximal coefficient is too large, recompute the activities */
   if( consdata->maxabsval > MAXVALRECOMP )
   {
      consdataRecomputeMinactivity(scip, consdata);
      consdataRecomputeMaxactivity(scip, consdata);
   }

   /* get the minimal and maximal activity of the constraint */
   consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

   minleftactivity = 0.0;
   maxleftactivity = 0.0;

   /* try to tighten each coefficient */
   i = 0;
   while( i < consdata->nvars )
   {
      var = consdata->vars[i];

      /* get coefficient and variable's bounds */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      val = consdata->vals[i];
      assert(!SCIPisZero(scip, val));

      /* check sign of coefficient */
      if( val >= 0.0 )
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS &&
            SCIPisGE(scip, minactivity + val, consdata->lhs) && SCIPisLE(scip, maxactivity - val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := max(lhs - minact, maxact - rhs)
             *   lhs' := lhs - (ai - ai')*li
             *   rhs' := rhs - (ai - ai')*ui
             */

            lval = consdata->lhs - minactivity;
            rval = maxactivity - consdata->rhs;

            /* Try to avoid cancellation, if there are only two variables */
            if( consdata->nvars == 2 )
            {
               SCIP_Real otherval;
               otherval = consdata->vals[1-i];

               if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->minactivityneginf + consdata->minactivityneginf == 0 )
               {
                  lval = consdata->lhs - val*lb;
                  lval -= otherval > 0.0 ? otherval * SCIPvarGetLbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetUbLocal(consdata->vars[1-i]);
               }

               if( !SCIPisInfinity(scip,consdata->rhs) && consdata->maxactivityneginf + consdata->maxactivityneginf == 0 )
               {
                  rval = val*ub - consdata->rhs;
                  rval += otherval > 0.0 ? otherval * SCIPvarGetUbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetLbLocal(consdata->vars[1-i]);
               }
            }

            newval = MAX(lval, rval);
            assert(SCIPisSumRelLE(scip, newval, val));

            /* Try to avoid cancellation in computation of lhs/rhs */
            newlhs = consdata->lhs - val * lb;
            newlhs += newval * lb;
            newrhs = consdata->rhs - val * ub;
            newrhs += newval * ub;

            if( !SCIPisSumRelEQ(scip, newval, val) )
            {
               SCIPdebugMessage("linear constraint <%s>: change coefficient %+.15g<%s> to %+.15g<%s>, act=[%.15g,%.15g], side=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
               {
                  SCIPdebugMessage("linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

                  SCIP_CALL( chgLhs(scip, cons, newlhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->lhs, newlhs));
               }

               if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
               {
                  SCIPdebugMessage("linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

                  SCIP_CALL( chgRhs(scip, cons, newrhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->rhs, newrhs));
               }
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -minleftactivity) )
            {
               assert(!SCIPisInfinity(scip, val));
               assert(!SCIPisInfinity(scip, lb));
               if( SCIPisInfinity(scip, -lb) )
                  minleftactivity = -SCIPinfinity(scip);
               else
                  minleftactivity += val * lb;
            }

            if( !SCIPisInfinity(scip, maxleftactivity) )
            {
               assert(!SCIPisInfinity(scip, val));
               assert(!SCIPisInfinity(scip, -ub));
               if( SCIPisInfinity(scip,ub) )
                  maxleftactivity = SCIPinfinity(scip);
               else
                  maxleftactivity += val * ub;
            }
         }
      }
      else
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS &&
            SCIPisGE(scip, minactivity - val, consdata->lhs) && SCIPisLE(scip, maxactivity + val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := min(rhs - maxact, minact - lhs)
             *   lhs' := lhs - (ai - ai')*ui
             *   rhs' := rhs - (ai - ai')*li
             */

            lval = minactivity - consdata->lhs;
            rval = consdata->rhs - maxactivity;

            /* Try to avoid cancellation, if there are only two variables */
            if( consdata->nvars == 2 )
            {
               SCIP_Real otherval;
               otherval = consdata->vals[1-i];

               if( !SCIPisInfinity(scip,-consdata->lhs) && consdata->minactivityneginf + consdata->minactivityneginf == 0 )
               {
                  lval = val*ub - consdata->lhs;
                  lval += otherval > 0.0 ? otherval * SCIPvarGetLbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetUbLocal(consdata->vars[1-i]);
               }

               if( !SCIPisInfinity(scip,consdata->rhs) && consdata->maxactivityneginf + consdata->maxactivityneginf == 0 )
               {
                  rval = consdata->rhs - val*lb;
                  rval -= otherval > 0.0 ? otherval * SCIPvarGetUbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetLbLocal(consdata->vars[1-i]);
               }
            }

            newval = MIN(lval, rval);
            assert(SCIPisSumRelGE(scip, newval, val));

            /* Try to avoid cancellation in computation of lhs/rhs */
            newlhs = consdata->lhs - val * ub;
            newlhs += newval * ub;
            newrhs = consdata->rhs - val * lb;
            newrhs += newval * lb;

            if( !SCIPisSumRelEQ(scip, newval, val) )
            {
               SCIPdebugMessage("linear constraint <%s>: change coefficient %+.15g<%s> to %+.15g<%s>, act=[%.15g,%.15g], side=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
               {
                  SCIPdebugMessage("linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

                  SCIP_CALL( chgLhs(scip, cons, newlhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->lhs, newlhs));
               }

               if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
               {
                  SCIPdebugMessage("linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

                  SCIP_CALL( chgRhs(scip, cons, newrhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->rhs, newrhs));
               }
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -minleftactivity) )
            {
               assert(!SCIPisInfinity(scip, -val));
               assert(!SCIPisInfinity(scip, -ub));
               if( SCIPisInfinity(scip, ub) )
                  minleftactivity = -SCIPinfinity(scip);
               else
                  minleftactivity += val * ub;
            }

            if( !SCIPisInfinity(scip, maxleftactivity) )
            {
               assert(!SCIPisInfinity(scip, -val));
               assert(!SCIPisInfinity(scip, lb));
               if( SCIPisInfinity(scip, -lb) )
                  maxleftactivity = SCIPinfinity(scip);
               else
                  maxleftactivity += val * lb;
            }
         }
      }
      ++i;
   }

   SCIPdebugMessage("minleftactivity = %.15g, rhs = %.15g\n",
      minleftactivity, consdata->rhs);
   SCIPdebugMessage("maxleftactivity = %.15g, lhs = %.15g\n", 
      maxleftactivity, consdata->lhs);

   /* minleft == \infty  ==>  minactivity == \infty */
   assert(!SCIPisInfinity(scip, -minleftactivity) || SCIPisInfinity(scip, -minactivity));
   assert(!SCIPisInfinity(scip, maxleftactivity) || SCIPisInfinity(scip, maxactivity));

   /* if the lhs is finite, we will check in the following whether the not non-redundant variables can make lhs feasible;
    * this is not valid, if the minactivity is -\infty (aggrlhs would be minus infinity in the following computation)
    * or if huge values contributed to the minactivity, because the minactivity is then just a relaxation
    * (<= the exact minactivity), and we might falsely claim variables to be redundant in the following
    */
   assert(!SCIPisInfinity(scip, minactivity));
   if( !SCIPisInfinity(scip, -consdata->lhs) && (SCIPisInfinity(scip, -minactivity) || minactisrelax) )
      return SCIP_OKAY;

   /* if the rhs is finite, we will check in the following whether the not non-redundant variables can make rhs feasible;
    * this is not valid, if the maxactivity is \infty (aggrrhs would be infinity in the following computation)
    * or if huge values contributed to the maxactivity, because the maxactivity is then just a relaxation
    * (>= the exact maxactivity), and we might falsely claim variables to be redundant in the following
    */
   assert(!SCIPisInfinity(scip, -maxactivity));
   if( !SCIPisInfinity(scip, consdata->rhs) && (SCIPisInfinity(scip, maxactivity) || maxactisrelax) )
      return SCIP_OKAY;

   /* correct lhs and rhs by min/max activity of surely non-redundant variables 
    * surely non-redundant variables are all those where a deviation from the bound makes the lhs/rhs redundant
    */
   aggrlhs = consdata->lhs - minactivity + minleftactivity;
   aggrrhs = consdata->rhs - maxactivity + maxleftactivity;

   /* check if the constraint contains variables which are redundant. The reasoning is the following:
    * Each non-redundant variable can make the lhs/rhs feasible with a deviation of only one in the bound.
    * If _all_ variables which are not non-redundant together cannot make lhs/rhs feasible, 
    * they can be removed from the constraint.
    * aggrrhs may contain some near-infinity value, but only if rhs is infinity.
    */
   if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasLT(scip, maxleftactivity, aggrlhs))
      && (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasGT(scip, minleftactivity, aggrrhs)) )
   {
      SCIP_Real minleftactivitypart;
      SCIP_Real maxleftactivitypart;
         
      assert(!SCIPisInfinity(scip, -consdata->lhs) || !SCIPisInfinity(scip, consdata->rhs));
      
      /* try to remove redundant variables from constraint */
      i = 0;
      while( i < consdata->nvars )
      {
         var = consdata->vars[i];
         minleftactivitypart = 0.0;
         maxleftactivitypart = 0.0;
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         /* get coefficient and variable's bounds */
         val = consdata->vals[i];
         assert(!SCIPisZero(scip, val));
         
         /* check sign of coefficient */
         if( val >= 0.0 )
         {     
            /* negation of condition above in case of positive val */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || 
               SCIPisLT(scip, minactivity + val, consdata->lhs) || SCIPisGT(scip, maxactivity - val, consdata->rhs) )
            {
               SCIPdebugMessage("minactivity = %g\tval = %g\tlhs = %g\n", minactivity, val, consdata->lhs);
               SCIPdebugMessage("maxactivity = %g\tval = %g\trhs = %g\n", maxactivity, val, consdata->rhs);
               SCIPdebugMessage("linear constraint <%s>: remove variable <%s> with coefficient <%g> from constraint since it is redundant\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), val);

               minleftactivitypart = val * lb;
               maxleftactivitypart = val * ub;
               
               SCIP_CALL( delCoefPos(scip, cons, i) );
               i--;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               /* we return above if the condition does not hold and deleting a variable cannot increase the number of
                * huge contributions
                */
               assert(!minactisrelax || SCIPisInfinity(scip, -consdata->lhs));
               assert(!maxactisrelax || SCIPisInfinity(scip, consdata->rhs));
            }
         }
         else 
         {
            /* negation of condition above in case of negative val */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || 
               SCIPisLT(scip, minactivity - val, consdata->lhs) || SCIPisGT(scip, maxactivity + val, consdata->rhs) )
            {               
               SCIPdebugMessage("linear constraint <%s>: remove variable <%s> with coefficient <%g> from constraint since it is redundant\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), val);

               minleftactivitypart = val * ub;
               maxleftactivitypart = val * lb;
               
               SCIP_CALL( delCoefPos(scip, cons, i) );
               i--;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               /* we return above if the condition does not hold and deleting a variable cannot increase the number of
                * huge contributions
                */
               assert(!minactisrelax || SCIPisInfinity(scip, -consdata->lhs));
               assert(!maxactisrelax || SCIPisInfinity(scip, consdata->rhs));
            }
         }

         /* the following update step is needed in every iteration cause otherwise it is possible that the surely none-
          * redundant variables could get deleted, 
          * e.g. y_1 + 16y_2 >= 25, y1 with bounds [9,12], y2 with bounds [0,2], minactivity would be 9, it follows that
          * y_2 is surely not redundant and y_1 is redundant so we would first delete y1 and without updating the sides
          * we would also delete y2 and as a result we would have gotten infeasibility */
         /* adjust lhs and right hand side */
         newlhs = consdata->lhs - minleftactivitypart;
         newrhs = consdata->rhs - maxleftactivitypart;

         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisFeasEQ(scip, newlhs, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);         
            SCIP_CALL( chgLhs(scip, cons, newlhs) );
            ++(*nchgsides);
            assert(SCIPisEQ(scip, consdata->lhs, newlhs));
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisFeasEQ(scip, newrhs, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);
            SCIP_CALL( chgRhs(scip, cons, newrhs) );
            ++(*nchgsides);
            assert(SCIPisEQ(scip, consdata->rhs, newrhs));
         }
         ++i;
      }
   }
   
   return SCIP_OKAY;
}

/* processes equality with only one variable by fixing the variable and deleting the constraint */
static
SCIP_RETCODE convertUnaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real fixval;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 1);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   /* calculate the value to fix the variable to */
   var = consdata->vars[0];
   val = consdata->vals[0];
   assert(!SCIPisZero(scip, val));
   fixval = SCIPselectSimpleValue(consdata->lhs/val - 0.9 * SCIPepsilon(scip),
      consdata->rhs/val + 0.9 * SCIPepsilon(scip), MAXDNOM);
   SCIPdebugMessage("linear equality <%s>: fix <%s> == %.15g\n",
      SCIPconsGetName(cons), SCIPvarGetName(var), fixval);

   /* fix variable */
   SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
   if( infeasible )
   {
      SCIPdebugMessage(" -> infeasible fixing\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   if( fixed )
      (*nfixedvars)++;
   
   /* disable constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   if( !consdata->upgraded )
      (*ndelconss)++;

   return SCIP_OKAY;
}

/* processes equality with exactly two variables by aggregating one of the variables and deleting the constraint */
static
SCIP_RETCODE convertBinaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMessage("linear constraint <%s>: aggregate %.15g<%s> + %.15g<%s> == %.15g\n",
      SCIPconsGetName(cons), consdata->vals[0], SCIPvarGetName(consdata->vars[0]),
      consdata->vals[1], SCIPvarGetName(consdata->vars[1]), consdata->rhs);

   /* aggregate the equality */
   SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], consdata->vals[0], consdata->vals[1],
         consdata->rhs, &infeasible, &redundant, &aggregated) );

   /* check for infeasibility of aggregation */
   if( infeasible )
   {
      SCIPdebugMessage(" -> infeasible aggregation\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* count the aggregation */
   if( aggregated )
      (*naggrvars)++;

   /* delete the constraint, if it is redundant */
   if( redundant )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );

      if( !consdata->upgraded )
         (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/** calculates the new lhs and rhs of the constraint after the given variable is aggregated out */
static
void getNewSidesAfterAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             slackvar,           /**< variable to be aggregated out */
   SCIP_Real             slackcoef,          /**< coefficient of variable in constraint */
   SCIP_Real*            newlhs,             /**< pointer to store new lhs of constraint */
   SCIP_Real*            newrhs              /**< pointer to store new rhs of constraint */
   )
{
   SCIP_Real slackvarlb;
   SCIP_Real slackvarub;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(newlhs != NULL);
   assert(newrhs != NULL);
   assert(!SCIPisInfinity(scip, -consdata->lhs));
   assert(!SCIPisInfinity(scip, consdata->rhs));

   slackvarlb = SCIPvarGetLbGlobal(slackvar);
   slackvarub = SCIPvarGetUbGlobal(slackvar);
   if( slackcoef > 0.0 )
   {
      if( SCIPisInfinity(scip, -slackvarlb) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = consdata->rhs - slackcoef * slackvarlb;
      if( SCIPisInfinity(scip, slackvarub) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = consdata->lhs - slackcoef * slackvarub;
   }
   else
   {
      if( SCIPisInfinity(scip, -slackvarlb) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = consdata->rhs - slackcoef * slackvarlb;
      if( SCIPisInfinity(scip, slackvarub) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = consdata->lhs - slackcoef * slackvarub;
   }
   assert(SCIPisLE(scip, *newlhs, *newrhs));
}

/* processes equality with more than two variables by multi-aggregating one of the variables and converting the equality
 * into an inequality; if multi-aggregation is not possible, tries to identify one continuous or integer variable that is
 * implicitly integral by this constraint
 *
 * @todo Check whether a more clever way of avoiding aggregation of variables containing implicitly integer variables can help.
 */
static
SCIP_RETCODE convertLongEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_VARTYPE bestslacktype;
   SCIP_VARTYPE slacktype;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real bestslackdomrng;
   SCIP_Bool bestremovescons;
   SCIP_Bool coefszeroone;
   SCIP_Bool coefsintegral;
   SCIP_Bool varsintegral;
   SCIP_Bool supinf;                         /* might the supremum of the multi-aggregation be infinite? */
   SCIP_Bool infinf;                         /* might the infimum of the multi-aggregation be infinite? */
   SCIP_Bool infeasible;

   int maxnlocksstay;
   int maxnlocksremove;
   int bestslackpos;
   int bestnlocks;
   int ncontvars;
   int contvarpos;
   int nintvars;
   int nimplvars;
   int intvarpos;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMessage("linear constraint <%s>: try to multi-aggregate equality\n", SCIPconsGetName(cons));

   /* We do not want to increase the total number of non-zeros due to the multi-aggregation.
    * Therefore, we have to restrict the number of locks of a variable that is aggregated out.
    *   maxnlocksstay:   maximal sum of lock numbers if the constraint does not become redundant after the aggregation
    *   maxnlocksremove: maximal sum of lock numbers if the constraint can be deleted after the aggregation
    */
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   maxnlocksstay = 0;
   if( consdata->nvars == 3 )
   {
      /* If the constraint becomes redundant, 3 non-zeros are removed, and we get 1 additional non-zero for each
       * constraint the variable appears in. Thus, the variable must appear in at most 3 other constraints.
       */
      maxnlocksremove = 3;
   }
   else if( consdata->nvars == 4 )
   {
      /* If the constraint becomes redundant, 4 non-zeros are removed, and we get 2 additional non-zeros for each
       * constraint the variable appears in. Thus, the variable must appear in at most 2 other constraints.
       */
      maxnlocksremove = 2;
   }
   else
   {
      /* If the constraint is redundant but has more than 4 variables, we can only accept one other constraint. */
      maxnlocksremove = 1;
   }

   /* the locks on this constraint can be ignored */
   if( SCIPconsIsChecked(cons) )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         maxnlocksstay++;
         maxnlocksremove++;
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         maxnlocksstay++;
         maxnlocksremove++;
      }
   }

   /* look for a slack variable s to convert a*x + s == b into lhs <= a*x <= rhs */
   vars = consdata->vars;
   vals = consdata->vals;
   bestslackpos = -1;
   bestslacktype = SCIP_VARTYPE_BINARY;
   bestnlocks = INT_MAX;
   bestremovescons = FALSE;
   bestslackdomrng = 0.0;
   coefszeroone = TRUE;
   coefsintegral = TRUE;
   varsintegral = TRUE;
   ncontvars = 0;
   contvarpos = -1;
   nintvars = 0;
   nimplvars = 0;
   intvarpos = -1;
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real absval;
      SCIP_Real varlb;
      SCIP_Real varub;
      SCIP_Bool iscont;
      int nlocks;

      assert(vars != NULL);
      assert(vals != NULL);

      var = vars[v];
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksDown(var) >= 1); /* because variable is locked in this equality */
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksUp(var) >= 1);
      varlb = SCIPvarGetLbGlobal(var);
      varub = SCIPvarGetUbGlobal(var);

      val = vals[v];
      absval = REALABS(val);
      assert(SCIPisPositive(scip, absval));

      slacktype = SCIPvarGetType(var);
      coefszeroone = coefszeroone && SCIPisEQ(scip, absval, 1.0);
      coefsintegral = coefsintegral && SCIPisIntegral(scip, val);
      varsintegral = varsintegral && (slacktype != SCIP_VARTYPE_CONTINUOUS);
      iscont = (slacktype == SCIP_VARTYPE_CONTINUOUS || slacktype == SCIP_VARTYPE_IMPLINT);
      if ( slacktype == SCIP_VARTYPE_IMPLINT )
         ++nimplvars;

      /* update candidates for continuous -> implint and integer -> implint conversion */
      if( iscont )
      {
         ncontvars++;
         contvarpos = v;
      }
      else if( slacktype == SCIP_VARTYPE_INTEGER )
      {
         nintvars++;
         intvarpos = v;
      }

      /* check, if variable is already fixed or aggregated */
      if( !SCIPvarIsActive(var) )
         continue;

      /* check, if variable is used in too many other constraints, even if this constraint could be deleted */
      nlocks = SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var);
      if( nlocks > maxnlocksremove )
         continue;

      /* check, if variable can be used as a slack variable */
      if( (iscont || (coefsintegral && varsintegral && SCIPisEQ(scip, absval, 1.0))) &&
         !SCIPdoNotMultaggrVar(scip, var) )
      {
         SCIP_Bool better;
         SCIP_Bool equal;
         SCIP_Real slackdomrng;

         if( SCIPisInfinity(scip, varub) || SCIPisInfinity(scip, -varlb) )
            slackdomrng = SCIPinfinity(scip);
         else
         {
            slackdomrng = (varub - varlb)*absval;
            assert(!SCIPisInfinity(scip, slackdomrng));
         }
         equal = FALSE;
         better = (slacktype > bestslacktype) || (bestslackpos == -1);
         if( !better && slacktype == bestslacktype )
         {
            better = (nlocks < bestnlocks);
            if( nlocks == bestnlocks && !bestremovescons )
            {
               better = SCIPisGT(scip, slackdomrng, bestslackdomrng);
               equal = !better && SCIPisGE(scip, slackdomrng, bestslackdomrng);
            }
         }

         if( better || equal )
         {
            SCIP_Real minresactivity;
            SCIP_Real maxresactivity;
            SCIP_Real newlhs;
            SCIP_Real newrhs;
            SCIP_Bool removescons;
            SCIP_Bool minisrelax;
            SCIP_Bool maxisrelax;
            SCIP_Bool isminsettoinfinity;
            SCIP_Bool ismaxsettoinfinity;

            /* check if the constraint becomes redundant after multi-aggregation */
            consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
               &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);

            /* do not perform the multi-aggregation due to numerics, if we have huge contributions in the residual
             * activity
             */
            if( minisrelax || maxisrelax )
               continue;

            getNewSidesAfterAggregation(scip, consdata, var, val, &newlhs, &newrhs);
            removescons = (SCIPisFeasLE(scip, newlhs, minresactivity) && SCIPisFeasLE(scip, maxresactivity, newrhs));

            /* check resactivities for reliability */
            if( removescons )
            {
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);

               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity)
                  && SCIPisFeasLE(scip, newlhs, minresactivity))
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);

               removescons = (SCIPisFeasLE(scip, newlhs, minresactivity) && SCIPisFeasLE(scip, maxresactivity, newrhs));
            }

            /* prefer variables that make the constraints redundant */
            if( bestremovescons && !removescons )
               continue;

            /* if the constraint does not become redundant, only accept the variable if it does not appear in
             * other constraints
             */
            if( !removescons && nlocks > maxnlocksstay )
               continue;

            better = better || (!bestremovescons && removescons);
            if( better )
            {
               bestslackpos = v;
               bestslacktype = slacktype;
               bestnlocks = nlocks;
               bestslackdomrng = slackdomrng;
               bestremovescons = removescons;
            }
         }
      }
   }

   /* if all coefficients and variables are integral, the right hand side must also be integral */
   if( coefsintegral && varsintegral && !SCIPisFeasIntegral(scip, consdata->rhs) )
   {
      SCIPdebugMessage("linear equality <%s> is integer infeasible\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   supinf = FALSE;
   infinf = FALSE;

   /* check whether the the infimum and the supremum of the multi-aggregation can be get infinite */
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != bestslackpos )
      {
         if( SCIPisGT(scip, consdata->vals[v], 0.0) )
         {
            supinf = supinf || SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[v]));
            infinf = infinf || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[v]));
         }
         else if( SCIPisLT(scip, consdata->vals[v], 0.0) )
         {
            supinf = supinf || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[v]));
            infinf = infinf || SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[v]));
         }
      }
   }
 
   /* If the infimum and the supremum of a multi-aggregation are both infinite, then the multi-aggregation might not be resolvable.
    * E.g., consider the equality z = x-y. If x and y are both fixed to +infinity, the value for z is not determined */     
   if( supinf && infinf )
   {      
      SCIPdebugMessage("do not perform multi-aggregation: infimum and supremum are both infinite\n");     
      return SCIP_OKAY;
   }

   /* if the slack variable is of integer type, and the constraint itself may take fractional values,
    * we cannot aggregate the variable, because the integrality condition would get lost
    * Similarly, if there are implicitly integral variables we cannot aggregate, since we might
    * loose the integrality condition for this variable.
    */
   if( bestslackpos >= 0
      && (bestslacktype == SCIP_VARTYPE_CONTINUOUS || bestslacktype == SCIP_VARTYPE_IMPLINT
         || (coefsintegral && varsintegral && nimplvars == 0)) )
   {
      SCIP_VAR* slackvar;
      SCIP_Real* scalars;
      SCIP_Real slackcoef;
      SCIP_Real aggrconst;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
      SCIP_Bool aggregated;

      /* we found a slack variable that only occurs in at most one other constraint:
       *   a_1*x_1 + ... + a_k*x_k + a'*s == rhs  ->  s == rhs - a_1/a'*x_1 - ... - a_k/a'*x_k
       */
      assert(bestslackpos < consdata->nvars);

      /* do not multi aggregate binary variables */
      if( SCIPvarIsBinary(vars[bestslackpos]) )
         return SCIP_OKAY;

      /* convert equality into inequality by deleting the slack variable:
       *  x + a*s == b, l <= s <= u   ->  b - a*u <= x <= b - a*l
       */
      slackvar = vars[bestslackpos];
      slackcoef = vals[bestslackpos];
      assert(!SCIPisZero(scip, slackcoef));
      aggrconst = consdata->rhs/slackcoef;

      getNewSidesAfterAggregation(scip, consdata, slackvar, slackcoef, &newlhs, &newrhs);
      assert(SCIPisLE(scip, newlhs, newrhs));
      SCIP_CALL( chgLhs(scip, cons, newlhs) );
      SCIP_CALL( chgRhs(scip, cons, newrhs) );
      SCIP_CALL( delCoefPos(scip, cons, bestslackpos) );

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &scalars, consdata->nvars) );

      /* set up the multi-aggregation */
      SCIPdebugMessage("linear constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(slackvar));
      for( v = 0; v < consdata->nvars; ++v )
      {
         scalars[v] = -consdata->vals[v]/slackcoef;
         SCIPdebugPrintf(" %+.15g<%s>", scalars[v], SCIPvarGetName(vars[v]));
      }
      SCIPdebugPrintf(" %+.15g, bounds of <%s>: [%.15g,%.15g], nlocks=%d, maxnlocks=%d, removescons=%u\n",
         aggrconst, SCIPvarGetName(slackvar), SCIPvarGetLbGlobal(slackvar), SCIPvarGetUbGlobal(slackvar),
         bestnlocks, bestremovescons ? maxnlocksremove : maxnlocksstay, bestremovescons);

      /* perform the multi-aggregation */
      SCIP_CALL( SCIPmultiaggregateVar(scip, slackvar, consdata->nvars, vars, scalars, aggrconst,
            &infeasible, &aggregated) );
      assert(aggregated);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &scalars);

      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMessage("linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      (*naggrvars)++;

      /* delete the constraint if it became redundant */
      if( bestremovescons )
      {
         SCIPdebugMessage("linear constraint <%s>: redundant after multi-aggregation\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         
         if( !consdata->upgraded )
            (*ndelconss)++;
      }
   }
   else if( ncontvars == 1 )
   {
      SCIP_VAR* var;

      assert(0 <= contvarpos && contvarpos < consdata->nvars);
      var = vars[contvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);

      if( coefsintegral
         && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         && SCIPisEQ(scip, REALABS(vals[contvarpos]), 1.0)
         && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* convert the continuous variable with coefficient 1.0 into an implicit integer variable */
         SCIPdebugMessage("linear constraint <%s>: converting continuous variable <%s> to implicit integer variable\n",
            SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );
         if( infeasible )
         {
            SCIPdebugMessage("infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
            *cutoff = TRUE;

            return SCIP_OKAY;
         }
      }
   }
   else if( ncontvars == 0 && nintvars == 1 && !coefszeroone )
   {
      SCIP_VAR* var;

      /* this seems to help for rococo instances, but does not for rout (where all coefficients are +/- 1.0)
       *  -> we don't convert integers into implints if the row is a 0/1-row
       */
      assert(varsintegral);
      assert(0 <= intvarpos && intvarpos < consdata->nvars);
      var = vars[intvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      if( coefsintegral
         && SCIPisEQ(scip, REALABS(vals[intvarpos]), 1.0)
         && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* convert the integer variable with coefficient 1.0 into an implicit integer variable */
         SCIPdebugMessage("linear constraint <%s>: converting integer variable <%s> to implicit integer variable\n",
            SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );
         if( infeasible )
         {
            SCIPdebugMessage("infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
            *cutoff = TRUE;

            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** checks if the given variables and their coefficient are equal (w.r.t. scaling factor) to the objective function */
static
SCIP_Bool checkEqualObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real*            scale,              /**< pointer to store the scaling factor between the constraint and the
					      *   objective function
					      */
   SCIP_Real*            offset              /**< pointer to store the offset of the objective function resulting by
					      *   this constraint
					      */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real objval;
   SCIP_Bool negated;
   int nvars;
   int v;

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(vars != NULL);

   for( v = 0; v < nvars; ++v )
   {
      negated = FALSE;
      var = vars[v];
      assert(vars != NULL);

      if( SCIPvarIsNegated(var) )
      {
         negated = TRUE;
         var = SCIPvarGetNegatedVar(var);
         assert(var != NULL);
      }

      objval = SCIPvarGetObj(var);

      /* if a variable has a zero objective coefficient the linear constraint is not a subset of the objective
       * function
       */
      if( SCIPisZero(scip, objval) )
         return FALSE;
      else
      {
         SCIP_Real val;

         val = consdata->vals[v];

         if( negated )
         {
            if( v == 0 )
            {
               /* the first variable defines the scale */
               (*scale) = val / -objval;

               (*offset) += val;
            }
            else if( SCIPisEQ(scip, -objval * (*scale), val) )
               (*offset) += val;
            else
               return FALSE;
         }
         else if( v == 0 )
         {
            /* the first variable defines the scale */
            (*scale) = val / objval;
         }
         else if( !SCIPisEQ(scip, objval * (*scale), val) )
            return FALSE;
      }
   }

   return TRUE;
}

/** check if the linear equality constraint is equal to a subset of the objective function; if so we can remove the
 *  objective coefficients and add an objective offset
 */
static
SCIP_RETCODE checkPartialObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear equation constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real offset;
   SCIP_Real scale;
   SCIP_Bool applicable;
   int nobjvars;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   nvars = consdata->nvars;
   nobjvars = SCIPgetNObjVars(scip);

   /* check if the linear equality constraints does not have more variables than the objective function */
   if( nvars > nobjvars || nvars == 0 )
      return SCIP_OKAY;

   offset = consdata->rhs;
   scale = 1.0;

   /* checks if the variables and their coefficients are equal (w.r.t. scaling factor) to the objective function */
   applicable = checkEqualObjective(scip, consdata, &scale, &offset);

   if( applicable )
   {
      SCIP_VAR** vars;

      vars = consdata->vars;
      assert(vars != NULL);

      offset /= scale;

      SCIPdebugMessage("linear equality constraint <%s> == %g (offset %g) is a subset of the objective function\n",
         SCIPconsGetName(cons), consdata->rhs, offset);

      /* set all objective coefficient to zero */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[v], 0.0) );
      }

      /* add an objective offset */
      SCIP_CALL( SCIPaddObjoffset(scip, offset) );
   }

   return SCIP_OKAY;
}

/** updates the cutoff if the given primal bound  (which is implied by the given constraint) is better */
static
SCIP_RETCODE updateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             primalbound         /**< feasible primal bound */
   )
{
   SCIP_Real cutoffbound;

   /* increase the cutoff bound value by an epsilon to ensue that solution with the value of the cutoff bound are still
    * accepted
    */
   cutoffbound = primalbound + SCIPcutoffbounddelta(scip);

   if( cutoffbound < SCIPgetCutoffbound(scip) )
   {
      SCIPdebugMessage("update cutoff bound <%g>\n", cutoffbound);

      SCIP_CALL( SCIPupdateCutoffbound(scip, cutoffbound) );
   }
   else
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* we cannot disable the enforcement and propagation on ranged rows, because the cutoffbound could only have
       * resulted from one side
       */
      if( SCIPisInfinity(scip, -consdata->lhs) || SCIPisInfinity(scip, consdata->rhs) )
      {
         /* in case the cutoff bound is worse then the currently known one, we additionally avoid enforcement and
          * propagation
          */
         SCIP_CALL( SCIPsetConsEnforced(scip, cons, FALSE) );
         SCIP_CALL( SCIPsetConsPropagated(scip, cons, FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** check if the linear constraint is parallel to objective function; if so update the cutoff bound and avoid that the
 *  constraint enters the LP by setting the initial and separated flag to FALSE
 */
static
SCIP_RETCODE checkParallelObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real offset;
   SCIP_Real scale;
   SCIP_Bool applicable;
   int nobjvars;
   int nvars;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* ignore equalities since these are covert by the method checkPartialObjective() */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   nvars = consdata->nvars;
   nobjvars = SCIPgetNObjVars(scip);

   /* check if the linear inequality constraints has the same number of variables as the objective function and if the
    * initial and/or separated flag is set to FALSE
    */
   if( nvars != nobjvars || (!SCIPconsIsInitial(cons) && !SCIPconsIsSeparated(cons)) )
      return SCIP_OKAY;

   offset = 0.0;
   scale = 1.0;

   /* checks if the variables and their coefficients are equal (w.r.t. scaling factor) to the objective function */
   applicable = checkEqualObjective(scip, consdata, &scale, &offset);

   if( applicable )
   {
      /* avoid that the linear constraint enters the LP since it is parallel to the objective function */
      SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
      SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );

      if( SCIPisPositive(scip, scale) )
      {
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_Real primalbound;

            primalbound = (consdata->rhs - offset) / scale;

            SCIPdebugMessage("constraint <%s> is parallel to objective function and provides a cutoff bound <%g>\n",
               SCIPconsGetName(cons), primalbound);

            SCIP_CALL( updateCutoffbound(scip, cons, primalbound) );
         }

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_Real lowerbound;

            lowerbound = (consdata->lhs - offset) / scale;

            SCIPdebugMessage("constraint <%s> is parallel to objective function and provides a lower bound <%g>\n",
               SCIPconsGetName(cons), lowerbound);

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, lowerbound) );
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_Real lowerbound;

            lowerbound = (consdata->rhs - offset) / scale;

            SCIPdebugMessage("constraint <%s> is parallel to objective function and provides a lower bound <%g>\n",
               SCIPconsGetName(cons), lowerbound);

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, lowerbound) );
         }

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_Real primalbound;

            primalbound = (consdata->lhs - offset) / scale;

            SCIPdebugMessage("constraint <%s> is parallel to objective function and provides a cutoff bound <%g>\n",
               SCIPconsGetName(cons), primalbound);

            SCIP_CALL( updateCutoffbound(scip, cons, primalbound) );
         }
      }
   }

   return SCIP_OKAY;
}

/** converts special equalities */
static
SCIP_RETCODE convertEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->removedfixings);

   /* do nothing on inequalities */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   /* depending on the number of variables, call a special conversion method */
   if( consdata->nvars == 1 )
   {
      /* fix variable */
      SCIP_CALL( convertUnaryEquality(scip, cons, cutoff, nfixedvars, ndelconss) );
   }
   else if( consdata->nvars == 2 )
   {
      /* aggregate one of the variables */
      SCIP_CALL( convertBinaryEquality(scip, cons, cutoff, naggrvars, ndelconss) );
   }
   else
   {
      /* check if the equality is part of the objective function */
      SCIP_CALL( checkPartialObjective(scip, cons) );

      /* try to multi-aggregate one of the variables */
      SCIP_CALL( convertLongEquality(scip, cons, cutoff, naggrvars, ndelconss) );
   }

   return SCIP_OKAY;
}

/** returns whether the linear sum of all variables/coefficients except the given one divided by the given value is always
 *  integral
 */
static
SCIP_Bool consdataIsResidualIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   pos,                /**< position of variable to be left out */
   SCIP_Real             val                 /**< value to divide the coefficients by */
   )
{
   int v;
   
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != pos && (!SCIPvarIsIntegral(consdata->vars[v]) || !SCIPisIntegral(scip, consdata->vals[v]/val)) )
         return FALSE;
   }

   return TRUE;
}

/* check if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i 
 * check if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i 
 */
static
void calculateMinvalAndMaxval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             side,               /**< lhs or rhs */
   SCIP_Real             val,                /**< coefficient */
   SCIP_Real             minresactivity,     /**< minimal residual activity */
   SCIP_Real             maxresactivity,     /**< maximal residual activity */
   SCIP_Real*            minval,             /**< pointer to store calculated minval */
   SCIP_Real*            maxval              /**< pointer to store calculated maxval */
   )
{
   assert(scip != NULL);
   assert(minval != NULL);
   assert(maxval != NULL);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, ABS(maxresactivity)) )
         *minval = -maxresactivity;
      else
         *minval = (side - maxresactivity)/val;
      
      if( SCIPisInfinity(scip, ABS(minresactivity)) )
         *maxval = -minresactivity;
      else
         *maxval = (side - minresactivity)/val;
   }
   else
   {
      if( SCIPisInfinity(scip, ABS(minresactivity)) )
         *minval = minresactivity;
      else
         *minval = (side - minresactivity)/val;
      
      if( SCIPisInfinity(scip, ABS(maxresactivity)) )
         *maxval = maxresactivity;
      else
         *maxval = (side - maxresactivity)/val;
   }
}


/* applies dual presolving for variables that are locked only once in a direction, and this locking is due to a
 * linear inequality
 */
static
SCIP_RETCODE dualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsexists;
   SCIP_Bool rhsexists;
   SCIP_Bool bestisint;
   SCIP_Bool bestislhs;
   int bestpos;
   int i;
   int maxotherlocks;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   /* only process checked constraints (for which the locks are increased);
    * otherwise we would have to check for variables with nlocks == 0, and these are already processed by the
    * dualfix presolver
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;
      
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
   rhsexists = !SCIPisInfinity(scip, consdata->rhs);

   /* search for a single-locked variable which can be multi-aggregated; if a valid continuous variable was found, we
    * can use it safely for aggregation and break the search loop
    */
   bestpos = -1;
   bestisint = TRUE;
   bestislhs = FALSE;

   /* We only want to multi-aggregate variables, if they appear in maximal one additional constraint,
    * everything else would produce fill-in. Exceptions:
    * - If there are only two variables in the constraint from which the multi-aggregation arises, no fill-in will be
    *   produced.
    * - If there are three variables in the constraint, multi-aggregation in three additional constraints will remove
    *   six nonzeros (three from the constraint and the three entries of the multi-aggregated variable) and add
    *   six nonzeros (two variables per substitution).
    * - If there at most four variables in the constraint, multi-aggregation in two additional constraints will remove
    *   six nonzeros (four from the constraint and the two entries of the multi-aggregated variable) and add
    *   six nonzeros (three variables per substitution). God exists! 
    */
   if( consdata->nvars <= 2 )
      maxotherlocks = INT_MAX;
   else if( consdata->nvars == 3 )
      maxotherlocks = 3;
   else if( consdata->nvars == 4 )
      maxotherlocks = 2;
   else
      maxotherlocks = 1;

   /* if this constraint has both sides, it also provides a lock for the other side and thus we can allow one more lock */
   if( lhsexists && rhsexists )
      maxotherlocks++;

   for( i = 0; i < consdata->nvars && bestisint; ++i )
   {
      SCIP_VAR* var;
      SCIP_Bool isint;
      SCIP_Real val;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool agglhs;
      SCIP_Bool aggrhs;

      var = consdata->vars[i];
      isint = (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      /* if we already found a candidate, skip integers */
      if( bestpos >= 0 && isint )
         continue;

      /* better do not multi-aggregate binary variables, since most plugins rely on their binary variables to be either
       * active, fixed, or single-aggregated with another binary variable
       */
      if( SCIPvarIsBinary(var) && consdata->nvars > 2 )
         continue;

      if ( SCIPdoNotMultaggrVar(scip, var) )
         continue;

      val = consdata->vals[i];
      obj = SCIPvarGetObj(var);
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* lhs <= a_0 * x_0 + a_1 * x_1 + ... + a_{n-1} * x_{n-1} <= rhs
       *
       * a_i >= 0, c_i >= 0, lhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i <= 0, lhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i >= 0, c_i <= 0, rhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i >= 0, rhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *       
       * but: all this is only applicable, if the aggregated value is inside x_i's bounds for all possible values
       *      of all x_j
       * furthermore: we only want to apply this, if no fill-in will be produced
       */
      agglhs = lhsexists
         && ((val > 0.0 && !SCIPisNegative(scip, obj) && SCIPvarGetNLocksDown(var) == 1 
               && SCIPvarGetNLocksUp(var) <= maxotherlocks)
            || (val < 0.0 && !SCIPisPositive(scip, obj) && SCIPvarGetNLocksUp(var) == 1
               && SCIPvarGetNLocksDown(var) <= maxotherlocks));
      aggrhs = rhsexists
         && ((val > 0.0 && !SCIPisPositive(scip, obj) && SCIPvarGetNLocksUp(var) == 1 
               && SCIPvarGetNLocksDown(var) <= maxotherlocks)            
            || (val < 0.0 && !SCIPisNegative(scip, obj)  && SCIPvarGetNLocksDown(var) == 1 
               && SCIPvarGetNLocksUp(var) <= maxotherlocks));
      if( agglhs || aggrhs )
      {
         SCIP_Real minresactivity;
         SCIP_Real maxresactivity;
         SCIP_Real minval;
         SCIP_Real maxval;
         SCIP_Bool minisrelax;
         SCIP_Bool maxisrelax;
         SCIP_Bool isminsettoinfinity;
         SCIP_Bool ismaxsettoinfinity;

         /* calculate bounds for \sum_{j \neq i} a_j * x_j */
         consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
            &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);
         assert(SCIPisLE(scip, minresactivity, maxresactivity));

         /* We called consdataGetActivityResiduals() saying that we do not need a good relaxation,
          * so whenever we have a relaxed activity, it should be relaxed to +/- infinity.
          * This is needed, because we do not want to rely on relaxed finite resactivities.
          */
         assert((!minisrelax || isminsettoinfinity) && (!maxisrelax || ismaxsettoinfinity));

         if( agglhs )
         {
            /* check if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            calculateMinvalAndMaxval(scip, consdata->lhs, val, minresactivity, maxresactivity, &minval, &maxval);

            assert(SCIPisLE(scip, minval, maxval));
            if( (!SCIPisInfinity(scip, -minval) && SCIPisFeasGE(scip, minval, lb)) &&
               (!SCIPisInfinity(scip, maxval) && SCIPisFeasLE(scip, maxval, ub)) )
            {
               SCIP_Real oldmaxresactivity;
               SCIP_Real oldminresactivity;
               SCIP_Bool recalculated;

               recalculated = FALSE;
               oldmaxresactivity = maxresactivity;
               oldminresactivity = minresactivity;

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
                  recalculated = !SCIPisEQ(scip, oldminresactivity, minresactivity);
                  isminsettoinfinity = TRUE; /* here it means only that it was even calculated */
               }

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
                  recalculated = recalculated || !SCIPisEQ(scip, oldmaxresactivity, maxresactivity);
                  ismaxsettoinfinity = TRUE; /* here it means only that it was even calculated */
               }

               /* minresactivity or maxresactivity wasn't reliable so recalculate min- and maxval*/
               if( recalculated )
               {
                  assert(SCIPisLE(scip, minresactivity, maxresactivity));

                  /* check again if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
                  calculateMinvalAndMaxval(scip, consdata->lhs, val, minresactivity, maxresactivity, &minval, &maxval);

                  assert(SCIPisLE(scip, minval, maxval));
               }

               if( !recalculated || (SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub)) )
               {
                  /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                   * in the multi-aggregation
                   */
                  if( !isint || (SCIPisIntegral(scip, consdata->lhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
                  {
                     bestpos = i;
                     bestisint = isint;
                     bestislhs = TRUE;
                     continue; /* no need to also look at the right hand side */
                  }
               }
            }
         }

         if( aggrhs )
         {
            /* check if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            calculateMinvalAndMaxval(scip, consdata->rhs, val, minresactivity, maxresactivity, &minval, &maxval);

            assert(SCIPisLE(scip,minval,maxval));
            if( (!SCIPisInfinity(scip, -minval) && SCIPisFeasGE(scip, minval, lb)) &&
               (!SCIPisInfinity(scip, maxval) && SCIPisFeasLE(scip, maxval, ub)) )
            {
               SCIP_Real oldmaxresactivity;
               SCIP_Real oldminresactivity;
               SCIP_Bool recalculated;

               recalculated = FALSE;
               oldmaxresactivity = maxresactivity;
               oldminresactivity = minresactivity;

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
                  recalculated = !SCIPisEQ(scip, oldminresactivity, minresactivity);
               }

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
                  recalculated = recalculated || !SCIPisEQ(scip, oldmaxresactivity, maxresactivity);
               }

               /* minresactivity or maxresactivity wasn't reliable so recalculate min- and maxval*/
               if( recalculated )
               {
                  /* check again if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
                  calculateMinvalAndMaxval(scip, consdata->rhs, val, minresactivity, maxresactivity, &minval, &maxval);
                  assert(SCIPisLE(scip,minval,maxval));
               }

               if( !recalculated || (SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub)) )
               {
                  /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                   * in the multi-aggregation
                   */
                  if( !isint || (SCIPisIntegral(scip, consdata->rhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
                  {
                     bestpos = i;
                     bestisint = isint;
                     bestislhs = FALSE;
                  }
               }
            }
         }
      }
   }

   if( bestpos >= 0 )
   {
      SCIP_VAR** aggrvars;
      SCIP_Real* aggrcoefs;
      SCIP_Real aggrconst;
      SCIP_VAR* bestvar;
      SCIP_Real bestval;
      int naggrs;
      int j;
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;
      SCIP_Bool supinf;                      /* might the supremum of the multi-aggregation be infinite? */
      SCIP_Bool infinf;                      /* might the infimum of the multi-aggregation be infinite? */

      assert(!bestislhs || lhsexists);
      assert(bestislhs || rhsexists);

      bestvar = consdata->vars[bestpos];
      bestval = consdata->vals[bestpos];
      assert(bestisint ==
         (SCIPvarGetType(bestvar) == SCIP_VARTYPE_BINARY || SCIPvarGetType(bestvar) == SCIP_VARTYPE_INTEGER));

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, consdata->nvars-1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrcoefs, consdata->nvars-1) );

      /* set up the multi-aggregation */
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMessage("linear constraint <%s> (dual): multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(bestvar));
      naggrs = 0;
      supinf = FALSE;
      infinf = FALSE;

      for( j = 0; j < consdata->nvars; ++j )
      {
         if( j != bestpos )
         {
            aggrvars[naggrs] = consdata->vars[j];
            aggrcoefs[naggrs] = -consdata->vals[j]/consdata->vals[bestpos];
            SCIPdebugPrintf(" %+.15g<%s>", aggrcoefs[naggrs], SCIPvarGetName(aggrvars[naggrs]));
            if( bestisint )
            {
               /* coefficient must be integral: round it to exact integral value */
               assert(SCIPisIntegral(scip, aggrcoefs[naggrs]));
               aggrcoefs[naggrs] = SCIPfloor(scip, aggrcoefs[naggrs]+0.5);
            }

            if( SCIPisGT(scip, aggrcoefs[naggrs], 0.0) )
            {
               supinf = supinf || SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[j]));
               infinf = infinf || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[j]));
            }
            else if( SCIPisLT(scip, aggrcoefs[naggrs], 0.0) )
            {
               supinf = supinf || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[j]));
               infinf = infinf || SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[j]));
            }

            naggrs++;
         }
      }
      aggrconst = (bestislhs ? consdata->lhs/bestval : consdata->rhs/bestval);
      SCIPdebugPrintf(" %+.15g, bounds of <%s>: [%.15g,%.15g]\n", aggrconst, SCIPvarGetName(bestvar),
         SCIPvarGetLbGlobal(bestvar), SCIPvarGetUbGlobal(bestvar));
      assert(naggrs == consdata->nvars-1);

      /* right hand side must be integral: round it to exact integral value */
      if( bestisint )
      {
         assert(SCIPisIntegral(scip, aggrconst));
         aggrconst = SCIPfloor(scip, aggrconst+0.5);
      }

      aggregated = FALSE;
      infeasible = FALSE;

      /* perform the multi-aggregation */
      if( !supinf || !infinf )
      {
         SCIP_CALL( SCIPmultiaggregateVar(scip, bestvar, naggrs, aggrvars, aggrcoefs, aggrconst, &infeasible, &aggregated) );
      }
      else
      {
         /* If the infimum and the supremum of a multi-aggregation are both infinite, then the multi-aggregation might not be resolvable.
          * E.g., consider the equality z = x-y. If x and y are both fixed to +infinity, the value for z is not determined */
         SCIPdebugMessage("do not perform multi-aggregation: infimum and supremum are both infinite\n");
      }
      /* free temporary memory */
      SCIPfreeBufferArray(scip, &aggrcoefs);
      SCIPfreeBufferArray(scip, &aggrvars);
            
      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMessage("linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      /* delete the constraint, if the aggregation was successful */
      if( aggregated )
      {
         SCIP_CALL( SCIPdelCons(scip, cons) );
         
         if( !consdata->upgraded )
            (*ndelconss)++;
         (*naggrvars)++;
      }
      else
      {
         SCIPdebugMessage("aggregation non successful!\n");
      }
   }

   return SCIP_OKAY;
}

/** converts all variables with fixed domain into FIXED variables */
static
SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars          /**< pointer to count the total number of fixed variables */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VARSTATUS varstatus;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool fixed;
   SCIP_Bool infeasible;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars != NULL);
      var = consdata->vars[v];
      varstatus = SCIPvarGetStatus(var);

      if( varstatus != SCIP_VARSTATUS_FIXED )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         if( SCIPisEQ(scip, lb, ub) )
         {
            SCIP_Real fixval;

            fixval = SCIPselectSimpleValue(lb - 0.9 * SCIPepsilon(scip), ub + 0.9 * SCIPepsilon(scip), MAXDNOM);
            SCIPdebugMessage("converting variable <%s> with fixed bounds [%.15g,%.15g] into fixed variable fixed at %.15g\n",
               SCIPvarGetName(var), lb, ub, fixval);
            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMessage(" -> infeasible fixing\n");
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( fixed )
               (*nfixedvars)++;
         }
      }
   }

   SCIP_CALL( applyFixings(scip, cons, &infeasible) );

   if( infeasible )
   {
      SCIPdebugMessage(" -> infeasible fixing\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->removedfixings);

   return SCIP_OKAY;
}

#define BINWEIGHT  1
#define INTWEIGHT  4
#define CONTWEIGHT 8

/** gets weight for variable in a "weighted number of variables" sum */
static
int getVarWeight(
   SCIP_VAR*             var                 /**< variable to get weight for */
   )
{
   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      return BINWEIGHT;
   case SCIP_VARTYPE_INTEGER:
   case SCIP_VARTYPE_IMPLINT:
      return INTWEIGHT;
   case SCIP_VARTYPE_CONTINUOUS:
      return CONTWEIGHT;
   default:
      SCIPerrorMessage("invalid variable type\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** tries to aggregate variables in equations a^Tx = lhs
 *  in case there are at most two binary variables with an odd coefficient and all other
 *  variables are not continuous and have an even coefficient then:
 *  - exactly one odd binary variables 
 *    this binary variables y can be fixed to 0 if the lhs is even and to 1 if the lhs is odd    
 *     - lhs is odd ->  y = 1
 *     - lhs is even -> y = 0
 *  - exactly two odd binary variables 
 *    aggregate the two binary variables with odd coefficient 
 *     - lhs is odd -> exactly one of the variable has to be 1 -> var1 + var2 = 1
 *     - lhs is even -> both have to take the same value -> var1 - var2 = 0
 */
static
SCIP_RETCODE aggregateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( cons != NULL );
   
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* check if the linear constraint is an equation with integral right hand side */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) || !SCIPisIntegral(scip, consdata->lhs) )
      return SCIP_OKAY;
   
   /* try to fix and aggregated variables until nothing is possible anymore */
   do
   {
      int v;
      int nvars;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Bool lhsodd;
   
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      SCIP_Bool aggregated;
      SCIP_Bool redundant;

      SCIP_VAR* var1;
      SCIP_VAR* var2;
      int noddvars;

      success = FALSE;
    
      lhs = consdata->lhs;
      vars = consdata->vars;
      vals = consdata->vals;
      nvars = consdata->nvars;
    
      assert( !SCIPisInfinity(scip, ABS(lhs)) );
   
      var1 = NULL;
      var2 = NULL;
      noddvars = 0;
      
      /* search for binary variables with an odd coefficient */
      for( v = 0; v < nvars && noddvars < 3; ++v )
      {
         SCIP_Longint val;
         
         /* all coefficients and variables have to be integral */
         if( !SCIPisIntegral(scip, vals[v]) || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            return SCIP_OKAY;
         
         val = (SCIP_Longint)SCIPfeasFloor(scip, vals[v]);
         if( val % 2 != 0 )
         {
            /* the odd values have to belong to binary variables */
            if( !SCIPvarIsBinary(vars[v]) )
               return SCIP_OKAY;
            
            if( noddvars == 0 )
               var1 = vars[v];
            else
               var2 = vars[v];
            
            noddvars++;
         }
      }
    
      /* check lhs is odd or even */
      lhsodd = (((SCIP_Longint)SCIPfeasFloor(scip, lhs)) % 2 != 0);
      
      if( noddvars == 1 )
      {
         assert( var1 != NULL );
         
         SCIPdebugMessage("linear constraint <%s>: try fixing variable <%s> to <%g>\n", 
            SCIPconsGetName(cons), SCIPvarGetName(var1), lhsodd ? 1.0 : 0.0);
         
         SCIP_CALL( SCIPfixVar(scip, var1, lhsodd? 1.0 : 0.0, &infeasible, &fixed) );
         
         /* check for infeasibility of fixing */
         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible fixing\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         
         if( fixed )
         {
            SCIPdebugMessage(" -> feasible fixing\n");
            (*nfixedvars)++;
            success = TRUE;
         }
      }
      else if( noddvars == 2 )
      {
         assert( var1 != NULL );
         assert( var2 != NULL );
         
         /* aggregate the two variables with odd coefficient 
          * - lhs is odd -> exactly one of the variable has to be 1 -> var1 + var2 = 1
          * - lhs is even -> both have to take the same value -> var1 - var2 = 0
          */
         SCIPdebugMessage("linear constraint <%s>: try aggregation of variables <%s> and <%s>\n", 
            SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
         
         SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, lhsodd ? 1.0 : -1.0,
               lhsodd ? 1.0 : 0.0, &infeasible, &redundant, &aggregated) );
      
         /* check for infeasibility of aggregation */
         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible aggregation\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         
         /* count the aggregation */
         if( aggregated )
         {
            SCIPdebugMessage(" -> feasible aggregation\n");
            (*naggrvars)++;
            success = TRUE;
         }
      }
      
      if( success )
      {
         /* apply fixings and aggregation to successfully rerun this presolving step */
         SCIP_CALL( applyFixings(scip, cons, &infeasible) );

         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible fixing\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         
         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );
      }
   }
   while( success );
   
   return SCIP_OKAY;
}


/*  tries to simplify coefficients and delete variables in inequalities lhs <= a^Tx <= rhs
 *  in case there is only one binary variable with an odd coefficient, all other
 *  variables are not continuous and have an even coefficient, and only one of the left and right 
 *  hand-sides is odd, then:
 *  1. the right hand-side is odd and the left hand-side is even or minus infinity then:
 *    - if the odd coefficient is equal to 1, delete the variable and decrease rhs by 1
 *    - otherwise, decrease coefficient and rhs by 1
 *  2. the left hand-side is odd and the right hand-side is even or infinity then:
 *    - if the odd coefficient equal to -1, delete the variable and increase lhs by 1
 *    - otherwise, increase coefficient and lhs by 1
 *  Afterwards we us the normalize method to further simplify the inequality 
 */
static
SCIP_RETCODE simplifyInequalities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides           /**< pointer to store the amount of changed sides */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;
   int v;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real oddbinval;
   SCIP_Bool lhsodd;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real val;
   
   SCIP_VAR* oddbinvar;
   int noddvals = 0;
   int pos = 0;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nchgcoefs != NULL );
   assert( nchgsides != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* try to delete variables and simplify constraint */
   do
   {
      lhs = consdata->lhs;
      rhs = consdata->rhs;

      assert( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) );
      
      /* check if right and left hand-side are both integral */
      if( !SCIPisIntegral(scip, rhs) || !SCIPisIntegral(scip, lhs) )
         return SCIP_OKAY;

      /* check if both sides are not even or odd together */
      if( SCIPisIntegral(scip, rhs / 2.0) == SCIPisIntegral(scip, lhs / 2.0))    
         return SCIP_OKAY;

      /* in case the left hand side is minus infinity and the right hand side is even, then there nothing to simplify */
      if( SCIPisInfinity(scip, -lhs) && SCIPisIntegral(scip, rhs / 2.0))    
         return SCIP_OKAY;
      
      /* in case the right hand side is infinity and the left hand side is even, then there nothing to simplify */
      if( SCIPisInfinity(scip, rhs) && SCIPisIntegral(scip, lhs / 2.0))    
         return SCIP_OKAY;

      oddbinvar = NULL;
      oddbinval = 0;
      success = FALSE;
      noddvals = 0;
      vars = consdata->vars;
      vals = consdata->vals;
      nvars = consdata->nvars;
      
      /* search for binary variables with an odd coefficient */
      for( v = 0; v < nvars; ++v )
      {
         /* all coefficients have to be integral and all variables not of continuous type */
         if( !SCIPisIntegral(scip, vals[v]) || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            return SCIP_OKAY;
         
         /* check if the coefficient is odd */
         val = vals[v] / 2.0;
         if( !SCIPisIntegral(scip, val) )
         {
            /* the odd values have to belong to binary variables */
            if( !SCIPvarIsBinary(vars[v]) )
               return SCIP_OKAY;
            
            oddbinvar = vars[v];
            oddbinval = vals[v];
            pos = v;
            noddvals++;
            if (noddvals >= 2)
               return SCIP_OKAY;
         }
      }
      
      /* now we found exactly one binary variables with an odd coefficient and all other variables have even
       * coefficients and are not of continuous type; furthermore, only one side is odd */
      if( noddvals)
      {
         assert(noddvals == 1);
         assert( oddbinvar != NULL );
         lhsodd = FALSE;
         
         /* check if lhs is odd or even */
         if ( !SCIPisInfinity(scip, -lhs) )
         {
            lhsodd = !SCIPisIntegral(scip, lhs / 2.0);
         }
         
         if ( lhsodd )
         {
            oddbinval++;
            SCIP_CALL( chgLhs(scip, cons, lhs + 1) );
            SCIPdebugMessage("linear constraint <%s>: increasing coefficient for variable <%s> to <%g> and lhs to <%g>\n",
               SCIPconsGetName(cons), SCIPvarGetName(oddbinvar), oddbinval, lhs + 1);
         }
         else
         {
            oddbinval--;
            SCIP_CALL( chgRhs(scip, cons, rhs - 1) );
            SCIPdebugMessage("linear constraint <%s>: decreasing coefficient for variable <%s> to <%g> and rhs to <%g>\n", 
               SCIPconsGetName(cons), SCIPvarGetName(oddbinvar), oddbinval , rhs - 1);
         }
         
         if ( SCIPisZero(scip, oddbinval) )
         {
            SCIP_CALL( delCoefPos( scip, cons, pos ) );
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, pos, oddbinval) );
         }

         (*nchgcoefs)++;
         (*nchgsides)++;

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );
         
         /* if only one var is left abort loop */
         if( consdata->nvars == 1 )
            success = FALSE;
         else
            success = TRUE;
      }
   }
   while( success );
   
   return SCIP_OKAY;
}


/* tries to aggregate an (in)equality and an equality in order to decrease the number of variables in the (in)equality:
 *   cons0 := a * cons0 + b * cons1,
 * where a = val1[v] and b = -val0[v] for common variable v which removes most variable weight;
 * for numerical stability, we will only accept integral a and b;
 * the variable weight is a weighted sum over all included variables, where each binary variable weighs BINWEIGHT,
 * each integer or implicit integer variable weighs INTWEIGHT and each continuous variable weighs CONTWEIGHT
 */
static
SCIP_RETCODE aggregateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< (in)equality to modify */
   SCIP_CONS*            cons1,              /**< equality to use for aggregation of cons0 */
   int*                  commonidx0,         /**< array with indices of variables in cons0, that appear also in cons1 */
   int*                  commonidx1,         /**< array with indices of variables in cons1, that appear also in cons0 */
   int*                  diffidx0minus1,     /**< array with indices of variables in cons0, that don't appear in cons1 */
   int*                  diffidx1minus0,     /**< array with indices of variables in cons1, that don't appear in cons0 */
   int                   nvarscommon,        /**< number of variables, that appear in both constraints */
   int                   commonidxweight,    /**< variable weight sum of common variables */
   int                   diffidx0minus1weight, /**< variable weight sum of variables in cons0, that don't appear in cons1 */
   int                   diffidx1minus0weight, /**< variable weight sum of variables in cons1, that don't appear in cons0 */
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   SCIP_Bool*            aggregated          /**< pointer to store whether an aggregation was made */
   )
{
   SCIP_CONSDATA* consdata0;
   SCIP_CONSDATA* consdata1;
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real aggrcoef;
   SCIP_Real scalarsum;
   SCIP_Real bestscalarsum;
   SCIP_Bool betterscalarsum;
   SCIP_Bool commonvarlindependent;  /* indicates whether coefficient vector of common variables in linearly dependent */
   int varweight;
   int nvars;
   int bestvarweight;
   int bestnvars;
   int bestv;
   int v;
   int i;

   assert(scip != NULL);
   assert(cons0 != NULL);
   assert(cons1 != NULL);
   assert(commonidx0 != NULL);
   assert(commonidx1 != NULL);
   assert(diffidx0minus1 != NULL);
   assert(diffidx1minus0 != NULL);
   assert(nvarscommon >= 1);
   assert(commonidxweight >= nvarscommon);
   assert(nchgcoefs != NULL);
   assert(aggregated != NULL);

   assert(SCIPconsIsActive(cons0));
   assert(SCIPconsIsActive(cons1));

   SCIPdebugMessage("try aggregation of <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

   /* cons0 is an (in)equality */
   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   assert(SCIPisLE(scip, consdata0->lhs, consdata0->rhs));
   assert(diffidx0minus1weight >= consdata0->nvars - nvarscommon);

   /* cons1 is an equality */
   consdata1 = SCIPconsGetData(cons1);
   assert(consdata1 != NULL);
   assert(consdata1->nvars >= 1);
   assert(SCIPisEQ(scip, consdata1->lhs, consdata1->rhs));
   assert(diffidx1minus0weight >= consdata1->nvars - nvarscommon);

   *aggregated = FALSE;

   /* search for the best common variable such that
    *   val1[var] * consdata0 - val0[var] * consdata1
    * has least weighted number of variables
    */
   bestvarweight = commonidxweight + diffidx0minus1weight;
   bestnvars = consdata0->nvars;
   bestv = -1;
   bestscalarsum = 0.0;
   commonvarlindependent = TRUE;
   for( v = 0; v < nvarscommon; ++v )
   {
      assert(consdata0->vars[commonidx0[v]] == consdata1->vars[commonidx1[v]]);
      a = consdata1->vals[commonidx1[v]];
      b = -consdata0->vals[commonidx0[v]];

      /* only try aggregation, if coefficients are integral (numerical stability) */
      if( SCIPisIntegral(scip, a) && SCIPisIntegral(scip, b) )
      {
         /* count the number of variables in the potential new constraint  a * consdata0 + b * consdata1 */
         varweight = diffidx0minus1weight + diffidx1minus0weight;
         nvars = consdata0->nvars + consdata1->nvars - 2*nvarscommon;
         scalarsum = REALABS(a) + REALABS(b);
         betterscalarsum = (scalarsum < bestscalarsum);
         for( i = 0; i < nvarscommon
                 && (varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum)); ++i )
         {
            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            if( !SCIPisZero(scip, aggrcoef) )
            {
               varweight += getVarWeight(consdata0->vars[commonidx0[i]]);
               nvars++;
            }
         }
         if( varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum) )
         {
            bestv = v;
            bestvarweight = varweight;
            bestnvars = nvars;
            bestscalarsum = scalarsum;
         }
      }
      
      /* update commonvarlindependent flag, if still TRUE:
       * v's common coefficient in cons1 / v's common coefficient in cons0 should be constant, i.e., equal 0's common coefficient in cons1 / 0's common coefficient in cons0
       */
      if( commonvarlindependent && v > 0 )
         commonvarlindependent = SCIPisEQ(scip,
            consdata1->vals[commonidx1[v]] * consdata0->vals[commonidx0[0]],
            consdata1->vals[commonidx1[0]] * consdata0->vals[commonidx0[v]]);
   }

   /* if better aggregation was found, create new constraint and delete old one */
   if( bestv != -1 || commonvarlindependent )
   {
      SCIP_CONS* newcons;
      SCIP_CONSDATA* newconsdata;
      SCIP_VAR** newvars;
      SCIP_Real* newvals;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
      int newnvars;

      if( bestv != -1 )
      {
         /* choose multipliers such that the multiplier for the (in)equality cons0 is positive */
         if( consdata1->vals[commonidx1[bestv]] > 0.0 )
         {
            a = consdata1->vals[commonidx1[bestv]];
            b = -consdata0->vals[commonidx0[bestv]];
         }
         else
         {
            a = -consdata1->vals[commonidx1[bestv]];
            b = consdata0->vals[commonidx0[bestv]];
         }
         assert(SCIPisIntegral(scip, a));
         assert(SCIPisPositive(scip, a));
         assert(SCIPisIntegral(scip, b));
         assert(!SCIPisZero(scip, b));
      }
      else
      {
         assert(commonvarlindependent);
         if( consdata1->vals[commonidx1[0]] > 0.0 )
         {
            a =  consdata1->vals[commonidx1[0]];
            b = -consdata0->vals[commonidx0[0]];
         }
         else
         {
            a = -consdata1->vals[commonidx1[0]];
            b =  consdata0->vals[commonidx0[0]];
         }
         assert(SCIPisPositive(scip, a));
         assert(!SCIPisZero(scip, b));

         /* if a/b is integral, then we can easily choose integer multipliers */
         if( SCIPisIntegral(scip, a/b) )
         {
            if( a/b > 0 )
            {
               a /= b;
               b = 1.0;
            }
            else
            {
               a /= -b;
               b = -1.0;
            }
         }

         /* setup best* variables that were not setup above because we are in the commonvarlindependent case */
         bestvarweight = diffidx0minus1weight + diffidx1minus0weight;
         bestnvars = consdata0->nvars + consdata1->nvars - 2*nvarscommon;
      }

      SCIPdebugMessage("aggregate linear constraints <%s> := %.15g*<%s> + %.15g*<%s>  ->  nvars: %d -> %d, weight: %d -> %d\n",
         SCIPconsGetName(cons0), a, SCIPconsGetName(cons0), b, SCIPconsGetName(cons1),
         consdata0->nvars, bestnvars, commonidxweight + diffidx0minus1weight, bestvarweight);
      SCIPdebugPrintCons(scip, cons0, NULL);
      SCIPdebugPrintCons(scip, cons1, NULL);

      /* get temporary memory for creating the new linear constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &newvars, bestnvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newvals, bestnvars) );

      /* calculate the common coefficients, if we have not recognized linear dependency */
      newnvars = 0;
      if( !commonvarlindependent )
      {
         for( i = 0; i < nvarscommon; ++i )
         {
            assert(0 <= commonidx0[i] && commonidx0[i] < consdata0->nvars);
            assert(0 <= commonidx1[i] && commonidx1[i] < consdata1->nvars);

            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            if( !SCIPisZero(scip, aggrcoef) )
            {
               assert(newnvars < bestnvars);
               newvars[newnvars] = consdata0->vars[commonidx0[i]];
               newvals[newnvars] = aggrcoef;
               newnvars++;
            }
         }
      }
      else
      {
         /* if we recognized linear dependency of the common coefficients, then the aggregation coefficient should be 0.0 for every common variable */
#ifndef NDEBUG
         for( i = 0; i < nvarscommon; ++i )
         {
            assert(0 <= commonidx0[i] && commonidx0[i] < consdata0->nvars);
            assert(0 <= commonidx1[i] && commonidx1[i] < consdata1->nvars);

            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            assert(SCIPisZero(scip, aggrcoef));
         }
#endif
      }

      /* calculate the coefficients appearing in cons0 but not in cons1 */
      for( i = 0; i < consdata0->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx0minus1[i] && diffidx0minus1[i] < consdata0->nvars);

         aggrcoef = a * consdata0->vals[diffidx0minus1[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata0->vars[diffidx0minus1[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }

      /* calculate the coefficients appearing in cons1 but not in cons0 */
      for( i = 0; i < consdata1->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx1minus0[i] && diffidx1minus0[i] < consdata1->nvars);

         aggrcoef = b * consdata1->vals[diffidx1minus0[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata1->vars[diffidx1minus0[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }
      assert(newnvars == bestnvars);

      /* calculate the new left and right hand side of the (in)equality */
      assert(!SCIPisInfinity(scip, -consdata1->lhs));
      assert(!SCIPisInfinity(scip, consdata1->rhs));
      if( SCIPisInfinity(scip, -consdata0->lhs) )
         newlhs = -SCIPinfinity(scip);
      else
         newlhs = a * consdata0->lhs + b * consdata1->lhs;
      if( SCIPisInfinity(scip, consdata0->rhs) )
         newrhs = SCIPinfinity(scip);
      else
         newrhs = a * consdata0->rhs + b * consdata1->rhs;

      /* create the new linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons0), newnvars, newvars, newvals, newlhs, newrhs,
            SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0),
            SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0),
            SCIPconsIsLocal(cons0), SCIPconsIsModifiable(cons0),
            SCIPconsIsDynamic(cons0), SCIPconsIsRemovable(cons0), SCIPconsIsStickingAtNode(cons0)) );

      newconsdata = SCIPconsGetData(newcons);
      assert(newconsdata != NULL);

      /* copy the upgraded flag from the old cons0 to the new constraint */
      newconsdata->upgraded = consdata0->upgraded;
      newconsdata->donotupgrade = consdata0->donotupgrade;

      /* normalize the new constraint */
      SCIP_CALL( normalizeCons(scip, newcons) );

      /* check, if we really want to use the new constraint instead of the old one:
       * use the new one, if the maximum norm doesn't grow too much
       */
      if( consdataGetMaxAbsval(SCIPconsGetData(newcons)) <= maxaggrnormscale * consdataGetMaxAbsval(consdata0) )
      {
         SCIPdebugMessage(" -> aggregated to <%s>\n", SCIPconsGetName(newcons));
         SCIPdebugPrintCons(scip, newcons, NULL);

         /* update the statistics: we changed all coefficients */
         if( !consdata0->upgraded )
            (*nchgcoefs) += consdata0->nvars + consdata1->nvars - nvarscommon;
         *aggregated = TRUE;

         /* delete the old constraint, and add the new linear constraint to the problem */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
      }

      /* release the new constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &newvals);
      SCIPfreeBufferArray(scip, &newvars);
   }

   return SCIP_OKAY;
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyLinearcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables and the 
 * coefficients are either equal or negated
 */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqLinearcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Bool coefsequal;
   SCIP_Bool coefsnegated;
   int i;

   assert(key1 != NULL);
   assert(key2 != NULL);
   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);

   scip = (SCIP*)userptr;
   assert(scip != NULL);

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   coefsequal = TRUE;
   coefsnegated = TRUE;

   for( i = 0; i < consdata1->nvars && (coefsequal || coefsnegated); ++i )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 ||
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         coefsequal = FALSE;
         coefsnegated = FALSE;
         break;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0);

      /* tests if coefficients are either equal or negated */
      val1 = consdata1->vals[i];
      val2 = consdata2->vals[i];
      coefsequal = coefsequal && SCIPisEQ(scip, val1, val2);
      coefsnegated = coefsnegated && SCIPisEQ(scip, val1, -val2);
   }

   return (coefsequal || coefsnegated);
}

#define MULTIPLIER 2048
/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValLinearcons)
{
   SCIP_CONSDATA* consdata;
   SCIP_Real maxabsrealval;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
   int maxabsval;
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   assert(key != NULL);
   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   assert(consdata->sorted);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   maxabsrealval = consdataGetMaxAbsval(consdata);
   /* hash value depends on vectors of variable indices */
   if( maxabsrealval > (SCIP_Real) INT_MAX )
      maxabsval = 0;
   else if( maxabsrealval < 1.0 )
      maxabsval = (int) (MULTIPLIER * maxabsrealval);
   else
      maxabsval = (int) maxabsrealval;

   hashval = (consdata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx + maxabsval; /*lint !e701*/

   return hashval;
}

/** updates the flags of the first constraint according to the ones of the second constraint */
static
SCIP_RETCODE updateFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(scip != NULL);
   assert(blkmem != NULL);
   assert(conss != NULL);
   assert(firstchange != NULL);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = SCIPcalcHashtableSize(10*nconss);
   hashtablesize = MAX(hashtablesize, HASHSIZE_LINEARCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyLinearcons, hashKeyEqLinearcons, hashKeyValLinearcons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      /* sorts the constraint */
      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      SCIP_CALL( consdataSort(scip, consdata0) );
      assert(consdata0->sorted);

      /* get constraint from current hash table with same variables as cons0 and with coefficients either equal or negated
       * to the ones of cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));

      if( cons1 != NULL )
      {
         SCIP_CONS* consstay;
         SCIP_CONS* consdel;
         SCIP_CONSDATA* consdatastay;
         SCIP_CONSDATA* consdatadel;
         SCIP_CONSDATA* consdata1;

         SCIP_Real lhs;
         SCIP_Real rhs;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));

         /* constraint found: create a new constraint with same coefficients and best left and right hand side;
          * delete old constraints afterwards
          */
         consdata1 = SCIPconsGetData(cons1);

         assert(consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);

         assert(consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         if( SCIPisEQ(scip, consdata0->vals[0], consdata1->vals[0]) )
         {
            /* the coefficients of both constraints are equal */
            assert(consdata0->nvars < 2 || SCIPisEQ(scip, consdata0->vals[1], consdata1->vals[1]));
            SCIPdebugMessage("aggregate linear constraints <%s> and <%s> with equal coefficients into single ranged row\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);

            lhs = MAX(consdata1->lhs, consdata0->lhs);
            rhs = MIN(consdata1->rhs, consdata0->rhs);
         }
         else
         {
            /* the coefficients of both rows are negations */
            assert(SCIPisEQ(scip, consdata0->vals[0], -(consdata1->vals[0])));
            assert(consdata0->nvars < 2 || SCIPisEQ(scip, consdata0->vals[1], -(consdata1->vals[1])));
            SCIPdebugMessage("aggregate linear constraints <%s> and <%s> with negated coefficients into single ranged row\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugPrintCons(scip, cons1, NULL);

            lhs = MAX(consdata1->lhs, -consdata0->rhs);
            rhs = MIN(consdata1->rhs, -consdata0->lhs);
         }

         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMessage("aggregated linear constraint <%s> is infeasible\n", SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* ensure that lhs <= rhs holds without tolerances as we only allow such rows to enter the LP */
         if( lhs > rhs )
         {
            rhs = (lhs + rhs)/2;
            lhs = rhs;
         }

         /* check which constraint has to stay;
          * changes applied to an upgraded constraint will not be considered in the instance */
         if( consdata1->upgraded && !consdata0->upgraded )
         {
            consstay = cons0;
            consdatastay = consdata0;
            consdel = cons1;
            consdatadel = consdata1;

            /* exchange consdel with consstay in hashtable */
            SCIP_CALL( SCIPhashtableRemove(hashtable, (void*) consdel) );
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) consstay) );
         }
         else
         {
            consstay = cons1;
            consdatastay = consdata1;
            consdel = cons0;
            consdatadel = consdata0;
         }

         /* update lhs and rhs of consstay */
         SCIP_CALL( chgLhs(scip, consstay, lhs) );
         SCIP_CALL( chgRhs(scip, consstay, rhs) );

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, consstay, consdel) );

         /* delete consdel */
         assert(!consdatastay->upgraded || (consdatastay->upgraded && consdatadel->upgraded));
         SCIP_CALL( SCIPdelCons(scip, consdel) );
         if( !consdatadel->upgraded )
            (*ndelconss)++;

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdatastay->changed && SCIPconsGetPos(consstay) < *firstchange )
            *firstchange = SCIPconsGetPos(consstay);

         assert(SCIPconsIsActive(consstay));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides,          /**< pointer to count number of changed left/right hand sides */
   int*                  nchgcoefs           /**< pointer to count number of changed coefficients */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   int* commonidx0;
   int* commonidx1;
   int* diffidx0minus1;
   int* diffidx1minus0;
   SCIP_Longint possignature0;
   SCIP_Longint negsignature0;
   SCIP_Bool cons0changed;
   SCIP_Bool cons0isequality;
   int diffidx1minus0size;
   int c;
   SCIP_Real cons0lhs;
   SCIP_Real cons0rhs;
   SCIP_Bool cons0upgraded;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);
   assert(nchgcoefs != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(cons0 != NULL);
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   cons0isequality = SCIPisEQ(scip, consdata0->lhs, consdata0->rhs);

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata0) );

   /* calculate bit signatures of cons0 for potentially positive and negative coefficients */
   consdataCalcSignatures(consdata0);
   possignature0 = consdata0->possignature;
   negsignature0 = consdata0->negsignature;

   /* get temporary memory for indices of common variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx0, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx0minus1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx1minus0, consdata0->nvars) );
   diffidx1minus0size = consdata0->nvars;

   cons0lhs = consdata0->lhs;
   cons0rhs = consdata0->rhs;
   cons0upgraded = consdata0->upgraded;

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && conss[chkind] != NULL; ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Longint possignature1;
      SCIP_Longint negsignature1;
      SCIP_Bool cons0dominateslhs;
      SCIP_Bool cons1dominateslhs;
      SCIP_Bool cons0dominatesrhs;
      SCIP_Bool cons1dominatesrhs;
      SCIP_Bool cons1isequality;
      SCIP_Bool coefsequal;
      SCIP_Bool coefsnegated;
      SCIP_Bool tryaggregation;
      int nvarscommon;
      int nvars0minus1;
      int nvars1minus0;
      int commonidxweight;
      int diffidx0minus1weight;
      int diffidx1minus0weight;
      int v0;
      int v1;

      assert(cons0lhs == consdata0->lhs);  /*lint !e777*/
      assert(cons0rhs == consdata0->rhs);  /*lint !e777*/
      assert(cons0upgraded == consdata0->upgraded);

      cons1 = conss[c];

      /* cons1 has become inactive during presolving of constraint pairs */
      if( cons1 == NULL )
         continue;
      
      assert(SCIPconsIsActive(cons0) && !SCIPconsIsModifiable(cons0));
      assert(SCIPconsIsActive(cons1) && !SCIPconsIsModifiable(cons1));

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* SCIPdebugMessage("preprocess linear constraint pair <%s>[chgd:%d, upgd:%d] and <%s>[chgd:%d, upgd:%d]\n",
         SCIPconsGetName(cons0), cons0changed, cons0upgraded,
         SCIPconsGetName(cons1), consdata1->changed, consdata1->upgraded); */

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( !cons0changed && !consdata1->changed )
         continue;
      
      /* if both constraints are already upgraded, skip the pair; 
       * because changes on these constraints cannot be applied to the instance anymore */
      if( cons0upgraded && consdata1->upgraded )
         continue;

      assert(consdata1->nvars >= 1);

      /* sort the constraint */
      SCIP_CALL( consdataSort(scip, consdata1) );

      /* calculate bit signatures of cons1 for potentially positive and negative coefficients */
      consdataCalcSignatures(consdata1);
      possignature1 = consdata1->possignature;
      negsignature1 = consdata1->negsignature;

      /* the signatures give a quick test to check for domination and equality of coefficients */
      coefsequal = (possignature0 == possignature1) && (negsignature0 == negsignature1);
      coefsnegated = (possignature0 == negsignature1) && (negsignature0 == possignature1);
      cons0dominateslhs = SCIPisGE(scip, cons0lhs, consdata1->lhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1dominateslhs = SCIPisGE(scip, consdata1->lhs, cons0lhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons0dominatesrhs = SCIPisLE(scip, cons0rhs, consdata1->rhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons1dominatesrhs = SCIPisLE(scip, consdata1->rhs, cons0rhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1isequality = SCIPisEQ(scip, consdata1->lhs, consdata1->rhs);
      tryaggregation = (cons0isequality || cons1isequality) && (maxaggrnormscale > 0.0);
      if( !cons0dominateslhs && !cons1dominateslhs && !cons0dominatesrhs && !cons1dominatesrhs
         && !coefsequal && !coefsnegated && !tryaggregation )
         continue;

      /* make sure, we have enough memory for the index set of V_1 \ V_0 */
      if( tryaggregation && consdata1->nvars > diffidx1minus0size )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &diffidx1minus0, consdata1->nvars) );
         diffidx1minus0size = consdata1->nvars;
      }

      /* check consdata0 against consdata1:
       * - if lhs0 >= lhs1 and for each variable v and each solution value x_v val0[v]*x_v <= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. left hand side
       * - if rhs0 <= rhs1 and for each variable v and each solution value x_v val0[v]*x_v >= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. right hand side
       * - if val0[v] == -val1[v] for all variables v, the two inequalities can be replaced by a single
       *   ranged row (or equality)
       * - if at least one constraint is an equality, count the weighted number of common variables W_c
       *   and the weighted number of variable in the difference sets W_0 = w(V_0 \ V_1), W_1 = w(V_1 \ V_0),
       *   where the weight of each variable depends on its type, such that aggregations in order to remove the
       *   number of continuous and integer variables are preferred:
       *   - if W_c > W_1, try to aggregate  consdata0 := a * consdata0 + b * consdata1  in order to decrease the
       *     variable weight in consdata0, where a = +/- val1[v] and b = -/+ val0[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons0
       *   - if W_c > W_0, try to aggregate  consdata1 := a * consdata1 + b * consdata0  in order to decrease the
       *     variable weight in consdata1, where a = +/- val0[v] and b = -/+ val1[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons1
       */

      /* check consdata0 against consdata1 for redundancy, or ranged row accumulation */
      nvarscommon = 0;
      commonidxweight = 0;
      nvars0minus1 = 0;
      diffidx0minus1weight = 0;
      nvars1minus0 = 0;
      diffidx1minus0weight = 0;
      v0 = 0;
      v1 = 0;
      while( (v0 < consdata0->nvars || v1 < consdata1->nvars)
         && (cons0dominateslhs || cons1dominateslhs || cons0dominatesrhs || cons1dominatesrhs
            || coefsequal || coefsnegated || tryaggregation) )
      {
         SCIP_VAR* var;
         SCIP_Real val0;
         SCIP_Real val1;
         int varcmp;

         /* test, if variable appears in only one or in both constraints */
         if( v0 < consdata0->nvars && v1 < consdata1->nvars )
            varcmp = SCIPvarCompare(consdata0->vars[v0], consdata1->vars[v1]);
         else if( v0 < consdata0->nvars )
            varcmp = -1;
         else
            varcmp = +1;

         switch( varcmp )
         {
         case -1:
            /* variable doesn't appear in consdata1 */
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = 0.0;
            if( tryaggregation )
            {
               diffidx0minus1[nvars0minus1] = v0;
               nvars0minus1++;
               diffidx0minus1weight += getVarWeight(var);
            }
            v0++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case +1:
            /* variable doesn't appear in consdata0 */
            var = consdata1->vars[v1];
            val0 = 0.0;
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               diffidx1minus0[nvars1minus0] = v1;
               nvars1minus0++;
               diffidx1minus0weight += getVarWeight(var);
            }
            v1++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(consdata0->vars[v0] == consdata1->vars[v1]);
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               commonidx0[nvarscommon] = v0;
               commonidx1[nvarscommon] = v1;
               nvarscommon++;
               commonidxweight += getVarWeight(var);
            }
            v0++;
            v1++;
            coefsequal = coefsequal && (SCIPisEQ(scip, val0, val1));
            coefsnegated = coefsnegated && (SCIPisEQ(scip, val0, -val1));
            break;

         default:
            SCIPerrorMessage("invalid comparison result\n");
            var = NULL;
            val0 = 0.0;
            val1 = 0.0;
            SCIPABORT();
         }
         assert(var != NULL);

         /* update domination criteria w.r.t. the coefficient and the variable's bounds */
         if( SCIPisGT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
         }
         else if( SCIPisLT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
         }
      }

      /* check for disaggregated ranged rows */
      if( coefsequal || coefsnegated )
      {
         SCIP_CONS* consstay;
         SCIP_CONS* consdel;
#ifndef NDEBUG
         SCIP_CONSDATA* consdatastay;
#endif
         SCIP_CONSDATA* consdatadel;
         SCIP_Real lhs;
         SCIP_Real rhs;
         int consinddel;

         /* the coefficients in both rows are either equal or negated: create a new constraint with same coefficients and
          * best left and right hand sides; delete the old constraints afterwards
          */
         SCIPdebugMessage("aggregate linear constraints <%s> and <%s> with %s coefficients into single ranged row\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1), coefsequal ? "equal" : "negated");
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         if( coefsequal )
         {
            /* the coefficients of both rows are equal */
            lhs = MAX(consdata0->lhs, consdata1->lhs);
            rhs = MIN(consdata0->rhs, consdata1->rhs);
         }
         else
         {
            /* the coefficients of both rows are negations */
            lhs = MAX(consdata0->lhs, -consdata1->rhs);
            rhs = MIN(consdata0->rhs, -consdata1->lhs);
         }
         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMessage("aggregated linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *cutoff = TRUE;
            break;
         }

         /* check which constraint has to stay; 
          * changes applied to an upgraded constraint will not be considered in the instance */
         if( consdata0->upgraded )
         {
            assert(!consdata1->upgraded);
            consstay = cons1;
#ifndef NDEBUG
            consdatastay = consdata1;
#endif

            consdel = cons0;
            consdatadel = consdata0;
            consinddel = chkind;
         }
         else
         {
            consstay = cons0;
#ifndef NDEBUG
            consdatastay = consdata0;
#endif

            consdel = cons1;
            consdatadel = consdata1;
            consinddel = c;
         }

         /* update the sides of consstay */
         SCIP_CALL( chgLhs(scip, consstay, lhs) );
         SCIP_CALL( chgRhs(scip, consstay, rhs) );
         if( !consdata0->upgraded )
         {
            assert(consstay == cons0);
            cons0lhs = consdata0->lhs;
            cons0rhs = consdata0->rhs;
         }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, consstay, consdel) ); 

         assert( !consdatastay->upgraded );
         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, consdel) );
         conss[consinddel] = NULL;
         if( !consdatadel->upgraded )
            (*ndelconss)++;
         continue;
      }

      /* check for domination: remove dominated sides, but don't touch equalities as long as they are not totally
       * redundant
       */
      if( cons1dominateslhs && (!cons0isequality || cons1dominatesrhs || SCIPisInfinity(scip, consdata0->rhs) ) )
      {
         /* left hand side is dominated by consdata1: delete left hand side of consdata0 */
         SCIPdebugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata1->lhs, consdata0->rhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata0->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons0, -SCIPinfinity(scip)) );
            cons0lhs = consdata0->lhs;
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons1, cons0) ); 

               (*nchgsides)++;
            }
         }
      }
      else if( cons0dominateslhs && (!cons1isequality || cons0dominatesrhs || SCIPisInfinity(scip, consdata1->rhs)) )
      {
         /* left hand side is dominated by consdata0: delete left hand side of consdata1 */
         SCIPdebugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata0->lhs, consdata1->rhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata1->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons0, cons1) ); 

               (*nchgsides)++;
            }
         }
      }
      if( cons1dominatesrhs && (!cons0isequality || cons1dominateslhs || SCIPisInfinity(scip, -consdata0->lhs)) )
      {
         /* right hand side is dominated by consdata1: delete right hand side of consdata0 */
         SCIPdebugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata1->rhs, consdata0->lhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata0->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons0, SCIPinfinity(scip)) );
            cons0rhs = consdata0->rhs;
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons1, cons0) ); 

               (*nchgsides)++;
            }
         }
      }
      else if( cons0dominatesrhs && (!cons1isequality || cons0dominateslhs || SCIPisInfinity(scip, -consdata1->lhs)) )
      {
         /* right hand side is dominated by consdata0: delete right hand side of consdata1 */
         SCIPdebugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata0->rhs, consdata1->lhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata1->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons0, cons1) ); 

               (*nchgsides)++;
            }
         }
      }

      /* check for now redundant constraints */
      if( SCIPisInfinity(scip, -consdata0->lhs) && SCIPisInfinity(scip, consdata0->rhs) )
      {
         /* consdata0 became redundant */
         SCIPdebugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         conss[chkind] = NULL;
         if( !consdata0->upgraded )
         {
            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( updateFlags(scip, cons1, cons0) ); 
            
            (*ndelconss)++;
         }
         continue;
      }
      if( SCIPisInfinity(scip, -consdata1->lhs) && SCIPisInfinity(scip, consdata1->rhs) )
      {
         /* consdata1 became redundant */
         SCIPdebugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         SCIP_CALL( SCIPdelCons(scip, cons1) );
         conss[c] = NULL;
         if( !consdata1->upgraded )
         {
            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( updateFlags(scip, cons0, cons1) ); 

            (*ndelconss)++;
         }
         continue;
      }

      /* check, if we want to aggregate an (in)equality with an equality:
       *   consdata0 := a * consdata0 + b * consdata1  or  consdata1 := a * consdata1 + b * consdata0
       */
      if( tryaggregation )
      {
         SCIP_Bool aggregated;

         assert(consdata0->nvars == nvarscommon + nvars0minus1);
         assert(consdata1->nvars == nvarscommon + nvars1minus0);
         
         aggregated = FALSE;
         if( cons1isequality && !consdata0->upgraded && commonidxweight > diffidx1minus0weight )
         {
            /* W_c > W_1: try to aggregate  consdata0 := a * consdata0 + b * consdata1 */
            SCIP_CALL( aggregateConstraints(scip, cons0, cons1, commonidx0, commonidx1, diffidx0minus1, diffidx1minus0,
                  nvarscommon, commonidxweight, diffidx0minus1weight, diffidx1minus0weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );

            /* update array of active constraints */
            if( aggregated )
            {
               assert(!SCIPconsIsActive(cons0));
               assert(SCIPconsIsActive(cons1));
               conss[chkind] = NULL;
            }
         }
         if( !aggregated && cons0isequality && !consdata1->upgraded && commonidxweight > diffidx0minus1weight )
         {
            /* W_c > W_0: try to aggregate  consdata1 := a * consdata1 + b * consdata0 */
            SCIP_CALL( aggregateConstraints(scip, cons1, cons0, commonidx1, commonidx0, diffidx1minus0, diffidx0minus1,
                  nvarscommon, commonidxweight, diffidx1minus0weight, diffidx0minus1weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );

            /* update array of active constraints */
            if( aggregated )
            {
               assert(!SCIPconsIsActive(cons1));
               assert(SCIPconsIsActive(cons0));
               conss[c] = NULL;
            }
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &diffidx1minus0);
   SCIPfreeBufferArray(scip, &diffidx0minus1);
   SCIPfreeBufferArray(scip, &commonidx1);
   SCIPfreeBufferArray(scip, &commonidx0);

   return SCIP_OKAY;
}

/** applies full dual presolving on variables that only appear in linear constraints */
static
SCIP_RETCODE fullDualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds             /**< pointer to count the number of bound changes */
   )
{
   SCIP_Real* redlb;
   SCIP_Real* redub;
   int* nlocksdown;
   int* nlocksup;
   SCIP_Bool* isimplint;
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   SCIP_VAR** conscontvars;
   int nvars;
   int nbinvars;
   int nintvars;
   int ncontvars;
   int v;
   int c;

   /* we calculate redundancy bounds with the following meaning:
    *   redlb[v] == k : if x_v >= k, we can always round x_v down to x_v == k without violating any constraint
    *   redub[v] == k : if x_v <= k, we can always round x_v up to x_v == k without violating any constraint
    * then:
    *   c_v >= 0 : x_v <= redlb[v] is feasible due to optimality
    *   c_v <= 0 : x_v >= redub[v] is feasible due to optimality
    */

   /* Additionally, we detect continuous variables that are implicitly integral.
    * A continuous variable j is implicit integral if it only has only +/-1 coefficients,
    * and all constraints (including the bounds as trivial constraints) in which:
    *   c_j > 0: the variable is down-locked,
    *   c_j < 0: the variable is up-locked,
    *   c_j = 0: the variable appears
    * have, apart from j, only integer variables with integral coefficients and integral sides.
    * This is because then, the value of the variable is either determined by one of its bounds or
    * by one of these constraints, and in all cases, the value of the variable is integral.
    */

   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(nchgbds != NULL);
   assert(!SCIPinProbing(scip));

   /* get active variables */
   nvars = SCIPgetNVars(scip);
   origvars = SCIPgetVars(scip);

   /* if the problem is a pure binary program, nothing can be achieved by full dual presolve */
   nbinvars = SCIPgetNBinVars(scip);
   if( nbinvars == nvars )
      return SCIP_OKAY;

   /* get number of continuous variables */
   ncontvars = SCIPgetNContVars(scip);
   nintvars = nvars - ncontvars;

   /* copy the variable array since this array might change during the curse of this algorithm */
   nvars = nvars - nbinvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &(origvars[nbinvars]), nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &redlb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksdown, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksup, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isimplint, ncontvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conscontvars, ncontvars) );

   /* initialize redundancy bounds */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY);
      redlb[v] = SCIPvarGetLbGlobal(vars[v]);
      redub[v] = SCIPvarGetUbGlobal(vars[v]);
   }
   BMSclearMemoryArray(nlocksdown, nvars);
   BMSclearMemoryArray(nlocksup, nvars);

   /* Initialize isimplint array: variable may be implied integer if both bounds are integral.
    * We better not use SCIPisFeasIntegral() in these checks.
    */
   for( v = 0; v < ncontvars; v++ )
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      
      var = vars[v + nintvars - nbinvars];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
      isimplint[v] = (SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb))
         && (SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub));
   }

   /* scan all constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* we only need to consider constraints that have been locked (i.e., checked constraints or constraints that are
       * part of checked disjunctions)
       */
      if( SCIPconsIsLocked(conss[c]) )
      {
         SCIP_CONSDATA* consdata;
         SCIP_Bool lhsexists;
         SCIP_Bool rhsexists;
         SCIP_Bool hasimpliedpotential;
         SCIP_Bool integralcoefs;
         int nlockspos;
         int contvarpos;
         int nconscontvars;
         int i;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* get number of times the constraint was locked */
         nlockspos = SCIPconsGetNLocksPos(conss[c]);

         /* we do not want to include constraints with locked negation (this would be too weird) */
         if( SCIPconsGetNLocksNeg(conss[c]) > 0 )
         {
            /* mark all continuous variables as not being implicit integral */
            for( i = 0; i < consdata->nvars; ++i )
            {
               SCIP_VAR* var;

               var = consdata->vars[i];
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  int contv;
                  contv = SCIPvarGetProbindex(var) - nintvars;
                  assert(0 <= contv && contv < ncontvars); /* variable should be active due to applyFixings() */
                  isimplint[contv] = FALSE;
               }
            }
            continue;
         }

         /* check for existing sides */
         lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
         rhsexists = !SCIPisInfinity(scip, consdata->rhs);

         /* count locks and update redundancy bounds */
         contvarpos = -1;
         nconscontvars = 0;
         hasimpliedpotential = FALSE;
         integralcoefs = !SCIPconsIsModifiable(conss[c]);
         
         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_VAR* var;
            SCIP_Real val;
            SCIP_Real minresactivity;
            SCIP_Real maxresactivity;
            SCIP_Real newredlb;
            SCIP_Real newredub;
            SCIP_Bool minisrelax;
            SCIP_Bool maxisrelax;
            SCIP_Bool isminsettoinfinity;
            SCIP_Bool ismaxsettoinfinity;
            int arrayindex;

            var = consdata->vars[i];
            val = consdata->vals[i];

            /* check if still all integer variables have integral coefficients */
            if( SCIPvarIsIntegral(var) )
               integralcoefs = integralcoefs && SCIPisIntegral(scip, val);

            /* we do not need to process binary variables */
            if( SCIPvarIsBinary(var) )
               continue;

            if( SCIPconsIsModifiable(conss[c]) )
            {
               minresactivity = -SCIPinfinity(scip);
               maxresactivity =  SCIPinfinity(scip);
               isminsettoinfinity = TRUE;
               ismaxsettoinfinity = TRUE;
            }
            else
            {
               /* calculate residual activity bounds if variable would be fixed to zero */
               consdataGetGlbActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
                  &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);

               /* We called consdataGetGlbActivityResiduals() saying that we do not need a good relaxation,
                * so whenever we have a relaxed activity, it should be relaxed to +/- infinity.
                * This is needed, because we do not want to rely on relaxed finite resactivities.
                */
               assert((!minisrelax || isminsettoinfinity) && (!maxisrelax || ismaxsettoinfinity));

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastglbminactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, TRUE);

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastglbmaxactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, TRUE);
            }

            arrayindex = SCIPvarGetProbindex(var) - nbinvars;

            assert(0 <= arrayindex && arrayindex < nvars); /* variable should be active due to applyFixings() */

            newredlb = redlb[arrayindex];
            newredub = redub[arrayindex];
            if( val > 0.0 )
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d > 0  ->  redundant in y if  x >= (lhs - min{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (isminsettoinfinity ? SCIPinfinity(scip) : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d > 0  ->  redundant in y if  x <= (rhs - max{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (ismaxsettoinfinity ? -SCIPinfinity(scip) : (consdata->rhs - maxresactivity)/val);
               }
            }
            else
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d < 0  ->  redundant in y if  x <= (lhs - min{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (isminsettoinfinity ? -SCIPinfinity(scip) : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d < 0  ->  redundant in y if  x >= (rhs - max{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (ismaxsettoinfinity ? SCIPinfinity(scip) : (consdata->rhs - maxresactivity)/val);
               }
            }

            /* if the variable is integer, we have to round the value to the next integral value */
            if( SCIPvarIsIntegral(var) )
            {
               if( !SCIPisInfinity(scip, newredlb) )
                  newredlb = SCIPceil(scip, newredlb);
               if( !SCIPisInfinity(scip, -newredub) )
                  newredub = SCIPfloor(scip, newredub);
            }

            /* update redundancy bounds */
            redlb[arrayindex] = MAX(redlb[arrayindex], newredlb);
            redub[arrayindex] = MIN(redub[arrayindex], newredub);

            /* collect the continuous variables of the constraint */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               int contv;

               assert(nconscontvars < ncontvars);
               contvarpos = i;
               conscontvars[nconscontvars] = var;
               nconscontvars++;

               contv = SCIPvarGetProbindex(var) - nintvars;
               assert(0 <= contv && contv < ncontvars);
               hasimpliedpotential = hasimpliedpotential || isimplint[contv];
            }
         }

         /* update implied integer status of continuous variables */
         if( hasimpliedpotential )
         {
            if( nconscontvars > 1 || !integralcoefs )
            {
               /* there is more than one continuous variable or the integer variables have fractional coefficients:
                * none of the continuous variables is implied integer
                */
               for( i = 0; i < nconscontvars; i++ )
               {
                  int contv;
                  contv = SCIPvarGetProbindex(conscontvars[i]) - nintvars;
                  assert(0 <= contv && contv < ncontvars);
                  isimplint[contv] = FALSE;
               }
            }
            else
            {
               SCIP_VAR* var;
               SCIP_Real val;
               SCIP_Real absval;
               int contv;

               /* there is exactly one continuous variable and the integer variables have integral coefficients:
                * this is the interesting case, and we have to check whether the coefficient is +/-1 and the corresponding
                * side(s) of the constraint is integral
                */
               assert(nconscontvars == 1);
               assert(0 <= contvarpos && contvarpos < consdata->nvars);
               var = consdata->vars[contvarpos];
               val = consdata->vals[contvarpos];
               contv = SCIPvarGetProbindex(var) - nintvars;
               assert(0 <= contv && contv < ncontvars);
               assert(isimplint[contv]);
               
               absval = REALABS(val);
               if( !SCIPisEQ(scip, absval, 1.0) )
                  isimplint[contv] =  FALSE;
               else
               {
                  SCIP_Real obj;
                  
                  obj = SCIPvarGetObj(var);
                  if( obj * val >= 0.0 && lhsexists )
                  {
                     /* the variable may be blocked by the constraint's left hand side */
                     isimplint[contv] = isimplint[contv] && SCIPisIntegral(scip, consdata->lhs);
                  }
                  if( obj * val <= 0.0 && rhsexists )
                  {
                     /* the variable may be blocked by the constraint's left hand side */
                     isimplint[contv] = isimplint[contv] && SCIPisIntegral(scip, consdata->rhs);
                  }
               }
            }
         }
      }
   }

   /* check if any bounds can be tightened due to optimality */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY);
      assert(SCIPvarGetNLocksDown(vars[v]) >= nlocksdown[v]);
      assert(SCIPvarGetNLocksUp(vars[v]) >= nlocksup[v]);

      var = vars[v];
      obj = SCIPvarGetObj(var);
      if( obj >= 0.0 )
      {
         /* making the variable as small as possible does not increase the objective:
          * check if all down locks of the variables are due to linear constraints;
          * if largest bound to make constraints redundant is -infinity, we better do nothing for numerical reasons
          */
         if( SCIPvarGetNLocksDown(var) == nlocksdown[v]
            && !SCIPisInfinity(scip, -redlb[v])
            && redlb[v] < SCIPvarGetUbGlobal(var) )
         {
            SCIP_Real ub;

            /* if x_v >= redlb[v], we can always round x_v down to x_v == redlb[v] without violating any constraint
             *  -> tighten upper bound to x_v <= redlb[v]
             */
            SCIPdebugMessage("variable <%s> only locked down in linear constraints: dual presolve <%s>[%.15g,%.15g] <= %.15g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redlb[v]);
            SCIP_CALL( SCIPtightenVarUb(scip, var, redlb[v], FALSE, &infeasible, &tightened) );
            assert(!infeasible);

            ub = SCIPvarGetUbGlobal(var);
            redub[v] = MIN(redub[v], ub);
            if( tightened )
               (*nchgbds)++;
         }
      }
      if( obj <= 0.0 )
      {
         /* making the variable as large as possible does not increase the objective:
          * check if all up locks of the variables are due to linear constraints;
          * if smallest bound to make constraints redundant is +infinity, we better do nothing for numerical reasons
          */
         if( SCIPvarGetNLocksUp(var) == nlocksup[v]
            && !SCIPisInfinity(scip, redub[v])
            && redub[v] > SCIPvarGetLbGlobal(var) )
         {
            SCIP_Real lb;

            /* if x_v <= redub[v], we can always round x_v up to x_v == redub[v] without violating any constraint
             *  -> tighten lower bound to x_v >= redub[v]
             */
            SCIPdebugMessage("variable <%s> only locked up in linear constraints: dual presolve <%s>[%.15g,%.15g] >= %.15g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redub[v]);
            SCIP_CALL( SCIPtightenVarLb(scip, var, redub[v], FALSE, &infeasible, &tightened) );
            assert(!infeasible);

            lb = SCIPvarGetLbGlobal(var);
            redlb[v] = MAX(redlb[v], lb);
            if( tightened )
               (*nchgbds)++;
         }
      }
   }

   /* upgrade continuous variables to implied integers */
   for( v = nintvars - nbinvars; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool infeasible;

      var = vars[v];
      assert(var != NULL);

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
      assert(SCIPvarGetNLocksDown(var) >= nlocksdown[v]);
      assert(SCIPvarGetNLocksUp(var) >= nlocksup[v]);
      assert(0 <= v - nintvars + nbinvars && v - nintvars + nbinvars < ncontvars);

      /* we can only conclude implied integrality if the variable appears in no other constraint */
      if( isimplint[v - nintvars + nbinvars]
         && SCIPvarGetNLocksDown(var) == nlocksdown[v]
         && SCIPvarGetNLocksUp(var) == nlocksup[v] )
      {

         /* since we locally copied the variable array we can change the variable type immediately */
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );

         if( infeasible )
         {
            SCIPdebugMessage("infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
            *cutoff = TRUE;

            break;
         }
         
         SCIPdebugMessage("dual presolve: converting continuous variable <%s>[%g,%g] to implicit integer\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &conscontvars);
   SCIPfreeBufferArray(scip, &isimplint);
   SCIPfreeBufferArray(scip, &nlocksup);
   SCIPfreeBufferArray(scip, &nlocksdown);
   SCIPfreeBufferArray(scip, &redub);
   SCIPfreeBufferArray(scip, &redlb);

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLinear)
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
SCIP_DECL_CONSINIT(consInitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   assert(nconss == 0 || conss != NULL);

   /* catch events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->eventdatas == NULL);

      /* catch all events */
      SCIP_CALL( consdataCatchAllEvents(scip, consdata, conshdlrdata->eventhdlr) );
      assert(consdata->eventdatas != NULL);
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* drop events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->eventdatas != NULL )
      {
         /* drop all events */
         SCIP_CALL( consdataDropAllEvents(scip, consdata, conshdlrdata->eventhdlr) );
         assert(consdata->eventdatas == NULL);
      }
   }

   return SCIP_OKAY;

}

#ifdef WITH_PRINTORIGCONSTYPES

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isRangedRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return !(SCIPisEQ(scip, SCIPconsGetData(cons)->lhs, SCIPconsGetData(cons)->rhs)
      || SCIPisInfinity(scip, -SCIPconsGetData(cons)->lhs) || SCIPisInfinity(scip, SCIPconsGetData(cons)->rhs) );
}

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isFiniteNonnegativeIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x                   /**< value */
   )
{
   assert(scip != NULL);

   return (!SCIPisInfinity(scip, x) && !SCIPisNegative(scip, x) && SCIPisIntegral(scip, x));
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreLinear)
{  /*lint --e{715}*/
   int counter[SCIP_CONSTYPE_GENERAL + 1];
   int c;

   assert(scip != NULL);

   /* initialize counter for constraint types to zero */
   BMSclearMemoryArray(counter, SCIP_CONSTYPE_GENERAL + 1);

   /* loop through all constraints */
   for( c = 0; c < nconss; c++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      int i;

      /* get constraint */
      cons = conss[c];
      assert(cons != NULL);

      /* get constraint data */
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* merge multiples and delete variables with zero coefficient */
      SCIP_CALL( mergeMultiples(scip, cons) );
      for( i = 0; i < consdata->nvars; i++ )
      {
         assert(!SCIPisZero(scip, consdata->vals[i]));
      }

      /* is constraint of type SCIP_CONSTYPE_EMPTY? */
      if( consdata->nvars == 0 )
      {
         SCIPdebugMessage("classified as EMPTY: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         counter[SCIP_CONSTYPE_EMPTY]++;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_FREE? */
      if( SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIPdebugMessage("classified as FREE: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         counter[SCIP_CONSTYPE_FREE]++;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_SINGLETON? */
      if( consdata->nvars == 1 )
      {
         SCIPdebugMessage("classified as SINGLETON: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         counter[SCIP_CONSTYPE_SINGLETON] += isRangedRow(scip, cons) ? 2 : 1;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_AGGREGATION? */
      if( consdata->nvars == 2 && SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      {
         SCIPdebugMessage("classified as AGGREGATION: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         counter[SCIP_CONSTYPE_AGGREGATION]++;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{VARBOUND}? */
      if( consdata->nvars == 2 )
      {
         SCIPdebugMessage("classified as VARBOUND: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         counter[SCIP_CONSTYPE_VARBOUND] += isRangedRow(scip, cons) ? 2 : 1;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{SETPARTITION, SETPACKING, SETCOVERING, CARDINALITY, INVKNAPSACK}? */
      {
         SCIP_Real scale;
         SCIP_Real b;
         SCIP_Bool unmatched;
         int nnegbinvars;

         unmatched = FALSE;
         nnegbinvars = 0;

         scale = REALABS(consdata->vals[0]);
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0);
            unmatched = unmatched || !SCIPisEQ(scip, REALABS(consdata->vals[i]), scale);

            if( consdata->vals[i] < 0.0 )
               nnegbinvars++;
         }

         if( !unmatched )
         {
            if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
            {
               b = consdata->rhs/scale + nnegbinvars;
               if( SCIPisEQ(scip, 1.0, b) )
               {
                  SCIPdebugMessage("classified as SETPARTITION: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  counter[SCIP_CONSTYPE_SETPARTITION]++;
                  continue;
               }
               else if( SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
               {
                  SCIPdebugMessage("classified as CARDINALITY: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  counter[SCIP_CONSTYPE_CARDINALITY]++;
                  continue;
               }
            }

            b = consdata->rhs/scale + nnegbinvars;
            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMessage("classified as SETPACKING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               counter[SCIP_CONSTYPE_SETPACKING]++;
               consdata->rhs = SCIPinfinity(scip);
            }
            else if( SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
            {
               SCIPdebugMessage("classified as INVKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               counter[SCIP_CONSTYPE_INVKNAPSACK]++;
               consdata->rhs = SCIPinfinity(scip);
            }

            b = consdata->lhs/scale + nnegbinvars;
            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMessage("classified as SETCOVERING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               counter[SCIP_CONSTYPE_SETCOVERING]++;
               consdata->lhs = -SCIPinfinity(scip);
            }

            if( SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
               continue;
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{EQKNAPSACK, BINPACKING, KNAPSACK}? */
      /* @todo If coefficients or rhs are not integral, we currently do not check
       * if the constraint could be scaled (finitely), such that they are.
       */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         b = consdata->rhs;
         unmatched = FALSE;
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0);
            unmatched = unmatched || !SCIPisIntegral(scip, consdata->vals[i]);

            if( SCIPisNegative(scip, consdata->vals[i]) )
               b -= consdata->vals[i];
         }
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         if( !unmatched )
         {
            if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
            {
               SCIPdebugMessage("classified as EQKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               counter[SCIP_CONSTYPE_EQKNAPSACK]++;
               continue;
            }
            else
            {
               SCIP_Bool matched;

               matched = FALSE;
               for( i = 0; i < consdata->nvars && !matched; i++ )
               {
                  matched = matched || SCIPisEQ(scip, b, REALABS(consdata->vals[i]));
               }

               SCIPdebugMessage("classified as %s: ", matched ? "BINPACKING" : "KNAPSACK");
               SCIPdebugPrintCons(scip, cons, NULL);
               counter[matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK]++;
            }

            if( SCIPisInfinity(scip, -consdata->lhs) )
               continue;
            else
               consdata->rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{INTKNAPSACK}? */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         unmatched = FALSE;

         b = consdata->rhs;
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->vars[i]));
            unmatched = unmatched || !SCIPisIntegral(scip, consdata->vals[i]);
            unmatched = unmatched || SCIPisNegative(scip, consdata->vals[i]);
         }

         if( !unmatched )
         {
            SCIPdebugMessage("classified as INTKNAPSACK: ");
            SCIPdebugPrintCons(scip, cons, NULL);
            counter[SCIP_CONSTYPE_INTKNAPSACK]++;

            if( SCIPisInfinity(scip, -consdata->lhs) )
               continue;
            else
               consdata->rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{MIXEDBINARY}? */
      {
         SCIP_Bool unmatched;

         unmatched = FALSE;
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            if( SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS
               && (SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0)
                  || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0)) )
               unmatched = TRUE;
         }

         if( !unmatched )
         {
            SCIPdebugMessage("classified as MIXEDBINARY (%d): ", isRangedRow(scip, cons) ? 2 : 1);
            SCIPdebugPrintCons(scip, cons, NULL);
            counter[SCIP_CONSTYPE_MIXEDBINARY] += isRangedRow(scip, cons) ? 2 : 1;
            continue;
         }
      }

      /* no special structure detected */
      SCIPdebugMessage("classified as GENERAL: ");
      SCIPdebugPrintCons(scip, cons, NULL);
      counter[SCIP_CONSTYPE_GENERAL] += isRangedRow(scip, cons) ? 2 : 1;
   }

   /* print statistics */
   SCIPinfoMessage(scip, NULL, "\nNumber of constraints according to type:\n");
   SCIPinfoMessage(scip, NULL, "----------------------------------------\n");
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_EMPTY        %6d\n",  0, counter[ 0]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_FREE         %6d\n",  1, counter[ 1]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_SINGLETON    %6d\n",  2, counter[ 2]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_AGGREGATION  %6d\n",  3, counter[ 3]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_VARBOUND     %6d\n",  4, counter[ 4]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_SETPARTITION %6d\n",  5, counter[ 5]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_SETPACKING   %6d\n",  6, counter[ 6]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_SETCOVERING  %6d\n",  7, counter[ 7]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_CARDINALITY  %6d\n",  8, counter[ 8]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_INVKNAPSACK  %6d\n",  9, counter[ 9]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_EQKNAPSACK   %6d\n", 10, counter[10]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_BINPACKING   %6d\n", 11, counter[11]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_KNAPSACK     %6d\n", 12, counter[12]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_INTKNAPSACK  %6d\n", 13, counter[13]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_MIXEDBINARY  %6d\n", 14, counter[14]);
   SCIPinfoMessage(scip, NULL, "%2d SCIP_CONSTYPE_GENERAL      %6d\n", 15, counter[15]);
   SCIPinfoMessage(scip, NULL, "----------------------------------------\n\n");

   SCIPinfoMessage(scip, NULL, "    EMPTY");
   SCIPinfoMessage(scip, NULL, "     FREE");
   SCIPinfoMessage(scip, NULL, "     SING");
   SCIPinfoMessage(scip, NULL, "     AGGR");
   SCIPinfoMessage(scip, NULL, "    VARBD");
   SCIPinfoMessage(scip, NULL, "  SETPART");
   SCIPinfoMessage(scip, NULL, "  SETPACK");
   SCIPinfoMessage(scip, NULL, "   SETCOV");
   SCIPinfoMessage(scip, NULL, "     CARD");
   SCIPinfoMessage(scip, NULL, "  INVKNAP");
   SCIPinfoMessage(scip, NULL, "   EQKNAP");
   SCIPinfoMessage(scip, NULL, "  BINPACK");
   SCIPinfoMessage(scip, NULL, "     KNAP");
   SCIPinfoMessage(scip, NULL, "  INTKNAP");
   SCIPinfoMessage(scip, NULL, "   MIXBIN");
   SCIPinfoMessage(scip, NULL, "      GEN\n");
   for( c = 0; c <= SCIP_CONSTYPE_GENERAL; c++ )
   {
      SCIPinfoMessage(scip, NULL, "%9d", counter[c]);
   }

   SCIPinfoMessage(scip, NULL, "\n\n");
   SCIPinfoMessage(scip, NULL, "----------------------------------------\n\n");

   return SCIP_OKAY;
}
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreLinear)
{  /*lint --e{715}*/
   int c;

   /* delete all linear constraints that were upgraded to a more specific constraint type;
    * make sure, only active variables remain in the remaining constraints
    */
   assert(scip != NULL);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
      {
         /* this is no problem reduction, because the upgraded constraint was added to the problem before, and the
          * (redundant) linear constraint was only kept in order to support presolving the the linear constraint handler
          */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
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
SCIP_DECL_CONSEXITSOL(consExitsolLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* if this is a restart, convert cutpool rows into linear constraints */
   if( restart )
   {
      int ncutsadded;

      ncutsadded = 0;

      /* create out of all active cuts in cutpool linear constraints */
      SCIP_CALL( SCIPconvertCutsToConss(scip, NULL, NULL, TRUE, &ncutsadded) );

      if( ncutsadded > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "(restart) converted %d cuts from the global cut pool into linear constraints\n", ncutsadded);
         /* an extra blank line should be printed separately since the buffer message handler only handles up to one
          * line correctly
          */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* free linear constraint */
   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of linear constraints\n");*/

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->vals, sourcedata->lhs, sourcedata->rhs) );

   /* copy the donotupgrade mark */
   targetdata->donotupgrade = sourcedata->donotupgrade;

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
SCIP_DECL_CONSINITLP(consInitlpLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLinear)
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
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMessage("Sepa method of linear constraints\n");*/

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
   maxbound = glblowerbound + conshdlrdata->maxcardbounddist * (cutoffbound - glblowerbound);
   separatecards = SCIPisLE(scip, loclowerbound, maxbound);
   separatecards = separatecards && (SCIPgetNLPBranchCands(scip) > 0);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts && !cutoff; ++c )
   {
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], NULL, separatecards, conshdlrdata->separateall, &ncuts, &cutoff) );
   }
   
   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMessage("Sepa method of linear constraints\n");*/

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
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], sol, TRUE, conshdlrdata->separateall, &ncuts, &cutoff) );
   }

   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLinear)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*SCIPdebugMessage("Enfolp method of linear constraints\n");*/

   /* check for violated constraints
    * LP is processed at current node -> we can add violated linear constraints to the SCIP_LP
    */
   *result = SCIP_FEASIBLE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
         *result = SCIP_SEPARATED;
      }
   }

   /* check all obsolete linear constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
         *result = SCIP_SEPARATED;
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLinear)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* if the solution is infeasible anyway due to objective value, skip the enforcement */
   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, TRUE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLinear)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Check method of linear constraints\n");*/

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, checklprows, &violated) );
   }

   if( violated )
   {
      *result = SCIP_INFEASIBLE;

      if( printreason )
      {
         SCIP_Real activity;
         SCIP_CONSDATA* consdata;

         consdata = SCIPconsGetData(conss[c-1]);
         assert( consdata != NULL);

         activity = consdataGetActivity(scip, consdata, sol);

         SCIP_CALL( SCIPprintCons(scip, conss[c-1], NULL ) );
         SCIPinfoMessage(scip, NULL, ";\n");
         if( SCIPisFeasLT(scip, activity, consdata->lhs) )
            SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhs - activity);

         if( SCIPisFeasGT(scip, activity, consdata->rhs) )
            SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", activity - consdata->rhs);
      }
   }
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool tightenbounds;
   SCIP_Bool cutoff;
   int nchgbds;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /*debugMessage("Prop method of linear constraints\n");*/

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

   /* process useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], tightenbounds, conshdlrdata->sortvars, &cutoff, &nchgbds) );
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


#define MAXCONSPRESOLROUNDS 10
/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Bool minactisrelax;
   SCIP_Bool maxactisrelax;
   SCIP_Bool cutoff;
   SCIP_Bool delay;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int firstchange;
   int firstupgradetry;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Presol method of linear constraints\n");*/

   /* remember old preprocessing counters */
   cutoff = FALSE;
   delay = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process single constraints */
   firstchange = INT_MAX;
   firstupgradetry = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      int npresolrounds;
      SCIP_Bool infeasible;

      infeasible = FALSE;

      cons = conss[c];
      assert(SCIPconsIsActive(cons));
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* constraint should not be already presolved in the initial round */
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->propagated);
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->boundstightened);
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->presolved);
      assert(consdata->propagated || !consdata->presolved);

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, cons, &infeasible) );

      if( infeasible )
      {
         SCIPdebugMessage(" -> infeasible fixing\n");
         cutoff = TRUE;
         break;
      }

      assert(consdata->removedfixings);

      /* we can only presolve linear constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
      if( firstupgradetry == INT_MAX && !consdata->upgradetried )
         firstupgradetry = c;

      /* check, if constraint is already preprocessed */
      if( consdata->presolved )
         continue;

      assert(SCIPconsIsActive(cons));

      SCIPdebugMessage("presolving linear constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      /* apply presolving as long as possible on the single constraint (however, abort after a certain number of rounds
       * to avoid nearly infinite cycling due to very small bound changes)
       */
      npresolrounds = 0;
      while( !consdata->presolved && npresolrounds < MAXCONSPRESOLROUNDS && !SCIPisStopped(scip) )
      {
         assert(!cutoff);
         npresolrounds++;

         /* mark constraint being presolved and propagated */
         consdata->presolved = TRUE;
         consdata->propagated = TRUE;

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );

         /* tighten left and right hand side due to integrality */
         SCIP_CALL( tightenSides(scip, cons, nchgsides) );

         /* check bounds */
         if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible: sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }

         /* tighten variable's bounds */
         SCIP_CALL( tightenBounds(scip, cons, conshdlrdata->sortvars, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* check for fixed variables */
         SCIP_CALL( fixVariables(scip, cons, &cutoff, nfixedvars) );
         if( cutoff )
            break;

         /* check constraint for infeasibility and redundancy */
         consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);
         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) || SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }
         else if( SCIPisFeasGE(scip, minactivity, consdata->lhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( SCIPdelCons(scip, cons) );
            assert(!SCIPconsIsActive(cons));

            if( !consdata->upgraded )
               (*ndelconss)++;
            break;
         }
         else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasGE(scip, minactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> left hand side is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> right hand side is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         assert(consdata->nvars >= 1); /* otherwise, it should be redundant or infeasible */

         /* extract cliques from constraint */
         SCIP_CALL( extractCliques(scip, cons, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* reduce big-M coefficients, that make the constraint redundant if the variable is on a bound */
         SCIP_CALL( consdataTightenCoefs(scip, cons, nchgcoefs, nchgsides) );

         /* try to simplify inequalities */
         if( conshdlrdata->simplifyinequalities )
         {
            SCIP_CALL( simplifyInequalities(scip, cons, nchgcoefs, nchgsides) );
         }
         
         /* aggregation variable in equations */
         if( conshdlrdata->aggregatevariables )
         {
            SCIP_CALL( aggregateVariables(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
            if( cutoff )
               break;
         }
      }

      if( !SCIPisStopped(scip) )
      {
	 /* convert special equalities */
	 if( !cutoff && SCIPconsIsActive(cons) )
	 {
	    SCIP_CALL( convertEquality(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
	 }

	 /* apply dual presolving for variables that appear in only one constraint */
	 if( !cutoff && SCIPconsIsActive(cons) && conshdlrdata->dualpresolving )
	 {
	    SCIP_CALL( dualPresolve(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
	 }

	 /* check if an inequality is parallel to the objective function */
	 if( !cutoff && SCIPconsIsActive(cons) )
	 {
	    SCIP_CALL( checkParallelObjective(scip, cons) );
	 }

	 /* remember the first changed constraint to begin the next aggregation round with */
	 if( firstchange == INT_MAX && consdata->changed )
	    firstchange = c;

	 /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
	 if( firstupgradetry == INT_MAX && !consdata->upgradetried )
	    firstupgradetry = c;
      }
   }

   /* process pairs of constraints: check them for redundancy and try to aggregate them;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    * (otherwise, we delay the presolving to be called again next time)
    */
   if( !cutoff && (conshdlrdata->presolusehashing || conshdlrdata->presolpairwise) && !SCIPisStopped(scip) )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
         && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides )
      {
         assert(firstchange >= 0);

         if( firstchange < nconss && conshdlrdata->presolusehashing ) 
         {
            /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
            SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, &cutoff,
                  ndelconss, nchgsides) );
         }
         
         if( firstchange < nconss && conshdlrdata->presolpairwise )
         {
            SCIP_CONS** usefulconss;
            int nusefulconss;
            int firstchangenew;
            SCIP_Longint npaircomparisons;

            npaircomparisons = 0;
            oldndelconss = *ndelconss;
            oldnchgsides = *nchgsides;
            oldnchgcoefs = *nchgcoefs;

            /* allocate temporary memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &usefulconss, nconss) );
      
            nusefulconss = 0;
            firstchangenew = -1;
            for( c = 0; c < nconss; ++c )
            {
               /* update firstchange */
               if( c == firstchange )
                  firstchangenew = nusefulconss;

               /* ignore inactive and modifiable constraints */
               if( !SCIPconsIsActive(conss[c]) || SCIPconsIsModifiable(conss[c]) )
                  continue;

               usefulconss[nusefulconss] = conss[c];
               nusefulconss++;     
            }
            firstchange = firstchangenew;
            assert(firstchangenew >= 0 && firstchangenew <= nusefulconss);

            for( c = firstchange; c < nusefulconss && !cutoff && !SCIPisStopped(scip); ++c )
            {
               /* constraint has become inactive or modifiable during pairwise presolving */
               if( usefulconss[c] == NULL )
                  continue;
              
               npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? c : (c - firstchange); /*lint !e776*/

               assert(SCIPconsIsActive(usefulconss[c]) && !SCIPconsIsModifiable(usefulconss[c]));
               SCIP_CALL( preprocessConstraintPairs(scip, usefulconss, firstchange, c, conshdlrdata->maxaggrnormscale,
                     &cutoff, ndelconss, nchgsides, nchgcoefs) );

               if( npaircomparisons > conshdlrdata->nmincomparisons )
               {
                  assert(npaircomparisons > 0);
                  if( ((*ndelconss - oldndelconss) + (*nchgsides - oldnchgsides)/2.0 + (*nchgcoefs - oldnchgcoefs)/10.0) / ((SCIP_Real) npaircomparisons) < conshdlrdata->mingainpernmincomp )
                     break;
                  oldndelconss = *ndelconss;
                  oldnchgsides = *nchgsides;
                  oldnchgcoefs = *nchgcoefs;
                  npaircomparisons = 0;
               }
            }
            /* free temporary memory */
            SCIPfreeBufferArray(scip, &usefulconss);
         }
      }
      else
         delay = TRUE;
   }

   /* before upgrading, check whether we can apply some additional dual presolving, because a variable only appears
    * in linear constraints and we therefore have full information about it
    */
   if( !cutoff && firstupgradetry < nconss
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides
      )
   {
      if( conshdlrdata->dualpresolving && !SCIPisStopped(scip) )
      {
         SCIP_CALL( fullDualPresolve(scip, conss, nconss, &cutoff, nchgbds) );
      }
   }

   /* try to upgrade constraints into a more specific constraint type;
    * only upgrade constraints, if no reductions were found in this round (otherwise, the linear constraint handler
    * may find additional reductions before giving control away to other (less intelligent?) constraint handlers)
    */
   if( !cutoff
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides
      )
   {
      for( c = firstupgradetry; c < nconss && !SCIPisStopped(scip); ++c )
      {
         cons = conss[c];

         /* don't upgrade modifiable constraints */
         if( SCIPconsIsModifiable(cons) )
            continue;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         /* only upgrade completely presolved constraints, that changed since the last upgrading call */
         if( consdata->upgradetried )
            continue;
         if( !consdata->presolved )
         {
            delay = TRUE;
            continue;
         }

         consdata->upgradetried = TRUE;
         if( SCIPconsIsActive(cons) )
         {
            SCIP_CONS* upgdcons;

            SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
            if( upgdcons != NULL )
            {
               /* add the upgraded constraint to the problem */
               SCIP_CALL( SCIPaddCons(scip, upgdcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &upgdcons) );
               (*nupgdconss)++;

               /* mark the linear constraint being upgraded and to be removed after presolving;
                * don't delete it directly, because it may help to preprocess other linear constraints
                */
               assert(!consdata->upgraded);
               consdata->upgraded = TRUE;

               /* delete upgraded inequalities immediately;
                * delete upgraded equalities, if we don't need it anymore for aggregation and redundancy checking
                */
               if( SCIPisLT(scip, consdata->lhs, consdata->rhs)
                  || !conshdlrdata->presolpairwise
                  || (conshdlrdata->maxaggrnormscale == 0.0) )
               {
                  SCIP_CALL( SCIPdelCons(scip, cons) );
               }
            }
         }
      }
   }
   else
      delay = TRUE;

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nupgdconss > oldnupgdconss || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLinear)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   SCIP_CALL( resolvePropagation(scip, cons, infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* update rounding locks of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisPositive(scip, consdata->vals[i]) )
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsLinear)
{
   assert(scip != NULL);
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
SCIP_DECL_CONSPRINT(consPrintLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   
   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );
    
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLinear)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_Real* sourcecoefs;
   const char* consname;
   int nvars;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsLinear(sourcescip, sourcecons);
   sourcecoefs = SCIPgetValsLinear(sourcescip, sourcecons); 
   nvars = SCIPgetNVarsLinear(sourcescip, sourcecons);
   
   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, sourcecoefs,
         SCIPgetLhsLinear(sourcescip, sourcecons), SCIPgetRhsLinear(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL || *valid == FALSE);

   return SCIP_OKAY;
}   

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int        nvars;
   int        coefssize;
   int        requsize;
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   char*      endptr;

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

   /* initialize buffers for storing the variables and coefficients */
   coefssize = 100;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  coefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, coefssize) );

   /* parse linear sum to get variables and coefficients */
   SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );

   if( *success && requsize > coefssize )
   {
      /* realloc buffers and try again */
      coefssize = requsize;
      SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  coefssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, coefssize) );

      SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );
      assert(!*success || requsize <= coefssize); /* if successful, then should have had enough space now */
   }

   if( !*success )
   {
      SCIPerrorMessage("no luck in parsing linear sum '%s'\n", str);
   }
   else
   {
      (*success) = FALSE;
      str = endptr;

      /* check for left or right hand side */
      while( isspace((unsigned char)*str) )
         ++str;

      /* check for free constraint */
      if( strncmp(str, "[free]", 6) == 0 )
      {
         if( !SCIPisInfinity(scip, -lhs) )
         {
            SCIPerrorMessage("cannot have left hand side and [free] status \n");
            return SCIP_OKAY;
         }
         (*success) = TRUE;
      }
      else
      {
         switch( *str )
         {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, &endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, &endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, &endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            return SCIP_OKAY;
         }
      }

      if( *success )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, coefs, lhs, rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }
   }

   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsLinear)
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

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);

   if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      /* update the activity values */
      if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
         consdataUpdateActivitiesLb(scip, consdata, var, oldbound, newbound, val, TRUE);
      else
      {
         assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
         consdataUpdateActivitiesUb(scip, consdata, var, oldbound, newbound, val, TRUE);
      }

      consdata->presolved = FALSE;

      /* bound change can turn the constraint infeasible or redundant only if it was a tightening */
      if( (eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0 )
         consdata->propagated = FALSE;

      /* check whether bound tightening might now be successful (if the current bound was relaxed, it might be
       * that it can be tightened again)
       */
      if( consdata->boundstightened )
      {
         switch( eventtype )
         {
         case SCIP_EVENTTYPE_LBTIGHTENED:
         case SCIP_EVENTTYPE_UBRELAXED:
            consdata->boundstightened = (val > 0.0 && SCIPisInfinity(scip, consdata->rhs))
               || (val < 0.0 && SCIPisInfinity(scip, -consdata->lhs));
            break;
         case SCIP_EVENTTYPE_LBRELAXED:
         case SCIP_EVENTTYPE_UBTIGHTENED:
            consdata->boundstightened = (val > 0.0 && SCIPisInfinity(scip, -consdata->lhs))
               || (val < 0.0 && SCIPisInfinity(scip, consdata->rhs));
            break;
         default:
            SCIPerrorMessage("invalid event type %d\n", eventtype);
            return SCIP_INVALIDDATA;
         }
      }
   }
   else if( (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0 )
   {
      /* we want to remove the fixed variable */
      consdata->presolved = FALSE;
      consdata->removedfixings = FALSE;
   }

   else if( (eventtype & SCIP_EVENTTYPE_VARUNLOCKED) != 0 )
   {
      /* there is only one lock left: we may multi-aggregate the variable as slack of an equation */
      assert(SCIPvarGetNLocksDown(var) <= 1);
      assert(SCIPvarGetNLocksUp(var) <= 1);
      consdata->presolved = FALSE;
   }
   else if( (eventtype & SCIP_EVENTTYPE_GBDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      /* update the activity values */
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
         consdataUpdateActivitiesGlbLb(scip, consdata, oldbound, newbound, val, TRUE);
      else
      {
         assert((eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0);
         consdataUpdateActivitiesGlbUb(scip, consdata, oldbound, newbound, val, TRUE);
      }
   }
   else if( (eventtype & SCIP_EVENTTYPE_OBJCHANGED) != 0 )
   {
      /* old and new objective value */
      SCIP_Real oldobj;
      SCIP_Real newobj;

      /* get old and new objective value */
      oldobj = SCIPeventGetOldobj(event);
      newobj = SCIPeventGetNewobj(event);

      if( (oldobj >= 0.0 && newobj < 0.0) || (oldobj < 0.0 && newobj >= 0.0) )
         consdataInvalidateActivities(consdata);
   }
   else
   {
      assert((eventtype & SCIP_EVENTTYPE_VARDELETED) != 0);
      consdata->varsdeleted = TRUE;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   int i;

   assert(scip != NULL);
   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* create array of variables and coefficients: sum_{i \in P} x_i - sum_{i \in N} x_i >= 1 - |N| */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nbdchginfos) );
   lhs = 1.0;
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      /**@todo extend linear conflict constraints to some non-binary cases */
      if( !SCIPvarIsBinary(vars[i]) )
         break;

      /* check whether the variable is fixed to zero (P) or one (N) in the conflict set */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) < 0.5 )
         vals[i] = 1.0;
      else
      {
         vals[i] = -1.0;
         lhs -= 1.0;
      }
   }

   if( i == nbdchginfos )
   {
      SCIP_CONS* cons;
      SCIP_CONS* upgdcons;
      char consname[SCIP_MAXSTRLEN];

      /* create a constraint out of the conflict set */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cf%"SCIP_LONGINT_FORMAT, SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, nbdchginfos, vars, vals, lhs, SCIPinfinity(scip),
            FALSE, separate, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );

      /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
      SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         cons = upgdcons;
      }
      
      /* add constraint to SCIP */
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      
      *result = SCIP_CONSADDED;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * Quadratic constraint upgrading
 */


/** upgrades quadratic constraints with only and at least one linear variables into a linear constraint
 */
static
SCIP_DECL_QUADCONSUPGD(upgradeConsQuadratic)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMessage("upgradeConsQuadratic called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   if( SCIPgetNQuadVarTermsQuadratic(scip, cons) > 0 )
      return SCIP_OKAY;
   if( SCIPgetNLinearVarsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* signal that we need more memory */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsQuadratic(scip, cons),
         SCIPgetLinearVarsQuadratic(scip, cons),
         SCIPgetCoefsLinearVarsQuadratic(scip, cons),
         SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );
   SCIPdebugMessage("created linear constraint:\n");
   SCIPdebugPrintCons(scip, upgdconss[0], NULL);

   return SCIP_OKAY;
}

/** tries to upgrade a nonlinear constraint into a linear constraint */
static
SCIP_DECL_NONLINCONSUPGD(upgradeConsNonlinear)
{
   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   /* no interest in nonlinear constraints */
   if( SCIPgetExprgraphNodeNonlinear(scip, cons) != NULL )
      return SCIP_OKAY;

   /* no interest in constant constraints */
   if( SCIPgetNLinearVarsNonlinear(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
         SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_CONFLICTHDLR* conflicthdlr;

   assert(scip != NULL);

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecLinear, NULL) );

   /* create conflict handler for linear constraints */
   SCIP_CALL( SCIPincludeConflicthdlrBasic(scip, &conflicthdlr, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         conflictExecLinear, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpLinear, consEnfopsLinear, consCheckLinear, consLockLinear,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyLinear, consCopyLinear) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteLinear) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsLinear) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitLinear) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreLinear) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolLinear) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeLinear) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsLinear) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsLinear) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitLinear) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpLinear) );
#ifdef WITH_PRINTORIGCONSTYPES
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreLinear) );
#endif
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseLinear) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolLinear, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintLinear) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropLinear, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropLinear) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpLinear, consSepasolLinear, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransLinear) );

   if( SCIPfindConshdlr(scip, "quadratic") != NULL )
   {
      /* include function that upgrades quadratic constraint to linear constraints */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, upgradeConsQuadratic, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the linear constraint upgrade in the nonlinear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, upgradeConsNonlinear, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add linear constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/tightenboundsfreq",
         "multiplier on propagation frequency, how often the bounds are tightened (-1: never, 0: only at root)",
         &conshdlrdata->tightenboundsfreq, TRUE, DEFAULT_TIGHTENBOUNDSFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/nmincomparisons",
         "number for minimal pairwise presolve comparisons",
         &conshdlrdata->nmincomparisons, TRUE, DEFAULT_NMINCOMPARISONS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/linear/mingainpernmincomparisons",
         "minimal gain per minimal pairwise presolve comparisons to repeat pairwise comparison round",
         &conshdlrdata->mingainpernmincomp, TRUE, DEFAULT_MINGAINPERNMINCOMP, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/linear/maxaggrnormscale",
         "maximal allowed relative gain in maximum norm for constraint aggregation (0.0: disable constraint aggregation)",
         &conshdlrdata->maxaggrnormscale, TRUE, DEFAULT_MAXAGGRNORMSCALE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/linear/maxcardbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for separating knapsack cardinality cuts",
         &conshdlrdata->maxcardbounddist, TRUE, DEFAULT_MAXCARDBOUNDDIST, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/separateall",
         "should all constraints be subject to cardinality cut generation instead of only the ones with non-zero dual value?",
         &conshdlrdata->separateall, FALSE, DEFAULT_SEPARATEALL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/aggregatevariables",
         "should presolving search for aggregations in equations",
         &conshdlrdata->aggregatevariables, TRUE, DEFAULT_AGGREGATEVARIABLES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/simplifyinequalities",
         "should presolving try to simplify inequalities",
         &conshdlrdata->simplifyinequalities, TRUE, DEFAULT_SIMPLIFYINEQUALITIES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/linear/sortvars", "apply binaries sorting in decr. order of coeff abs value?",
         &conshdlrdata->sortvars, TRUE, DEFAULT_SORTVARS, NULL, NULL) );

   return SCIP_OKAY;
}

/** includes a linear constraint update method into the linear constraint handler */
SCIP_RETCODE SCIPincludeLinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority,           /**< priority of upgrading method */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LINCONSUPGRADE* linconsupgrade;
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];
   
   assert(scip != NULL);
   assert(linconsupgd != NULL);
   assert(conshdlrname != NULL );

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check if linear constraint update method already exists in constraint handler data */
   if( !conshdlrdataHasUpgrade(scip, conshdlrdata, linconsupgd, conshdlrname) )
   {
      /* create a linear constraint upgrade data object */
      SCIP_CALL( linconsupgradeCreate(scip, &linconsupgrade, linconsupgd, priority) );
      
      /* insert linear constraint update method into constraint handler data */
      SCIP_CALL( conshdlrdataIncludeUpgrade(scip, conshdlrdata, linconsupgrade) );
    
      /* adds parameter to turn on and off the upgrading step */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/linear/upgrade/%s", conshdlrname);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable linear upgrading for constraint handler <%s>", conshdlrname);
      SCIP_CALL( SCIPaddBoolParam(scip,
            paramname, paramdesc,
            &linconsupgrade->active, FALSE, TRUE, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
   SCIP_Bool             dynamic,            /**< Is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
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

   assert(scip != NULL);
   assert(cons != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE && nvars > 0 )
   {
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_Real constant = 0.0;
      int nconsvars;
      int requiredsize;

      nconsvars = nvars;
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvars, vars, nconsvars) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvals, vals, nconsvars) );

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, nconsvars, &constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      /* adjust sides and check that we do not subtract infinity values */
      if( SCIPisInfinity(scip, REALABS(constant)) )
      {
	 if( constant < 0.0 )
	 {
	    if( SCIPisInfinity(scip, lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = -SCIPinfinity(scip);
	    rhs = -SCIPinfinity(scip);
	 }
	 else
	 {
	    if( SCIPisInfinity(scip, -lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, -rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = SCIPinfinity(scip);
	    rhs = SCIPinfinity(scip);
	 }
      }
      else
      {
	 if( !SCIPisInfinity(scip, REALABS(lhs)) )
	    lhs -= constant;
	 if( !SCIPisInfinity(scip, REALABS(rhs)) )
	    rhs -= constant;

	 if( SCIPisInfinity(scip, -lhs) )
	    lhs = -SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, lhs) )
	    lhs = SCIPinfinity(scip);

	 if( SCIPisInfinity(scip, rhs) )
	    rhs = SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, -rhs) )
	    rhs = -SCIPinfinity(scip);
      }

      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nconsvars, consvars, consvals, lhs, rhs) );
      assert(consdata != NULL);

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, vals, lhs, rhs) );
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
SCIP_RETCODE SCIPcreateConsBasicLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, vals, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates by copying and captures a linear constraint */
SCIP_RETCODE SCIPcopyConsLinear(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in source variable array */
   SCIP_VAR**            sourcevars,         /**< source variables of the linear constraints */
   SCIP_Real*            sourcecoefs,        /**< coefficient array of the linear constraint, or NULL if all coefficients are one */
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

   if( SCIPisGT(scip, lhs, rhs) )
   {
      *valid = FALSE;
      return SCIP_OKAY;
   }

   (*valid) = TRUE;

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
      SCIP_CALL( SCIPduplicateBufferArray(scip, &coefs, sourcecoefs, nvars) );
   }
   else 
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
      for( v = 0; v < nvars; ++v )
         coefs[v] = 1.0;      
   }

   constant = 0.0;
   
   /* transform source variable to active variables of the source SCIP since only these can be mapped to variables of
    * the target SCIP */
   if( SCIPisTransformed(sourcescip) )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, nvars, &constant, &requiredsize, TRUE) );
      
      if( requiredsize > nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requiredsize) );
         
         SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nvars);
      }
   }
   else
   {
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &coefs[v], &constant) );
      }
   }
   
   /* map variables of the source constraint to variables of the target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_VAR* var;
      var = vars[v];
      
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
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
SCIP_RETCODE SCIPaddCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_Real constant = 0.0;
      SCIP_Real rhs;
      SCIP_Real lhs;
      int nconsvars;
      int requiredsize;
      int v;

      nconsvars = 1;
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
      consvars[0] = var;
      consvals[0] = val;

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, nconsvars, &constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      /* adjust sides and check that we do not subtract infinity values */
      /* constant is infinite */
      if( SCIPisInfinity(scip, REALABS(constant)) )
      {
	 if( constant < 0.0 )
	 {
	    if( SCIPisInfinity(scip, lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = -SCIPinfinity(scip);
	    rhs = -SCIPinfinity(scip);
	 }
	 else
	 {
	    if( SCIPisInfinity(scip, -lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, -rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = SCIPinfinity(scip);
	    rhs = SCIPinfinity(scip);
	 }
      }
      /* constant is not infinite */
      else
      {
	 if( !SCIPisInfinity(scip, REALABS(lhs)) )
	    lhs -= constant;
	 if( !SCIPisInfinity(scip, REALABS(rhs)) )
	    rhs -= constant;

	 if( SCIPisInfinity(scip, -lhs) )
	    lhs = -SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, lhs) )
	    lhs = SCIPinfinity(scip);

	 if( SCIPisInfinity(scip, rhs) )
	    rhs = SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, -rhs) )
	    rhs = -SCIPinfinity(scip);
      }

      /* add all active variables to constraint */
      for( v = nconsvars - 1; v >= 0; --v )
      {
	 SCIP_CALL( addCoef(scip, cons, consvars[v], consvals[v]) );
      }

      /* update left and right hand sides */
      SCIP_CALL( chgLhs(scip, cons, lhs));
      SCIP_CALL( chgRhs(scip, cons, rhs));

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      SCIP_CALL( addCoef(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
SCIP_Real SCIPgetLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of linear constraint */
SCIP_Real SCIPgetRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** changes left hand side of linear constraint */
SCIP_RETCODE SCIPchgLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
SCIP_RETCODE SCIPchgRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

/** gets the number of variables in the linear constraint */
int SCIPgetNVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the linear constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_Real* SCIPgetValsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vals;
}

/** gets the activity of the linear constraint in the given solution */
SCIP_Real SCIPgetActivityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolActivity(scip, consdata->row, sol);
   else
      return consdataGetActivity(scip, consdata, sol);
}

/** gets the feasibility of the linear constraint in the given solution */
SCIP_Real SCIPgetFeasibilityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolFeasibility(scip, consdata->row, sol);
   else
      return consdataGetFeasibility(scip, consdata, sol);
}

/** gets the dual solution of the linear constraint in the current LP */
SCIP_Real SCIPgetDualsolLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the linear constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
SCIP_RETCODE SCIPupgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real poscoeffsum;
   SCIP_Real negcoeffsum;
   SCIP_Bool integral;
   int nposbin;
   int nnegbin;
   int nposint;
   int nnegint;
   int nposimpl;
   int nnegimpl;
   int nposcont;
   int nnegcont;
   int ncoeffspone;
   int ncoeffsnone;
   int ncoeffspint;
   int ncoeffsnint;
   int ncoeffspfrac;
   int ncoeffsnfrac;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(upgdcons != NULL);

   *upgdcons = NULL;

   /* we cannot upgrade a modifiable linear constraint, since we don't know what additional coefficients to expect */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get the constraint handler and check, if it's really a linear constraint */
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* get constraint handler data and constraint data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint was already upgraded and will be deleted anyway after preprocessing */
   if( consdata->upgraded )
      return SCIP_OKAY;

   /* do not upgrade marked constraints */
   if ( consdata->donotupgrade )
      return SCIP_OKAY;

   /* check, if the constraint is already stored as LP row */
   if( consdata->row != NULL )
   {
      if( SCIProwIsInLP(consdata->row) )
      {
         SCIPerrorMessage("cannot upgrade linear constraint that is already stored as row in the LP\n");
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* normalize constraint */
   SCIP_CALL( normalizeCons(scip, cons) );


   /*
    * calculate some statistics on linear constraint
    */

   nposbin = 0;
   nnegbin = 0;
   nposint = 0;
   nnegint = 0;
   nposimpl = 0;
   nnegimpl = 0;
   nposcont = 0;
   nnegcont = 0;
   ncoeffspone = 0;
   ncoeffsnone = 0;
   ncoeffspint = 0;
   ncoeffsnint = 0;
   ncoeffspfrac = 0;
   ncoeffsnfrac = 0;
   integral = TRUE;
   poscoeffsum = 0.0;
   negcoeffsum = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      var = consdata->vars[i];
      val = consdata->vals[i];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, val));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposbin++;
         else
            nnegbin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposint++;
         else
            nnegint++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposimpl++;
         else
            nnegimpl++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisEQ(scip, lb, ub) && SCIPisIntegral(scip, val * lb);
         if( val >= 0.0 )
            nposcont++;
         else
            nnegcont++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
      if( SCIPisEQ(scip, val, 1.0) )
         ncoeffspone++;
      else if( SCIPisEQ(scip, val, -1.0) )
         ncoeffsnone++;
      else if( SCIPisIntegral(scip, val) )
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspint++;
         else
            ncoeffsnint++;
      }
      else
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspfrac++;
         else
            ncoeffsnfrac++;
      }
      if( SCIPisPositive(scip, val) )
         poscoeffsum += val;
      else
         negcoeffsum += val;
   }


   /*
    * call the upgrading methods
    */

   SCIPdebugMessage("upgrading linear constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nlinconsupgrades);
   SCIPdebugMessage(" +bin=%d -bin=%d +int=%d -int=%d +impl=%d -impl=%d +cont=%d -cont=%d +1=%d -1=%d +I=%d -I=%d +F=%d -F=%d possum=%.15g negsum=%.15g integral=%u\n",
      nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
      ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
      poscoeffsum, negcoeffsum, integral);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nlinconsupgrades && *upgdcons == NULL; ++i )
   {
      if( conshdlrdata->linconsupgrades[i]->active )
      {
         SCIP_CALL( conshdlrdata->linconsupgrades[i]->linconsupgd(scip, cons, consdata->nvars,
               consdata->vars, consdata->vals, consdata->lhs, consdata->rhs,
               nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
               ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
               poscoeffsum, negcoeffsum, integral,
               upgdcons) );
      }
   }

#ifdef SCIP_DEBUG
   if( *upgdcons != NULL )
   {
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMessage(" -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(SCIPconsGetHdlr(*upgdcons)));
      SCIPdebugPrintCons(scip, *upgdcons, NULL);
   }
#endif

   return SCIP_OKAY;
}


/** forbids upgrading of constraint */
SCIP_RETCODE SCIPmarkDoNotUpgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint to mark */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* get the constraint handler and check, if it's really a linear constraint */
   conshdlr = SCIPconsGetHdlr(cons);
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* get data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->donotupgrade = TRUE;

   return SCIP_OKAY;
}

/** sets upgrading flag of linear constraint 
 *
 *  @note the donotupgrade flag should only be changed from TRUE to FALSE, by the caller who set it to TRUE
 */
SCIP_RETCODE SCIPsetUpgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint to mark */
   SCIP_Bool             upgradeallowed      /**< allow upgrading? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   /* get the constraint handler and check, if it's really a linear constraint */
   conshdlr = SCIPconsGetHdlr(cons);
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* get data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!consdata->upgraded);

   if( upgradeallowed && consdata->donotupgrade )
   {
      consdata->donotupgrade = FALSE;

      /* update the upgrade flag to try again */
      consdata->upgradetried = FALSE;
   }
   else if( !upgradeallowed )
   {
      if( consdata->donotupgrade )
      {
         /* @todo: change donotupgrade flag to a counter */
         SCIPwarningMessage(scip, "constraint is already marked not to be upgraded\n");
      }
      consdata->donotupgrade = TRUE;
   }
   return SCIP_OKAY;
}
