/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr.c
 * @brief  constraint handler for expression constraints (in particular, nonlinear constraints)
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef SCIP_DEBUG
#define ENFO_LOGGING
#endif

/* enable to get log output for enforcement */
/* #define ENFO_LOGGING */
/* define to get enforcement logging into file */
/* #define ENFOLOGFILE "consexpr_enfo.log" */

/* define to get more debug output from domain propagation */
/* #define DEBUG_PROP */

/*lint -e528*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_expr.h"
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/struct_cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_log.h"
#include "scip/cons_expr_abs.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_entropy.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_cos.h"
#include "scip/cons_expr_nlhdlr_bilinear.h"
#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr_nlhdlr_default.h"
#include "scip/cons_expr_nlhdlr_quadratic.h"
#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr_iterator.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "expr"
#define CONSHDLR_DESC          "constraint handler for expressions"
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */

#define VERTEXPOLY_MAXPERTURBATION      1e-3 /**< maximum perturbation */
#define VERTEXPOLY_USEDUALSIMPLEX       TRUE /**< use dual or primal simplex algorithm? */
#define VERTEXPOLY_RANDNUMINITSEED  20181029 /**< seed for random number generator, which is used to move points away from the boundary */
#define VERTEXPOLY_ADJUSTFACETFACTOR     1e1 /**< adjust resulting facets in checkRikun() up to a violation of this value times lpfeastol */

#define BRANCH_RANDNUMINITSEED      20191229 /**< seed for random number generator, which is used to select from several similar good branching candidates */

/* properties of the expression constraint handler statistics table */
#define TABLE_NAME_EXPR                          "expression"
#define TABLE_DESC_EXPR                          "expression constraint handler statistics"
#define TABLE_POSITION_EXPR                      12500                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_EXPR                SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

/* enable nonlinear constraint upgrading */
#include "scip/cons_nonlinear.h"
#define NONLINCONSUPGD_PRIORITY   600000 /**< priority of the constraint handler for upgrading of nonlinear constraints */

/* enable quadratic constraint upgrading */
#include "scip/cons_quadratic.h"
#define QUADCONSUPGD_PRIORITY     600000 /**< priority of the constraint handler for upgrading of quadratic constraints */

/** ensures that a block memory array has at least a given size
 *
 *  if cursize is 0, then *array1 can be NULL
 */
#define ENSUREBLOCKMEMORYARRAYSIZE(scip, array1, cursize, minsize)      \
   do {                                                                 \
      int __newsize;                                                    \
      assert((scip)  != NULL);                                          \
      if( (cursize) >= (minsize) )                                      \
         break;                                                         \
      __newsize = SCIPcalcMemGrowSize(scip, minsize);                   \
      assert(__newsize >= (minsize));                                   \
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(array1), cursize, __newsize) ); \
      (cursize) = __newsize;                                            \
   } while( FALSE )

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/** translates x to 2^x for non-negative integer x */
#define POWEROFTWO(x) (0x1u << (x))

#ifdef ENFO_LOGGING
#define ENFOLOG(x) if( SCIPgetSubscipDepth(scip) == 0 && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL ) { x }

FILE* enfologfile = NULL;

#else
#define ENFOLOG(x)
#endif

/*
 * Data structures
 */

/** expression constraint update method */
struct SCIP_ExprConsUpgrade
{
   SCIP_DECL_EXPRCONSUPGD((*exprconsupgd));  /**< method to call for upgrading expression constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};
typedef struct SCIP_ExprConsUpgrade SCIP_EXPRCONSUPGRADE;

/** constraint data for expr constraints */
struct SCIP_ConsData
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression that represents this constraint */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   SCIP_CONSEXPR_EXPR**  varexprs;           /**< array containing all variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   SCIP_Bool             catchedevents;      /**< do we catch events on variables? */

   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution (used temporarily inside constraint handler) */
   SCIP_Real             gradnorm;           /**< norm of gradient of constraint function in current solution (if evaluated) */
   unsigned int          gradnormsoltag;     /**< tag of solution used that gradnorm corresponds to */

   unsigned int          ispropagated:1;     /**< did we propagate the current bounds already? */
   unsigned int          issimplified:1;     /**< did we simplify the expression tree already? */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   int                   nlockspos;          /**< number of positive locks */
   int                   nlocksneg;          /**< number of negative locks */

   /* repair infeasible solutions */
   SCIP_VAR*             linvardecr;         /**< variable that may be decreased without making any other constraint infeasible, or NULL if none */
   SCIP_VAR*             linvarincr;         /**< variable that may be increased without making any other constraint infeasible, or NULL if none */
   SCIP_Real             linvardecrcoef;     /**< linear coefficient of linvardecr */
   SCIP_Real             linvarincrcoef;     /**< linear coefficient of linvarincr */

   int                   consindex;          /**< an index of the constraint that is unique among all expr-constraints in this SCIP instance and is constant */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* expression handler */
   SCIP_CONSEXPR_EXPRHDLR** exprhdlrs;       /**< expression handlers */
   int                      nexprhdlrs;      /**< number of expression handlers */
   int                      exprhdlrssize;   /**< size of exprhdlrs array */

   SCIP_CONSEXPR_EXPRHDLR*  exprvarhdlr;     /**< variable expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprvalhdlr;     /**< value expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprsumhdlr;     /**< summation expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprprodhdlr;    /**< product expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprpowhdlr;     /**< power expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprsignpowhdlr; /**< signed power expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprexphdlr;     /**< exponential expression handler */

   /* nonlinear handler */
   SCIP_CONSEXPR_NLHDLR**   nlhdlrs;         /**< nonlinear handlers */
   int                      nnlhdlrs;        /**< number of nonlinear handlers */
   int                      nlhdlrssize;     /**< size of nlhdlrs array */

   /* constraint upgrades */
   SCIP_EXPRCONSUPGRADE**   exprconsupgrades;     /**< nonlinear constraint upgrade methods for specializing expression constraints */
   int                      exprconsupgradessize; /**< size of exprconsupgrades array */
   int                      nexprconsupgrades;    /**< number of expression constraint upgrade methods */

   /* other plugins */
   SCIP_EVENTHDLR*          eventhdlr;       /**< handler for variable bound change events */
   SCIP_HEUR*               subnlpheur;      /**< a pointer to the subnlp heuristic, if available */
   SCIP_HEUR*               trysolheur;      /**< a pointer to the trysol heuristic, if available */

   /* expression iterator */
   int                      nactiveiter;     /**< number of currently active iterators */
   unsigned int             lastvisitedtag;  /**< last visited tag used by iterators */

   /* tags and counters */
   int                      auxvarid;        /**< unique id for the next auxiliary variable */

   unsigned int             lastsoltag;      /**< last solution tag used to evaluate current solution */
   unsigned int             curboundstag;    /**< tag indicating current variable bounds */
   unsigned int             lastboundrelax;  /**< tag when bounds where most recently relaxed */
   unsigned int             lastdifftag;     /**< last tag used for computing gradients */
   unsigned int             enforound;       /**< total number of enforcement calls, including current one */

   int                      lastconsindex;   /**< last used consindex, plus one */

   /* parameters */
   int                      maxproprounds;   /**< limit on number of propagation rounds for a set of constraints within one round of SCIP propagation */
   char                     varboundrelax;   /**< strategy on how to relax variable bounds during bound tightening */
   SCIP_Real                varboundrelaxamount; /**< by how much to relax variable bounds during bound tightening */
   SCIP_Real                conssiderelaxamount; /**< by how much to relax constraint sides during bound tightening */
   SCIP_Real                vp_maxperturb;   /**< maximal relative perturbation of reference point */
   SCIP_Real                vp_adjfacetthreshold; /**< adjust computed facet up to a violation of this value times lpfeastol */
   SCIP_Bool                vp_dualsimplex;  /**< whether to use dual simplex instead of primal simplex for facet computing LP */
   SCIP_Bool                reformbinprods;  /**< whether to reformulate products of binary variables during presolving */
   SCIP_Bool                reformbinprodsand;/**< whether to use the AND constraint handler for reformulating binary products */
   int                      reformbinprodsfac; /**< minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled) */
   SCIP_Bool                tightenlpfeastol;/**< whether to tighten LP feasibility tolerance during enforcement, if it seems useful */
   SCIP_Bool                propinenforce;   /**< whether to (re)run propagation in enforcement */
   SCIP_Real                weakcutthreshold;/**< threshold for when to regard a cut from an estimator as weak */
   SCIP_Real                strongcutmaxcoef;/**< "strong" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef] */
   SCIP_Bool                strongcutefficacy;/**< consider efficacy requirement when deciding whether a cut is "strong" */
   SCIP_Bool                forcestrongcut;  /**< whether to force "strong" cuts in enforcement */
   SCIP_Real                enfoauxviolfactor;/**< an expression will be enforced if the "auxiliary" violation is at least enfoauxviolfactor times the "original" violation */
   SCIP_Real                weakcutminviolfactor; /**< retry with weak cuts for constraints with violation at least this factor of maximal violated constraints */
   char                     violscale;       /**< method how to scale violations to make them comparable (not used for feasibility check) */
   char                     checkvarlocks;   /**< whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction) */
   int                      branchauxmindepth; /**< from which depth on to allow branching on auxiliary variables */
   SCIP_Bool                branchexternal;  /**< whether to use external branching candidates for branching */
   SCIP_Real                branchhighviolfactor; /**< consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints */
   SCIP_Real                branchhighscorefactor; /**< consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables */
   SCIP_Real                branchviolweight;/**< weight by how much to consider the violation assigned to a variable for its branching score */
   SCIP_Real                branchdualweight;/**< weight by how much to consider the dual values of rows that contain a variable for its branching score */
   SCIP_Real                branchpscostweight;/**< weight by how much to consider the pseudo cost of a variable for its branching score */
   SCIP_Real                branchdomainweight; /**< weight by how much to consider the domain width in branching score */
   SCIP_Real                branchvartypeweight;/**< weight by how much to consider variable type in branching score */
   char                     branchscoreagg;  /**< how to aggregate several branching scores given for the same expression ('a'verage, 'm'aximum, or 's'um) */
   char                     branchviolsplit; /**< method used to split violation in expression onto variables ('e'venly, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width) */
   SCIP_Real                branchpscostreliable; /**< minimum pseudo-cost update count required to consider pseudo-costs reliable */

   /* statistics */
   SCIP_Longint             nweaksepa;       /**< number of times we used "weak" cuts for enforcement */
   SCIP_Longint             ntightenlp;      /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing */
   SCIP_Longint             ndesperatetightenlp; /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing because we didn't know anything better */
   SCIP_Longint             ndesperatebranch;/**< number of times we branched on some variable because normal enforcement was not successful */
   SCIP_Longint             ndesperatecutoff;/**< number of times we cut off a node in enforcement because no branching candidate could be found */
   SCIP_Longint             nforcelp;        /**< number of times we forced solving the LP when enforcing a pseudo solution */
   SCIP_CLOCK*              canonicalizetime;/**< time spend for canonicalization */
   SCIP_Longint             ncanonicalizecalls; /**< number of times we called canonicalization */

   /* facets of envelops of vertex-polyhedral functions */
   SCIP_RANDNUMGEN*         vp_randnumgen;   /**< random number generator used to perturb reference point */
   SCIP_LPI*                vp_lp[SCIP_MAXVERTEXPOLYDIM+1];  /**< LPs used to compute facets for functions of different dimension */

   /* hashing of bilinear terms */
   SCIP_HASHTABLE*          bilinhashtable;  /**< hash table for bilinear terms */
   SCIP_CONSEXPR_BILINTERM* bilinterms;      /**< bilinear terms */
   int                      nbilinterms;     /**< total number of bilinear terms */
   int                      bilintermssize;  /**< size of bilinterms array */

   /* branching */
   SCIP_RANDNUMGEN*         branchrandnumgen;/**< random number generated used in branching variable selection */
   char                     branchpscostupdatestrategy; /**< value of parameter branching/lpgainnormalize */

   /* misc */
   SCIP_Bool                checkedvarlocks; /**< whether variables contained in a single constraint have been already considered */
};

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*           varmap;           /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*           consmap;          /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool               global;           /**< should a global or a local copy be created */
   SCIP_Bool               valid;            /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** printing to dot file data */
struct SCIP_ConsExpr_PrintDotData
{
   FILE*                   file;             /**< file to print to */
   SCIP_CONSEXPR_ITERATOR* iterator;         /**< iterator to use */
   SCIP_Bool               closefile;        /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*           leaveexprs;       /**< hashmap storing leave (no children) expressions */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint;  /**< flags that indicate what to print for each expression */
};

/** branching candidate with various scores */
typedef struct
{
   SCIP_CONSEXPR_EXPR*     expr;             /**< expression that holds branching candidate */
   SCIP_Real               auxviol;          /**< aux-violation score of candidate */
   SCIP_Real               domain;           /**< domain score of candidate */
   SCIP_Real               dual;             /**< dual score of candidate */
   SCIP_Real               pscost;           /**< pseudo-cost score of candidate */
   SCIP_Real               vartype;          /**< variable type score of candidate */
   SCIP_Real               weighted;         /**< weighted sum of other scores, see scoreBranchingCandidates() */
} BRANCHCAND;

/*
 * Local methods
 */

/** creates an expression */
static
SCIP_RETCODE createExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_CONSEXPR_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   )
{
   int c;

   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);
   assert(exprdata == NULL || exprhdlr->copydata != NULL); /* copydata must be available if there is expression data */
   assert(exprdata == NULL || exprhdlr->freedata != NULL); /* freedata must be available if there is expression data */

   SCIP_CALL( SCIPallocClearBlockMemory(scip, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;
   (*expr)->curvature = SCIP_EXPRCURV_UNKNOWN;
   (*expr)->auxfilterpos = -1;

   /* initialize an empty interval for interval evaluation */
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &(*expr)->activity);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureConsExprExpr((*expr)->children[c]);
   }

   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** frees an expression */
static
SCIP_RETCODE freeExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to free the expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(*expr != NULL);
   assert((*expr)->nuses == 1);

   /* free children array, if any */
   SCIPfreeBlockMemoryArrayNull(scip, &(*expr)->children, (*expr)->childrensize);

   /* expression should not be locked anymore */
   assert((*expr)->nlockspos == 0);
   assert((*expr)->nlocksneg == 0);

   SCIPfreeBlockMemory(scip, expr);
   assert(*expr == NULL);

   return SCIP_OKAY;
}

/** frees auxiliary variables of expression, if any */
static
SCIP_RETCODE freeAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler, can be NULL */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which auxvar to free, if any */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);

   if( expr->auxvar == NULL )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "remove auxiliary variable %s for expression %p\n", SCIPvarGetName(expr->auxvar), (void*)expr);

   /* check that if not finishing up, noone else is still using the auxvar
    * once we release it, noone else would take care of unlocking it
    * note that we do not free auxvars in exitsolve if we are restarting
    * TODO: this doesn't work when run from unittests, so should find a safer way to define exceptions than checking the stage
    */
   /* assert(SCIPvarGetNUses(expr->auxvar) == 2 || SCIPgetStage(scip) >= SCIP_STAGE_EXITSOLVE); */

   /* remove variable locks; we assume that no other plugin is still using the auxvar for deducing any type of
    * reductions or cutting planes; unfortunately, we cannot check the auxvar->nuses here because the auxiliary
    * variable might be still captured by SCIP's NLP or some other plugin
    */
   SCIP_CALL( SCIPaddVarLocks(scip, expr->auxvar, -1, -1) );

   if( expr->auxfilterpos >= 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      if( conshdlr == NULL )
         conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPdropVarEvent(scip, expr->auxvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)expr, expr->auxfilterpos) );
      expr->auxfilterpos = -1;
   }

   /* release auxiliary variable */
   SCIP_CALL( SCIPreleaseVar(scip, &expr->auxvar) );
   assert(expr->auxvar == NULL);

   return SCIP_OKAY;
}

/** frees data used for enforcement, that is, nonlinear handlers and auxiliary variables */
static
SCIP_RETCODE freeEnfoData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler, can be NULL */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression whose enforcement data will be released */
   SCIP_Bool             freeauxvar          /**< whether aux var should be released */
   )
{
   int e;

   /* free auxiliary variable */
   if( freeauxvar )
   {
      SCIP_CALL( freeAuxVar(scip, conshdlr, expr) );
      assert(expr->auxvar == NULL);
   }

   /* free data stored by nonlinear handlers */
   for( e = 0; e < expr->nenfos; ++e )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr;

      assert(expr->enfos[e] != NULL);

      nlhdlr = expr->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      if( expr->enfos[e]->issepainit )
      {
         /* call the separation deinitialization callback of the nonlinear handler */
         SCIP_CALL( SCIPexitsepaConsExprNlhdlr(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata) );
         expr->enfos[e]->issepainit = FALSE;
      }

      /* free nlhdlr exprdata, if there is any and there is a method to free this data */
      if( expr->enfos[e]->nlhdlrexprdata != NULL && nlhdlr->freeexprdata != NULL )
      {
         SCIP_CALL( (*nlhdlr->freeexprdata)(scip, nlhdlr, expr, &expr->enfos[e]->nlhdlrexprdata) );
         assert(expr->enfos[e]->nlhdlrexprdata == NULL);
      }

      /* free enfo data */
      SCIPfreeBlockMemory(scip, &expr->enfos[e]); /*lint !e866 */
   }

   /* free array with enfo data */
   SCIPfreeBlockMemoryArrayNull(scip, &expr->enfos, expr->nenfos);
   expr->nenfos = 0;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_MAPVAR(transformVar)
{   /*lint --e{715}*/
   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(sourcescip == targetscip);

   /* transform variable (does not capture target variable) */
   SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcevar, targetvar) );
   assert(*targetvar != NULL);

   /* caller assumes that target variable has been captured */
   SCIP_CALL( SCIPcaptureVar(sourcescip, *targetvar) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_MAPVAR(copyVar)
{
   COPY_MAPVAR_DATA* data;
   SCIP_Bool valid;

   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(mapvardata != NULL);

   data = (COPY_MAPVAR_DATA*)mapvardata;

   SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevar, targetvar, data->varmap, data->consmap, data->global, &valid) );
   assert(*targetvar != NULL);

   /* if copy was not valid, store so in mapvar data */
   if( !valid )
      data->valid = FALSE;

   /* caller assumes that target variable has been captured */
   SCIP_CALL( SCIPcaptureVar(targetscip, *targetvar) );

   return SCIP_OKAY;
}

/** copies an expression including subexpressions
 *
 * @note If copying fails due to an expression handler not being available in the targetscip, then *targetexpr will be set to NULL.
 */
static
SCIP_RETCODE copyExpr(
   SCIP*                 sourcescip,         /**< SCIP data structure corresponding to source expression */
   SCIP*                 targetscip,         /**< SCIP data structure where target expression will live */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   sourceexpr,         /**< expression to be copied */
   SCIP_CONSEXPR_EXPR**  targetexpr,         /**< buffer to store pointer to copy of source expression */
   SCIP_DECL_CONSEXPR_MAPVAR((*mapvar)),  /**< variable mapping function, or NULL for identity mapping */
   void*                 mapvardata          /**< data of variable mapping function */
   )
{
   SCIP_CONSHDLR* targetconsexprhdlr = NULL;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPRITERATOR_USERDATA expriteruserdata;
   SCIP_CONSEXPR_EXPR* expr;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(consexprhdlr != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);

   if( targetscip != sourcescip )
      targetconsexprhdlr = SCIPfindConshdlr(targetscip, CONSHDLR_NAME);
   else
      targetconsexprhdlr = consexprhdlr;
   assert(targetconsexprhdlr != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, consexprhdlr, SCIPblkmem(sourcescip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, sourceexpr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );  /*TODO use FALSE, i.e., don't duplicate common subexpr? */
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR | SCIP_CONSEXPRITERATOR_VISITEDCHILD);

   expr = sourceexpr;
   while( !SCIPexpriteratorIsEnd(it) )
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_ENTEREXPR :
         {
            /* create expr that will hold the copy */
            SCIP_CONSEXPR_EXPR* exprcopy;

            /* if the source is a variable expression create a variable expression directly; otherwise copy the expression data */
            if( SCIPisConsExprExprVar(expr) )
            {
               SCIP_VAR* sourcevar;
               SCIP_VAR* targetvar;

               sourcevar = SCIPgetConsExprExprVarVar(expr);
               assert(sourcevar != NULL);
               targetvar = NULL;

               /* get the corresponding variable in the target SCIP */
               if( mapvar != NULL )
               {
                  SCIP_CALL( mapvar(targetscip, &targetvar, sourcescip, sourcevar, mapvardata) );
                  SCIP_CALL( SCIPcreateConsExprExprVar(targetscip, targetconsexprhdlr, &exprcopy, targetvar) );

                  /* we need to release once since it has been captured by the mapvar() and SCIPcreateConsExprExprVar() call */
                  SCIP_CALL( SCIPreleaseVar(targetscip, &targetvar) );
               }
               else
               {
                  targetvar = sourcevar;
                  SCIP_CALL( SCIPcreateConsExprExprVar(targetscip, targetconsexprhdlr, &exprcopy, targetvar) );
               }
            }
            else
            {
               SCIP_CONSEXPR_EXPRHDLR* targetexprhdlr;
               SCIP_CONSEXPR_EXPRDATA* targetexprdata;

               /* get the exprhdlr of the target scip */
               if( targetscip != sourcescip )
               {
                  assert(targetconsexprhdlr != NULL);

                  targetexprhdlr = SCIPfindConsExprExprHdlr(targetconsexprhdlr,
                     SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));

                  if( targetexprhdlr == NULL )
                  {
                     /* expression handler not in target scip (probably did not have a copy callback) -> abort */
                     expriteruserdata.ptrval = NULL;
                     SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

                     expr = SCIPexpriteratorSkipDFS(it);
                     continue;
                  }
               }
               else
               {
                  targetexprhdlr = SCIPgetConsExprExprHdlr(expr);
               }
               assert(targetexprhdlr != NULL);

               /* copy expression data */
               if( expr->exprdata != NULL )
               {
                  assert(expr->exprhdlr->copydata != NULL);
                  SCIP_CALL( expr->exprhdlr->copydata(targetscip, targetexprhdlr, &targetexprdata, sourcescip, expr, mapvar, mapvardata) );
               }
               else
               {
                  targetexprdata = NULL;
               }

               /* create in targetexpr an expression of the same type as expr, but without children for now */
               SCIP_CALL( SCIPcreateConsExprExpr(targetscip, &exprcopy, targetexprhdlr, targetexprdata, 0, NULL) );
            }

            /* store targetexpr */
            expriteruserdata.ptrval = exprcopy;
            SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITEDCHILD :
         {
            /* just visited child so a copy of himself should be available; append it */
            SCIP_CONSEXPR_EXPR* exprcopy;
            SCIP_CONSEXPR_EXPR* childcopy;

            exprcopy = (SCIP_CONSEXPR_EXPR*)SCIPexpriteratorGetCurrentUserData(it).ptrval;

            /* get copy of child */
            childcopy = (SCIP_CONSEXPR_EXPR*)SCIPexpriteratorGetChildUserDataDFS(it).ptrval;
            if( childcopy == NULL )
            {
               /* abort */
               /* release exprcopy (should free also the already copied children) */
               SCIP_CALL( SCIPreleaseConsExprExpr(targetscip, (SCIP_CONSEXPR_EXPR**)&exprcopy) );

               expriteruserdata.ptrval = NULL;
               SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }

            /* append child to exprcopy */
            SCIP_CALL( SCIPappendConsExprExpr(targetscip, exprcopy, childcopy) );

            /* release childcopy (still captured by exprcopy) */
            SCIP_CALL( SCIPreleaseConsExprExpr(targetscip, &childcopy) );

            break;
         }

         default:
            /* we should never be called in this stage */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

   /* the target expression should be stored in the userdata of the sourceexpr (can be NULL if aborted) */
   *targetexpr = (SCIP_CONSEXPR_EXPR*)SCIPexpriteratorGetExprUserData(it, sourceexpr).ptrval;

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(SCIP* scip);

/** copy expression and nonlinear handlers from sourceconshdlr to (target's) scip consexprhdlr */
static
SCIP_RETCODE copyConshdlrExprExprHdlr(
   SCIP*                 scip,               /**< (target) SCIP data structure */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint expression handler */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   int                i;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLRDATA* sourceconshdlrdata;

   assert(strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   assert(conshdlr != sourceconshdlr);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   sourceconshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert(sourceconshdlrdata != NULL);

   /* copy expression handlers */
   *valid = TRUE;
   for( i = 0; i < sourceconshdlrdata->nexprhdlrs; i++ )
   {
      SCIP_Bool localvalid;
      SCIP_CONSEXPR_EXPRHDLR* sourceexprhdlr;

      sourceexprhdlr = sourceconshdlrdata->exprhdlrs[i];

      if( sourceexprhdlr->copyhdlr != NULL )
      {
         SCIP_CALL( sourceexprhdlr->copyhdlr(scip, conshdlr, sourceconshdlr, sourceexprhdlr, &localvalid) );
         *valid &= localvalid;
      }
      else
      {
         *valid = FALSE;
      }
   }

   /* set pointer to important expression handlers in conshdlr of target SCIP */
   conshdlrdata->exprvarhdlr = SCIPfindConsExprExprHdlr(conshdlr, "var");
   conshdlrdata->exprvalhdlr = SCIPfindConsExprExprHdlr(conshdlr, "val");
   conshdlrdata->exprsumhdlr = SCIPfindConsExprExprHdlr(conshdlr, "sum");
   conshdlrdata->exprprodhdlr = SCIPfindConsExprExprHdlr(conshdlr, "prod");
   conshdlrdata->exprpowhdlr = SCIPfindConsExprExprHdlr(conshdlr, "pow");
   conshdlrdata->exprsignpowhdlr = SCIPfindConsExprExprHdlr(conshdlr, "signpower");
   conshdlrdata->exprexphdlr = SCIPfindConsExprExprHdlr(conshdlr, "exp");

   /* copy nonlinear handlers */
   for( i = 0; i < sourceconshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CONSEXPR_NLHDLR* sourcenlhdlr;

      /* TODO for now just don't copy disabled nlhdlr, a clean way would probably be to first copy and disable then */
      sourcenlhdlr = sourceconshdlrdata->nlhdlrs[i];
      if( sourcenlhdlr->copyhdlr != NULL && sourcenlhdlr->enabled )
      {
         SCIP_CALL( sourcenlhdlr->copyhdlr(scip, conshdlr, sourceconshdlr, sourcenlhdlr) );
      }
   }

   return SCIP_OKAY;
}

/** tries to automatically convert an expression constraint into a more specific and more specialized constraint */
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
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;
   int i;

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

   /* if there are no upgrade methods, we can stop */
   if( conshdlrdata->nexprconsupgrades == 0 )
      return SCIP_OKAY;

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   /* call the upgrading methods */
   SCIPdebugMsg(scip, "upgrading expression constraint <%s> (up to %d upgrade methods): ",
      SCIPconsGetName(cons), conshdlrdata->nexprconsupgrades);
   SCIPdebugPrintCons(scip, cons, NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nexprconsupgrades; ++i )
   {
      if( !conshdlrdata->exprconsupgrades[i]->active )
         continue;

      assert(conshdlrdata->exprconsupgrades[i]->exprconsupgd != NULL);

      SCIP_CALL( conshdlrdata->exprconsupgrades[i]->exprconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->exprconsupgrades[i]->exprconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         int j;

         SCIPdebugMsg(scip, " -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugMsgPrint(scip, "\t");
            SCIPdebugPrintCons(scip, upgdconss[j], NULL);

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMsg(scip, "delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** interval evaluation of variables as used in bound tightening
 *
 * Returns slightly relaxed local variable bounds of a variable as interval.
 * Does not relax beyond integer values, thus does not relax bounds on integer variables at all.
 */
static
SCIP_DECL_CONSEXPR_INTEVALVAR(intEvalVarBoundTightening)
{
   SCIP_INTERVAL interval;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(lb <= ub);  /* can SCIP ensure by now that variable bounds are not contradicting? */

   /* implicit integer variables may have non-integer bounds, apparently (run space25a) */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
   {
      lb = EPSROUND(lb, 0.0); /*lint !e835*/
      ub = EPSROUND(ub, 0.0); /*lint !e835*/
   }

   /* integer variables should always have integral bounds in SCIP */
   assert(EPSFRAC(lb, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/
   assert(EPSFRAC(ub, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/

   switch( conshdlrdata->varboundrelax )
   {
      case 'n' : /* no relaxation */
         break;

      case 'a' : /* relax by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
         {
            /* reduce lb by epsilon, or to the next integer value, which ever is larger */
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - conshdlrdata->varboundrelaxamount);
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            /* increase ub by epsilon, or to the next integer value, which ever is smaller */
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + conshdlrdata->varboundrelaxamount);
         }

         break;
      }

      case 'b' : /* relax always by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
            lb -= conshdlrdata->varboundrelaxamount;

         if( !SCIPisInfinity(scip, ub) )
            ub += conshdlrdata->varboundrelaxamount;

         break;
      }

      case 'r' : /* relax by relative value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         /* relax bounds by epsilon*max(1,|bnd|), instead of just epsilon as in case 'a', thus we trust the first log(epsilon) digits
          * however, when domains get small, relaxing can excessively weaken bound tightening, thus do only fraction of |ub-lb| if that is smaller
          * further, do not relax beyond next integer value
          */
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(lb)), 0.001 * REALABS(ub-lb)));  /*lint !e666*/
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(ub)), 0.001 * REALABS(ub-lb)));  /*lint !e666*/
         }

         break;
      }

      default :
      {
         SCIPerrorMessage("Unsupported value '%c' for varboundrelax option.\n");
         SCIPABORT();
         break;
      }
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** interval evaluation of variables as used in redundancy check
 *
 * Returns local variable bounds of a variable, relaxed by feastol, as interval.
 */
static
SCIP_DECL_CONSEXPR_INTEVALVAR(intEvalVarRedundancyCheck)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(lb <= ub);  /* can SCIP ensure by now that variable bounds are not contradicting? */

   /* relax variable bounds, if there are bounds and variable is not fixed
    * (actually some assert complains if trying SCIPisRelEQ if both bounds are at different infinity)
    */
   if( !(SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub)) && !SCIPisRelEQ(scip, lb, ub) )
   {
      if( !SCIPisInfinity(scip, -lb) )
         lb -= SCIPfeastol(scip);

      if( !SCIPisInfinity(scip, ub) )
         ub += SCIPfeastol(scip);
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** returns whether intersecting oldinterval with newinterval would provide a properly smaller interval
 *
 * If subsetsufficient is TRUE, then the intersection being smaller than oldinterval is sufficient.
 * If subsetsufficient is FALSE, then we require
 *  - a change from an unbounded interval to a bounded one, or
 *  - or a change from an unfixed (width > epsilon) to a fixed interval, or
 *  - a minimal tightening of one of the interval bounds as defined by SCIPis{Lb,Ub}Better.
 */
static
SCIP_Bool isIntervalBetter(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_Bool               subsetsufficient, /**< whether the intersection being a proper subset of oldinterval is sufficient */
   SCIP_INTERVAL           newinterval,      /**< new interval */
   SCIP_INTERVAL           oldinterval       /**< old interval */
   )
{
   assert(scip != NULL);
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newinterval));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, oldinterval));

   if( subsetsufficient )
      /* oldinterval \cap newinterval < oldinterval iff not oldinterval is subset of newinterval */
      return !SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, oldinterval, newinterval);

   /* check whether lower bound of interval becomes finite */
   if( oldinterval.inf <= -SCIP_INTERVAL_INFINITY && newinterval.inf > -SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether upper bound of interval becomes finite */
   if( oldinterval.sup >=  SCIP_INTERVAL_INFINITY && newinterval.sup >  SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether intersection will have width <= epsilon, if oldinterval doesn't have yet */
   if( !SCIPisEQ(scip, oldinterval.inf, oldinterval.sup) && SCIPisEQ(scip, MAX(oldinterval.inf, newinterval.inf), MIN(oldinterval.sup, newinterval.sup)) ) /*lint !e666*/
      return TRUE;

   /* check whether lower bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisLbBetter(scip, newinterval.inf, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   /* check whether upper bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisUbBetter(scip, newinterval.sup, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   return FALSE;
}


/** tightens the bounds of the auxiliary variable associated with an expression (or original variable if being a variable-expression) according to its activity
 *
 *  Nothing will happen if SCIP is not in presolve or solve.
 */
static
SCIP_RETCODE tightenAuxVarBounds(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be tightened */
   SCIP_Bool               force,            /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a cutoff was detected */
   int*                    ntightenings      /**< buffer to add the total number of tightenings, or NULL */
   )
{
   SCIP_VAR* var;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   /* the current activity must be valid and should not be empty */
   assert(expr->activitytag >= SCIPconshdlrGetData(conshdlr)->lastboundrelax || SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, expr->activity));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity));

   *cutoff = FALSE;

   /* do not tighten variable in problem stage (important for unittests)
    * TODO put some kind of #ifdef UNITTEST around this once the unittest are modified to include the .c file (again)?
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   var = SCIPgetConsExprExprAuxVar(expr);
   if( var == NULL )
      return SCIP_OKAY;

   /* force tightening if it would mean fixing the variable */
   force = force || SCIPisEQ(scip, expr->activity.inf, expr->activity.sup);

   /* try to tighten lower bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarLb(scip, var, expr->activity.inf, force, cutoff, &tightened) );
   if( tightened )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened lb on auxvar <%s> to %.15g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var));
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening lb on auxvar <%s> to %.15g\n", SCIPvarGetName(var), expr->activity.inf);
      return SCIP_OKAY;
   }

   /* try to tighten upper bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarUb(scip, var, expr->activity.sup, force, cutoff, &tightened) );
   if( tightened )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened ub on auxvar <%s> to %.15g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening ub on auxvar <%s> to %.15g\n", SCIPvarGetName(var), expr->activity.sup);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** propagate bounds of the expressions in a given expression tree and tries to tighten the bounds of the auxiliary
 *  variables accordingly
 */
static
SCIP_RETCODE forwardPropExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression */
   SCIP_Bool               force,            /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool               tightenauxvars,   /**< should the bounds of auxiliary variables be tightened? */
   SCIP_DECL_CONSEXPR_INTEVALVAR((*intevalvar)), /**< function to call to evaluate interval of variable, or NULL */
   void*                   intevalvardata,   /**< data to be passed to intevalvar call */
   SCIP_QUEUE*             reversepropqueue, /**< queue to add candidates for reverse propagation, or NULL */
   SCIP_Bool*              infeasible,       /**< buffer to store whether the problem is infeasible (NULL if not needed) */
   int*                    ntightenings      /**< buffer to store the number of auxiliary variable tightenings (NULL if not needed) */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(rootexpr != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;
   if( ntightenings != NULL )
      *ntightenings = 0;

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* if value is valid and empty, then we cannot improve, so do nothing */
   if( rootexpr->activitytag >= conshdlrdata->lastboundrelax && SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rootexpr->activity) )
   {
      SCIPdebugMsg(scip, "stored activity of root expr is empty and valid (activitytag >= lastboundrelax (%u)), skip forwardPropExpr -> cutoff\n", conshdlrdata->lastboundrelax);

      if( infeasible != NULL )
         *infeasible = TRUE;

      return SCIP_OKAY;
   }

   /* if value is up-to-date, then nothing to do */
   if( rootexpr->activitytag == conshdlrdata->curboundstag )
   {
      SCIPdebugMsg(scip, "activitytag of root expr equals curboundstag (%u), skip forwardPropExpr\n", conshdlrdata->curboundstag);

      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rootexpr->activity)); /* handled in previous if() */

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexpriteratorCreate(&it, consexprhdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, rootexpr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it);  )
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
         {
            /* skip child if it has been evaluated already */
            SCIP_CONSEXPR_EXPR* child;

            child = SCIPexpriteratorGetChildExprDFS(it);
            if( conshdlrdata->curboundstag == child->activitytag )
            {
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, child->activity) )
               {
                  if( infeasible != NULL )
                     *infeasible = TRUE;
               }

#ifndef NDEBUG
               /* We do not check here whether the child should be added to the reversepropqueue.
                * This should have happened when the activitytag of the child was set to curboundstag, I believe.
                * So we can assert that the child activity should be a subset of the auxiliary variable bounds.
                * Since SCIP sometimes moves variable bounds slightly and we also use relaxed variable bounds below,
                * we have to add some epsilons here.
                */
               else if( child->auxvar != NULL )
               {
                  SCIP_INTERVAL auxvarbounds;
                  auxvarbounds = intevalvar(scip, child->auxvar, intevalvardata);
                  assert(reversepropqueue == NULL || child->inqueue ||
                     ((auxvarbounds.inf <= -SCIP_INTERVAL_INFINITY || SCIPisRelGE(scip, child->activity.inf, auxvarbounds.inf)) &&
                      (auxvarbounds.sup >=  SCIP_INTERVAL_INFINITY || SCIPisRelLE(scip, child->activity.sup, auxvarbounds.sup))));
               }
#endif

               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_CONSEXPRITERATOR_LEAVEEXPR :
         {
            SCIP_INTERVAL prevactivity;
            SCIP_INTERVAL auxvarbounds = { -SCIP_DEFAULT_INFINITY, SCIP_DEFAULT_INFINITY };  /* init just for lint and scan-build */

            /* we should not have entered this expression if its activity was already uptodate */
            assert(expr->activitytag < conshdlrdata->curboundstag);

            if( expr->activitytag < conshdlrdata->lastboundrelax )
            {
               /* reset activity to infinity if invalid, because SCIPtightenConsExprExprInterval seems to assume valid activity in expr */
               SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &expr->activity);
            }
            else if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               /* if already empty, then don't try to compute even better activity
                * we should have noted that we are infeasible, though (if not remove the assert and enable below code)
                */
               assert(infeasible == NULL || *infeasible);
               /* if( infeasible != NULL )
                  *infeasible = TRUE; */
               break;
            }

#ifdef DEBUG_PROP
            SCIPdebugMsg(scip, "interval evaluation of expr %p ", (void*)expr);
            SCIP_CALL( SCIPprintConsExprExpr(scip, consexprhdlr, expr, NULL) );
            SCIPdebugMsgPrint(scip, ", current activity = [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif

            /* store activity that we start with */
            prevactivity = expr->activity;

            /* reset existing activity of expression if we are collecting expressions for reverse propagation
             * the reason for this is that expr->activity might currently store bounds from the previous
             * reverse propagation and we want to collect those expressions for the next reverse propagation
             * where the forward propagation does not already provide as good activity as those given by
             * previous reverse propagation (i.e., expressions where there is potential for reverse propagation
             * because we know tighter bounds on the expression than what is given by forward propagation)
             */
            if( reversepropqueue != NULL )
               SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &expr->activity);

            /* run interval eval of nonlinear handlers or expression handler */
            if( expr->nenfos > 0 )
            {
               SCIP_CONSEXPR_NLHDLR* nlhdlr;
               SCIP_INTERVAL nlhdlrinterval;
               int e;

               /* for nodes with enforcement (having auxvar, thus during solve), nlhdlrs take care of interval evaluation */
               for( e = 0; e < expr->nenfos && !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity); ++e )
               {
                  nlhdlr = expr->enfos[e]->nlhdlr;
                  assert(nlhdlr != NULL);

                  /* skip nlhdlr if it does not provide interval evaluation */
                  if( !SCIPhasConsExprNlhdlrInteval(nlhdlr) )
                     continue;

                  /* let nlhdlr evaluate current expression */
                  nlhdlrinterval = expr->activity;
                  SCIP_CALL( SCIPintevalConsExprNlhdlr(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &nlhdlrinterval, intevalvar, intevalvardata) );
#ifdef DEBUG_PROP
                  SCIPdebugMsg(scip, " nlhdlr <%s>::inteval = [%.20g, %.20g]", nlhdlr->name, nlhdlrinterval.inf, nlhdlrinterval.sup);
#endif

                  /* update expr->activity by intersecting with computed activity */
                  SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, nlhdlrinterval);
#ifdef DEBUG_PROP
                  SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
               }
            }
            else
            {
               /* for node without enforcement (no auxvar, maybe in presolve), call the callback of the exprhdlr directly */
               SCIP_INTERVAL exprhdlrinterval = expr->activity;
               SCIP_CALL( SCIPintevalConsExprExprHdlr(scip, expr, &exprhdlrinterval, intevalvar, intevalvardata) );
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " exprhdlr <%s>::inteval = [%.20g, %.20g]", expr->exprhdlr->name, exprhdlrinterval.inf, exprhdlrinterval.sup);
#endif

               /* update expr->activity by intersecting with computed activity */
               SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, exprhdlrinterval);
#ifdef DEBUG_PROP
               SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
            }

            /* if expression is integral, then we try to tighten the interval bounds a bit
             * this should undo the addition of some unnecessary safety added by use of nextafter() in interval arithmetics, e.g., when doing pow()
             * it would be ok to use ceil() and floor(), but for safety we use SCIPceil and SCIPfloor for now
             */
            if( expr->isintegral )
            {
               if( expr->activity.inf > -SCIP_INTERVAL_INFINITY )
                  expr->activity.inf = SCIPceil(scip, expr->activity.inf);
               if( expr->activity.sup <  SCIP_INTERVAL_INFINITY )
                  expr->activity.sup = SCIPfloor(scip, expr->activity.sup);
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " applying integrality: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
            }

            /* mark the current node to be infeasible if either the lower/upper bound is above/below +/- SCIPinfinity() */
            if( SCIPisInfinity(scip, expr->activity.inf) || SCIPisInfinity(scip, -expr->activity.sup) )
            {
               SCIPdebugMsg(scip, "cut off due to activity [%g,%g] beyond infinity\n", expr->activity.inf, expr->activity.sup);
               SCIPintervalSetEmpty(&expr->activity);
            }

            /* get interval with auxiliary variable bounds */
            if( expr->auxvar != NULL )
            {
               auxvarbounds = intevalvar(scip, expr->auxvar, intevalvardata);
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " auxvar <%s> bounds: [%.15g,%.15g]\n", SCIPvarGetName(expr->auxvar), auxvarbounds.inf, auxvarbounds.sup);
#endif

               /* it would be odd if the domain of an auxiliary variable were empty */
               assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, auxvarbounds));
            }

            if( reversepropqueue != NULL && !expr->inqueue && !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               /* check whether previously valid activity (maybe set by reverse-propagation) or auxiliary variable bounds (if any)
                * would result in a tightening; if so, add expr to reversepropqueue
                */
               SCIP_INTERVAL compareinterval;

               if( expr->auxvar != NULL )
                  SCIPintervalIntersectEps(&compareinterval, SCIPepsilon(scip), prevactivity, auxvarbounds);
               else
                  compareinterval = prevactivity;

               /* if compareinterval allow a further tightening, then it may provide tighter bounds for children
                * thus add this expression to the reversepropqueue
                */
               if( !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, compareinterval) && isIntervalBetter(scip, force, compareinterval, expr->activity) )
               {
#ifdef DEBUG_PROP
                  SCIPdebugMsg(scip, " insert expr <%p> (%s) into reversepropqueue, new activity = [%.15g,%.15g] is not subset of previous one = [%.15g,%.15g]\n", (void*)expr, SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), expr->activity.inf, expr->activity.sup, compareinterval.inf, compareinterval.sup);
#endif
                  SCIP_CALL( SCIPqueueInsert(reversepropqueue, expr) );
                  expr->inqueue = TRUE;
               }
            }

            if( !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               /* tighten activity from inteval with prevactivity and auxiliary variable bounds (if any) */
               SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, prevactivity);
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " apply previous activity = [%.20g, %.20g]: [%.20g, %.20g]\n", prevactivity.inf, prevactivity.sup, expr->activity.inf, expr->activity.sup);
#endif
               if( expr->auxvar != NULL )
               {
                  SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, auxvarbounds);
#ifdef DEBUG_PROP
                  SCIPdebugMsg(scip, " apply auxvar <%s> bounds = [%.20g, %.20g]: [%.20g, %.20g]\n", SCIPvarGetName(expr->auxvar), auxvarbounds.inf, auxvarbounds.sup);
#endif
               }
            }

            /* remember that activity is uptodate now */
            expr->activitytag = conshdlrdata->curboundstag;

            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               if( infeasible != NULL )
                  *infeasible = TRUE;
            }
            else if( tightenauxvars )
            {
               SCIP_Bool tighteninfeasible;

               SCIP_CALL( tightenAuxVarBounds(scip, consexprhdlr, expr, force, &tighteninfeasible, ntightenings) );
               if( tighteninfeasible )
               {
                  if( infeasible != NULL )
                     *infeasible = TRUE;
                  SCIPintervalSetEmpty(&expr->activity);
               }
            }

            break;
         }

         default:
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** propagates bounds for each sub-expression in a given queue by starting from the root expressions
 *
 *  the expression will be traversed in breadth first search by using this queue
 *
 *  @note calling this function requires feasible intervals for each sub-expression; this is guaranteed by calling
 *  forwardPropExpr() before calling this function
 */
static
SCIP_RETCODE reversePropQueue(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_QUEUE*             queue,            /**< queue of expression to propagate */
   SCIP_Bool               force,            /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool               allexprs,         /**< whether reverseprop should be called for all expressions, regardless of whether their interval was tightened */
   SCIP_Bool*              infeasible,       /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   int*                    ntightenings      /**< buffer to store the number of (variable) tightenings */
   )
{
   assert(queue != NULL);
   assert(infeasible != NULL);
   assert(ntightenings != NULL);

   *infeasible = FALSE;
   *ntightenings = 0;

   /* main loop that calls reverse propagation for expressions on the queue
    * when reverseprop finds a tightening for an expression, then that expression is added to the queue (within the reverseprop call)
    */
   while( !SCIPqueueIsEmpty(queue) && !(*infeasible) )
   {
      SCIP_CONSEXPR_EXPR* expr;
      int e;

      expr = (SCIP_CONSEXPR_EXPR*) SCIPqueueRemove(queue);
      assert(expr != NULL);

      /* mark that the expression is not in the queue anymore */
      expr->inqueue = FALSE;

      if( expr->nenfos > 0 )
      {
         /* for nodes with enforcement, call reverse propagation callbacks of nlhdlrs */
         for( e = 0; e < expr->nenfos && !*infeasible; ++e )
         {
            SCIP_CONSEXPR_NLHDLR* nlhdlr;
            int nreds;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* call the reverseprop of the nlhdlr */
#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, "call reverse propagation for ");
            SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
            SCIPdebugMsgPrint(scip, " in [%g,%g] using nlhdlr <%s>\n", expr->activity.inf, expr->activity.sup, nlhdlr->name);
#endif

            nreds = 0;
            SCIP_CALL( SCIPreversepropConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, queue, infeasible, &nreds, force) );
            assert(nreds >= 0);
            *ntightenings += nreds;
         }
      }
      else
      {
         /* if node without enforcement (no auxvar or in presolve), call reverse propagation callback of exprhdlr directly */
         int nreds = 0;

#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "call reverse propagation for ");
         SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
         SCIPdebugMsgPrint(scip, " in [%g,%g] using exprhdlr <%s>\n", expr->activity.inf, expr->activity.sup, expr->exprhdlr->name);
#endif

         /* call the reverseprop of the exprhdlr */
         SCIP_CALL( SCIPreversepropConsExprExprHdlr(scip, conshdlr, expr, queue, infeasible, &nreds, force) );
         assert(nreds >= 0);
         *ntightenings += nreds;
      }

      /* if allexprs is set, then make sure that all children of expr with children are in the queue
       * SCIPtightenConsExprExpr only adds children to the queue which have reverseprop capability
       */
      if( allexprs )
      {
         int i;
         for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
         {
            SCIP_CONSEXPR_EXPR* child;

            child = SCIPgetConsExprExprChildren(expr)[i];

            if( !child->inqueue && SCIPgetConsExprExprNChildren(child) > 0 )
            {
               /* SCIPdebugMsg(scip, "allexprs: insert expr <%p> (%s) into reversepropqueue\n", (void*)child, SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))); */
               SCIP_CALL( SCIPqueueInsert(queue, (void*) child) );
               child->inqueue = TRUE;
            }
         }
      }

      /* stop propagation if the problem is infeasible */
      if( *infeasible )
         break;
   }

   /* reset expr->inqueue for all remaining expr's in queue (can happen in case of early stop due to infeasibility) */
   while( !SCIPqueueIsEmpty(queue) )
   {
      SCIP_CONSEXPR_EXPR* expr;

      expr = (SCIP_CONSEXPR_EXPR*) SCIPqueueRemove(queue);

      /* mark that the expression is not in the queue anymore */
      expr->inqueue = FALSE;
   }

   return SCIP_OKAY;
}

/** calls domain propagation for a given set of constraints; the algorithm alternates calls of forward and reverse
 *  propagation; the latter only for nodes which have been tightened during the propagation loop;
 *
 *  the propagation algorithm works as follows:
 *
 *   1.) apply forward propagation (update activities) and collect expressions for which auxiliary variables (during solve)
 *       or constraint sides (during presolve) provide tighter bounds
 *
 *   2.) apply reverse propagation to all collected expressions; don't explore
 *       sub-expressions which have not changed since the beginning of the propagation loop
 *
 *   3.) if we have found enough tightenings go to 1.) otherwise leave propagation loop
 *
 *  @note after calling forward propagation for a constraint we mark this constraint as propagated; this flag might be
 *  reset during the reverse propagation when we find a bound tightening of a variable expression contained in the
 *  constraint; resetting this flag is done in the EVENTEXEC callback of the event handler
 *
 *  @note when using forward and reverse propagation alternatingly we reuse expression activites computed in previous
 *  iterations
 */
static
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds,            /**< buffer to add the number of changed bounds */
   int*                  ndelconss           /**< buffer to add the number of deleted constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool success = FALSE;
   SCIP_Bool allexprs;
   int ntightenings;
   int roundnr;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);
   assert(ndelconss != NULL);

   /* no constraints to propagate */
   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* TODO maybe only do this if first call or simplify or someone else changed the expression graph */
   allexprs = (SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;
   cutoff = FALSE;

   /* main propagation loop */
   do
   {
      SCIP_QUEUE* queue;

      SCIPdebugMsg(scip, "start propagation round %d\n", roundnr);

      /* create queue */
      SCIP_CALL( SCIPqueueCreate(&queue, SCIPgetNVars(scip), 2.0) );

      /* apply forward propagation (update expression activities)
       * and add promising root expressions into queue for reversepropagation
       */
      for( i = 0; i < nconss; ++i )
      {
         consdata = SCIPconsGetData(conss[i]);
         assert(consdata != NULL);

         /* skip deleted, non-active, or propagation-disabled constraints */
         if( SCIPconsIsDeleted(conss[i]) || !SCIPconsIsActive(conss[i]) )
            continue;

         /* in the first round, we reevaluate all bounds to remove some possible leftovers that could be in this
          * expression from a reverse propagation in a previous propagation round
          * (TODO: do we still need this since we have the tag's???
          * this means that we propagate all constraints even if there was only very few boundchanges that related to only a few constraints)
          * in other rounds, we skip already propagated constraints
          */
         if( (consdata->ispropagated && roundnr > 0) || !SCIPconsIsPropagationEnabled(conss[i]) )
            continue;

         /* update activities in expression and collect initial candidates for reverse propagation */
         SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s> (round %d): ", SCIPconsGetName(conss[i]), roundnr);
         SCIPdebugPrintCons(scip, conss[i], NULL);

         ntightenings = 0;
         SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, force, TRUE, intEvalVarBoundTightening, (void*)SCIPconshdlrGetData(conshdlr), allexprs ? NULL : queue, &cutoff, &ntightenings) );
         assert(cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, consdata->expr->activity));

#ifdef SCIP_DEBUG
         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> found empty bound for an expression during forward propagation of constraint %s\n",
               SCIPconsGetName(conss[i]));
         }
#endif

         if( !cutoff && consdata->expr->auxvar == NULL )
         {
            /* intersect activity with constraint sides (relaxed by epsilon)
             * if we have auxvar (not in presolve), then bounds of the auxvar are initially set to [lhs,rhs], so doing this again is useless
             * otherwise, SCIPtightenConsExprExprInterval will take care of adding the constraint expr to the queue if the sides implied a tightening
             */
            SCIP_INTERVAL conssides;

            /* relax sides by SCIPepsilon() and handle infinite sides */
            SCIP_Real lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - conshdlrdata->conssiderelaxamount;
            SCIP_Real rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + conshdlrdata->conssiderelaxamount;
            SCIPintervalSetBounds(&conssides, lhs, rhs);

            SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, consdata->expr, conssides, force, allexprs ? NULL : queue, &cutoff, &ntightenings) );

            if( cutoff )
            {
               SCIPdebugMsg(scip, " -> cutoff after intersect with conssides\n");
               break;
            }
         }

         /* mark constraint as propagated; this will be reset via the event system when we find a variable tightening */
         consdata->ispropagated = TRUE;

         assert(ntightenings >= 0);
         *nchgbds += ntightenings;

         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> cutoff\n");
            break;
         }

         if( ntightenings > 0 )
            *result = SCIP_REDUCEDDOM;

         if( allexprs && !consdata->expr->inqueue )
         {
            /* SCIPdebugMsg(scip, "allexprs: insert expr <%p> (%s) into reversepropqueue\n", (void*)consdata->expr, SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(consdata->expr))); */
            SCIP_CALL( SCIPqueueInsert(queue, consdata->expr) );
            consdata->expr->inqueue = TRUE;
         }
      }

      if( !cutoff )
      {
         /* apply backward propagation */
         SCIP_CALL( reversePropQueue(scip, conshdlr, queue, force, allexprs, &cutoff, &ntightenings) );

         /* @todo add parameter for the minimum number of tightenings to trigger a new propagation round */
         success = ntightenings > 0;

         if( nchgbds != NULL )
            *nchgbds += ntightenings;

         if( success )
            *result = SCIP_REDUCEDDOM;
      }
      else
      {
         while( !SCIPqueueIsEmpty(queue) )
         {
            SCIP_CONSEXPR_EXPR* expr;
            expr = (SCIP_CONSEXPR_EXPR*) SCIPqueueRemove(queue);
            expr->inqueue = FALSE;
         }
      }

      assert(SCIPqueueIsEmpty(queue));
      SCIPqueueFree(&queue);

      if( cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* do this only for the first round */
      allexprs = FALSE;

#if 0
      for( i = 0; i < nconss; ++i )
      {
         SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
         SCIPprintConsExprExprDotInit(scip, conshdlr, &dotdata, NULL, SCIP_CONSEXPR_PRINTDOT_ACTIVITY | SCIP_CONSEXPR_PRINTDOT_EXPRSTRING);
         SCIPprintConsExprExprDot(scip, dotdata, SCIPconsGetData(conss[i])->expr);
         SCIPprintConsExprExprDotFinal(scip, &dotdata);
      }
#endif
   }
   while( success && ++roundnr < conshdlrdata->maxproprounds );

   return SCIP_OKAY;
}

/** checks constraints for redundancy
 *
 * Checks whether the activity of constraint functions is a subset of the constraint sides (relaxed by feastol).
 * To compute the activity, we use forwardPropExpr(), but relax variable bounds by feastol, because solutions to be checked
 * might violate variable bounds by up to feastol, too.
 * This is the main reason why the redundancy check is not done in propConss(), which relaxes variable bounds by epsilon only.
 *
 * Also removes constraints of the form lhs <= variable <= rhs.
 *
 * @todo it would be sufficient to check constraints for which we know that they are not currently violated by a valid solution
 *
 * @note This could should not run during solving, because the forwardProp takes the bounds of auxiliary variables into account.
 * For the root expression, these bounds are already set to the constraint sides, so that the activity of every expression
 * would appear as if the constraint is redundant.
 */
static
SCIP_RETCODE checkRedundancyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store whether infeasibility has been identified */
   int*                  ndelconss,          /**< buffer to add the number of deleted constraints */
   int*                  nchgbds             /**< buffer to add the number of variable bound tightenings */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL sides;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgbds != NULL);

   /* no constraints to check */
   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase curboundstag and pretend bounds were relaxed
    * we do this here to trigger a reevaluation of all bounds, since we will relax variable bounds
    * for the redundancy check differently than for domain propagation
    */
   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;

   SCIPdebugMsg(scip, "checking %d constraints for redundancy\n", nconss);

   *cutoff = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      if( !SCIPconsIsActive(conss[i]) || SCIPconsIsDeleted(conss[i]) )
         continue;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* handle constant expressions separately: either the problem is infeasible or the constraint is redundant */
      if( consdata->expr->exprhdlr == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         SCIP_Real value = SCIPgetConsExprExprValueValue(consdata->expr);

         if(  (!SCIPisInfinity(scip, -consdata->lhs) && value < consdata->lhs - SCIPfeastol(scip))
            || (!SCIPisInfinity(scip, consdata->rhs) && value > consdata->rhs + SCIPfeastol(scip)) )
         {
            SCIPdebugMsg(scip, "constant constraint <%s> is infeasible: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);
            *cutoff = TRUE;

            return SCIP_OKAY;
         }

         SCIPdebugMsg(scip, "constant constraint <%s> is redundant: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* handle variable expressions separately: tighten variable bounds to constraint sides, then remove constraint (now redundant) */
      if( consdata->expr->exprhdlr == SCIPgetConsExprExprHdlrVar(conshdlr) )
      {
         SCIP_VAR* var;
         SCIP_Bool tightened;

         var = SCIPgetConsExprExprVarVar(consdata->expr);
         assert(var != NULL);

         SCIPdebugMsg(scip, "variable constraint <%s> can be made redundant: <%s>[%g,%g] in [%g,%g] ", SCIPconsGetName(conss[i]), SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), consdata->lhs, consdata->rhs);

         /* ensure that variable bounds are within constraint sides */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, consdata->lhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               return SCIP_OKAY;
         }

         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, consdata->rhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               return SCIP_OKAY;
         }

         /* delete the (now) redundant constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* reevaluate all bounds to remove some possible leftovers that could be in this
       * expression from a reverse propagation in a previous propagation round
       *
       * we relax variable bounds by feastol here, as solutions that are checked later can also violate
       * variable bounds by up to feastol
       * (relaxing fixed variables seems to be too much, but they would be removed by presolve soon anyway)
       */
      SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s>: ", SCIPconsGetName(conss[i]));
      SCIPdebugPrintCons(scip, conss[i], NULL);

      SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, FALSE, FALSE, intEvalVarRedundancyCheck, NULL, NULL, cutoff, NULL) );
      assert(*cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, consdata->expr->activity));

      /* it is unlikely that we detect infeasibility by doing forward propagation */
      if( *cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         return SCIP_OKAY;
      }

      assert(consdata->expr->activitytag == conshdlrdata->curboundstag);
      activity = consdata->expr->activity;

      /* relax sides by feastol
       * we could accept every solution that violates constraints up to feastol as redundant, so this is the most permissive we can be
       */
      SCIPintervalSetBounds(&sides,
         SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - SCIPfeastol(scip),
         SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + SCIPfeastol(scip));

      if( SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, activity, sides) )
      {
         SCIPdebugMsg(scip, " -> redundant: activity [%g,%g] within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, " -> not redundant: activity [%g,%g] not within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);
   }

   /* make sure bounds are reevaluated again, since we relaxed bounds in a different way */
   ++conshdlrdata->curboundstag;
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;

   return SCIP_OKAY;
}

/** install nlhdlrs in one expression */
static
SCIP_RETCODE detectNlhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression for which to run detection routines */
   SCIP_CONS*            cons,               /**< constraint for which expr == consdata->expr, otherwise NULL */
   SCIP_CONSEXPR_NLHDLR**  nlhdlrssuccess,   /**< buffer for nlhdlrs that had success detecting structure at expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrssuccessexprdata, /**< buffer for exprdata of nlhdlrs */
   SCIP_Bool*            infeasible          /**< buffer to indicate whether infeasibility has been detected */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool enforcedbelow;
   SCIP_Bool enforcedabove;
   SCIP_CONSEXPR_EXPRENFO_METHOD enforcemethods;
   SCIP_Bool nlhdlrenforcedbelow;
   SCIP_Bool nlhdlrenforcedabove;
   SCIP_CONSEXPR_EXPRENFO_METHOD nlhdlrenforcemethods;
   SCIP_Bool success;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;
   int ntightenings;
   int nsuccess;
   int e, h;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrssuccess != NULL);
   assert(nlhdlrssuccessexprdata != NULL);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);

   /* there should be no enforcer yet, i.e., detection should not have been run already */
   assert(expr->nenfos == 0);
   assert(expr->enfos == NULL);

   /* analyze expression with nonlinear handlers
    * if nobody positively (up) locks expr -> only need to enforce expr >= auxvar -> no need for underestimation
    * if nobody negatively (down) locks expr -> only need to enforce expr <= auxvar -> no need for overestimation
    */
   nsuccess = 0;
   enforcemethods = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcedbelow = (SCIPgetConsExprExprNLocksPos(expr) == 0); /* no need for underestimation */
   enforcedabove = (SCIPgetConsExprExprNLocksNeg(expr) == 0); /* no need for overestimation */

   SCIPdebugMsg(scip, "detecting nlhdlrs for %s expression %p (%s); start with below %d above %d\n",
      cons != NULL ? "root" : "non-root", (void*)expr, SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), enforcedbelow, enforcedabove);

   for( h = 0; h < conshdlrdata->nnlhdlrs && !*infeasible; ++h )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[h];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlrs */
      if( !nlhdlr->enabled )
         continue;

      /* call detect routine of nlhdlr */
      nlhdlrexprdata = NULL;
      success = FALSE;
      nlhdlrenforcemethods = enforcemethods;
      nlhdlrenforcedbelow = enforcedbelow;
      nlhdlrenforcedabove = enforcedabove;
      SCIP_CALL( SCIPdetectConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, cons, &nlhdlrenforcemethods, &nlhdlrenforcedbelow, &nlhdlrenforcedabove, &success, &nlhdlrexprdata) );

      /* detection is only allowed to augment to the various parameters (enforce "more", add "more" methods) */
      assert(nlhdlrenforcemethods >= enforcemethods);
      assert(nlhdlrenforcedbelow >= enforcedbelow);
      assert(nlhdlrenforcedabove >= enforcedabove);

      if( !success )
      {
         /* nlhdlrexprdata can only be non-NULL if it provided some functionality */
         assert(nlhdlrexprdata == NULL);
         assert(nlhdlrenforcemethods == enforcemethods);
         assert(nlhdlrenforcedbelow == enforcedbelow);
         assert(nlhdlrenforcedabove == enforcedabove);

         continue;
      }

      SCIPdebugMsg(scip, "nlhdlr <%s> detect successful; now enforced below: %d above: %d methods: %d\n",
         SCIPgetConsExprNlhdlrName(nlhdlr), nlhdlrenforcedbelow, nlhdlrenforcedabove, nlhdlrenforcemethods);

      /* if the nlhdlr enforces, then it must have added at least one enforcement method */
      assert(nlhdlrenforcemethods > enforcemethods || (nlhdlrenforcedbelow == enforcedbelow && nlhdlrenforcedabove == enforcedabove));

      /* remember nlhdlr and its data */
      nlhdlrssuccess[nsuccess] = nlhdlr;
      nlhdlrssuccessexprdata[nsuccess] = nlhdlrexprdata;
      ++nsuccess;

      /* update enforcement flags */
      enforcemethods = nlhdlrenforcemethods;
      enforcedbelow = nlhdlrenforcedbelow;
      enforcedabove = nlhdlrenforcedabove;

      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         /* call reverse propagation of nlhdlr
          * This can ensure that just created auxiliary variables take only values that are within the domain of functions that use them,
          * e.g., sqrt(x) in [-infty,infty] will ensure x >= 0, thus regardless of [-infty,infty] being pretty useless.
          * Another reason to do this already here is that LP solving and separation will be called next, which could already profit
          * from the tighter bounds (or: cons_expr_pow spits out a warning in separation if the child can be negative and exponent not integral).
          * NOTE: This assumes that reverseprop of the nlhdlr can be called before a preceding inteval call.
          */
         SCIP_CALL( SCIPreversepropConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, NULL, infeasible, &ntightenings, FALSE) );
      }
   }

   /* stop if the expression cannot be enforced but we are already in solving stage
    * (as long as the expression provides its callbacks, the default nlhdlr should have provided all enforcement methods)
    */
   if( (!enforcedbelow || !enforcedabove) && !*infeasible && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIPerrorMessage("no nonlinear handler provided enforcement for %s expression %s auxvar\n",
         SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)),
         (!enforcedbelow && !enforcedabove) ? "==" : (!enforcedbelow ? "<=" : ">="));
      return SCIP_ERROR;
   }

   if( nsuccess == 0 )
      return SCIP_OKAY;

   /* copy collected nlhdlrs into expr->enfos */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->enfos, nsuccess) );
   for( e = 0; e < nsuccess; ++e )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[e]) );  /*lint !e866 */
      expr->enfos[e]->nlhdlr = nlhdlrssuccess[e];
      expr->enfos[e]->nlhdlrexprdata = nlhdlrssuccessexprdata[e];
      expr->enfos[e]->issepainit = FALSE;
   }
   expr->nenfos = nsuccess;

   return SCIP_OKAY;
}

/** detect nlhdlrs that can handle the expressions */
static
SCIP_RETCODE detectNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected while creating the auxiliary vars */
   )
{
   SCIP_CONSEXPR_NLHDLR** nlhdlrssuccess;   /* buffer for nlhdlrs that had success detecting structure at expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrssuccessexprdata; /* buffer for exprdata of nlhdlrs */
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_INTERVAL activity;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(infeasible != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);  /* should only be called in presolve or initlp */

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );

   /* allocate some buffer for temporary storage of nlhdlr detect result */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlhdlrssuccess, conshdlrdata->nnlhdlrs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlhdlrssuccessexprdata, conshdlrdata->nnlhdlrs) );

   *infeasible = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL && conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         /* make sure activities in expression are uptodate
          * we do this here to have bounds for the auxiliary variables and for a reverseprop call at the end
          * we don't do auxiliary variables if in presolve, so do only in solving
          */
         SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, consdata->expr, &activity, TRUE) );
         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
         {
            SCIPdebugMsg(scip, "infeasibility detected in activity calculation of constraint <%s>\n", SCIPconsGetName(conss[i]));
            *infeasible = TRUE;
            break;
         }

#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIP_SOL* debugsol;

            SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

            if( debugsol != NULL ) /* it can be compiled WITH_DEBUG_SOLUTION, but still no solution given */
            {
               /* evaluate expression in debug solution, so we can set the solution value of created auxiliary variables
                * in SCIPcreateConsExprExprAuxVar()
                */
               SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, consdata->expr, debugsol, 0) );
            }
         }
#endif
      }

      /* compute integrality information for all subexpressions */
      SCIP_CALL( SCIPcomputeConsExprExprIntegral(scip, conshdlr, consdata->expr) );

      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         /* ensure auxiliary variable for root expression exists (not always necessary, but it simplifies things) */
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, consdata->expr, NULL) );
         assert(consdata->expr->auxvar != NULL);  /* couldn't this fail if the expression is only a variable? */

         /* change the bounds of the auxiliary variable of the root node to [lhs,rhs] */
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->expr->auxvar, consdata->lhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while creating vars: lhs of constraint (%g) > ub of node (%g)\n",
               consdata->lhs, SCIPvarGetUbLocal(consdata->expr->auxvar));
            break;
         }

         SCIP_CALL( SCIPtightenVarUb(scip, consdata->expr->auxvar, consdata->rhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while creating vars: rhs of constraint (%g) < lb of node (%g)\n",
               consdata->rhs, SCIPvarGetLbLocal(consdata->expr->auxvar));
            break;
         }
      }

      SCIP_CALL( SCIPexpriteratorInit(it, consdata->expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );  /* TODO init once for all conss */
      expr = SCIPexpriteratorGetCurrent(it);
      while( !SCIPexpriteratorIsEnd(it) )
      {
         if( expr->nenfos > 0 )
         {
            /* because of common sub-expressions it might happen that we already detected a nonlinear handler and added it to the expr
             * then also the subtree has been investigated already and we can stop iterating further down
             * HOWEVER: most likely we have been running DETECT with cons == NULL, which may interest less nlhdlrs
             * thus, if expr is the root expression, then rerun DETECT
             */
            if( expr == consdata->expr )
            {
               SCIP_CALL( freeEnfoData(scip, conshdlr, expr, FALSE) );
            }
            else
            {
               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }
         }

         /* during solve: if there is an auxiliary variable here, then there is some-one requiring that
          *   an auxvar equals (or approximates) to value of this expression or we are at the root expression (expr==consdata->expr)
          *   thus, we need to find nlhdlrs
          * during presolve: we do detect for all expressions for now, expecting the handler to only become
          *   active if they want to contribute in propagation
          */
         if( expr->auxvar != NULL || SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
         {
            SCIP_CALL( detectNlhdlr(scip, conshdlr, expr, expr == consdata->expr ? conss[i] : NULL, nlhdlrssuccess, nlhdlrssuccessexprdata, infeasible) );

            if( *infeasible )
               break;
         }

         expr = SCIPexpriteratorGetNext(it);
      }

      if( *infeasible )
      {
         SCIPdebugMsg(scip, "infeasibility detected while detecting nlhdlr\n");
         break;
      }
   }

   SCIPexpriteratorFree(&it);
   SCIPfreeBufferArray(scip, &nlhdlrssuccessexprdata);
   SCIPfreeBufferArray(scip, &nlhdlrssuccess);

   return SCIP_OKAY;
}

/** stores all variable expressions into a given constraint */
static
SCIP_RETCODE storeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already */
   if( consdata->varexprs != NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs == NULL);
   assert(consdata->nvarexprs == 0);

   /* create array to store all variable expressions; the number of variable expressions is bounded by SCIPgetNTotalVars() */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip)) );

   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, consdata->expr, consdata->varexprs, &(consdata->nvarexprs)) );
   assert(SCIPgetNTotalVars(scip) >= consdata->nvarexprs);

   /* realloc array if there are less variable expression than variables */
   if( SCIPgetNTotalVars(scip) > consdata->nvarexprs )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip), consdata->nvarexprs) );
   }

   return SCIP_OKAY;
}

/** frees all variable expression stored in storeVarExprs() */
static
SCIP_RETCODE freeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   int i;

   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* release variable expressions */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->varexprs[i]) );
      assert(consdata->varexprs[i] == NULL);
   }

   /* free variable expressions */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->varexprs, consdata->nvarexprs);
   consdata->varexprs = NULL;
   consdata->nvarexprs = 0;

   return SCIP_OKAY;
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPevalConsExprExpr(scip, SCIPconsGetHdlr(cons), consdata->expr, sol, soltag) );
   activity = SCIPgetConsExprExprValue(consdata->expr);

   /* consider constraint as violated if it is undefined in the current point */
   if( activity == SCIP_INVALID ) /*lint !e777*/
   {
      consdata->lhsviol = SCIPinfinity(scip);
      consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* compute violations */
   consdata->lhsviol = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : consdata->lhs  - activity;
   consdata->rhsviol = SCIPisInfinity(scip,  consdata->rhs) ? -SCIPinfinity(scip) : activity - consdata->rhs;

   return SCIP_OKAY;
}

/** returns absolute violation of a constraint
 *
 * @note This does not reevaluate the violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_Real getConsAbsViolation(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return MAX3(0.0, consdata->lhsviol, consdata->rhsviol);
}

/** computes relative violation of a constraint
 *
 * @note This does not reevaluate the absolute violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_RETCODE getConsRelViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            viol,               /**< buffer to store violation */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real scale;

   assert(cons != NULL);
   assert(viol != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *viol = getConsAbsViolation(cons);

   if( conshdlrdata->violscale == 'n' )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, *viol) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( conshdlrdata->violscale == 'a' )
   {
      scale = MAX(1.0, REALABS(SCIPgetConsExprExprValue(consdata->expr)));  /*lint !e666*/

      /* consider value of side that is violated for scaling, too */
      if( consdata->lhsviol > 0.0 && REALABS(consdata->lhs) > scale )
      {
         assert(!SCIPisInfinity(scip, -consdata->lhs));
         scale = REALABS(consdata->lhs);
      }
      else if( consdata->rhsviol > 0.0 && REALABS(consdata->rhs) > scale )
      {
         assert(!SCIPisInfinity(scip,  consdata->rhs));
         scale = REALABS(consdata->rhs);
      }

      *viol /= scale;
      return SCIP_OKAY;
   }

   /* if not 'n' or 'a', then it has to be 'g' at the moment */
   assert(conshdlrdata->violscale == 'g');
   if( soltag == 0 || consdata->gradnormsoltag != soltag )
   {
      /* we need the varexprs to conveniently access the gradient */
      SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

      /* update cached value of norm of gradient */
      consdata->gradnorm = 0.0;

      /* compute gradient */
      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, consdata->expr, sol, soltag) );

      /* gradient evaluation error -> no scaling */
      if( SCIPgetConsExprExprDerivative(consdata->expr) != SCIP_INVALID ) /*lint !e777*/
      {
         int i;
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real deriv;

            assert(consdata->expr->difftag == consdata->varexprs[i]->difftag);
            deriv = SCIPgetConsExprExprDerivative(consdata->varexprs[i]);
            if( deriv == SCIP_INVALID ) /*lint !e777*/
            {
               /* SCIPdebugMsg(scip, "gradient evaluation error for component %d\n", i); */
               consdata->gradnorm = 0.0;
               break;
            }

            consdata->gradnorm += deriv*deriv;
         }
      }
      consdata->gradnorm = sqrt(consdata->gradnorm);
      consdata->gradnormsoltag = soltag;
   }

   *viol /= MAX(1.0, consdata->gradnorm);

   return SCIP_OKAY;
}

/** returns whether constraint is currently violated
 *
 * @note This does not reevaluate the violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_Bool isConsViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   return getConsAbsViolation(cons) > SCIPfeastol(scip);
}

/** returns absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(x) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(x) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(x).
 * If f could not be evaluated, then return SCIPinfinity and set both violover and violunder to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that the expression has been evaluated
 */
static
SCIP_Real getExprAbsOrigViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   SCIP_Real auxvarvalue;

   assert(expr != NULL);
   assert(expr->auxvar != NULL);

   if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);

   if( SCIPgetConsExprExprNLocksNeg(expr) > 0 && auxvarvalue > expr->evalvalue )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - expr->evalvalue;
   }

   if( SCIPgetConsExprExprNLocksPos(expr) > 0 && expr->evalvalue > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return expr->evalvalue - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;
   return 0.0;
}

/** returns absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(w) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(w) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(w).
 * If f could not be evaluated, then return SCIPinfinity and set both violover and violunder to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that f(w) is passed in with auxvalue.
 */
static
SCIP_Real getExprAbsAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   SCIP_Real auxvarvalue;

   assert(expr != NULL);
   assert(expr->auxvar != NULL);

   if( auxvalue == SCIP_INVALID )  /*lint !e777*/
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);

   if( SCIPgetConsExprExprNLocksNeg(expr) > 0 && auxvarvalue > auxvalue )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - auxvalue;
   }

   if( SCIPgetConsExprExprNLocksPos(expr) > 0 && auxvalue > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return auxvalue - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;
   return 0.0;
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

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* check if we have catched variable events already */
   if( consdata->catchedevents )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "catchVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(consdata->varexprs[i]));

      SCIP_CALL( SCIPcatchConsExprExprVarEvent(scip, consdata->varexprs[i], eventhdlr, cons) );
   }

   consdata->catchedevents = TRUE;

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

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we have catched variable events already */
   if( !consdata->catchedevents )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   SCIPdebugMsg(scip, "dropVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = consdata->nvarexprs - 1; i >= 0; --i )
   {
      assert(consdata->varexprs[i] != NULL);

      SCIP_CALL( SCIPdropConsExprExprVarEvent(scip, consdata->varexprs[i], eventhdlr, cons) );
   }

   consdata->catchedevents = FALSE;

   return SCIP_OKAY;
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSEXPR_EXPR* expr;

   eventtype = SCIPeventGetType(event);
   assert((eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0 || (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0);

   assert(eventdata != NULL);
   expr = (SCIP_CONSEXPR_EXPR*) eventdata;

   SCIPdebugMsg(scip, "  exec event %#x for variable <%s>\n", eventtype, SCIPvarGetName(SCIPeventGetVar(event)));

   /* for real variables notify constraints to repropagate and possibly resimplify */
   if( SCIPisConsExprExprVar(expr) )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS** conss;
      int nconss;
      int c;

      nconss = SCIPgetConsExprExprVarNConss(expr);
      conss = SCIPgetConsExprExprVarConss(expr);
      assert(conss != NULL || nconss == 0);

      for( c = 0; c < nconss; ++c )
      {
         assert(conss[c] != NULL);  /*lint !e613*/
         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/

         /* if boundchange, then mark constraints to be propagated again */
         if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != (unsigned int) 0 )
         {
            consdata->ispropagated = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for propagate and simplify\n", SCIPconsGetName(conss[c]));  /*lint !e613*/

            /* store handler for below */
            conshdlr = SCIPconsGetHdlr(conss[c]);  /*lint !e613*/
         }

         /* if still in presolve, then mark constraints to be simplified again */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
         {
            consdata->issimplified = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for simplify\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
         }
      }
   }

   /* update curboundstag and lastboundrelax */
   if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != (unsigned int) 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      if( conshdlr == NULL )
         conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* increase tag on bounds */
      /* TODO maybe do not increase if we did not use the new tag yet, e.g., when there is a sequence of bound changes,
       * so we do not run out of numbers so soon
       */
      ++conshdlrdata->curboundstag;
      assert(conshdlrdata->curboundstag > 0);

      /* remember also if we relaxed bounds now */
      if( eventtype & SCIP_EVENTTYPE_BOUNDRELAXED )
         conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
   }

   return SCIP_OKAY;
}

/** propagates variable locks through expression and adds lock to variables */
static
SCIP_RETCODE propagateLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   nlockspos,          /**< number of positive locks */
   int                   nlocksneg           /**< number of negative locks */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPRITERATOR_USERDATA ituserdata;

   assert(expr != NULL);

   /* if no locks, then nothing to do, then do nothing */
   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR | SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR);
   assert(SCIPexpriteratorGetCurrent(it) == expr); /* iterator should not have moved */

   /* store locks in root node */
   ituserdata.intvals[0] = nlockspos;
   ituserdata.intvals[1] = nlocksneg;
   SCIPexpriteratorSetCurrentUserData(it, ituserdata);

   while( !SCIPexpriteratorIsEnd(it) )
   {
      /* collect locks */
      ituserdata = SCIPexpriteratorGetCurrentUserData(it);
      nlockspos = ituserdata.intvals[0];
      nlocksneg = ituserdata.intvals[1];

      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_ENTEREXPR:
         {
            if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrVar(conshdlr) )
            {
               /* if a variable, then also add nlocksneg/nlockspos via SCIPaddVarLocks() */
               SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetConsExprExprVarVar(expr), nlocksneg, nlockspos) );
            }

            /* add locks to expression */
            expr->nlockspos += nlockspos;
            expr->nlocksneg += nlocksneg;

            /* add monotonicity information if expression has been locked for the first time */
            if( expr->nlockspos == nlockspos && expr->nlocksneg == nlocksneg && expr->nchildren > 0
               && expr->exprhdlr->monotonicity != NULL )
            {
               int i;

               assert(expr->monotonicity == NULL);
               assert(expr->monotonicitysize == 0);

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->monotonicity, expr->nchildren) );
               expr->monotonicitysize = expr->nchildren;

               /* store the monotonicity for each child */
               for( i = 0; i < expr->nchildren; ++i )
               {
                  SCIP_CALL( (*expr->exprhdlr->monotonicity)(scip, expr, i, &expr->monotonicity[i]) );
               }
            }
            break;
         }

         case SCIP_CONSEXPRITERATOR_LEAVEEXPR :
         {
            /* remove monotonicity information if expression has been unlocked */
            if( expr->nlockspos == 0 && expr->nlocksneg == 0 && expr->monotonicity != NULL )
            {
               assert(expr->monotonicitysize > 0);
               /* keep this assert for checking whether someone changed an expression without updating locks properly */
               assert(expr->monotonicitysize == expr->nchildren);

               SCIPfreeBlockMemoryArray(scip, &expr->monotonicity, expr->monotonicitysize);
               expr->monotonicitysize = 0;
            }
            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
         {
            SCIP_MONOTONE monotonicity;

            assert(expr->monotonicity != NULL || expr->exprhdlr->monotonicity == NULL);

            /* get monotonicity of child */
            /* NOTE: the monotonicity stored in an expression might be different from the result obtained by
             * SCIPgetConsExprExprMonotonicity
             */
            monotonicity = expr->monotonicity != NULL ? expr->monotonicity[SCIPexpriteratorGetChildIdxDFS(it)] : SCIP_MONOTONE_UNKNOWN;

            /* compute resulting locks of the child expression */
            switch( monotonicity )
            {
               case SCIP_MONOTONE_INC:
                  ituserdata.intvals[0] = nlockspos;
                  ituserdata.intvals[1] = nlocksneg;
                  break;
               case SCIP_MONOTONE_DEC:
                  ituserdata.intvals[0] = nlocksneg;
                  ituserdata.intvals[1] = nlockspos;
                  break;
               case SCIP_MONOTONE_UNKNOWN:
                  ituserdata.intvals[0] = nlockspos + nlocksneg;
                  ituserdata.intvals[1] = nlockspos + nlocksneg;
                  break;
               case SCIP_MONOTONE_CONST:
                  ituserdata.intvals[0] = 0;
                  ituserdata.intvals[1] = 0;
                  break;
            }
            /* set locks in child expression */
            SCIPexpriteratorSetChildUserData(it, ituserdata);

            break;
         }

         default :
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** main function for adding locks to expressions and variables; locks for an expression constraint are used to update
 *  locks for all sub-expressions and variables; locks of expressions depend on the monotonicity of expressions
 *  w.r.t. their children, e.g., consider the constraint x^2 <= 1 with x in [-2,-1] implies an up-lock for the root
 *  expression (pow) and a down-lock for its child x because x^2 is decreasing on [-2,-1]; since the monotonicity (and thus
 *  the locks) might also depend on variable bounds, the function remembers the computed monotonicity information ofcan
 *  each expression until all locks of an expression have been removed, which implies that updating the monotonicity
 *  information during the next locking of this expression does not break existing locks
 *
 *  @note when modifying the structure of an expression, e.g., during simplification, it is necessary to remove all
 *        locks from an expression and repropagating them after the structural changes have been applied; because of
 *        existing common sub-expressions, it might be necessary to remove the locks of all constraints to ensure
 *        that an expression is unlocked (see canonicalizeConstraints() for an example)
 */
static
SCIP_RETCODE addLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< expression constraint */
   int                   nlockspos,          /**< number of positive rounding locks */
   int                   nlocksneg           /**< number of negative rounding locks */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* no constraint sides -> nothing to lock */
   if( SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, -consdata->lhs) )
      return SCIP_OKAY;

   /* make sure activities are uptodate when root expression is locked for the first time */
   if( consdata->expr->nlockspos == 0 && consdata->expr->nlocksneg == 0 )
   {
      SCIP_INTERVAL activity;
      SCIP_CALL( SCIPevalConsExprExprActivity(scip, SCIPconsGetHdlr(cons), consdata->expr, &activity, TRUE) );
   }

   /* remember locks */
   consdata->nlockspos += nlockspos;
   consdata->nlocksneg += nlocksneg;

   assert(consdata->nlockspos >= 0);
   assert(consdata->nlocksneg >= 0);

   /* compute locks for lock propagation */
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos + nlocksneg, nlockspos + nlocksneg));
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos, nlocksneg));
   }
   else
   {
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlocksneg, nlockspos));
   }

   return SCIP_OKAY;
}

/** returns an equivalent expression for a given expression if possible; it adds the expression to key2expr if the map
 *  does not contain the key
 */
static
SCIP_RETCODE findEqualExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR *  expr,               /**< expression to replace */
   SCIP_MULTIHASH*       key2expr,           /**< mapping of hashes to expressions */
   SCIP_CONSEXPR_EXPR**  newexpr             /**< pointer to store an equivalent expression (NULL if there is none) */
   )
{  /*lint --e{438}*/
   SCIP_MULTIHASHLIST* multihashlist;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(key2expr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;
   multihashlist = NULL;

   do
   {
      /* search for an equivalent expression */
      *newexpr = (SCIP_CONSEXPR_EXPR*)(SCIPmultihashRetrieveNext(key2expr, &multihashlist, (void*)expr));

      if( *newexpr == NULL )
      {
         /* processed all expressions like expr from hash table, so insert expr */
         SCIP_CALL( SCIPmultihashInsert(key2expr, (void*) expr) );
         break;
      }
      else if( expr != *newexpr )
      {
         assert(SCIPcompareConsExprExprs(expr, *newexpr) == 0);
         break;
      }
      else
      {
         /* can not replace expr since it is already contained in the hashtablelist */
         assert(expr == *newexpr);
         *newexpr = NULL;
         break;
      }
   }
   while( TRUE ); /*lint !e506*/

   return SCIP_OKAY;
}

/** hashes an expression using an already existing iterator
 *
 * The iterator must by of type DFS with allowrevisit=FALSE and the only leaveexpr stage enabled.
 * The hashes of all visited expressions will be stored in the iterators expression data.
 */
static
SCIP_RETCODE hashExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to hash */
   SCIP_CONSEXPR_ITERATOR* hashiterator,     /**< iterator to use for hashing */
   int*                  nvisitedexprs       /**< counter to increment by the number of expressions visited, or NULL */
   )
{
   SCIP_CONSEXPRITERATOR_USERDATA iterdata;
   unsigned int* childrenhashes;
   int childrenhashessize;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashiterator != NULL);

   childrenhashessize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &childrenhashes, childrenhashessize) );

   for( expr = SCIPexpriteratorRestartDFS(hashiterator, expr); !SCIPexpriteratorIsEnd(hashiterator); expr = SCIPexpriteratorGetNext(hashiterator) ) /*lint !e441*/
   {
      assert(SCIPexpriteratorGetStageDFS(hashiterator) == SCIP_CONSEXPRITERATOR_LEAVEEXPR);

      if( nvisitedexprs != NULL )
         ++*nvisitedexprs;

      /* collect hashes of children */
      if( childrenhashessize < expr->nchildren )
      {
         childrenhashessize = SCIPcalcMemGrowSize(scip, expr->nchildren);
         SCIP_CALL( SCIPreallocBufferArray(scip, &childrenhashes, childrenhashessize) );
      }
      for( i = 0; i < expr->nchildren; ++i )
         childrenhashes[i] = SCIPexpriteratorGetExprUserData(hashiterator, expr->children[i]).uintval;

      SCIP_CALL( SCIPhashConsExprExprHdlr(scip, expr, &iterdata.uintval, childrenhashes) );

      SCIPexpriteratorSetCurrentUserData(hashiterator, iterdata);
   }

   SCIPfreeBufferArray(scip, &childrenhashes);

   return SCIP_OKAY;
}

/** get key of hash element */
static
SCIP_DECL_HASHGETKEY(hashCommonSubexprGetKey)
{
   return elem;
}  /*lint !e715*/

/** checks if two expressions are structurally the same */
static
SCIP_DECL_HASHKEYEQ(hashCommonSubexprEq)
{
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;

   expr1 = (SCIP_CONSEXPR_EXPR*)key1;
   expr2 = (SCIP_CONSEXPR_EXPR*)key2;
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   return expr1 == expr2 || SCIPcompareConsExprExprs(expr1, expr2) == 0;
}  /*lint !e715*/

/** get value of hash element when comparing with another expression */
static
SCIP_DECL_HASHKEYVAL(hashCommonSubexprKeyval)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_ITERATOR* hashiterator;

   expr = (SCIP_CONSEXPR_EXPR*) key;
   assert(expr != NULL);

   hashiterator = (SCIP_CONSEXPR_ITERATOR*) userptr;
   assert(hashiterator != NULL);

   return SCIPexpriteratorGetExprUserData(hashiterator, expr).uintval;
}  /*lint !e715*/

/** replaces common sub-expressions in the current expression graph by using a hash key for each expression; the
 *  algorithm consists of two steps:
 *
 *  1. traverse through all expressions trees of given constraints and compute for each of them a (not necessarily
 *     unique) hash
 *
 *  2. initialize an empty hash table and traverse through all expression; check for each of them if we can find a
 *     structural equivalent expression in the hash table; if yes we replace the expression by the expression inside the
 *     hash table, otherwise we add it to the hash table
 *
 *  @note the hash keys of the expressions are used for the hashing inside the hash table; to compute if two expressions
 *  (with the same hash) are structurally the same we use the function SCIPcompareConsExprExprs()
 */
static
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_CONSEXPR_ITERATOR* hashiterator;
   SCIP_CONSEXPR_ITERATOR* repliterator;
   SCIP_MULTIHASH* key2expr;
   SCIP_CONSDATA* consdata;
   int i;
   int nexprs = 0;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);

   if( nconss == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexpriteratorCreate(&hashiterator, SCIPconsGetHdlr(conss[0]), SCIPblkmem(scip)) );

   /* compute all hashes for each sub-expression */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr == NULL )
         continue;

      if( !SCIPexpriteratorIsInit(hashiterator) )
      {
         /* first constraint with non-NULL expr: initialize iterator (set type and stopstage) */
         SCIP_CALL( SCIPexpriteratorInit(hashiterator, consdata->expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
         SCIPexpriteratorSetStagesDFS(hashiterator, SCIP_CONSEXPRITERATOR_LEAVEEXPR);
      }

      SCIP_CALL( hashExpr(scip, consdata->expr, hashiterator, &nexprs) );
   }

   /* replace equivalent sub-expressions */
   SCIP_CALL( SCIPmultihashCreate(&key2expr, SCIPblkmem(scip), nexprs,
         hashCommonSubexprGetKey, hashCommonSubexprEq, hashCommonSubexprKeyval, (void*)hashiterator) );

   SCIP_CALL( SCIPexpriteratorCreate(&repliterator, SCIPconsGetHdlr(conss[0]), SCIPblkmem(scip)) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSEXPR_EXPR* newroot;
      SCIP_CONSEXPR_EXPR* newchild;
      SCIP_CONSEXPR_EXPR* child;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr == NULL )
         continue;

      /* check the root for equivalence separately first */
      SCIP_CALL( findEqualExpr(scip, consdata->expr, key2expr, &newroot) );

      if( newroot != NULL )
      {
         assert(newroot != consdata->expr);
         assert(SCIPcompareConsExprExprs(consdata->expr, newroot) == 0);

         SCIPdebugMsg(scip, "replacing common root expression of constraint <%s>: %p -> %p\n", SCIPconsGetName(conss[i]), (void*)consdata->expr, (void*)newroot);

         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

         consdata->expr = newroot;
         SCIPcaptureConsExprExpr(newroot);

         continue;
      }

      /* replace equivalent sub-expressions in the tree */
      SCIP_CALL( SCIPexpriteratorInit(repliterator, consdata->expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
      SCIPexpriteratorSetStagesDFS(repliterator, SCIP_CONSEXPRITERATOR_VISITINGCHILD);

      while( !SCIPexpriteratorIsEnd(repliterator) )
      {
         child = SCIPexpriteratorGetChildExprDFS(repliterator);
         assert(child != NULL);

         /* try to find an equivalent expression */
         SCIP_CALL( findEqualExpr(scip, child, key2expr, &newchild) );

         /* replace child with newchild */
         if( newchild != NULL )
         {
            assert(child != newchild);
            assert(SCIPcompareConsExprExprs(child, newchild) == 0);

            SCIPdebugMsg(scip, "replacing common child expression %p -> %p\n", (void*)child, (void*)newchild);

            SCIP_CALL( SCIPreplaceConsExprExprChild(scip, SCIPexpriteratorGetCurrent(repliterator), SCIPexpriteratorGetChildIdxDFS(repliterator), newchild) );

            (void) SCIPexpriteratorSkipDFS(repliterator);
         }
         else
         {
            (void) SCIPexpriteratorGetNext(repliterator);
         }
      }
   }

   /* free memory */
   SCIPexpriteratorFree(&repliterator);
   SCIPmultihashFree(&key2expr);
   SCIPexpriteratorFree(&hashiterator);

   return SCIP_OKAY;
}

/** helper function to either simplify or reformulate an expression and its subexpressions */
static
SCIP_RETCODE reformulateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_Bool               simplify,         /**< should the expression be simplified or reformulated? */
   SCIP_CONSEXPR_EXPR**    simplified,       /**< buffer to store simplified expression */
   SCIP_Bool*              changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*              infeasible        /**< buffer to store whether infeasibility has been detected */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_ITERATOR* it;

   assert(scip != NULL);
   assert(rootexpr != NULL);
   assert(simplified != NULL);
   assert(changed != NULL);
   assert(infeasible != NULL);

   /* simplify bottom up
    * when leaving an expression it simplifies it and stores the simplified expr in its iterators expression data
    * after the child was visited, it is replaced with the simplified expr
    */
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, rootexpr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );  /* TODO can we set allowrevisited to FALSE?*/
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITEDCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   *changed = FALSE;
   *infeasible = FALSE;
   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_VISITEDCHILD:
         {
            SCIP_CONSEXPR_EXPR* newchild;
            SCIP_CONSEXPR_EXPR* child;

            newchild = (SCIP_CONSEXPR_EXPR*)SCIPexpriteratorGetChildUserDataDFS(it).ptrval;
            child = SCIPexpriteratorGetChildExprDFS(it);
            assert(newchild != NULL);

            /* if child got simplified, replace it with the new child */
            if( newchild != child )
            {
               SCIP_CALL( SCIPreplaceConsExprExprChild(scip, expr, SCIPexpriteratorGetChildIdxDFS(it), newchild) );
            }

            /* we do not need to hold newchild anymore */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &newchild) );

            break;
         }

         case SCIP_CONSEXPRITERATOR_LEAVEEXPR:
         {
            SCIP_CONSEXPR_EXPR* refexpr = NULL;
            SCIP_CONSEXPRITERATOR_USERDATA iterdata;

            /* use simplification of expression handlers */
            if( simplify )
            {
               if( SCIPhasConsExprExprHdlrSimplify(expr->exprhdlr) )
               {
                  SCIP_CALL( SCIPsimplifyConsExprExprHdlr(scip, conshdlr, expr, &refexpr) );
                  if( expr != refexpr )
                  {
                     SCIP_INTERVAL activity;

                     *changed = TRUE;

                     /* make sure valid activities are available for the new expr (and its children)
                      * we might expect them to be present in nlhdlr detect later
                      */
                     SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, refexpr, &activity, TRUE) );

                     if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
                        *infeasible = TRUE;
                  }
               }
               else
               {
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sum")  != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "var") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "abs") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "log") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "exp") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "pow") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sin") != 0);
                  assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "cos") != 0);

                  /* if an expression handler doesn't implement simplify, we assume all those type of expressions are simplified
                   * we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created
                   */
                  refexpr = expr;
                  SCIPcaptureConsExprExpr(refexpr);
               }
               assert(refexpr != NULL);
            }
            else /* use nonlinear handler to reformulate the expression */
            {
               SCIP_CONSHDLRDATA* conshdlrdata;
               int k;

               conshdlrdata = SCIPconshdlrGetData(conshdlr);
               assert(conshdlrdata != NULL);

               /* iterate through nonlinear handlers and call reformulation callbacks;
                *
                * TODO store nonlinear handlers that implement the reformulation callback separately
                * TODO sort nonlinear handlers according to their priorities
                */
               for( k = 0; k < conshdlrdata->nnlhdlrs; ++k )
               {
                  assert(conshdlrdata->nlhdlrs[k] != NULL);

                  if( SCIPhasConsExprNlhdlrReformulate(conshdlrdata->nlhdlrs[k]) )
                  {
                     SCIP_CALL( SCIPreformulateConsExprNlhdlr(scip, conshdlr, conshdlrdata->nlhdlrs[k], expr, &refexpr) );

                     /* stop calling other nonlinear handlers as soon as the reformulation was successful */
                     if( refexpr != NULL && refexpr != expr )
                     {
                        SCIP_INTERVAL activity;

                        *changed = TRUE;

                        /* make sure valid activities are available for the new expr (and its children)
                         * we might expect them to be present in nlhdlr detect later
                         */
                        SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, refexpr, &activity, TRUE) );

                        if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
                           *infeasible = TRUE;

                        break;
                     }
                  }
               }

               /* no nonlinear handlers implements the reformulation callback -> capture expression manually */
               if( refexpr == NULL )
               {
                  refexpr = expr;
                  SCIPcaptureConsExprExpr(refexpr);
               }
            }

            iterdata.ptrval = (void*) refexpr;
            SCIPexpriteratorSetCurrentUserData(it, iterdata);

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }
   }

   *simplified = (SCIP_CONSEXPR_EXPR*)SCIPexpriteratorGetExprUserData(it, rootexpr).ptrval;
   assert(*simplified != NULL);

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}


/** scales the sides of the constraint l <= sum_i c_i f_i(x) <= r according to the following rules:
 *
 *  let n_+ the number of positive coefficients c_i and n_- be the number of negative coefficients
 *
 *   i. scale by -1 if n_+ < n_-
 *
 *  ii. scale by -1 if n_+ = n_- & r = INF
 */
static
SCIP_RETCODE scaleConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_Bool*            changed             /**< buffer to store if the expression of cons changed */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPgetConsExprExprHdlr(consdata->expr) == SCIPgetConsExprExprHdlrSum(conshdlr) )
   {
      SCIP_Real* coefs;
      SCIP_Real constant;
      int nchildren;
      int counter = 0;

      coefs = SCIPgetConsExprExprSumCoefs(consdata->expr);
      constant = SCIPgetConsExprExprSumConstant(consdata->expr);
      nchildren = SCIPgetConsExprExprNChildren(consdata->expr);

      /* handle special case when constraint is l <= -f(x) <= r and f(x) not a sum: simplfy ensures f is not a sum */
      if( nchildren == 1 && constant == 0.0 && coefs[0] == -1.0 )
      {
         SCIP_CONSEXPR_EXPR* expr;
         expr = consdata->expr;

         consdata->expr = SCIPgetConsExprExprChildren(expr)[0];
         assert(SCIPgetConsExprExprHdlr(consdata->expr) != SCIPgetConsExprExprHdlrSum(conshdlr));

         SCIPcaptureConsExprExpr(consdata->expr);

         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
         *changed = TRUE;
         return SCIP_OKAY;
      }

      /* compute n_+ - n_i */
      for( i = 0; i < nchildren; ++i )
         counter += coefs[i] > 0 ? 1 : -1;

      if( counter < 0 || (counter == 0 && SCIPisInfinity(scip, consdata->rhs)) )
      {
         SCIP_CONSEXPR_EXPR* expr;
         SCIP_Real* newcoefs;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &newcoefs, nchildren) );

         for( i = 0; i < nchildren; ++i )
            newcoefs[i] = -coefs[i];

         /* create a new sum expression */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, nchildren, SCIPgetConsExprExprChildren(consdata->expr), newcoefs, -constant) );

         /* replace expression in constraint data and scale sides */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );
         consdata->expr = expr;
         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         /* free memory */
         SCIPfreeBufferArray(scip, &newcoefs);

         *changed = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** helper method to decide whether a given expression is product of at least two binary variables */
static
SCIP_Bool isBinaryProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   int nchildren;
   int i;

   assert(expr != NULL);

   /* check whether the expression is a product */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      return FALSE;

   nchildren = SCIPgetConsExprExprNChildren(expr);

   /* don't consider products with a coefficient != 1 and products with a single child; simplification will take care
    * of this expression later
    */
   if( nchildren <= 1 || SCIPgetConsExprExprProductCoef(expr) != 1.0 )
      return FALSE;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CONSEXPR_EXPR* child;
      SCIP_VAR* var;
      SCIP_Real ub;
      SCIP_Real lb;

      child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      if( !SCIPisConsExprExprVar(child) )
         return FALSE;

      var = SCIPgetConsExprExprVarVar(child);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      /* check whether variable is integer and has [0,1] as variable bounds */
      if( !SCIPvarIsIntegral(var) || !SCIPisEQ(scip, lb, 0.0) || !SCIPisEQ(scip, ub, 1.0) )
         return FALSE;
   }

   return TRUE;
}

/** helper method to collect all bilinear binary product terms */
static
SCIP_RETCODE getBilinearBinaryTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   sumexpr,            /**< sum expression */
   SCIP_VAR**            xs,                 /**< array to collect first variable of each bilinear binary product */
   SCIP_VAR**            ys,                 /**< array to collect second variable of each bilinear binary product */
   int*                  childidxs,          /**< array to store the index of the child of each stored bilinear binary product */
   int*                  nterms              /**< pointer to store the total number of bilinear binary terms */
   )
{
   int i;

   assert(sumexpr != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(childidxs != NULL);
   assert(nterms != NULL);

   *nterms = 0;

   for( i = 0; i < SCIPgetConsExprExprNChildren(sumexpr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child;

      child = SCIPgetConsExprExprChildren(sumexpr)[i];
      assert(child != NULL);

      if( SCIPgetConsExprExprNChildren(child) == 2 && isBinaryProduct(scip, conshdlr, child) )
      {
         SCIP_VAR* x = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]);
         SCIP_VAR* y = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[1]);

         assert(x != NULL);
         assert(y != NULL);

         if( x != y )
         {
            xs[*nterms] = x;
            ys[*nterms] = y;
            childidxs[*nterms] = i;
            ++(*nterms);
         }
      }
   }

   return SCIP_OKAY;
}

/** helper method to reformulate x_i * sum_j c_ij x_j */
static
SCIP_RETCODE reformulateFactorizedBinaryQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR*             facvar,             /**< variable that has been factorized */
   SCIP_VAR**            vars,               /**< variables of sum_j c_ij x_j */
   SCIP_Real*            coefs,              /**< coefficients of sum_j c_ij x_j */
   int                   nvars,              /**< total number of variables in sum_j c_ij x_j */
   SCIP_CONSEXPR_EXPR**  newexpr,            /**< pointer to store the new expression */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_VAR* auxvar;
   SCIP_CONS* newcons;
   SCIP_Real minact = 0.0;
   SCIP_Real maxact = 0.0;
   SCIP_Bool integral = TRUE;
   char name [SCIP_MAXSTRLEN];
   int i;

   assert(facvar != NULL);
   assert(vars != NULL);
   assert(nvars > 1);
   assert(newexpr != NULL);

   /* compute minimum and maximum activity of sum_j c_ij x_j */
   /* TODO could compute minact and maxact for facvar=0 and facvar=1 separately, taking implied bounds into account, allowing for possibly tighter big-M's below */
   for( i = 0; i < nvars; ++i )
   {
      minact += MIN(coefs[i], 0.0);
      maxact += MAX(coefs[i], 0.0);
      integral = integral && SCIPisIntegral(scip, coefs[i]);
   }
   assert(minact <= maxact);

   /* create and add auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateVarBasic(scip, &auxvar, name, minact, maxact, 0.0, integral ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );

   /* create and add z - maxact x <= 0 */
   if( !SCIPisZero(scip, maxact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -maxact, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add  0 <= z - minact x */
   if( !SCIPisZero(scip, minact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -minact, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add minact <= sum_j c_j x_j - z + minact x_i */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, minact, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, minact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, minact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create and add sum_j c_j x_j - z + maxact x_i <= maxact */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_4", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, -SCIPinfinity(scip), maxact) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, maxact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, maxact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create variable expression */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, newexpr, auxvar) );

   /* release auxvar */
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   return SCIP_OKAY;
}

/** helper method to generate an expression for a sum of product of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getFactorizedBinaryQuadraticExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_CONSEXPR_EXPR*   sumexpr,            /**< sum expression */
   int                   minterms,           /**< minimum number of terms in a the sum of x_i sum_j c_j x_j */
   SCIP_CONSEXPR_EXPR**  newexpr,            /**< pointer to store the expression that represents the binary quadratic */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_CONSEXPR_EXPR** exprs = NULL;
   SCIP_VAR** tmpvars = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_VAR** xs = NULL;
   SCIP_VAR** ys = NULL;
   SCIP_Real* exprcoefs = NULL;
   SCIP_Real* tmpcoefs = NULL;
   SCIP_Real* sumcoefs;
   SCIP_Bool* isused  = NULL;
   int* childidxs = NULL;
   int* count = NULL;
   int nchildren;
   int nexprs = 0;
   int nterms;
   int nvars;
   int ntotalvars;
   int i;

   assert(sumexpr != NULL);
   assert(minterms > 1);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* check whether sumexpr is indeed a sum */
   if( SCIPgetConsExprExprHdlr(sumexpr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   nchildren = SCIPgetConsExprExprNChildren(sumexpr);
   sumcoefs = SCIPgetConsExprExprSumCoefs(sumexpr);
   nvars = SCIPgetNVars(scip);
   ntotalvars = SCIPgetNTotalVars(scip);

   /* check whether there are enough terms available */
   if( nchildren < minterms )
      return SCIP_OKAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &childidxs, nchildren) );

   /* collect all bilinear binary product terms */
   SCIP_CALL( getBilinearBinaryTerms(scip, conshdlr, sumexpr, xs, ys, childidxs, &nterms) );

   /* check whether there are enough terms available */
   if( nterms < minterms )
      goto TERMINATE;

   /* store how often each variable appears in a bilinear binary product */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &count, ntotalvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &isused, nchildren) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprcoefs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, MIN(nterms, nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, MIN(nterms, nvars)) );

   for( i = 0; i < nterms; ++i )
   {
      int xidx;
      int yidx;

      assert(xs[i] != NULL);
      assert(ys[i] != NULL);

      xidx = SCIPvarGetIndex(xs[i]);
      assert(xidx < ntotalvars);
      yidx = SCIPvarGetIndex(ys[i]);
      assert(yidx < ntotalvars);

      ++count[xidx];
      ++count[yidx];

      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(xs[i]), count[xidx]);
      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(ys[i]), count[yidx]);
   }

   /* sort variables; don't change order of count array because it depends on problem indices */
   {
      int* tmpcount;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpcount, count, nvars) );
      SCIPsortDownIntPtr(tmpcount, (void**)vars, nvars);
      SCIPfreeBufferArray(scip, &tmpcount);
   }

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* facvar = vars[i];
      int ntmpvars = 0;
      int j;

      /* skip candidate if there are not enough terms left */
      if( count[SCIPvarGetIndex(vars[i])] < minterms )
         continue;

      SCIPdebugMsg(scip, "consider facvar = %s with count = %d\n", SCIPvarGetName(facvar), count[SCIPvarGetIndex(vars[i])]);

      /* collect variables for x_i * sum_j c_ij x_j */
      for( j = 0; j < nterms; ++j )
      {
         int childidx = childidxs[j];
         assert(childidx >= 0 && childidx < nchildren);

         if( !isused[childidx] && (xs[j] == facvar || ys[j] == facvar) )
         {
            SCIP_Real coef;
            int xidx;
            int yidx;

            coef = sumcoefs[childidx];
            assert(coef != 0.0);

            /* collect corresponding variable */
            tmpvars[ntmpvars] = (xs[j] == facvar) ? ys[j] : xs[j];
            tmpcoefs[ntmpvars] = coef;
            ++ntmpvars;

            /* update counters */
            xidx = SCIPvarGetIndex(xs[j]);
            assert(xidx < ntotalvars);
            yidx = SCIPvarGetIndex(ys[j]);
            assert(yidx < ntotalvars);
            --count[xidx];
            --count[yidx];
            assert(count[xidx] >= 0);
            assert(count[yidx] >= 0);

            /* mark term to be used */
            isused[childidx] = TRUE;
         }
      }
      assert(ntmpvars >= minterms);
      assert(SCIPvarGetIndex(facvar) < ntotalvars);
      assert(count[SCIPvarGetIndex(facvar)] == 0); /* facvar should not appear in any other bilinear term */

      /* create required constraints and store the generated expression */
      SCIP_CALL( reformulateFactorizedBinaryQuadratic(scip, conshdlr, cons, facvar, tmpvars, tmpcoefs, ntmpvars, &exprs[nexprs], naddconss) );
      exprcoefs[nexprs] = 1.0;
      ++nexprs;
   }

   /* factorization was only successful if at least one expression has been generated */
   if( nexprs > 0 )
   {
      int nexprsold = nexprs;

      /* add all children of the sum that have not been used */
      for( i = 0; i < nchildren; ++i )
      {
         if( !isused[i] )
         {
            exprs[nexprs] = SCIPgetConsExprExprChildren(sumexpr)[i];
            exprcoefs[nexprs] = sumcoefs[i];
            ++nexprs;
         }
      }

      /* create a new sum expression */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, newexpr, nexprs, exprs, exprcoefs, SCIPgetConsExprExprSumConstant(sumexpr)) );

      /* release all expressions that have been generated by reformulateFactorizedBinaryQuadratic() */
      for( i = 0; i < nexprsold; ++i )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[i]) );
      }
   }

TERMINATE:
   /* free memory */
   SCIPfreeBufferArrayNull(scip, &tmpcoefs);
   SCIPfreeBufferArrayNull(scip, &tmpvars);
   SCIPfreeBufferArrayNull(scip, &exprcoefs);
   SCIPfreeBufferArrayNull(scip, &exprs);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &isused);
   SCIPfreeBufferArrayNull(scip, &count);
   SCIPfreeBufferArray(scip, &childidxs);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);

   return SCIP_OKAY;
}

/** helper method to create an AND constraint or varbound constraints for a given binary product expression */
static
SCIP_RETCODE getBinaryProductExprDo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   prodexpr,           /**< product expression */
   SCIP_CONSEXPR_EXPR**  newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   SCIP_Bool             empathy4and         /**< whether to use an AND constraint, if possible */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS* cons;
   SCIP_Real* coefs;
   SCIP_VAR* w;
   char name[SCIP_MAXSTRLEN];
   int nchildren;
   int i;

   assert(conshdlr != NULL);
   assert(prodexpr != NULL);
   assert(newexpr != NULL);

   nchildren = SCIPgetConsExprExprNChildren(prodexpr);
   assert(nchildren >= 2);

   /* memory to store the variables of the variable expressions (+1 for w) */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nchildren + 1) );

   /* prepare the names of the variable and the constraints */
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform");
   for( i = 0; i < nchildren; ++i )
   {
      vars[i] = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(prodexpr)[i]);
      coefs[i] = 1.0;
      assert(vars[i] != NULL);
      (void) strcat(name, "_");
      (void) strcat(name, SCIPvarGetName(vars[i]));
   }

   /* create and add variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIPdebugMsg(scip, "  created auxiliary variable %s\n", name);

   /* use variable bound constraints if it is a bilinear product and there is no empathy for an AND constraint */
   if( nchildren == 2 && !empathy4and )
   {
      SCIP_VAR* x = vars[0];
      SCIP_VAR* y = vars[1];

      assert(x != NULL);
      assert(y != NULL);
      assert(x != y);

      /* create and add x - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, x, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add y - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, y, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add x + y - w <= 1 */
      vars[2] = w;
      coefs[2] = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 3, vars, coefs, -SCIPinfinity(scip), 1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 3;
   }
   else
   {
      /* create, add, and release AND constraint */
      SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, name, w, nchildren, vars) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPdebugMsg(scip, "  create AND constraint\n");

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 1;
   }

   /* create variable expression */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, newexpr, w) );

   /* release created variable */
   SCIP_CALL( SCIPreleaseVar(scip, &w) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** helper method to generate an expression for the product of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getBinaryProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_CONSEXPR_EXPR*   prodexpr,           /**< product expression */
   SCIP_CONSEXPR_EXPR**  newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nchildren;

   assert(prodexpr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* only consider products of binary variables */
   if( !isBinaryProduct(scip, conshdlr, prodexpr) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   nchildren = SCIPgetConsExprExprNChildren(prodexpr);
   assert(nchildren >= 2);

   /* check whether there is already an expression that represents the product */
   if( SCIPhashmapExists(exprmap, (void*)prodexpr) )
   {
      *newexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapGetImage(exprmap, (void*)prodexpr);
      assert(*newexpr != NULL);

      /* capture expression */
      SCIPcaptureConsExprExpr(*newexpr);
   }
   else
   {
      SCIPdebugMsg(scip, "  product expression %p has been considered for the first time\n", (void*)prodexpr);

      if( nchildren == 2 )
      {
         SCIP_CLIQUE** xcliques;
         SCIP_VAR* x;
         SCIP_VAR* y;
         SCIP_Bool found_clique = FALSE;
         int c;

         /* get variables from the product expression */
         x = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(prodexpr)[0]);
         assert(x != NULL);
         y = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(prodexpr)[1]);
         assert(y != NULL);
         assert(x != y);

         /* first try to find a clique containing both variables */
         xcliques = SCIPvarGetCliques(x, TRUE);

         /* look in cliques containing x */
         for( c = 0; c < SCIPvarGetNCliques(x, TRUE); ++c )
         {
            if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* x + y <= 1 => x*y = 0 */
            {
               /* create zero value expression */
               SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, newexpr, 0.0) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 1;

               found_clique = TRUE;
               break;
            }

            if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* x + (1-y) <= 1 => x*y = x */
            {
               /* create variable expression for x */
               SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, newexpr, x) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 2;

               found_clique = TRUE;
               break;
            }
         }

         if( !found_clique )
         {
            xcliques = SCIPvarGetCliques(x, FALSE);

            /* look in cliques containing complement of x */
            for( c = 0; c < SCIPvarGetNCliques(x, FALSE); ++c )
            {
               if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* (1-x) + y <= 1 => x*y = y */
               {
                  /* create variable expression for y */
                  SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, newexpr, y) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 1;

                  found_clique = TRUE;
                  break;
               }

               if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* (1-x) + (1-y) <= 1 => x*y = x + y - 1 */
               {
                  /* create sum expression */
                  SCIP_CONSEXPR_EXPR* sum_children[2];
                  SCIP_Real sum_coefs[2];
                  SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &sum_children[0], x) );
                  SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &sum_children[1], y) );
                  sum_coefs[0] = 1.0;
                  sum_coefs[1] = 1.0;
                  SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, newexpr, 2, sum_children, sum_coefs, -1.0) );

                  SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sum_children[0]) );
                  SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sum_children[1]) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 3;

                  found_clique = TRUE;
                  break;
               }
            }
         }

         /* if the variables are not in a clique, do standard linearization */
         if( !found_clique )
         {
            SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss,
               conshdlrdata->reformbinprodsand) );
         }
      }
      else
      {
         /* linearize binary product using an AND constraint because nchildren > 2 */
         SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss,
            conshdlrdata->reformbinprodsand) );
      }

      /* hash variable expression */
      SCIP_CALL( SCIPhashmapInsert(exprmap, (void*)prodexpr, *newexpr) );
   }

   return SCIP_OKAY;
}

/** helper function to replace binary products in a given expression constraints */
static
SCIP_RETCODE replaceBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_CONSEXPR_ITERATOR* it,               /**< expression iterator */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(exprmap != NULL);
   assert(it != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   SCIPdebugMsg(scip, "  check constraint %s\n", SCIPconsGetName(cons));

   for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      SCIP_CONSEXPR_EXPR* newexpr = NULL;
      SCIP_CONSEXPR_EXPR* childexpr;
      int childexpridx;

      childexpridx = SCIPexpriteratorGetChildIdxDFS(it);
      assert(childexpridx >= 0 && childexpridx < SCIPgetConsExprExprNChildren(expr));
      childexpr = SCIPexpriteratorGetChildExprDFS(it);
      assert(childexpr != NULL);

      /* try to factorize variables in a sum expression that contains several products of binary variables */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, childexpr,
            conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* try to create an expression that represents a product of binary variables */
      if( newexpr == NULL )
      {
         SCIP_CALL( getBinaryProductExpr(scip, conshdlr, exprmap, childexpr, &newexpr, naddconss, nchgcoefs) );
      }

      if( newexpr != NULL )
      {
         assert(naddconss == NULL || *naddconss > 0 || nchgcoefs == NULL || *nchgcoefs > 0);

         /* replace product expression */
         SCIP_CALL( SCIPreplaceConsExprExprChild(scip, expr, childexpridx, newexpr) );

         /* note that the expression has been captured by getBinaryProductExpr and SCIPreplaceConsExprExprChild */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &newexpr) );

         /* mark the constraint to not be simplied anymore */
         consdata->issimplified = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** reformulates products of binary variables during presolving in the following way:
 *
 * Let sum_{i,j} Q_ij x_i x_j be a subexpression that only contains binary variables. Each term x_i x_j is
 * reformulated with the help of an extra (implicit integer) variable z_ij in {0,1}:
 *
 *    z_ij <= x_i, z_ij <= x_j, x_i + x_j - z_ij <= 1
 *
 * Before reformulating x_i x_j in this way, it is checked whether there is a clique that contains x_i and x_j. These
 * cliques allows for a better reformulation. There are four cases:
 *
 *    1. x_i + x_j <= 1 implies that x_i x_j = 0
 *
 *    2. x_i + (1 - x_j) <= 1 implies x_i x_j = x_i
 *
 *    3. (1 - x_i) + x_j <= 1 implies x_i x_j = x_j
 *
 *    4. (1 - x_i) + (1 - x_j) <= 1 implies x_i x_j = x_i + x_j - 1
 *
 * The reformulation using z_ij or the cliques is implemented in getBinaryProductExpr().
 *
 * Introducing too many extra variables and constraints can have a negative impact on the performance (e.g., due to
 * slow probing). For this reason, it is checked in getFactorizedBinaryQuadraticExpr() whether sum_{i,j} Q_ij x_i x_j
 * contains large (>= reformbinprodsfac parameter) lower sums of the form x_i sum_{j} Q_ij x_j. Such a lower sum is
 * reformulated with only one extra variable w_i:
 *
 *    maxact := sum_j max{0, Q_ij}, minact := sum_j min{0, Q_ij}
 *    minact x_i <= w_i, w_i <= maxact x_i
 *    minact <= sum_j Q_ij x_j - w_i + minact x_i
 *    maxact >= sum_j Q_ij x_j - w_i + maxact x_i
 *
 * We mark w_i to be implicit integer if all Q_ij are integer. After each replacment of a lower sum, it
 * is checked whether there are enough terms left to factorize other binary variables. Lower sums with a larger number
 * of terms are prioritized.
 */
static
SCIP_RETCODE presolveBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< expression constraints */
   int                   nconss,             /**< total number of expression constraints */
   int*                  naddconss,          /**< pointer to store the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to store the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_HASHMAP* exprmap;
   SCIP_CONSEXPR_ITERATOR* it;
   int c;

   assert(conshdlr != NULL);

   /* no expression constraints or binary variables -> skip */
   if( nconss == 0 || SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;
   assert(conss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create expression hash map */
   SCIP_CALL( SCIPhashmapCreate(&exprmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* create expression iterator */
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD);

   SCIPdebugMsg(scip, "call presolveBinaryProducts()\n");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONSEXPR_EXPR* newexpr = NULL;

      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* try to reformulate the root expression */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, conss[c], consdata->expr,
            conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* release the root node if another expression has been found */
      if( newexpr != NULL )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );
         consdata->expr = newexpr;

         /* mark constraint to be not simplified anymore */
         consdata->issimplified = FALSE;
      }

      /* replace each product of binary variables separately */
      SCIP_CALL( replaceBinaryProducts(scip, conshdlr, conss[c], exprmap, it, naddconss, nchgcoefs) );
   }

   /* free memory */
   SCIPhashmapFree(&exprmap);
   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** simplifies expressions and replaces common subexpressions for a set of constraints
 * @todo put the constant to the constraint sides
 */
static
SCIP_RETCODE canonicalizeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< total number of constraints */
   SCIP_PRESOLTIMING     presoltiming,       /**< presolve timing (SCIP_PRESOLTIMING_ALWAYS if not in presolving) */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   int*                  ndelconss,          /**< counter to add number of deleted constraints, or NULL */
   int*                  naddconss,          /**< counter to add number of added constraints, or NULL */
   int*                  nchgcoefs           /**< counter to add number of changed coefficients, or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it;
   int* nlockspos;
   int* nlocksneg;
   SCIP_Bool havechange;
   SCIP_Bool reformulate = FALSE;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* update number of canonizalize calls */
   ++(conshdlrdata->ncanonicalizecalls);

   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->canonicalizetime) );

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );

   *infeasible = FALSE;

   /* check whether at least one nonlinear handler implements the reformulation callback */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      if( SCIPhasConsExprNlhdlrReformulate(conshdlrdata->nlhdlrs[i]) )
      {
         reformulate = TRUE;
         break;
      }
   }

   /* set havechange to TRUE in the first call of canonicalize; otherwise we might not replace common subexpressions */
   havechange = conshdlrdata->ncanonicalizecalls == 1;

   /* free nonlinear handlers information from expressions */  /* TODO can skip this in first presolve round */
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSEXPR_EXPR* expr;

      assert(conss != NULL);
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIPdebugMsg(scip, "free nonlinear handler data for expression %p\n", (void*)expr);

         assert(expr->auxvar == NULL);  /* should not have been created yet or have been removed in INITPRE (if restart) */

         /* remove nonlinear handlers in expression and their data */
         SCIP_CALL( freeEnfoData(scip, conshdlr, expr, FALSE) );
      }
   }

   /* allocate memory for storing locks of each constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   /* unlock all constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* remember locks */
      nlockspos[i] = consdata->nlockspos;
      nlocksneg[i] = consdata->nlocksneg;

      /* remove locks */
      SCIP_CALL( addLocks(scip, conss[i], -consdata->nlockspos, -consdata->nlocksneg) );
      assert(consdata->nlockspos == 0);
      assert(consdata->nlocksneg == 0);
   }

#ifndef NDEBUG
   /* check whether all locks of each expression have been removed */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSEXPR_EXPR* expr;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      SCIP_CALL( SCIPexpriteratorInit(it, consdata->expr, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, TRUE) );
      for( expr = consdata->expr; !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         assert(expr != NULL);
         assert(expr->nlocksneg == 0);
         assert(expr->nlockspos == 0);
      }
   }
#endif

   /* reformulate products of binary variables */
   if( conshdlrdata->reformbinprods && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING
      && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      int tmpnaddconss = 0;
      int tmpnchgcoefs = 0;

      /* call this function before simplification because expressions might not be simplified after reformulating
       * binary products; the detection of some nonlinear handlers might assume that expressions are simplified
       */
      SCIP_CALL( presolveBinaryProducts(scip, conshdlr, conss, nconss, &tmpnaddconss, &tmpnchgcoefs) );

      /* update counters */
      if( naddconss != NULL )
         *naddconss = tmpnaddconss;
      if( nchgcoefs != NULL )
         *nchgcoefs = tmpnchgcoefs;

      /* check whether at least one expression has changed */
      if( tmpnaddconss + tmpnchgcoefs > 0 )
         havechange = TRUE;
   }

   for( i = 0; i < nconss; ++i )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* call simplify for each expression */
      if( !consdata->issimplified && consdata->expr != NULL )
      {
         SCIP_CONSEXPR_EXPR* simplified;
         SCIP_Bool changed;

         changed = FALSE;
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, consdata->expr, &simplified, &changed, infeasible) );
         consdata->issimplified = TRUE;

         if( changed )
            havechange = TRUE;

         /* If root expression changed, then we need to take care updating the locks as well (the consdata is the one holding consdata->expr "as a child").
          * If root expression did not change, some subexpression may still have changed, but the locks were taking care of in the corresponding SCIPreplaceConsExprExprChild() call.
          */
         if( simplified != consdata->expr )
         {
            assert(changed);

            /* release old expression */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

            /* store simplified expression */
            consdata->expr = simplified;
         }
         else
         {
            /* The simplify captures simplified in any case, also if nothing has changed.
             * Therefore, we have to release it here.
             */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );
         }

         if( *infeasible )
            break;

         /* scale constraint sides */
         SCIP_CALL( scaleConsSides(scip, conshdlr, conss[i], &changed) );

         if( changed )
            havechange = TRUE;

         /* handle constant root expression; either the problem is infeasible or the constraint is redundant */
         if( consdata->expr->exprhdlr == SCIPgetConsExprExprHdlrValue(conshdlr) )
         {
            SCIP_Real value = SCIPgetConsExprExprValueValue(consdata->expr);
            if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasNegative(scip, value - consdata->lhs)) ||
                (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasPositive(scip, value - consdata->rhs)) )
            {
               SCIPdebugMsg(scip, "<%s> with constant expression found infeasible\n", SCIPconsGetName(conss[i]));
               SCIPdebugPrintCons(scip, conss[i], NULL);
               *infeasible = TRUE;
               break;
            }
            else
            {
               SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
               SCIP_CALL( SCIPdelCons(scip, conss[i]) );
               if( ndelconss != NULL )
                  ++*ndelconss;
               havechange = TRUE;
            }
         }
      }

      /* call reformulation callback of nonlinear handlers for each expression */
      if( reformulate && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      {
         SCIP_CONSEXPR_EXPR* refexpr;
         SCIP_Bool changed;

         if( consdata->expr != NULL )
         {
            SCIP_CALL( SCIPreformulateConsExprExpr(scip, conshdlr, consdata->expr, &refexpr, &changed, infeasible) );

            if( changed )
               havechange = TRUE;

            if( refexpr != consdata->expr )
            {
               assert(changed);

               /* release old expression */
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

               /* store reformulated expression */
               consdata->expr = refexpr;
            }
            else
            {
               /* The reformulation captures refexpr in any case, also if nothing has changed.
                * Therefore, we have to release it here.
                */
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &refexpr) );
            }
         }

         if( *infeasible )
            break;
      }
   }

   /* replace common subexpressions */
   if( havechange && !*infeasible )
   {
      SCIP_CONS** consssorted;

      SCIP_CALL( replaceCommonSubexpressions(scip, conss, nconss) );

      /* FIXME: this is a dirty hack for updating the variable expressions stored inside an expression which might have
       * been changed after simplification; now we completely recollect all variable expression and variable events
       */

      /* Each variable stores the constraints for which it catched varbound events sorted by the constraint index.
       * Thus, for performance reasons, it is better to call dropVarEvents in descending order of constraint index.
       */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortPtr((void**)consssorted, SCIPcompareConsExprIndex, nconss);

      for( i = nconss-1; i >= 0; --i )
      {
         assert(i == 0 || SCIPcompareConsExprIndex((void*)consssorted[i-1], (void*)consssorted[i]) < 0);
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }
      for( i = 0; i < nconss; ++i )
      {
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(consssorted[i])) );
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
      }

      SCIPfreeBufferArray(scip, &consssorted);
   }

   /* restore locks */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsDeleted(conss[i]) )
         continue;

      SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
   }

   /* run nlhdlr detect if in presolving stage (that is, not in exitpre) */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !*infeasible )
   {
      SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss, infeasible) );
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);
   SCIPexpriteratorFree(&it);

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->canonicalizetime) );

   return SCIP_OKAY;
}

/** @name Parsing methods
 * @{
 * Here is an attempt at defining the grammar of an expression.
 * We use upper case names for variables (in the grammar sense) and terminals are between "".
 * Loosely speaking, a Base will be any "block", a Factor is a Base to a power, a Term is a product of Factors
 * and an Expression is a sum of terms.
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * parse(Expr|Term|Base) returns an SCIP_CONSEXPR_EXPR
 *
 * @todo We can change the grammar so that Factor becomes base and we allow a Term to be
 *       <pre> Term       -> Factor { ("*" | "/" | "^") Factor } </pre>
 */

#ifdef PARSE_DEBUG
#define debugParse                      printf
#else
#define debugParse                      while( FALSE ) printf
#endif
static
SCIP_RETCODE parseExpr(SCIP*, SCIP_CONSHDLR*, SCIP_HASHMAP*, const char*, const char**, SCIP_CONSEXPR_EXPR**);

/** Parses base to build a value, variable, sum, or function-like ("func(...)") expression.
 * <pre>
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 */
static
SCIP_RETCODE parseBase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between SCIP vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  basetree            /**< buffer to store the expr parsed by Base */
   )
{
   SCIP_VAR* var;

   debugParse("parsing base from %s\n", expr); /*lint !e506 !e681*/

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string\n");
      return SCIP_READERROR;
   }

   if( *expr == '<' )
   {
      /* parse a variable */
      SCIP_CALL( SCIPparseVarName(scip, expr, &var, (char**)newpos) );

      if( var == NULL )
      {
         SCIPerrorMessage("Could not find variable with name '%s'\n", expr);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* check if we have already created an expression out of this var */
      if( SCIPhashmapExists(vartoexprvarmap, (void *)var) )
      {
         debugParse("Variable %s has been parsed, capturing its expression\n", SCIPvarGetName(var)); /*lint !e506 !e681*/
         *basetree = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(vartoexprvarmap, (void *)var);
         SCIPcaptureConsExprExpr(*basetree);
      }
      else
      {
         debugParse("First time parsing variable %s, creating varexpr and adding it to hashmap\n", SCIPvarGetName(var)); /*lint !e506 !e681*/
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, basetree, var) );
         SCIP_CALL( SCIPhashmapInsert(vartoexprvarmap, (void*)var, (void*)(*basetree)) );
      }
   }
   else if( *expr == '(' )
   {
      /* parse expression */
      SCIP_CALL( parseExpr(scip, conshdlr, vartoexprvarmap, ++expr, newpos, basetree) );
      expr = *newpos;

      /* expect ')' */
      if( *expr != ')' )
      {
         SCIPerrorMessage("Read a '(', parsed expression inside --> expecting closing ')'. Got <%c>: rest of string <%s>\n", *expr, expr);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, basetree) );
         return SCIP_READERROR;
      }
      ++expr;
      debugParse("Done parsing expression, continue with <%s>\n", expr); /*lint !e506 !e681*/
   }
   else if( isdigit(*expr) )
   {
      /* parse number */
      SCIP_Real value;
      if( !SCIPstrToRealValue(expr, &value, (char**)&expr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", expr);
         return SCIP_READERROR;
      }
      debugParse("Parsed value %g, creating a value-expression.\n", value); /*lint !e506 !e681*/
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, basetree, value) );
   }
   else if( isalpha(*expr) )
   {
      /* a (function) name is coming, should find exprhandler with such name */
      int i;
      char operatorname[SCIP_MAXSTRLEN];
      SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
      SCIP_Bool success;

      /* get name */
      i = 0;
      while( *expr != '(' && !isspace((unsigned char)*expr) && *expr != '\0' )
      {
         operatorname[i] = *expr;
         ++expr;
         ++i;
      }
      operatorname[i] = '\0';

      /* after name we must see a '(' */
      if( *expr != '(' )
      {
         SCIPerrorMessage("Expected '(' after operator name <%s>, but got %s.\n", operatorname, expr);
         return SCIP_READERROR;
      }

      /* search for expression handler */
      exprhdlr = SCIPfindConsExprExprHdlr(conshdlr, operatorname);

      /* check expression handler exists and has a parsing method */
      if( exprhdlr == NULL )
      {
         SCIPerrorMessage("No expression handler with name <%s> found.\n", operatorname);
         return SCIP_READERROR;
      }

      ++expr;
      SCIP_CALL( SCIPparseConsExprExprHdlr(scip, conshdlr, exprhdlr, expr, newpos, basetree, &success) );

      if( !success )
      {
         SCIPerrorMessage("Error while expression handler <%s> was parsing %s\n", operatorname, expr);
         assert(*basetree == NULL);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* we should see the ')' of Op "(" OpExpression ") */
      assert(*expr == ')');

      /* move one character forward */
      ++expr;
   }
   else
   {
      /* Base -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ") */
      SCIPerrorMessage("Expected a number, (expression), <varname>, Opname(Opexpr), instead got <%c> from %s\n", *expr, expr);
      return SCIP_READERROR;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses a factor and builds a product-expression if there is an exponent, otherwise returns the base expression.
 * <pre>
 * Factor -> Base [ "^" "number" | "^(" "number" ")" ]
 * </pre>
 */
static
SCIP_RETCODE parseFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             isdenominator,      /**< whether factor is in the denominator */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  factortree          /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_CONSEXPR_EXPR*  basetree;
   SCIP_Real exponent;

   debugParse("parsing factor from %s\n", expr); /*lint !e506 !e681*/

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string.\n");
      return SCIP_READERROR;
   }

   /* parse Base */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseBase(scip, conshdlr, vartoexprvarmap, expr, newpos, &basetree) );
   expr = *newpos;

   /* check if there is an exponent */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '^' )
   {

      ++expr;
      while( isspace((unsigned char)*expr) )
         ++expr;

      if( *expr == '\0' )
      {
         SCIPerrorMessage("Unexpected end of expression string after '^'.\n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
         return SCIP_READERROR;
      }

      if( *expr == '(' )
      {
         ++expr;

         /* it is exponent with parenthesis; expect number possibly starting with + or - */
         if( !SCIPstrToRealValue(expr, &exponent, (char**)&expr) )
         {
            SCIPerrorMessage("error parsing number from <%s>\n", expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }

         /* expect the ')' */
         while( isspace((unsigned char)*expr) )
            ++expr;
         if( *expr != ')' )
         {
            SCIPerrorMessage("error in parsing exponent: expected ')', received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
         ++expr;
      }
      else
      {
         /* no parenthesis, we should see just a positive number */

         /* expect a digit */
         if( isdigit(*expr) )
         {
            if( !SCIPstrToRealValue(expr, &exponent, (char**)&expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected a digit, received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
      }

      debugParse("parsed the exponent %g\n", exponent); /*lint !e506 !e681*/
   }
   else
   {
      /* there is no explicit exponent */
      exponent = 1.0;
   }
   *newpos = expr;

   /* multiply with -1 when we are in the denominator */
   if( isdenominator )
      exponent *= -1.0;

   /* create power */
   if( exponent != 1.0 )
   {
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, factortree, basetree, exponent) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
   }
   else
      /* Factor consists of this unique Base */
      *factortree = basetree;

   return SCIP_OKAY;
}

/** Parses a term and builds a product-expression, where each factor is a child.
 * <pre>
 * Term -> Factor { ("*" | "/" ) Factor }
 * </pre>
 */
static
SCIP_RETCODE parseTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  termtree            /**< buffer to store the expr parsed by Term */
   )
{
   SCIP_CONSEXPR_EXPR* factortree;

   debugParse("parsing term from %s\n", expr); /*lint !e506 !e681*/

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(scip, conshdlr, FALSE, vartoexprvarmap, expr, newpos, &factortree) );
   expr = *newpos;

   debugParse("back to parsing Term, continue parsing from %s\n", expr); /*lint !e506 !e681*/

   /* check if Terms has another Factor incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '*' || *expr == '/' )
   {
      /* initialize termtree as a product expression with a single term, so we can append the extra Factors */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 1, &factortree, 1.0) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factortree) );

      /* loop: parse Factor, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Bool isdivision;

         isdivision = (*expr == '/') ? TRUE : FALSE;

         debugParse("while parsing term, read char %c\n", *expr); /*lint !e506 !e681*/

         ++expr;
         retcode = parseFactor(scip, conshdlr, isdivision, vartoexprvarmap, expr, newpos, &factortree);

         /* release termtree, if parseFactor fails with a read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, termtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created factor */
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, *termtree, factortree) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factortree) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '*' || *expr == '/' );
   }
   else
   {
      /* Term consists of this unique factor */
      *termtree = factortree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses an expression and builds a sum-expression with children.
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * </pre>
 */
static
SCIP_RETCODE parseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  exprtree            /**< buffer to store the expr parsed by Expr */
   )
{
   SCIP_Real sign;
   SCIP_CONSEXPR_EXPR* termtree;

   debugParse("parsing expression %s\n", expr); /*lint !e506 !e681*/

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   /* if '+' or '-', store it */
   sign = 1.0;
   if( *expr == '+' || *expr == '-' )
   {
      debugParse("while parsing expression, read char %c\n", *expr); /*lint !e506 !e681*/
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree) );
   expr = *newpos;

   debugParse("back to parsing expression (we have the following term), continue parsing from %s\n", expr); /*lint !e506 !e681*/

   /* check if Expr has another Term incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '+' || *expr == '-' )
   {
      if( SCIPgetConsExprExprHdlr(termtree) == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         /* initialize exprtree as a sum expression with a constant only, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 0, NULL, NULL, sign * SCIPgetConsExprExprValueValue(termtree)) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }
      else
      {
         /* initialize exprtree as a sum expression with a single term, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }

      /* loop: parse Term, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Real coef;

         /* check if we have a "coef * <term>" */
         if( SCIPstrToRealValue(expr, &coef, (char**)newpos) )
         {
            while( isspace((unsigned char)**newpos) )
               ++(*newpos);

            if( **newpos != '*' )
            {
               /* no '*', so fall back to parsing term after sign */
               coef = (*expr == '+') ? 1.0 : -1.0;
               ++expr;
            }
            else
            {
               /* keep coefficient in coef and continue parsing term after coefficient */
               expr = (*newpos)+1;

               while( isspace((unsigned char)*expr) )
                  ++expr;
            }
         }
         else
         {
            coef = (*expr == '+') ? 1.0 : -1.0;
            ++expr;
         }

         debugParse("while parsing expression, read coefficient %g\n", coef); /*lint !e506 !e681*/

         retcode = parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree);

         /* release exprtree if parseTerm fails with an read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, exprtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created term */
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *exprtree, termtree, coef) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '+' || *expr == '-' );
   }
   else
   {
      /* Expr consists of this unique ['+' | '-'] Term */
      if( sign  < 0.0 )
      {
         assert(sign == -1.0);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }
      else
         *exprtree = termtree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** @} */  /* end of parsing methods */

/** given a cons_expr expression, creates an equivalent classic (nlpi-) expression */
static
SCIP_RETCODE makeClassicExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   sourceexpr,         /**< expression to convert */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to created expression */
   SCIP_CONSEXPR_EXPR**  varexprs,           /**< variable expressions that might occur in expr, their position in this array determines the varidx */
   int                   nvarexprs           /**< number of variable expressions */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_EXPR** children = NULL;
   int nchildren;
   int c;

   assert(scip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(sourceexpr);
   nchildren = SCIPgetConsExprExprNChildren(sourceexpr);

   /* collect children expressions from children, if any */
   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );
      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( makeClassicExpr(scip, SCIPgetConsExprExprChildren(sourceexpr)[c], &children[c], varexprs, nvarexprs) );
         assert(children[c] != NULL);
      }
   }

   /* create target expression */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "var") == 0 )
   {
      int varidx;

      /* find variable expression in varexprs array
       * the position in the array determines the index of the variable in the classic expression
       * TODO if varexprs are sorted, then can do this more efficient
       */
      for( varidx = 0; varidx < nvarexprs; ++varidx )
         if( varexprs[varidx] == sourceexpr )
            break;
      assert(varidx < nvarexprs);

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_VARIDX, varidx) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "val") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_CONST, SCIPgetConsExprExprValueValue(sourceexpr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "sum") == 0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, nchildren, children, SCIPgetConsExprExprSumCoefs(sourceexpr), SCIPgetConsExprExprSumConstant(sourceexpr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "pow") == 0 )
   {
      SCIP_Real exponent;

      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);

      exponent = SCIPgetConsExprExprPowExponent(sourceexpr);
      if( EPSISINT(exponent, 0.0) )  /*lint !e835*/
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_INTPOWER, *children, (int)exponent) );
      }
      else
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_REALPOWER, *children, exponent) );
      }
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "signpower") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_SIGNPOWER, *children,
         SCIPgetConsExprExprPowExponent(sourceexpr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "prod") == 0 )
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, SCIPgetConsExprExprProductCoef(sourceexpr), nchildren, NULL, NULL) );
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), targetexpr, nchildren, children, 1, &monomial, 0.0, FALSE) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "abs") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_ABS, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "exp") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_EXP, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "log") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_LOG, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "sin") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_SIN, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "cos") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_COS, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "entropy") == 0 )
   {
      SCIP_EXPR* childcopy;
      SCIP_Real minusone = -1.0;

      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &childcopy, children[0]) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &childcopy, SCIP_EXPR_LOG, childcopy) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_MUL, children[0], childcopy) );
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, 1, targetexpr, &minusone, 0.0) );
   }
   else
   {
      SCIPerrorMessage("unsupported expression handler <%s>, cannot convert to classical expression\n", SCIPgetConsExprExprHdlrName(exprhdlr));
      return SCIP_ERROR;
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** create a nonlinear row representation of an expr constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< expression constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* classicexpr = NULL;
   SCIP_VAR** nlvars = NULL;
   int nnlvars = 0;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* @todo pass correct curvature */
   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
         0, NULL, NULL, 0, NULL, 0, NULL, NULL, consdata->lhs, consdata->rhs, SCIP_EXPRCURV_UNKNOWN) );

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   if( SCIPgetConsExprExprHdlr(consdata->expr) == conshdlrdata->exprsumhdlr )
   {
      /* if root is a sum, then split into linear, quadratic, and expression */
      SCIP_CONSEXPR_EXPR* child;
      SCIP_Real* coefs;

      /* constant term of sum */
      SCIP_CALL( SCIPchgNlRowConstant(scip, consdata->nlrow, SCIPgetConsExprExprSumConstant(consdata->expr)) );

      coefs = SCIPgetConsExprExprSumCoefs(consdata->expr);

      for( i = 0; i < SCIPgetConsExprExprNChildren(consdata->expr); ++i )
      {
         child = SCIPgetConsExprExprChildren(consdata->expr)[i];

         if( SCIPisConsExprExprVar(child) )
         {
            /* linear term */
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(child), coefs[i]) );
         }
         else if( SCIPgetConsExprExprHdlr(child) == conshdlrdata->exprpowhdlr &&
            SCIPgetConsExprExprPowExponent(child) == 2.0 &&
            SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(child)[0]) )
         {
            /* square term  */
            SCIP_QUADELEM quadelem;

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0])) );
               quadelem.idx1 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }
            quadelem.idx2 = quadelem.idx1;
            quadelem.coef = coefs[i];

            SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
         }
         else if( SCIPgetConsExprExprHdlr(child) == conshdlrdata->exprprodhdlr &&
            SCIPgetConsExprExprNChildren(child) == 2 &&
            SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(child)[0]) &&
            SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(child)[1]) )
         {
            /* bilinear term */
            SCIP_QUADELEM quadelem;

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0])) );
               quadelem.idx1 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }

            quadelem.idx2 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[1]));
            if( quadelem.idx2 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[1])) );
               quadelem.idx2 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }

            quadelem.coef = coefs[i];

            SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
         }
         else
         {
            /* general nonlinear term */
            SCIP_EXPR* classicchild;

            /* make classic expression of child i */
            SCIP_CALL( makeClassicExpr(scip, child, &classicchild, consdata->varexprs, consdata->nvarexprs) );

            /* create or extend classicexpr */
            if( classicexpr == NULL )
            {
               SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &classicexpr, 1, &classicchild, coefs + i, 0.0) );
            }
            else
            {
               SCIP_CALL( SCIPexprAddToLinear(SCIPblkmem(scip), classicexpr, 1, &coefs[i], &classicchild, 0.0) );
            }
         }
      }

      if( classicexpr != NULL )
      {
         /* reindex variables in classicexpr so that only used variables are left */
         int* varsusage;
         int* reindexvars;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &nlvars, consdata->nvarexprs) );
         SCIP_CALL( SCIPallocBufferArray(scip, &reindexvars, consdata->nvarexprs) );
         SCIP_CALL( SCIPallocClearBufferArray(scip, &varsusage, consdata->nvarexprs) );

         /* get count how often variables are used in expr */
         SCIPexprGetVarsUsage(classicexpr, varsusage);

         /* sort out unused variables and collect and reindex remaining variables */
         nnlvars = 0;
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            if( varsusage[i] == 0 )
            {
               reindexvars[i] = -1;
            }
            else
            {
               reindexvars[i] = nnlvars;
               nlvars[nnlvars] = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
               ++nnlvars;
            }
         }

         SCIPexprReindexVars(classicexpr, reindexvars);

         SCIPfreeBufferArray(scip, &varsusage);
         SCIPfreeBufferArray(scip, &reindexvars);
      }
   }
   else if( SCIPgetConsExprExprHdlr(consdata->expr) == conshdlrdata->exprpowhdlr &&
      SCIPgetConsExprExprPowExponent(consdata->expr) == 2.0 &&
      SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(consdata->expr)[0]) )
   {
      /* if root is a x^2, then set the quadratic part of the nlrow */
      SCIP_QUADELEM quadelem;

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(consdata->expr)[0])) );
      quadelem.idx1 = 0;
      quadelem.idx2 = 0;
      quadelem.coef = 1.0;

      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
   }
   else if( SCIPgetConsExprExprHdlr(consdata->expr) == conshdlrdata->exprprodhdlr &&
      SCIPgetConsExprExprNChildren(consdata->expr) == 2 &&
      SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(consdata->expr)[0]) &&
      SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(consdata->expr)[1]) )
   {
      /* if root is a bilinear term x*y, then set the quadratic part of the nlrow */
      SCIP_QUADELEM quadelem;

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(consdata->expr)[0])) );
      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(consdata->expr)[1])) );

      quadelem.idx1 = 0;
      quadelem.idx2 = 1;
      quadelem.coef = 1.0;

      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
   }
   else
   {
      /* make classic expression */
      SCIP_CALL( makeClassicExpr(scip, consdata->expr, &classicexpr, consdata->varexprs, consdata->nvarexprs) );

      /* collect variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &nlvars, consdata->nvarexprs) );

      nnlvars = consdata->nvarexprs;
      for( i = 0; i < consdata->nvarexprs; ++i )
         nlvars[i] = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
   }
   assert((classicexpr != NULL) == (nlvars != NULL));

   if( classicexpr != NULL )
   {
      /* make classic expression tree */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, classicexpr, nnlvars, 0, NULL) );

      /* set variables in expression tree */
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, nnlvars, nlvars) );
      SCIPfreeBufferArray(scip, &nlvars);

      /* add expression tree in nlrow (this will make a copy) */
      SCIP_CALL( SCIPsetNlRowExprtree(scip, consdata->nlrow, exprtree) );

      /* free exprtree */
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** compares nonlinear handler by priority
 *
 * if handlers have same priority, then compare by name
 */
static
int nlhdlrCmp(
   void*                 hdlr1,              /**< first handler */
   void*                 hdlr2               /**< second handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* h1;
   SCIP_CONSEXPR_NLHDLR* h2;

   assert(hdlr1 != NULL);
   assert(hdlr2 != NULL);

   h1 = (SCIP_CONSEXPR_NLHDLR*)hdlr1;
   h2 = (SCIP_CONSEXPR_NLHDLR*)hdlr2;

   if( h1->priority != h2->priority )
      return (int)(h1->priority - h2->priority);

   return strcmp(h1->name, h2->name);
}

/** frees auxiliary variables which have been added to compute an outer approximation */
static
SCIP_RETCODE freeAuxVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSDATA* consdata;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL);
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr == NULL )
         continue;

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIP_CALL( freeAuxVar(scip, conshdlr, expr) );
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** calls separation initialization callback for each expression */
static
SCIP_RETCODE initSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether the problem is infeasible or not */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   int c, e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(infeasible != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   *infeasible = FALSE;
   for( c = 0; c < nconss && !*infeasible; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      /* call separation initialization callback for 'initial' constraints only */
      if( !SCIPconsIsInitial(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it) && !*infeasible; expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         /* call initsepa of all nlhdlrs in expr */
         for( e = 0; e < expr->nenfos; ++e )
         {
            SCIP_CONSEXPR_NLHDLR* nlhdlr;
            SCIP_Bool underestimate;
            SCIP_Bool overestimate;
            assert(expr->enfos[e] != NULL);

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* only init sepa if there is an initsepa callback */
            if( !SCIPhasConsExprNlhdlrInitSepa(nlhdlr) )
               continue;

            assert(!expr->enfos[e]->issepainit);

            /* check whether expression needs to be under- or overestimated */
            overestimate = SCIPgetConsExprExprNLocksNeg(expr) > 0;
            underestimate = SCIPgetConsExprExprNLocksPos(expr) > 0;
            assert(underestimate || overestimate);

            /* call the separation initialization callback of the nonlinear handler */
            SCIP_CALL( SCIPinitsepaConsExprNlhdlr(scip, conshdlr, conss[c], nlhdlr, expr,
               expr->enfos[e]->nlhdlrexprdata, overestimate, underestimate, infeasible) );
            expr->enfos[e]->issepainit = TRUE;

            if( *infeasible )
            {
               /* stop everything if we detected infeasibility */
               SCIPdebugMsg(scip, "detect infeasibility for constraint %s during initsepa()\n", SCIPconsGetName(conss[c]));
               break;
            }
         }
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** gets weight of variable when splitting violation score onto several variables in an expression */
static
SCIP_Real getViolSplitWeight(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_VAR*               var,              /**< variable */
   SCIP_SOL*               sol               /**< current solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   switch( conshdlrdata->branchviolsplit )
   {
      case 'e' :  /* evenly: everyone gets the same score */
         return 1.0;

      case 'm' :  /* midness of solution: 0.5 if in middle of domain, 0.05 if close to lower or upper bound */
      {
         SCIP_Real weight;
         weight = MIN(SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var) - SCIPgetSolVal(scip, sol, var)) / (SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var)); /*lint !e666*/
         return MAX(0.05, weight);
      }

      case 'd' :  /* domain width */
         return SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

      case 'l' :  /* logarithmic domain width: log-scale if width is below 0.1 or above 10, otherwise actual width */
      {
         SCIP_Real width = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
         assert(width > 0.0);
         if( width > 10.0 )
            return 10.0*log10(width);
         if( width < 0.1 )
            return 0.1/(-log10(width));
         return width;
      }

      default :
         SCIPerrorMessage("invalid value for parameter constraints/expr/branching/violsplit");
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** adds violation-branching score to a set of expressions, thereby distributing the score
 *
 * Each expression must either be a variable expression or have an aux-variable.
 *
 * If unbounded variables are present, each unbounded var gets an even score.
 * If no unbounded variables, then parameter constraints/expr/branching/violsplit decides weight for each var.
 */
static
void addConsExprExprsViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation-branching score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_VAR* var;
   SCIP_Real weight;
   SCIP_Real weightsum = 0.0; /* sum of weights over all candidates with bounded domain */
   int nunbounded = 0;  /* number of candidates with unbounded domain */
   int i;

   assert(exprs != NULL);
   assert(nexprs > 0);
   assert(success != NULL);

   if( nexprs == 1 )
   {
      SCIPaddConsExprExprViolScore(scip, conshdlr, exprs[0], violscore);
      *success = TRUE;
      return;
   }

   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetConsExprExprAuxVar(exprs[i]);
      assert(var != NULL);

      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         ++nunbounded;
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         weightsum += getViolSplitWeight(scip, conshdlr, var, sol);
   }

   *success = FALSE;
   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetConsExprExprAuxVar(exprs[i]);
      assert(var != NULL);

      if( nunbounded > 0 )
      {
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIPaddConsExprExprViolScore(scip, conshdlr, exprs[i], violscore / nunbounded);
            *success = TRUE;
         }
      }
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         assert(weightsum > 0.0);

         weight = getViolSplitWeight(scip, conshdlr, var, sol);
         SCIPaddConsExprExprViolScore(scip, conshdlr, exprs[i], violscore * weight / weightsum);
         SCIPdebugMsg(scip, "add score %g (%g%% of %g) to <%s>[%g,%g]\n", violscore * weight / weightsum,
            100*weight / weightsum, violscore,
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         *success = TRUE;
      }
   }
}

/** adds violation-branching score to children of expression for given auxiliary variables
 *
 * Iterates over the successors of expr to find expressions that are associated with one of the given auxiliary variables.
 * Adds violatoin-branching scores to all found exprs by means of addConsExprExprsViolScore().
 *
 * @note This method may modify the given auxvars array by means of sorting.
 */
static
SCIP_RETCODE addConsExprExprViolScoresAuxVars(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression where to start searching */
   SCIP_Real               violscore,        /**< violation score to add to expression */
   SCIP_VAR**              auxvars,          /**< auxiliary variables for which to find expression */
   int                     nauxvars,         /**< number of auxiliary variables */
   SCIP_SOL*               sol,              /**< current solution (NULL for the LP solution) */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_VAR* auxvar;
   SCIP_CONSEXPR_EXPR** exprs;
   int nexprs;
   int pos;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvars != NULL);
   assert(success != NULL);

   /* sort variables to make lookup below faster */
   SCIPsortPtr((void**)auxvars, SCIPvarComp, nauxvars);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_BFS, FALSE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nauxvars) );
   nexprs = 0;

   for( expr = SCIPexpriteratorGetNext(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) )  /*lint !e441*/
   {
      auxvar = SCIPgetConsExprExprAuxVar(expr);
      if( auxvar == NULL )
         continue;

      /* if auxvar of expr is contained in auxvars array, add branching score to expr */
      if( SCIPsortedvecFindPtr((void**)auxvars, SCIPvarComp, auxvar, nauxvars, &pos) )
      {
         assert(auxvars[pos] == auxvar);

         SCIPdebugMsg(scip, "adding branchingscore for expr %p with auxvar <%s>\n", expr, SCIPvarGetName(auxvar));
         exprs[nexprs++] = expr;

         if( nexprs == nauxvars )
            break;
      }
   }

   SCIPexpriteratorFree(&it);

   if( nexprs > 0 )
   {
      SCIP_CALL( SCIPaddConsExprExprsViolScore(scip, conshdlr, exprs, nexprs, violscore, sol, success) );
   }
   else
      *success = FALSE;

   SCIPfreeBufferArray(scip, &exprs);

   return SCIP_OKAY;
}

/** registers all unfixed variables in violated constraints as branching candidates */
static
SCIP_RETCODE registerBranchingCandidatesAllUnfixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nnotify != NULL);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      /* register all variables that have not been fixed yet */
      assert(consdata->varexprs != NULL);
      for( i = 0; i < consdata->nvarexprs; ++i )
      {
         var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
         assert(var != NULL);

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, getConsAbsViolation(conss[c]), SCIP_INVALID) );
            ++(*nnotify);
         }
      }
   }

   return SCIP_OKAY;
}

/** registers all variables in violated constraints with branching scores as external branching candidates */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< buffer to store whether at least one branching candidate was added */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it = NULL;
   int c;

   assert(conshdlr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPgetConsExprBranchAux(scip, conshdlr) )
   {
      SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   }

   /* register external branching candidates */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->varexprs != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      if( !SCIPgetConsExprBranchAux(scip, conshdlr) )
      {
         int i;

         /* if not branching on auxvars, then violation-branching scores will have been added to original variables
          * only, so we can loop over variable expressions
          */
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real violscore;
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_VAR* var;

            violscore = SCIPgetConsExprExprViolScore(conshdlr, consdata->varexprs[i]);

            /* skip variable expressions that do not have a violation score */
            if( violscore == 0.0 )
               continue;

            var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate "\
                        "with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )

               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }

            /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
             * several times as external branching candidate, see SCIPgetConsExprExprViolScore()
             */
            consdata->varexprs[i]->violscoretag = 0;
         }
      }
      else
      {
         SCIP_CONSEXPR_EXPR* expr;
         SCIP_VAR* var;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real violscore;

         for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
         {
            violscore = SCIPgetConsExprExprViolScore(conshdlr, expr);
            if( violscore == 0.0 )
               continue;

            /* if some nlhdlr added a branching score for this expression, then it considered this expression as a
             * variable, so this expression should either be an original variable or have an auxiliary variable
             */
            var = SCIPgetConsExprExprAuxVar(expr);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate "\
                        "with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )

               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }
         }
      }
   }

   if( SCIPgetConsExprBranchAux(scip, conshdlr) )
      SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** collect branching candidates from violated constraints
 *
 * Fills array with expressions that serve as branching candidates.
 * Collects those expressions that have a branching score assigned and stores the score in the auxviol field of the
 * branching candidate.
 *
 * If branching on aux-variables is allowed, then iterate through expressions of violated constraints, otherwise iterate
 * through variable-expressions only.
 */
static
SCIP_RETCODE collectBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   BRANCHCAND*           cands,              /**< array where to store candidates, must be at least SCIPgetNVars() long */
   int*                  ncands              /**< number of candidates found */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it = NULL;
   int c;
   int attempt;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetConsExprBranchAux(scip, conshdlr) )
   {
      SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   }

   *ncands = 0;
   for( attempt = 0; attempt < 2; ++attempt )
   {
      /* collect branching candidates from violated constraints
       * in the first attempt, consider only constraints with large violation
       * in the second attempt, consider all remaining violated constraints
       */
      for( c = 0; c < nconss; ++c )
      {
         SCIP_Real consviol;

         assert(conss != NULL && conss[c] != NULL);

         /* consider only violated constraints */
         if( !isConsViolated(scip, conss[c]) )
            continue;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         assert(consdata->varexprs != NULL);

         SCIP_CALL( getConsRelViolation(scip, conss[c], &consviol, sol, soltag) );

         if( attempt == 0 && consviol < conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;
         else if( attempt == 1 && consviol >= conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;

         if( !SCIPgetConsExprBranchAux(scip, conshdlr) )
         {
            int i;

            /* if not branching on auxvars, then violation-branching scores will be available for original variables
             * only, so we can loop over variable expressions
             * unfortunately, we don't know anymore which constraint contributed the violation-branching score to the
             * variable, therefore we invalidate the score of a variable after processing it.
             */
            for( i = 0; i < consdata->nvarexprs; ++i )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               /* skip variable expressions that do not have a valid violation score */
               if( conshdlrdata->enforound != consdata->varexprs[i]->violscoretag )
                  continue;

               var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = consdata->varexprs[i];
               cands[*ncands].auxviol = SCIPgetConsExprExprViolScore(conshdlr, consdata->varexprs[i]);
               ++(*ncands);

               /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
                * several times as external branching candidate */
               consdata->varexprs[i]->violscoretag = 0;
            }
         }
         else
         {
            SCIP_CONSEXPR_EXPR* expr;
            SCIP_Real lb;
            SCIP_Real ub;

            for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
            {
               if( expr->violscoretag != conshdlrdata->enforound )
                  continue;

               /* if some nlhdlr added a branching score for this expression, then it considered this expression as
                * variables, so this expression should either be an original variable or have an auxiliary variable
                */
               var = SCIPgetConsExprExprAuxVar(expr);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = expr;
               cands[*ncands].auxviol = SCIPgetConsExprExprViolScore(conshdlr, expr);
               ++(*ncands);
            }
         }
      }

      /* if we have branching candidates, then we don't need another attempt */
      if( *ncands > 0 )
         break;
   }

   if( SCIPgetConsExprBranchAux(scip, conshdlr) )
      SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** computes a branching score for a variable that reflects how important branching on this variable would be for
 * improving the dual bound from the LP relaxation
 *
 * Assume the Lagrangian for the current LP is something of the form
 *   L(x,z,lambda) = c'x + sum_i lambda_i (a_i'x - z_i + b_i) + ...
 * where x are the original variables, z the auxiliary variables,
 * and a_i'x - z_i + b_i <= 0 are the rows of the LP.
 *
 * Assume that a_i'x + b_i <= z_i was derived from some nonlinear constraint f(x) <= z and drop index i.
 * If we could have used not only an estimator, but the actual function f(x), then this would
 * have contributed lambda*(f(x) - z) to the Lagrangian function (though the value of z would be different).
 * Using a lot of handwaving, we claim that
 *   lambda_i * (f(x) - a_i'x + b_i)
 * is a value that can be used to quantity how much improving the estimator a'x + b <= z could change the dual bound.
 * If an estimator depended on local bounds, then it could be improved by branching.
 * We use row-is-local as proxy for estimator-depending-on-lower-bounds.
 *
 * To score a variable, we then sum the values lambda_i * (f(x) - a_i'x + b_i) for all rows in which the variable appears.
 * To scale, we divide by the LP objective value (if >1).
 *
 * TODO if we branch only on original variables, we neglect here estimators that are build on auxiliary variables
 *     these are affected by the bounds on original variables indirectly (through forward-propagation)
 * TODO if we branch also on auxiliary variables, then separating z from the x-variables in the row a'x+b <= z should happen
 *     in effect, we should go from the row to the expression for which it was generated and consider only variables that
 *     would also be branching candidates
 */
static
SCIP_Real getDualBranchscore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraints handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows;
   int nrows;
   int r;
   SCIP_Real dualscore;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(var != NULL);

   /* if LP not solved, then the dual branching score is not available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return 0.0;

   /* if var is not in the LP, then the dual branching score is not available */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0.0;

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
      return 0.0;

   nrows = SCIPcolGetNLPNonz(col);  /* TODO there is a big warning on when not to use this method; is the check for SCIPcolIsInLP sufficient? */
   rows = SCIPcolGetRows(col);

   /* SCIPinfoMessage(scip, enfologfile, " dualscoring <%s>\n", SCIPvarGetName(var)); */

   /* aggregate duals from all rows from consexpr with non-zero dual
    * TODO: this is a quick-and-dirty implementation, and not used by default
    *   in the long run, this should be either removed or replaced by a proper implementation
    */
   dualscore = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      SCIP_Real estimategap;
      const char* estimategapstr;

      /* rows from cuts that may be replaced by tighter ones after branching are the interesting ones
       * these would typically be local, unless they are created at the root node
       * so not check for local now, but trust that estimators that do not improve after branching will have an estimategap of 0
      if( !SCIProwIsLocal(rows[r]) )
         continue;
       */
      if( SCIProwGetOriginConshdlr(rows[r]) != conshdlr )
         continue;
      if( SCIPisZero(scip, SCIProwGetDualsol(rows[r])) )
         continue;

      estimategapstr = strstr(SCIProwGetName(rows[r]), "_estimategap=");
      if( estimategapstr == NULL ) /* gap not stored, maybe because it was 0 */
         continue;
      estimategap = atof(estimategapstr + 13);
      assert(estimategap >= 0.0);
      if( !SCIPisFinite(estimategap) || SCIPisHugeValue(scip, estimategap) )
         estimategap = SCIPgetHugeValue(scip);

      /* SCIPinfoMessage(scip, enfologfile, "  row <%s> contributes %g*|%g|: ", SCIProwGetName(rows[r]), estimategap, SCIProwGetDualsol(rows[r]));
      SCIP_CALL( SCIPprintRow(scip, rows[r], enfologfile) ); */

      dualscore += estimategap * REALABS(SCIProwGetDualsol(rows[r]));
   }

   /* divide by optimal value of LP for scaling */
   dualscore /= MAX(1.0, REALABS(SCIPgetLPObjval(scip)));  /*lint !e666*/

   return dualscore;
}

/** computes branching scores (including weighted score) for a set of candidates
 *
 * For each candidate in the array, compute and store the various branching scores (violation, pseudo-costs, vartype, domainwidth).
 * For pseudo-costs, it's possible that the score is not available, in which case cands[c].pscost will be set to SCIP_INVALID.
 *
 * For each score, compute the maximum over all candidates.
 *
 * Then compute for each candidate a "weighted" score using the weights as specified by parameters
 * and the scores as previously computed, but scale each score to be in [0,1], i.e., divide each score by the maximum
 * score all candidate.
 * Further divide by the sum of all weights where a score was available (even if the score was 0).
 *
 * For example:
 * - Let variable x have violation-score 10.0 and pseudo-cost-score 5.0.
 * - Let variable y have violation-score 12.0 but no pseudo-cost-score (because it hasn't yet been branched on sufficiently often).
 * - Assuming violation is weighted by 2.0 and pseudo-costs are weighted by 3.0.
 * - Then the weighted scores for x will be (2.0 * 10.0/12.0 + 3.0 * 5.0/5.0) / (2.0 + 3.0) = 0.9333.
 *   The weighted score for y will be (2.0 * 12.0/12.0) / 2.0 = 1.0.
 */
static
void scoreBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BRANCHCAND*           cands,              /**< branching candidates */
   int                   ncands,             /**< number of candidates */
   SCIP_SOL*             sol                 /**< solution to enforce (NULL for the LP solution) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND maxscore;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initialize counts to 0 */
   memset(&maxscore, 0, sizeof(BRANCHCAND));

   for( c = 0; c < ncands; ++c )
   {
      if( conshdlrdata->branchviolweight > 0.0 )
      {
         /* cands[c].auxviol was set in collectBranchingCandidates, so only update maxscore here */
         maxscore.auxviol = MAX(maxscore.auxviol, cands[c].auxviol);
      }

      if( conshdlrdata->branchdomainweight > 0.0 )
      {
         SCIP_Real domainwidth;
         SCIP_VAR* var;

         var = SCIPgetConsExprExprAuxVar(cands[c].expr);
         assert(var != NULL);

         /* get domain width, taking infinity at 1e20 on purpose */
         domainwidth = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

         /* domain-score is going to be log(2*infinity / domainwidth) if domain width >= 1
          * and log(2 * infinity *  MAX(epsilon, domainwidth)) for domain width < 1
          * the idea is to penalize very large and very small domains
          */
         if( domainwidth >= 1.0 )
            cands[c].domain = log10(2 * SCIPinfinity(scip) / domainwidth);
         else
            cands[c].domain = log10(2 * SCIPinfinity(scip) * MAX(SCIPepsilon(scip), domainwidth));  /*lint !e666*/

         maxscore.domain = MAX(cands[c].domain, maxscore.domain);
      }
      else
         cands[c].domain = 0.0;

      if( conshdlrdata->branchdualweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetConsExprExprAuxVar(cands[c].expr);
         assert(var != NULL);

         cands[c].dual = getDualBranchscore(scip, conshdlr, var);
         maxscore.dual = MAX(cands[c].dual, maxscore.dual);
      }

      if( conshdlrdata->branchpscostweight > 0.0 && SCIPgetNObjVars(scip) > 0 )
      {
         SCIP_VAR* var;

         var = SCIPgetConsExprExprAuxVar(cands[c].expr);
         assert(var != NULL);

         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            cands[c].pscost = SCIP_INVALID;
         else
         {
            SCIP_Real brpoint;
            SCIP_Real pscostdown;
            SCIP_Real pscostup;
            char strategy;

            /* decide how to compute pseudo-cost scores
             * this should be consistent with the way how pseudo-costs are updated in the core, which is decided by
             * branching/lpgainnormalize for continuous variables and move in LP-value for non-continuous variables
             */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               strategy = conshdlrdata->branchpscostupdatestrategy;
            else
               strategy = 'l';

            brpoint = SCIPgetBranchingPoint(scip, var, SCIP_INVALID);

            /* branch_relpscost deems pscosts as reliable, if the pseudo-count is at least something between 1 and 4
             * or it uses some statistical tests involving SCIPisVarPscostRelerrorReliable
             * For here, I use a simple #counts >= branchpscostreliable.
             * TODO use SCIPgetVarPseudocostCount() instead?
             */
            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_DOWNWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint)));
                     break;
                  case 'd' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var)));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, SCIPgetSolVal(scip, sol, var)) )
                        pscostdown = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, sol, var) <= SCIPadjustedVarUb(scip, var, brpoint) )
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPgetSolVal(scip, NULL, var) - SCIPadjustedVarUb(scip, var, brpoint)));
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostdown = SCIP_INVALID;
               }
            }
            else
               pscostdown = SCIP_INVALID;

            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_UPWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var));
                     break;
                  case 'd' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, var)) )
                        pscostup = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, NULL, var) >= SCIPadjustedVarLb(scip, var, brpoint) )
                        pscostup = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarLb(scip, var, brpoint) - SCIPgetSolVal(scip, NULL, var) );
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostup = SCIP_INVALID;
               }
            }
            else
               pscostup = SCIP_INVALID;

            if( pscostdown == SCIP_INVALID && pscostup == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = SCIP_INVALID;
            else if( pscostdown == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = pscostup;
            else if( pscostup == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = pscostdown;
            else
               cands[c].pscost = SCIPgetBranchScore(scip, NULL, pscostdown, pscostup);  /* pass NULL for var to avoid multiplication with branch-factor */
         }

         if( cands[c].pscost != SCIP_INVALID )  /*lint !e777*/
            maxscore.pscost = MAX(cands[c].pscost, maxscore.pscost);
      }

      if( conshdlrdata->branchvartypeweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetConsExprExprAuxVar(cands[c].expr);
         assert(var != NULL);

         switch( SCIPvarGetType(var) )
         {
            case SCIP_VARTYPE_BINARY :
               cands[c].vartype = 1.0;
               break;
            case SCIP_VARTYPE_INTEGER :
               cands[c].vartype = 0.1;
               break;
            case SCIP_VARTYPE_IMPLINT :
               cands[c].vartype = 0.01;
               break;
            case SCIP_VARTYPE_CONTINUOUS :
            default:
               cands[c].vartype = 0.0;
         }
         maxscore.vartype = MAX(cands[c].vartype, maxscore.vartype);
      }
   }

   /* now computed a weighted score for each candidate from the single scores
    * the single scores are scaled to be in [0,1] for this
    */
   for( c = 0; c < ncands; ++c )
   {
      SCIP_Real weightsum;

      ENFOLOG(
         SCIP_VAR* var;
         var = SCIPgetConsExprExprAuxVar(cands[c].expr);
         SCIPinfoMessage(scip, enfologfile, " scoring <%8s>[%7.1g,%7.1g]:(", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         )

      cands[c].weighted = 0.0;
      weightsum = 0.0;

      if( maxscore.auxviol > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchviolweight * cands[c].auxviol / maxscore.auxviol;
         weightsum += conshdlrdata->branchviolweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(viol)", conshdlrdata->branchviolweight, cands[c].auxviol / maxscore.auxviol); )
      }

      if( maxscore.domain > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdomainweight * cands[c].domain / maxscore.domain;
         weightsum += conshdlrdata->branchdomainweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(domain)", conshdlrdata->branchdomainweight, cands[c].domain / maxscore.domain); )
      }

      if( maxscore.dual > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdualweight * cands[c].dual / maxscore.dual;
         weightsum += conshdlrdata->branchdualweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(dual)", conshdlrdata->branchdualweight, cands[c].dual / maxscore.dual); )
      }

      /* use pseudo-costs, if we have some for at least half the candidates */
      if( maxscore.pscost > 0.0 )
      {
         if( cands[c].pscost != SCIP_INVALID )  /*lint !e777*/
         {
            cands[c].weighted += conshdlrdata->branchpscostweight * cands[c].pscost / maxscore.pscost;
            weightsum += conshdlrdata->branchpscostweight;

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(pscost)", conshdlrdata->branchpscostweight, cands[c].pscost / maxscore.pscost); )
         }
         else
         {
            /* do not add pscostscore, if not available, also do not add into weightsum */
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " +0.0*    n/a(pscost)"); )
         }
      }

      if( maxscore.vartype > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchvartypeweight * cands[c].vartype / maxscore.vartype;
         weightsum += conshdlrdata->branchvartypeweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%6.2g(vartype)", conshdlrdata->branchvartypeweight, cands[c].vartype / maxscore.vartype); )
      }
      assert(weightsum > 0.0);  /* we should have got at least one valid score */
      cands[c].weighted /= weightsum;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " ) / %g = %g\n", weightsum, cands[c].weighted); )
   }
}

/** compare two branching candidates by their weighted score
 *
 * if weighted score is equal, use variable index of (aux)var
 */
static
SCIP_DECL_SORTINDCOMP(branchcandCompare)
{
   BRANCHCAND* cands = (BRANCHCAND*)dataptr;

   if( cands[ind1].weighted != cands[ind2].weighted )  /*lint !e777*/
      return cands[ind1].weighted < cands[ind2].weighted ? -1 : 1;
   else
      return SCIPvarGetIndex(SCIPgetConsExprExprAuxVar(cands[ind1].expr)) - SCIPvarGetIndex(SCIPgetConsExprExprAuxVar(cands[ind2].expr));
}

/** do branching or register branching candidates */
static
SCIP_RETCODE branching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_RESULT*          result              /**< pointer to store the result of branching */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND* cands;
   int ncands;
   SCIP_VAR* var;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->branchexternal )
   {
      /* just register branching candidates as external */
      SCIP_Bool success;

      SCIP_CALL( registerBranchingCandidates(scip, conshdlr, conss, nconss, &success) );
      if( success )
         *result = SCIP_INFEASIBLE;

      return SCIP_OKAY;
   }

   /* collect branching candidates and their auxviol-score */
   SCIP_CALL( SCIPallocBufferArray(scip, &cands, SCIPgetNVars(scip)) );
   SCIP_CALL( collectBranchingCandidates(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, cands, &ncands) );

   /* if no unfixed branching candidate in all violated constraint, then it's probably numerics that prevented us to separate or decide a cutoff
    * we will return here and let the fallbacks in consEnfo() decide how to proceed
    */
   if( ncands == 0 )
      goto TERMINATE;

   if( ncands > 1 )
   {
      /* if there are more than one candidate, then compute scores and select */
      int* perm;
      int c;
      int left;
      int right;
      SCIP_Real threshold;

      /* compute additional scores on branching candidates and weighted score */
      scoreBranchingCandidates(scip, conshdlr, cands, ncands, sol);

      /* sort candidates by weighted score */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, ncands) );
      SCIPsortDown(perm, branchcandCompare, (void*)cands, ncands);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g)\n", ncands,
         SCIPvarGetName(SCIPgetConsExprExprAuxVar(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetConsExprExprAuxVar(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      /* binary search to find first low-scored (score below branchhighscorefactor * maximal-score)  candidate */
      left = 0;
      right = ncands - 1;
      threshold = conshdlrdata->branchhighscorefactor * cands[perm[0]].weighted;
      while( left < right )
      {
         int mid = (left + right) / 2;
         if( cands[perm[mid]].weighted >= threshold )
            left = mid + 1;
         else
            right = mid;
      }
      assert(left <= ncands);

      if( left < ncands )
      {
         if( cands[perm[left]].weighted >= threshold )
         {
            assert(left + 1 == ncands || cands[perm[left + 1]].weighted < threshold);
            ncands = left + 1;
         }
         else
         {
            assert(cands[perm[left]].weighted < threshold);
            ncands = left;
         }
      }
      assert(ncands > 0);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g) after removing low scores\n", ncands,
         SCIPvarGetName(SCIPgetConsExprExprAuxVar(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetConsExprExprAuxVar(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      if( ncands > 1 )
      {
         /* choose at random from candidates 0..ncands-1 */
         if( conshdlrdata->branchrandnumgen == NULL )
         {
            SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->branchrandnumgen, BRANCH_RANDNUMINITSEED, TRUE) );
         }
         c = SCIPrandomGetInt(conshdlrdata->branchrandnumgen, 0, ncands - 1);
         var = SCIPgetConsExprExprAuxVar(cands[perm[c]].expr);
      }
      else
         var = SCIPgetConsExprExprAuxVar(cands[perm[0]].expr);

      SCIPfreeBufferArray(scip, &perm);
   }
   else
   {
      var = SCIPgetConsExprExprAuxVar(cands[0].expr);
   }
   assert(var != NULL);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " branching on variable <%s>[%g,%g]\n", SCIPvarGetName(var),
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)); )

   SCIP_CALL( SCIPbranchVarVal(scip, var, SCIPgetBranchingPoint(scip, var, SCIP_INVALID), &downchild, &eqchild,
            &upchild) );
   if( downchild != NULL || eqchild != NULL || upchild != NULL )
      *result = SCIP_BRANCHED;
   else
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      *result = SCIP_REDUCEDDOM;

 TERMINATE:
   SCIPfreeBufferArray(scip, &cands);

   return SCIP_OKAY;
}

/** call enforcement or estimate callback of nonlinear handler
 *
 * Calls the enforcement callback, if available.
 * Otherwise, calls the estimate callback, if available, and constructs a cut from the estimator.
 *
 * If cut is weak, but estimator is not tight, tries to add branching candidates.
 */
static
SCIP_RETCODE enforceExprNlhdlr(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_CONSEXPR_NLHDLR* nlhdlr,             /**< nonlinear handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler data of expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Bool             separated,          /**< whether another nonlinear handler already added a cut for this expression */
   SCIP_Bool             allowweakcuts,      /**< whether we allow for weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real mincutviolation;

   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* decide on minimal violation of cut */
   if( sol == NULL )
      mincutviolation = SCIPgetLPFeastol(scip);  /* we enforce an LP solution */
   else
      mincutviolation = SCIPfeastol(scip);

   /* call enforcement callback of the nlhdlr */
   SCIP_CALL( SCIPenfoConsExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
            allowweakcuts, separated, inenforcement, result) );

   /* if it was not running (e.g., because it was not available) or did not find anything, then try with estimator callback */
   if( *result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    sepa of nlhdlr %s succeeded with result %d\n",
               SCIPgetConsExprNlhdlrName(nlhdlr), *result); )
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* now call the estimator callback of the nlhdlr */
   if( SCIPhasConsExprNlhdlrEstimate(nlhdlr) )
   {
      SCIP_ROWPREP* rowprep;
      SCIP_VAR* auxvar;
      SCIP_Real auxvarvalue;
      SCIP_Real cutviol;
      SCIP_Real estimateval = SCIP_INVALID;
      SCIP_Bool sepasuccess = FALSE;
      SCIP_Bool branchscoresuccess = FALSE;

      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

      auxvar = SCIPgetConsExprExprAuxVar(expr);
      assert(auxvar != NULL);

      SCIP_CALL( SCIPestimateConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
               SCIPgetSolVal(scip, sol, auxvar), rowprep, &sepasuccess, inenforcement, &branchscoresuccess) );

      if( sepasuccess )
      {
         /* complete estimator to cut */
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

         cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, NULL);
         if( cutviol > 0.0 )
         {
            auxvarvalue = SCIPgetSolVal(scip, sol, auxvar);

            /* check whether cut is weak (if f(x) not defined, then it's never weak) */
            if( !allowweakcuts && auxvalue != SCIP_INVALID )  /*lint !e777*/
            {
               /* let the estimator be c'x-b, the auxvar is z (=auxvarvalue), and the expression is f(x) (=auxvalue)
                * then if we are underestimating and since the cut is violated, we should have z <= c'x-b <= f(x)
                * cutviol is c'x-b - z, so estimator value is c'x-b = z + cutviol
                * if the estimator value (c'x-b) is too close to z (auxvarvalue), when compared to f(x) (auxvalue), then
                * let's call this a weak cut that is, it's a weak cut if c'x-b <= z + weakcutthreshold * (f(x)-z)
                *   <->   c'x-b - z <= weakcutthreshold * (f(x)-z)
                *
                * if we are overestimating, we have z >= c'x-b >= f(x)
                * cutviol is z - (c'x-b), so estimator value is c'x-b = z - cutviol
                * it's weak if c'x-b >= f(x) + (1-weakcutthreshold) * (z - f(x))
                *   <->   c'x-b - z >= weakcutthreshold * (f(x)-z)
                *
                * when linearizing convex expressions, then we should have c'x-b = f(x), so they would never be weak
                */
               if( (!overestimate && ( cutviol <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
                   ( overestimate && (-cutviol >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut is too "\
                           "weak: auxvarvalue %g estimateval %g auxvalue %g (over %d)\n",
                           SCIPgetConsExprNlhdlrName(nlhdlr), auxvarvalue,
                           auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate); )
                  sepasuccess = FALSE;
               }
            }

            /* save estimator value for later, see long comment above why this gives the value for c'x-b */
            estimateval = auxvarvalue + (!overestimate ? cutviol : -cutviol);
         }
         else
         {
            sepasuccess = FALSE;
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut does not "\
                     "separate\n", SCIPgetConsExprNlhdlrName(nlhdlr)); )
         }
      }
      else
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s failed\n",
                  SCIPgetConsExprNlhdlrName(nlhdlr)); )
      }

      /* clean up estimator */
      if( sepasuccess )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded: auxvarvalue %g estimateval "\
                  "%g auxvalue %g (over %d)\n    ", SCIPgetConsExprNlhdlrName(nlhdlr), auxvarvalue,
                  auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate);
                  SCIPprintRowprep(scip, rowprep, enfologfile); )

         /* if not allowweakcuts, then do not attempt to get cuts more violated by scaling them up,
          * instead, may even scale them down, that is, scale so that max coef is close to 1
          */
         if( !allowweakcuts )
         {
            SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, conshdlrdata->strongcutmaxcoef, &sepasuccess) );

            if( !sepasuccess )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup cut failed due to bad numerics\n"); )
            }
            else
            {
               cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, &sepasuccess);
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup succeeded, violation = %g and %sreliable, min "\
                        "requ viol = %g\n", cutviol, sepasuccess ? "" : "not ", mincutviolation); )
               if( sepasuccess )
                  sepasuccess = cutviol > mincutviolation;
            }

            if( sepasuccess && auxvalue != SCIP_INVALID ) /*lint !e777*/
            {
               /* check whether cut is weak now
                * auxvar z may now have a coefficient due to scaling (down) in cleanup - take this into account when
                * reconstructing estimateval from cutviol (TODO improve or remove?)
                */
               SCIP_Real auxvarcoef = 0.0;
               int i;

               /* get absolute value of coef of auxvar in row; this makes the check more expensive than it should be */
               for( i = 0; i < rowprep->nvars; ++i )
               {
                  if( rowprep->vars[i] == auxvar )
                  {
                     auxvarcoef = REALABS(rowprep->coefs[i]);
                     break;
                  }
               }

               if( auxvarcoef == 0.0 ||
                   (!overestimate && ( cutviol / auxvarcoef <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
                   ( overestimate && (-cutviol / auxvarcoef >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )  /*lint !e644*/
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut is too weak after cleanup: auxvarvalue %g" \
                           "estimateval %g auxvalue %g (over %d)\n", auxvalue, auxvarvalue,
                           auxvarvalue + (overestimate ? -cutviol : cutviol) / auxvarcoef, auxvalue, overestimate); )
                  sepasuccess = FALSE;
               }
            }
         }
         else
         {
            /* TODO if violations are really tiny, then maybe handle special (decrease LP feastol, for example) */

            /* if estimate didn't report branchscores explicitly, then consider branching on those children for which
             * the following cleanup changes coefficients (we had/have this in cons_expr_sum this way)
             */
            if( !branchscoresuccess )
               rowprep->recordmodifications = TRUE;

            SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, mincutviolation, &cutviol, &sepasuccess) );

            if( !sepasuccess )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup failed, %d coefs modified, cutviol %g\n",
                        rowprep->nmodifiedvars, cutviol); )
            }

            /* if cleanup left us with a useless cut, then consider branching on variables for which coef were changed */
            if( !sepasuccess && !branchscoresuccess && rowprep->nmodifiedvars > 0 )
            {
               SCIP_Real violscore;

#ifdef BRSCORE_ABSVIOL
               violscore = getExprAbsAuxViolation(scip, expr, auxvalue, sol, NULL, NULL);
#else
               SCIP_CALL( SCIPgetConsExprExprRelAuxViolation(scip, conshdlr, expr, auxvalue, sol, &violscore, NULL,
                        NULL) );
#endif
               SCIP_CALL( addConsExprExprViolScoresAuxVars(scip, conshdlr, expr, violscore, rowprep->modifiedvars,
                        rowprep->nmodifiedvars, sol, &branchscoresuccess) );

               /* addConsExprExprBranchScoresAuxVars can fail if the only var for which the coef was changed is this
                * expr's auxvar
                * I don't think it makes sense to branch on that one (would it?)
                * it can also fail if everything is fixed
                */
               assert(branchscoresuccess || (rowprep->nmodifiedvars == 1 && rowprep->modifiedvars[0] == auxvar));
            }
         }
      }

      /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
      if( sepasuccess )
      {
         SCIP_ROW* row;

         if( conshdlrdata->branchdualweight > 0.0 )
         {
            /* store remaining gap |f(x)-estimateval| in row name, which could be used in getDualBranchscore
             * skip if gap is zero
             */
            if( auxvalue == SCIP_INVALID )  /*lint !e777*/
               strcat(rowprep->name, "_estimategap=inf");
            else if( !SCIPisEQ(scip, auxvalue, estimateval) )
            {
               char gap[40];
               sprintf(gap, "_estimategap=%g", REALABS(auxvalue - estimateval));
               strcat(rowprep->name, gap);
            }
         }

         SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

         if( !allowweakcuts && conshdlrdata->strongcutefficacy && !SCIPisCutEfficacious(scip, sol, row) )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut efficacy %g is too low (minefficacy=%g)\n",
               SCIPgetCutEfficacy(scip, sol, row), SCIPgetSepaMinEfficacy(scip)); )
         }
         else
         {
            SCIP_Bool infeasible;

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    adding cut ");
            SCIP_CALL( SCIPprintRow(scip, row, enfologfile) ); )

            /* I take !allowweakcuts as equivalent for having a strong cut (we usually have allowweakcuts=TRUE only if we haven't found strong cuts before)
             */
            SCIP_CALL( SCIPaddRow(scip, row, conshdlrdata->forcestrongcut && !allowweakcuts && inenforcement, &infeasible) );

            if( infeasible )
            {
               *result = SCIP_CUTOFF;
               ++nlhdlr->ncutoffs;
            }
            else
            {
               *result = SCIP_SEPARATED;
               ++nlhdlr->nseparated;
            }
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
      else if( branchscoresuccess )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed, but branching" \
                  "candidates added\n", SCIPgetConsExprNlhdlrName(nlhdlr)); )

         /* well, not branched, but addConsExprExprViolScoresAuxVars() added scores to (aux)variables and that makes the
          * expressions eligible for branching candidate, see enforceConstraints() and branching()
          */
         *result = SCIP_BRANCHED;
      }
      else
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed and no " \
                  "branching candidates%s\n", SCIPgetConsExprNlhdlrName(nlhdlr), (allowweakcuts && inenforcement) ?
                  " (!)" : ""); )
      }

      SCIPfreeRowprep(scip, &rowprep);
   }

   return SCIP_OKAY;
}

/** tries to enforce violation in an expression by separation, bound tightening, or finding a branching candidate
 *
 * if not inenforcement, then we should be called by consSepa, and thus only try separation
 */
static
SCIP_RETCODE enforceExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraints handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Bool             allowweakcuts,      /**< whether we allow weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_Real origviol;
   SCIP_Bool underestimate;
   SCIP_Bool overestimate;
   SCIP_Real auxviol;
   SCIP_Bool auxunderestimate;
   SCIP_Bool auxoverestimate;
   SCIP_RESULT hdlrresult;
   int e;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr->auxvar != NULL);  /* there must be a variable attached to the expression in order to construct a cut here */
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* make sure that this expression has been evaluated */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol, soltag) );

   /* decide whether under- or overestimate is required and get amount of violation */
   origviol = getExprAbsOrigViolation(scip, expr, sol, &underestimate, &overestimate);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* no sufficient violation w.r.t. the original variables -> skip expression */
   if( !overestimate && !underestimate )
      return SCIP_OKAY;

   /* check aux-violation w.r.t. each nonlinear handlers and try to enforce when there is a decent violation */
   for( e = 0; e < expr->nenfos; ++e )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr;

      nlhdlr = expr->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
      SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );
      ENFOLOG(
         SCIPinfoMessage(scip, enfologfile, "  expr ");
         SCIPprintConsExprExpr(scip, conshdlr, expr, enfologfile);
         SCIPinfoMessage(scip, enfologfile, " (%p): evalvalue %.15g auxvarvalue %.15g [%.15g,%.15g], nlhdlr <%s>" \
            "auxvalue: %.15g\n", (void*)expr, expr->evalvalue, SCIPgetSolVal(scip, sol, expr->auxvar),
            expr->activity.inf, expr->activity.sup, nlhdlr->name, expr->enfos[e]->auxvalue);
      )

      /* TODO if expr is root of constraint (consdata->expr == expr),
       * then compare auxvalue with constraint sides instead of auxvarvalue, as the former is what actually matters
       * that is, if auxvalue is good enough for the constraint to be satisfied, but when looking at evalvalue we see
       * the the constraint is violated, then some of the auxvars that nlhdlr uses is not having a good enough value,
       * so we should enforce in these auxiliaries first
       * if changing this here, we must also adapt analyzeViolation
       */

      auxviol = getExprAbsAuxViolation(scip, expr, expr->enfos[e]->auxvalue, sol, &auxunderestimate, &auxoverestimate);
      assert(auxviol >= 0.0);

      /* if aux-violation is much smaller than orig-violation, then better enforce further down in the expression first */
      if( !SCIPisInfinity(scip, auxviol) && auxviol < conshdlrdata->enfoauxviolfactor * origviol )  /*lint !e777*/
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with" \
                  "auxviolation %g << origviolation %g under:%d over:%d\n", nlhdlr->name, (void*)expr,
                  expr->exprhdlr->name, auxviol, origviol, underestimate, overestimate); )

         /* TODO expr->lastenforced = conshdlrdata->enforound;  ??? */
         continue;
      }

      /* if aux-violation is small (below feastol) and we look only for strong cuts, then it's unlikely to give a strong cut, so skip it */
      if( !allowweakcuts && auxviol < SCIPfeastol(scip) )  /*lint !e777*/
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with tiny " \
                  "auxviolation %g under:%d over:%d\n", nlhdlr->name, (void*)expr, expr->exprhdlr->name, auxviol,
                  underestimate, overestimate); )

         /* TODO expr->lastenforced = conshdlrdata->enforound;  ??? */
         continue;
      }

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   enforce using nlhdlr <%s> for expr %p (%s) with auxviolation " \
               "%g origviolation %g under:%d over:%d weak:%d\n", nlhdlr->name, (void*)expr, expr->exprhdlr->name,
               auxviol, origviol, underestimate, overestimate, allowweakcuts); )

      /* if we want overestimation and violation w.r.t. auxiliary variables is also present on this side, then call separation of nlhdlr */
      if( overestimate && auxoverestimate )  /*lint !e777*/
      {
         /* call the separation or estimation callback of the nonlinear handler for overestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, sol,
            expr->enfos[e]->auxvalue, TRUE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            expr->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", nlhdlr->name); )
            *result = SCIP_SEPARATED;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", nlhdlr->name); )
            *result = SCIP_REDUCEDDOM;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", nlhdlr->name); )
            assert(inenforcement);

            /* separation and domain reduction takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            expr->lastenforced = conshdlrdata->enforound;
         }
      }

      /* if we want underestimation and violation w.r.t. auxiliary variables is also present on this side, then call separation of nlhdlr */
      if( underestimate && auxunderestimate )  /*lint !e777*/
      {
         /* call the separation or estimation callback of the nonlinear handler for underestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, sol,
            expr->enfos[e]->auxvalue, FALSE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            expr->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", nlhdlr->name); )
            *result = SCIP_SEPARATED;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", nlhdlr->name); )
            *result = SCIP_REDUCEDDOM;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", nlhdlr->name); )
            assert(inenforcement);

            /* separation takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            expr->lastenforced = conshdlrdata->enforound;
         }
      }
   }

   return SCIP_OKAY;
}

/** helper function to enforce a single constraint */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_CONSEXPR_ITERATOR* it,               /**< expression iterator that we can just use here */
   SCIP_Bool             allowweakcuts,      /**< whether to allow weak cuts in this round */
   SCIP_Bool             inenforcement,      /**< whether to we are in enforcement, and not just separation */
   SCIP_RESULT*          result,             /**< pointer to update with result of the enforcing call */
   SCIP_Bool*            success             /**< buffer to store whether some enforcement took place */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPR* expr;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(it != NULL);
   assert(result != NULL);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *success = FALSE;

   for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      SCIP_RESULT resultexpr;

      /* we can only enforce if there is an auxvar to compare with */
      if( expr->auxvar == NULL )
         continue;

      assert(expr->lastenforced <= conshdlrdata->enforound);
      if( expr->lastenforced == conshdlrdata->enforound )
      {
         ENFOLOG(
            SCIPinfoMessage(scip, enfologfile, "  skip expr ");
            SCIPprintConsExprExpr(scip, conshdlr, expr, enfologfile);
            SCIPinfoMessage(scip, enfologfile, " as already enforced in this enforound\n");
         )
         *success = TRUE;
         continue;
      }

      SCIP_CALL( enforceExpr(scip, conshdlr, cons, expr, sol, soltag, allowweakcuts, inenforcement, &resultexpr) );

      /* if not enforced, then we must not have found a cutoff, cut, domain reduction, or branchscore */
      assert((expr->lastenforced == conshdlrdata->enforound) == (resultexpr != SCIP_DIDNOTFIND));
      if( expr->lastenforced == conshdlrdata->enforound )
         *success = TRUE;

      if( resultexpr == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if( resultexpr == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;

      if( resultexpr == SCIP_REDUCEDDOM && *result != SCIP_SEPARATED )
         *result = SCIP_REDUCEDDOM;

      if( resultexpr == SCIP_BRANCHED && *result != SCIP_SEPARATED && *result != SCIP_REDUCEDDOM )
         *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** try to separate violated constraints and, if in enforcement, register branching scores
 *
 * Sets result to
 * - SCIP_DIDNOTFIND, if nothing of the below has been done
 * - SCIP_CUTOFF, if node can be cutoff,
 * - SCIP_SEPARATED, if a cut has been added,
 * - SCIP_REDUCEDDOM, if a domain reduction has been found,
 * - SCIP_BRANCHED, if branching has been done,
 * - SCIP_REDUCEDDOM, if a variable got fixed (in an attempt to branch on it),
 * - SCIP_INFEASIBLE, if external branching candidates were registered
 */
static
SCIP_RETCODE enforceConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, and not just separation */
   SCIP_Real             maxrelconsviol,     /**< largest scaled violation among all violated expr-constraints, only used if in enforcement */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_Bool consenforced;  /* whether any expression in constraint could be enforced */
   int c;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase tag to tell whether branching scores in expression belong to this sweep
    * and which expressions have already been enforced in this sweep
    * (we also want to distinguish sepa rounds, so this need to be here and not in consEnfo)
    */
   ++(conshdlrdata->enforound);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, TRUE) );

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      /* skip constraints that are not enabled or deleted */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      /* skip constraints that have separation disabled if we are only in separation */
      if( !inenforcement && !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      ENFOLOG(
      {
         SCIP_CONSDATA* consdata;
         int i;
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         SCIPinfoMessage(scip, enfologfile, " constraint ");
         SCIP_CALL( SCIPprintCons(scip, conss[c], enfologfile) );
         SCIPinfoMessage(scip, enfologfile, "\n with viol %g and point\n", getConsAbsViolation(conss[c]));
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_VAR* var;
            var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
            SCIPinfoMessage(scip, enfologfile, "  %-10s = %15g bounds: [%15g,%15g]\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      })

      SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, FALSE, inenforcement, result, &consenforced) );

      if( *result == SCIP_CUTOFF )
         break;

      if( !consenforced && inenforcement )
      {
         SCIP_Real viol;

         SCIP_CALL( getConsRelViolation(scip, conss[c], &viol, sol, soltag) );
         if( viol > conshdlrdata->weakcutminviolfactor * maxrelconsviol )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " constraint <%s> could not be enforced, try again with weak "\
                     "cuts allowed\n", SCIPconsGetName(conss[c])); )

            SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, TRUE, inenforcement, result, &consenforced) );

            if( consenforced )
               ++conshdlrdata->nweaksepa;  /* TODO maybe this should not be counted per constraint, but per enforcement round? */

            if( *result == SCIP_CUTOFF )
               break;
         }
      }
   }

   SCIPexpriteratorFree(&it);

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   /* if having branching scores, then propagate them from expressions with children to variable expressions */
   if( *result == SCIP_BRANCHED )
   {
      /* having result set to branched here means only that we have branching candidates, we still need to do the actual
       * branching
       */
      SCIP_CALL( branching(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, result) );

      /* branching should either have branched: result == SCIP_BRANCHED,
       * or fixed a variable: result == SCIP_REDUCEDDOM,
       * or have registered external branching candidates: result == SCIP_INFEASIBLE,
       * or have not done anything: result == SCIP_DIDNOTFIND
       */
      assert(*result == SCIP_BRANCHED || *result == SCIP_REDUCEDDOM || *result == SCIP_INFEASIBLE || *result == SCIP_DIDNOTFIND);
   }

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   return SCIP_OKAY;
}

/** collect (and print (if debugging enfo)) information on violation in expressions
 *
 * assumes that constraint violations have been computed
 */
static
SCIP_RETCODE analyzeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Real*            maxabsconsviol,     /**< buffer to store maximal absolute violation of constraints */
   SCIP_Real*            maxrelconsviol,     /**< buffer to store maximal relative violation of constraints */
   SCIP_Real*            minauxviol,         /**< buffer to store minimal (nonzero) violation of auxiliaries */
   SCIP_Real*            maxauxviol,         /**< buffer to store maximal violation of auxiliaries (violation in "extended formulation") */
   SCIP_Real*            maxvarboundviol     /**< buffer to store maximal violation of variable bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real v;
   int c;

   assert(conss != NULL || nconss == 0);
   assert(maxabsconsviol != NULL);
   assert(maxrelconsviol != NULL);
   assert(maxauxviol != NULL);
   assert(maxvarboundviol != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   *maxabsconsviol = 0.0;
   *maxrelconsviol = 0.0;
   *minauxviol = SCIPinfinity(scip);
   *maxauxviol = 0.0;
   *maxvarboundviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      v = getConsAbsViolation(conss[c]);
      *maxabsconsviol = MAX(*maxabsconsviol, v);

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      SCIP_CALL( getConsRelViolation(scip, conss[c], &v, sol, soltag) );
      *maxrelconsviol = MAX(*maxrelconsviol, v);

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIP_Real auxvarvalue;
         SCIP_Real auxvarlb;
         SCIP_Real auxvarub;
         SCIP_Bool violunder;
         SCIP_Bool violover;
         SCIP_Real origviol;
         SCIP_Real auxviol;
         int e;

         if( expr->auxvar == NULL )
         {
            /* check violation of variable bounds of original variable */
            if( SCIPisConsExprExprVar(expr) )
            {
               SCIP_VAR* var;
               var = SCIPgetConsExprExprVarVar(expr);
               auxvarvalue = SCIPgetSolVal(scip, sol, var);
               auxvarlb = SCIPvarGetLbLocal(var);
               auxvarub = SCIPvarGetUbLocal(var);

               origviol = 0.0;
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  origviol = auxvarlb - auxvarvalue;
               else if( auxvarub < auxvarvalue && !SCIPisInfinity(scip, auxvarub) )
                  origviol = auxvarvalue - auxvarub;
               if( origviol <= 0.0 )
                  continue;

               *maxvarboundviol = MAX(*maxvarboundviol, origviol);

               ENFOLOG(
               SCIPinfoMessage(scip, enfologfile, "var <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(var), auxvarlb, auxvarub, auxvarvalue);
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  SCIPinfoMessage(scip, enfologfile, " var >= lb violated by %g", auxvarlb - auxvarvalue);
               if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
                  SCIPinfoMessage(scip, enfologfile, " var <= ub violated by %g", auxvarvalue - auxvarub);
               SCIPinfoMessage(scip, enfologfile, "\n");
               )
            }

            continue;
         }

         auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);
         auxvarlb = SCIPvarGetLbLocal(expr->auxvar);
         auxvarub = SCIPvarGetUbLocal(expr->auxvar);

         /* check violation of variable bounds of auxiliary variable */
         if( auxvarlb - auxvarvalue > *maxvarboundviol && !SCIPisInfinity(scip, -auxvarlb) )
            *maxvarboundviol = auxvarlb - auxvarvalue;
         else if( auxvarvalue - auxvarub > *maxvarboundviol && !SCIPisInfinity(scip,  auxvarub) )
            *maxvarboundviol = auxvarvalue - auxvarub;

         origviol = getExprAbsOrigViolation(scip, expr, sol, &violunder, &violover);

         ENFOLOG(
         if( origviol > 0.0 || auxvarlb > auxvarvalue || auxvarub < auxvarvalue )
         {
            SCIPinfoMessage(scip, enfologfile, "expr ");
            SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, enfologfile) );
            SCIPinfoMessage(scip, enfologfile, " (%p)[%.15g,%.15g] = %.15g\n", (void*)expr, expr->activity.inf, expr->activity.sup, expr->evalvalue);

            SCIPinfoMessage(scip, enfologfile, "  auxvar <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(expr->auxvar), auxvarlb, auxvarub, auxvarvalue);
            if( origviol > 0.0 )
               SCIPinfoMessage(scip, enfologfile, " auxvar %s expr violated by %g", violunder ? ">=" : "<=", origviol);
            if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
               SCIPinfoMessage(scip, enfologfile, " auxvar >= auxvar's lb violated by %g", auxvarlb - auxvarvalue);
            if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
               SCIPinfoMessage(scip, enfologfile, " auxvar <= auxvar's ub violated by %g", auxvarvalue - auxvarub);
            SCIPinfoMessage(scip, enfologfile, "\n");
         }
         )

         /* no violation w.r.t. the original variables -> skip expression */
         if( origviol == 0.0 )
            continue;

         /* TODO remove? origviol shouldn't be mixed up with auxviol */
         *maxauxviol = MAX(*maxauxviol, origviol);  /*lint !e666*/
         *minauxviol = MIN(*minauxviol, origviol);  /*lint !e666*/

         /* compute aux-violation for each nonlinear handlers */
         for( e = 0; e < expr->nenfos; ++e )
         {
            SCIP_CONSEXPR_NLHDLR* nlhdlr;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
            SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "  nlhdlr <%s> = %.15g", nlhdlr->name, expr->enfos[e]->auxvalue); )

            auxviol = getExprAbsAuxViolation(scip, expr, expr->enfos[e]->auxvalue, sol, &violunder, &violover);

            if( auxviol > 0.0 )  /*lint !e777*/
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " auxvar %s nlhdlr-expr violated by %g", violover ? "<=" : ">=", auxviol); )
               *maxauxviol = MAX(*maxauxviol, auxviol);
               *minauxviol = MIN(*minauxviol, auxviol);
            }
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "\n"); )
         }
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
} /*lint !e715*/

/** enforcement of constraints called by enfolp and enforelax */
static
SCIP_RETCODE consEnfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real maxabsconsviol;
   SCIP_Real maxrelconsviol;
   SCIP_Real minauxviol;
   SCIP_Real maxauxviol;
   SCIP_Real maxvarboundviol;
   unsigned int soltag;
   int nnotify;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   soltag = ++conshdlrdata->lastsoltag;

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: all expr-constraints feasible, skip enforcing\n",
               SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   SCIP_CALL( analyzeViolation(scip, conshdlr, conss, nconss, sol, soltag, &maxabsconsviol, &maxrelconsviol,
            &minauxviol, &maxauxviol, &maxvarboundviol) );

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: enforcing constraints with max conssviol=%e (rel=%e), "\
            "auxviolations in %g..%g, variable bounds violated by at most %g\n",
            SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), maxabsconsviol, maxrelconsviol, minauxviol, maxauxviol,
            maxvarboundviol); )

   assert(maxvarboundviol <= SCIPgetLPFeastol(scip));

   /* try to propagate */
   if( conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;
      int ndelconss = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds, &ndelconss) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* tighten the LP tolerance if violation in variables bounds is larger than aux-violation (max |expr - auxvar| over
    * all violated expr/auxvar in violated constraints)
    */
   if( conshdlrdata->tightenlpfeastol && maxvarboundviol > maxauxviol && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) &&
         sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));  /*lint !e666*/
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bound violation %g larger than auxiliary violation %g, "\
               "reducing LP feastol to %g\n", maxvarboundviol, maxauxviol, SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, TRUE, maxrelconsviol, result) );

   if( *result == SCIP_CUTOFF || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED ||
         *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   assert(*result == SCIP_DIDNOTFIND);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " could not enforce violation %g in regular ways, LP feastol=%g, "\
            "becoming desperate now...\n", maxabsconsviol, SCIPgetLPFeastol(scip)); )

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxvarboundviol) && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));  /*lint !e666*/
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bounds are violated by more than eps, reduced LP "\
               "feasibility tolerance to %g\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxauxviol) && SCIPisPositive(scip,
            SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      /* try whether tighten the LP feasibility tolerance could help
       * maybe it is just some cut that hasn't been taken into account sufficiently
       * in the next enforcement round, we would then also allow even weaker cuts, as we want a minimal cut violation of LP's feastol
       * unfortunately, we do not know the current LP solution primal infeasibility, so sometimes this just repeats without effect
       * until the LP feastol reaches epsilon
       */
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxauxviol / 2.0, SCIPgetLPFeastol(scip) / 10.0)));  /*lint !e666*/
      ++conshdlrdata->ndesperatetightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " reduced LP feasibility tolerance to %g and hope\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   /* try to propagate, if not tried above TODO(?) allow to disable this as well */
   if( !conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;
      int ndelconss = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds, &ndelconss) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* could not find branching candidates even when looking at minimal violated (>eps) expressions
    * now look if we find any unfixed variable that we could still branch on
    */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify > 0 )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " registered %d unfixed variables as branching candidates\n", nnotify); )
      ++conshdlrdata->ndesperatebranch;

      *result = SCIP_INFEASIBLE;  /* enforceConstraints may have changed it to SCIP_DIDNOTFIND */

      return SCIP_OKAY;
   }

   /* if everything is fixed in violated constraints, then let's cut off the node
    * - bound tightening with all vars fixed should prove cutoff, but interval arithmetic overestimates and so the
    *   result may not be conclusive (when constraint violations are small)
    * - if tightenlpfeastol=FALSE, then the LP solution that we try to enforce here may just not be within bounds
    *   sufficiently (see st_e40)
    * - but if the LP solution is really within bounds and since variables are fixed, cutting off the node is actually
    *   not "desperate", but a pretty obvious thing to do
    */
   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " enforcement with max. violation %g failed; cutting off node\n", maxabsconsviol); )
   *result = SCIP_CUTOFF;

   /* it's only "desperate" if the LP solution does not coincide with variable fixings (should we use something tighter than epsilon here?) */
   if( !SCIPisZero(scip, maxvarboundviol) )
      ++conshdlrdata->ndesperatecutoff;

   return SCIP_OKAY;
}

/** separation for all violated constraints to be used by SEPA callbacks */
static
SCIP_RETCODE consSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   unsigned int soltag;
   SCIP_Bool haveviol = FALSE;
   int c;

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   soltag = ++conshdlrdata->lastsoltag;

   /* compute violations */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         haveviol = TRUE;
   }

   /* if none of our constraints are violated, don't attempt separation */
   if( !haveviol )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: skip separation of non-violated constraints\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: separation\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )

   /* call separation */
   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, FALSE, SCIP_INVALID, result) );

   return SCIP_OKAY;
}

/** checks for a linear variable that can be increased or decreased without harming feasibility */
static
void consdataFindUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int poslock;
   int neglock;
   int i;

   assert(conshdlr != NULL);
   assert(consdata != NULL);

   consdata->linvarincr = NULL;
   consdata->linvardecr = NULL;
   consdata->linvarincrcoef = 0.0;
   consdata->linvardecrcoef = 0.0;

   /* root expression is not a sum -> no unlocked linear variable available */
   if( SCIPgetConsExprExprHdlr(consdata->expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return;

   for( i = 0; i < SCIPgetConsExprExprNChildren(consdata->expr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child;

      child = SCIPgetConsExprExprChildren(consdata->expr)[i];
      assert(child != NULL);

      /* check whether the child is a variable expression */
      if( SCIPisConsExprExprVar(child) )
      {
         SCIP_VAR* var = SCIPgetConsExprExprVarVar(child);
         SCIP_Real coef = SCIPgetConsExprExprSumCoefs(consdata->expr)[i];

         if( coef > 0.0 )
         {
            poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         }
         else
         {
            poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         }

         if( SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) - neglock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can decrease x without harming other constraints */
            /* if we have already one candidate, then take the one where the loss in the objective function is less */
            if( (consdata->linvardecr == NULL) ||
               (SCIPvarGetObj(consdata->linvardecr) / consdata->linvardecrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvardecr = var;
               consdata->linvardecrcoef = coef;
            }
         }

         if( SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) - poslock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can increase x without harm */
            /* if we have already one candidate, then take the one where the loss in the objective function is less */
            if( (consdata->linvarincr == NULL) ||
               (SCIPvarGetObj(consdata->linvarincr) / consdata->linvarincrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvarincr = var;
               consdata->linvarincrcoef = coef;
            }
         }
      }
   }

   assert(consdata->linvarincr == NULL || consdata->linvarincrcoef != 0.0);
   assert(consdata->linvardecr == NULL || consdata->linvardecrcoef != 0.0);

#ifdef SCIP_DEBUG
   if( consdata->linvarincr != NULL )
   {
      SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvarincr));
   }
   if( consdata->linvardecr != NULL )
   {
      SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvardecr));
   }
#endif
}

/** Given a solution where every expression constraint is either feasible or can be made feasible by
 *  moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 *
 *  The method assumes that this is always possible and that not all constraints are feasible already.
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* newsol;
   int c;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   *success = FALSE;

   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );
   SCIPdebugMsg(scip, "attempt to make solution from <%s> feasible by shifting linear variable\n",
      sol != NULL ? (SCIPsolGetHeur(sol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(sol)) : "tree") : "LP");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      SCIP_Real viol = 0.0;
      SCIP_Real delta;
      SCIP_Real gap;

      assert(consdata != NULL);

      /* get absolute violation and sign */
      if( consdata->lhsviol > SCIPfeastol(scip) )
         viol = consdata->lhsviol; /* lhs - activity */
      else if( consdata->rhsviol > SCIPfeastol(scip) )
         viol = -consdata->rhsviol; /* rhs - activity */
      else
         continue; /* constraint is satisfied */

      if( consdata->linvarincr != NULL &&
         ((viol > 0.0 && consdata->linvarincrcoef > 0.0) || (viol < 0.0 && consdata->linvarincrcoef < 0.0)) )
      {
         SCIP_VAR* var = consdata->linvarincr;

         /* compute how much we would like to increase var */
         delta = viol / consdata->linvarincrcoef;
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
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy lhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));  /*lint !e613*/

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvarincrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvardecr != NULL &&
         ((viol > 0.0 && consdata->linvardecrcoef < 0.0) || (viol < 0.0 && consdata->linvardecrcoef > 0.0)) )
      {
         SCIP_VAR* var = consdata->linvardecr;

         /* compute how much we would like to decrease var */
         delta = viol / consdata->linvardecrcoef;
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
            SCIP_CALL( SCIPincSolVal(scip, newsol, consdata->linvardecr, delta) );
            /*lint --e{613} */
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy rhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvardecrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so probably we could not make constraint feasible due to variable bounds, thus give up */
      break;
   }

   /* if we have a solution that should satisfy all quadratic constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic
    */
   if( c == nconss && (SCIPisInfinity(scip, SCIPgetUpperbound(scip)) || SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip))) )
   {
      SCIPdebugMsg(scip, "pass solution with objective val %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      assert(conshdlrdata->trysolheur != NULL);
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );

      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** merges constraints that have the same root expression */
static
SCIP_RETCODE presolMergeConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< pointer to store whether at least one constraint could be deleted */
   )
{
   SCIP_HASHMAP* expr2cons;
   SCIP_Bool* updatelocks;
   int* nlockspos;
   int* nlocksneg;
   int c;

   assert(success != NULL);

   *success = FALSE;

   /* not enough constraints available */
   if( nconss <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&expr2cons, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &updatelocks, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      /* ignore deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add expression to the hash map if not seen so far */
      if( !SCIPhashmapExists(expr2cons, (void*)consdata->expr) )
      {
         SCIP_CALL( SCIPhashmapInsertInt(expr2cons, (void*)consdata->expr, c) );
      }
      else
      {
         SCIP_CONSDATA* imgconsdata;
         int idx;

         idx = SCIPhashmapGetImageInt(expr2cons, (void*)consdata->expr);
         assert(idx >= 0 && idx < nconss);

         imgconsdata = SCIPconsGetData(conss[idx]);
         assert(imgconsdata != NULL);
         assert(imgconsdata->expr == consdata->expr);

         SCIPdebugMsg(scip, "merge constraint %g <= %s <= %g with %g <= %s <= %g\n", consdata->lhs,
            SCIPconsGetName(conss[c]), consdata->rhs, imgconsdata->lhs, SCIPconsGetName(conss[idx]), imgconsdata->rhs);

         /* check whether locks need to be updated */
         if( !updatelocks[idx] && ((SCIPisInfinity(scip, -imgconsdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs))
            || (SCIPisInfinity(scip, imgconsdata->rhs) && !SCIPisInfinity(scip, consdata->rhs))) )
         {
            nlockspos[idx] = imgconsdata->nlockspos;
            nlocksneg[idx] = imgconsdata->nlocksneg;
            SCIP_CALL( addLocks(scip, conss[idx], -imgconsdata->nlockspos, -imgconsdata->nlocksneg) );
            updatelocks[idx] = TRUE;
         }

         /* update constraint sides */
         imgconsdata->lhs = MAX(imgconsdata->lhs, consdata->lhs);
         imgconsdata->rhs = MIN(imgconsdata->rhs, consdata->rhs);

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         *success = TRUE;
      }
   }

   /* restore locks of updated constraints */
   if( *success )
   {
      for( c = 0; c < nconss; ++c )
      {
         if( updatelocks[c] )
         {
            SCIP_CALL( addLocks(scip, conss[c], nlockspos[c], nlocksneg[c]) );
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);
   SCIPfreeBufferArray(scip, &updatelocks);
   SCIPhashmapFree(&expr2cons);

   return SCIP_OKAY;
}

/** print statistics for expression handlers */
static
void printExprHdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "Expression Handlers: %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "SimplCalls", "Simplified", "EstimCalls", "#IntEval", "PropCalls", "Cuts", "Cutoffs", "DomReds", "BranchScor", "EstimTime", "PropTime", "IntEvalTi", "SimplifyTi");

   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      SCIP_CONSEXPR_EXPRHDLR* exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      SCIPinfoMessage(scip, file, "  %-17s:", exprhdlr->name);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->nsimplifycalls);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->nsimplified);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->nestimatecalls);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->nintevalcalls);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->npropcalls);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->ncutsfound);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->ncutoffs);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->ndomreds);
      SCIPinfoMessage(scip, file, " %10lld", exprhdlr->nbranchscores);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, exprhdlr->estimatetime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, exprhdlr->proptime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, exprhdlr->intevaltime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, exprhdlr->simplifytime));
      SCIPinfoMessage(scip, file, "\n");
   }
}

/** print statistics for nonlinear handlers */
static
void printNlhdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "Nlhdlrs            : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "EnfoCalls", "#IntEval", "PropCalls", "Detects", "Separated", "Cutoffs", "DomReds", "BranchScor", "Reforms", "DetectTime", "EnfoTime", "PropTime", "IntEvalTi", "ReformTi");

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlr */
      if( !nlhdlr->enabled )
         continue;

      SCIPinfoMessage(scip, file, "  %-17s:", nlhdlr->name);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nenfocalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nintevalcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->npropcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndetections);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nseparated);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ncutoffs);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndomreds);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nbranchscores);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nreformulates);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->detecttime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->enfotime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->proptime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->intevaltime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->reformulatetime));
      SCIPinfoMessage(scip, file, "\n");
   }
}

/** print statistics for constraint handlers */
static
void printConshdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "ConsExpr Enforce   : %10s %10s %10s %10s %10s %10s\n", "WeakSepa", "TightenLP", "DespTghtLP", "DespBranch", "DespCutoff", "ForceLP");
   SCIPinfoMessage(scip, file, "  %-18s", "");
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nweaksepa);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ntightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatetightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatebranch);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatecutoff);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nforcelp);
   SCIPinfoMessage(scip, file, "\n");
   SCIPinfoMessage(scip, file, "ConsExpr Presolve  : %10s\n", "CanonTime");
   SCIPinfoMessage(scip, file, "  %-18s", "");
   SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, conshdlrdata->canonicalizetime));
   SCIPinfoMessage(scip, file, "\n");
}


/*
 * vertex polyhedral separation
 */

/** builds LP used to compute facets of the convex envelope of vertex-polyhedral functions */
static
SCIP_RETCODE buildVertexPolyhedralSeparationLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of (unfixed) variables in vertex-polyhedral functions */
   SCIP_LPI**            lp                  /**< pointer to store created LP */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   int* beg;
   int* ind;
   unsigned int nnonz;
   unsigned int ncols;
   unsigned int nrows;
   unsigned int i;
   unsigned int k;

   assert(scip != NULL);
   assert(lp != NULL);
   assert(nvars > 0);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);

   SCIPdebugMsg(scip, "Building LP for computing facets of convex envelope of vertex-polyhedral function\n");

   /* create lpi to store the LP */
   SCIP_CALL( SCIPlpiCreate(lp, SCIPgetMessagehdlr(scip), "facet finding LP", SCIP_OBJSEN_MINIMIZE) );

   nrows = (unsigned int)nvars + 1;
   ncols = POWEROFTWO((unsigned int)nvars);
   nnonz = (ncols * (nrows + 1)) / 2;

   /* allocate necessary memory; set obj, lb, and ub to zero */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );

   /* calculate nonzero entries in the LP */
   for( i = 0, k = 0; i < ncols; ++i )
   {
      int row;
      unsigned int a;

      /* an upper bound of 1.0 is implied by the last row, but I presume that LP solvers prefer unbounded variables */
      ub[i] = SCIPlpiInfinity(*lp);

      SCIPdebugMsg(scip, "col %i starts at position %d\n", i, k);
      beg[i] = (int)k;
      row = 0;

      /* iterate through the bit representation of i */
      a = 1;
      while( a <= i )
      {
         if( (a & i) != 0 )
         {
            val[k] = 1.0;
            ind[k] = row;

            SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", row, i, k);

            ++k;
         }

         a <<= 1;
         ++row;
         assert(0 <= row && row <= SCIP_MAXVERTEXPOLYDIM);
         assert(POWEROFTWO(row) == a);
      }

      /* put 1 as a coefficient for sum_{i} \lambda_i = 1 row (last row) */
      val[k] = 1.0;
      ind[k] = (int)nrows - 1;
      ++k;
      SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", nrows - 1, i, k);
   }
   assert(k == nnonz);

   /* load all data into LP interface
    * we can assume nrows (=nvars+1) <= ncols (=2^nvars), so we can pass lb as dummy lhs and rhs
    */
   assert(nrows <= ncols);
   SCIP_CALL( SCIPlpiLoadColLP(*lp, SCIP_OBJSEN_MINIMIZE,
      (int)ncols, obj, lb, ub, NULL,
      (int)nrows, lb, lb, NULL,
      (int)nnonz, beg, ind, val) );

   /* for the last row, we can set the rhs to 1.0 already */
   ind[0] = (int)nrows - 1;
   val[0] = 1.0;
   SCIP_CALL( SCIPlpiChgSides(*lp, 1, ind, val, val) );

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** the given facet might not be a valid under(over)estimator, because of numerics and bad fixings; we compute \f$
 * \max_{v \in V} f(v) - (\alpha v + \beta) \f$ (\f$\max_{v \in V} \alpha v + \beta - f(v) \f$) where \f$ V \f$ is the
 * set of vertices of the domain
 */
static
SCIP_Real computeVertexPolyhedralMaxFacetError(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether we check for an over or underestimator */
   SCIP_Real*            funvals,            /**< array containing the evaluation of the function at all corners, length: 2^nvars */
   SCIP_Real*            box,                /**< box for which facet was computed, length: 2*nallvars */
   int                   nallvars,           /**< number of all variables */
   int                   nvars,              /**< number of unfixed variables */
   int*                  nonfixedpos,        /**< indices of unfixed variables, length: nvars */
   SCIP_Real*            facetcoefs,         /**< current facet candidate's coefficients, length: nallvars */
   SCIP_Real             facetconstant       /**< current facet candidate's constant, length: nallvars */
   )
{
   SCIP_Real maxerror;
   SCIP_Real facetval;
   SCIP_Real funval;
   SCIP_Real error;
   unsigned int i;
   unsigned int ncorners;
   unsigned int prev;

   assert(scip != NULL);
   assert(funvals != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(facetcoefs != NULL);

   ncorners = POWEROFTWO(nvars);
   maxerror = 0.0;

   /* check the origin (all variables at lower bound) */
   facetval = facetconstant;
   for( i = 0; i < (unsigned int) nallvars; ++i )
      facetval += facetcoefs[i] * box[2*i];

   /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
   funval = funvals[0];
   if( overestimate )
      error = funval - facetval;
   else
      error = facetval - funval;

   /* update maximum error */
   maxerror = MAX(error, maxerror);

   prev = 0;
   for( i = 1; i < ncorners; ++i )
   {
      unsigned int gray;
      unsigned int diff;
      unsigned int pos;
      int origpos;

      gray = i ^ (i >> 1);
      diff = gray ^ prev;

      /* compute position of unique 1 of diff */
      pos = 0;
      while( (diff >>= 1) != 0 )
         ++pos;
      assert(pos < (unsigned int)nvars);

      origpos = nonfixedpos[pos];

      if( gray > prev )
         facetval += facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);
      else
         facetval -= facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);

      /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
      funval = funvals[gray];
      if( overestimate )
         error = funval - facetval;
      else
         error = facetval - funval;

      /* update  maximum error */
      maxerror = MAX(error, maxerror);

      prev = gray;
   }

   SCIPdebugMsg(scip, "maximum error of facet: %2.8e\n", maxerror);

   return maxerror;
}

/** computes a facet of the convex or concave envelope of a vertex polyhedral function using by solving an LP */
static
SCIP_RETCODE computeVertexPolyhedralFacetLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   int*                  nonfixedpos,        /**< indices of nonfixed variables */
   SCIP_Real*            funvals,            /**< values of function in all corner points (w.r.t. nonfixed variables) */
   int                   nvars,              /**< number of nonfixed variables */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an zero'ed array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LPI* lp;
   SCIP_Real* aux; /* used to transform x^* and then to store LP solution */
   int* inds;
   int ncols;
   int nrows;
   int i;
   SCIP_Real facetvalue;
   SCIP_Real mindomwidth;
   SCIP_RETCODE lpsolveretcode;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(funvals != NULL);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->vp_randnumgen == NULL && conshdlrdata->vp_maxperturb > 0.0 )
   {
      SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->vp_randnumgen, VERTEXPOLY_RANDNUMINITSEED, TRUE) );
   }

   /* construct an LP for this size, if not having one already */
   if( conshdlrdata->vp_lp[nvars] == NULL )
   {
      SCIP_CALL( buildVertexPolyhedralSeparationLP(scip, nvars, &conshdlrdata->vp_lp[nvars]) );
   }
   lp = conshdlrdata->vp_lp[nvars];
   assert(lp != NULL);

   /* get number of cols and rows of separation lp */
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );

   /* number of columns should equal the number of corners = 2^nvars */
   assert(ncols == (int)POWEROFTWO(nvars));

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &aux, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );

   /*
    * set up the described LP on the transformed space
    */

   for( i = 0; i < ncols; ++i )
      inds[i] = i;

   /* compute T^-1(x^*), i.e. T^-1(x^*)_i = (x^*_i - lb_i)/(ub_i - lb_i) */
   mindomwidth = 2*SCIPinfinity(scip);
   for( i = 0; i < nrows-1; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      assert(i < nvars);

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      solval = xstar[varpos];

      if( ub - lb < mindomwidth )
         mindomwidth = ub - lb;

      /* explicitly handle solution which violate bounds of variables (this can happen because of tolerances) */
      if( solval <= lb )
         aux[i] = 0.0;
      else if( solval >= ub )
         aux[i] = 1.0;
      else
         aux[i] = (solval - lb) / (ub - lb);

      /* perturb point to hopefully obtain a facet of the convex envelope */
      if( conshdlrdata->vp_maxperturb > 0.0 )
      {
         assert(conshdlrdata->vp_randnumgen != NULL);

         if( aux[i] == 1.0 )
            aux[i] -= SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else if( aux[i] == 0.0 )
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else
         {
            SCIP_Real perturbation;

            perturbation = MIN( aux[i], 1.0 - aux[i] ) / 2.0;
            perturbation = MIN( perturbation, conshdlrdata->vp_maxperturb );
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, -perturbation, perturbation);
         }
         assert(0.0 < aux[i] && aux[i] < 1.0);
      }

      SCIPdebugMsg(scip, "LP row %d in [%e, %e]\n", i, aux[i], aux[i]);
   }

   /* update LP */
   SCIP_CALL( SCIPlpiChgObj(lp, ncols, inds, funvals) );
   SCIP_CALL( SCIPlpiChgSides(lp, nrows-1, inds, aux, aux) );
   SCIP_CALL( SCIPlpiChgObjsen(lp, overestimate ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE) );

   /* we can stop the LP solve if will not meet the target value anyway, but only if xstar hasn't been perturbed */
   if( conshdlrdata->vp_maxperturb == 0.0 && !SCIPisInfinity(scip, REALABS(targetvalue)) )
   {
      SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_OBJLIM, targetvalue) );
   }
   /* set an iteration limit so we do not run forever */
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPITLIM, 100*ncols) );
   /* since we work with the dual of the LP, primal feastol determines how much we want the computed facet to be the best possible one */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_FEASTOL, SCIPfeastol(scip)) );
   /* since we work with the dual of the LP, dual feastol determines validity of the facet
    * if some ub-lb is small, we need higher accuracy, since below we divide coefs by ub-lb (we moved and scaled the box)
    * thus, we set the dual feastol to be between SCIPepsilon and SCIPfeastol
    */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, MIN(SCIPfeastol(scip), MAX(SCIPepsilon(scip), mindomwidth * SCIPfeastol(scip)))) ); /*lint !e666*/

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPINFO, 1) );
#endif
   /* SCIP_CALL( SCIPlpiWriteLP(lp, "lp.lp") ); */

   /*
    * solve the LP and store the resulting facet for the transformed space
    */
   if( conshdlrdata->vp_dualsimplex )
   {
      lpsolveretcode = SCIPlpiSolveDual(lp);
   }
   else
   {
      lpsolveretcode = SCIPlpiSolvePrimal(lp);
   }
   if( lpsolveretcode == SCIP_LPERROR )
   {
      SCIPdebugMsg(scip, "LP error, aborting.\n");
      goto CLEANUP;
   }
   SCIP_CALL( lpsolveretcode );

   /* any dual feasible solution should provide a valid estimator (and a dual optimal one a facet) */
   if( !SCIPlpiIsDualFeasible(lp) )
   {
      SCIPdebugMsg(scip, "LP not solved to dual feasibility, aborting.\n");
      goto CLEANUP;
   }

   /* get dual solution (facet of convex envelope); again, we have to be careful since the LP can have more rows and
    * columns than needed, in particular, \bar \beta is the last dual multiplier
    */
   SCIP_CALL( SCIPlpiGetSol(lp, NULL, NULL, aux, NULL, NULL) );

   for( i = 0; i < nvars; ++i )
      facetcoefs[nonfixedpos[i]] = aux[i];
   /* last dual multiplier is the constant */
   *facetconstant = aux[nrows - 1];


#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "facet for the transformed problem: ");
   for( i = 0; i < nallvars; ++i )
   {
      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[i], i);
   }
   SCIPdebugMsgPrint(scip, "%3.4e\n", *facetconstant);
#endif

   /*
    * transform the facet to original space and compute value at x^*, i.e., alpha x + beta
    */

   SCIPdebugMsg(scip, "facet in orig. space: ");

   facetvalue = 0.0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      assert(!SCIPisEQ(scip, lb, ub));

      /* alpha_i := alpha_bar_i / (ub_i - lb_i) */
      facetcoefs[varpos] = facetcoefs[varpos] / (ub - lb);

      /* beta = beta_bar - sum_i alpha_i * lb_i */
      *facetconstant -= facetcoefs[varpos] * lb;

      /* evaluate */
      facetvalue += facetcoefs[varpos] * xstar[varpos];

      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[varpos], varpos);
   }
   SCIPdebugMsgPrint(scip, "%3.4e ", *facetconstant);

   /* add beta to the facetvalue: at this point in the code, facetvalue = g(x^*) */
   facetvalue += *facetconstant;

   SCIPdebugMsgPrint(scip, "has value %g, target = %g\n", facetvalue, targetvalue);

    /* if overestimate, then we want facetvalue < targetvalue
    * if underestimate, then we want facetvalue > targetvalue
    * if none holds, give up
    * so maybe here we should check against the minimal violation
    */
   if( overestimate == (facetvalue > targetvalue) )
   {
      SCIPdebugMsg(scip, "missed the target, facetvalue %g targetvalue %g, overestimate=%d\n", facetvalue, targetvalue, overestimate);
      goto CLEANUP;
   }

   /* if we made it until here, then we have a nice facet */
   *success = TRUE;

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &aux);

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a univariant vertex polyhedral function
 *
 * In other words, compute the line that passes through two given points.
 */
static
SCIP_RETCODE computeVertexPolyhedralFacetUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             left,               /**< left coordinate */
   SCIP_Real             right,              /**< right coordinate */
   SCIP_Real             funleft,            /**< value of function in left coordinate */
   SCIP_Real             funright,           /**< value of function in right coordinate */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoef,          /**< buffer to store coefficient of facet defining inequality */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   assert(scip != NULL);
   assert(SCIPisLE(scip, left, right));
   assert(!SCIPisInfinity(scip, -left));
   assert(!SCIPisInfinity(scip, right));
   assert(SCIPisFinite(funleft) && funleft != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(funright) && funright != SCIP_INVALID);  /*lint !e777*/
   assert(success != NULL);
   assert(facetcoef != NULL);
   assert(facetconstant != NULL);

   *facetcoef = (funright - funleft) / (right - left);
   *facetconstant = funleft - *facetcoef * left;

   *success = TRUE;

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a bivariate vertex polyhedral function */
static
SCIP_RETCODE computeVertexPolyhedralFacetBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real             p1[2],              /**< first vertex of box */
   SCIP_Real             p2[2],              /**< second vertex of box */
   SCIP_Real             p3[2],              /**< third vertex of box */
   SCIP_Real             p4[2],              /**< forth vertex of box */
   SCIP_Real             p1val,              /**< value in p1 */
   SCIP_Real             p2val,              /**< value in p2 */
   SCIP_Real             p3val,              /**< value in p3 */
   SCIP_Real             p4val,              /**< value in p4 */
   SCIP_Real             xstar[2],           /**< point to be separated */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least 2 */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   SCIP_Real alpha, beta, gamma_, delta;
   SCIP_Real xstarval, candxstarval = 0.0;
   int leaveout;

   assert(scip != NULL);
   assert(success != NULL);
   assert(SCIPisFinite(p1val) && p1val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p2val) && p2val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p3val) && p3val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p4val) && p4val != SCIP_INVALID);  /*lint !e777*/
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* if we want an underestimator, flip f(x,y), i.e., do as if we compute an overestimator for -f(x,y) */
   if( !overestimate )
   {
      p1val = -p1val;
      p2val = -p2val;
      p3val = -p3val;
      p4val = -p4val;
      targetvalue = -targetvalue;
   }

   SCIPdebugMsg(scip, "p1 = (%g, %g), f(p1) = %g\n", p1[0], p1[1], p1val);
   SCIPdebugMsg(scip, "p2 = (%g, %g), f(p2) = %g\n", p2[0], p2[1], p2val);
   SCIPdebugMsg(scip, "p3 = (%g, %g), f(p3) = %g\n", p3[0], p3[1], p3val);
   SCIPdebugMsg(scip, "p4 = (%g, %g), f(p4) = %g\n", p4[0], p4[1], p4val);

   /* Compute coefficients alpha, beta, gamma (>0), delta such that
    *   alpha*x + beta*y + gamma*z = delta
    * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
    * the fourth corner point lies below this hyperplane.
    * Since we assume that f is vertex-polyhedral, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
    *    alpha*x + beta*y - delta <= -gamma * f(x,y),
    * or, equivalently,
    *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
    */
   for( leaveout = 1; leaveout <= 4; ++leaveout )
   {
      switch( leaveout)
      {
         case 1 :
            /* get hyperplane through p2, p3, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p1, then go to next candidate */
            if( alpha * p1[0] + beta * p1[1] + gamma_ * p1val - delta > 0.0 )
               continue;
            break;

         case 2 :
            /* get hyperplane through p1, p3, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p2, then go to next candidate */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val - delta > 0.0 )
               continue;
            break;

         case 3 :
            /* get hyperplane through p1, p2, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p3, then go to next candidate */
            if( alpha * p3[0] + beta * p3[1] + gamma_ * p3val - delta > 0.0 )
               continue;
            break;

         case 4 :
            /* get hyperplane through p1, p2, p3 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p4, then stop */
            if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val - delta > 0.0 )
               continue;
            break;

         default: /* only for lint */
            alpha = SCIP_INVALID;
            beta = SCIP_INVALID;
            gamma_ =  SCIP_INVALID;
            delta = SCIP_INVALID;
            break;
      }

      /* check if bad luck: should not happen if numerics are fine */
      if( SCIPisZero(scip, gamma_) )
         continue;
      assert(!SCIPisNegative(scip, gamma_));

      /* if coefficients become tiny because division by gamma makes them < SCIPepsilon(scip), then skip, too */
      if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
         ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta/gamma_)) )
         continue;

      SCIPdebugMsg(scip, "alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

      /* value of hyperplane candidate in xstar */
      xstarval = -alpha/gamma_ * xstar[0] -beta/gamma_ * xstar[1] + delta/gamma_;

      /* if reaching target and first or better than previous candidate, then update */
      if( xstarval <= targetvalue && (!*success || xstarval < candxstarval) )
      {
         /* flip hyperplane */
         if( !overestimate )
            gamma_ = -gamma_;

         facetcoefs[0] = -alpha/gamma_;
         facetcoefs[1] = -beta/gamma_;
         *facetconstant = delta/gamma_;

         *success = TRUE;
         candxstarval = xstarval;
      }
   }

   return SCIP_OKAY;
}

/** hash key retrieval function for bilinear term entries */
static
SCIP_DECL_HASHGETKEY(bilinearTermsGetHashkey)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int idx;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userptr;
   assert(conshdlrdata != NULL);

   idx = ((int)(size_t)elem) - 1;
   assert(idx >= 0 && idx < conshdlrdata->nbilinterms);

   return (void*)&conshdlrdata->bilinterms[idx];
}

/** returns TRUE iff the bilinear term entries are equal */
static
SCIP_DECL_HASHKEYEQ(bilinearTermsIsHashkeyEq)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_BILINTERM* entry1;
   SCIP_CONSEXPR_BILINTERM* entry2;

   /* get corresponding entries */
   entry1 = (SCIP_CONSEXPR_BILINTERM*)key1;
   entry2 = (SCIP_CONSEXPR_BILINTERM*)key2;
   assert(entry1->x != NULL && entry1->y != NULL);
   assert(entry2->x != NULL && entry2->y != NULL);
   assert(SCIPvarCompare(entry1->x, entry1->y) < 1);
   assert(SCIPvarCompare(entry2->x, entry2->y) < 1);

   return entry1->x == entry2->x && entry1->y == entry2->y;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(bilinearTermsGetHashkeyVal)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_BILINTERM* entry;

   entry = (SCIP_CONSEXPR_BILINTERM*)key;
   assert(entry->x != NULL && entry->y != NULL);
   assert(SCIPvarCompare(entry->x, entry->y) < 1);

   return SCIPhashTwo(SCIPvarGetIndex(entry->x), SCIPvarGetIndex(entry->y));
}

/** resizes array of bilinear terms */
static
SCIP_RETCODE bilinearTermsResize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int                   reqsize             /**< required size */
   )
{
   int newsize;

   assert(conshdlrdata != NULL);

   /* check whether array is large enough */
   if( reqsize <= conshdlrdata->bilintermssize )
      return SCIP_OKAY;

   /* compute new size */
   newsize = SCIPcalcMemGrowSize(scip, reqsize);
   assert(reqsize <= newsize);

   /* realloc array */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bilinterms, conshdlrdata->bilintermssize,
      newsize) );
   conshdlrdata->bilintermssize = newsize;

   return SCIP_OKAY;
}

/** stores the variables of a bilinear term in the data of the constraint handler */
static
SCIP_RETCODE bilinearTermsInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   )
{
   SCIP_CONSEXPR_BILINTERM* term;

   assert(conshdlrdata != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(nlockspos >= 0);
   assert(nlocksneg >= 0);

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   /* ensure size of bilinterms array */
   SCIP_CALL( bilinearTermsResize(scip, conshdlrdata, conshdlrdata->nbilinterms + 1) );

   /* set values in the created bilinear term */
   term = &conshdlrdata->bilinterms[conshdlrdata->nbilinterms];
   assert(term != NULL);
   term->x = x;
   term->y = y;
   term->auxvar = auxvar;
   term->nlockspos = nlockspos;
   term->nlocksneg = nlocksneg;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, x) );
   SCIP_CALL( SCIPcaptureVar(scip, y) );
   if( auxvar != NULL )
   {
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   /* increase the total number of bilinear terms */
   ++(conshdlrdata->nbilinterms);

   return SCIP_OKAY;
}

/** iterates through all expressions of all expression constraints and adds the corresponding bilinear terms to the
 *  hash table
 */
static
SCIP_RETCODE bilinearTermsInsertAll(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< expression constraints */
   int                   nconss              /**< total number of expression constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPRHDLR* producthdlr;
   SCIP_CONSEXPR_EXPRHDLR* powhdlr;
   int c;

   assert(conss != NULL || nconss == 0);

   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether the bilinear terms have been stored already */
   if( conshdlrdata->bilinterms != NULL )
      return SCIP_OKAY;

   /* create and initialize iterator */
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR);

   /* get product and pow expression handlers */
   producthdlr = SCIPgetConsExprExprHdlrProduct(conshdlr);
   powhdlr = SCIPgetConsExprExprHdlrPower(conshdlr);

   /* iterate through all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONSEXPR_EXPR* expr;

      assert(conss != NULL && conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* iterate through all expressions */
      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIP_CONSEXPR_EXPR** children = SCIPgetConsExprExprChildren(expr);
         SCIP_VAR* x = NULL;
         SCIP_VAR* y = NULL;

         /* check whether the expression is of the form f(..)^2 */
         if( SCIPgetConsExprExprHdlr(expr) == powhdlr && SCIPgetConsExprExprPowExponent(expr) == 2.0 )
         {
            x = SCIPgetConsExprExprAuxVar(children[0]);
            y = x;
         }
         /* check whether the expression is of the form f(..) * g(..) */
         else if( SCIPgetConsExprExprHdlr(expr) == producthdlr && SCIPgetConsExprExprNChildren(expr) == 2 )
         {
            x = SCIPgetConsExprExprAuxVar(children[0]);
            y = SCIPgetConsExprExprAuxVar(children[1]);
         }

         /* add variables to the hash table */
         if( x != NULL && y != NULL )
         {
            SCIP_CALL( bilinearTermsInsert(scip, conshdlrdata, x, y, SCIPgetConsExprExprAuxVar(expr),
               SCIPgetConsExprExprNLocksPos(expr), SCIPgetConsExprExprNLocksNeg(expr)) );
         }
      }
   }

   /* release iterator */
   SCIPexpriteratorFree(&it);

   /* create hash table and insert stored bilinear terms */
   if( conshdlrdata->nbilinterms > 0 )
   {
      int i;

      assert(conshdlrdata->bilinhashtable == NULL);

      SCIP_CALL( SCIPhashtableCreate(&conshdlrdata->bilinhashtable, SCIPblkmem(scip), conshdlrdata->nbilinterms,
         bilinearTermsGetHashkey, bilinearTermsIsHashkeyEq, bilinearTermsGetHashkeyVal,
         (void*)conshdlrdata) );

      for( i = 0; i < conshdlrdata->nbilinterms; ++i )
      {
         /* insert the index of the bilinear term into the hash table; note that the index of the i-th element is (i+1)
          * because zero can not be inserted into hash table
          */
         SCIP_CALL( SCIPhashtableInsert(conshdlrdata->bilinhashtable, (void*)(size_t)(i+1)) );/*lint !e571 !e776*/
      }
   }

   return SCIP_OKAY;
}

/** frees array of bilinear terms and hash table */
static
SCIP_RETCODE bilinearTermsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int i;

   assert(conshdlrdata != NULL);

   /* check whether bilinear terms have been stored */
   if( conshdlrdata->bilinterms == NULL )
   {
      assert(conshdlrdata->bilinterms == NULL);
      assert(conshdlrdata->nbilinterms == 0);
      assert(conshdlrdata->bilintermssize == 0);

      return SCIP_OKAY;
   }

   /* release variables */
   for( i = 0; i < conshdlrdata->nbilinterms; ++i )
   {
      /* it might be that there is a bilinear term without a corresponding auxiliary variable */
      if( conshdlrdata->bilinterms[i].auxvar != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].auxvar) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].y) );
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].x) );
   }

   /* free hash table */
   if( conshdlrdata->bilinhashtable != NULL )
   {
      SCIPhashtableFree(&conshdlrdata->bilinhashtable);
   }

   /* free bilinterms array; reset counters */
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bilinterms, conshdlrdata->bilintermssize);
   conshdlrdata->nbilinterms = 0;
   conshdlrdata->bilintermssize = 0;

   return SCIP_OKAY;
}

/** returns whether the variable of a given variable expression is a candidate for presolSingleLockedVars(), i.e.,
 *  the variable is only contained in a single expression constraint, has no objective coefficient, has finite
 *  variable bounds, and is not binary
 */
static
SCIP_Bool isSingleLockedCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_VAR* var;

   assert(SCIPisConsExprExprVar(expr));

   var = SCIPgetConsExprExprVarVar(expr);
   assert(var != NULL);

   return SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) == SCIPgetConsExprExprNLocksNeg(expr)
      && SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) == SCIPgetConsExprExprNLocksPos(expr)
      && SCIPgetConsExprExprVarNConss(expr) == 1 && SCIPisZero(scip, SCIPvarGetObj(var))
      && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))
      && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY;
}

/** removes all variable expressions that are contained in a given expression from a hash map */
static
SCIP_RETCODE removeSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_ITERATOR* it,               /**< expression iterator */
   SCIP_HASHMAP*         exprcands           /**< map to hash variable expressions */
   )
{
   SCIP_CONSEXPR_EXPR* e;

   for( e = SCIPexpriteratorRestartDFS(it, expr); !SCIPexpriteratorIsEnd(it); e = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      if( SCIPisConsExprExprVar(e) && SCIPhashmapExists(exprcands, (void*)e) )
      {
         SCIP_CALL( SCIPhashmapRemove(exprcands, (void*)e) );
      }
   }

   return SCIP_OKAY;
}

/** presolving method to fix a variable x_i to one of its bounds if the variable is only contained in a single
 *  expression contraint g(x) <= rhs (>= lhs) if g is concave (convex) in x_i;  if a continuous variable has bounds
 *  [0,1], then the variable type is changed to be binary; otherwise a bound disjunction constraint is added
 *
 *  @todo the same reduction can be applied if g(x) is not concave, but monotone in x_i for g(x) <= rhs
 *  @todo extend this to cases where a variable can appear in a monomial with an exponent, essentially relax
 *    g(x) to sum_i [a_i,b_i] x^{p_i} for a single variable x and try to conclude montonicity or convexity/concavity
 *    on this (probably have one or two flags per variable and update this whenever another x^{p_i} is found)
 */
static
SCIP_RETCODE presolSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   int*                  nchgvartypes,       /**< pointer to store the total number of changed variable types */
   int*                  naddconss,          /**< pointer to store the total number of added constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_CONSEXPR_EXPR** singlelocked;
   SCIP_CONSEXPR_EXPRHDLR* prodhdlr;
   SCIP_CONSEXPR_EXPRHDLR* powhdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* exprcands;
   SCIP_Bool hasbounddisj;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int nsinglelocked = 0;
   int i;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(nchgvartypes != NULL);
   assert(naddconss != NULL);
   assert(infeasible != NULL);

   *nchgvartypes = 0;
   *naddconss = 0;
   *infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* only consider constraints with one finite side */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
      return SCIP_OKAY;

   /* only consider sum expressions */
   if( consdata->expr == NULL || SCIPgetConsExprExprHdlr(consdata->expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   /* remember which side is finite */
   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* get product and power handlers */
   prodhdlr = SCIPgetConsExprExprHdlrProduct(conshdlr);
   powhdlr = SCIPgetConsExprExprHdlrPower(conshdlr);

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(&exprcands, SCIPblkmem(scip), consdata->nvarexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &singlelocked, consdata->nvarexprs) );

   /* check all variable expressions for single locked variables */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);

      if( isSingleLockedCand(scip, consdata->varexprs[i]) )
      {
         SCIP_CALL( SCIPhashmapInsert(exprcands, (void*)consdata->varexprs[i], NULL) );
         singlelocked[nsinglelocked++] = consdata->varexprs[i];
      }
   }
   SCIPdebugMsg(scip, "found %d single locked variables for constraint %s\n", nsinglelocked, SCIPconsGetName(cons));

   if( nsinglelocked > 0 )
   {
      SCIP_CONSEXPR_EXPR** children;
      SCIP_CONSEXPR_ITERATOR* it;
      int nchildren;

      children = SCIPgetConsExprExprChildren(consdata->expr);
      nchildren = SCIPgetConsExprExprNChildren(consdata->expr);

      /* create iterator */
      SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
      SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR);

      for( i = 0; i < nchildren; ++i )
      {
         SCIP_CONSEXPR_EXPR* child;
         SCIP_Real coef;

         child = children[i];
         assert(child != NULL);
         coef = SCIPgetConsExprExprSumCoefs(consdata->expr)[i];

         /* ignore linear terms */
         if( SCIPisConsExprExprVar(child) )
            continue;

         /* consider products prod_j f_j(x); ignore f_j(x) if it is a single variable, otherwise iterate through the
          * expression that represents f_j and remove each variable expression from exprcands
          */
         else if( SCIPgetConsExprExprHdlr(child) == prodhdlr )
         {
            int j;

            for( j = 0; j < SCIPgetConsExprExprNChildren(child); ++j )
            {
               SCIP_CONSEXPR_EXPR* grandchild = SCIPgetConsExprExprChildren(child)[j];

               if( !SCIPisConsExprExprVar(grandchild) )
               {
                  /* mark all variable expressions that are contained in the expression */
                  SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
               }
            }
         }
         /* fixing a variable x to one of its bounds is only valid for ... +x^p >= lhs or ... -x^p <= rhs if p = 2k
          * for an integer k > 1
          */
         else if( SCIPgetConsExprExprHdlr(child) == powhdlr )
         {
            SCIP_CONSEXPR_EXPR* grandchild = SCIPgetConsExprExprChildren(child)[0];
            SCIP_Real exponent = SCIPgetConsExprExprPowExponent(child);
            SCIP_Bool valid;

            /* check for even integral exponent */
            valid = exponent > 1.0 && fmod(exponent, 2.0) == 0.0;

            if( !valid || !SCIPisConsExprExprVar(grandchild) || (hasrhs && coef > 0.0) || (haslhs && coef < 0.0) )
            {
               /* mark all variable expressions that are contained in the expression */
               SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
            }
         }
         /* all other cases cannot be handled */
         else
         {
            /* mark all variable expressions that are contained in the expression */
            SCIP_CALL( removeSingleLockedVars(scip, child, it, exprcands) );
         }
      }

      /* free expression iterator */
      SCIPexpriteratorFree(&it);
   }

   /* check whether the bound disjunction constraint handler is available */
   hasbounddisj = SCIPfindConshdlr(scip, "bounddisjunction") != NULL;

   /* fix variable to one of its bounds by either changing its variable type or adding a disjunction constraint */
   for( i = 0; i < nsinglelocked; ++i )
   {
      /* only consider expressions that are still contained in the exprcands map */
      if( SCIPhashmapExists(exprcands, (void*)singlelocked[i]) )
      {
         SCIP_CONS* newcons;
         SCIP_VAR* vars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];
         SCIP_VAR* var;

         var = SCIPgetConsExprExprVarVar(singlelocked[i]);
         assert(var != NULL);
         SCIPdebugMsg(scip, "found single locked variable %s in [%g,%g] that can be fixed to one of its bounds\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

         /* try to change the variable type to binary */
         if( conshdlrdata->checkvarlocks == 't' && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0) && SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0) )
         {
            assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            ++(*nchgvartypes);

            if( *infeasible )
            {
               SCIPdebugMsg(scip, "detect infeasibility after changing variable type of <%s>\n", SCIPvarGetName(var));
               break;
            }
         }
         /* add bound disjunction constraint if bounds of the variable are finite */
         else if( hasbounddisj && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            vars[0] = var;
            vars[1] = var;
            boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
            boundtypes[1] = SCIP_BOUNDTYPE_UPPER;
            bounds[0] = SCIPvarGetUbGlobal(var);
            bounds[1] = SCIPvarGetLbGlobal(var);

            SCIPdebugMsg(scip, "add bound disjunction constraint for %s\n", SCIPvarGetName(var));

            /* create, add, and release bound disjunction constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quadvarbnddisj_%s", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &newcons, name, 2, vars, boundtypes, bounds, TRUE, TRUE,
               TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            ++(*naddconss);
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &singlelocked);
   SCIPhashmapFree(&exprcands);

   return SCIP_OKAY;
}

/** @} */

/*
 * Callback methods of constraint handler
 */

/** upgrades quadratic constraint to expr constraint */
static
SCIP_DECL_QUADCONSUPGD(quadconsUpgdExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* consexprhdlr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* varexpr;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* powexpr;
   SCIP_CONSEXPR_EXPR* twoexprs[2];
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_BILINTERM* bilinterm;
   int pos;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMsg(scip, "quadconsUpgdExpr called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* no interest in linear constraints */
   if( SCIPgetNQuadVarTermsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* signal that we need more memory */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   if( SCIPgetNBilinTermsQuadratic(scip, cons) > 0 )
   {
      /* we will need SCIPfindQuadVarTermQuadratic later, so ensure now that quad var terms are sorted */
      SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, cons) );
   }

   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &expr, 0, NULL, NULL, 0.0) );

   /* append linear terms */
   for( i = 0; i < SCIPgetNLinearVarsQuadratic(scip, cons); ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, SCIPgetLinearVarsQuadratic(scip, cons)[i]) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexpr, SCIPgetCoefsLinearVarsQuadratic(scip, cons)[i]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
   }

   /* array to store variable expression for each quadratic variable */
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNQuadVarTermsQuadratic(scip, cons)) );

   /* create var exprs for quadratic vars; append linear and square part of quadratic terms */
   for( i = 0; i < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++i )
   {
      quadvarterm = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexprs[i], quadvarterm->var) );

      if( quadvarterm->lincoef != 0.0 )
      {
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexprs[i], quadvarterm->lincoef) );
      }

      if( quadvarterm->sqrcoef != 0.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr, varexprs[i], 2.0) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, powexpr, quadvarterm->sqrcoef) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr) );
      }
   }

   /* append bilinear terms */
   for( i = 0; i < SCIPgetNBilinTermsQuadratic(scip, cons); ++i)
   {
      bilinterm = &SCIPgetBilinTermsQuadratic(scip, cons)[i];

      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var1, &pos) );
      assert(pos >= 0);
      assert(pos < SCIPgetNQuadVarTermsQuadratic(scip, cons));
      assert(SCIPgetQuadVarTermsQuadratic(scip, cons)[pos].var == bilinterm->var1);
      twoexprs[0] = varexprs[pos];

      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var2, &pos) );
      assert(pos >= 0);
      assert(pos < SCIPgetNQuadVarTermsQuadratic(scip, cons));
      assert(SCIPgetQuadVarTermsQuadratic(scip, cons)[pos].var == bilinterm->var2);
      twoexprs[1] = varexprs[pos];

      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 2, twoexprs, 1.0) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, prodexpr, bilinterm->coef) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   }

   /* release variable expressions */
   for( i = 0; i < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++i )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[i]) );
   }

   SCIPfreeBufferArray(scip, &varexprs);

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsExpr(scip, upgdconss, SCIPconsGetName(cons),
      expr, SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons),
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

   SCIPdebugMsg(scip, "created expr constraint:\n");
   SCIPdebugPrintCons(scip, *upgdconss, NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** upgrades nonlinear constraint to expr constraint */
static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdExpr)
{
   SCIP_CONSHDLR* consexprhdlr;
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_CONSEXPR_EXPR* expr;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));
   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   SCIPdebugMsg(scip, "nonlinconsUpgdExpr called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   /* try to create a cons_expr expression from an expression graph node */
   SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &expr, exprgraph, node) );

   /* if that didn't work, then because we do not support a certain expression type yet -> no upgrade */
   if( expr == NULL )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   if( SCIPgetNLinearVarsNonlinear(scip, cons) > 0 )
   {
      /* add linear terms */
      SCIP_CONSEXPR_EXPR* varexpr;
      int i;

      /* ensure expr is a sum expression */
      if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(consexprhdlr) )
      {
         SCIP_CONSEXPR_EXPR* sumexpr;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr, 1, &expr, NULL, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

         expr = sumexpr;
      }

      for( i = 0; i < SCIPgetNLinearVarsNonlinear(scip, cons); ++i )
      {
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, SCIPgetLinearVarsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexpr, SCIPgetLinearCoefsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
      }
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsExpr(scip, upgdconss, SCIPconsGetName(cons),
      expr, SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

   SCIPdebugMsg(scip, "created expr constraint:\n");
   SCIPdebugPrintCons(scip, *upgdconss, NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExpr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(valid != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* create basic data of constraint handler and include it to scip */
   SCIP_CALL( includeConshdlrExprBasic(scip) );

   /* copy expression and nonlinear handlers */
   SCIP_CALL( copyConshdlrExprExprHdlr(scip, conshdlr, valid) );

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      if( exprhdlr->freehdlr != NULL )
      {
         SCIP_CALL( (*exprhdlr->freehdlr)(scip, conshdlr, exprhdlr, &exprhdlr->data) );
      }

      /* free clocks */
      SCIP_CALL( SCIPfreeClock(scip, &(exprhdlr)->simplifytime) );
      SCIP_CALL( SCIPfreeClock(scip, &(exprhdlr)->intevaltime) );
      SCIP_CALL( SCIPfreeClock(scip, &(exprhdlr)->proptime) );
      SCIP_CALL( SCIPfreeClock(scip, &(exprhdlr)->estimatetime) );

      SCIPfreeMemory(scip, &exprhdlr->name);
      SCIPfreeMemoryNull(scip, &exprhdlr->desc);

      SCIPfreeMemory(scip, &exprhdlr);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize);

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      if( nlhdlr->freehdlrdata != NULL )
      {
         SCIP_CALL( (*nlhdlr->freehdlrdata)(scip, nlhdlr, &nlhdlr->data) );
      }

      /* free clocks */
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->detecttime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->enfotime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->proptime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->intevaltime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->reformulatetime) );

      SCIPfreeMemory(scip, &nlhdlr->name);
      SCIPfreeMemoryNull(scip, &nlhdlr->desc);

      SCIPfreeMemory(scip, &nlhdlr);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize);
   conshdlrdata->nlhdlrssize = 0;

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nexprconsupgrades; ++i )
   {
      assert(conshdlrdata->exprconsupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->exprconsupgrades[i]);  /*lint !e866*/
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->exprconsupgrades, conshdlrdata->exprconsupgradessize);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->canonicalizetime) );

   assert(conshdlrdata->vp_randnumgen == NULL);
#ifndef NDEBUG
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
      assert(conshdlrdata->vp_lp[i] == NULL);
#endif

   assert(conshdlrdata->branchrandnumgen == NULL);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* make sure current activity tags in expressions are invalid, because we start catching variable events only now */
   conshdlrdata->lastboundrelax = ++conshdlrdata->curboundstag;
   /* set to 1 so it is larger than initial value of lastenforound in exprs */
   conshdlrdata->enforound = 1;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(conss[i])) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
   }

   /* sort nonlinear handlers by priority, in decreasing order */
   if( conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, nlhdlrCmp, conshdlrdata->nnlhdlrs);

   /* get heuristics for later use */
   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* reset statistics in expression handlers */
   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      exprhdlr->nestimatecalls = 0;
      exprhdlr->nintevalcalls = 0;
      exprhdlr->npropcalls = 0;
      exprhdlr->ncutsfound = 0;
      exprhdlr->ncutoffs = 0;
      exprhdlr->ndomreds = 0;
      exprhdlr->nbranchscores = 0;
      exprhdlr->nsimplifycalls = 0;
      exprhdlr->nsimplified = 0;

      SCIP_CALL( SCIPresetClock(scip, exprhdlr->estimatetime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->proptime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->intevaltime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->simplifytime) );
   }

   /* reset statistics in nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      nlhdlr->nenfocalls = 0;
      nlhdlr->nintevalcalls = 0;
      nlhdlr->npropcalls = 0;
      nlhdlr->nseparated = 0;
      nlhdlr->ncutoffs = 0;
      nlhdlr->ndomreds = 0;
      nlhdlr->nbranchscores = 0;
      nlhdlr->ndetections = 0;

      SCIP_CALL( SCIPresetClock(scip, nlhdlr->detecttime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->proptime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->intevaltime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->reformulatetime) );
   }

   /* reset statistics in constraint handler */
   conshdlrdata->nweaksepa = 0;
   conshdlrdata->ntightenlp = 0;
   conshdlrdata->ndesperatebranch = 0;
   conshdlrdata->ndesperatecutoff = 0;
   conshdlrdata->ndesperatetightenlp = 0;
   conshdlrdata->nforcelp = 0;
   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->canonicalizetime) );

#ifdef ENFOLOGFILE
   ENFOLOG( enfologfile = fopen(ENFOLOGFILE, "w"); )
#endif

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS** consssorted;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( nconss > 0 )
   {
      /* for better performance of dropVarEvents, we sort by index, descending */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortDownPtr((void**)consssorted, SCIPcompareConsExprIndex, nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }

      SCIPfreeBufferArray(scip, &consssorted);
   }

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   if( conshdlrdata->vp_randnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->vp_randnumgen);

   /* free LPs used to construct facets of envelops of vertex-polyhedral functions */
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
   {
      if( conshdlrdata->vp_lp[i] != NULL )
      {
         SCIP_CALL( SCIPlpiFree(&conshdlrdata->vp_lp[i]) );
      }
   }

   if( conshdlrdata->branchrandnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->branchrandnumgen);

   ENFOLOG(
      if( enfologfile != NULL )
      {
         fclose(enfologfile);
         enfologfile = NULL;
      }
   )

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreExpr)
{  /*lint --e{715}*/

   /* remove auxiliary variables when a restart has happened; this ensures that the previous branch-and-bound tree
    * removed all of his captures on variables
    */
   if( SCIPgetNRuns(scip) > 1 )
   {
      SCIP_CALL( freeAuxVars(scip, conshdlr, conss, nconss) );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreExpr)
{  /*lint --e{715}*/
   SCIP_Bool infeasible;

   if( nconss == 0 )
      return SCIP_OKAY;

   /* skip some extra work if already known to be infeasible */
   if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
      return SCIP_OKAY;

   /* simplify constraints and replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   /* currently SCIP does not offer to communicate this,
    * but at the moment this can only become true if canonicalizeConstraints called detectNlhdlrs (which it doesn't do in EXITPRESOLVE stage)
    * or if a constraint expression became constant
    */
   assert(!infeasible);

   /* tell SCIP that we have something nonlinear */
   SCIPenableNLP(scip);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   /* skip a number of initializations if we are already infeasible
    * if infeasibility was found by our boundtightening, then curvature check may also fail as some exprhdlr (e.g., pow)
    * assumes nonempty activities in expressions
    */
   if( SCIPgetStatus(scip) != SCIP_STATUS_INFEASIBLE )
   {
      SCIP_CONSDATA* consdata;
      int c;

      for( c = 0; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
         assert(consdata != NULL);
         assert(consdata->expr != NULL);

         /* check for a linear variable that can be increase or decreased without harming feasibility */
         consdataFindUnlockedLinearVar(scip, conshdlr, consdata);

         /* call curvature detection of expression handlers; TODO do we really need this? */
         SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, consdata->expr) );

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
      }
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initialize nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[i];
      if( nlhdlr->init != NULL )
      {
         SCIP_CALL( (*nlhdlr->init)(scip, nlhdlr) );
      }
   }

   if( conshdlrdata->branchpscostweight > 0.0 )
   {
      SCIP_CALL( SCIPgetCharParam(scip, "branching/lpgainnormalize", &(conshdlrdata->branchpscostupdatestrategy)) );
      if( strchr("lds", conshdlrdata->branchpscostupdatestrategy) == NULL )
      {
         SCIPerrorMessage("branching/lpgainnormalize strategy %c unknown\n", conshdlrdata->branchpscostupdatestrategy);
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   int c;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   /* call deinitialization callbacks of expression and nonlinear handlers
    * free nonlinear handlers information from expressions
    * remove auxiliary variables from expressions, if not restarting; otherwise do so in CONSINITPRE
    */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIPdebugMsg(scip, "exitsepa and free nonlinear handler data for expression %p\n", (void*)expr);

         /* remove nonlinear handlers in expression and their data and auxiliary variables if not restarting */
         SCIP_CALL( freeEnfoData(scip, conshdlr, expr, !restart) );
      }
   }

   SCIPexpriteratorFree(&it);

   /* deinitialize nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CONSEXPR_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[i];
      if( nlhdlr->exit != NULL )
      {
         SCIP_CALL( (*nlhdlr->exit)(scip, nlhdlr) );
      }
   }

   /* free nonlinear row representations */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   /* free hash table for bilinear terms */
   SCIP_CALL( bilinearTermsFree(scip, conshdlrdata) );

   /* reset flag to allow another call of presolSingleLockedVars() after a restart */
   conshdlrdata->checkedvarlocks = FALSE;

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);

   /* constraint locks should have been removed */
   assert((*consdata)->nlockspos == 0);
   assert((*consdata)->nlocksneg == 0);

   /* free variable expressions */
   SCIP_CALL( freeVarExprs(scip, *consdata) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*consdata)->expr) );

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransExpr)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( copyExpr(scip, scip, conshdlr, sourcedata->expr, &targetexpr, transformVar, NULL) );
   assert(targetexpr != NULL);  /* copyExpr cannot fail if source and target scip are the same */

   /* create transformed cons (captures targetexpr) */
   SCIP_CALL( SCIPcreateConsExpr(scip, targetcons, SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
      SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
      SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons)) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpExpr)
{
   /* register non linear handlers TODO: do we want this here? */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss, infeasible) );

   /* if creating auxiliary variables detected an infeasible (because of bounds), stop initing lp */
   if( *infeasible )
      return SCIP_OKAY;

   /* call seaparation initialization callbacks of the expression handlers */
   SCIP_CALL( initSepa(scip, conshdlr, conss, nconss, infeasible) );

   /* collect all bilinear terms */
   SCIP_CALL( bilinearTermsInsertAll(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIP_RESULT propresult;
   unsigned int soltag;
   int nchgbds;
   int ndelconss;
   int nnotify;
   int c;

   soltag = ++conshdlrdata->lastsoltag;

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], NULL, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* try to propagate
    * TODO obey propinenfo parameter, but we need something to recognize cutoff
    */
   nchgbds = 0;
   ndelconss = 0;
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds, &ndelconss) );

   if( (propresult == SCIP_CUTOFF) || (propresult == SCIP_REDUCEDDOM) )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* register all unfixed variables in all violated constraints as branching candidates */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );
   if( nnotify > 0 )
   {
      SCIPdebugMsg(scip, "registered %d external branching candidates\n", nnotify);

      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "could not find branching candidates, forcing to solve LP\n");
   *result = SCIP_SOLVELP;
   ++conshdlrdata->nforcelp;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   SCIP_Bool          maypropfeasible;
   unsigned int soltag;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   soltag = ++(conshdlrdata->lastsoltag);
   maxviol = 0.0;
   maypropfeasible = conshdlrdata->trysolheur != NULL && SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED
      && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;

   /* check nonlinear constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
      {
         *result = SCIP_INFEASIBLE;
         maxviol = MAX(maxviol, getConsAbsViolation(conss[c]));  /*lint !e666*/

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* print reason for infeasibility */
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhsviol);
            }
            if( consdata->rhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", consdata->rhsviol);
            }
         }
         else if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible && !completely )
         {
            /* if we don't want to pass to subnlp heuristic and don't need to print reasons, then can stop checking here */
            return SCIP_OKAY;
         }

         /* do not try to shift linear variables if violation is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, getConsAbsViolation(conss[c])) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            /* update information on linear variables that may be in- or decreased, if initsolve has not done so yet */
            if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
               consdataFindUnlockedLinearVar(scip, conshdlr, consdata);

            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef > 0.0) &&
                  !(consdata->linvardecr != NULL && consdata->linvardecrcoef < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(consdata->rhsviol > SCIPfeastol(scip));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef < 0.0) &&
                  !(consdata->linvardecr != NULL && consdata->linvardecrcoef > 0.0) )
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

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL && !SCIPisInfinity(scip, maxviol) )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   int nchgbds = 0;
   int ndelconss = 0;

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, SCIPgetDepth(scip) == 0, result, &nchgbds, &ndelconss) );
   assert(nchgbds >= 0);

   /* TODO would it make sense to check for redundant constraints? */

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int c;

   *result = SCIP_DIDNOTFIND;

   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* simplify constraints and replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, presoltiming, &infeasible, ndelconss, naddconss, nchgcoefs) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* merge constraints with the same root expression */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      SCIP_Bool success;

      SCIP_CALL( presolMergeConss(scip, conss, nconss, &success) );
      if( success )
         *result = SCIP_SUCCESS;
   }

   /* propagate constraints */
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, (presoltiming & (SCIP_PRESOLTIMING_MEDIUM | SCIP_PRESOLTIMING_EXHAUSTIVE)) != 0, result, nchgbds, ndelconss) );
   if( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   /* check for redundant constraints, remove constraints that are a value expression */
   SCIP_CALL( checkRedundancyConss(scip, conshdlr, conss, nconss, &infeasible, ndelconss, nchgbds) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* try to upgrade constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool upgraded;

      /* skip inactive and deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsActive(conss[c]) )
         continue;

      SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );  /*lint !e794*/
   }

   /* fix variables that are contained in only one expression constraint to their upper or lower bounds, if possible */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && SCIPisPresolveFinished(scip)
      && !conshdlrdata->checkedvarlocks && conshdlrdata->checkvarlocks != 'd' )
   {
      /* run this presolving technique only once because we don't want to generate identical bound disjunction
       * constraints multiple times
       */
      conshdlrdata->checkedvarlocks = TRUE;

      for( c = 0; c < nconss; ++c )
      {
         int tmpnchgvartypes = 0;
         int tmpnaddconss = 0;

         SCIP_CALL( presolSingleLockedVars(scip, conshdlr, conss[c], &tmpnchgvartypes, &tmpnaddconss, &infeasible) );
         SCIPdebugMsg(scip, "presolSingleLockedVars() for %s: nchgvartypes=%d naddconss=%d infeas=%u\n",
            SCIPconsGetName(conss[c]), tmpnchgvartypes, tmpnaddconss, infeasible);

         if( infeasible )
         {
            SCIPdebugMsg(scip, "presolSingleLockedVars() detected infeasibility\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         (*nchgvartypes) += tmpnchgvartypes;
         (*naddconss) += tmpnaddconss;
      }
   }

   if( *ndelconss > 0 || *nchgbds > 0 || *nupgdconss > 0 || *naddconss > 0 || *nchgvartypes > 0 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExpr NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   /* add locks */
   SCIP_CALL( addLocks(scip, cons, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveExpr)
{  /*lint --e{715}*/

   /* store variable expressions */
   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(cons)) );
   }

   /* simplify root expression if the constraint has been added after presolving */
   if( SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( !consdata->issimplified )
      {
         SCIP_CONSEXPR_EXPR* simplified;
         SCIP_Bool infeasible;
         SCIP_Bool changed;

         /* simplify constraint */
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, consdata->expr, &simplified, &changed, &infeasible) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );
         assert(simplified != NULL);
         consdata->expr = simplified;
         consdata->issimplified = TRUE;
      }
   }

   /* add manually locks to constraints that are not checked for feasibility */
   if( !SCIPconsIsChecked(cons) )
   {
      assert(SCIPconsGetData(cons)->nlockspos == 0);
      assert(SCIPconsGetData(cons)->nlocksneg == 0);

      SCIP_CALL( addLocks(scip, cons, 1, 0) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(cons)) );
   }

#ifdef SCIP_DISABLED_CODE   /* we probably don't need this (?) */
   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_CONSEXPR_ITERATOR* it;

      consdata = SCIPconsGetData(cons);
      SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

      for( expr = SCIPexpriteratorRestartDFS(it, consdata->expr); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIPdebugMsg(scip, "consdeactivate: free nonlinear handler data for expression %p\n", (void*)expr);

         /* remove nonlinear handlers in expression and their data; keep auxiliary variable */
         SCIP_CALL( freeEnfoData(scip, conshdlr, expr, FALSE) );
      }

      SCIPexpriteratorFree(&it);
   }
#endif

   /* remove locks that have been added in consActiveExpr() */
   if( !SCIPconsIsChecked(cons) )
   {
      SCIP_CALL( addLocks(scip, cons, -1, 0) );

      assert(SCIPconsGetData(cons)->nlockspos == 0);
      assert(SCIPconsGetData(cons)->nlocksneg == 0);
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsExpr NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExpr)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print expression */
   if( consdata->expr != NULL )
   {
      SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, consdata->expr, file) );
   }
   else
   {
      SCIPinfoMessage(scip, file, "0");
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


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyExpr)
{  /*lint --e{715}*/
   COPY_MAPVAR_DATA mapvardata;
   SCIP_CONSEXPR_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   assert(cons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   mapvardata.varmap = varmap;
   mapvardata.consmap = consmap;
   mapvardata.global = global;
   mapvardata.valid = TRUE; /* hope the best */

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( copyExpr(sourcescip, scip, sourceconshdlr, sourcedata->expr, &targetexpr, copyVar, &mapvardata) );

   if( targetexpr == NULL )
   {
      *cons = NULL;
      *valid = FALSE;

      return SCIP_OKAY;
   }

   /* validity depends only on the SCIPgetVarCopy() returns from copyVar, which are accumulated in mapvardata.valid */
   *valid = mapvardata.valid;

   /* create copy (captures targetexpr) */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name != NULL ? name : SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   SCIP_CONSEXPR_EXPR* consexprtree;

   SCIPdebugMsg(scip, "cons_expr::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = FALSE;

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, (char**)&endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the beginning of the expression and not a left-hand-side */
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

   debugParse("str should start at beginning of expr: %s\n", str); /*lint !e506 !e681*/

   /* parse expression: so far we did not allocate memory, so can just return in case of readerror */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, str, &str, &consexprtree) );

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
         return SCIP_OKAY;
      }
      *success = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, (char**)&endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
            return SCIP_OKAY;
      }
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name,
      consexprtree, lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   assert(*cons != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );

   debugParse("created expression constraint: <%s>\n", SCIPconsGetName(*cons)); /*lint !e506 !e681*/

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   /* check whether array is too small in order to store all variables */
   if( varssize < consdata->nvarexprs )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      vars[i] = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
      assert(vars[i] != NULL);
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   *nvars = consdata->nvarexprs;
   *success = TRUE;

   return SCIP_OKAY;
}


/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0 /* TODO? */
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsExpr NULL
#endif


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputExpr)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* print statistics for expression handlers */
   printExprHdlrStatistics(scip, conshdlr, file);

   /* print statistics for nonlinear handlers */
   printNlhdlrStatistics(scip, conshdlr, file);

   /* print statistics for constraint handler */
   printConshdlrStatistics(scip, conshdlr, file);

   return SCIP_OKAY;
}

/** creates the handler for an expression handler and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR**    exprhdlr,     /**< buffer where to store expression handler */
   const char*                 name,         /**< name of expression handler (must not be NULL) */
   const char*                 desc,         /**< description of expression handler (can be NULL) */
   unsigned int                precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval)),     /**< point evaluation callback (can not be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA* data          /**< data of expression handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(name != NULL);
   assert(exprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, exprhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*exprhdlr)->precedence = precedence;
   (*exprhdlr)->eval = eval;
   (*exprhdlr)->data = data;

   /* create clocks */
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->estimatetime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->proptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->intevaltime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->simplifytime) );

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize, conshdlrdata->nexprhdlrs+1);

   conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs] = *exprhdlr;
   ++conshdlrdata->nexprhdlrs;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeHdlr(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr)), /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr))  /**< handler free callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copyhdlr = copyhdlr;
   exprhdlr->freehdlr = freehdlr;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free expression data */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata)), /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata))  /**< expression data free callback (can be NULL if data does not need to be freed) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;

   return SCIP_OKAY;
}

/** set the print callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print))    /**< print callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->print = print;

   return SCIP_OKAY;
}

/** set the parse callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrParse(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse))    /**< parse callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->parse = parse;

   return SCIP_OKAY;
}

/** set the curvature detection callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCurvature(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCURVATURE((*curvature)) /**< curvature detection callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->curvature = curvature;

   return SCIP_OKAY;
}

/** set the monotonicity detection callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrMonotonicity(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->monotonicity = monotonicity;

   return SCIP_OKAY;
}

/** set the integrality detection callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrIntegrality(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
   )
{ /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->integrality = integrality;

   return SCIP_OKAY;
}

/** set the hash callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrHash(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRHASH((*hash))      /**< hash callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->hash = hash;

   return SCIP_OKAY;
}

/** set the compare callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCompare(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOMPARE((*compare))/**< compare callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->compare = compare;

   return SCIP_OKAY;
}

/** set the derivative evaluation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrBwdiff(
            SCIP*                      scip,          /**< SCIP data structure */
            SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
            SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
            SCIP_DECL_CONSEXPR_EXPRBWDIFF((*bwdiff))  /**< derivative evaluation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->bwdiff = bwdiff;

   return SCIP_OKAY;
}

/** set the interval evaluation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrIntEval(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval))/**< interval evaluation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;

   return SCIP_OKAY;
}

/** set the simplify callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrSimplify(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRSIMPLIFY((*simplify))  /**< simplify callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->simplify = simplify;

   return SCIP_OKAY;
}

/** set the reverse propagation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrReverseProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRREVERSEPROP((*reverseprop))/**< reverse propagation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->reverseprop = reverseprop;

   return SCIP_OKAY;
}

/** set the separation and estimation callbacks of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPREXITSEPA((*exitsepa)), /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRESTIMATE((*estimate))  /**< estimator callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->initsepa = initsepa;
   exprhdlr->exitsepa = exitsepa;
   exprhdlr->estimate = estimate;

   return SCIP_OKAY;
}

/** gives expression handlers */
SCIP_CONSEXPR_EXPRHDLR** SCIPgetConsExprExprHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprhdlrs;
}

/** gives number of expression handlers */
int SCIPgetConsExprExprNHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->nexprhdlrs;
}

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_CONSEXPR_EXPRHDLR* SCIPfindConsExprExprHdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of expression handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nexprhdlrs; ++h )
      if( strcmp(SCIPgetConsExprExprHdlrName(conshdlrdata->exprhdlrs[h]), name) == 0 )
         return conshdlrdata->exprhdlrs[h];

   return NULL;
}

/** returns expression handler for variable expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrVar(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvarhdlr;
}

/** returns expression handler for constant value expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrValue(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvalhdlr;
}

/** returns expression handler for sum expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSum(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprsumhdlr;
}

/** returns expression handler for product expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrProduct(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprprodhdlr;
}

/** returns expression handler for power expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprpowhdlr;
}

/** returns expression handler for signed power expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSignPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprsignpowhdlr;
}

/** returns expression handler for exponential expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrExponential(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprexphdlr;
}

/** gives the name of an expression handler */
const char* SCIPgetConsExprExprHdlrName(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->name;
}

/** gives the description of an expression handler (can be NULL) */
const char* SCIPgetConsExprExprHdlrDescription(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->desc;
}

/** gives the precedence of an expression handler */
unsigned int SCIPgetConsExprExprHdlrPrecedence(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->precedence;
}

/** gives the data of an expression handler */
SCIP_CONSEXPR_EXPRHDLRDATA* SCIPgetConsExprExprHdlrData(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr      /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->data;
}
/** returns whether expression handler implements the print callback */
SCIP_Bool SCIPhasConsExprExprHdlrPrint(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->print != NULL;
}

/** returns whether expression handler implements the backward differentiation callback */
SCIP_Bool SCIPhasConsExprExprHdlrBwdiff(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->bwdiff != NULL;
}

/** returns whether expression handler implements the interval evaluation callback */
SCIP_Bool SCIPhasConsExprExprHdlrIntEval(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->inteval != NULL;
}

/** returns whether expression handler implements the estimator callback */
SCIP_Bool SCIPhasConsExprExprHdlrEstimate(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->estimate != NULL;
}

/** returns whether expression handler implements the simplification callback */
SCIP_Bool SCIPhasConsExprExprHdlrSimplify(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->simplify != NULL;
}

/** returns whether expression handler implements the curvature callback */
SCIP_Bool SCIPhasConsExprExprHdlrCurvature(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->curvature != NULL;
}

/** returns whether expression handler implements the reverse propagation callback */
SCIP_Bool SCIPhasConsExprExprHdlrReverseProp(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->reverseprop != NULL;
}

/** returns whether expression handler implements the initialization callback */
SCIP_Bool SCIPhasConsExprExprHdlrInitSepa(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->initsepa != NULL;
}

/** returns whether expression handler implements the deinitialization callback */
SCIP_Bool SCIPhasConsExprExprHdlrExitSepa(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->exitsepa != NULL;
}

/** calls the print callback of an expression handler */
SCIP_DECL_CONSEXPR_EXPRPRINT(SCIPprintConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);

   if( SCIPhasConsExprExprHdlrPrint(expr->exprhdlr) )
   {
      SCIP_CALL( (*expr->exprhdlr->print)(scip, expr, stage, currentchild, parentprecedence, file) );
   }
   else
   {
      /* default: <hdlrname>(<child1>, <child2>, ...) */
      switch( stage )
      {
         case SCIP_CONSEXPRITERATOR_ENTEREXPR :
         {
            SCIPinfoMessage(scip, file, SCIPgetConsExprExprHdlrName(expr->exprhdlr));
            if( SCIPgetConsExprExprNChildren(expr) > 0 )
            {
               SCIPinfoMessage(scip, file, "(");
            }
            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITEDCHILD :
         {
            if( currentchild < SCIPgetConsExprExprNChildren(expr)-1 )
            {
               SCIPinfoMessage(scip, file, ", ");
            }
            else
            {
               SCIPinfoMessage(scip, file, ")");
            }

            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
         case SCIP_CONSEXPRITERATOR_LEAVEEXPR :
         default:
            break;
      }
   }

   return SCIP_OKAY;
}

/** calls the parse callback of an expression handler */
SCIP_DECL_CONSEXPR_EXPRPARSE(SCIPparseConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(exprhdlr != NULL);
   assert(expr != NULL);
   assert(success != NULL);

   *expr = NULL;

   if( exprhdlr->parse == NULL )
   {
      /* TODO we could just look for a comma separated list of operands and try to initialize the expr with this one?
       * That would be sufficient for sin, cos, exp, log, abs, for example.
       */
      SCIPdebugMessage("Expression handler <%s> has no parsing method.\n", SCIPgetConsExprExprHdlrName(exprhdlr));
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* give control to exprhdlr's parser */
   SCIP_CALL( exprhdlr->parse(scip, consexprhdlr, exprhdlr, string, endstring, expr, success) );

   assert(*success || (*expr == NULL));

   return SCIP_OKAY;
}

/** calls the expression hash callback */
SCIP_DECL_CONSEXPR_EXPRHASH(SCIPhashConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashkey != NULL);

   if( expr->exprhdlr->hash != NULL )
   {
      SCIP_CALL( (*expr->exprhdlr->hash)(scip, expr, hashkey, childrenhashes) );
   }
   else
   {
      int i;

      /* compute initial hash from expression handler name if callback is not implemented
       * this can lead to more collisions and thus a larger number of expensive expression compare calls
       */
      *hashkey = 0;
      for( i = 0; expr->exprhdlr->name[i] != '\0'; i++ )
         *hashkey += (unsigned int) expr->exprhdlr->name[i]; /*lint !e571*/

      *hashkey = SCIPcalcFibHash((SCIP_Real)*hashkey);

      /* now make use of the hashkeys of the children */
      for( i = 0; i < expr->nchildren; ++i )
         *hashkey ^= childrenhashes[i];
   }

   return SCIP_OKAY;
}

/** calls the expression compare callback */
SCIP_DECL_CONSEXPR_EXPRCOMPARE(SCIPcompareConsExprExprHdlr)
{
   assert(expr1 != NULL);
   assert(expr2 != NULL);
   assert(expr1->exprhdlr == expr2->exprhdlr);

   if( expr1->exprhdlr->compare != NULL )
   {
      /* enforces OR1-OR4 */
      return expr1->exprhdlr->compare(expr1, expr2);
   }
   else
   {
      /* enforces OR5: default comparison method of expressions of the same type:
       * expr1 < expr2 if and only if expr1_i = expr2_i for all i < k and expr1_k < expr2_k.
       * if there is no such k, use number of children to decide
       * if number of children is equal, both expressions are equal
       * @note: Warning, this method doesn't know about expression data. So if your expressions have special data,
       * you must implement the compare callback: SCIP_DECL_CONSEXPR_EXPRCMP
       */
      int i;
      int nchildren1;
      int nchildren2;
      int compareresult;

      nchildren1 = SCIPgetConsExprExprNChildren(expr1);
      nchildren2 = SCIPgetConsExprExprNChildren(expr2);

      for( i = 0; i < nchildren1 && i < nchildren2; ++i )
      {
         compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[i]);
         if( compareresult != 0 )
            return compareresult;
      }

      return nchildren1 == nchildren2 ? 0 : nchildren1 < nchildren2 ? -1 : 1;
   }
}

/** calls the backward-differentiation callback of an expression handler
 *
 * further, allows to different w.r.t. given expression and children values
 */
SCIP_RETCODE SCIPbwdiffConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*        expr,         /**< expression */
   int                        childidx,     /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*                 derivative,   /**< buffer to store value of derivative */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real                  exprval       /**< value for expression, used only if childrenvals is not NULL */
)
{
   SCIP_Real* origchildrenvals;
   SCIP_Real origexprval;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);
   assert(derivative != NULL);

   if( expr->exprhdlr->bwdiff == NULL )
   {
      *derivative = SCIP_INVALID;
      return SCIP_OKAY;
   }

   /* temporarily overwrite the evalvalue in all children and expr with values from childrenvals and exprval, resp. */
   if( childrenvals != NULL )
   {
      if( expr->nchildren > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &origchildrenvals, expr->nchildren) );

         for( c = 0; c < expr->nchildren; ++c )
         {
            origchildrenvals[c] = expr->children[c]->evalvalue;
            expr->children[c]->evalvalue = childrenvals[c];
         }
      }

      origexprval = expr->evalvalue;
      expr->evalvalue = exprval;
   }

   SCIP_CALL( expr->exprhdlr->bwdiff(scip, expr, childidx, derivative) );

   /* restore original evalvalues in children */
   if( childrenvals != NULL )
   {
      if( expr->nchildren > 0 )
      {
         for( c = 0; c < expr->nchildren; ++c )
            expr->children[c]->evalvalue = origchildrenvals[c];  /*lint !e644*/

         SCIPfreeBufferArray(scip, &origchildrenvals);
      }

      expr->evalvalue = origexprval;   /*lint !e644*/
   }

   return SCIP_OKAY;
}

/** calls the evaluation callback of an expression handler
 *
 * further, allows to evaluate w.r.t. given children values
 */
SCIP_RETCODE SCIPevalConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*        expr,         /**< expression */
   SCIP_Real*                 val,          /**< buffer store value of expression */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*                  sol           /**< solution that is evaluated (used by the var-expression) */
)
{
   SCIP_Real* origvals = NULL;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);
   assert(expr->exprhdlr->eval != NULL);
   assert(val != NULL);

   /* temporarily overwrite the evalvalue in all children with values from childrenvals */
   if( childrenvals != NULL && expr->nchildren > 0 )
   {
      int c;

      SCIP_CALL( SCIPallocBufferArray(scip, &origvals, expr->nchildren) );

      for( c = 0; c < expr->nchildren; ++c )
      {
         origvals[c] = expr->children[c]->evalvalue;
         expr->children[c]->evalvalue = childrenvals[c];
      }
   }

   /* call expression eval callback */
   SCIP_CALL( expr->exprhdlr->eval(scip, expr, val, sol) );

   /* restore original evalvalues in children */
   if( origvals != NULL )
   {
      int c;
      for( c = 0; c < expr->nchildren; ++c )
         expr->children[c]->evalvalue = origvals[c];

      SCIPfreeBufferArray(scip, &origvals);
   }

   return SCIP_OKAY;
}

/** calls the expression interval evaluation callback */
SCIP_DECL_CONSEXPR_EXPRINTEVAL(SCIPintevalConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(interval != NULL);

   if( SCIPhasConsExprExprHdlrIntEval(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->intevaltime) );
      SCIP_CALL( expr->exprhdlr->inteval(scip, expr, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->intevaltime) );

      ++expr->exprhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls estimator method of expression handler */
SCIP_DECL_CONSEXPR_EXPRESTIMATE(SCIPestimateConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(coefs != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPhasConsExprExprHdlrEstimate(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->estimatetime) );
      SCIP_CALL( expr->exprhdlr->estimate(scip, conshdlr, expr, sol, overestimate, targetvalue, coefs, constant, islocal, success, branchcand) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->estimatetime) );

      /* update statistics */
      ++expr->exprhdlr->nestimatecalls;
   }

   return SCIP_OKAY;
}

/** calls the simplification method of an expression handler */
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(SCIPsimplifyConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);

   if( SCIPhasConsExprExprHdlrSimplify(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->simplifytime) );
      SCIP_CALL( expr->exprhdlr->simplify(scip, conshdlr, expr, simplifiedexpr) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->simplifytime) );

      /* update statistics */
      ++(expr->exprhdlr->nsimplifycalls);
      if( expr != *simplifiedexpr )
         ++(expr->exprhdlr->nsimplified);
   }

   return SCIP_OKAY;
}

/** calls the curvature check method of an expression handler */
SCIP_DECL_CONSEXPR_EXPRCURVATURE(SCIPcurvatureConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPhasConsExprExprHdlrCurvature(expr->exprhdlr) )
   {
      SCIP_CALL( expr->exprhdlr->curvature(scip, conshdlr, expr, exprcurvature, success, childcurv) );
   }

   return SCIP_OKAY;
}


/** calls the expression callback for reverse propagation */
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(SCIPreversepropConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   *infeasible = FALSE;
   *nreductions = 0;

   if( SCIPhasConsExprExprHdlrReverseProp(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->proptime) );
      SCIP_CALL( expr->exprhdlr->reverseprop(scip, conshdlr, expr, reversepropqueue, infeasible, nreductions, force) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->proptime) );

      /* update statistics */
      expr->exprhdlr->ndomreds += *nreductions;
      if( *infeasible )
         ++(expr->exprhdlr->ncutoffs);
      ++(expr->exprhdlr->npropcalls);
   }

   return SCIP_OKAY;
}

/** calls the separation initialization method of an expression handler */
SCIP_DECL_CONSEXPR_EXPRINITSEPA(SCIPinitsepaConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   if( SCIPhasConsExprExprHdlrInitSepa(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->estimatetime) );
      SCIP_CALL( expr->exprhdlr->initsepa(scip, conshdlr, cons, expr, overestimate, underestimate, infeasible) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->estimatetime) );

      /* update statistics */
      if( *infeasible )
         ++(expr->exprhdlr->ncutoffs);
      ++(expr->exprhdlr->nestimatecalls);
   }

   return SCIP_OKAY;
}

/** calls the separation deinitialization method of an expression handler */
SCIP_DECL_CONSEXPR_EXPREXITSEPA(SCIPexitsepaConsExprExprHdlr)
{
   assert(scip != NULL);
   assert(expr != NULL);

   if( SCIPhasConsExprExprHdlrExitSepa(expr->exprhdlr) )
   {
      SCIP_CALL( SCIPstartClock(scip, expr->exprhdlr->estimatetime) );
      SCIP_CALL( expr->exprhdlr->exitsepa(scip, expr) );
      SCIP_CALL( SCIPstopClock(scip, expr->exprhdlr->estimatetime) );
   }

   return SCIP_OKAY;
}

/** increments the branching score count of an expression handler */
void SCIPincrementConsExprExprHdlrNBranchScore(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr
   )
{
   assert(exprhdlr != NULL);

   ++exprhdlr->nbranchscores;
}

/** returns whether we are ok to branch on auxiliary variables
 *
 * Currently returns whether depth of node in B&B tree is at least value of constraints/expr/branching/aux parameter.
 */
SCIP_Bool SCIPgetConsExprBranchAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->branchauxmindepth <= SCIPgetDepth(scip);
}

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPcreateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_CONSEXPR_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   )
{
   SCIP_CALL( createExpr(scip, expr, exprhdlr, exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** creates and captures an expression with up to two children */
SCIP_RETCODE SCIPcreateConsExprExpr2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data */
   SCIP_CONSEXPR_EXPR*     child1,           /**< first child (can be NULL) */
   SCIP_CONSEXPR_EXPR*     child2            /**< second child (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(exprhdlr != NULL);

   if( child1 != NULL && child2 != NULL )
   {
      SCIP_CONSEXPR_EXPR* pair[2];
      pair[0] = child1;
      pair[1] = child2;

      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 2, pair) );
   }
   else if( child2 == NULL )
   {
      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, child1 == NULL ? 0 : 1, &child1) );
   }
   else
   {
      /* child2 != NULL, child1 == NULL */
      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 1, &child2) );
   }

   return SCIP_OKAY;
}

/** creates and captures an expression from a node in an (old-style) expression graph */
SCIP_RETCODE SCIPcreateConsExprExpr3(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRGRAPH*         exprgraph,        /**< expression graph */
   SCIP_EXPRGRAPHNODE*     node              /**< expression graph node */
   )
{
   SCIP_CONSEXPR_EXPR** children = NULL;
   int nchildren;
   int c = 0;

   assert(expr != NULL);
   assert(node != NULL);

   *expr = NULL;
   nchildren = SCIPexprgraphGetNodeNChildren(node);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &children[c], exprgraph, SCIPexprgraphGetNodeChildren(node)[c]) );
         if( children[c] == NULL )
            goto TERMINATE;
      }

   }

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
      case SCIP_EXPR_CONST :
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, consexprhdlr, expr, SCIPexprgraphGetNodeOperatorReal(node)) );
         break;

      case SCIP_EXPR_VARIDX :
      {
         int varidx;

         varidx = SCIPexprgraphGetNodeOperatorIndex(node);
         assert(varidx >= 0);
         assert(varidx < SCIPexprgraphGetNVars(exprgraph));

         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, expr, (SCIP_VAR*)SCIPexprgraphGetVars(exprgraph)[varidx]) );

         break;
      }

      case SCIP_EXPR_PLUS:
      {
         assert(nchildren == 2);
         assert(children != NULL && children[0] != NULL && children[1] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_MINUS:
      {
         SCIP_Real coefs[2] = {1.0, -1.0};

         assert(nchildren == 2);
         assert(children != NULL && children[0] != NULL && children[1] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, coefs, 0.0) );

         break;
      }

      case SCIP_EXPR_MUL:
      {
         assert(nchildren == 2);
         assert(children != NULL && children[0] != NULL && children[1] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, children, 1.0) );

         break;
      }

      case SCIP_EXPR_DIV:
      {
         SCIP_CONSEXPR_EXPR* factors[2];

         assert(nchildren == 2);
         assert(children != NULL && children[0] != NULL && children[1] != NULL);

         factors[0] = children[0];
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &factors[1], children[1], -1.0) );
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, factors, 1.0) );

         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factors[1]) );

         break;
      }

      case SCIP_EXPR_SQUARE:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, expr, *children, 2.0) );

         break;
      }

      case SCIP_EXPR_SQRT:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, expr, *children, 0.5) );

         break;
      }

      case SCIP_EXPR_REALPOWER:
      {
         SCIP_Real exponent;

         exponent = SCIPexprgraphGetNodeRealPowerExponent(node);

         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, expr, *children, exponent) );

         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         SCIP_Real exponent;

         exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);

         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, expr, *children, exponent) );

         break;
      }

      case SCIP_EXPR_SIGNPOWER:
      {
         SCIP_Real exponent;

         exponent = (SCIP_Real)SCIPexprgraphGetNodeSignPowerExponent(node);

         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprSignPower(scip, consexprhdlr, expr, *children, exponent) );

         break;
      }

      case SCIP_EXPR_SUM:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_PRODUCT:
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, nchildren, children, 1.0) );

         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, SCIPexprgraphGetNodeLinearCoefs(node), SCIPexprgraphGetNodeLinearConstant(node)) );

         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_QUADELEM quadelem;
         SCIP_CONSEXPR_EXPR* prod;
         int i;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodeQuadraticConstant(node)) );

         /* append linear terms */
         if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node) != NULL )
         {
            for( i = 0; i < nchildren; ++i )
            {
               if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i] != 0.0 )
               {
                  assert(children != NULL);
                  SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[i], SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i]) );
               }
            }
         }

         /* append quadratic terms */
         for( i = 0; i < SCIPexprgraphGetNodeQuadraticNQuadElements(node); ++i )
         {
            quadelem = SCIPexprgraphGetNodeQuadraticQuadElements(node)[i];

            if( quadelem.idx1 == quadelem.idx2 )
            {
               assert(children != NULL);
               SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &prod, children[quadelem.idx1], 2.0) );
            }
            else
            {
               SCIP_CONSEXPR_EXPR* prodchildren[2];

               assert(children != NULL);

               prodchildren[0] = children[quadelem.idx1];
               prodchildren[1] = children[quadelem.idx2];

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 2, prodchildren, 1.0) );
            }

            SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, prod, quadelem.coef) );

            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prod) );
         }

         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPRDATA_MONOMIAL* monom;
         int m;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodePolynomialConstant(node)) );

         /* append monomials */
         for( m = 0; m < SCIPexprgraphGetNodePolynomialNMonomials(node); ++m )
         {
            SCIP_Real* exponents;

            monom = SCIPexprgraphGetNodePolynomialMonomials(node)[m];
            exponents = SCIPexprGetMonomialExponents(monom);

            if( SCIPexprGetMonomialNFactors(monom) == 1 && (exponents == NULL || exponents[0] == 1.0) )
            {
               assert(children != NULL && children[SCIPexprGetMonomialChildIndices(monom)[0]] != NULL);

               /* monom is linear in child -> append child itself */
               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[SCIPexprGetMonomialChildIndices(monom)[0]], SCIPexprGetMonomialCoef(monom)) );
            }
            else
            {
               /* monom is nonlinear -> translate into a product expression */
               SCIP_CONSEXPR_EXPR* monomial;
               int f;

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &monomial, 0, NULL, 1.0) );

               for( f = 0; f < SCIPexprGetMonomialNFactors(monom); ++f )
               {
                  assert(children != NULL && children[SCIPexprGetMonomialChildIndices(monom)[f]] != NULL);
                  if( exponents == NULL || exponents[f] == 1.0 )
                  {
                     SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, monomial, children[SCIPexprGetMonomialChildIndices(monom)[f]]) );
                  }
                  else
                  {
                     SCIP_CONSEXPR_EXPR* powexpr;

                     SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr, children[SCIPexprGetMonomialChildIndices(monom)[f]], exponents[f]) );
                     SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, monomial, powexpr) );
                     SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr) );
                  }
               }

               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, monomial, SCIPexprGetMonomialCoef(monom)) );
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &monomial) );
            }
         }

         break;
      }

      case SCIP_EXPR_EXP:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprExp(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_LOG:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprLog(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_ABS:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprAbs(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_SIN:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprSin(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_COS:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprCos(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_TAN:
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_SIGN:
      case SCIP_EXPR_USER:
      case SCIP_EXPR_PARAM:
      case SCIP_EXPR_LAST:
      default:
         goto TERMINATE;
   }


TERMINATE:
   /* release all created children expressions (c-1...0) */
   for( --c; c >= 0; --c )
   {
      assert(children != NULL && children[c] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[c]) );
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** creates and captures an expression representing a monomial */
SCIP_RETCODE SCIPcreateConsExprExpr4(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   int                     nfactors,         /**< number of factors in monomial */
   SCIP_VAR**              vars,             /**< variables in the monomial */
   SCIP_Real*              exponents         /**< exponent in each factor, or NULL if all 1.0 */
   )
{
   assert(consexprhdlr != NULL);
   assert(expr != NULL);
   assert(nfactors >= 0);

   /* return 1 as constant expression if there are no factors */
   if( nfactors == 0 )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, consexprhdlr, expr, 1.0) );
   }
   else if( nfactors == 1 )
   {
      /* only one factor and exponent is 1 => return factors[0] */
      if( exponents == NULL || exponents[0] == 1.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, expr, vars[0]) );
      }
      else
      {
         SCIP_CONSEXPR_EXPR* varexpr;

         /* create variable and power expression */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, vars[0]) );
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, expr, varexpr, exponents[0]) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
      }
   }
   else
   {
      SCIP_CONSEXPR_EXPR** children;
      int i;

      /* allocate memory to store the children */
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nfactors) );

      /* create children */
      for( i = 0; i < nfactors; ++i )
      {
         /* check whether to create a power expression or not, i.e., exponent == 1 */
         if( exponents == NULL || exponents[i] == 1.0 )
         {
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &children[i], vars[i]) );
         }
         else
         {
            SCIP_CONSEXPR_EXPR* varexpr;

            /* create variable and pow expression */
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, vars[i]) );
            SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &children[i], varexpr, exponents[i]) );
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
         }
      }

      /* create product expression */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, nfactors, children, 1.0) );

      /* release children */
      for( i = 0; i < nfactors; ++i )
      {
         assert(children[i] != NULL);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[i]) );
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &children);
   }

   return SCIP_OKAY;
}

/** appends child to the children list of expr */
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(child != NULL);
   assert(expr->monotonicitysize == 0);  /* should not append child while mononoticity is stored in expr (not updated here) */
   assert(expr->nlocksneg == 0);  /* should not append child while expression is locked (not updated here) */
   assert(expr->nlockspos == 0);  /* should not append child while expression is locked (not updated here) */

   ENSUREBLOCKMEMORYARRAYSIZE(scip, expr->children, expr->childrensize, expr->nchildren + 1);

   expr->children[expr->nchildren] = child;
   ++expr->nchildren;

   /* capture child */
   SCIPcaptureConsExprExpr(child);

   return SCIP_OKAY;
}

/** remove all children of expr
 *
 * only use if you really know what you are doing
 */
SCIP_RETCODE SCIPremoveConsExprExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   for( c = 0; c < expr->nchildren; ++c )
   {
      assert(expr->children[c] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(expr->children[c])) );
   }

   expr->nchildren = 0;

   return SCIP_OKAY;
}

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_RETCODE SCIPreplaceConsExprExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_CONSEXPR_EXPR*     newchild          /**< the new child */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(newchild != NULL);
   assert(childidx < SCIPgetConsExprExprNChildren(expr));
   assert(expr->monotonicitysize == 0);  /* should not append child while mononoticity is stored in expr (not updated here) */
   assert(expr->nlocksneg == 0);  /* should not append child while expression is locked (not updated here) */
   assert(expr->nlockspos == 0);  /* should not append child while expression is locked (not updated here) */

   /* do nothing if child is not changing */
   if( newchild == expr->children[childidx] )
      return SCIP_OKAY;

   /* capture new child (do this before releasing the old child in case there are equal */
   SCIPcaptureConsExprExpr(newchild);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(expr->children[childidx])) );
   expr->children[childidx] = newchild;

   return SCIP_OKAY;
}

/** duplicates the given expression */
SCIP_RETCODE SCIPduplicateConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< original expression */
   SCIP_CONSEXPR_EXPR**  copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_Bool             copychildren        /**< whether children (and all successors) should be copied, too */
   )
{
   if( copychildren )
   {
      SCIP_CALL( copyExpr(scip, scip, consexprhdlr, expr, copyexpr, NULL, NULL) );
   }
   else
   {
      /* copy expression data */
      SCIP_CONSEXPR_EXPRDATA* exprdatacopy = NULL;
      if( SCIPgetConsExprExprData(expr) != NULL )
      {
         assert(expr->exprhdlr->copydata != NULL);
         SCIP_CALL( expr->exprhdlr->copydata(scip, expr->exprhdlr, &exprdatacopy, scip, expr, NULL, NULL) );
      }

      /* create expression with same handler and copied data, but without children */
      SCIP_CALL( SCIPcreateConsExprExpr(scip, copyexpr, expr->exprhdlr, exprdatacopy, 0, NULL) );
   }

   assert(*copyexpr != NULL);

   return SCIP_OKAY;
}

/** gets the number of times the expression is currently captured */
int SCIPgetConsExprExprNUses(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nuses;
}

/** captures an expression (increments usage count) */
void SCIPcaptureConsExprExpr(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   ++expr->nuses;
}

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_RETCODE SCIPreleaseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  rootexpr            /**< pointer to expression to be released */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;

   assert(rootexpr != NULL);
   assert(*rootexpr != NULL);
   assert((*rootexpr)->nuses > 0);

   if( (*rootexpr)->nuses > 1 )
   {
      --(*rootexpr)->nuses;
      *rootexpr = NULL;

      return SCIP_OKAY;
   }

   /* handle the root expr separately: free enfodata and expression data here */
   SCIP_CALL( freeEnfoData(scip, NULL, *rootexpr, TRUE) );

   if( (*rootexpr)->exprdata != NULL )
   {
      assert((*rootexpr)->exprhdlr->freedata != NULL);
      SCIP_CALL( (*rootexpr)->exprhdlr->freedata(scip, *rootexpr) );
   }

   SCIP_CALL( SCIPexpriteratorCreate(&it, SCIPfindConshdlr(scip, CONSHDLR_NAME), SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, *rootexpr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_VISITEDCHILD);
   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it) ; )
   {
      /* expression should be used by its parent and maybe by the iterator (only the root!)
       * in VISITEDCHILD we assert that expression is only used by its parent
       */
      assert(expr != NULL);
      assert(0 <= expr->nuses && expr->nuses <= 2);

      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
         {
            /* check whether a child needs to be visited (nuses == 1)
             * if not, then we still have to release it
             */
            SCIP_CONSEXPR_EXPR* child;

            child = SCIPexpriteratorGetChildExprDFS(it);
            if( child->nuses > 1 )
            {
               /* child is not going to be freed: just release it */
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }

            assert(child->nuses == 1);

            /* free child's enfodata and expression data when entering child */
            SCIP_CALL( freeEnfoData(scip, NULL, child, TRUE) );

            if( child->exprdata != NULL )
            {
               assert(child->exprhdlr->freedata != NULL);
               SCIP_CALL( child->exprhdlr->freedata(scip, child) );
               assert(child->exprdata == NULL);
            }

            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITEDCHILD :
         {
            /* free child after visiting it */
            SCIP_CONSEXPR_EXPR* child;

            child = SCIPexpriteratorGetChildExprDFS(it);
            /* child should only be used by its parent */
            assert(child->nuses == 1);

            /* child should have no data associated */
            assert(child->exprdata == NULL);

            /* free child expression */
            SCIP_CALL( freeExpr(scip, &child) );
            expr->children[SCIPexpriteratorGetChildIdxDFS(it)] = NULL;

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

   SCIPexpriteratorFree(&it);

   /* handle the root expr separately: free its children and itself here */
   SCIP_CALL( freeExpr(scip, rootexpr) );

   return SCIP_OKAY;
}

/** gives the number of children of an expression */
int SCIPgetConsExprExprNChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nchildren;
}

/** gives the children of an expression (can be NULL if no children) */
SCIP_CONSEXPR_EXPR** SCIPgetConsExprExprChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->children;
}

/** gets the handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlr(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprhdlr;
}

/** gets the expression data of an expression */
SCIP_CONSEXPR_EXPRDATA* SCIPgetConsExprExprData(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprdata;
}

/** returns whether an expression is a variable expression */
SCIP_Bool SCIPisConsExprExprVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return strcmp(expr->exprhdlr->name, "var") == 0;
}

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_VAR* SCIPgetConsExprExprAuxVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return SCIPisConsExprExprVar(expr) ? SCIPgetConsExprExprVarVar(expr) : expr->auxvar;
}

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
void SCIPsetConsExprExprData(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPRDATA* exprdata          /**< expression data to be set (can be NULL) */
   )
{
   assert(expr != NULL);
   assert(exprdata == NULL || expr->exprhdlr->copydata != NULL);  /* copydata must be available if there is expression data */
   assert(exprdata == NULL || expr->exprhdlr->freedata != NULL);  /* freedata must be available if there is expression data */

   expr->exprdata = exprdata;
}

/** print an expression as info-message */
SCIP_RETCODE SCIPprintConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be printed */
   FILE*                   file              /**< file to print to, or NULL for stdout */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPRITERATOR_STAGE stage;
   int currentchild;
   unsigned int parentprecedence;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, consexprhdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ALLSTAGES);

   while( !SCIPexpriteratorIsEnd(it) )
   {
      assert(expr->exprhdlr != NULL);
      stage = SCIPexpriteratorGetStageDFS(it);

      if( stage == SCIP_CONSEXPRITERATOR_VISITEDCHILD || stage == SCIP_CONSEXPRITERATOR_VISITINGCHILD )
         currentchild = SCIPexpriteratorGetChildIdxDFS(it);
      else
         currentchild = -1;

      if( SCIPexpriteratorGetParentDFS(it) != NULL )
         parentprecedence = SCIPgetConsExprExprHdlrPrecedence(SCIPgetConsExprExprHdlr(SCIPexpriteratorGetParentDFS(it)));
      else
         parentprecedence = 0;

      SCIP_CALL( SCIPprintConsExprExprHdlr(scip, expr, stage, currentchild, parentprecedence, file) );

      expr = SCIPexpriteratorGetNext(it);
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format */
SCIP_RETCODE SCIPprintConsExprExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   )
{
   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(dotdata != NULL);

   if( file == NULL )
      file = stdout;

   SCIP_CALL( SCIPallocBlockMemory(scip, dotdata) );

   (*dotdata)->file = file;
   SCIP_CALL( SCIPexpriteratorCreate(&(*dotdata)->iterator, consexprhdlr, SCIPblkmem(scip)) );
   (*dotdata)->closefile = FALSE;
   (*dotdata)->whattoprint = whattoprint;
   SCIP_CALL( SCIPhashmapCreate(&(*dotdata)->leaveexprs, SCIPblkmem(scip), 100) );

   SCIPinfoMessage(scip, file, "strict digraph exprgraph {\n");
   SCIPinfoMessage(scip, file, "node [fontcolor=white, style=filled, rankdir=LR]\n");

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_RETCODE SCIPprintConsExprExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   )
{
   FILE* f;

   assert(dotdata != NULL);
   assert(filename != NULL);

   f = fopen(filename, "w");
   if( f == NULL )
   {
      SCIPerrorMessage("could not open file <%s> for writing\n", filename);  /* error code would be in errno */
      return SCIP_FILECREATEERROR;
   }

   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, consexprhdlr, dotdata, f, whattoprint) );
   (*dotdata)->closefile = TRUE;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPprintConsExprExprDot(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata,      /**< data as initialized by \ref SCIPprintConsExprExprDotInit() */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   )
{
   SCIP_Real color;
   int c;

   assert(scip != NULL);
   assert(dotdata != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);

   SCIP_CALL( SCIPexpriteratorInit(dotdata->iterator, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   while( !SCIPexpriteratorIsEnd(dotdata->iterator) )
   {
      /* print expression as dot node */

      if( SCIPgetConsExprExprNChildren(expr) == 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(dotdata->leaveexprs, (void*)expr, NULL) );
      }

      /* make up some color from the expression type (it's name) */
      color = 0.0;
      for( c = 0; expr->exprhdlr->name[c] != '\0'; ++c )
         color += (tolower(expr->exprhdlr->name[c]) - 'a') / 26.0;
      color = SCIPfrac(scip, color);
      SCIPinfoMessage(scip, dotdata->file, "n%p [fillcolor=\"%g,%g,%g\", label=\"", expr, color, color, color);

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EXPRHDLR )
      {
         SCIPinfoMessage(scip, dotdata->file, "%s\\n", SCIPgetConsExprExprHdlrName(expr->exprhdlr));
      }

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EXPRSTRING )
      {
         SCIP_CALL( SCIPprintConsExprExprHdlr(scip, expr, SCIP_CONSEXPRITERATOR_ENTEREXPR, -1, 0, dotdata->file) );
         for( c = 0; c < expr->nchildren; ++c )
         {
            SCIP_CALL( SCIPprintConsExprExprHdlr(scip, expr, SCIP_CONSEXPRITERATOR_VISITINGCHILD, c, 0, dotdata->file) );
            SCIPinfoMessage(scip, dotdata->file, "c%d", c);
            SCIP_CALL( SCIPprintConsExprExprHdlr(scip, expr, SCIP_CONSEXPRITERATOR_VISITEDCHILD, c, 0, dotdata->file) );
         }
         SCIP_CALL( SCIPprintConsExprExprHdlr(scip, expr, SCIP_CONSEXPRITERATOR_LEAVEEXPR, -1, 0, dotdata->file) );

         SCIPinfoMessage(scip, dotdata->file, "\\n");
      }

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_NUSES )
      {
         /* print number of uses */
         SCIPinfoMessage(scip, dotdata->file, "%d uses\\n", expr->nuses);
      }

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_NUSES )
      {
         /* print number of locks */
         SCIPinfoMessage(scip, dotdata->file, "%d,%d +,-locks\\n", expr->nlockspos, expr->nlocksneg);
      }

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EVALVALUE )
      {
         /* print eval value */
         SCIPinfoMessage(scip, dotdata->file, "val=%g", expr->evalvalue);

         if( (dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EVALTAG) == SCIP_CONSEXPR_PRINTDOT_EVALTAG )
         {
            /* print also eval tag */
            SCIPinfoMessage(scip, dotdata->file, " (%u)", expr->evaltag);
         }
         SCIPinfoMessage(scip, dotdata->file, "\\n");
      }

      if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_ACTIVITY )
      {
         /* print activity */
         SCIPinfoMessage(scip, dotdata->file, "[%g,%g]", expr->activity.inf, expr->activity.sup);

         if( (dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_ACTIVITYTAG) == SCIP_CONSEXPR_PRINTDOT_ACTIVITYTAG )
         {
            /* print also activity eval tag */
            SCIPinfoMessage(scip, dotdata->file, " (%u)", expr->activitytag);
         }
         SCIPinfoMessage(scip, dotdata->file, "\\n");
      }

      SCIPinfoMessage(scip, dotdata->file, "\"]\n");  /* end of label and end of node */

      /* add edges from expr to its children */
      for( c = 0; c < expr->nchildren; ++c )
         SCIPinfoMessage(scip, dotdata->file, "n%p -> n%p [label=\"c%d\"]\n", (void*)expr, (void*)expr->children[c], c);

      expr = SCIPexpriteratorGetNext(dotdata->iterator);
   }

   return SCIP_OKAY;
}

/** finishes printing of expressions in dot format */
SCIP_RETCODE SCIPprintConsExprExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata      /**< buffer where dot printing data has been stored */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_HASHMAPENTRY* entry;
   FILE* file;
   int i;

   assert(dotdata != NULL);
   assert(*dotdata != NULL);

   file = (*dotdata)->file;
   assert(file != NULL);

   /* iterate through all entries of the map */
   SCIPinfoMessage(scip, file, "{rank=same;");
   for( i = 0; i < SCIPhashmapGetNEntries((*dotdata)->leaveexprs); ++i )
   {
      entry = SCIPhashmapGetEntry((*dotdata)->leaveexprs, i);

      if( entry != NULL )
      {
         expr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetOrigin(entry);
         assert(expr != NULL);
         assert(SCIPgetConsExprExprNChildren(expr) == 0);

         SCIPinfoMessage(scip, file, " n%p", expr);
      }
   }
   SCIPinfoMessage(scip, file, "}\n");

   SCIPinfoMessage(scip, file, "}\n");

   SCIPhashmapFree(&(*dotdata)->leaveexprs);

   SCIPexpriteratorFree(&(*dotdata)->iterator);

   if( (*dotdata)->closefile )
      fclose((*dotdata)->file);

   SCIPfreeBlockMemory(scip, dotdata);

   return SCIP_OKAY;
}

/** shows a single expression by use of dot and gv
 *
 * This function is meant for debugging purposes.
 * It prints the expression into a temporary file in dot format, then calls dot to create a postscript file, then calls ghostview (gv) to show the file.
 * SCIP will hold until ghostscript is closed.
 */
SCIP_RETCODE SCIPshowConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   )
{
   /* this function is for developers, so don't bother with C variants that don't have popen() */
#if _POSIX_C_SOURCE < 2
   SCIPerrorMessage("No POSIX version 2. Try http://distrowatch.com/.");
   return SCIP_ERROR;
#else
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
   FILE* f;

   assert(expr != NULL);

   /* call dot to generate postscript output and show it via ghostview */
   f = popen("dot -Tps | gv --media=a3 -", "w");
   if( f == NULL )
   {
      SCIPerrorMessage("Calling popen() failed");
      return SCIP_FILECREATEERROR;
   }

   /* print all of the expression into the pipe */
   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), &dotdata, f, SCIP_CONSEXPR_PRINTDOT_ALL) );
   SCIP_CALL( SCIPprintConsExprExprDot(scip, dotdata, expr) );
   SCIP_CALL( SCIPprintConsExprExprDotFinal(scip, &dotdata) );

   /* close the pipe */
   (void) pclose(f);

   return SCIP_OKAY;
#endif
}

/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPdismantleConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to dismantle */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   int depth = -1;

   SCIP_CALL( SCIPexpriteratorCreate(&it, SCIPfindConshdlr(scip, CONSHDLR_NAME), SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR | SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   for( ; !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_ENTEREXPR:
         {
            int nspaces;
            const char* type;

            ++depth;
            nspaces = 3 * depth;
            type = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr));

            /* use depth of expression to align output */
            SCIPinfoMessage(scip, NULL, "%*s[%s]: ", nspaces, "", type);

            if( strcmp(type, "var") == 0 )
            {
               SCIP_VAR* var;

               var = SCIPgetConsExprExprVarVar(expr);
               SCIPinfoMessage(scip, NULL, "%s in [%g, %g]", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
                  SCIPvarGetUbLocal(var));
            }
            else if(strcmp(type, "sum") == 0)
               SCIPinfoMessage(scip, NULL, "%g", SCIPgetConsExprExprSumConstant(expr));
            else if(strcmp(type, "prod") == 0)
               SCIPinfoMessage(scip, NULL, "%g", SCIPgetConsExprExprProductCoef(expr));
            else if(strcmp(type, "val") == 0)
               SCIPinfoMessage(scip, NULL, "%g", SCIPgetConsExprExprValueValue(expr));
            else if(strcmp(type, "pow") == 0 || strcmp(type, "signpower") == 0)
               SCIPinfoMessage(scip, NULL, "%g", SCIPgetConsExprExprPowExponent(expr));
            else if(strcmp(type, "exp") == 0)
               SCIPinfoMessage(scip, NULL, "\n");
            else if(strcmp(type, "log") == 0)
               SCIPinfoMessage(scip, NULL, "\n");
            else if(strcmp(type, "abs") == 0)
               SCIPinfoMessage(scip, NULL, "\n");
            else
               SCIPinfoMessage(scip, NULL, "NOT IMPLEMENTED YET\n");

            if(expr->nenfos > 0 )
            {
               int i;
               SCIPinfoMessage(scip, NULL, "   {");

               for( i = 0; i < expr->nenfos - 1; ++i )
                  SCIPinfoMessage(scip, NULL, "%s, ", expr->enfos[i]->nlhdlr->name);

               SCIPinfoMessage(scip, NULL, "%s}", expr->enfos[i]->nlhdlr->name);
            }
            SCIPinfoMessage(scip, NULL, "\n");

            break;
         }

         case SCIP_CONSEXPRITERATOR_VISITINGCHILD:
         {
            int nspaces;
            const char* type;

            nspaces = 3 * depth;
            type = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr));

            if( strcmp(type, "sum") == 0 )
            {
               SCIPinfoMessage(scip, NULL, "%*s   ", nspaces, "");
               SCIPinfoMessage(scip, NULL, "[coef]: %g\n", SCIPgetConsExprExprSumCoefs(expr)[SCIPexpriteratorGetChildIdxDFS(it)]);
            }

            break;
         }

         case SCIP_CONSEXPRITERATOR_LEAVEEXPR:
         {
            --depth;
            break;
         }

         default:
            /* shouldn't be here */
            SCIPABORT();
            break;
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** Creates an expression from a string.
 * We specify the grammar that defines the syntax of an expression. Loosely speaking, a Base will be any "block",
 * a Factor is a Base to a power, a Term is a product of Factors and an Expression is a sum of terms
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * See also @ref parseExpr.
 */
SCIP_RETCODE SCIPparseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos,           /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the expr parsed */
   )
{
   const char* finalpos_;
   SCIP_RETCODE retcode;
   SCIP_HASHMAP* vartoexprvarmap;

   SCIP_CALL( SCIPhashmapCreate(&vartoexprvarmap, SCIPblkmem(scip), 5 * SCIPgetNVars(scip)) );

   /* if parseExpr fails, we still want to free hashmap */
   retcode = parseExpr(scip, consexprhdlr, vartoexprvarmap, exprstr, &finalpos_, expr);

   SCIPhashmapFree(&vartoexprvarmap);

   if( finalpos != NULL )
      *finalpos = finalpos_;

   return retcode;
}

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPgetConsExprExprEvalTag().
 */
SCIP_RETCODE SCIPevalConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(expr != NULL);

   /* if value is up-to-date, then nothing to do */
   if( soltag != 0 && expr->evaltag == soltag )
      return SCIP_OKAY;

   /* assume we'll get a domain error, so we don't have to get this expr back if we abort the iteration
    * if there is no domain error, then we will overwrite the evalvalue in the last leaveexpr stage
    */
   expr->evalvalue = SCIP_INVALID;
   expr->evaltag = soltag;

   SCIP_CALL( SCIPexpriteratorCreate(&it, consexprhdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   while( !SCIPexpriteratorIsEnd(it) )
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
         {
            SCIP_CONSEXPR_EXPR* child;

            if( soltag == 0 )
               break;

            /* check whether child has been evaluated for that solution already */
            child = SCIPexpriteratorGetChildExprDFS(it);
            if( soltag == child->evaltag )
            {
               if( child->evalvalue == SCIP_INVALID ) /*lint !e777*/
                  goto TERMINATE;

               /* skip this child
                * this already returns the next one, so continue with loop
                */
               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_CONSEXPRITERATOR_LEAVEEXPR :
         {
            SCIP_CALL( SCIPevalConsExprExprHdlr(scip, expr, &expr->evalvalue, NULL, sol) );
            expr->evaltag = soltag;

            if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
               goto TERMINATE;

            break;
         }

         default :
            /* we should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

TERMINATE:
   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evalvalue;
}

/** sets the evaluation value */
void SCIPsetConsExprExprEvalValue(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_Real               value,            /**< value to set */
   unsigned int            tag               /**< tag of solution that was evaluated, or 0 */
   )
{
   assert(expr != NULL);

   expr->evalvalue = value;
   expr->evaltag = tag;
}

/** gives the evaluation tag from the last evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evaltag;
}

/** computes the gradient for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprPartialDiff().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_RETCODE SCIPcomputeConsExprExprGradient(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real derivative;
   unsigned int difftag;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(rootexpr != NULL);

   /* ensure expression is evaluated */
   SCIP_CALL( SCIPevalConsExprExpr(scip, consexprhdlr, rootexpr, sol, soltag) );

   /* check if expression could not be evaluated */
   if( SCIPgetConsExprExprValue(rootexpr) == SCIP_INVALID ) /*lint !e777*/
   {
      rootexpr->derivative = SCIP_INVALID;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   if( rootexpr->exprhdlr == conshdlrdata->exprvalhdlr )
   {
      rootexpr->derivative = 0.0;
      return SCIP_OKAY;
   }

   difftag = ++(conshdlrdata->lastdifftag);

   rootexpr->derivative = 1.0;
   rootexpr->difftag = difftag;

   SCIP_CALL( SCIPexpriteratorCreate(&it, consexprhdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, rootexpr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      assert(expr->evalvalue != SCIP_INVALID); /*lint !e777*/

      if( expr->exprhdlr->bwdiff == NULL )
      {
         rootexpr->derivative = SCIP_INVALID;
         break;
      }

      child = SCIPexpriteratorGetChildExprDFS(it);
      assert(child != NULL);

      /* reset the value of the partial derivative w.r.t. a variable expression if we see it for the first time */
      if( child->difftag != difftag && SCIPisConsExprExprVar(child) )
         child->derivative = 0.0;

      /* update differentiation tag of the child */
      child->difftag = difftag;

      /* call backward differentiation callback */
      if( child->exprhdlr == conshdlrdata->exprvalhdlr )
      {
         derivative = 0.0;
      }
      else
      {
         derivative = SCIP_INVALID;
         SCIP_CALL( SCIPbwdiffConsExprExprHdlr(scip, expr, SCIPexpriteratorGetChildIdxDFS(it), &derivative, NULL, 0.0) );

         if( derivative == SCIP_INVALID ) /*lint !e777*/
         {
            rootexpr->derivative = SCIP_INVALID;
            break;
         }
      }

      /* update partial derivative stored in the child expression
       * for a variable, we have to sum up the partial derivatives of the root w.r.t. this variable over all parents
       * for other intermediate expressions, we only store the partial derivative of the root w.r.t. this expression
       */
      if( !SCIPisConsExprExprVar(child) )
         child->derivative = expr->derivative * derivative;
      else
         child->derivative += expr->derivative * derivative;
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error) */
SCIP_Real SCIPgetConsExprExprPartialDiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression which has been used in the last SCIPcomputeConsExprExprGradient() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   )
{
   SCIP_CONSEXPR_EXPR* varexpr;
   SCIP_HASHMAP* var2expr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(consexprhdlr), CONSHDLR_NAME) == 0);
   assert(expr != NULL);
   assert(var != NULL);
   assert(expr->exprhdlr != SCIPgetConsExprExprHdlrValue(consexprhdlr) || expr->derivative == 0.0);

   /* return 0.0 for value expression */
   if( strcmp(expr->exprhdlr->name, "val") == 0 )
      return 0.0;

   /* check if an error occurred during the last SCIPcomputeConsExprExprGradient() call */
   if( expr->derivative == SCIP_INVALID ) /*lint !e777*/
      return SCIP_INVALID;

   /* use variable to expressions mapping which is stored as the expression handler data */
   var2expr = (SCIP_HASHMAP*)SCIPgetConsExprExprHdlrData(SCIPgetConsExprExprHdlrVar(consexprhdlr));
   assert(var2expr != NULL);
   assert(SCIPhashmapExists(var2expr, var));

   varexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(var2expr, var);
   assert(varexpr != NULL);
   assert(SCIPisConsExprExprVar(varexpr));

   /* use difftag to decide whether the variable belongs to the expression */
   return (expr->difftag != varexpr->difftag) ? 0.0 : varexpr->derivative;
}

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error) */
SCIP_Real SCIPgetConsExprExprDerivative(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->derivative;
}

/** returns the activity of the expression
 *
 * The caller needs to make sure that the activity is valid.
 * For expression and nonlinear handlers, this is made sure when the following callbacks are called:
 * - interval evaluation (intervals for children only)
 * - reverse propagation
 * - monotonicity computation
 * - convexity detection
 * - structure detection
 */
SCIP_INTERVAL SCIPgetConsExprExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
#ifndef NDEBUG
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);

   /* check whether activity is valid */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(expr->activitytag >= conshdlrdata->lastboundrelax);
#endif

   return expr->activity;
}

/** returns the tag associated with the activity of the expression
 *
 * Can be compared with SCIPgetConsExprCurBoundsTag() and SCIPgetConsExprLastBoundRelaxTag()
 * to check whether the activity currently stored in this expression is current and valid, respectively.
 */
unsigned int SCIPgetConsExprExprActivityTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   return expr->activitytag;
}

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is not valid (some bound was relaxed since last evaluation).
 * If validsufficient is set to FALSE, then it will also reevaluate activity if a bound tightening was happening
 * since last evaluation.
 */
SCIP_RETCODE SCIPevalConsExprExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler, or NULL */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_INTERVAL*          activity,         /**< interval where to store expression */
   SCIP_Bool               validsufficient   /**< whether any valid activity is sufficient */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(activity != NULL);

   if( consexprhdlr == NULL )
      consexprhdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   if( expr->activitytag < conshdlrdata->lastboundrelax ||
      (!validsufficient && expr->activitytag < conshdlrdata->curboundstag) )
   {
      /* update activity of expression */
      SCIP_CALL( forwardPropExpr(scip, consexprhdlr, expr, FALSE, FALSE, intEvalVarBoundTightening, conshdlrdata, NULL, NULL, NULL) );

      assert(expr->activitytag == conshdlrdata->curboundstag);
   }

   *activity = expr->activity;

   return SCIP_OKAY;
}

/** tightens the activity of an expression and bounds of corresponding (auxiliary) variable (if any)
 *
 *  If a reversepropqueue is given, then the expression will be added to the queue if the tightening is sufficient.
 */
SCIP_RETCODE SCIPtightenConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be tightened */
   SCIP_INTERVAL           newbounds,        /**< new bounds for the expression */
   SCIP_Bool               force,            /**< force tightening of variable bound and add to reversepropqueue even if tightening is below bound strengthening tolerance */
   SCIP_QUEUE*             reversepropqueue, /**< reverse propagation queue, or NULL if not in reverse propagation */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a cutoff was detected */
   int*                    ntightenings      /**< buffer to add the total number of tightenings, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   /* the code below assumes that current activity is valid
    * if it turns out that we cannot ensure that, then we should change code
    */
   assert(expr->activitytag >= SCIPconshdlrGetData(conshdlr)->lastboundrelax || SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, expr->activity));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity));

   *cutoff = FALSE;

#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, "Trying to tighten bounds of expr ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
   SCIPdebugMsgPrint(scip, " from [%.15g,%.15g] to [%.15g,%.15g] (force=%d)\n", expr->activity.inf, expr->activity.sup, newbounds.inf, newbounds.sup, force);
#endif

   if( expr->isintegral )
   {
      /* apply integrality to new bounds
       * it should be ok to use normal ceil() and floor(), but for safety, we use SCIPceil and SCIPfloor for now
       */
      if( newbounds.inf > -SCIP_INTERVAL_INFINITY )
         newbounds.inf = SCIPceil(scip, newbounds.inf);
      if( newbounds.sup <  SCIP_INTERVAL_INFINITY )
         newbounds.sup = SCIPfloor(scip, newbounds.sup);
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "applied integrality: [%.15g,%.15g]\n", newbounds.inf, newbounds.sup);
#endif
   }

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      SCIPdebugMsg(scip, "cut off due to new bounds being empty\n");

      SCIPintervalSetEmpty(&expr->activity);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* treat the new bounds as empty if either the lower/upper bound is above/below +/- SCIPinfinity() */
   if( SCIPisInfinity(scip, newbounds.inf) || SCIPisInfinity(scip, -newbounds.sup) )
   {
      SCIPdebugMsg(scip, "cut off due to new bounds being beyond infinity\n");

      SCIPintervalSetEmpty(&expr->activity);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if newbounds do not allow a sufficient tightening, then do not store new activity and do not add to queue for reverse propagation */
   if( !isIntervalBetter(scip, force, newbounds, expr->activity) )
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "new activity [%g,%g] for expr %p not sufficiently tighter than existing one -- rejecting update\n", newbounds.inf, newbounds.sup, (void*)expr);
#endif
      return SCIP_OKAY;
   }

   /* apply new bounds to expr activity, but ensure it's not-empty if very close disjoint intervals */
   SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, newbounds);

#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, "new activity [%g,%g] for expr %p\n", expr->activity.inf, expr->activity.sup, (void*)expr);
#endif

   /* check if the new bounds lead to an empty interval */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
   {
      SCIPdebugMsg(scip, "cut off due to empty intersection with previous activity\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if a reversepropagation queue is given, then add expression to that queue if it has at least one child and could have a reverseprop callback and we see a tightening */
   if( reversepropqueue != NULL && !expr->inqueue && (expr->nenfos > 0 || SCIPhasConsExprExprHdlrReverseProp(expr->exprhdlr)) )
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, " insert expr <%p> (%s) into reversepropqueue\n", (void*)expr, SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));
#endif
      SCIP_CALL( SCIPqueueInsert(reversepropqueue, expr) );
      expr->inqueue = TRUE;
   }

   /* update bounds on auxiliary variable */
   SCIP_CALL( tightenAuxVarBounds(scip, conshdlr, expr, force, cutoff, ntightenings) );

   return SCIP_OKAY;
}

/** mark constraints that include this expression to be propagated again
 *
 * This can be used by, e.g., nlhdlrs, to trigger a new propagation of constraints without
 * a change of variable bounds, e.g., because new information on the expression is available
 * that could potentially lead to tighter expression activity values.
 *
 * Note, that this call marks also constraints for propagation which only share some variable
 * with this expression.
 */
SCIP_RETCODE SCIPmarkConsExprExprPropagate(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to propagate again */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   for( ; !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) )  /*lint !e441*/
   {
      if( !SCIPisConsExprExprVar(expr) )
         continue;

      conss = SCIPgetConsExprExprVarConss(expr);
      nconss = SCIPgetConsExprExprVarNConss(expr);

      for( c = 0; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         consdata->ispropagated = FALSE;
      }
   }

   SCIPexpriteratorFree(&it);

   SCIPincrementConsExprCurBoundsTag(conshdlr, FALSE);

   return SCIP_OKAY;
}

/** increments the curboundstag and resets lastboundrelax in constraint handler data
 *
 * @note This method is not intended for normal use.
 *   These tags are maintained by the event handler for variable bound change events.
 *   This method is used by some unittests.
 */
void SCIPincrementConsExprCurBoundsTag(
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_Bool               boundrelax        /**< indicates whether a bound was relaxed, i.e., lastboundrelax should be set too */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);

   if( boundrelax )
      conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
}

/** adds violation-branching score to an expression
 *
 * Adds a score to the expression-specific violation-branching score, thereby marking it as branching candidate.
 * The expression must either be a variable expression or have an aux-variable.
 * In the latter case, branching on auxiliary variables must have been enabled.
 * In case of doubt, use SCIPaddConsExprExprsViolScore(). Roughly, the difference between these functions is that the current
 * function adds the violscore to the expression directly, while SCIPaddConsExprExprsViolScore() will split the
 * violation score among all the given expressions according to constraints/expr/branching/violsplit. See
 * SCIPaddConsExprExprsViolScore() for more details.
 */
void SCIPaddConsExprExprViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression where to add branching score */
   SCIP_Real               violscore         /**< violation score to add to expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(violscore >= 0.0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if not allowing to branch on auxvars, then expr must be a var-expr */
   assert(SCIPgetConsExprBranchAux(scip, conshdlr) || expr->exprhdlr == conshdlrdata->exprvarhdlr);
   /* if allowing to branch on auxvars, then expr must be a var-expr or have an auxvar */
   assert(!SCIPgetConsExprBranchAux(scip, conshdlr) || (expr->exprhdlr == conshdlrdata->exprvarhdlr || expr->auxvar != NULL));

   /* reset branching score if we are in a different enfo round */
   if( expr->violscoretag != conshdlrdata->enforound )
   {
      expr->violscoresum = violscore;
      expr->violscoremax = violscore;
      expr->nviolscores = 1;
      expr->violscoretag = conshdlrdata->enforound;
      return;
   }

   expr->violscoresum += violscore;
   if( violscore > expr->violscoremax )
      expr->violscoremax = violscore;
   ++(expr->nviolscores);
}

/** adds violation-branching score to a set of expressions, distributing the score among all the expressions.
 *
 * Each expression must either be a variable expression or have an aux-variable.
 * If branching on aux-variables is disabled, then the violation branching score will be distributed among all among the
 * variables present in exprs
 */
SCIP_RETCODE SCIPaddConsExprExprsViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation-branching score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one branchscore was added */
   )
{
   /* distribute violation as branching score to original variables in children of expr that are marked in branchcand */
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_CONSEXPR_EXPR* e;
   int nvars;
   int varssize;
   int i;

   assert(exprs != NULL || nexprs == 0);
   assert(success != NULL);

   if( nexprs == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* if allowing to branch on auxiliary variables, then call internal addConsExprExprsViolScore immediately */
   if( SCIPgetConsExprBranchAux(scip, conshdlr) )
   {
      addConsExprExprsViolScore(scip, conshdlr, exprs, nexprs, violscore, sol, success);
      return SCIP_OKAY;
   }

   /* if not allowing to branch on aux vars, then create new array containing var expressions that exprs depend on */
   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, varssize) );

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, NULL, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   for( i = 0; i < nexprs; ++i )
   {
      for( e = SCIPexpriteratorRestartDFS(it, exprs[i]); !SCIPexpriteratorIsEnd(it); e = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         assert(e != NULL);

         if( SCIPisConsExprExprVar(e) )
         {
            /* add variable expression to vars array */
            if( varssize == nvars )
            {
               varssize = SCIPcalcMemGrowSize(scip, nvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &varexprs, varssize) );
            }
            assert(varssize > nvars);

            varexprs[nvars++] = e;
         }
      }
   }

   SCIPexpriteratorFree(&it);

   addConsExprExprsViolScore(scip, conshdlr, varexprs, nvars, violscore, sol, success);

   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_Real SCIPgetConsExprExprViolScore(
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(expr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->enforound != expr->violscoretag )
      return 0.0;

   if( expr->nviolscores == 0 )
      return 0.0;

   switch( conshdlrdata->branchscoreagg )
   {
      case 'a' :
         /* average */
         return expr->violscoresum / expr->nviolscores;

      case 'm' :
         /* maximum */
         return expr->violscoremax;

      case 's' :
         /* sum */
         return expr->violscoresum;

      default:
         SCIPerrorMessage("Invalid value %c for branchscoreagg parameter\n", conshdlrdata->branchscoreagg);
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** returns the hash value of an expression */
SCIP_RETCODE SCIPgetConsExprExprHash(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   unsigned int*           hashval           /**< pointer to store the hash value */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashval != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, SCIPfindConshdlr(scip, CONSHDLR_NAME), SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   SCIP_CALL( hashExpr(scip, expr, it, NULL) );

   *hashval = SCIPexpriteratorGetExprUserData(it, expr).uintval;

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}


/** creates and gives the auxiliary variable for a given expression
 *
 * @note if auxiliary variable already present for that expression, then only returns this variable
 * @note for a variable expression it returns the corresponding variable
 */
SCIP_RETCODE SCIPcreateConsExprExprAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR**            auxvar              /**< buffer to store pointer to auxiliary variable, or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VARTYPE vartype;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(expr != NULL);

   /* if we already have auxvar, then just return it */
   if( expr->auxvar != NULL )
   {
      if( auxvar != NULL )
         *auxvar = expr->auxvar;
      return SCIP_OKAY;
   }

   /* if expression is a variable-expression, then return that variable */
   if( expr->exprhdlr == SCIPgetConsExprExprHdlrVar(conshdlr) )
   {
      if( auxvar != NULL )
         *auxvar = SCIPgetConsExprExprVarVar(expr);
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);

   /* it doesn't harm much to have an auxvar for a constant, as this can be handled well by the default hdlr,
    * but it usually indicates a missing simplify
    * if we find situations where we need to have an auxvar for a constant, then remove this assert
    */
   assert(expr->exprhdlr != SCIPgetConsExprExprHdlrValue(conshdlr));

   /* create and capture auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxvar_%s_%d", expr->exprhdlr->name, conshdlrdata->auxvarid);
   ++conshdlrdata->auxvarid;

   /* type of auxiliary variable depends on integrality information of the expression */
   vartype = SCIPisConsExprExprIntegral(expr) ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS;

   SCIP_CALL( SCIPcreateVarBasic(scip, &expr->auxvar, name, MAX( -SCIPinfinity(scip), expr->activity.inf ),
      MIN( SCIPinfinity(scip), expr->activity.sup ), 0.0, vartype) ); /*lint !e666*/
   SCIP_CALL( SCIPaddVar(scip, expr->auxvar) );

   /* mark the auxiliary variable to be added for the relaxation only
    * this prevents SCIP to create linear constraints from cuts or conflicts that contain auxiliary variables,
    * or to copy the variable to a subscip
    */
   SCIPvarMarkRelaxationOnly(expr->auxvar);

   SCIPdebugMsg(scip, "added auxiliary variable %s [%g,%g] for expression %p\n", SCIPvarGetName(expr->auxvar), SCIPvarGetLbGlobal(expr->auxvar), SCIPvarGetUbGlobal(expr->auxvar), (void*)expr);

   /* add variable locks in both directions */
   SCIP_CALL( SCIPaddVarLocks(scip, expr->auxvar, 1, 1) );

#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      /* store debug solution value of auxiliary variable
       * assumes that expression has been evaluated in debug solution before
       */
      SCIP_CALL( SCIPdebugAddSolVal(scip, expr->auxvar, SCIPgetConsExprExprValue(expr)) );
   }
#endif

   /* catch bound change events on this variable, since bounds on this variable take part of activity computation */
   assert(expr->auxfilterpos == -1);
   SCIP_CALL( SCIPcatchVarEvent(scip, expr->auxvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)expr, &expr->auxfilterpos) );

   if( auxvar != NULL )
      *auxvar = expr->auxvar;

   return SCIP_OKAY;
}

/** @name Simplifying methods
 *
 * This is largely inspired in Joel Cohen's
 * Computer algebra and symbolic computation: Mathematical methods
 * In particular Chapter 3
 * The other fountain of inspiration is the current simplifying methods in expr.c.
 *
 * Note: The things to keep in mind when adding simplification rules are the following.
 * I will be using the product expressions as an example.
 * There are mainly 3 parts of the simplification process. You need to decide
 * at which stage the simplification rule makes sense.
 * 1. Simplify each factor (simplifyFactor): At this stage we got the children of the product expression.
 * At this point, each child is simplified when viewed as a stand-alone
 * expression, but not necessarily when viewed as child of a product
 * expression. Rules like SP2, SP7, etc are enforced at this point.
 * 2. Multiply the factors (mergeProductExprlist): At this point rules like SP4, SP5 and SP14 are enforced.
 * 3. Build the actual simplified product expression (buildSimplifiedProduct):
 * At this point rules like SP10, SP11, etc are enforced.
 *
 * **During step 1. and 2. do not forget to set the flag changed to TRUE when something actually changes**
 *
 * Definition of simplified expressions
 * ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * An expression is simplified if it
 * - is a value expression
 * - is a var expression
 * - is a product expression such that
 *    SP1:  every child is simplified
 *    SP2:  no child is a product
 *    SP4:  no two children are the same expression (those should be multiplied)
 *    SP5:  the children are sorted [commutative rule]
 *    SP7:  no child is a value
 *    SP8:  its coefficient is 1.0 (otherwise should be written as sum)
 *    SP10: it has at least two children
 *    ? at most one child is an abs
 *    SP11: no two children are expr*log(expr)
 *    (TODO: we could handle more complicated stuff like x*y*log(x) -> - y * entropy(x), but I am not sure this should
 *    happen at the simplifcation level, or (x*y) * log(x*y), which currently simplifies to x * y * log(x*y))
 *    SP12: if it has two children, then neither of them is a sum (expand sums)
 *    SP13: no child is a sum with a single term
 *    SP14: at most one child is an exp
 * - is a (signed)power expression such that
 *   TODO: Some of these criteria are too restrictive for signed powers; for example, the exponent does not need to be
 *   an integer for signedpower to distribute over a product (POW5, POW6, POW8). Others can also be improved
 *    POW1: exponent is not 0
 *    POW2: exponent is not 1
 *    POW3: its child is not a value
 *    POW4: its child is simplified
 *    POW5: if exponent is integer, its child is not a product
 *    POW6: if exponent is integer, its child is not a sum with a single term ((2*x)^2 -> 4*x^2)
 *    POW7: if exponent is 2, its child is not a sum (expand sums)
 *    POW8: if exponent is integer, its child is not a power
 *    POW9: its child is not a sum with a single term with a positive coefficient: (25*x)^0.5 -> 5 x^0.5
 *    POW10: its child is not a binary variable: b^e and e > 0 --> b, b^e and e < 0 --> fix b to 1
 *    POW11: its child is not an exponential: exp(expr)^e --> exp(e * expr)
 * - is a signedpower expression such that
 *   TODO: Some of these criteria are too restrictive for signed powers; for example, the exponent does not need to be
 *   an integer for signedpower to distribute over a product (SPOW5, SPOW6, SPOW8). Others can also be improved
 *    SPOW1: exponent is not 0
 *    SPOW2: exponent is not 1
 *    SPOW3: its child is not a value
 *    SPOW4: its child is simplified
 *    SPOW5: (TODO) do we want to distribute signpowers over products like we do powers?
 *    SPOW6: exponent is not an odd integer: (signpow odd expr) -> (pow odd expr)
 *    SPOW8: if exponent is integer, its child is not a power
 *    SPOW9: its child is not a sum with a single term: (25*x)^0.5 -> 5 x^0.5
 *    SPOW10: its child is not a binary variable: b^e and e > 0 --> b, b^e and e < 0 --> fix b to 1
 *    SPOW11: its child is not an exponential: exp(expr)^e --> exp(e * expr)
 *    SPOW?: TODO: what happens when child is another signed power?
 *    SPOW?: if child >= 0 -> transform to normal power; if child < 0 -> transform to - normal power
 * - is a sum expression such that
 *    SS1: every child is simplified
 *    SS2: no child is a sum
 *    SS3: no child is a value (values should go in the constant of the sum)
 *    SS4: no two children are the same expression (those should be summed up)
 *    SS5: the children are sorted [commutative rule]
 *    SS6: it has at least one child
 *    SS7: if it consists of a single child, then either constant is != 0.0 or coef != 1
 *    SS8: no child has coefficient 0
 *    SS9: if a child c is a product that has an exponential expression as one of its factors, then the coefficient of c is +/-1.0
 *    SS10: if a child c is an exponential, then the coefficient of c is +/-1.0 (TODO)
 *    x if it consists of a single child, then its constant != 0.0 (otherwise, should be written as a product)
 * - it is a function with simplified arguments, but not all of them can be values
 * ? a logarithm doesn't have a product as a child
 * ? the exponent of an exponential is always 1
 *
 * ORDERING RULES
 * ^^^^^^^^^^^^^^
 * These rules define a total order on *simplified* expressions.
 * There are two groups of rules, when comparing equal type expressions and different type expressions
 * Equal type expressions:
 * OR1: u,v value expressions: u < v <=> val(u) < val(v)
 * OR2: u,v var expressions: u < v <=> SCIPvarGetIndex(var(u)) < SCIPvarGetIndex(var(v))
 * OR3: u,v are both sum or product expression: < is a lexicographical order on the terms
 * OR4: u,v are both pow: u < v <=> base(u) < base(v) or, base(u) == base(v) and expo(u) < expo(v)
 * OR5: u,v are u = FUN(u_1, ..., u_n), v = FUN(v_1, ..., v_m): u < v <=> For the first k such that u_k != v_k, u_k < v_k,
 *      or if such a k doesn't exist, then n < m.
 *
 * Different type expressions:
 * OR6: u value, v other: u < v always
 * OR7: u sum, v var or func: u < v <=> u < 0+v
 *      In other words, u = sum_{i = 1}^n alpha_i u_i, then u < v <=> u_n < v or if u_n = v and alpha_n < 1
 * OR8: u product, v pow, sum, var or func: u < v <=> u < 1*v
 *      In other words, u = Pi_{i = 1}^n u_i,  then u < v <=> u_n < v
 *      @note: since this applies only to simplified expressions, the form of the product is correct. Simplified products
 *             do *not* have constant coefficients
 * OR9: u pow, v sum, var or func: u < v <=> u < v^1
 * OR10: u var, v func: u < v always
 * OR11: u func, v other type of func: u < v <=> name(type(u)) < name(type(v))
 * OR12: none of the rules apply: u < v <=> ! v < u
 * Examples:
 * OR12: x < x^2 ?:  x is var and x^2 product, so none applies.
 *       Hence, we try to answer x^2 < x ?: x^2 < x <=> x < x or if x = x and 2 < 1 <=> 2 < 1 <=> False, so x < x^2 is True
 *       x < x^-1 --OR12--> ~(x^-1 < x) --OR9--> ~(x^-1 < x^1) --OR4--> ~(x < x or -1 < 1) --> ~True --> False
 *       x*y < x --OR8--> x*y < 1*x --OR3--> y < x --OR2--> False
 *       x*y < y --OR8--> x*y < 1*y --OR3--> y < x --OR2--> False
 *
 * Algorithm
 * ^^^^^^^^^
 * The recursive version of the algorithm is
 *
 * EXPR simplify(expr)
 *    for c in 1..expr->nchildren
 *       expr->children[c] = simplify(expr->children[c])
 *    end
 *    return expr->exprhdlr->simplify(expr)
 * end
 *
 * Important: Whatever is returned by a simplify callback **has** to be simplified.
 * Also, all children of the given expression **are** already simplified
 *
 * @{
 */

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note: The given expressions are assumed to be simplified.
 */
int SCIPcompareConsExprExprs(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr1;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr2;
   int retval;

   exprhdlr1 = SCIPgetConsExprExprHdlr(expr1);
   exprhdlr2 = SCIPgetConsExprExprHdlr(expr2);

   /* expressions are of the same kind/type; use compare callback or default method */
   if( exprhdlr1 == exprhdlr2 )
   {
      return SCIPcompareConsExprExprHdlr(expr1, expr2);
   }

   /* expressions are of different kind/type */
   /* enforces OR6 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "val") == 0 )
   {
      return -1;
   }
   /* enforces OR12 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "val") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR7 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "sum") == 0 )
   {
      int compareresult;
      int nchildren;

      nchildren = SCIPgetConsExprExprNChildren(expr1);
      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* "base" of the largest expression of the sum is equal to expr2, coefficient might tell us that expr2 is larger */
      if( SCIPgetConsExprExprSumCoefs(expr1)[nchildren-1] < 1.0 )
         return -1;

      /* largest expression of sum is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "sum") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR8 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "prod") == 0 )
   {
      int compareresult;
      int nchildren;

      nchildren = SCIPgetConsExprExprNChildren(expr1);
      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* largest expression of product is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "prod") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR9 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "pow") == 0 )
   {
      int compareresult;

      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[0], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* base equal to expr2, exponent might tell us that expr2 is larger */
      if( SCIPgetConsExprExprPowExponent(expr1) < 1.0 )
         return -1;

      /* power expression is larger => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "pow") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR10 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "var") == 0 )
      return -1;
   /* enforces OR12 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "var") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR11 */
   retval = strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), SCIPgetConsExprExprHdlrName(exprhdlr2));
   return retval == 0 ? 0 : retval < 0 ? -1 : 1;
}

/** simplifies an expression
 *
 * The given expression will be released and overwritten with the simplified expression.
 * To keep the expression, duplicate it via SCIPduplicateConsExprExpr before calling this method.
 */
SCIP_RETCODE SCIPsimplifyConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**    simplified,       /**< buffer to store simplified expression */
   SCIP_Bool*              changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*              infeasible        /**< buffer to store whether infeasibility has been detected */
)
{
   assert(rootexpr != NULL);
   assert(simplified != NULL);

   SCIP_CALL( reformulateConsExprExpr(scip, conshdlr, rootexpr, TRUE, simplified, changed, infeasible) );

   return SCIP_OKAY;
}

/**@} */  /* end of simplifying methods */

/** reformulate an expression; this functions works similar as SCIPsimplifyConsExprExpr() but instead of calling the
 *  simplify callback of an expression handler it iterates through all nonlinear handlers and uses the reformulation
 *  callback
 */
SCIP_RETCODE SCIPreformulateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**    refrootexpr,      /**< buffer to store reformulated expression */
   SCIP_Bool*              changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*              infeasible        /**< buffer to store whether infeasibility has been detected */
   )
{
   assert(rootexpr != NULL);
   assert(refrootexpr != NULL);

   SCIP_CALL( reformulateConsExprExpr(scip, conshdlr, rootexpr, FALSE, refrootexpr, changed, infeasible) );

   return SCIP_OKAY;
}

/** sets the curvature of an expression */
void SCIPsetConsExprExprCurvature(
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   )
{
   assert(expr != NULL);
   expr->curvature = curvature;
}

/** returns the curvature of an expression */
SCIP_EXPRCURV SCIPgetConsExprExprCurvature(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->curvature;
}

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 */
SCIP_RETCODE SCIPcomputeConsExprExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPRCURV curv;
   SCIP_INTERVAL activity;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   SCIP_Bool success;
   SCIP_EXPRCURV trialcurv[3] = { SCIP_EXPRCURV_LINEAR, SCIP_EXPRCURV_CONVEX, SCIP_EXPRCURV_CONCAVE };
   int i, c;

   assert(scip != NULL);
   assert(expr != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* ensure activities are uptodate */
   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &activity, TRUE) );

   childcurvsize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) )  /*lint !e441*/
   {
      curv = SCIP_EXPRCURV_UNKNOWN;

      if( expr->exprhdlr->curvature == NULL )
      {
         /* set curvature in expression */
         SCIPsetConsExprExprCurvature(expr, curv);
         continue;
      }

      if( SCIPgetConsExprExprNChildren(expr) > childcurvsize )
      {
         childcurvsize = SCIPcalcMemGrowSize(scip, SCIPgetConsExprExprNChildren(expr));
         SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
      }

      /* SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
      SCIPinfoMessage(scip, NULL, " (%p)", expr); */
      for( i = 0; i < 3; ++i )
      {
         /* check if expression can have a curvature trialcurv[i] */
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, expr, trialcurv[i], &success, childcurv) );
         /* SCIPinfoMessage(scip, NULL, " %s? %d", SCIPexprcurvGetName(trialcurv[i]), success); */
         if( !success )
            continue;

         /* check if conditions on children are satisfied */
         for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
         {
            if( (childcurv[c] & SCIPgetConsExprExprCurvature(SCIPgetConsExprExprChildren(expr)[c])) != childcurv[c] )
            {
               success = FALSE;
               break;
            }
         }

         if( success )
         {
            curv = trialcurv[i];
            break;
         }
      }

      /* set curvature in expression */
      SCIPsetConsExprExprCurvature(expr, curv);
      /* SCIPinfoMessage(scip, NULL, " -> curv = %s\n", SCIPexprcurvGetName(curv)); */
   }

   SCIPexpriteratorFree(&it);

   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

/** returns the monotonicity of an expression w.r.t. to a given child */
SCIP_MONOTONE SCIPgetConsExprExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   childidx            /**< index of child */
   )
{
   SCIP_MONOTONE monotonicity = SCIP_MONOTONE_UNKNOWN;

   assert(expr != NULL);
   assert(childidx >= 0 || expr->nchildren == 0);
   assert(childidx < expr->nchildren);

   /* check whether the expression handler implements the monotonicity callback */
   if( expr->exprhdlr->monotonicity != NULL )
   {
      SCIP_CALL_ABORT( (*expr->exprhdlr->monotonicity)(scip, expr, childidx, &monotonicity) );
   }

   return monotonicity;
}

/** returns the number of positive rounding locks of an expression */
int SCIPgetConsExprExprNLocksPos(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nlockspos;
}

/** returns the number of negative rounding locks of an expression */
int SCIPgetConsExprExprNLocksNeg(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nlocksneg;
}

/** computes integrality information of a given expression and all its subexpressions; the integrality information can
 * be accessed via SCIPisConsExprExprIntegral()
 */
SCIP_RETCODE SCIPcomputeConsExprExprIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_LEAVEEXPR);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      /* compute integrality information */
      expr->isintegral = FALSE;

      if( expr->exprhdlr->integrality != NULL )
      {
         /* get curvature from expression handler */
         SCIP_CALL( (*expr->exprhdlr->integrality)(scip, expr, &expr->isintegral) );
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** returns whether an expression is integral */
SCIP_Bool SCIPisConsExprExprIntegral(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->isintegral;
}

/** returns the total number of variables in an expression
 *
 * The function counts variables in common sub-expressions only once.
 */
SCIP_RETCODE SCIPgetConsExprExprNVars(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   int*                    nvars             /**< buffer to store the total number of variables */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nvars != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   *nvars = 0;
   for( ; !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      if( SCIPisConsExprExprVar(expr) )
         ++(*nvars);

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** returns all variable expressions contained in a given expression; the array to store all variable expressions needs
 * to be at least of size the number of unique variables in the expression which is given by SCIPgetConsExprExprNVars()
 * and can be bounded by SCIPgetNVars().
 *
 * @note function captures variable expressions
 */
SCIP_RETCODE SCIPgetConsExprExprVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPR**    varexprs,         /**< array to store all variable expressions */
   int*                    nvarexprs         /**< buffer to store the total number of variable expressions */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(expr != NULL);
   assert(varexprs != NULL);
   assert(nvarexprs != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );

   *nvarexprs = 0;
   for( ; !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      assert(expr != NULL);

      if( SCIPisConsExprExprVar(expr) )
      {
         /* add variable expression to array and capture expr */
         assert(SCIPgetNTotalVars(scip) >= *nvarexprs + 1);

         varexprs[(*nvarexprs)++] = expr;

         /* capture expression */
         SCIPcaptureConsExprExpr(expr);
      }
   }

   /* @todo sort variable expressions here? */

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(x) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(x) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(x).
 *
 * If necessary, f is evaluated in the given solution. If that fails (domain error),
 * then viol is set to SCIPinfinity and both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetConsExprExprAbsOrigViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* make sure expression has been evaluated */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol, soltag) );

   /* get violation from internal method */
   *viol = getExprAbsOrigViolation(scip, expr, sol, violunder, violover);

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(w) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(w) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(w).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetConsExprExprAbsAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   return SCIP_OKAY;
}

/** computes relative violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * Taking the absolute violation from SCIPgetConsExprExprAbsAuxViolation, this function returns
 * the absolute violation divided by max(1,|f(w)|).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetConsExprExprRelAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   if( !SCIPisInfinity(scip, *viol) )
   {
      assert(auxvalue != SCIP_INVALID);  /*lint !e777*/
      *viol /= MAX(1.0, REALABS(auxvalue));  /*lint !e666*/
   }

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** gets the index an expression iterator can use to store iterator specific data in an expression */
SCIP_RETCODE SCIPactivateConsExprExprHdlrIterator(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   int*                       iterindex       /**< buffer to store iteration index */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   assert(iterindex != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->nactiveiter + 1 >= SCIP_CONSEXPRITERATOR_MAXNACTIVE )
   {
      SCIPerrorMessage("Maximal number of active expression iterators reached.\n");
      return SCIP_MAXDEPTHLEVEL;
   }

   *iterindex = conshdlrdata->nactiveiter++;

   return SCIP_OKAY;
}

/** returns the index that an expression iterator used to store iterator specific data in an expression */
void SCIPdeactivateConsExprExprHdlrIterator(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   int                        iterindex       /**< iteration index that is not used anymore */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   assert(iterindex >= 0);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* the iterindex must be the one of the last initialized iterator */
   assert(iterindex == conshdlrdata->nactiveiter-1);

   --conshdlrdata->nactiveiter;
}

/** get a new tag that can be used to mark an expression as visited */
unsigned int SCIPgetConsExprExprHdlrNewVisitedTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   return ++conshdlrdata->lastvisitedtag;
}

/** gets tag indicating current local variable bounds */
unsigned int SCIPgetConsExprCurBoundsTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);

   return conshdlrdata->curboundstag;
}

/** gets the curboundstag at the last time where variable bounds were relaxed */
unsigned int SCIPgetConsExprLastBoundRelaxTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);

   return conshdlrdata->lastboundrelax;
}

/** collects all bilinear terms for a given set of constraints
 *
 * @note This method should only be used for unit tests that depend on SCIPgetConsExprBilinTerms()
 *       or SCIPgetConsExprBilinTerm().
 */
SCIP_RETCODE SCIPcollectConsExprBilinTerms(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_CONS**                conss,          /**< expression constraints */
   int                        nconss          /**< total number of expression constraints */
   )
{
   assert(consexprhdlr != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( bilinearTermsInsertAll(scip, consexprhdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** returns the total number of bilinear terms that are contained in all expression constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
int SCIPgetConsExprNBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nbilinterms;
}

/** returns all bilinear terms that are contained in all expression constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->bilinterms;
}

/** returns the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerm(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_BILINTERM entry;
   int idx;

   assert(consexprhdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   /* use a new entry to find the image in the bilinear hash table */
   entry.x = x;
   entry.y = y;
   idx = (int)(size_t)SCIPhashtableRetrieve(conshdlrdata->bilinhashtable, (void*)&entry) - 1;
   assert(idx >= -1 && idx < conshdlrdata->nbilinterms);

   if( idx >= 0 )
   {
      assert(conshdlrdata->bilinterms[idx].x == x);
      assert(conshdlrdata->bilinterms[idx].y == y);
      return &conshdlrdata->bilinterms[idx];
   }

   return NULL;
}

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create expr constraint handler data */
   SCIP_CALL( SCIPallocClearMemory(scip, &conshdlrdata) );
   conshdlrdata->lastsoltag = 1;
   conshdlrdata->curboundstag = 1;
   conshdlrdata->lastboundrelax = 1;
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->canonicalizetime) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyExpr,
         consFreeExpr, consInitExpr, consExitExpr,
         consInitpreExpr, consExitpreExpr, consInitsolExpr, consExitsolExpr,
         consDeleteExpr, consTransExpr, consInitlpExpr,
         consSepalpExpr, consSepasolExpr, consEnfolpExpr, consEnforelaxExpr, consEnfopsExpr, consCheckExpr,
         consPropExpr, consPresolExpr, consRespropExpr, consLockExpr,
         consActiveExpr, consDeactiveExpr,
         consEnableExpr, consDisableExpr, consDelvarsExpr,
         consPrintExpr, consCopyExpr, consParseExpr,
         consGetVarsExpr, consGetNVarsExpr, consGetDiveBdChgsExpr, conshdlrdata) );

   if( SCIPfindConshdlr(scip, "quadratic") != NULL )
   {
      /* include function that upgrades quadratic constraint to expr constraints */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, quadconsUpgdExpr, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdExpr, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add expr constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelax",
         "strategy on how to relax variable bounds during bound tightening: relax (n)ot, relax by (a)bsolute value, relax always by a(b)solute value, relax by (r)relative value",
         &conshdlrdata->varboundrelax, TRUE, 'r', "nabr", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelaxamount",
         "by how much to relax variable bounds during bound tightening if strategy 'a', 'b', or 'r'",
         &conshdlrdata->varboundrelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/conssiderelaxamount",
         "by how much to relax constraint sides during bound tightening",
         &conshdlrdata->conssiderelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpmaxperturb",
         "maximal relative perturbation of reference point when computing facet of envelope of vertex-polyhedral function (dim>2)",
         &conshdlrdata->vp_maxperturb, TRUE, VERTEXPOLY_MAXPERTURBATION, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpadjfacetthresh",
         "adjust computed facet of envelope of vertex-polyhedral function up to a violation of this value times LP feasibility tolerance",
         &conshdlrdata->vp_adjfacetthreshold, TRUE, VERTEXPOLY_ADJUSTFACETFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/vpdualsimplex",
         "whether to use dual simplex instead of primal simplex for LP that computes facet of vertex-polyhedral function",
         &conshdlrdata->vp_dualsimplex, TRUE, VERTEXPOLY_USEDUALSIMPLEX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprods",
         "whether to reformulate products of binary variables during presolving",
         &conshdlrdata->reformbinprods, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsand",
         "whether to use the AND constraint handler for reformulating binary products",
         &conshdlrdata->reformbinprodsand, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsfac",
         "minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled)",
         &conshdlrdata->reformbinprodsfac, FALSE, 50, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/tightenlpfeastol",
         "whether to tighten LP feasibility tolerance during enforcement, if it seems useful",
         &conshdlrdata->tightenlpfeastol, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propinenforce",
         "whether to (re)run propagation in enforcement",
         &conshdlrdata->propinenforce, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutthreshold",
         "threshold for when to regard a cut from an estimator as weak (lower values allow more weak cuts)",
         &conshdlrdata->weakcutthreshold, TRUE, 0.2, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/strongcutmaxcoef",
         "\"strong\" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef]",
         &conshdlrdata->strongcutmaxcoef, TRUE, 1000.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/strongcutefficacy",
         "consider efficacy requirement when deciding whether a cut is \"strong\"",
         &conshdlrdata->strongcutefficacy, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forcestrongcut",
         "whether to force \"strong\" cuts in enforcement",
         &conshdlrdata->forcestrongcut, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/enfoauxviolfactor",
         "an expression will be enforced if the \"auxiliary\" violation is at least this factor times the \"original\" violation",
         &conshdlrdata->enfoauxviolfactor, TRUE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutminviolfactor",
         "retry enfo of constraint with weak cuts if violation is least this factor of maximal violated constraints",
         &conshdlrdata->weakcutminviolfactor, TRUE, 0.5, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/violscale",
         "method how to scale violations to make them comparable (not used for feasibility check): (n)one, (a)ctivity and side, norm of (g)radient",
         &conshdlrdata->violscale, TRUE, 'n', "nag", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/checkvarlocks",
         "whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction)",
         &conshdlrdata->checkvarlocks, TRUE, 't', "bdt", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/branching/aux",
         "from which depth on in the tree to allow branching on auxiliary variables (variables added for extended formulation)",
         &conshdlrdata->branchauxmindepth, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branching/external",
         "whether to use external branching candidates and branching rules for branching",
         &conshdlrdata->branchexternal, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highviolfactor",
         "consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints",
         &conshdlrdata->branchhighviolfactor, FALSE, 0.0, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highscorefactor",
         "consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables",
         &conshdlrdata->branchhighscorefactor, FALSE, 0.9, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/violweight",
         "weight by how much to consider the violation assigned to a variable for its branching score",
         &conshdlrdata->branchviolweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/dualweight",
         "weight by how much to consider the dual values of rows that contain a variable for its branching score",
         &conshdlrdata->branchdualweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostweight",
         "weight by how much to consider the pseudo cost of a variable for its branching score",
         &conshdlrdata->branchpscostweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/domainweight",
         "weight by how much to consider the domain width in branching score",
         &conshdlrdata->branchdomainweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/vartypeweight",
         "weight by how much to consider variable type (continuous: 0, binary: 1, integer: 0.1, impl-integer: 0.01) in branching score",
         &conshdlrdata->branchvartypeweight, FALSE, 0.5, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/scoreagg",
         "how to aggregate several branching scores given for the same expression: 'a'verage, 'm'aximum, 's'um",
         &conshdlrdata->branchscoreagg, TRUE, 's', "ams", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/violsplit",
         "method used to split violation in expression onto variables: 'e'venly, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width",
         &conshdlrdata->branchviolsplit, TRUE, 'm', "emdl", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostreliable",
         "minimum pseudo-cost update count required to consider pseudo-costs reliable",
         &conshdlrdata->branchpscostreliable, FALSE, 2.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   /* include handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, CONSHDLR_NAME "_boundchange",
         "signals a bound change to an expression constraint", processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   /* include table for statistics */
   assert(SCIPfindTable(scip, TABLE_NAME_EXPR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_EXPR, TABLE_DESC_EXPR, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputExpr,
         NULL, TABLE_POSITION_EXPR, TABLE_EARLIEST_STAGE_EXPR) );

   return SCIP_OKAY;
}

/** creates the handler for expr constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( includeConshdlrExprBasic(scip) );

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* include and remember handler for variable expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "var") == 0);
   conshdlrdata->exprvarhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for constant value expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrValue(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "val") == 0);
   conshdlrdata->exprvalhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for sum expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sum") == 0);
   conshdlrdata->exprsumhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for product expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "prod") == 0);
   conshdlrdata->exprprodhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for exponential expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrExp(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "exp") == 0);
   conshdlrdata->exprexphdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for logarithmic expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrLog(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "log") == 0);

   /* include handler for absolute expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrAbs(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "abs") == 0);

   /* include handler for power expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrPow(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "pow") == 0);
   conshdlrdata->exprpowhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for signed power expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSignpower(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "signpower") == 0);
   conshdlrdata->exprsignpowhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for entropy expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrEntropy(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "entropy") == 0);

   /* include handler for sine expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSin(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sin") == 0);

   /* include handler for cosine expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrCos(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "cos") == 0);

   /* include default nonlinear handler */
   SCIP_CALL( SCIPincludeConsExprNlhdlrDefault(scip, conshdlr) );

   /* include nonlinear handler for quadratics */
   SCIP_CALL( SCIPincludeConsExprNlhdlrQuadratic(scip, conshdlr) );

   /* include nonlinear handler for convex expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrConvex(scip, conshdlr) );

   /* include nonlinear handler for concave expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrConcave(scip, conshdlr) );

   /* include nonlinear handler for bilinear expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrBilinear(scip, conshdlr) );

   /* include nonlinear handler for perspective reformulations */
   SCIP_CALL( SCIPincludeConsExprNlhdlrPerspective(scip, conshdlr) );

   return SCIP_OKAY;
}

/** includes an expression constraint upgrade method into the expression constraint handler */
SCIP_RETCODE SCIPincludeExprconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_EXPRCONSUPGD((*exprconsupgd)),  /**< method to call for upgrading expression constraint, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_EXPRCONSUPGRADE* exprconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(conshdlrname != NULL );

   /* ignore empty upgrade functions */
   if( exprconsupgd == NULL )
      return SCIP_OKAY;

   /* find the expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether upgrade method exists already */
   for( i = conshdlrdata->nexprconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->exprconsupgrades[i]->exprconsupgd == exprconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade method %p for constraint handler <%s>.\n", exprconsupgd, conshdlrname); /*lint !e611*/
#endif
         return SCIP_OKAY;
      }
   }

   /* create an expression constraint upgrade data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &exprconsupgrade) );
   exprconsupgrade->exprconsupgd = exprconsupgd;
   exprconsupgrade->priority   = priority;
   exprconsupgrade->active     = active;

   /* insert expression constraint upgrade method into constraint handler data */
   assert(conshdlrdata->nexprconsupgrades <= conshdlrdata->exprconsupgradessize);
   if( conshdlrdata->nexprconsupgrades+1 > conshdlrdata->exprconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nexprconsupgrades+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->exprconsupgrades, conshdlrdata->nexprconsupgrades, newsize) );
      conshdlrdata->exprconsupgradessize = newsize;
   }
   assert(conshdlrdata->nexprconsupgrades+1 <= conshdlrdata->exprconsupgradessize);

   for( i = conshdlrdata->nexprconsupgrades; i > 0 && conshdlrdata->exprconsupgrades[i-1]->priority < exprconsupgrade->priority; --i )
      conshdlrdata->exprconsupgrades[i] = conshdlrdata->exprconsupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nexprconsupgrades);
   conshdlrdata->exprconsupgrades[i] = exprconsupgrade;
   conshdlrdata->nexprconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable expression upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &exprconsupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsExpr() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(expr != NULL);

   /* find the expr constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("expr constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* TODO remove this once we allow for local expression constraints */
   if( local && SCIPgetDepth(scip) != 0 )
   {
      SCIPerrorMessage("Locally valid expression constraints are not supported, yet.\n");
      return SCIP_INVALIDCALL;
   }

   /* TODO remove this once we allow for non-initial expression constraints */
   if( !initial )
   {
      SCIPerrorMessage("Non-initial expression constraints are not supported, yet.\n");
      return SCIP_INVALIDCALL;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &consdata) );
   consdata->expr = expr;
   consdata->lhs = lhs;
   consdata->rhs = rhs;
   consdata->consindex = conshdlrdata->lastconsindex++;

   /* capture expression */
   SCIPcaptureConsExprExpr(consdata->expr);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic expression constraint */
SCIP_RETCODE SCIPcreateConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
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
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_Real* coefs;
   int i;

   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadterms == 0 || (quadvars1 != NULL && quadvars2 != NULL && quadcoefs != NULL));

   /* get expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadterms + nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nquadterms + nlinvars) );

   /* create children for quadratic terms */
   for( i = 0; i < nquadterms; ++i )
   {
      assert(quadvars1 != NULL && quadvars1[i] != NULL);
      assert(quadvars2 != NULL && quadvars2[i] != NULL);

      /* quadratic term */
      if( quadvars1[i] == quadvars2[i] )
      {
         SCIP_CONSEXPR_EXPR* xexpr;

         /* create variable expression */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, quadvars1[i]) );

         /* create pow expression */
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &children[i], xexpr, 2.0) );

         /* release variable expression; note that the variable expression is still captured by children[i] */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
      }
      else /* bilinear term */
      {
         SCIP_CONSEXPR_EXPR* exprs[2];

         /* create variable expressions */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[0], quadvars1[i]) );
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[1], quadvars2[i]) );

         /* create product expression */
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &children[i], 2, exprs, 1.0) );

         /* release variable expressions; note that the variable expressions are still captured by children[i] */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[1]) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[0]) );
      }

      /* store coefficient */
      coefs[i] = quadcoefs[i];
   }

   /* create children for linear terms */
   for( i = 0; i < nlinvars; ++i )
   {
      assert(linvars != NULL && linvars[i] != NULL);

      /* create variable expression; release variable expression after the constraint has been created */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &children[nquadterms + i], linvars[i]) );

      /* store coefficient */
      coefs[nquadterms + i] = lincoefs[i];
   }

   /* create sum expression */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, nquadterms + nlinvars, children, coefs, 0.0) );

   /* create expression constraint */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name, sumexpr, lhs, rhs, initial, separate, enforce, check, propagate,
      local, modifiable, dynamic, removable) );

   /* release sum expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );

   /* release children */
   for( i = 0; i < nquadterms + nlinvars; ++i )
   {
      assert(children[i] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[i]) );
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &children);

   return SCIP_OKAY;
}

/** returns the expression of the given expression constraint */
SCIP_CONSEXPR_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not expression\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** gets the left hand side of an expression constraint */
SCIP_Real SCIPgetLhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets the right hand side of an expression constraint */
SCIP_Real SCIPgetRhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets absolute violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(scip), the constraint would be considered feasible.
 */
SCIP_RETCODE SCIPgetAbsViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0) );
   *viol = getConsAbsViolation(cons);

   return SCIP_OKAY;
}

/** gets scaled violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/expr/violscale.
 */
SCIP_RETCODE SCIPgetRelViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0) );
   SCIP_CALL( getConsRelViolation(scip, cons, viol, sol, 0) );

   return SCIP_OKAY;
}

/** gives the unique index of an expression constraint
 *
 * Each expression constraint gets an index assigned when it is created.
 * This index never changes and is unique among all expression constraints
 * within the same SCIP instance.
 * Thus, it can be used to sort a set of expression constraints.
 */
int SCIPgetConsExprIndex(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->consindex;
}

/** compares two expression constraints by its index
 *
 * Usable as compare operator in array sort functions.
 */
int SCIPcompareConsExprIndex(
   void*                 cons1,
   void*                 cons2
   )
{
   return SCIPgetConsExprIndex((SCIP_CONS*)cons1) - SCIPgetConsExprIndex((SCIP_CONS*)cons2);
}

/** returns an equivalent linear constraint if possible */
SCIP_RETCODE SCIPgetLinearConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_CONS**           lincons             /**< buffer to store linear constraint data */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;
   SCIP_CONSEXPR_EXPRHDLR* varhdlr;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR** vars;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(lincons != NULL);

   *lincons = NULL;
   expr = SCIPgetExprConsExpr(scip, cons);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   assert(sumhdlr != NULL);
   varhdlr = SCIPgetConsExprExprHdlrVar(conshdlr);
   assert(varhdlr != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* not a linear constraint if the root expression is not a sum */
   if( expr == NULL || expr->exprhdlr != sumhdlr )
      return SCIP_OKAY;

   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];

      /* at least one child is not a variable -> not a linear constraint */
      if( child->exprhdlr != varhdlr )
         return SCIP_OKAY;
   }

   /* collect all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPgetConsExprExprNChildren(expr)) );
   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];

      assert(child->exprhdlr == varhdlr);
      vars[i] = SCIPgetConsExprExprVarVar(child);
   }

   /* consider constant part of the sum expression */
   lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : (consdata->lhs - SCIPgetConsExprExprSumConstant(expr));
   rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : (consdata->rhs - SCIPgetConsExprExprSumConstant(expr));

   SCIP_CALL( SCIPcreateConsLinear(scip, lincons, SCIPconsGetName(cons),
         SCIPgetConsExprExprNChildren(expr), vars, SCIPgetConsExprExprSumCoefs(expr),
         lhs, rhs,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   /* free memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
SCIP_RETCODE SCIPgetLinvarMayDecreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   consdataFindUnlockedLinearVar(scip, conshdlr, consdata);

   *var = consdata->linvardecr;
   *coef = consdata->linvardecrcoef;

   return SCIP_OKAY;
}

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
SCIP_RETCODE SCIPgetLinvarMayIncreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   consdataFindUnlockedLinearVar(scip, conshdlr, consdata);

   *var = consdata->linvarincr;
   *coef = consdata->linvarincrcoef;

   return SCIP_OKAY;
}

/** detects nonlinear handlers that can handle the expressions and creates needed auxiliary variables
 *
 *  @note this method is only used for testing purposes
 */
SCIP_RETCODE SCIPdetectConsExprNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected while creating the auxiliary vars */
   )
{
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(infeasible != NULL);

   SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss, infeasible) );

   return SCIP_OKAY;
}

/** creates the nonlinearity handler and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprNlhdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLR**      nlhdlr,       /**< buffer where to store nonlinear handler */
   const char*                 name,         /**< name of nonlinear handler (must not be NULL) */
   const char*                 desc,         /**< description of nonlinear handler (can be NULL) */
   int                         priority,     /**< priority of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRDETECT((*detect)), /**< structure detection callback of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLREVALAUX((*evalaux)), /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_CONSEXPR_NLHDLRDATA*   data          /**< data of nonlinear handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nlhdlr != NULL);
   assert(name != NULL);
   assert(detect != NULL);
   assert(evalaux != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, nlhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*nlhdlr)->priority = priority;
   (*nlhdlr)->data = data;
   (*nlhdlr)->detect = detect;
   (*nlhdlr)->evalaux = evalaux;

   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->detecttime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->enfotime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->proptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->intevaltime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->reformulatetime) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/expr/nlhdlr/%s/enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname, "should this nonlinear handler be used",
      &(*nlhdlr)->enabled, FALSE, TRUE, NULL, NULL) );

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize, conshdlrdata->nnlhdlrs+1);

   conshdlrdata->nlhdlrs[conshdlrdata->nnlhdlrs] = *nlhdlr;
   ++conshdlrdata->nnlhdlrs;

   /* sort nonlinear handlers by priority, in decreasing order
    * will happen in INIT, so only do when called late
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_INIT && conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, nlhdlrCmp, conshdlrdata->nnlhdlrs);

   return SCIP_OKAY;
}

/** set the copy handler callback of a nonlinear handler */
void SCIPsetConsExprNlhdlrCopyHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR((*copy)) /**< copy callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->copyhdlr = copy;
}

/** set the nonlinear handler callback to free the nonlinear handler data */
void SCIPsetConsExprNlhdlrFreeHdlrData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freehdlrdata = freehdlrdata;
}

/** set the expression handler callback to free expression specific data of nonlinear handler */
void SCIPsetConsExprNlhdlrFreeExprData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freeexprdata = freeexprdata;
}

/** set the initialization and deinitialization callback of a nonlinear handler */
void SCIPsetConsExprNlhdlrInitExit(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINIT((*init)),   /**< initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXIT((*exit_))    /**< deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->init = init;
   nlhdlr->exit = exit_;
}

/** set the reformulate callback of a nonlinear handler */
void SCIPsetConsExprNlhdlrReformulate(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE((*reformulate)) /**< reformulation callback */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->reformulate = reformulate;
}

/** set the propagation callbacks of a nonlinear handler */
void SCIPsetConsExprNlhdlrProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINTEVAL((*inteval)), /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->inteval = inteval;
   nlhdlr->reverseprop = reverseprop;
}

/** set the separation callbacks of a nonlinear handler */
void SCIPsetConsExprNlhdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRENFO((*enfo)),         /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRESTIMATE((*estimate)), /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);
   assert(enfo != NULL || estimate != NULL);

   nlhdlr->initsepa = initsepa;
   nlhdlr->enfo = enfo;
   nlhdlr->estimate = estimate;
   nlhdlr->exitsepa = exitsepa;
}

/** gives name of nonlinear handler */
const char* SCIPgetConsExprNlhdlrName(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->name;
}

/** gives description of nonlinear handler, can be NULL */
const char* SCIPgetConsExprNlhdlrDesc(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->desc;
}

/** gives priority of nonlinear handler */
int SCIPgetConsExprNlhdlrPriority(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->priority;
}

/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_CONSEXPR_NLHDLR* SCIPfindConsExprNlhdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of nonlinear handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), name) == 0 )
         return conshdlrdata->nlhdlrs[h];

   return NULL;
}

/** gives handler data of nonlinear handler */
SCIP_CONSEXPR_NLHDLRDATA* SCIPgetConsExprNlhdlrData(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->data;
}

/** gives nonlinear handler expression data
 *
 * @return NULL if expr has not been detected by nlhdlr or nlhdlr did not store data
 */
SCIP_CONSEXPR_NLHDLREXPRDATA* SCIPgetConsExprNlhdlrExprData(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_CONSEXPR_EXPR*        expr           /**< expression */
   )
{
   int e;

   for( e = 0; e < expr->nenfos; ++e )
   {
      if( expr->enfos[e]->nlhdlr == nlhdlr )
         return expr->enfos[e]->nlhdlrexprdata;
   }

   return NULL;
}

/** returns whether nonlinear handler implements the reformulation callback */
SCIP_Bool SCIPhasConsExprNlhdlrReformulate(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->reformulate != NULL;
}

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_Bool SCIPhasConsExprNlhdlrInteval(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->inteval != NULL;
}

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_Bool SCIPhasConsExprNlhdlrReverseProp(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->reverseprop != NULL;
}

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_Bool SCIPhasConsExprNlhdlrInitSepa(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->initsepa != NULL;
}

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_Bool SCIPhasConsExprNlhdlrExitSepa(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->exitsepa != NULL;
}

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_Bool SCIPhasConsExprNlhdlrEnfo(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->enfo != NULL;
}

/** returns whether nonlinear handler implements the estimator callback */
SCIP_Bool SCIPhasConsExprNlhdlrEstimate(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->estimate != NULL;
}

/** call the detect callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRDETECT(SCIPdetectConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->detect != NULL);
   assert(nlhdlr->detecttime != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->detecttime) );
   SCIP_CALL( nlhdlr->detect(scip, conshdlr, nlhdlr, expr, cons, enforcemethods, enforcedbelow, enforcedabove, success, nlhdlrexprdata) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->detecttime) );

   if( *success )
      ++nlhdlr->ndetections;

   return SCIP_OKAY;
}

/** calls the reformulation callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(SCIPreformulateConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->reformulatetime != NULL);

   *refexpr = NULL;

   if( nlhdlr->reformulate == NULL )
      return SCIP_OKAY;

   /* call reformulation callback */
   SCIP_CALL( SCIPstartClock(scip, nlhdlr->reformulatetime) );
   SCIP_CALL( nlhdlr->reformulate(scip, conshdlr, nlhdlr, expr, refexpr) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->reformulatetime) );

   /* check whether reformulation was successful */
   if( *refexpr != NULL && *refexpr != expr )
      ++nlhdlr->nreformulates;

   return SCIP_OKAY;
}

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(SCIPevalauxConsExprNlhdlr)
{
   assert(nlhdlr != NULL);
   assert(nlhdlr->evalaux != NULL);

   SCIP_CALL( nlhdlr->evalaux(scip, nlhdlr, expr, nlhdlrexprdata, auxvalue, sol) );

   return SCIP_OKAY;
}

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(SCIPintevalConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->intevaltime != NULL);

   if( nlhdlr->inteval != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->intevaltime) );
      SCIP_CALL( nlhdlr->inteval(scip, nlhdlr, expr, nlhdlrexprdata, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->intevaltime) );

      ++nlhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(SCIPreversepropConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->proptime != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   if( nlhdlr->reverseprop == NULL )
   {
      *infeasible = FALSE;
      *nreductions = 0;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->proptime) );
   SCIP_CALL( nlhdlr->reverseprop(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, reversepropqueue, infeasible, nreductions, force) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->proptime) );

   /* update statistics */
   nlhdlr->ndomreds += *nreductions;
   if( *infeasible )
      ++nlhdlr->ncutoffs;
   ++nlhdlr->npropcalls;

   return SCIP_OKAY;
}

/** calls the separation initialization callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(SCIPinitsepaConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(infeasible != NULL);

   if( nlhdlr->initsepa == NULL )
   {
      *infeasible = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->initsepa(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, overestimate, underestimate, infeasible) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   ++nlhdlr->nenfocalls;
   if( *infeasible )
      ++nlhdlr->ncutoffs;

   return SCIP_OKAY;
}

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(SCIPexitsepaConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);

   if( nlhdlr->exitsepa != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( nlhdlr->exitsepa(scip, nlhdlr, expr, nlhdlrexprdata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );
   }

   return SCIP_OKAY;
}

/** calls the enforcement callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRENFO(SCIPenfoConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(result != NULL);

   if( nlhdlr->enfo == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->enfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, allowweakcuts, separated, addbranchscores, result) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;
   switch( *result )
   {
      case SCIP_SEPARATED :
         ++nlhdlr->nseparated;
         break;
      case SCIP_BRANCHED:
         ++nlhdlr->nbranchscores;
         break;
      case SCIP_CUTOFF:
         ++nlhdlr->ncutoffs;
         break;
      default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** calls the estimator callback of a nonlinear handler */
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(SCIPestimateConsExprNlhdlr)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(success != NULL);
   assert(addedbranchscores != NULL);

   if( nlhdlr->estimate == NULL )
   {
      *success = FALSE;
      *addedbranchscores = FALSE;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->estimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, targetvalue, rowprep, success, addbranchscores, addedbranchscores) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;

   return SCIP_OKAY;
}

/* computes a facet of the convex or concave envelope of a vertex polyhedral function
 * see (doxygen-)comment of this function in cons_expr.h
 * (this is by intention not a doxygen comment)
 */
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   SCIP_Real* corner;
   SCIP_Real* funvals;
   int* nonfixedpos;
   SCIP_Real maxfaceterror;
   int nvars; /* number of nonfixed variables */
   unsigned int ncorners;
   unsigned int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(function != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* identify fixed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &nonfixedpos, nallvars) );
   nvars = 0;
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         continue;
      nonfixedpos[nvars] = j;
      nvars++;
   }

   /* if all variables are fixed, then we could provide something trivial, but that wouldn't be the job of separation
    * if too many variables are not fixed, then we do nothing currently
    */
   if( nvars == 0 || nvars > SCIP_MAXVERTEXPOLYDIM )
   {
      SCIPwarningMessage(scip, "SCIPcomputeFacetVertexPolyhedral() called with %d nonfixed variables. Must be between [1,%d].\n", nvars, SCIP_MAXVERTEXPOLYDIM);
      SCIPfreeBufferArray(scip, &nonfixedpos);
      return SCIP_OKAY;
   }

   /* compute f(v^i) for each corner v^i of [l,u] */
   ncorners = POWEROFTWO(nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &funvals, ncorners) );
   SCIP_CALL( SCIPallocBufferArray(scip, &corner, nallvars) );
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         corner[j] = (box[2 * j] + box[2 * j + 1]) / 2.0;
   }
   for( i = 0; i < ncorners; ++i )
   {
      SCIPdebugMsg(scip, "corner %d: ", i);
      for( j = 0; j < nvars; ++j )
      {
         int varpos = nonfixedpos[j];
         /* if j'th bit of row index i is set, then take upper bound on var j, otherwise lower bound var j
          * we check this by shifting i for j positions to the right and checking whether the last bit is set
          */
         if( (i >> j) & 0x1 )
            corner[varpos] = box[2 * varpos + 1]; /* ub of var */
         else
            corner[varpos] = box[2 * varpos ]; /* lb of var */
         SCIPdebugMsgPrint(scip, "%g, ", corner[varpos]);
         assert(!SCIPisInfinity(scip, REALABS(corner[varpos])));
      }

      funvals[i] = function(corner, nallvars, fundata);

      SCIPdebugMsgPrint(scip, "obj = %e\n", funvals[i]);

      if( funvals[i] == SCIP_INVALID || SCIPisInfinity(scip, REALABS(funvals[i])) )  /*lint !e777*/
      {
         SCIPdebugMsg(scip, "cannot compute underestimator; function value at corner is too large %g\n", funvals[i]);
         goto CLEANUP;
      }
   }

   /* clear coefs array; below we only fill in coefs for nonfixed variables */
   BMSclearMemoryArray(facetcoefs, nallvars);

   if( nvars == 1 )
   {
      SCIP_CALL( computeVertexPolyhedralFacetUnivariate(scip, box[2 * nonfixedpos[0]], box[2 * nonfixedpos[0] + 1], funvals[0], funvals[1], success, &facetcoefs[nonfixedpos[0]], facetconstant) );

      /* check whether target has been missed */
      if( *success && overestimate == (*facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]] > targetvalue) )
      {
         SCIPdebugMsg(scip, "computed secant, but missed target %g (facetvalue=%g, overestimate=%d)\n", targetvalue, *facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]], overestimate);
         *success = FALSE;
      }
   }
   else if( nvars == 2 )
   {
      int idx1 = nonfixedpos[0];
      int idx2 = nonfixedpos[1];
      double p1[2] = { box[2*idx1],   box[2*idx2]   }; /* corner 0: 0>>0 & 0x1 = 0, 0>>1 & 0x1 = 0 */
      double p2[2] = { box[2*idx1+1], box[2*idx2]   }; /* corner 1: 1>>0 & 0x1 = 1, 1>>1 & 0x1 = 0 */
      double p3[2] = { box[2*idx1],   box[2*idx2+1] }; /* corner 2: 2>>0 & 0x1 = 0, 2>>1 & 0x1 = 1 */
      double p4[2] = { box[2*idx1+1], box[2*idx2+1] }; /* corner 3: 3>>0 & 0x1 = 1, 3>>1 & 0x1 = 1 */
      double xstar2[2] = { xstar[idx1], xstar[idx2] };
      double coefs[2] = { 0.0, 0.0 };

      SCIP_CALL( computeVertexPolyhedralFacetBivariate(scip, overestimate, p1, p2, p3, p4, funvals[0], funvals[1], funvals[2], funvals[3], xstar2, targetvalue, success, coefs, facetconstant) );

      facetcoefs[idx1] = coefs[0];
      facetcoefs[idx2] = coefs[1];
   }
   else
   {
      SCIP_CALL( computeVertexPolyhedralFacetLP(scip, conshdlr, overestimate, xstar, box, nallvars, nonfixedpos, funvals, nvars, targetvalue, success, facetcoefs, facetconstant) );
   }
   if( !*success )
   {
      SCIPdebugMsg(scip, "no success computing facet, %d vars\n", nvars);
      goto CLEANUP;
   }

   /*
    *  check and adjust facet with the algorithm of Rikun et al.
    */

   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, overestimate, funvals, box, nallvars, nvars, nonfixedpos, facetcoefs, *facetconstant);

   /* adjust constant part of the facet by maxerror to make it a valid over/underestimator (not facet though) */
   if( maxfaceterror > 0.0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real midval;
      SCIP_Real feastol;

      feastol = SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetLPFeastol(scip) : SCIPfeastol(scip);

      /* evaluate function in middle point to get some idea for a scaling */
      for( j = 0; j < nvars; ++j )
         corner[nonfixedpos[j]] = (box[2 * nonfixedpos[j]] + box[2 * nonfixedpos[j] + 1]) / 2.0;
      midval = function(corner, nallvars, fundata);
      if( midval == SCIP_INVALID )  /*lint !e777*/
         midval = 1.0;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* there seem to be numerical problems if the error is too large; in this case we reject the facet */
      if( maxfaceterror > conshdlrdata->vp_adjfacetthreshold * feastol * fabs(midval) )
      {
         SCIPdebugMsg(scip, "ignoring facet due to instability, it cuts off a vertex by %g (midval=%g).\n", maxfaceterror, midval);
         *success = FALSE;
         goto CLEANUP;
      }

      SCIPdebugMsg(scip, "maximum facet error %g (midval=%g), adjust constant to make cut valid!\n", maxfaceterror, midval);

      if( overestimate )
         *facetconstant += maxfaceterror;
      else
         *facetconstant -= maxfaceterror;
   }

   /* if we made it until here, then we have a nice facet */
   assert(*success);

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &funvals);
   SCIPfreeBufferArray(scip, &corner);
   SCIPfreeBufferArray(scip, &nonfixedpos);

   return SCIP_OKAY;
}
